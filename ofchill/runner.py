"""Parallel runner with result caching.

A job is (run_id, ChilldownModel, spec_metadata).  The config is built in the
parent and shipped (as JSON) to a worker process; each worker keeps ONE
openflume session (the default session) for its lifetime.  Results are cached
by a hash of the full NetworkConfig (+ RESULT_VERSION + openflume version) in
results/_cache/<hash>.{npz,json}; results/<tag>/index.json maps run ids to
hashes and spec metadata.
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
import threading
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import numpy as np

from .data import ROOT

RESULTS = ROOT / "results"
CACHE = RESULTS / "_cache"
RESULT_VERSION = "ofchill-result-1"
STALL_SCALED = 0.01  # stepResidualsScaled above this = stalled step (schema.ts)


def _of_version() -> str:
    try:
        import openflume

        return getattr(openflume, "__version__", "?")
    except Exception:  # pragma: no cover
        return "?"


def config_hash(config: dict) -> str:
    s = json.dumps(config, sort_keys=True, separators=(",", ":")) + RESULT_VERSION + _of_version()
    return hashlib.sha1(s.encode()).hexdigest()[:16]


# ----------------------------------------------------------------------------
# wall material cp (energy diagnostic only; same NIST fits as the solver presets)
# ----------------------------------------------------------------------------
_NIST = {
    "stainless-steel-304": [22.0061, -127.5528, 303.647, -381.0098, 274.0328, -112.9212, 24.7593, -2.239153, 0],
    "ofhc-copper": [-1.91844, -0.15973, 8.61013, -18.996, 21.9661, -12.7328, 3.54322, -0.3797, 0],
}


def wall_cp(material: str | None, T: np.ndarray) -> np.ndarray | None:
    c = _NIST.get(material or "")
    if c is None:
        return None
    x = np.log10(np.asarray(T, float))
    return 10 ** np.polynomial.polynomial.polyval(x, c)


def wall_energy(material, mass, T1, T0) -> float | None:
    """mass * int_{T1}^{T0} cp dT (J), positive when cooled from T0 to T1."""
    Ts = np.linspace(T1, T0, 401)
    cp = wall_cp(material, Ts)
    if cp is None:
        return None
    return float(mass * np.trapezoid(cp, Ts))


# ----------------------------------------------------------------------------
# result extraction
# ----------------------------------------------------------------------------
def _arr(v) -> np.ndarray:
    return np.array([np.nan if (x is None or (isinstance(x, float) and not np.isfinite(x))) else x for x in v], float)


def extract(raw: dict, layout: dict) -> tuple[dict[str, np.ndarray], dict]:
    t = np.asarray(raw["times"], float)
    nodes, br, sol, cond = raw["nodes"], raw["branches"], raw.get("solidNodes") or {}, raw.get("conductors") or {}
    wx = np.asarray(layout["wall_x"])
    W = np.array([sol[s]["temperature"] for s in layout["wall_ids"]], float)
    arrs: dict[str, np.ndarray] = {
        "t": t, "wall_x": wx, "wall_T": W,
        "h": np.array([_arr(cond[c].get("heatTransferCoeff", [])) for c in layout["conv_ids"]]),
        "q": np.array([_arr(cond[c].get("heatFlux", [])) for c in layout["conv_ids"]]),
        "Q": np.array([_arr(cond[c]["heatRate"]) for c in layout["conv_ids"]]),
        "fluid_P": np.array([_arr(nodes[f]["pressure"]) for f in layout["fluid_ids"]]),
        "fluid_T": np.array([_arr(nodes[f]["temperature"]) for f in layout["fluid_ids"]]),
        "fluid_x": np.array([_arr(nodes[f].get("quality", [None] * len(t))) for f in layout["fluid_ids"]]),
        "fluid_rho": np.array([_arr(nodes[f]["density"]) for f in layout["fluid_ids"]]),
        "mdot_in": _arr(br[layout["inlet_branch"]]["mdot"]),
        "mdot_out": _arr(br[layout["outlet_branch"]]["mdot"]),
        "step_res_scaled": _arr(raw.get("stepResidualsScaled") or []),
    }
    for name, xs in (layout.get("stations") or {}).items():
        j = int(np.clip(np.searchsorted(wx, xs) - 1, 0, len(wx) - 2))
        w = (xs - wx[j]) / (wx[j + 1] - wx[j])
        arrs[f"T_{name}"] = (1 - w) * W[j] + w * W[j + 1]
    arrs["phase_code"] = np.array(
        [[{"liquid": 0, "twoPhase": 1, "vapor": 2}.get(p, 3) for p in nodes[f].get("phase", [])] for f in layout["fluid_ids"]],
        dtype=np.int8,
    )
    diag = energy_diag(raw, layout, arrs)
    return arrs, diag


def energy_diag(raw: dict, layout: dict, a: dict) -> dict:
    t = a["t"]
    n = len(t)
    if n < 2:
        return {}
    dt = np.diff(t)
    Q = a["Q"]  # fluid -> wall positive
    wall_heat = float(-np.sum(Q[:, 1:] * dt))
    W = a["wall_T"]
    mat = layout.get("wall_material")
    we = [wall_energy(mat, m, W[i, -1], W[i, 0]) for i, m in enumerate(layout["wall_mass"])]
    wall_e = None if any(v is None for v in we) else float(sum(we))
    nodes = raw["nodes"]
    mIn, mOut = a["mdot_in"], a["mdot_out"]
    hIn = _arr(nodes[layout["inlet_node"]]["enthalpy"])
    on = nodes[raw_outlet_from(layout)]
    Af = layout["flow_area"]
    vout = mOut / (_arr(on["density"]) * Af)
    hOut = _arr(on["enthalpy"]) + 0.5 * vout**2
    enth = float(np.sum((mOut[1:] * hOut[1:] - mIn[1:] * hIn[1:]) * dt))
    dU = 0.0
    m0 = m1 = 0.0
    for f, V in zip(layout["fluid_ids"], layout["fluid_volumes"]):
        rho = _arr(nodes[f]["density"])
        u = _arr(nodes[f]["internalEnergy"])
        dU += (rho[-1] * u[-1] - rho[0] * u[0]) * V
        m0 += rho[0] * V
        m1 += rho[-1] * V
    d = {
        "wallHeatJ": wall_heat, "wallEnergyRemovedJ": wall_e, "fluidEnthalpyGainJ": enth + dU,
        "massInKg": float(np.sum(mIn[1:] * dt)), "massOutKg": float(np.sum(mOut[1:] * dt)),
        "fluidInventoryChangeKg": m1 - m0,
    }
    if wall_e:
        d["wallHeat_over_wallEnergy_m1"] = wall_heat / wall_e - 1
        d["fluidGain_over_wallEnergy_m1"] = d["fluidEnthalpyGainJ"] / wall_e - 1
    d["mass_imbalance_kg"] = d["massInKg"] - d["massOutKg"] - d["fluidInventoryChangeKg"]
    return d


def raw_outlet_from(layout: dict) -> str:
    # the outlet branch runs from the last internal node
    return layout["fluid_ids"][-1]


# ----------------------------------------------------------------------------
# segmented solve (workaround for the openflume 0.3.1 CoolProp-wasm exception leak)
# ----------------------------------------------------------------------------
def _shift_sched(sched, t0):
    if not sched:
        return sched
    ts = [p[0] for p in sched]
    vs = [p[1] for p in sched]
    out = [[0.0, float(np.interp(t0, ts, vs))]] + [[p[0] - t0, p[1]] for p in sched if p[0] > t0]
    return out


def _segment_config(cfg: dict, t0: float, t1: float, state: dict | None) -> dict:
    c = json.loads(json.dumps(cfg))
    c["settings"]["endTime"] = t1 - t0
    for n in c["nodes"]:
        for k in ("pressureSchedule", "temperatureSchedule"):
            if k in n:
                n[k] = _shift_sched(n[k], t0)
                if k == "pressureSchedule":
                    n["pressure"] = n[k][0][1]
        if state and n["type"] == "internal" and n["id"] in state["nodes"]:
            st = state["nodes"][n["id"]]
            n["pressure"] = st["P"]
            n.pop("quality", None)
            if st.get("x") is not None:
                n["quality"] = st["x"]
                n.pop("temperature", None)
            else:
                n["temperature"] = st["T"]
    for sn in c.get("solidNodes") or []:
        if "temperatureSchedule" in sn:
            sn["temperatureSchedule"] = _shift_sched(sn["temperatureSchedule"], t0)
        if state and sn["id"] in state["solids"]:
            sn["temperature"] = state["solids"][sn["id"]]
    for b in c["branches"]:
        comp = b["component"]
        if "massFlowSchedule" in comp:
            comp["massFlowSchedule"] = _shift_sched(comp["massFlowSchedule"], t0)
            comp["massFlow"] = comp["massFlowSchedule"][0][1]
    return c


def _end_state(raw: dict) -> dict:
    st = {"nodes": {}, "solids": {}}
    for nid, nd in raw["nodes"].items():
        x = (nd.get("quality") or [None])[-1]
        ph = (nd.get("phase") or [None])[-1]
        st["nodes"][nid] = {"P": nd["pressure"][-1], "T": nd["temperature"][-1],
                            "x": x if (ph == "twoPhase" and x is not None and 0 < x < 1) else None}
    for sid, sd in (raw.get("solidNodes") or {}).items():
        st["solids"][sid] = sd["temperature"][-1]
    return st


def _merge(acc: dict | None, raw: dict, t0: float) -> dict:
    if acc is None:
        return raw
    acc["times"] += [t + t0 for t in raw["times"][1:]]
    for table in ("nodes", "branches", "solidNodes", "conductors"):
        for eid, vals in (raw.get(table) or {}).items():
            for q, series in vals.items():
                if isinstance(series, list) and q in acc[table][eid]:
                    acc[table][eid][q] += series[1:]
    for k in ("stepResiduals", "stepResidualsScaled"):
        if k in raw and k in acc:
            acc[k] += raw[k]
    acc["converged"] = bool(acc.get("converged")) and bool(raw.get("converged"))
    return acc


def solve_segmented(of, cfg: dict, seg: float, sess_factory, log=None) -> tuple[dict, int]:
    """Solve in segments of length seg, each in a FRESH session; on a core error retry the segment
    with seg/4 (twice).  Returns (merged raw, number of segments)."""
    tend = cfg["settings"]["endTime"]
    t0, acc, state, nseg = 0.0, None, None, 0
    while t0 < tend - 1e-9:
        length = seg
        for attempt in range(3):
            t1 = min(tend, t0 + length)
            sess = sess_factory()
            try:
                raw = of.Network.from_dict(_segment_config(cfg, t0, t1, state), session=sess).solve().raw
                break
            except Exception as e:  # noqa: BLE001
                if log:
                    log(f"segment {t0:.1f}-{t1:.1f} failed ({str(e)[:80]}); retry shorter")
                length /= 4
                if attempt == 2:
                    raise
            finally:
                try:
                    sess._kill()
                except Exception:
                    pass
        acc = _merge(acc, raw, t0)
        state = _end_state(raw)
        t0 = t1
        nseg += 1
    return acc, nseg


# ----------------------------------------------------------------------------
# worker
# ----------------------------------------------------------------------------
def _worker(job: dict) -> dict:
    """Runs in a worker process.  Never raises: failures are recorded."""
    import openflume as of
    from openflume._bridge import default_session

    rid, h = job["run_id"], job["hash"]
    out_npz, out_json = Path(job["npz"]), Path(job["json"])
    meta = {"run_id": rid, "hash": h, "spec": job["spec"], "openflume": _of_version(), "pid": os.getpid()}
    sess = default_session()
    timed_out = threading.Event()

    def _kill():
        timed_out.set()
        try:
            sess._kill()
        except Exception:
            pass

    t0 = time.time()
    raw = None
    try:
        sess.start()
        t_start = time.time()
        timer = threading.Timer(job["timeout"], _kill)
        timer.daemon = True
        timer.start()
        try:
            prog = {"t": 0.0}
            meta["progress"] = prog
            if job.get("segment_s"):
                from openflume._bridge import Session

                raw_m, nseg = solve_segmented(of, json.loads(job["config"]), job["segment_s"], Session)
                meta["n_segments"] = nseg

                class _R:
                    raw = raw_m

                res = _R()
            else:
                net = of.Network.from_dict(json.loads(job["config"]), session=sess)

                def _prog(p):
                    prog["t"] = p.get("time", prog["t"])

                res = net.solve(progress=_prog)
            raw = res.raw
        finally:
            timer.cancel()
        meta["solve_s"] = time.time() - t_start
        arrs, diag = extract(raw, job["layout"])
        sr = arrs["step_res_scaled"]
        stalled = np.where(sr > STALL_SCALED)[0].tolist() if sr.size else []
        meta.update(
            status="ok", converged=bool(raw.get("converged")), aborted=bool(raw.get("aborted", False)),
            n_steps=int(len(arrs["t"]) - 1), t_end=float(arrs["t"][-1]),
            stalled_steps=stalled, n_stalled=len(stalled),
            max_scaled_residual=float(np.nanmax(sr)) if sr.size else None,
            warnings=list(raw.get("warnings") or [])[:20], diag=diag,
            h_at_floor_frac=float(np.mean(arrs["h"] <= 5.0 + 1e-9)) if arrs["h"].size else None,
        )
        np.savez_compressed(out_npz, **arrs)
    except BaseException as e:  # noqa: BLE001
        # A CoolProp-wasm abort ("memory access out of bounds", "Aborted") poisons the Node module:
        # every later solve in this session would fail instantly.  Always restart the session.
        try:
            sess._kill()
        except Exception:
            pass
        meta.update(
            status="timeout" if timed_out.is_set() else "failed",
            error=f"{type(e).__name__}: {e}"[:2000],
            traceback=traceback.format_exc()[-4000:],
            solve_s=time.time() - t0,
        )
    if isinstance(meta.get("progress"), dict):
        meta["last_sim_t"] = meta.pop("progress")["t"]
    meta["wall_s"] = time.time() - t0
    out_json.write_text(json.dumps(meta, indent=1, default=float))
    return meta


# ----------------------------------------------------------------------------
# driver
# ----------------------------------------------------------------------------
def make_job(run_id: str, model, spec: dict, timeout: float = 900.0, segment_s: float | None = None) -> dict:
    cfg_json = json.dumps(model.config, sort_keys=True)
    h = config_hash(model.config) + (f"_seg{segment_s:g}" if segment_s else "")
    CACHE.mkdir(parents=True, exist_ok=True)
    return {
        "run_id": run_id, "hash": h, "config": cfg_json, "layout": model.layout(), "spec": spec,
        "timeout": timeout, "npz": str(CACHE / f"{h}.npz"), "json": str(CACHE / f"{h}.json"),
        "segment_s": segment_s,
    }


def is_cached(job: dict, retry_failed: bool = False) -> bool:
    j = Path(job["json"])
    if not j.exists():
        return False
    m = json.loads(j.read_text())
    if m.get("status") == "ok":
        return Path(job["npz"]).exists()
    return not retry_failed


def run_jobs(jobs: list[dict], tag: str, max_workers: int = 6, force: bool = False, retry_failed: bool = False,
             log=print) -> dict[str, dict]:
    """Run jobs (skipping cached ones) and write results/<tag>/index.json."""
    tagdir = RESULTS / tag
    tagdir.mkdir(parents=True, exist_ok=True)
    index_f = tagdir / "index.json"
    index = json.loads(index_f.read_text()) if index_f.exists() else {}
    todo = [j for j in jobs if force or not is_cached(j, retry_failed)]
    for j in jobs:
        index[j["run_id"]] = {"hash": j["hash"], "spec": j["spec"]}
    index_f.write_text(json.dumps(index, indent=1))
    log(f"[ofchill] {len(jobs)} jobs, {len(jobs) - len(todo)} cached, running {len(todo)} on {max_workers} workers")
    metas: dict[str, dict] = {}
    if todo:
        with ProcessPoolExecutor(max_workers=min(max_workers, len(todo))) as ex:
            futs = {ex.submit(_worker, j): j for j in todo}
            for f in as_completed(futs):
                j = futs[f]
                try:
                    m = f.result()
                except Exception as e:  # worker process died
                    m = {"run_id": j["run_id"], "hash": j["hash"], "spec": j["spec"], "status": "failed",
                         "error": f"worker crashed: {e}"}
                    Path(j["json"]).write_text(json.dumps(m, indent=1))
                metas[j["run_id"]] = m
                log(f"[ofchill] {m['run_id']}: {m.get('status')} conv={m.get('converged')} "
                    f"stalled={m.get('n_stalled')} solve={m.get('solve_s', 0):.1f}s "
                    f"{m.get('error', '')[:200]}")
                sys.stdout.flush()
    for j in jobs:
        if j["run_id"] not in metas:
            metas[j["run_id"]] = json.loads(Path(j["json"]).read_text())
    return metas


def load_run(tag: str, run_id: str) -> tuple[dict, dict[str, np.ndarray] | None]:
    index = json.loads((RESULTS / tag / "index.json").read_text())
    h = index[run_id]["hash"]
    meta = json.loads((CACHE / f"{h}.json").read_text())
    npz = CACHE / f"{h}.npz"
    arrs = dict(np.load(npz)) if npz.exists() else None
    return meta, arrs


def load_index(tag: str) -> dict:
    return json.loads((RESULTS / tag / "index.json").read_text())
