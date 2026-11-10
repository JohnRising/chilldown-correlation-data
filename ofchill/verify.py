"""Verification helpers: expression-vs-numpy closure check, fit errors,
reproduction of prior runs."""

from __future__ import annotations

import json
import math

import numpy as np

from .closure import FALLBACK_H_FLOOR, CustomClosure


def tiny_network(closure, *, inlet: dict, P_in=3.0e5, P_out=2.95e5, mdot=0.01, Tf0=200.0, Tw0=250.0,
                 D=0.01168, steps=8, dt=0.05, fluid="Nitrogen") -> dict:
    """1 internal fluid node + 1 solid wall node + 1 convection conductor."""
    A = math.pi / 4 * D * D
    inl = {"id": "in", "type": "boundary", "x": 0, "y": 0, "pressure": P_in, **inlet}
    return {
        "meta": {"name": "tiny closure test", "version": 2},
        "settings": {"mode": "transient", "tolerance": 1e-5, "maxIterations": 200, "relaxation": 0.7,
                     "endTime": steps * dt, "dt": dt, "timeStepping": "fixed"},
        "fluid": {"model": "realFluid", "params": {"fluidName": fluid}},
        "nodes": [
            inl,
            {"id": "f", "type": "internal", "x": 100, "y": 0, "pressure": 0.5 * (P_in + P_out), "temperature": Tf0,
             "volume": 2e-6},
            {"id": "out", "type": "boundary", "x": 200, "y": 0, "pressure": P_out, "temperature": 300.0},
        ],
        "solidNodes": [{"id": "w", "type": "solid", "x": 100, "y": 80, "temperature": Tw0, "mass": 0.02,
                        "cp": {"material": "stainless-steel-304"}}],
        "conductors": [{"id": "c", "from": "f", "to": "w",
                        "type": {"kind": "convection", "area": 2e-3, "correlation": closure.correlation(D, A)}}],
        "branches": [
            {"id": "src", "from": "in", "to": "f", "component": {"type": "flowSource", "massFlow": mdot}},
            {"id": "ex", "from": "f", "to": "out",
             "component": {"type": "pipe", "length": 0.05, "diameter": D, "roughness": 1.5e-6}},
        ],
    }


def scope_from_result(raw: dict, fluid_node: str, wall_node: str, branch_ids: list[str], flow_area: float, D: float,
                      area: float, fluid_name: str | None = None) -> dict[str, np.ndarray]:
    """Rebuild the custom-h scope at each recorded step from the result."""
    nd = raw["nodes"][fluid_node]
    f = lambda v: np.array([np.nan if x is None else x for x in v], float)  # noqa: E731
    G = 0.5 * sum(np.abs(f(raw["branches"][b]["mdot"])) for b in branch_ids) / flow_area
    sc = dict(
        t=np.asarray(raw["times"], float), Tf=f(nd["temperature"]), P=f(nd["pressure"]), rho=f(nd["density"]),
        mu=f(nd["viscosity"]), k=f(nd["thermalConductivity"]), cp=f(nd["specificHeat"]),
        quality=f(nd.get("quality", [None] * len(raw["times"]))), G=G, D=np.full_like(G, D),
        flowArea=np.full_like(G, flow_area), area=np.full_like(G, area),
    )
    if wall_node:
        sc["Tw"] = f(raw["solidNodes"][wall_node]["temperature"])
    # In the dome the RESULT reports thermalConductivity = NaN, but the custom
    # scope carries the mixture k = 1/(x/kg + (1-x)/kf) (realFluid.domePH).
    ph = nd.get("phase") or []
    fluid = fluid_name or "Nitrogen"
    import CoolProp.CoolProp as CP

    for i, p in enumerate(ph):
        if p == "twoPhase" and not np.isfinite(sc["k"][i]):
            x = sc["quality"][i]
            kf = CP.PropsSI("L", "P", sc["P"][i], "Q", 0, fluid)
            kg = CP.PropsSI("L", "P", sc["P"][i], "Q", 1, fluid)
            sc["k"][i] = 1.0 / (x / kg + (1 - x) / kf)
    sc["Pr"] = sc["cp"] * sc["mu"] / sc["k"]
    sc["Re"] = sc["G"] * D / sc["mu"]
    return sc


def closure_vs_numpy(closure: CustomClosure, inlet: dict, **kw) -> dict:
    """Solve the tiny network with ``closure`` and compare the reported
    heatTransferCoeff (a FRESH, un-relaxed evaluation at each recorded state,
    floored at 5 W/m^2K) with closure.evaluate on the reported state."""
    import openflume as of

    cfg = tiny_network(closure, inlet=inlet, **kw)
    res = of.Network.from_dict(cfg).solve()
    raw = res.raw
    D = cfg["conductors"][0]["type"]["correlation"]["diameter"]
    A = cfg["conductors"][0]["type"]["correlation"]["flowArea"]
    sc = scope_from_result(raw, "f", "w", ["src", "ex"], A, D, cfg["conductors"][0]["type"]["area"], closure.fluid)
    h_of = np.asarray(raw["conductors"]["c"]["heatTransferCoeff"], float)
    h_np = np.maximum(closure.evaluate(**sc), FALLBACK_H_FLOOR)
    rel = np.abs(h_np / h_of - 1)
    return dict(converged=res.converged, phase=raw["nodes"]["f"].get("phase"), h_openflume=h_of.tolist(),
                h_numpy=h_np.tolist(), max_rel_err=float(np.max(rel)), expr_len=len(closure.expression))


def closure_vs_numpy_uf(closure: CustomClosure, case_no=8, end_time=2.0, N=10, dt=0.1, **opts) -> dict:
    """Same check on a short UF network run (all conductors, all steps;
    covers vapour, two-phase and liquid nodes)."""
    import openflume as of

    from .builder import build_uf_case

    m = build_uf_case(case_no, N=N, dt=dt, closure=closure, end_time=end_time, **opts)
    raw = of.Network.from_dict(m.config).solve().raw
    errs, phases = [], set()
    for j, (fid, sid, cid) in enumerate(zip(m.fluid_ids, m.wall_ids, m.conv_ids)):
        brs = [b["id"] for b in m.config["branches"] if fid in (b["from"], b["to"])]
        cond = next(c for c in m.config["conductors"] if c["id"] == cid)
        sc = scope_from_result(raw, fid, sid, brs, m.flow_area, m.D, cond["type"]["area"])
        h_of = np.asarray(raw["conductors"][cid]["heatTransferCoeff"], float)
        zp = {"z": float(m.wall_x[j])} if "z" in closure.params else None
        h_np = np.maximum(closure.evaluate(zp, **sc), FALLBACK_H_FLOOR)
        errs.append(np.abs(h_np / h_of - 1))
        phases |= set(raw["nodes"][fid].get("phase", []))
    e = np.concatenate(errs)
    hs = np.concatenate([np.asarray(raw["conductors"][c]["heatTransferCoeff"], float) for c in m.conv_ids])
    return dict(h_range=[float(hs.min()), float(hs.max())], max_rel_err=float(e.max()), p99_rel_err=float(np.quantile(e, 0.99)), n=int(e.size),
                phases=sorted(phases))
