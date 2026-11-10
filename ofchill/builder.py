"""Network builders: generic pipe chilldown + the UF 2016 case builder.

Topology ('internal' layout with half-lumped end nodes; identical to the
reference TypeScript builder OpenFLUME-uf/scripts/uf/ufChilldown.ts):

  f0 (boundary: inlet P(t), quality or T)
   └─ pipe0  flowSource ṁ(t)       (or any component given as inlet_component)
  f1 .. f{N+1}  internal fluid nodes at x_j = j L/N  (end nodes half volume)
   └─ pipe1..pipeN   pipe segments of length L/N
  f{N+2} (boundary: exit P(t)) via pipe{N+1} (exit_length long)
  s0..sN    wall nodes at x_j (end nodes half mass), conv_j: f{j+1} <-> s_j,
            cond_j: s_j <-> s_{j+1} (axial conduction)
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Callable, Sequence

import numpy as np

from .closure import Closure
from .data import UF_GEOM, boxcar_schedule, load_uf_case


@dataclass
class ChilldownModel:
    config: dict
    N: int
    L: float
    D: float
    wall_x: np.ndarray
    wall_mass: np.ndarray
    wall_ids: list[str]
    fluid_ids: list[str]  # internal fluid node ids aligned with wall_ids
    conv_ids: list[str]
    inlet_branch: str
    outlet_branch: str
    inlet_node: str
    outlet_node: str
    flow_area: float
    stations: dict[str, float] = field(default_factory=dict)
    wall_cp: Callable[[np.ndarray], np.ndarray] | None = None
    info: dict = field(default_factory=dict)

    def layout(self) -> dict:
        """JSON-serialisable layout (for the runner)."""
        return dict(
            N=self.N, L=self.L, D=self.D, wall_x=self.wall_x.tolist(), wall_mass=self.wall_mass.tolist(),
            wall_ids=self.wall_ids, fluid_ids=self.fluid_ids, conv_ids=self.conv_ids,
            inlet_branch=self.inlet_branch, outlet_branch=self.outlet_branch, inlet_node=self.inlet_node,
            outlet_node=self.outlet_node, flow_area=self.flow_area, stations=self.stations,
            fluid_volumes=self.info.get("fluid_volumes"), wall_material=self.info.get("wall_material"), info=self.info,
        )


def build_pipe_chilldown(
    *,
    L: float,
    D: float,
    closure: Closure,
    fluid: str,
    wall_T0: float | Callable[[float], float] | Sequence[float],
    p_in_schedule: list[list[float]],
    p_out_schedule: list[list[float]],
    mdot_schedule: list[list[float]] | None = None,
    OD: float | None = None,
    wall_area: float | None = None,
    wall_rho: float = 8027.0,
    wall_material: str | dict = "stainless-steel-304",
    wall_k: dict | None = None,
    N: int = 10,
    dt: float = 0.1,
    end_time: float | None = None,
    inlet_quality: float | None = 0.0,
    inlet_temperature: float | None = None,
    inertia: bool = True,
    roughness: float = 1.5e-6,
    inlet_length: float = 0.05,
    exit_length: float = 0.05,
    exit_temperature: float = 300.0,
    fluid_T0: float | Callable[[float], float] | None = None,
    inlet_component: dict | None = None,
    stations: dict[str, float] | None = None,
    settings: dict | None = None,
    name: str = "pipe chilldown",
    time_stepping: str = "fixed",
    parasitic: tuple[float, float, float] | None = None,
    axial_conduction: bool = True,
    friction_factor: float | None = None,
) -> ChilldownModel:
    """Generic uniform-tube chilldown network (see module docstring).

    wall_material: OpenFLUME material preset name ('stainless-steel-304',
    'ofhc-copper', 'aluminum-6061-t6', ...) or a SolidPropertySpec dict used
    for both cp and k (override k with ``wall_k``).  Wall cross-section from
    OD or ``wall_area``.  wall_T0 / fluid_T0: constant, f(x) or per-node list.
    Schedules are [[t, value], ...] in SI (Pa, kg/s).  Without mdot_schedule
    pass ``inlet_component`` (any branch component dict) instead.
    parasitic = (q0, h_par, T_amb): wall heat leak q0 + h_par (T_amb - Tw) [W/m2 of inner area]
    (constant heatInput + a conduction link to an ambient node; CLOSURE_SPEC §9: (98.6, 2.2, 293)).
    """
    A_f = math.pi / 4 * D * D
    if wall_area is None:
        if OD is None:
            raise ValueError("give OD or wall_area")
        wall_area = math.pi / 4 * (OD * OD - D * D)
    segL = L / N
    vol = A_f * segL
    mass_seg = wall_rho * wall_area * segL
    conv_area = math.pi * D * segL
    if end_time is None:
        end_time = float(p_in_schedule[-1][0])
    wall_cp = {"material": wall_material} if isinstance(wall_material, str) else wall_material
    wall_kspec = wall_k or wall_cp
    pos = lambda x: {"x": x, "y": 0, "z": 0}  # noqa: E731

    def _prof(f, j, x):
        if f is None:
            return None
        if callable(f):
            return float(f(x))
        if np.ndim(f) == 0:
            return float(f)
        return float(f[j])

    P0in = p_in_schedule[0][1]
    P0out = p_out_schedule[0][1]
    inlet = {
        "id": "f0", "type": "boundary", "x": 0, "y": 0, "position": pos(-inlet_length),
        "pressure": P0in, "pressureSchedule": p_in_schedule, "label": "Inlet",
    }
    if inlet_temperature is not None:
        inlet["temperature"] = inlet_temperature
    else:
        inlet["quality"] = inlet_quality if inlet_quality is not None else 0.0
    nodes: list[dict] = [inlet]
    solids, conductors, branches = [], [], []
    wall_x, wall_mass, wall_ids, fluid_ids, conv_ids, fvols = [], [], [], [], [], []
    for j in range(N + 1):
        frac = 0.5 if j in (0, N) else 1.0
        x = j * segL
        Tw0 = _prof(wall_T0, j, x)
        Tf0 = _prof(fluid_T0, j, x) if fluid_T0 is not None else Tw0
        fid, sid, cid = f"f{j + 1}", f"s{j}", f"conv{j}"
        nodes.append({
            "id": fid, "type": "internal", "x": (j + 1) * 170, "y": 0, "position": pos(x),
            "pressure": P0in + (P0out - P0in) * j / N, "temperature": Tf0, "volume": vol * frac,
            "label": f"Node {j}",
        })
        solids.append({
            "id": sid, "type": "solid", "x": (j + 1) * 170, "y": 80, "position": pos(x),
            "temperature": Tw0, "mass": mass_seg * frac, "cp": wall_cp, "label": f"Wall {j}",
        })
        conductors.append({
            "id": cid, "from": fid, "to": sid,
            "type": {"kind": "convection", "area": conv_area * frac,
                     "correlation": closure.correlation(D, A_f, x=x, seg_len=segL * frac)},
            "label": f"Conv {j}",
        })
        wall_x.append(x); wall_mass.append(mass_seg * frac); wall_ids.append(sid); fluid_ids.append(fid)
        conv_ids.append(cid); fvols.append(vol * frac)
    out_id = f"f{N + 2}"
    nodes.append({
        "id": out_id, "type": "boundary", "x": (N + 2) * 170, "y": 0, "position": pos(L + exit_length),
        "pressure": P0out, "temperature": exit_temperature, "pressureSchedule": p_out_schedule, "label": "Exit",
    })
    if parasitic is not None:
        q0, hpar, Tamb = parasitic
        solids.append({"id": "amb", "type": "ambient", "x": 0, "y": 200, "temperature": Tamb,
                       "label": "Ambient (parasitic heat)"})
        for j, sn in enumerate(solids[:N + 1]):
            a_in = conv_area * (0.5 if j in (0, N) else 1.0)
            sn["heatInput"] = q0 * a_in
            conductors.append({"id": f"par{j}", "from": "amb", "to": sn["id"],
                               "type": {"kind": "conduction", "k": hpar, "area": a_in, "length": 1.0},
                               "label": f"Parasitic {j}"})
    for j in range(N if axial_conduction else 0):
        conductors.append({
            "id": f"cond{j}", "from": f"s{j}", "to": f"s{j + 1}",
            "type": {"kind": "conduction", "k": wall_kspec, "area": wall_area, "length": segL},
            "label": f"Cond {j}",
        })
    def pipe(ln):
        c = {"type": "pipe", "length": ln, "diameter": D, "roughness": roughness}
        if friction_factor is not None:
            c["frictionFactor"] = friction_factor
        return c
    if inlet_component is None:
        if mdot_schedule is None:
            raise ValueError("give mdot_schedule or inlet_component")
        inlet_component = {"type": "flowSource", "massFlow": mdot_schedule[0][1], "massFlowSchedule": mdot_schedule}
    branches.append({"id": "pipe0", "from": "f0", "to": "f1", "component": inlet_component, "label": "Inlet"})
    for j in range(1, N + 1):
        branches.append({"id": f"pipe{j}", "from": f"f{j}", "to": f"f{j + 1}", "component": pipe(segL),
                         "label": f"Pipe {j}"})
    branches.append({"id": f"pipe{N + 1}", "from": f"f{N + 1}", "to": out_id, "component": pipe(exit_length),
                     "label": "Exit branch"})
    st = {
        "mode": "transient", "tolerance": 1e-5, "maxIterations": 200, "relaxation": 0.7,
        "endTime": end_time, "dt": dt, "timeStepping": time_stepping,
    }
    if inertia:
        st.update(momentumFlux=True, kineticEnergy=True)
    if settings:
        st.update(settings)
    config = {
        "meta": {"name": name, "version": 2},
        "settings": st,
        "fluid": {"model": "realFluid", "params": {"fluidName": fluid}},
        "nodes": nodes,
        "solidNodes": solids,
        "conductors": conductors,
        "branches": branches,
    }
    return ChilldownModel(
        config=config, N=N, L=L, D=D, wall_x=np.array(wall_x), wall_mass=np.array(wall_mass), wall_ids=wall_ids,
        fluid_ids=fluid_ids, conv_ids=conv_ids, inlet_branch="pipe0", outlet_branch=f"pipe{N + 1}",
        inlet_node="f0", outlet_node=out_id, flow_area=A_f, stations=dict(stations or {}),
        info=dict(fluid=fluid, fluid_volumes=fvols, wall_material=wall_material if isinstance(wall_material, str) else None,
                  closure=closure.to_spec(), name=name),
    )


# Case 4: P1/P2 collapse from ~130 kPa to 6-30 kPa at t = 18.8 s (end of the
# test; below the N2 triple point, CoolProp fails at the boundary) -> stop before.
UF_END_OVERRIDE = {4: 18.7}


def build_uf_case(
    case_no: int,
    N: int = 10,
    dt: float = 0.1,
    closure: Closure | None = None,
    *,
    g_shift_s: float = 0.0,
    inlet_quality: float = 0.0,
    exit_p_scale: float = 1.0,
    end_time: float | None = None,
    inertia: bool = True,
    knot: float = 0.1,
    half_window: float = 0.1,
    window_tol: float = 1e-6,
    data_json: str | None = None,
    settings: dict | None = None,
    **extra: Any,
) -> ChilldownModel:
    """UF 2016 low-g LN2 chilldown, case 1..10 (Darr et al. 2016).

    Inlet ṁ = A·G(t + g_shift_s) (flowSource), inlet boundary at P1(t) with
    quality ``inlet_quality``, exit boundary at exit_p_scale·P2(t);
    schedules are 0.2-s boxcar means at 0.1-s knots; initial wall/fluid T
    linear through the two station TC means at t = 0.
    """
    from .closures import get_closure

    if closure is None:
        closure = get_closure("miropolskii")
    if data_json is not None:
        # Exact-reproduction path: the prior agent's uf_case<N>.json (full-
        # precision xls times/values) with the TS boxcar (tol = 0).
        import json as _json

        dj = _json.loads(open(data_json).read())
        t = np.asarray(dj["t"], float)
        cols = {"G": np.asarray(dj["G"]), "P1": np.asarray(dj["P1"]), "P2": np.asarray(dj["P2"])}
        T0u, T0d = float(dj["T0up"]), float(dj["T0dn"])
        window_tol = 0.0
    else:
        d = load_uf_case(case_no)
        t = d.t
        cols = {"G": d.col("G_kgm2s"), "P1": d.col("P1_kPa"), "P2": d.col("P2_kPa")}
        T0u, T0d = float(d.station_mean("up")[0]), float(d.station_mean("dn")[0])
    if end_time is None:
        end_time = UF_END_OVERRIDE.get(case_no, float(t[-1]))
    tEnd = float(end_time)
    g = UF_GEOM
    A_f = math.pi / 4 * g["D"] ** 2
    kw = dict(knot=knot, half=half_window, tol=window_tol)
    mdot = boxcar_schedule(t, cols["G"], tEnd, shift=g_shift_s, scale=A_f, **kw)
    pin = boxcar_schedule(t, cols["P1"], tEnd, scale=1e3, **kw)
    pout = boxcar_schedule(t, cols["P2"], tEnd, scale=1e3 * exit_p_scale, **kw)
    xu, xd = g["stations"]["up"], g["stations"]["dn"]
    wall_T0 = lambda x: T0u + (T0d - T0u) * (x - xu) / (xd - xu)  # noqa: E731
    return build_pipe_chilldown(
        L=g["L"], D=g["D"], OD=g["OD"], wall_rho=g["rho_wall"], wall_material="stainless-steel-304",
        closure=closure, fluid="Nitrogen", wall_T0=wall_T0, p_in_schedule=pin, p_out_schedule=pout,
        mdot_schedule=mdot, N=N, dt=dt, end_time=tEnd, inlet_quality=inlet_quality, inertia=inertia,
        roughness=g["roughness"], stations=g["stations"], settings=settings, name=f"UF chilldown case {case_no}",
        **extra,
    )
