"""Closure plug-ins for convection conductors.

Two kinds:

* ``NamedClosure('miropolskii' | 'dittusBoelter' | 'darrHartwig' | 'ttWf', **extra)``
  -> correlation {model, diameter, flowArea, ...}.  darrHartwig/ttWf get the
  per-conductor ``axialPosition`` (and ttWf ``segmentLength``) from the builder.

* ``CustomClosure(fn, params, fluid)`` -> correlation {model:'custom',
  expression, params, diameter, flowArea}.  ``fn(s)`` is ordinary Python
  written with ``ofchill.expr`` functions over the scope object ``s``:

      s.t s.Tf s.Tw s.P s.G s.D s.area s.flowArea s.rho s.mu s.k s.cp s.Pr s.Re s.quality
      s.param('name')      -> params.name in the emitted expression
      s.fp                 -> FluidProps fits (s.fp.Tsat(s.P), s.fp.vap('k', T, s.P) ...)

  The same ``fn`` runs (a) symbolically -> ``closure.expression`` (string for
  OpenFLUME), (b) numerically -> ``closure.evaluate(Tw=..., P=..., ...)``.

  IMPORTANT (v0.3.1 scope semantics, see NOTES.md): an identifier absent from
  the runtime scope makes the WHOLE expression fail -> h = fallback floor
  (5 W/m^2K).  ``cp``/``Pr`` are absent inside the two-phase dome and
  ``quality`` is absent outside it; only reference them inside a ``where``
  branch that is guaranteed not to be taken otherwise (the ternary
  short-circuits).  ``s.x_eq`` / ``s.in_dome`` give a robust phase test from
  rho, P and the saturation-density fits that never touches absent names.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from typing import Any, Callable

import numpy as np

from . import expr as E
from .props import FluidProps, fluid_props

SCOPE_NAMES = ["t", "Tf", "Tw", "P", "G", "D", "area", "flowArea", "rho", "mu", "k", "cp", "Pr", "Re", "quality"]
NAMED_MODELS = ("miropolskii", "dittusBoelter", "darrHartwig", "ttWf")
FALLBACK_H_FLOOR = 5.0  # src/core/correlations/types.ts
DOME_EPS = 2e-4  # relative margin on the saturation densities for the dome test


class Closure:
    name: str = "closure"

    def correlation(self, D: float, flow_area: float, *, x: float | None = None, seg_len: float | None = None) -> dict:
        raise NotImplementedError

    def to_spec(self) -> dict:
        raise NotImplementedError


@dataclass
class NamedClosure(Closure):
    model: str
    extra: dict = field(default_factory=dict)
    label: str | None = None

    def __post_init__(self):
        if self.model not in NAMED_MODELS:
            raise ValueError(f"unknown named model {self.model!r}; one of {NAMED_MODELS}")

    @property
    def name(self) -> str:
        return self.label or self.model

    def correlation(self, D, flow_area, *, x=None, seg_len=None):
        c = {"model": self.model, "diameter": D, "flowArea": flow_area, **self.extra}
        if self.model in ("darrHartwig", "ttWf"):
            c["axialPosition"] = float(x if x is not None else 0.0)
        if self.model == "ttWf":
            if seg_len is None:
                raise ValueError("ttWf needs segmentLength")
            c["segmentLength"] = float(seg_len)
        return c

    def to_spec(self):
        return {"model": self.model, "extra": self.extra, "label": self.name}


class _Scope:
    """Scope object handed to closure functions."""

    def __init__(self, values: dict[str, Any], params: dict[str, float], fp: FluidProps | None, symbolic: bool):
        self._v = values
        self._params = params
        self._symbolic = symbolic
        self.fp = fp

    def __getattr__(self, name):
        if name.startswith("_"):
            raise AttributeError(name)
        if name in SCOPE_NAMES:
            if name not in self._v:
                raise KeyError(f"scope variable {name!r} not provided")
            return self._v[name]
        raise AttributeError(name)

    def param(self, name: str):
        if name not in self._params:
            raise KeyError(f"unknown param {name!r}")
        return E.Var(f"params.{name}") if self._symbolic else self._params[name]

    # robust phase helpers (never reference cp / quality) --------------------
    def x_eq(self):
        """Equilibrium quality from the node density and the saturation-density
        fits (HEM: 1/rho = x/rho_v + (1-x)/rho_l), clipped to [0, 1]."""
        P, rho = self.P, self.rho
        vl, vv = 1.0 / self.fp.rho_l(P), 1.0 / self.fp.rho_v(P)
        return E.clip((1.0 / rho - vl) / (vv - vl), 0.0, 1.0)

    def in_dome(self, eps: float = DOME_EPS):
        """True when rho lies within the saturation densities widened by eps
        (so a real dome node is never classified single-phase)."""
        P, rho = self.P, self.rho
        return E.logical_and(rho < self.fp.rho_l(P) * (1 + eps), rho > self.fp.rho_v(P) * (1 - eps))


@dataclass
class CustomClosure(Closure):
    fn: Callable[[_Scope], Any]
    params: dict = field(default_factory=dict)
    fluid: str | None = "Nitrogen"
    label: str | None = None

    @property
    def name(self) -> str:
        return self.label or getattr(self.fn, "__name__", "custom")

    @property
    def fp(self) -> FluidProps | None:
        return fluid_props(self.fluid) if self.fluid else None

    @property
    def expression(self) -> str:
        if not hasattr(self, "_expr_cache"):
            scope = _Scope({n: E.Var(n) for n in SCOPE_NAMES}, self.params, self.fp, True)
            out = self.fn(scope)
            self._expr_cache = E._emit(out)
        return self._expr_cache

    def evaluate(self, _params: dict | None = None, **values) -> np.ndarray:
        """Numerical evaluation of the RAW h (no floor/relaxation).  Pass the
        scope quantities as scalars or arrays (only those the closure uses);
        ``_params`` overrides params (e.g. the per-conductor ``z``)."""
        vals = {k: np.asarray(v, dtype=float) for k, v in values.items()}
        scope = _Scope(vals, {**self.params, **(_params or {})}, self.fp, False)
        with np.errstate(all="ignore"):
            return np.asarray(self.fn(scope), dtype=float)

    def correlation(self, D, flow_area, *, x=None, seg_len=None):
        params = {k: float(v) for k, v in self.params.items()}
        if "z" in params and x is not None:  # per-conductor axial position (reserved param name)
            params["z"] = float(x)
        return {
            "model": "custom",
            "expression": self.expression,
            "params": params,
            "diameter": D,
            "flowArea": flow_area,
        }

    def to_spec(self):
        h = hashlib.sha1(self.expression.encode()).hexdigest()[:10]
        return {"model": "custom", "label": self.name, "params": self.params, "fluid": self.fluid, "expr_sha1": h}

    def with_params(self, **kw) -> "CustomClosure":
        return CustomClosure(self.fn, {**self.params, **kw}, self.fluid, self.label)

    def with_label(self, label: str) -> "CustomClosure":
        return CustomClosure(self.fn, dict(self.params), self.fluid, label)


def closure(params: dict | None = None, fluid: str | None = "Nitrogen", label: str | None = None):
    """Decorator: ``@closure(params={'C': 0.023})`` turns ``fn(s)`` into a CustomClosure."""

    def deco(fn):
        return CustomClosure(fn, dict(params or {}), fluid, label or fn.__name__)

    return deco
