"""Closure compiler + property-fit tests.

Run:  .venv/bin/python -m ofchill.tests.test_closure_expr   (or pytest, if installed)
"""

from __future__ import annotations

import json

import numpy as np

from ofchill import expr as E
from ofchill.closure import closure
from ofchill.closures import db_vapor, dittus_boelter, miropolskii_expr, probe_all
from ofchill.props import verify_fits
from ofchill.verify import closure_vs_numpy

INLETS = ({"temperature": 150.0}, {"quality": 0.3}, {"quality": 0.0})


def test_emit_basic():
    x = E.Var("x")
    assert E._emit(-x**2) == "(-(x ^ 2.0))"
    assert E._emit(E.where(x < 1, 2 * x, E.minimum(x, 3))) == "((x < 1.0) ? (2.0 * x) : min(x, 3.0))"
    assert E._emit(E.clip(x, -1.5, 2)) == "clamp(x, (-1.5), 2.0)"


def test_numpy_mode_matches_python():
    x = np.linspace(0.1, 3, 7)
    np.testing.assert_allclose(E.where(x < 1, 2 * x, E.minimum(x, 3)), np.where(x < 1, 2 * x, np.minimum(x, 3)))
    np.testing.assert_allclose(dittus_boelter(100.0, 0.01, 1e-5, 0.02, 1000.0), 0.023 * 1e5**0.8 * 0.5**0.4 * 2)


def test_expression_vs_numpy_tiny_network():
    """OpenFLUME-reported h (fresh evaluation at each recorded state, floored
    at 5) == numpy evaluation of the same closure on the reported state."""
    for cl in (probe_all, db_vapor, miropolskii_expr):
        for inlet in INLETS:
            r = closure_vs_numpy(cl, inlet)
            assert r["converged"], (cl.name, inlet)
            assert r["max_rel_err"] < 1e-9, (cl.name, inlet, r["max_rel_err"])


def test_absent_identifier_falls_back_to_floor():
    """`quality` is absent outside the dome -> whole expression -> 5 W/m2K."""

    @closure(params={}, fluid=None, label="uses_quality")
    def uses_quality(s):
        return 1000.0 + 0 * s.quality

    import openflume as of

    from ofchill.verify import tiny_network

    raw = of.Network.from_dict(tiny_network(uses_quality, inlet={"temperature": 150.0})).solve().raw
    h = raw["conductors"]["c"]["heatTransferCoeff"]
    assert all(abs(v - 5.0) < 1e-12 for v in h), h
    raw = of.Network.from_dict(tiny_network(uses_quality, inlet={"quality": 0.3})).solve().raw
    ph = raw["nodes"]["f"]["phase"]
    h = raw["conductors"]["c"]["heatTransferCoeff"]
    assert all((v == 1000.0) == (p == "twoPhase") for v, p in zip(h, ph)), (h, ph)


def test_property_fits():
    for fl in ("Nitrogen", "Oxygen", "Hydrogen"):
        r = verify_fits(fl, n_sat=300, n_vap=1500)
        assert r["max_sat"] < 0.01 and r["max_vap"] < 0.01, json.dumps(r)


if __name__ == "__main__":
    import sys
    import time

    fails = 0
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            t0 = time.time()
            try:
                fn()
                print(f"PASS {name} ({time.time() - t0:.1f}s)")
            except Exception as e:  # noqa: BLE001
                fails += 1
                print(f"FAIL {name}: {type(e).__name__}: {e}")
    sys.exit(1 if fails else 0)
