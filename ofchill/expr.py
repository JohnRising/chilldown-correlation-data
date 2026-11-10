"""Tiny dual-mode expression DSL: the same Python code evaluates numerically
(numpy) or builds an expression tree that is emitted in OpenFLUME's safe
expression language (src/core/usercode/expression.ts, v0.3.1).

Use the functions in this module (``exp``, ``log``, ``where``, ``minimum``...)
instead of numpy's inside closures; they dispatch on argument type.

Emission rules (all verified against the v0.3.1 parser):
  * every binary operation is fully parenthesised (no reliance on precedence);
  * ``**`` -> ``^`` (Math.pow), ``where`` -> ternary ``c ? a : b`` (short-circuit),
    ``minimum/maximum`` -> ``min/max``, ``clip`` -> ``clamp``;
  * ``&`` / ``|`` on boolean Exprs -> ``&&`` / ``||``;
  * numbers are emitted with ``repr`` (round-trip exact doubles).
The language has no ``let``: shared subtrees are duplicated in the string.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

Number = (int, float, np.floating, np.integer)


class Expr:
    __array_priority__ = 1000  # make numpy defer to our reflected operators

    # --- arithmetic -------------------------------------------------------
    def __add__(self, o):
        return Bin("+", self, o)

    def __radd__(self, o):
        return Bin("+", o, self)

    def __sub__(self, o):
        return Bin("-", self, o)

    def __rsub__(self, o):
        return Bin("-", o, self)

    def __mul__(self, o):
        return Bin("*", self, o)

    def __rmul__(self, o):
        return Bin("*", o, self)

    def __truediv__(self, o):
        return Bin("/", self, o)

    def __rtruediv__(self, o):
        return Bin("/", o, self)

    def __pow__(self, o):
        return Bin("^", self, o)

    def __rpow__(self, o):
        return Bin("^", o, self)

    def __neg__(self):
        return Un("-", self)

    def __pos__(self):
        return self

    # --- comparisons / logic (return boolean Exprs) -------------------------
    def __lt__(self, o):
        return Bin("<", self, o)

    def __le__(self, o):
        return Bin("<=", self, o)

    def __gt__(self, o):
        return Bin(">", self, o)

    def __ge__(self, o):
        return Bin(">=", self, o)

    def __and__(self, o):
        return Bin("&&", self, o)

    def __rand__(self, o):
        return Bin("&&", o, self)

    def __or__(self, o):
        return Bin("||", self, o)

    def __ror__(self, o):
        return Bin("||", o, self)

    def __invert__(self):
        return Un("!", self)

    def __bool__(self):
        raise TypeError("symbolic Expr has no truth value; use where(cond, a, b)")

    # --- emission -------------------------------------------------------
    def emit(self) -> str:  # pragma: no cover - abstract
        raise NotImplementedError

    def __str__(self) -> str:
        return self.emit()

    def __repr__(self) -> str:
        s = self.emit()
        return f"Expr({s[:120]}{'...' if len(s) > 120 else ''})"


def _emit(x: Any) -> str:
    if isinstance(x, Expr):
        return x.emit()
    if isinstance(x, (bool, np.bool_)):
        raise TypeError("boolean literals are not supported in emitted expressions")
    if isinstance(x, Number):
        v = float(x)
        if not math.isfinite(v):
            raise ValueError(f"non-finite literal {v}")
        s = repr(v)
        return f"({s})" if v < 0 else s
    raise TypeError(f"cannot emit {type(x).__name__}")


class Var(Expr):
    def __init__(self, name: str):
        self.name = name

    def emit(self) -> str:
        return self.name


class Bin(Expr):
    def __init__(self, op: str, a, b):
        self.op, self.a, self.b = op, a, b

    def emit(self) -> str:
        return f"({_emit(self.a)} {self.op} {_emit(self.b)})"


class Un(Expr):
    def __init__(self, op: str, a):
        self.op, self.a = op, a

    def emit(self) -> str:
        return f"({self.op}{_emit(self.a)})"


class Call(Expr):
    def __init__(self, fn: str, *args):
        self.fn, self.args = fn, args

    def emit(self) -> str:
        return f"{self.fn}({', '.join(_emit(a) for a in self.args)})"


class Cond(Expr):
    def __init__(self, c, a, b):
        self.c, self.a, self.b = c, a, b

    def emit(self) -> str:
        return f"({_emit(self.c)} ? {_emit(self.a)} : {_emit(self.b)})"


def is_sym(*xs) -> bool:
    return any(isinstance(x, Expr) for x in xs)


# --- dual-mode functions ----------------------------------------------------
def _f1(name: str, npf):
    def f(x):
        if isinstance(x, Expr):
            return Call(name, x)
        with np.errstate(all="ignore"):
            return npf(x)

    f.__name__ = name
    return f


exp = _f1("exp", np.exp)
log = _f1("log", np.log)
sqrt = _f1("sqrt", np.sqrt)
tanh = _f1("tanh", np.tanh)
fabs = _f1("abs", np.abs)
sin = _f1("sin", np.sin)
cos = _f1("cos", np.cos)


def minimum(*xs):
    if is_sym(*xs):
        return Call("min", *xs)
    out = xs[0]
    for x in xs[1:]:
        out = np.minimum(out, x)
    return out


def maximum(*xs):
    if is_sym(*xs):
        return Call("max", *xs)
    out = xs[0]
    for x in xs[1:]:
        out = np.maximum(out, x)
    return out


def clip(x, lo, hi):
    if is_sym(x, lo, hi):
        return Call("clamp", x, lo, hi)
    return np.minimum(np.maximum(x, lo), hi)


def smoothstep(e0, e1, x):
    if is_sym(e0, e1, x):
        return Call("smoothstep", e0, e1, x)
    x = np.asarray(x, float)
    e0, e1 = np.asarray(e0, float), np.asarray(e1, float)
    with np.errstate(all="ignore"):
        t = np.clip((x - e0) / (e1 - e0), 0.0, 1.0)
        out = t * t * (3 - 2 * t)
    return np.where(e1 == e0, np.where(x < e0, 0.0, 1.0), out)


def where(c, a, b):
    """Ternary.  Symbolic -> ``(c ? a : b)`` (short-circuit in OpenFLUME, so the
    untaken branch may reference identifiers absent from the scope)."""
    if is_sym(c, a, b):
        return Cond(c, a, b)
    with np.errstate(all="ignore"):
        return np.where(c, a, b)


def logical_and(a, b):
    if is_sym(a, b):
        return Bin("&&", a, b)
    return np.logical_and(a, b)


def logical_or(a, b):
    if is_sym(a, b):
        return Bin("||", a, b)
    return np.logical_or(a, b)


def poly(coefs, x):
    """Horner polynomial sum(c_i x^i) (works for Expr and arrays)."""
    acc = coefs[-1]
    for c in reversed(coefs[:-1]):
        acc = acc * x + c
    return acc


pi = math.pi
