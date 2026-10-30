"""Columnar evaluation of ROOT/C++ style cut and weight strings.

The ROOT macros steer everything with ``TTreeFormula`` strings such as::

    (l1_eta<1.37 || l1_eta>1.52) && fabs(el_eta[0])<2.47 && l1_IFF==3

This module parses that C++ expression subset into a small AST and evaluates it
on whole arrays (numpy, or awkward for jagged branches), so the existing
selection strings can be reused verbatim in the YAML configs.

Supported grammar (C++ precedence):
    ternary ?:  ||  &&  |  &  == !=  < <= > >=  + -  * / %  unary ! - +
    postfix: ``x[i]`` (constant index into a jagged branch), ``f(a, b)``
Functions: fabs/abs, sqrt, log, log10, exp, pow, sin, cos, tan, sinh, cosh,
tanh, atan2, min, max, floor, ceil and their ``TMath::`` / ``std::`` spellings.

Semantics follow ``TTreeFormula``: logical operators return 0/1, any non-zero
value is true, and an event where a ``x[i]`` index is out of range evaluates to
0 for the whole expression (ROOT skips such entries).
"""

from __future__ import annotations

import re
from functools import lru_cache

import numpy as np

try:  # awkward is only needed for jagged (vector) branches
    import awkward as ak
except ImportError:  # pragma: no cover
    ak = None


class ExprError(ValueError):
    """Raised for malformed expressions or missing variables."""


_TOKEN_RE = re.compile(
    r"""
    (?P<ws>\s+)
  | (?P<num>(?:\d+\.\d*|\.\d+|\d+)(?:[eE][+-]?\d+)?[fFlLuU]*)
  | (?P<name>[A-Za-z_][A-Za-z0-9_]*(?:::[A-Za-z_][A-Za-z0-9_]*)*(?:\$)?)
  | (?P<op>&&|\|\||==|!=|<=|>=|[-+*/%<>!()\[\],?:&|])
    """,
    re.VERBOSE,
)

_FUNCS = {
    "fabs": np.abs, "abs": np.abs,
    "sqrt": np.sqrt, "log": np.log, "log10": np.log10, "exp": np.exp,
    "pow": np.power, "power": np.power,
    "sin": np.sin, "cos": np.cos, "tan": np.tan,
    "sinh": np.sinh, "cosh": np.cosh, "tanh": np.tanh,
    "atan2": np.arctan2, "min": np.minimum, "max": np.maximum,
    "floor": np.floor, "ceil": np.ceil,
}


def _canonical_func(name):
    for prefix in ("TMath::", "std::"):
        if name.startswith(prefix):
            name = name[len(prefix):]
    return name.lower()


def _tokenize(text):
    pos, out = 0, []
    while pos < len(text):
        m = _TOKEN_RE.match(text, pos)
        if not m:
            raise ExprError(f"Unexpected character {text[pos]!r} at {pos} in: {text}")
        pos = m.end()
        kind = m.lastgroup
        if kind == "ws":
            continue
        val = m.group(kind)
        if kind == "num":
            val = val.rstrip("fFlLuU")
            val = float(val) if any(c in val for c in ".eE") else int(val)
        out.append((kind, val))
    out.append(("end", None))
    return out


class _Parser:
    # binary operators, lowest precedence first
    _LEVELS = [("||",), ("&&",), ("|",), ("&",), ("==", "!="), ("<", "<=", ">", ">="), ("+", "-"), ("*", "/", "%")]

    def __init__(self, text):
        self.text = text
        self.toks = _tokenize(text)
        self.i = 0

    def peek(self):
        return self.toks[self.i]

    def take(self, value=None):
        tok = self.toks[self.i]
        if value is not None and tok[1] != value:
            raise ExprError(f"Expected {value!r} but found {tok[1]!r} in: {self.text}")
        self.i += 1
        return tok

    def parse(self):
        node = self.ternary()
        if self.peek()[0] != "end":
            raise ExprError(f"Trailing tokens from {self.peek()[1]!r} in: {self.text}")
        return node

    def ternary(self):
        cond = self.binary(0)
        if self.peek() == ("op", "?"):
            self.take("?")
            a = self.ternary()
            self.take(":")
            b = self.ternary()
            return ("?", cond, a, b)
        return cond

    def binary(self, level):
        if level == len(self._LEVELS):
            return self.unary()
        node = self.binary(level + 1)
        while self.peek()[0] == "op" and self.peek()[1] in self._LEVELS[level]:
            op = self.take()[1]
            node = ("bin", op, node, self.binary(level + 1))
        return node

    def unary(self):
        tok = self.peek()
        if tok[0] == "op" and tok[1] in ("!", "-", "+"):
            self.take()
            return ("un", tok[1], self.unary())
        return self.postfix()

    def postfix(self):
        node = self.atom()
        while self.peek() == ("op", "["):
            self.take("[")
            idx = self.ternary()
            self.take("]")
            if idx[0] != "num" or not isinstance(idx[1], int):
                raise ExprError(f"Only constant integer indices are supported in: {self.text}")
            node = ("idx", node, idx[1])
        return node

    def atom(self):
        kind, val = self.take()
        if kind == "num":
            return ("num", val)
        if kind == "name":
            if val in ("true", "kTRUE"):
                return ("num", 1)
            if val in ("false", "kFALSE"):
                return ("num", 0)
            if self.peek() == ("op", "("):
                fname = _canonical_func(val)
                if fname not in _FUNCS:
                    raise ExprError(f"Unknown function {val!r} in: {self.text}")
                self.take("(")
                args = []
                if self.peek() != ("op", ")"):
                    args.append(self.ternary())
                    while self.peek() == ("op", ","):
                        self.take(",")
                        args.append(self.ternary())
                self.take(")")
                return ("call", fname, args)
            return ("var", val)
        if (kind, val) == ("op", "("):
            node = self.ternary()
            self.take(")")
            return node
        raise ExprError(f"Unexpected token {val!r} in: {self.text}")


@lru_cache(maxsize=None)
def parse(text):
    """Parse an expression string into an AST (cached)."""
    if text is None or str(text).strip() == "":
        return ("num", 1)
    return _Parser(str(text)).parse()


def variables(text):
    """Return the set of branch names an expression depends on."""
    out = set()

    def walk(node):
        tag = node[0]
        if tag == "var":
            out.add(node[1])
        elif tag == "idx":
            walk(node[1])
        elif tag == "un":
            walk(node[2])
        elif tag == "bin":
            walk(node[2]); walk(node[3])
        elif tag == "?":
            walk(node[1]); walk(node[2]); walk(node[3])
        elif tag == "call":
            for a in node[2]:
                walk(a)

    walk(parse(text))
    return out


def _to_numpy(x):
    """numpy view of a branch; floats are promoted to double, as TTreeFormula
    evaluates everything in double precision."""
    if ak is not None and isinstance(x, ak.Array):
        x = ak.to_numpy(x)
    x = np.asarray(x)
    if x.dtype.kind == "f" and x.dtype != np.float64:
        x = x.astype(np.float64)
    return x


def _truth(x):
    x = np.asarray(x)
    return x if x.dtype == np.bool_ else x != 0


class _Evaluator:
    def __init__(self, arrays, n):
        self.arrays = arrays
        self.n = n
        self.valid = None  # entries where every x[i] index was in range

    def _flag_invalid(self, ok):
        self.valid = ok if self.valid is None else (self.valid & ok)

    def ev(self, node):
        tag = node[0]
        if tag == "num":
            return node[1]
        if tag == "var":
            name = node[1]
            if name not in self.arrays:
                raise ExprError(f"Variable {name!r} not available (missing branch?)")
            arr = self.arrays[name]
            if ak is not None and isinstance(arr, ak.Array) and arr.ndim > 1:
                return arr  # jagged: must be indexed
            return _to_numpy(arr)
        if tag == "idx":
            arr = self.arrays.get(node[1][1]) if node[1][0] == "var" else None
            if arr is None:
                arr = self.ev(node[1])
            k = node[2]
            if ak is not None and isinstance(arr, ak.Array) and arr.ndim > 1:
                ok = _to_numpy(ak.num(arr, axis=1) > k)
                self._flag_invalid(ok)
                padded = ak.pad_none(arr, k + 1, axis=1)[:, k]
                vals = _to_numpy(ak.fill_none(padded, 0))
                return vals
            arr = np.asarray(arr)
            if arr.ndim == 2:
                if arr.shape[1] <= k:
                    self._flag_invalid(np.zeros(len(arr), bool))
                    return np.zeros(len(arr))
                return arr[:, k]
            raise ExprError(f"Cannot index scalar branch {node[1]} with [{k}]")
        if tag == "un":
            v = self.ev(node[2])
            if node[1] == "!":
                return ~_truth(v) if np.ndim(v) else int(not v)
            if node[1] == "-":
                return -_np(v)
            return v
        if tag == "bin":
            op = node[1]
            a = self.ev(node[2])
            b = self.ev(node[3])
            if op == "&&":
                return _truth(a) & _truth(b)
            if op == "||":
                return _truth(a) | _truth(b)
            if op in ("&", "|"):
                a, b = np.asarray(a), np.asarray(b)
                if a.dtype.kind in "iu" and b.dtype.kind in "iu":
                    return a & b if op == "&" else a | b
                return (_truth(a) & _truth(b)) if op == "&" else (_truth(a) | _truth(b))
            a, b = np.asarray(a), np.asarray(b)
            if op == "==":
                return a == b
            if op == "!=":
                return a != b
            if op == "<":
                return a < b
            if op == "<=":
                return a <= b
            if op == ">":
                return a > b
            if op == ">=":
                return a >= b
            a, b = _np(a), _np(b)
            if op == "+":
                return a + b
            if op == "-":
                return a - b
            if op == "*":
                return a * b
            if op in ("/", "%"):
                ints = np.asarray(a).dtype.kind in "iub" and np.asarray(b).dtype.kind in "iub"
                with np.errstate(divide="ignore", invalid="ignore"):
                    if op == "%":
                        return np.fmod(a, b)
                    q = np.true_divide(a, b)
                    return np.trunc(q) if ints else q
        if tag == "?":
            c = _truth(self.ev(node[1]))
            return np.where(c, self.ev(node[2]), self.ev(node[3]))
        if tag == "call":
            args = [_np(self.ev(a)) for a in node[2]]
            with np.errstate(divide="ignore", invalid="ignore"):
                return _FUNCS[node[1]](*args)
        raise ExprError(f"Bad node {node}")  # pragma: no cover


def _np(v):
    """Operand for arithmetic: booleans count as 0/1 integers (C++)."""
    if np.ndim(v) == 0:
        return v
    a = np.asarray(v)
    return a.astype(np.int64) if a.dtype == np.bool_ else a


def evaluate(text, arrays, n=None, dtype=None):
    """Evaluate ``text`` on a mapping of branch arrays.

    Returns a numpy array of length ``n`` (number of entries in the chunk), of
    the natural result type (bool for cuts) unless ``dtype`` is given.
    Entries for which an ``x[i]`` index was out of range are set to 0.
    """
    if n is None:
        n = len(next(iter(arrays.values())))
    ev = _Evaluator(arrays, n)
    out = ev.ev(parse(text))
    if np.ndim(out) == 0:
        out = np.full(n, out, dtype=dtype or (np.float64 if isinstance(out, float) else np.int64))
    else:
        out = np.asarray(out)
        if out.dtype == np.bool_ and dtype is None:
            pass
        elif dtype is not None and out.dtype != dtype:
            out = out.astype(dtype)
    if ev.valid is not None:
        out = np.where(ev.valid, out, np.zeros((), out.dtype))
    return out


def evaluate_mask(text, arrays, n=None):
    """Evaluate a cut string as a boolean mask."""
    return _truth(evaluate(text, arrays, n))
