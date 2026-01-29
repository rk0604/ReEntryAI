"""
interval_math.py

Minimal interval arithmetic utilities for interval propagation (interval Euler, etc.).

Design goals:
- Works with floats OR intervals via promote()
- Includes inclusion functions for common math: sin, cos, sqrt, exp
- Includes vector ("box") helpers for re-entry state propagation

Important:
- Interval arithmetic is conservative (guaranteed enclosures), not tight.
- For division, if the denominator interval contains 0, raise ValueError by default.

Usage pattern:
    from interval_math import Interval, promote, box_add, box_scalar_mul

    # Example: interval Euler
    x_box_next = box_add(x_box, box_scalar_mul(dt, f_inclusion(x_box)))
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Sequence, Tuple, Union, Optional
import math


Number = Union[int, float]
MaybeInterval = Union["Interval", Number]
Box = List["Interval"]  # vector of intervals


@dataclass(frozen=True)
class Interval:
    """
    Closed interval [lo, hi].
    Invariant: lo <= hi, unless you construct via empty interval intentionally (not recommended).
    """
    lo: float
    hi: float

    def __post_init__(self) -> None:
        if self.lo > self.hi:
            raise ValueError(f"Invalid Interval: lo ({self.lo}) > hi ({self.hi})")

    # ---------------------------
    # Basic properties
    # ---------------------------
    def width(self) -> float:
        return self.hi - self.lo

    def mid(self) -> float:
        return 0.5 * (self.lo + self.hi)

    def contains(self, x: float) -> bool:
        return self.lo <= x <= self.hi

    def is_punctual(self) -> bool:
        return self.lo == self.hi

    # ---------------------------
    # Set hull / union wrapper
    # ---------------------------
    def hull(self, other: "Interval") -> "Interval":
        """Smallest interval containing self and other."""
        return Interval(min(self.lo, other.lo), max(self.hi, other.hi))

    # ---------------------------
    # Arithmetic operators
    # ---------------------------
    def __neg__(self) -> "Interval":
        return Interval(-self.hi, -self.lo)

    def __add__(self, other: MaybeInterval) -> "Interval":
        o = promote(other)
        return Interval(self.lo + o.lo, self.hi + o.hi)

    def __radd__(self, other: MaybeInterval) -> "Interval":
        return self.__add__(other)

    def __sub__(self, other: MaybeInterval) -> "Interval":
        o = promote(other)
        # [a,b] - [c,d] = [a-d, b-c]
        return Interval(self.lo - o.hi, self.hi - o.lo)

    def __rsub__(self, other: MaybeInterval) -> "Interval":
        # other - self
        o = promote(other)
        return Interval(o.lo - self.hi, o.hi - self.lo)

    def __mul__(self, other: MaybeInterval) -> "Interval":
        o = promote(other)
        a, b, c, d = self.lo, self.hi, o.lo, o.hi
        candidates = (a * c, a * d, b * c, b * d)
        return Interval(min(candidates), max(candidates))

    def __rmul__(self, other: MaybeInterval) -> "Interval":
        return self.__mul__(other)

    def reciprocal(self) -> "Interval":
        """1 / [a,b]. Raises if 0 in [a,b]."""
        if self.lo <= 0.0 <= self.hi:
            raise ValueError(f"Cannot take reciprocal: interval contains 0: {self}")
        # For positive-only or negative-only intervals:
        return Interval(1.0 / self.hi, 1.0 / self.lo)

    def __truediv__(self, other: MaybeInterval) -> "Interval":
        o = promote(other)
        return self * o.reciprocal()

    def __rtruediv__(self, other: MaybeInterval) -> "Interval":
        # other / self
        o = promote(other)
        return o * self.reciprocal()

    # ---------------------------
    # Useful extras
    # ---------------------------
    def abs(self) -> "Interval":
        """Interval absolute value enclosure."""
        if self.lo >= 0:
            return Interval(self.lo, self.hi)
        if self.hi <= 0:
            return Interval(-self.hi, -self.lo)
        # crosses 0
        return Interval(0.0, max(-self.lo, self.hi))
    
    def log(self) -> "Interval":
        """Natural log of a positive interval"""
        if self.lo <= 0.0:
            raise ValueError(f"log undefined for interval containing non-positive values: {self}")
        return Interval(math.log(self.lo), math.log(self.hi))

    def sqrt(self) -> "Interval":
        """sqrt([a,b]) with a >= 0 required."""
        if self.lo < 0:
            raise ValueError(f"sqrt undefined for interval with negative part: {self}")
        return Interval(math.sqrt(self.lo), math.sqrt(self.hi))

    def exp(self) -> "Interval":
        """exp is monotone increasing."""
        return Interval(math.exp(self.lo), math.exp(self.hi))

    def pow_int(self, n: int) -> "Interval":
        """Integer power enclosure (supports n >= 0)."""
        if n < 0:
            raise ValueError("pow_int expects n >= 0")
        if n == 0:
            return Interval(1.0, 1.0)
        if n == 1:
            return self
        if n % 2 == 1:
            # odd power is monotone increasing
            return Interval(self.lo ** n, self.hi ** n)
        # even power: could have min at 0
        a, b = self.lo, self.hi
        candidates = (a ** n, b ** n)
        if a <= 0.0 <= b:
            return Interval(0.0, max(candidates))
        return Interval(min(candidates), max(candidates))

    def sin(self) -> "Interval":
        return interval_sin(self)

    def cos(self) -> "Interval":
        return interval_cos(self)

    def __repr__(self) -> str:
        return f"Interval({self.lo}, {self.hi})"


def promote(x: MaybeInterval) -> Interval:
    """Convert float/int to punctual interval, leave Interval unchanged."""
    if isinstance(x, Interval):
        return x
    return Interval(float(x), float(x))


# ---------------------------
# Scalar helpers (explicit)
# ---------------------------
def scalar_times_interval(alpha: float, iv: Interval) -> Interval:
    """alpha * [lo,hi] with bound flip when alpha < 0."""
    if alpha >= 0:
        return Interval(alpha * iv.lo, alpha * iv.hi)
    return Interval(alpha * iv.hi, alpha * iv.lo)


# ---------------------------
# Trig inclusion functions
# ---------------------------
_TWO_PI = 2.0 * math.pi
_HALF_PI = 0.5 * math.pi


def _contains_critical(a: float, b: float, c: float, period: float) -> bool:
    """
    Returns True if there exists integer k such that c + k*period is in [a,b].
    """
    # Find k range where c + k*period in [a,b]
    k_min = math.ceil((a - c) / period)
    k_max = math.floor((b - c) / period)
    return k_min <= k_max


def interval_sin(x: Interval) -> Interval:
    """
    Inclusion for sin([a,b]).
    Strategy:
    - If width >= 2π -> [-1,1]
    - Else evaluate endpoints and include ±1 if critical points fall inside.
      sin reaches:
        +1 at  π/2 + 2kπ
        -1 at -π/2 + 2kπ  (equiv 3π/2 + 2kπ)
    """
    a, b = x.lo, x.hi
    if b - a >= _TWO_PI:
        return Interval(-1.0, 1.0)

    sa = math.sin(a)
    sb = math.sin(b)
    lo = min(sa, sb)
    hi = max(sa, sb)

    # Check if +1 is attained within [a,b]
    if _contains_critical(a, b, _HALF_PI, _TWO_PI):
        hi = 1.0
    # Check if -1 is attained within [a,b]
    if _contains_critical(a, b, -_HALF_PI, _TWO_PI):
        lo = -1.0

    return Interval(lo, hi)


def interval_cos(x: Interval) -> Interval:
    """
    Inclusion for cos([a,b]).
    cos reaches:
      +1 at 0 + 2kπ
      -1 at π + 2kπ
    """
    a, b = x.lo, x.hi
    if b - a >= _TWO_PI:
        return Interval(-1.0, 1.0)

    ca = math.cos(a)
    cb = math.cos(b)
    lo = min(ca, cb)
    hi = max(ca, cb)

    # +1 critical points
    if _contains_critical(a, b, 0.0, _TWO_PI):
        hi = 1.0
    # -1 critical points
    if _contains_critical(a, b, math.pi, _TWO_PI):
        lo = -1.0

    return Interval(lo, hi)


# ---------------------------
# Box (vector of intervals) helpers
# ---------------------------
def box_from_numbers(vals: Sequence[Number]) -> Box:
    return [Interval(float(v), float(v)) for v in vals]


def box_width(box: Box) -> float:
    """Max component width (∞-norm width)."""
    return max(iv.width() for iv in box) if box else 0.0


def box_mid(box: Box) -> List[float]:
    return [iv.mid() for iv in box]


def box_hull(a: Box, b: Box) -> Box:
    if len(a) != len(b):
        raise ValueError("box_hull: mismatched dimensions")
    return [a[i].hull(b[i]) for i in range(len(a))]


def box_add(a: Box, b: Box) -> Box:
    if len(a) != len(b):
        raise ValueError("box_add: mismatched dimensions")
    return [a[i] + b[i] for i in range(len(a))]


def box_sub(a: Box, b: Box) -> Box:
    if len(a) != len(b):
        raise ValueError("box_sub: mismatched dimensions")
    return [a[i] - b[i] for i in range(len(a))]


def box_scalar_mul(alpha: float, box: Box) -> Box:
    return [scalar_times_interval(alpha, iv) for iv in box]


def box_apply_unary(fn, box: Box) -> Box:
    """Apply a unary Interval->Interval function to each component."""
    return [fn(iv) for iv in box]


def box_contains(box: Box, point: Sequence[float]) -> bool:
    if len(box) != len(point):
        return False
    return all(box[i].contains(float(point[i])) for i in range(len(box)))


def box_split(box: Box, idx: Optional[int] = None) -> Tuple[Box, Box]:
    """
    Split a box into two boxes by bisecting one interval (default: widest dimension).
    This is your main tool later to reduce overestimation.
    """
    if not box:
        raise ValueError("box_split: empty box")

    if idx is None:
        # choose widest interval
        widths = [iv.width() for iv in box]
        idx = max(range(len(widths)), key=lambda i: widths[i])

    iv = box[idx]
    m = iv.mid()
    left_iv = Interval(iv.lo, m)
    right_iv = Interval(m, iv.hi)

    left = box[:]
    right = box[:]
    left[idx] = left_iv
    right[idx] = right_iv
    return left, right


# ---------------------------
# Convenience: "safe" helpers used a lot in dynamics
# ---------------------------
def dynamic_pressure(rho: MaybeInterval, V: MaybeInterval) -> Interval:
    """
    q = 0.5 * rho * V^2
    Works for rho and V being floats or Intervals.
    """
    rho_i = promote(rho)
    V_i = promote(V)
    return Interval(0.5, 0.5) * rho_i * V_i.pow_int(2)


def clamp_interval(iv: Interval, lo: float, hi: float) -> Interval:
    """
    Clamp interval to [lo, hi] (intersection with [lo,hi]).
    Raises if intersection is empty.
    """
    new_lo = max(iv.lo, lo)
    new_hi = min(iv.hi, hi)
    if new_lo > new_hi:
        raise ValueError(f"clamp_interval empty intersection: {iv} ∩ [{lo},{hi}]")
    return Interval(new_lo, new_hi)


# ---------------------------
# Small self-test (optional)
# ---------------------------
if __name__ == "__main__":
    # Quick sanity checks
    x = Interval(0.9, 1.1)
    fx = -x
    print("x:", x, "  -x:", fx)

    # One interval Euler step for x_dot = -x with h=0.1
    h = 0.1
    x1 = x + scalar_times_interval(h, -x)
    print("Euler step:", x1)

    # sin/cos checks
    print("sin([0, pi]):", interval_sin(Interval(0.0, math.pi)))
    print("cos([0, pi]):", interval_cos(Interval(0.0, math.pi)))

    # dynamic pressure check
    q = dynamic_pressure(Interval(1.0, 1.2), Interval(7000.0, 7100.0))
    print("q:", q)
