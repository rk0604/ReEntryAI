"""
interval_math.py

Minimal interval arithmetic utilities for interval propagation, such as the
interval Euler step.

Design goals:
    Works with floats or intervals through promote().
    Includes inclusion functions for common math: sin, cos, sqrt, exp.
    Includes vector helpers, called boxes, for re-entry state propagation.

Notes:
    Interval arithmetic is conservative. It returns guaranteed enclosures, not
    tight ranges.
    For division, an interval denominator that contains zero raises ValueError.

Usage pattern:
    from interval_math import Interval, promote, box_add, box_scalar_mul

    # One interval Euler step:
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
    A closed interval [lo, hi].

    Invariant: lo is less than or equal to hi. Operators are defined so that the
    result always encloses the true result of the same operation on every pair
    of points drawn from the operands.

    Fields:
        lo: lower bound.
        hi: upper bound.
    """
    lo: float
    hi: float

    def __post_init__(self) -> None:
        """Reject an inverted interval where lo is greater than hi."""
        if self.lo > self.hi:
            raise ValueError(f"Invalid Interval: lo ({self.lo}) > hi ({self.hi})")

    # Basic properties
    def width(self) -> float:
        """Return the width hi minus lo."""
        return self.hi - self.lo

    def mid(self) -> float:
        """Return the midpoint, the average of lo and hi."""
        return 0.5 * (self.lo + self.hi)

    def contains(self, x: float) -> bool:
        """Return True when x lies within the closed interval."""
        return self.lo <= x <= self.hi

    def is_punctual(self) -> bool:
        """Return True when the interval is a single point, lo equal to hi."""
        return self.lo == self.hi

    # Set hull
    def hull(self, other: "Interval") -> "Interval":
        """Return the smallest interval that contains both this and other."""
        return Interval(min(self.lo, other.lo), max(self.hi, other.hi))

    # Arithmetic operators
    def __neg__(self) -> "Interval":
        """Return the negated interval [minus hi, minus lo]."""
        return Interval(-self.hi, -self.lo)

    def __add__(self, other: MaybeInterval) -> "Interval":
        """Return the sum interval. The other operand may be a number."""
        o = promote(other)
        return Interval(self.lo + o.lo, self.hi + o.hi)

    def __radd__(self, other: MaybeInterval) -> "Interval":
        """Support number plus interval."""
        return self.__add__(other)

    def __sub__(self, other: MaybeInterval) -> "Interval":
        """Return the difference interval [lo minus other.hi, hi minus other.lo]."""
        o = promote(other)
        # [a,b] - [c,d] = [a-d, b-c]
        return Interval(self.lo - o.hi, self.hi - o.lo)

    def __rsub__(self, other: MaybeInterval) -> "Interval":
        """Support number minus interval, computed as other minus self."""
        o = promote(other)
        return Interval(o.lo - self.hi, o.hi - self.lo)

    def __mul__(self, other: MaybeInterval) -> "Interval":
        """Return the product interval, the min and max of the four corner products."""
        o = promote(other)
        a, b, c, d = self.lo, self.hi, o.lo, o.hi
        candidates = (a * c, a * d, b * c, b * d)
        return Interval(min(candidates), max(candidates))

    def __rmul__(self, other: MaybeInterval) -> "Interval":
        """Support number times interval."""
        return self.__mul__(other)

    def reciprocal(self) -> "Interval":
        """
        Return 1 divided by the interval.

        Raises ValueError when the interval contains zero, since the reciprocal
        is unbounded there.
        """
        if self.lo <= 0.0 <= self.hi:
            raise ValueError(f"Cannot take reciprocal: interval contains 0: {self}")
        # Valid for a strictly positive or strictly negative interval.
        return Interval(1.0 / self.hi, 1.0 / self.lo)

    def __truediv__(self, other: MaybeInterval) -> "Interval":
        """Return self divided by other, as self times the reciprocal of other."""
        o = promote(other)
        return self * o.reciprocal()

    def __rtruediv__(self, other: MaybeInterval) -> "Interval":
        """Support number divided by interval, computed as other times reciprocal."""
        o = promote(other)
        return o * self.reciprocal()

    # Useful extras
    def abs(self) -> "Interval":
        """
        Return the interval enclosure of the absolute value.

        When the interval crosses zero the lower bound becomes zero and the
        upper bound is the larger magnitude of the two ends.
        """
        if self.lo >= 0:
            return Interval(self.lo, self.hi)
        if self.hi <= 0:
            return Interval(-self.hi, -self.lo)
        # The interval crosses zero.
        return Interval(0.0, max(-self.lo, self.hi))

    def log(self) -> "Interval":
        """Return the natural log of a strictly positive interval."""
        if self.lo <= 0.0:
            raise ValueError(f"log undefined for interval containing non-positive values: {self}")
        return Interval(math.log(self.lo), math.log(self.hi))

    def sqrt(self) -> "Interval":
        """Return the square root. Requires the lower bound to be non negative."""
        if self.lo < 0:
            raise ValueError(f"sqrt undefined for interval with negative part: {self}")
        return Interval(math.sqrt(self.lo), math.sqrt(self.hi))

    def exp(self) -> "Interval":
        """Return the exponential. The function is increasing, so ends map directly."""
        return Interval(math.exp(self.lo), math.exp(self.hi))

    def pow_int(self, n: int) -> "Interval":
        """
        Return the interval raised to a non negative integer power.

        Odd powers are increasing, so the ends map directly. Even powers can
        reach zero when the interval crosses zero, so the lower bound becomes
        zero in that case.

        Parameters:
            n: the integer exponent, which must be zero or greater.
        """
        if n < 0:
            raise ValueError("pow_int expects n >= 0")
        if n == 0:
            return Interval(1.0, 1.0)
        if n == 1:
            return self
        if n % 2 == 1:
            # Odd power is monotone increasing.
            return Interval(self.lo ** n, self.hi ** n)
        # Even power can attain its minimum at zero.
        a, b = self.lo, self.hi
        candidates = (a ** n, b ** n)
        if a <= 0.0 <= b:
            return Interval(0.0, max(candidates))
        return Interval(min(candidates), max(candidates))

    def sin(self) -> "Interval":
        """Return the interval enclosure of sin over this interval."""
        return interval_sin(self)

    def cos(self) -> "Interval":
        """Return the interval enclosure of cos over this interval."""
        return interval_cos(self)

    def __repr__(self) -> str:
        """Return a readable representation, for example Interval(0.9, 1.1)."""
        return f"Interval({self.lo}, {self.hi})"


def promote(x: MaybeInterval) -> Interval:
    """
    Convert a number into a point interval, or return an Interval unchanged.

    Parameters:
        x: an Interval or a number.

    Returns:
        x itself when it is already an Interval, otherwise Interval(x, x).
    """
    if isinstance(x, Interval):
        return x
    return Interval(float(x), float(x))


# Scalar helpers
def scalar_times_interval(alpha: float, iv: Interval) -> Interval:
    """
    Multiply an interval by a scalar, flipping the bounds when the scalar is
    negative.

    Parameters:
        alpha: the scalar multiplier.
        iv:    the interval to scale.

    Returns:
        The scaled interval.
    """
    if alpha >= 0:
        return Interval(alpha * iv.lo, alpha * iv.hi)
    return Interval(alpha * iv.hi, alpha * iv.lo)

# Trig inclusion functions
_TWO_PI = 2.0 * math.pi
_HALF_PI = 0.5 * math.pi


def _contains_critical(a: float, b: float, c: float, period: float) -> bool:
    """
    Report whether a critical point lies inside an interval.

    Returns True when some integer k exists such that c plus k times period
    falls inside [a, b].

    Parameters:
        a, b:   interval ends.
        c:      a base critical point of the function.
        period: the spacing between critical points.
    """
    # Find the range of k for which c + k*period lands in [a, b].
    k_min = math.ceil((a - c) / period)
    k_max = math.floor((b - c) / period)
    return k_min <= k_max


def interval_sin(x: Interval) -> Interval:
    """
    Return the interval enclosure of sin over the interval x.

    Strategy: if the interval spans a full period or more, the result is the
    whole range [-1, 1]. Otherwise evaluate the ends and widen to plus one or
    minus one when a peak or trough falls inside. Sine reaches plus one at
    pi/2 plus multiples of two pi, and minus one at minus pi/2 plus multiples of
    two pi.
    """
    a, b = x.lo, x.hi
    if b - a >= _TWO_PI:
        return Interval(-1.0, 1.0)

    sa = math.sin(a)
    sb = math.sin(b)
    lo = min(sa, sb)
    hi = max(sa, sb)

    # Widen to plus one when a peak is inside the interval.
    if _contains_critical(a, b, _HALF_PI, _TWO_PI):
        hi = 1.0
    # Widen to minus one when a trough is inside the interval.
    if _contains_critical(a, b, -_HALF_PI, _TWO_PI):
        lo = -1.0

    return Interval(lo, hi)


def interval_cos(x: Interval) -> Interval:
    """
    Return the interval enclosure of cos over the interval x.

    Cosine reaches plus one at multiples of two pi and minus one at pi plus
    multiples of two pi. The same widen at the ends approach as interval_sin is
    used.
    """
    a, b = x.lo, x.hi
    if b - a >= _TWO_PI:
        return Interval(-1.0, 1.0)

    ca = math.cos(a)
    cb = math.cos(b)
    lo = min(ca, cb)
    hi = max(ca, cb)

    # Widen to plus one at a peak.
    if _contains_critical(a, b, 0.0, _TWO_PI):
        hi = 1.0
    # Widen to minus one at a trough.
    if _contains_critical(a, b, math.pi, _TWO_PI):
        lo = -1.0

    return Interval(lo, hi)


# Box (vector of intervals) helpers
def box_from_numbers(vals: Sequence[Number]) -> Box:
    """Build a box of point intervals from a sequence of numbers."""
    return [Interval(float(v), float(v)) for v in vals]


def box_width(box: Box) -> float:
    """Return the largest component width across the box, or zero when empty."""
    return max(iv.width() for iv in box) if box else 0.0


def box_mid(box: Box) -> List[float]:
    """Return the midpoint of each component as a list of floats."""
    return [iv.mid() for iv in box]


def box_hull(a: Box, b: Box) -> Box:
    """Return the component wise hull of two boxes of equal length."""
    if len(a) != len(b):
        raise ValueError("box_hull: mismatched dimensions")
    return [a[i].hull(b[i]) for i in range(len(a))]


def box_add(a: Box, b: Box) -> Box:
    """Return the component wise sum of two boxes of equal length."""
    if len(a) != len(b):
        raise ValueError("box_add: mismatched dimensions")
    return [a[i] + b[i] for i in range(len(a))]


def box_sub(a: Box, b: Box) -> Box:
    """Return the component wise difference of two boxes of equal length."""
    if len(a) != len(b):
        raise ValueError("box_sub: mismatched dimensions")
    return [a[i] - b[i] for i in range(len(a))]


def box_scalar_mul(alpha: float, box: Box) -> Box:
    """Scale every component of the box by the scalar alpha."""
    return [scalar_times_interval(alpha, iv) for iv in box]


def box_apply_unary(fn, box: Box) -> Box:
    """Apply a unary interval function to each component and return the new box."""
    return [fn(iv) for iv in box]


def box_contains(box: Box, point: Sequence[float]) -> bool:
    """Return True when every component of the point lies inside the matching interval."""
    if len(box) != len(point):
        return False
    return all(box[i].contains(float(point[i])) for i in range(len(box)))


def box_split(box: Box, idx: Optional[int] = None) -> Tuple[Box, Box]:
    """
    Split a box into two by bisecting one component.

    Parameters:
        box: the box to split.
        idx: index of the component to bisect. When None the widest component
             is chosen.

    Returns:
        A pair of boxes, the left holding the lower half of the chosen component
        and the right holding the upper half. This is the main tool for
        reducing overestimation.
    """
    if not box:
        raise ValueError("box_split: empty box")

    if idx is None:
        # Choose the widest component when no index is given.
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

# Convenience helpers used commonly in dynamics
def dynamic_pressure(rho: MaybeInterval, V: MaybeInterval) -> Interval:
    """
    Return the dynamic pressure interval, one half rho times V squared.

    Parameters:
        rho: density, as a float or interval.
        V:   speed, as a float or interval.
    """
    rho_i = promote(rho)
    V_i = promote(V)
    return Interval(0.5, 0.5) * rho_i * V_i.pow_int(2)


def clamp_interval(iv: Interval, lo: float, hi: float) -> Interval:
    """
    Clamp an interval to the range [lo, hi], the intersection of the two.

    Parameters:
        iv: the interval to clamp.
        lo: lower clamp bound.
        hi: upper clamp bound.

    Returns:
        The intersection interval. Raises ValueError when the intersection is
        empty.
    """
    new_lo = max(iv.lo, lo)
    new_hi = min(iv.hi, hi)
    if new_lo > new_hi:
        raise ValueError(f"clamp_interval empty intersection: {iv} ∩ [{lo},{hi}]")
    return Interval(new_lo, new_hi)

# interval euler time step for the ODEs
# def interval_euler_step(X, dt, f):
#     """
#     X : list[Interval]  (state box)
#     dt: float
#     f : function(t, X) -> list[Interval] (derivatives)
#     """
#     Xdot = f(X)
#     return box_add(X, box_scalar_mul(dt, Xdot))

def interval_euler_step(X, dt, f, t=None):
    """
    Take one interval Euler step.

    Parameters:
        X:  the state box, a list of intervals.
        dt: the time step.
        f:  the derivative function. It is called as f(t, X) when t is given,
            otherwise as f(X).
        t:  optional time passed to f.

    Returns:
        The next state box, X plus dt times the derivative box.
    """
    Xdot = f(t, X) if t is not None else f(X)
    return box_add(X, box_scalar_mul(dt, Xdot))



# self-test
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
