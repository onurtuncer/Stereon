#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""Generate corpus/planar-curves/: pairs of plane curves with a known answer.

Every unordered pair of the built-in plane curves, Line2d, Circle2d,
Ellipse2d, Hyperbola2d, Parabola2d and BSplineCurve2d (polynomial and
rational, counted as two kinds), in configurations whose every meeting
follows from a closed-form argument: crossings, tangencies of order two and
four, separations 3 2^-40 either side of a tangency (relative to the case's
coordinates), and coincident or overlapping arcs. The expected result never
comes from the kernel.

The argument. Every curve but a line is symmetric about the local x axis,
with its vertex on it, and is one of two shapes in x and Y = y^2:

    Y = q(x), q of degree <= 2       circles, ellipses, hyperbolas, parabolas,
                                     the rational arcs (exact circles,
                                     ellipses and hyperbolas) and the
                                     quadratic spline, which is a parabola
    x = p(Y), p of degree <= 2       the quartic spline x = A t^2 + B t^4,
                                     y = w t, and the quadratic spline again

so two of them meet where a polynomial of degree two in x (or in Y)
vanishes, and a line where one in its own parameter does. The roots are
taken in Q(sqrt D) with Fraction coefficients, the input doubles' exact
values, so every sign (a discriminant, Y at a root) is decided exactly:
Y > 0 gives two mirror points of the root's multiplicity m, Y = 0 one point
on the axis of multiplicity 2 m, Y < 0 none. Multiplicity 1 is a crossing,
2 and 4 tangencies (odd crosses, even does not). Branches, arcs and knots
are kept apart by float margins that the generator asserts.

Placements. "origin" and "far" (a translation by dyadic integers near
ADR-0001's size box) are exact on every input. "rotated" is the similarity
with matrix [[3, -4], [4, 3]], a rotation by atan(4/3) times 5, exact on
every point, radius and direction; an ellipse's, hyperbola's or parabola's
unit frame (3/5, 4/5) rounds, so regimes whose answer a rounding can change
(a tangency or coincidence with such a conic) are generated at the exact
placements only, as in generate_surface_corpus.py. A circle keeps the
identity frame everywhere: its point set does not depend on it. Each
placement is generated at two scales, 1 and 100.

Usage: generate_planar_corpus.py [--root <repository root>] [--check]
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass, field
from decimal import Decimal, getcontext
from fractions import Fraction as Fr
from pathlib import Path
from typing import Callable

getcontext().prec = 80

ADDED = "2026-11-12"
KINDS = ("line", "circle", "ellipse", "hyperbola", "parabola", "spline", "rational")
SCALES = (1, 100)


@dataclass(frozen=True)
class Placement:
    name: str
    matrix: tuple  # ((m00, m01), (m10, m11)), integers
    stretch: int  # the matrix's scale: |M v| = stretch |v|
    translation: tuple  # Fractions


PLACEMENTS = (
    Placement("origin", ((1, 0), (0, 1)), 1, (Fr(0), Fr(0))),
    Placement("far", ((1, 0), (0, 1)), 1, (Fr(393216), Fr(-262144))),
    Placement("rotated", ((3, -4), (4, 3)), 5, (Fr(12), Fr(-7))),
)
EXACT_ONLY = ("origin", "far")

SPLINE_TANGENCY = ("a tangency of two B-splines fails closed (v0.2 ledger, 2D curve-curve, "
                   "Remaining)")
SPLINE_ALONG_CONIC = "a B-spline along a conic fails closed (v0.2 ledger, 2D curve-curve, Remaining)"
SPLINE_ALONG_SPLINE = ("a B-spline along another B-spline other than an identical one fails "
                       "closed (v0.2 ledger, 2D curve-curve, Remaining)")
SPLINE_ALONG_CUT_LINE = ("a B-spline along a line that the requested domains cut fails closed "
                         "(v0.2 ledger, 2D curve-curve, Remaining)")

# Pair kinds whose every decision the ledger says is exact: their points
# must come Decided.
EXACT_KINDS = {("line", "line"), ("line", "circle"), ("line", "ellipse"), ("line", "hyperbola"),
               ("line", "parabola"), ("circle", "circle")}


# --- Exact numbers in Q(sqrt d) ----------------------------------------------


def sign(x: Fr) -> int:
    return (x > 0) - (x < 0)


@dataclass(frozen=True)
class Surd:
    """a + b sqrt(d), with Fraction a, b and d >= 0."""

    a: Fr
    b: Fr = Fr(0)
    d: Fr = Fr(0)

    def __add__(self, o):
        o = lift(o, self.d)
        return Surd(self.a + o.a, self.b + o.b, self.d if self.b else o.d)

    __radd__ = __add__

    def __neg__(self):
        return Surd(-self.a, -self.b, self.d)

    def __sub__(self, o):
        return self + (-lift(o, self.d))

    def __rsub__(self, o):
        return lift(o, self.d) - self

    def __mul__(self, o):
        o = lift(o, self.d)
        d = self.d if self.b else o.d
        if self.b and o.b:
            assert self.d == o.d
        return Surd(self.a * o.a + self.b * o.b * d, self.a * o.b + self.b * o.a, d)

    __rmul__ = __mul__

    def sign(self) -> int:
        if self.b == 0 or self.d == 0:
            return sign(self.a)
        sa, sb = sign(self.a), sign(self.b)
        if sa == 0 or sa == sb:
            return sb
        # Opposite signs: compare a^2 with b^2 d.
        c = sign(self.a * self.a - self.b * self.b * self.d)
        return sa if c > 0 else (sb if c < 0 else 0)

    def dec(self) -> Decimal:
        value = Decimal(self.a.numerator) / Decimal(self.a.denominator)
        if self.b and self.d:
            root = (Decimal(self.d.numerator) / Decimal(self.d.denominator)).sqrt()
            value += Decimal(self.b.numerator) / Decimal(self.b.denominator) * root
        return value


def lift(x, d=Fr(0)) -> Surd:
    if isinstance(x, Surd):
        return x
    return Surd(Fr(x), Fr(0), d)


def evaluate(coefficients, x):
    """c0 + c1 x + c2 x^2 at a Fraction or a Surd."""
    value = lift(Fr(0))
    for c in reversed(coefficients):
        value = value * x + c
    return value


def roots(coefficients):
    """The real roots of c0 + c1 x + c2 x^2, as (Surd, multiplicity), or
    None if the polynomial vanishes identically."""
    c0, c1, c2 = (list(coefficients) + [Fr(0)] * 3)[:3]
    if c2 == 0:
        if c1 == 0:
            return None if c0 == 0 else []
        return [(lift(-c0 / c1), 1)]
    disc = c1 * c1 - 4 * c2 * c0
    centre = -c1 / (2 * c2)
    if disc < 0:
        return []
    if disc == 0:
        return [(lift(centre), 2)]
    half = 1 / (2 * c2)
    return [(Surd(centre, -abs(half), disc), 1), (Surd(centre, abs(half), disc), 1)]


def dsqrt(x: Decimal) -> Decimal:
    return x.sqrt() if x > 0 else Decimal(0)


def _series(x: Decimal, odd: bool, alternating: bool) -> Decimal:
    """sin/cos (alternating) or sinh/cosh of x by their Taylor series, at the
    80-digit context: platform-independent, unlike the C library's."""
    term = x if odd else Decimal(1)
    total, k, x2 = term, 1 if odd else 0, x * x
    while True:
        term = term * x2 / ((k + 1) * (k + 2))
        if alternating:
            term = -term
        k += 2
        if abs(term) < Decimal(10) ** -70 * max(abs(total), Decimal(1)):
            return total + term
        total += term


def trig(name: str, t: float) -> float:
    """cos, sin, cosh or sinh of the double t, rounded once to the nearest
    double. math.cos and friends differ in the last bit between Windows,
    glibc and macOS, which made the corpus platform-dependent."""
    x = Decimal(t)
    odd = name in ("sin", "sinh")
    return float(_series(x, odd, alternating=name in ("sin", "cos")))


# --- Curves in the local frame ----------------------------------------------


@dataclass
class Curve:
    """A curve in the local frame: its exact relation, the data the kernel
    reads, and what keeps a point of the relation off the curve."""

    kind: str  # one of KINDS
    rel: str  # "Y" (Y = q(x)), "X" (x = p(Y)) or "line"
    coeffs: tuple = ()  # q or p, Fractions
    line: tuple = ()  # (origin, direction), Fraction pairs
    data: dict = field(default_factory=dict)  # local kernel data, see emit()
    keep: Callable | None = None  # margin(x, y) > 0 on the curve, < 0 off it
    knots_y: tuple = ()  # |y| of the spline knots, kept away from points
    frame: bool = False  # a unit frame that a rotation rounds
    domain: tuple | None = None
    scale: Fr = Fr(1)  # a local size, for the feature scale

    def tangent(self, x: float, y: float):
        if self.rel == "line":
            return (float(self.line[1][0]), float(self.line[1][1]))
        c = [float(v) for v in self.coeffs]
        if self.rel == "Y":
            # Gradient of y^2 - q(x) is (-q'(x), 2y): the tangent is (2y, q'(x)).
            return (2.0 * y, c[1] + 2.0 * c[2] * x)
        big = c[1] + 2.0 * c[2] * y * y
        return (2.0 * y * big, 1.0)


def as_y(c: Curve):
    """q with Y = q(x) for a curve of either shape, if it has one."""
    if c.rel == "Y":
        return c.coeffs
    p0, p1, p2 = c.coeffs
    assert p2 == 0 and p1 != 0, "a quartic only meets lines and curves x = p(Y)"
    return (-p0 / p1, 1 / p1, Fr(0))


def as_x(c: Curve):
    if c.rel == "X":
        return c.coeffs
    q0, q1, q2 = c.coeffs
    assert q2 == 0 and q1 != 0
    return (-q0 / q1, 1 / q1, Fr(0))


def frac(v) -> Fr:
    return v if isinstance(v, Fr) else Fr(v)


def pt(x, y):
    return (frac(x), frac(y))


def line(origin, direction, domain=None) -> Curve:
    o, d = pt(*origin), pt(*direction)
    return Curve("line", "line", line=(o, d), data={"type": "line", "origin": o, "direction": d},
                 domain=domain, scale=max(abs(o[0]), abs(o[1]), 1))


def vertical(v) -> Curve:
    """The line x = v, upwards from (v, -1) at speed 2."""
    return line((v, -1), (0, 2))


def circle(v, s, r, domain=None) -> Curve:
    r = frac(r)
    c = frac(v) + s * r
    return Curve("circle", "Y", (r * r - c * c, 2 * c, Fr(-1)),
                 data={"type": "circle", "center": (c, Fr(0)), "radius": r},
                 domain=domain, scale=abs(c) + r)


def ellipse(v, s, ax, ay, flip=False, domain=None, round_frame=True) -> Curve:
    """Semi-axes ax along x and ay along y, the vertex at x = v and the
    curve on side s of it. flip turns the frame half a turn."""
    ax, ay = frac(ax), frac(ay)
    c = frac(v) + s * ax
    k = ay * ay / (ax * ax)
    if ax >= ay:
        x, y = (1, 0), (0, 1)
        major, minor = ax, ay
    else:
        x, y = (0, 1), (-1, 0)
        major, minor = ay, ax
    if flip:
        x, y = (-x[0], -x[1]), (-y[0], -y[1])
    return Curve("ellipse", "Y", (ay * ay - k * c * c, 2 * k * c, -k),
                 data={"type": "ellipse", "center": (c, Fr(0)), "major": major, "minor": minor,
                       "x_axis": x, "y_axis": y},
                 frame=round_frame and ax != ay, domain=domain, scale=abs(c) + max(ax, ay))


def hyperbola(v, s, a, b, flip=False, domain=None) -> Curve:
    """The branch with its vertex at x = v, opening to side s; flip
    reverses its parameter (y_axis negated)."""
    a, b = frac(a), frac(b)
    c = frac(v) - s * a
    k = b * b / (a * a)
    side = (s, 0)
    return Curve("hyperbola", "Y", (k * c * c - b * b, -2 * k * c, k),
                 data={"type": "hyperbola", "center": (c, Fr(0)), "major": a, "minor": b,
                       "x_axis": side, "y_axis": (0, -1 if flip else 1)},
                 keep=lambda x, y, c=float(c), s=s: s * (x - c), frame=True, domain=domain,
                 scale=abs(c) + a + b)


def parabola(v, s, f, flip=False, domain=None) -> Curve:
    f, v = frac(f), frac(v)
    return Curve("parabola", "Y", (-4 * f * s * v, 4 * f * s, Fr(0)),
                 data={"type": "parabola", "vertex": (v, Fr(0)), "focal": f, "x_axis": (s, 0),
                       "y_axis": (0, -1 if flip else 1)},
                 frame=True, domain=domain, scale=abs(v) + 4 * f)


# Bernstein and B-spline arithmetic, exact.


def bernstein(power):
    """Bernstein coefficients on [0, 1] of sum power[j] s^j."""
    n = len(power) - 1
    return [sum(Fr(math.comb(i, j), math.comb(n, j)) * power[j] for j in range(i + 1))
            for i in range(n + 1)]


def insert_knot(degree, knots, poles, weights, u):
    """Boehm's knot insertion, on homogeneous poles, exact."""
    k = max(i for i in range(len(knots) - 1) if knots[i] <= u < knots[i + 1])
    hom = [(p[0] * w, p[1] * w, w) for p, w in zip(poles, weights)]
    new = []
    for i in range(len(hom) + 1):
        if i <= k - degree:
            new.append(hom[i])
        elif i > k:
            new.append(hom[i - 1])
        else:
            alpha = (u - knots[i]) / (knots[i + degree] - knots[i])
            new.append(tuple(alpha * hom[i][j] + (1 - alpha) * hom[i - 1][j] for j in range(3)))
    knots = sorted(knots + [u])
    return knots, [(h[0] / h[2], h[1] / h[2]) for h in new], [h[2] for h in new]


def de_boor(degree, knots, poles, weights, u):
    """The point at u, exact."""
    k = max(i for i in range(degree, len(knots) - degree - 1) if knots[i] <= u)
    hom = [(poles[j][0] * weights[j], poles[j][1] * weights[j], weights[j])
           for j in range(k - degree, k + 1)]
    for r in range(1, degree + 1):
        for j in range(degree, r - 1, -1):
            i = j + k - degree
            alpha = (u - knots[i]) / (knots[i + degree - r + 1] - knots[i])
            hom[j] = tuple(alpha * hom[j][m] + (1 - alpha) * hom[j - 1][m] for m in range(3))
    h = hom[degree]
    return (h[0] / h[2], h[1] / h[2])


def bspline(kind, rel, coeffs, degree, knots, poles, weights=None, keep=None, knots_y=(),
            domain=None) -> Curve:
    weights = weights or [Fr(1)] * len(poles)
    data = {"type": "bspline", "degree": degree, "knots": [frac(k) for k in knots],
            "poles": [pt(*p) for p in poles]}
    if kind == "rational":
        data["weights"] = [frac(w) for w in weights]
    size = max(max(abs(p[0]), abs(p[1])) for p in data["poles"])
    return Curve(kind, rel, tuple(coeffs), data=data, keep=keep, knots_y=tuple(knots_y),
                 domain=domain, scale=size)


def spline(v, s, a, w, extra_knot=False, domain=None) -> Curve:
    """The parabola x = v + s a t^2, y = w t, t = 2u - 1 in [-1, 1], as a
    quadratic B-spline with knots inserted at 1/4 and 3/4 (|y| = w/2)."""
    v, a, w = frac(v), frac(a), frac(w)
    poles = [(v + s * a, -w), (v - s * a, Fr(0)), (v + s * a, w)]
    knots = [Fr(0)] * 3 + [Fr(1)] * 3
    ones = [Fr(1)] * 3
    knots, poles, ones = insert_knot(2, knots, poles, ones, Fr(1, 4))
    knots, poles, ones = insert_knot(2, knots, poles, ones, Fr(3, 4))
    kys = [w / 2]
    if extra_knot:
        knots, poles, ones = insert_knot(2, knots, poles, ones, Fr(1, 2))
        kys.append(Fr(0))
    return bspline("spline", "X", (v, s * a / (w * w), Fr(0)), 2, knots, poles,
                   keep=lambda x, y, w=float(w): w - abs(y), knots_y=kys, domain=domain)


def quartic(v, s, a, b, w) -> Curve:
    """x = v + s (a t^2 + b t^4), y = w t, t = 3 (2u - 1) in [-3, 3]: a
    quartic Bezier. Over [-1, 1] the poles of t^2 would have thirds; over
    [-3, 3] they are 9 (1, 0, -1/3, 0, 1), integers."""
    v, a, b, w = frac(v), frac(a), frac(b), frac(w)
    t2 = [Fr(9), Fr(-36), Fr(36), Fr(0), Fr(0)]
    t4 = [sum(t2[i] * t2[k - i] for i in range(k + 1)) for k in range(5)]
    xs = bernstein([(v if j == 0 else 0) + s * (a * t2[j] + b * t4[j]) for j in range(5)])
    ys = bernstein([-3 * w, 6 * w, Fr(0), Fr(0), Fr(0)])
    knots = [Fr(0)] * 5 + [Fr(1)] * 5
    return bspline("spline", "X", (v, s * a / (w * w), s * b / w ** 4), 4, knots,
                   list(zip(xs, ys)), keep=lambda x, y, h=float(3 * w): h - abs(y))


def straight_spline(v) -> Curve:
    v = frac(v)
    poles = [(v, Fr(-4)), (v, Fr(-1)), (v, Fr(1)), (v, Fr(4))]
    knots = [Fr(0)] * 3 + [Fr(1, 2)] + [Fr(1)] * 3
    c = bspline("spline", "line", (), 2, knots, poles)
    c.line = ((v, Fr(-4)), (Fr(0), Fr(1)))
    return c


def straight_rational(v) -> Curve:
    v = frac(v)
    poles = [(v, Fr(-4)), (v, Fr(0)), (v, Fr(4))]
    c = bspline("rational", "line", (), 2, [Fr(0)] * 3 + [Fr(1)] * 3, poles,
                [Fr(1), Fr(2), Fr(1)])
    c.line = ((v, Fr(-4)), (Fr(0), Fr(1)))
    return c


def rational(v, s, k, sx=1, sy=1, domain=None) -> Curve:
    """The arc |y| <= 12 k sy of the ellipse with semi-axes 15 k sx, 15 k sy
    (a circle when sx = sy) whose vertex is at x = v, on side s: a rational
    quadratic Bezier with weights 5, 3, 5, exact since cos(theta) = 3/5 for
    the half angle with sin(theta) = 4/5."""
    v, k, sx, sy = frac(v), frac(k), frac(sx), frac(sy)
    ax, ay = 15 * k * sx, 15 * k * sy
    c = v + s * ax
    poles = [(c - s * 9 * k * sx, -12 * k * sy), (c - s * 25 * k * sx, Fr(0)),
             (c - s * 9 * k * sx, 12 * k * sy)]
    q = ay * ay / (ax * ax)
    return bspline("rational", "Y", (ay * ay - q * c * c, 2 * q * c, -q), 2,
                   [Fr(0)] * 3 + [Fr(1)] * 3, poles, [Fr(5), Fr(3), Fr(5)],
                   keep=lambda x, y, h=float(12 * k * sy): h - abs(y), domain=domain)


def rational_hyperbola(v, s, m) -> Curve:
    """The arc |y| <= 15 m of the branch x^2 - y^2 = (20 m)^2 with its
    vertex at x = v, on side s: poles (25, -15), (16, 0), (25, 15) times m
    about the centre, weights 1, 5/4, 1, the conic through both ends, tangent
    there, and through the vertex."""
    v, m = frac(v), frac(m)
    r = 20 * m
    c = v - s * r
    poles = [(c + s * 25 * m, -15 * m), (c + s * 16 * m, Fr(0)), (c + s * 25 * m, 15 * m)]
    return bspline("rational", "Y", (c * c - r * r, -2 * c, Fr(1)), 2, [Fr(0)] * 3 + [Fr(1)] * 3,
                   poles, [Fr(1), Fr(5, 4), Fr(1)],
                   keep=lambda x, y, h=float(15 * m): h - abs(y))


# --- The oracle ---------------------------------------------------------------


class Coincident(Exception):
    pass


def meet(a: Curve, b: Curve):
    """Every point of both relations, as (x Surd, y Decimal, multiplicity),
    before the curves' own restrictions."""
    if a.rel == "line" and b.rel == "line":
        (o1, d1), (o2, d2) = a.line, b.line
        cross = d1[0] * d2[1] - d1[1] * d2[0]
        off = (o2[0] - o1[0], o2[1] - o1[1])
        if cross == 0:
            if off[0] * d1[1] - off[1] * d1[0] == 0:
                raise Coincident
            return []
        t = (off[0] * d2[1] - off[1] * d2[0]) / cross
        return [(lift(o1[0] + t * d1[0]), lift(o1[1] + t * d1[1]).dec(), 1)]
    if b.rel == "line":
        a, b = b, a
    if a.rel == "line":
        (px, py), (dx, dy) = a.line
        if b.rel == "Y" or b.coeffs[2] == 0:
            q0, q1, q2 = as_y(b)
            poly = (py * py - q0 - q1 * px - q2 * px * px,
                    2 * py * dy - q1 * dx - 2 * q2 * px * dx,
                    dy * dy - q2 * dx * dx)
            rs = roots(poly)
            if rs is None:
                raise Coincident
            return [(px + t * dx, (py + t * dy).dec(), m) for t, m in rs]
        p = b.coeffs
        if dx == 0:
            return mirror_y(roots((p[0] - px, p[1], p[2])), p)
        assert dy == 0, "a quartic only meets vertical and horizontal lines"
        y2 = py * py
        return [(evaluate(p, lift(y2)), Decimal(py.numerator) / Decimal(py.denominator), 1)]
    as_q = lambda c: c.rel == "Y" or c.coeffs[2] == 0  # noqa: E731
    as_p = lambda c: c.rel == "X" or c.coeffs[2] == 0  # noqa: E731
    if as_q(a) and as_q(b):
        qa, qb = as_y(a), as_y(b)
        rs = roots(tuple(u - v for u, v in zip(qa, qb)))
        if rs is None:
            raise Coincident
        out = []
        for x, m in rs:
            big_y = evaluate(qa, x)
            out += mirror(x, big_y, m)
        return out
    assert as_p(a) and as_p(b), "a quartic only meets lines and curves x = p(Y)"
    pa, pb = as_x(a), as_x(b)
    p = tuple(u - v for u, v in zip(pa, pb))
    return mirror_y(roots(p), pa)


def mirror(x: Surd, big_y: Surd, m: int):
    s = big_y.sign()
    if s < 0:
        return []
    if s == 0:
        return [(x, Decimal(0), 2 * m)]
    y = dsqrt(big_y.dec())
    return [(x, -y, m), (x, y, m)]


def mirror_y(rs, p):
    if rs is None:
        raise Coincident
    out = []
    for big_y, m in rs:
        s = big_y.sign()
        if s < 0:
            continue
        x = evaluate(p, big_y)
        if s == 0:
            out.append((x, Decimal(0), 2 * m))
        else:
            y = dsqrt(big_y.dec())
            out += [(x, -y, m), (x, y, m)]
    return out


def solve(a: Curve, b: Curve, margin: float):
    """The points of both curves, (x, y, multiplicity, sin of the angle),
    each decided away from every restriction by `margin`."""
    out = []
    for x, y, m in meet(a, b):
        xf, yf = float(x.dec()), float(y)
        inside = True
        for c in (a, b):
            if c.keep is not None:
                k = c.keep(xf, yf)
                assert abs(k) > margin, f"a point at ({xf}, {yf}) is within {margin} of an end"
                inside = inside and k > 0
            for ky in c.knots_y:
                assert abs(abs(yf) - float(ky)) > margin, f"a point at y = {yf} is near a knot"
        if not inside:
            continue
        ta, tb = a.tangent(xf, yf), b.tangent(xf, yf)
        cross = abs(ta[0] * tb[1] - ta[1] * tb[0]) / (math.hypot(*ta) * math.hypot(*tb))
        out.append((x.dec(), y, m, cross))
    return out


# --- Placement -------------------------------------------------------------


@dataclass
class Frame:
    placement: Placement
    scale: int

    def point(self, p):
        (m00, m01), (m10, m11) = self.placement.matrix
        t = self.placement.translation
        s = self.scale
        return (t[0] + s * (m00 * p[0] + m01 * p[1]), t[1] + s * (m10 * p[0] + m11 * p[1]))

    def point_dec(self, x: Decimal, y: Decimal):
        (m00, m01), (m10, m11) = self.placement.matrix
        t = self.placement.translation
        s = self.scale
        tx = Decimal(t[0].numerator) / Decimal(t[0].denominator)
        ty = Decimal(t[1].numerator) / Decimal(t[1].denominator)
        return (tx + s * (m00 * x + m01 * y), ty + s * (m10 * x + m11 * y))

    def vector(self, d):
        (m00, m01), (m10, m11) = self.placement.matrix
        return (self.scale * (m00 * d[0] + m01 * d[1]), self.scale * (m10 * d[0] + m11 * d[1]))

    def length(self, r):
        return self.scale * self.placement.stretch * r

    def unit(self, d):
        (m00, m01), (m10, m11) = self.placement.matrix
        k = self.placement.stretch
        return (Fr(m00 * d[0] + m01 * d[1], k), Fr(m10 * d[0] + m11 * d[1], k))


def exact(v: Fr) -> float:
    f = float(v)
    assert Fr(f) == v, f"{v} is not a double"
    return f


def nearest(v: Fr) -> float:
    return float(v)


def emit(c: Curve, name: str, frame: Frame) -> dict:
    d = c.data
    out = {"name": name, "type": d["type"]}
    pair = lambda p: [exact(u) for u in frame.point(p)]  # noqa: E731
    if d["type"] == "line":
        out["origin"] = pair(d["origin"])
        out["direction"] = [exact(u) for u in frame.vector(d["direction"])]
    elif d["type"] == "bspline":
        out["degree"] = d["degree"]
        out["knots"] = [exact(k) for k in d["knots"]]
        out["poles"] = [pair(p) for p in d["poles"]]
        if "weights" in d:
            out["weights"] = [exact(w) for w in d["weights"]]
    else:
        key = "vertex" if d["type"] == "parabola" else "center"
        out[key] = pair(d[key])
        if d["type"] == "circle":
            out["radius"] = exact(frame.length(d["radius"]))
            out["x_axis"], out["y_axis"] = [1.0, 0.0], [0.0, 1.0]
        else:
            if d["type"] == "parabola":
                out["focal"] = exact(frame.length(d["focal"]))
            else:
                out["major"] = exact(frame.length(d["major"]))
                out["minor"] = exact(frame.length(d["minor"]))
            if c.frame:
                out["x_axis"] = [nearest(u) for u in frame.unit(d["x_axis"])]
                out["y_axis"] = [nearest(u) for u in frame.unit(d["y_axis"])]
            else:
                out["x_axis"] = [float(u) for u in d["x_axis"]]
                out["y_axis"] = [float(u) for u in d["y_axis"]]
    if c.domain is not None:
        out["domain"] = [exact(frac(u)) for u in c.domain]
    return out


def placed_eval(c: dict, t: float):
    """A placed analytic curve's or line's point at t, in floats."""
    kind = c["type"]
    if kind == "line":
        return [c["origin"][i] + t * c["direction"][i] for i in range(2)]
    x, y = c["x_axis"], c["y_axis"]
    if kind == "circle":
        o, a, b = c["center"], c["radius"] * trig("cos", t), c["radius"] * trig("sin", t)
    elif kind == "ellipse":
        o, a, b = c["center"], c["major"] * trig("cos", t), c["minor"] * trig("sin", t)
    elif kind == "hyperbola":
        o, a, b = c["center"], c["major"] * trig("cosh", t), c["minor"] * trig("sinh", t)
    else:
        o, a, b = c["vertex"], t * t / (4.0 * c["focal"]), t
    return [o[i] + a * x[i] + b * y[i] for i in range(2)]


# --- Regimes ---------------------------------------------------------------


@dataclass
class Regime:
    name: str
    a: Curve
    b: Curve
    overlap: Callable | None = None  # (placed A, placed B, frame) -> [(from, to)]
    fail: str | None = None  # may_fail_closed reason
    placements: tuple = ("origin", "far", "rotated")


def regimes_for(ka: str, kb: str, delta: Fr):
    """The regimes of the unordered pair (ka, kb), ka before kb in KINDS."""
    out = []
    atom = ATOMS
    pair = (ka, kb)
    splines = sum(k in ("spline", "rational") for k in pair)
    tangent_fail = SPLINE_TANGENCY if splines == 2 else None

    def exact_only(*curves):
        return EXACT_ONLY if any(c.frame for c in curves) else ("origin", "far", "rotated")

    if pair == ("line", "line"):
        a = vertical(0)
        out.append(Regime("crossing", a, line((1, -1), (1, 3))))
        out.append(Regime("apart", a, vertical(delta)))
        a = line((0, -1), (0, 2), domain=(0, 4))
        b = line((0, 5), (0, -4), domain=(0, 2))
        out.append(Regime("overlapping-segments", a, b,
                          overlap=lambda A, B, f: [(placed_eval(A, 0.0), placed_eval(A, 3.0))]))
        return out

    # Crossing: B's vertex moved into A, or a slanted line through A.
    a = atom[ka](0, 1)
    if ka == "line":
        b = atom[kb](0, 1)
        out.append(Regime("crossing", line((1, -1), (1, 3)), b))
    else:
        b = atom[kb](SHIFT.get(pair, Fr(5, 4)), -1)
        out.append(Regime("crossing", a, b, fail=None))
    # Tangency at the common vertex, B on the other side.
    b0 = atom[kb](0, -1 if ka != "line" else 1)
    out.append(Regime("tangent", a, b0, fail=tangent_fail, placements=exact_only(a, b0)))
    # 3 2^-40 either side of it: apart, and two crossings.
    away = -delta if ka != "line" else delta
    out.append(Regime("apart", a, atom[kb](away, b0_side(ka)),
                      fail=None, placements=("origin", "far", "rotated")))
    out.append(Regime("near", a, atom[kb](-away, b0_side(ka)),
                      fail=None, placements=("origin", "far", "rotated")))
    # Order four: the same side, the same curvature at the vertex.
    if pair in OSCULATING:
        ca, cb = OSCULATING[pair]()
        out.append(Regime("osculating", ca, cb, fail=tangent_fail,
                          placements=exact_only(ca, cb)))
    # Internal tangency of order two, where the pair has no other.
    if pair in INTERNAL:
        ca, cb = INTERNAL[pair]()
        out.append(Regime("tangent-internal", ca, cb, fail=tangent_fail,
                          placements=exact_only(ca, cb)))
    if pair in COINCIDENT:
        for r in COINCIDENT[pair]():
            if r.placements == ("origin", "far", "rotated"):
                r.placements = exact_only(r.a, r.b)
            out.append(r)
    return out


def b0_side(ka):
    return 1 if ka == "line" else -1


ATOMS = {
    "line": lambda v, s: vertical(v),
    "circle": lambda v, s: circle(v, s, 4),
    "ellipse": lambda v, s: ellipse(v, s, 6, 3),
    "hyperbola": lambda v, s: hyperbola(v, s, 2, 3),
    "parabola": lambda v, s: parabola(v, s, 1),
    "spline": lambda v, s: spline(v, s, 2, 4),
    "rational": lambda v, s: rational(v, s, Fr(1, 2)),
}

# Where B's vertex goes for a crossing, when 5/4 puts a point near a knot or an
# arc's end.
SHIFT: dict = {}

# Curvature 1/8 at the vertex: C(8), E(32, 16), E(128, 32), H(2, 4), H(8, 8),
# P(4), the spline a = 1, w = 4 (curvature 2a / w^2), the quartic with the
# same a; curvature 2/15: C(7.5), E(30, 15), H(7.5, 7.5), P(3.75), the
# rational circle k = 1/2 (radius 7.5) and the rational ellipse k = 1/2,
# sx = 4, sy = 2 (30 by 15).
OSCULATING = {
    ("line", "spline"): lambda: (vertical(0), quartic(0, 1, 0, Fr(1, 16), 1)),
    ("circle", "ellipse"): lambda: (circle(0, 1, 8), ellipse(0, 1, 32, 16)),
    ("circle", "hyperbola"): lambda: (circle(0, 1, 8), hyperbola(0, 1, 2, 4)),
    ("circle", "parabola"): lambda: (circle(0, 1, 8), parabola(0, 1, 4)),
    ("circle", "spline"): lambda: (circle(0, 1, 8), spline(0, 1, 1, 4)),
    ("circle", "rational"): lambda: (circle(0, 1, Fr(15, 2)), rational(0, 1, Fr(1, 2), 4, 2)),
    ("ellipse", "ellipse"): lambda: (ellipse(0, 1, 32, 16), ellipse(0, 1, 128, 32)),
    ("ellipse", "hyperbola"): lambda: (ellipse(0, 1, 32, 16), hyperbola(0, 1, 2, 4)),
    ("ellipse", "parabola"): lambda: (ellipse(0, 1, 32, 16), parabola(0, 1, 4)),
    ("ellipse", "spline"): lambda: (ellipse(0, 1, 32, 16), spline(0, 1, 1, 4)),
    ("ellipse", "rational"): lambda: (ellipse(0, 1, 30, 15), rational(0, 1, Fr(1, 2))),
    ("hyperbola", "hyperbola"): lambda: (hyperbola(0, 1, 2, 4), hyperbola(0, 1, 8, 8)),
    ("hyperbola", "parabola"): lambda: (hyperbola(0, 1, 2, 4), parabola(0, 1, 4)),
    ("hyperbola", "spline"): lambda: (hyperbola(0, 1, 2, 4), spline(0, 1, 1, 4)),
    ("hyperbola", "rational"): lambda: (hyperbola(0, 1, Fr(15, 2), Fr(15, 2)),
                                        rational(0, 1, Fr(1, 2))),
    ("parabola", "spline"): lambda: (parabola(0, 1, 4), quartic(0, 1, 1, Fr(1, 128), 4)),
    ("parabola", "rational"): lambda: (parabola(0, 1, Fr(15, 4)), rational(0, 1, Fr(1, 2))),
    ("spline", "spline"): lambda: (spline(0, 1, 1, 4), quartic(0, 1, 1, Fr(1, 128), 4)),
    ("rational", "rational"): lambda: (rational(0, 1, Fr(1, 2)), rational(0, 1, Fr(1, 2), 4, 2)),
}

INTERNAL = {
    ("circle", "circle"): lambda: (circle(0, 1, 4), circle(0, 1, 2)),
    ("parabola", "parabola"): lambda: (parabola(0, 1, 1), parabola(0, 1, 2)),
}


def ends_of(curve_key):
    """The overlap's ends: B's first and last poles (the whole spline)."""
    def ends(A, B, f):
        c = A if A["type"] == "bspline" and curve_key == "A" else B
        return [(c["poles"][0], c["poles"][-1])]
    return ends


def at(curve_key, lo, hi):
    def ends(A, B, f):
        c = A if curve_key == "A" else B
        return [(placed_eval(c, lo), placed_eval(c, hi))]
    return ends


def spline_at(curve_key, lo, hi, local):
    """Ends of an identical spline's overlap: the exact points at lo and hi
    of the local curve, placed."""
    def ends(A, B, f):
        d = local.data
        w = d.get("weights", [Fr(1)] * len(d["poles"]))
        p = [f.point(de_boor(d["degree"], d["knots"], d["poles"], w, frac(u))) for u in (lo, hi)]
        return [([float(p[0][0]), float(p[0][1])], [float(p[1][0]), float(p[1][1])])]
    return ends


def coincident_spline_spline():
    a = spline(0, 1, 2, 4, domain=(Fr(3, 8), Fr(7, 8)))
    b = spline(0, 1, 2, 4, domain=(Fr(1, 8), Fr(5, 8)))
    same = spline(0, 1, 2, 4)
    refined = spline(0, 1, 2, 4, extra_knot=True)
    return [Regime("identical-arcs", a, b, overlap=spline_at("A", Fr(3, 8), Fr(5, 8), same)),
            Regime("refined-copy", same, refined, overlap=ends_of("A"), fail=SPLINE_ALONG_SPLINE)]


def coincident_rational_rational():
    a = rational(0, 1, Fr(1, 2), domain=(Fr(1, 8), Fr(5, 8)))
    b = rational(0, 1, Fr(1, 2), domain=(Fr(3, 8), Fr(7, 8)))
    same = rational(0, 1, Fr(1, 2))
    return [Regime("identical-arcs", a, b, overlap=spline_at("A", Fr(3, 8), Fr(5, 8), same))]


COINCIDENT = {
    ("line", "spline"): lambda: [
        Regime("along", vertical(0), straight_spline(0), overlap=ends_of("B")),
        Regime("along-cut", line((0, -1), (0, 2), domain=(0, 1)), straight_spline(0),
               overlap=lambda A, B, f: [(placed_eval(A, 0.0), placed_eval(A, 1.0))],
               fail=SPLINE_ALONG_CUT_LINE)],
    ("line", "rational"): lambda: [
        Regime("along", vertical(0), straight_rational(0), overlap=ends_of("B"))],
    ("circle", "circle"): lambda: [
        Regime("overlapping-arcs", circle(0, 1, 4, domain=(1, 3)), circle(0, 1, 4, domain=(2, 4)),
               overlap=at("A", 2.0, 3.0))],
    ("circle", "ellipse"): lambda: [
        Regime("round-ellipse", circle(0, 1, 4),
               ellipse(0, 1, 4, 4, domain=(1, Fr(5, 2)), round_frame=False),
               overlap=at("B", 1.0, 2.5))],
    ("circle", "rational"): lambda: [
        Regime("along", circle(0, 1, Fr(15, 2)), rational(0, 1, Fr(1, 2)), overlap=ends_of("B"),
               fail=SPLINE_ALONG_CONIC)],
    ("ellipse", "ellipse"): lambda: [
        Regime("half-turn-arcs", ellipse(0, 1, 6, 3, domain=(Fr(1, 2), Fr(5, 2))),
               ellipse(0, 1, 6, 3, flip=True, domain=(4, 5)), overlap=at("B", 4.0, 5.0))],
    ("ellipse", "rational"): lambda: [
        Regime("along", ellipse(0, 1, 30, 15), rational(0, 1, Fr(1, 2), 4, 2),
               overlap=ends_of("B"), fail=SPLINE_ALONG_CONIC)],
    ("hyperbola", "hyperbola"): lambda: [
        Regime("reversed-arcs", hyperbola(0, 1, 2, 3, domain=(-1, 2)),
               hyperbola(0, 1, 2, 3, flip=True, domain=(0, 5)), overlap=at("A", -1.0, 0.0))],
    ("hyperbola", "rational"): lambda: [
        Regime("along", hyperbola(0, 1, 5, 5), rational_hyperbola(0, 1, Fr(1, 4)),
               overlap=ends_of("B"), fail=SPLINE_ALONG_CONIC)],
    ("parabola", "parabola"): lambda: [
        Regime("reversed-arcs", parabola(0, 1, 1, domain=(-2, 3)),
               parabola(0, 1, 1, flip=True, domain=(-1, 4)), overlap=at("A", -2.0, 1.0))],
    ("parabola", "spline"): lambda: [
        Regime("along", parabola(0, 1, 2), spline(0, 1, 2, 4), overlap=ends_of("B"),
               fail=SPLINE_ALONG_CONIC)],
    ("spline", "spline"): coincident_spline_spline,
    ("rational", "rational"): coincident_rational_rational,
}


# --- Cases -----------------------------------------------------------------


def power_of_two_above(x: Fr) -> Fr:
    e = 0
    while Fr(2) ** e < x:
        e += 1
    while e > -60 and Fr(2) ** (e - 1) >= x:
        e -= 1
    return Fr(2) ** e


def magnitude(curves, frame):
    """The power of two above every placed coordinate and length."""
    values = [Fr(1)]
    for c in curves:
        d = c.data
        for key in ("origin", "center", "vertex"):
            if key in d:
                values += [abs(u) for u in frame.point(d[key])]
        for key in ("radius", "major", "minor", "focal"):
            if key in d:
                values.append(abs(frame.length(d[key])))
        for p in d.get("poles", []):
            values += [abs(u) for u in frame.point(p)]
    return power_of_two_above(max(values))


def local_size(curves):
    return power_of_two_above(max(Fr(abs(c.scale)) for c in curves))


def contact_of(m):
    return "crossing" if m % 2 else "tangential"


def build(ka, kb, r: Regime, frame: Frame, big: Fr):
    feature = float(local_size([r.a, r.b]) * frame.scale * frame.placement.stretch)
    ulp = float(big) * 2.0 ** -52
    names = ("A", "B")
    A, B = emit(r.a, "A", frame), emit(r.b, "B", frame)
    expected = {"points": [], "overlaps": []}
    try:
        pts = solve(r.a, r.b, 2.0 ** -12 * float(local_size([r.a, r.b])))
        assert r.overlap is None, f"{r.name}: an overlap was expected, the relations differ"
    except Coincident:
        assert r.overlap is not None, f"{r.name}: the curves coincide"
        pts = []
        within = 2.0 ** -30 * feature + 2.0 ** 8 * ulp
        for start, end in r.overlap(A, B, frame):
            expected["overlaps"].append({"from": [float(u) for u in start],
                                         "to": [float(u) for u in end], "within": within})
    for x, y, m, cross in pts:
        if m == 1:
            within = 2.0 ** 8 * ulp / max(cross, 1e-300) + 2.0 ** -40 * feature
        elif m == 2:
            within = 2.0 ** 8 * math.sqrt(ulp * feature)
        else:
            within = 2.0 ** 4 * (ulp * feature ** 3) ** 0.25
        gx, gy = frame.point_dec(x, y)
        expected["points"].append({"position": [float(gx), float(gy)], "contact": contact_of(m),
                                   "order": m, "within": within})
    ps = expected["points"]
    for i, p in enumerate(ps):
        for q in ps[i + 1:]:
            gap = math.dist(p["position"], q["position"])
            assert 4.0 * max(p["within"], q["within"]) < gap, \
                f"{r.name}: points {p['position']} and {q['position']} within reach of each other"
    if (ka, kb) in EXACT_KINDS:
        expected["decided"] = True
    if r.fail:
        expected["may_fail_closed"] = r.fail
    return [A, B], expected


def describe(ka, kb, r):
    what = {
        "crossing": "cross",
        "tangent": "touch at a common vertex from either side (order two)",
        "tangent-internal": "touch at a common vertex from one side (order two)",
        "osculating": "touch at a common vertex with equal curvature (order four)",
        "apart": "are 3 2^-40 (relative) apart from a tangency",
        "near": "cross twice 3 2^-40 (relative) past a tangency",
    }.get(r.name, f"coincide ({r.name})")
    return f"A {ka} and a {kb} that {what}"


def case_for(ka, kb, regime_name, placement, scale):
    frame = Frame(placement, scale)
    # The size first, with no shift, then the shift 3 2^-40 of it.
    r0 = next(r for r in regimes_for(ka, kb, Fr(0)) if r.name == regime_name)
    big = magnitude([r0.a, r0.b], frame)
    delta = 3 * Fr(2) ** -40 * power_of_two_above(big / (scale * placement.stretch))
    r = next(r for r in regimes_for(ka, kb, delta) if r.name == regime_name)
    if placement.name not in r.placements:
        return None
    curves, expected = build(ka, kb, r, frame, big)
    family = f"{ka}-{kb}"
    tags = ["generated", "planar", r.name]
    if expected.get("may_fail_closed"):
        tags.append("may-fail-closed")
    return {
        "format": "stereon-corpus-case/1",
        "id": f"planar-curves/{family}/{r.name}-{placement.name}-x{scale}",
        "kind": "reference" if r.name == "crossing" else "hard",
        "summary": f"{describe(ka, kb, r)}, placed at {placement.name}, scaled by {scale}.",
        "added": ADDED,
        "source": "Generated by tools/generate_planar_corpus.py",
        "licence": "MPL-2.0",
        "references": ["ADR-0029", "ADR-0053"],
        "tags": tags,
        "inputs": {"curves": curves, "scale": scale, "magnitude": float(big)},
        "operation": {"name": "intersect", "operands": ["A", "B"]},
        "expected": expected,
    }


def all_cases():
    for i, ka in enumerate(KINDS):
        for kb in KINDS[i:]:
            for r in regimes_for(ka, kb, Fr(0)):
                for placement in PLACEMENTS:
                    for scale in SCALES:
                        case = case_for(ka, kb, r.name, placement, scale)
                        if case is not None:
                            yield case


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument("--check", action="store_true", help="fail if a file differs")
    args = parser.parse_args()
    stale = []
    count = 0
    expected_paths = set()
    for case in all_cases():
        count += 1
        path = args.root / "corpus" / (case["id"] + ".json")
        expected_paths.add(path.resolve())
        text = json.dumps(case, indent=2) + "\n"
        if args.check:
            if not path.exists() or path.read_text(encoding="utf-8") != text:
                stale.append(case["id"])
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(text, encoding="utf-8", newline="\n")
    suite = args.root / "corpus" / "planar-curves"
    if suite.exists():
        for path in sorted(suite.glob("*/*.json")):
            if path.resolve() not in expected_paths:
                if args.check:
                    stale.append(str(path.relative_to(args.root)))
                else:
                    path.unlink()
    if stale:
        print(f"generate_planar_corpus: {len(stale)} stale case(s), first {stale[0]}", file=sys.stderr)
        return 1
    print(f"generate_planar_corpus: {count} cases {'match' if args.check else 'written'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
