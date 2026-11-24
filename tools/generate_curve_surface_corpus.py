#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""Generate corpus/curve-surface/: curves crossing surfaces at a known point.

For every pair of a curve kind and a surface kind, a point P of the surface
is taken at chart parameters away from its seams and poles, with its unit
normal n, both from the surface's defining formula evaluated here. The curve
is moved rigidly so that its point at a chosen parameter t* is P and its
tangent there is d = sin(b) n + cos(b) e, e a unit tangent of the surface:
the curve crosses the surface at P at the angle b to the tangent plane. Three
angles are generated: 90 degrees (along the normal), 20 degrees and 1 degree
(grazing, still transverse).

The claim of a case is local: the result has a decided crossing at P whose
curve parameter is t*. The pair may also meet elsewhere.

Usage: generate_curve_surface_corpus.py [--root <repository root>] [--check]
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

ANGLES = {"normal": 90.0, "oblique": 20.0, "grazing": 1.0}
CURVES = ("line", "circle", "ellipse", "parabola", "hyperbola", "bspline", "nurbs")
SURFACES = ("plane", "cylinder", "cone", "sphere", "torus", "ellipsoid", "bspline", "nurbs")

# The frame every analytic surface is given in: orthonormal, its axis
# X x Y = (0, -0.8, 0.6).
SX = (1.0, 0.0, 0.0)
SY = (0.0, 0.6, 0.8)
CENTER = (0.5, -0.25, 1.0)
# The direction the curves' planes are tilted towards, across d.
TILT = (0.3, -0.5, 0.8)
# How the surface tangent e is turned from the first chart direction.
TURN = 0.7

# The curves in their own plane, small beside the surfaces.
BEZIER = ((-1.0, 0.0), (-0.5, 0.75), (0.5, -0.75), (1.0, 0.0))
QUARTER = ((1.0, 0.0), (1.0, 1.0), (0.0, 1.0))
QUARTER_WEIGHTS = (1.0, math.sqrt(0.5), 1.0)
CURVE_PARAMETER = {
    "line": 0.5,
    "circle": 0.8,
    "ellipse": 2.0,
    "parabola": -0.6,
    "hyperbola": 0.4,
    "bspline": 0.4,
    "nurbs": 0.3,
}


def add(p, q, s=1.0):
    return tuple(a + s * b for a, b in zip(p, q))


def scale(v, s):
    return tuple(s * c for c in v)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def unit(v):
    n = math.sqrt(dot(v, v))
    return tuple(c / n for c in v)


AXIS = cross(SX, SY)


def frame(x, y, z):
    """The point C + x SX + y SY + z AXIS."""
    return add(add(add(CENTER, SX, x), SY, y), AXIS, z)


def framed(v):
    return add(add(scale(SX, v[0]), SY, v[1]), AXIS, v[2])


# --- Surfaces ----------------------------------------------------------------


def bernstein(n, t):
    return [math.comb(n, i) * t**i * (1.0 - t) ** (n - i) for i in range(n + 1)]


def bernstein_derivative(n, t):
    lower = bernstein(n - 1, t)
    return [n * ((lower[i - 1] if i > 0 else 0.0) - (lower[i] if i < n else 0.0)) for i in range(n + 1)]


def patch_poles(rational):
    """A bicubic (or, rational, biquadratic) net over [0, 1]^2, row by row
    with u along the rows, wavy in the surface frame's axis."""
    if rational:
        poles, weights = [], []
        for i in range(3):
            for j in range(3):
                z = 0.3 * math.sin(1.1 * i + 0.6 * j) + (0.5 if (i, j) == (1, 1) else 0.0)
                poles.append(frame(1.5 * i - 1.5, 1.5 * j - 1.5, z))
                weights.append(0.75 if 1 in (i, j) and (i, j) != (1, 1) else (0.6 if (i, j) == (1, 1) else 1.0))
        return 2, poles, weights
    poles = []
    for i in range(4):
        for j in range(4):
            z = 0.4 * math.sin(1.3 * i + 0.7 * j) - 0.1 * i * j
            poles.append(frame(i - 1.5, j - 1.5, z))
    return 3, poles, None


def patch_eval(rational, u, v):
    """The patch's point and its partials in u and v."""
    degree, poles, weights = patch_poles(rational)
    bu, bv = bernstein(degree, u), bernstein(degree, v)
    du, dv = bernstein_derivative(degree, u), bernstein_derivative(degree, v)
    n = degree + 1
    h = [0.0] * 4
    hu = [0.0] * 4
    hv = [0.0] * 4
    for i in range(n):
        for j in range(n):
            w = weights[i * n + j] if weights else 1.0
            p = poles[i * n + j]
            homogeneous = (w * p[0], w * p[1], w * p[2], w)
            for k in range(4):
                h[k] += bu[i] * bv[j] * homogeneous[k]
                hu[k] += du[i] * bv[j] * homogeneous[k]
                hv[k] += bu[i] * dv[j] * homogeneous[k]
    point = tuple(h[k] / h[3] for k in range(3))
    su = tuple((hu[k] - point[k] * hu[3]) / h[3] for k in range(3))
    sv = tuple((hv[k] - point[k] * hv[3]) / h[3] for k in range(3))
    return point, su, sv


def surface_point(kind):
    """A point of the surface, two independent tangents there and the unit
    normal, from the defining formula, and the surface as the case states it."""
    u, v = 0.9, 0.6
    if kind == "plane":
        point, su, sv = frame(0.7, -0.4, 0.0), SX, SY
        data = {"origin": CENTER}
    elif kind == "cylinder":
        r = 1.5
        point = frame(r * math.cos(u), r * math.sin(u), v)
        su, sv = framed((-r * math.sin(u), r * math.cos(u), 0.0)), AXIS
        data = {"origin": CENTER, "radius": r}
    elif kind == "cone":
        r, a = 1.0, 0.5
        rho = r + v * math.sin(a)
        point = frame(rho * math.cos(u), rho * math.sin(u), v * math.cos(a))
        su = framed((-rho * math.sin(u), rho * math.cos(u), 0.0))
        sv = framed((math.sin(a) * math.cos(u), math.sin(a) * math.sin(u), math.cos(a)))
        data = {"origin": CENTER, "radius": r, "semi_angle": a}
    elif kind == "sphere":
        r = 2.0
        point = frame(r * math.cos(v) * math.cos(u), r * math.cos(v) * math.sin(u), r * math.sin(v))
        su = framed((-math.cos(v) * math.sin(u), math.cos(v) * math.cos(u), 0.0))
        sv = framed((-math.sin(v) * math.cos(u), -math.sin(v) * math.sin(u), math.cos(v)))
        data = {"center": CENTER, "radius": r}
    elif kind == "torus":
        big, small = 3.0, 1.0
        rho = big + small * math.cos(v)
        point = frame(rho * math.cos(u), rho * math.sin(u), small * math.sin(v))
        su = framed((-math.sin(u), math.cos(u), 0.0))
        sv = framed((-math.sin(v) * math.cos(u), -math.sin(v) * math.sin(u), math.cos(v)))
        data = {"center": CENTER, "major": big, "minor": small}
    elif kind == "ellipsoid":
        a, b, c = 2.0, 1.5, 1.0
        point = frame(a * math.cos(v) * math.cos(u), b * math.cos(v) * math.sin(u), c * math.sin(v))
        su = framed((-a * math.cos(v) * math.sin(u), b * math.cos(v) * math.cos(u), 0.0))
        sv = framed((-a * math.sin(v) * math.cos(u), -b * math.sin(v) * math.sin(u), c * math.cos(v)))
        data = {"center": CENTER, "a": a, "b": b, "c": c}
    else:
        rational = kind == "nurbs"
        point, su, sv = patch_eval(rational, 0.4, 0.55)
        degree, poles, weights = patch_poles(rational)
        knots = [0.0] * (degree + 1) + [1.0] * (degree + 1)
        data = {
            "degree_u": degree,
            "degree_v": degree,
            "knots_u": knots,
            "knots_v": knots,
            "poles": [list(p) for p in poles],
        }
        if weights:
            data["weights"] = list(weights)
        surface = {"name": "S", "type": "bspline", **data}
        return point, su, sv, unit(cross(su, sv)), surface
    surface = {"name": "S", "type": kind, **data, "x_axis": SX, "y_axis": SY}
    return point, su, sv, unit(cross(su, sv)), surface


# --- Curves ------------------------------------------------------------------


def bezier(poles, weights, t):
    """A rational Bezier's point and derivative, in its plane."""
    n = len(poles) - 1
    b = bernstein(n, t)
    db = bernstein_derivative(n, t)
    w = weights or (1.0,) * len(poles)
    hx = sum(bi * wi * p[0] for bi, wi, p in zip(b, w, poles))
    hy = sum(bi * wi * p[1] for bi, wi, p in zip(b, w, poles))
    hw = sum(bi * wi for bi, wi in zip(b, w))
    dx = sum(bi * wi * p[0] for bi, wi, p in zip(db, w, poles))
    dy = sum(bi * wi * p[1] for bi, wi, p in zip(db, w, poles))
    dw = sum(bi * wi for bi, wi in zip(db, w))
    x, y = hx / hw, hy / hw
    return (x, y), ((dx - x * dw) / hw, (dy - y * dw) / hw)


def local(kind, t):
    """A curve's point and tangent in its own plane."""
    if kind == "line":
        return (t, 0.0), (1.0, 0.0)
    if kind == "circle":
        return (0.8 * math.cos(t), 0.8 * math.sin(t)), (-0.8 * math.sin(t), 0.8 * math.cos(t))
    if kind == "ellipse":
        return (1.2 * math.cos(t), 0.6 * math.sin(t)), (-1.2 * math.sin(t), 0.6 * math.cos(t))
    if kind == "parabola":
        return (t * t, t), (2.0 * t, 1.0)
    if kind == "hyperbola":
        return (0.8 * math.cosh(t), 0.5 * math.sinh(t)), (0.8 * math.sinh(t), 0.5 * math.cosh(t))
    if kind == "bspline":
        return bezier(BEZIER, None, t)
    return bezier(QUARTER, QUARTER_WEIGHTS, t)


def placement(kind, t, point, direction):
    """The rigid motion taking the curve's point at t to `point` and its
    tangent there along `direction`: the plane's origin and axes."""
    p, d = local(kind, t)
    tangent = unit(d)
    e1 = unit(direction)
    e2 = unit(add(TILT, e1, -dot(TILT, e1)))
    x = add(scale(e1, tangent[0]), e2, -tangent[1])
    y = add(scale(e1, tangent[1]), e2, tangent[0])
    origin = add(add(point, x, -p[0]), y, -p[1])
    return origin, x, y


def round_all(v):
    # Rounded to a grid far below the claims' reach: the last digits of a
    # platform's cos and cosh would otherwise change the files between
    # machines.
    return [round(c, 12) for c in v]


def describe(kind, origin, x, y):
    """The curve as the case file states it."""
    curve = {"name": "C", "type": kind if kind not in ("bspline", "nurbs") else "bspline"}
    x, y = round_all(x), round_all(y)
    if kind == "line":
        curve.update(origin=origin, direction=x, domain=[-2.0, 2.0])
    elif kind == "circle":
        curve.update(center=origin, radius=0.8, x_axis=x, y_axis=y, domain=[0.0, 2.0 * math.pi])
    elif kind == "ellipse":
        curve.update(center=origin, major=1.2, minor=0.6, x_axis=x, y_axis=y, domain=[0.0, 2.0 * math.pi])
    elif kind == "parabola":
        curve.update(vertex=origin, focal=0.25, x_axis=x, y_axis=y, domain=[-1.5, 1.5])
    elif kind == "hyperbola":
        curve.update(center=origin, major=0.8, minor=0.5, x_axis=x, y_axis=y, domain=[-1.5, 1.5])
    else:
        local_poles = BEZIER if kind == "bspline" else QUARTER
        poles = [round_all(add(add(origin, x, p[0]), y, p[1])) for p in local_poles]
        degree = len(local_poles) - 1
        curve.update(degree=degree, knots=[0.0] * (degree + 1) + [1.0] * (degree + 1), poles=poles, domain=[0.0, 1.0])
        if kind == "nurbs":
            curve["weights"] = list(QUARTER_WEIGHTS)
    for key in ("origin", "center", "vertex"):
        if key in curve:
            curve[key] = round_all(curve[key])
    return curve


def case_for(curve_kind, surface_kind, angle_name, degrees):
    point, su, _, normal, surface = surface_point(surface_kind)
    along = unit(su)
    across = cross(normal, along)
    e = add(scale(along, math.cos(TURN)), across, math.sin(TURN))
    b = math.radians(degrees)
    direction = add(scale(normal, math.sin(b)), e, math.cos(b))
    t = CURVE_PARAMETER[curve_kind]
    origin, x, y = placement(curve_kind, t, point, direction)
    for key in ("origin", "center", "x_axis", "y_axis"):
        if key in surface:
            surface[key] = list(surface[key])
    case_id = f"curve-surface/{curve_kind}-{surface_kind}/{angle_name}"
    expected = {"point": round_all(point), "parameter": t, "angle_degrees": degrees}
    limit = KNOWN_LIMITS.get((curve_kind, surface_kind, angle_name))
    if limit:
        expected["may_fail_closed"] = True
        expected["known_limit"] = limit
    return {
        "format": "stereon-corpus-case/1",
        "id": case_id,
        "kind": "reference",
        "summary": f"A {curve_kind} crossing a {surface_kind} at a known point at {degrees:g} degrees to its tangent plane.",
        "added": "2026-11-22",
        "source": "Generated by tools/generate_curve_surface_corpus.py",
        "licence": "MPL-2.0",
        "references": ["ADR-0029", "ADR-0053"],
        "tags": ["generated", angle_name],
        "inputs": {"curve": describe(curve_kind, origin, x, y), "surface": surface},
        "operation": {"name": "intersect", "operands": ["C", "S"]},
        "expected": expected,
    }


# Cases the kernel refuses today (IntersectionFailed), never answering
# wrongly, as measured on 2026-11-22 (main at 898ed69). A case that starts
# passing keeps its mark until the generator drops it.
SPLINE_SURFACE_GRAZING = (
    "Krawczyk's test on C(s) - S(u, v) does not settle a crossing at 1 degree within the "
    "subdivision budget of the B-spline surface search (CurveSplineSurface.cpp)"
)
KNOWN_LIMITS = {
    **{(c, s, "grazing"): SPLINE_SURFACE_GRAZING for c in CURVES for s in ("bspline", "nurbs")},
    ("line", "cone", "grazing"): (
        "the line-cone handler leaves a root's existence or multiplicity unresolved at 1 degree "
        "(LineSurface.cpp)"
    ),
}


def all_cases():
    for c in CURVES:
        for s in SURFACES:
            for angle_name, degrees in ANGLES.items():
                yield case_for(c, s, angle_name, degrees)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument("--check", action="store_true", help="fail if a file differs")
    args = parser.parse_args()
    stale = []
    count = 0
    for case in all_cases():
        count += 1
        path = args.root / "corpus" / (case["id"] + ".json")
        text = json.dumps(case, indent=2) + "\n"
        if args.check:
            if not path.exists() or path.read_text(encoding="utf-8") != text:
                stale.append(case["id"])
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(text, encoding="utf-8", newline="\n")
    if stale:
        print(f"generate_curve_surface_corpus: {len(stale)} stale case(s), first {stale[0]}", file=sys.stderr)
        return 1
    print(f"generate_curve_surface_corpus: {count} cases {'match' if args.check else 'written'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
