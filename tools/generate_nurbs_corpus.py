#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""Generate corpus/nurbs-ssi/: pairs of B-spline surfaces with a known answer.

Every surface is a B-spline: the exact NURBS form of a plane, a cylinder or
a sphere (rational quadratic arcs), a bilinear saddle that is exactly
z = x y, or a bicubic bump that is exactly z = 9 H u (1 - u) v (1 - v). The
expected topology follows from the analytic surface each one is, and the
witnesses from closed forms; the expected result never comes from the
kernel. The configurations cover what ADR-0010's NURBS corpus asks for:
transverse crossings, separation and near-coincidence a few tolerances
either side, small loops near a maximum, and contact at a sphere's pole.
Each is placed three ways and at two scales, by moving the poles, which a
NURBS surface follows exactly.

Usage: generate_nurbs_corpus.py [--root <repository root>] [--check]
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

EX, EY, EZ = (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)
PLACEMENTS = {
    "origin": (EX, EY, EZ, (0.0, 0.0, 0.0)),
    # A rational rotation, so that the frame is orthonormal to rounding.
    "tilted": ((0.6, 0.8, 0.0), (-0.48, 0.36, 0.8), (0.64, -0.48, 0.6), (12.0, -7.0, 5.0)),
    # Near the edge of the size box of ADR-0001, as the other corpora.
    "far": (EX, EY, EZ, (480000.0, -300000.0, 210000.0)),
}
SCALES = (1.0, 100.0)
TOLERANCE = 1e-5
HALF_ROOT = math.sqrt(0.5)
# The rational quadratic full circle: nine poles, four arcs.
CIRCLE_KNOTS = [0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0]
CIRCLE_POINTS = [(1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1), (1, 0)]
CIRCLE_WEIGHTS = [1.0, HALF_ROOT, 1.0, HALF_ROOT, 1.0, HALF_ROOT, 1.0, HALF_ROOT, 1.0]
# The rational quadratic half circle from the south pole to the north one.
HALF_KNOTS = [0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0]
HALF_POINTS = [(0, -1), (1, -1), (1, 0), (1, 1), (0, 1)]
HALF_WEIGHTS = [1.0, HALF_ROOT, 1.0, HALF_ROOT, 1.0]
LINEAR = [0.0, 0.0, 1.0, 1.0]
CUBIC = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]

TANGENT_REASON = "a tangency along a curve or at a point of rational surfaces; the generic path fails closed on it"
NODE_REASON = ("two branches cross at a node of rational surfaces; the generic path fails closed at "
               "a branch point it cannot prove a transverse node of polynomial surfaces")
COINCIDENT_REASON = "the surfaces coincide within the tolerance; the generic path fails closed on an overlap"
TINY_REASON = "a loop a few thousandths of the patch across; subdivision may exhaust its budget first"


def surface(degree_u, degree_v, knots_u, knots_v, poles, weights=None):
    """A B-spline surface: poles row by row, i along u."""
    s = {"type": "bspline", "degree_u": degree_u, "degree_v": degree_v,
         "knots_u": knots_u, "knots_v": knots_v, "poles": poles}
    if weights:
        s["weights"] = weights
    return s


def plane_patch(origin, x, y, lx, ly):
    """The rectangle origin + u lx x + v ly y over the unit square, bilinear."""
    poles = []
    for i in range(2):
        for j in range(2):
            poles.append(tuple(origin[k] + i * lx * x[k] + j * ly * y[k] for k in range(3)))
    return surface(1, 1, LINEAR, LINEAR, poles)


def cylinder(center, radius, half_height):
    """The cylinder of the given radius about the z axis through `center`,
    |z - center z| <= half_height: a circle in u, a line in v."""
    poles, weights = [], []
    for (cx, cy), w in zip(CIRCLE_POINTS, CIRCLE_WEIGHTS):
        for h in (-half_height, half_height):
            poles.append((center[0] + radius * cx, center[1] + radius * cy, center[2] + h))
            weights.append(w)
    return surface(2, 1, CIRCLE_KNOTS, LINEAR, poles, weights)


def cylinder_x(center, radius, half_length):
    """The same about the x axis."""
    poles, weights = [], []
    for (cy, cz), w in zip(CIRCLE_POINTS, CIRCLE_WEIGHTS):
        for h in (-half_length, half_length):
            poles.append((center[0] + h, center[1] + radius * cy, center[2] + radius * cz))
            weights.append(w)
    return surface(2, 1, CIRCLE_KNOTS, LINEAR, poles, weights)


def sphere(center, radius):
    """The sphere: a full circle in u, a half circle from pole to pole in v.
    The first and last rows of v collapse to the poles."""
    poles, weights = [], []
    for (cx, cy), wu in zip(CIRCLE_POINTS, CIRCLE_WEIGHTS):
        for (rho, z), wv in zip(HALF_POINTS, HALF_WEIGHTS):
            poles.append((center[0] + radius * rho * cx, center[1] + radius * rho * cy,
                          center[2] + radius * z))
            weights.append(wu * wv)
    return surface(2, 2, CIRCLE_KNOTS, HALF_KNOTS, poles, weights)


def saddle():
    """The bilinear patch through (+-1, +-1, +-1): exactly z = x y over
    [-1, 1]^2, with x = 2u - 1 and y = 2v - 1."""
    poles = [(-1.0, -1.0, 1.0), (-1.0, 1.0, -1.0), (1.0, -1.0, -1.0), (1.0, 1.0, 1.0)]
    return surface(1, 1, LINEAR, LINEAR, poles)


def bump(height):
    """The bicubic over [0, 3]^2 whose four inner poles stand at `height`:
    x = 3u, y = 3v (evenly spaced poles) and z = 9 H u (1 - u) v (1 - v),
    since B1 + B2 = 3 u (1 - u). Its top, 9 H / 16, is at the centre."""
    poles = []
    for i in range(4):
        for j in range(4):
            inner = i in (1, 2) and j in (1, 2)
            poles.append((float(i), float(j), height if inner else 0.0))
    return surface(3, 3, CUBIC, CUBIC, poles)


def regime(name, a, b, curves=0, points=0, contact="crossing", witnesses=(), overlap=False, hard=None):
    return {"name": name, "a": a, "b": b, "curves": curves, "points": points, "contact": contact,
            "witnesses": list(witnesses), "overlap": overlap, "hard": hard}


def plane_plane():
    flat = plane_patch((-2.0, -2.0, 0.0), EX, EY, 4.0, 4.0)
    upright = plane_patch((0.5, -2.0, -2.0), EY, EZ, 4.0, 4.0)  # x = 0.5
    return [
        regime("crossing", flat, upright, curves=1, witnesses=[(0.5, 0.0, 0.0), (0.5, 1.5, 0.0)]),
        regime("apart-10tol", flat, plane_patch((-2.0, -2.0, 10 * TOLERANCE), EX, EY, 4.0, 4.0)),
        regime("coincident-half-tol", flat, plane_patch((-2.0, -2.0, 0.5 * TOLERANCE), EX, EY, 4.0, 4.0),
               overlap=True, hard=COINCIDENT_REASON),
    ]


def saddle_plane():
    flat = lambda z: plane_patch((-2.0, -2.0, z), EX, EY, 4.0, 4.0)  # noqa: E731
    c = 0.5
    return [
        # x y = 1/2 inside [-1, 1]^2: two branches, x and y of one sign.
        regime("hyperbola", saddle(), flat(c), curves=2,
               witnesses=[(1.0, c, c), (c, 1.0, c), (-1.0, -c, c), (-c, -1.0, c), (0.75, c / 0.75, c)]),
        # x y = 0: two lines crossing at a transverse node of polynomial
        # surfaces, which the generic path splits there.
        regime("node", saddle(), flat(0.0), curves=1, witnesses=[(0.5, 0.0, 0.0), (0.0, -0.5, 0.0)]),
    ]


def bump_plane():
    height = 1.0
    top = 9.0 * height / 16.0
    out = []
    for name, depth, hard in (("loop-1e-3", 1e-3, None), ("loop-1e-5", 1e-5, TINY_REASON)):
        level = top - depth
        # Along v = 1/2: (9 H / 4) u (1 - u) = level, so u (1 - u) = q.
        q = level / (9.0 * height / 4.0)
        u = (1.0 - math.sqrt(1.0 - 4.0 * q)) / 2.0
        flat = plane_patch((-1.0, -1.0, level), EX, EY, 5.0, 5.0)
        out.append(regime(name, bump(height), flat, curves=1,
                          witnesses=[(3.0 * u, 1.5, level), (3.0 * (1.0 - u), 1.5, level)], hard=hard))
    return out


def cylinder_plane():
    tube = cylinder((0.0, 0.0, 0.0), 1.0, 3.0)
    return [
        regime("circle", tube, plane_patch((-2.0, -2.0, 0.5), EX, EY, 4.0, 4.0), curves=1,
               witnesses=[(1.0, 0.0, 0.5), (0.0, -1.0, 0.5)]),
        # The plane x = 1/2 along the axis: the lines x = 1/2, y = +-sqrt(3)/2.
        regime("two-lines", tube, plane_patch((0.5, -2.0, -4.0), EY, EZ, 4.0, 8.0), curves=2,
               witnesses=[(0.5, math.sqrt(0.75), 1.0), (0.5, -math.sqrt(0.75), -1.0)]),
        regime("tangent-line", tube, plane_patch((1.0, -2.0, -4.0), EY, EZ, 4.0, 8.0), curves=1,
               contact="tangential", witnesses=[(1.0, 0.0, 0.0)], hard=TANGENT_REASON),
    ]


def sphere_plane():
    ball = sphere((0.0, 0.0, 0.0), 1.0)
    flat = lambda z: plane_patch((-2.0, -2.0, z), EX, EY, 4.0, 4.0)  # noqa: E731
    near = 1.0 - 1e-3
    rho = math.sqrt(1.0 - near * near)
    return [
        regime("circle", ball, flat(0.25), curves=1, witnesses=[(math.sqrt(1.0 - 0.0625), 0.0, 0.25)]),
        regime("near-pole", ball, flat(near), curves=1, witnesses=[(rho, 0.0, near), (0.0, rho, near)]),
        regime("pole-contact", ball, flat(1.0), points=1, contact="tangential", witnesses=[(0.0, 0.0, 1.0)],
               hard=TANGENT_REASON),
    ]


def cylinder_cylinder():
    small = cylinder((0.0, 0.0, 0.0), 1.0, 3.0)
    return [
        # Radii 1 about z and 2 about x: the smaller passes through the
        # larger, two loops, through (0, +-1, +-sqrt 3).
        regime("two-loops", small, cylinder_x((0.0, 0.0, 0.0), 2.0, 3.0), curves=2,
               witnesses=[(0.0, 1.0, math.sqrt(3.0)), (0.0, 1.0, -math.sqrt(3.0))]),
        regime("node", small, cylinder_x((0.0, 0.0, 0.0), 1.0, 3.0), curves=1,
               witnesses=[(0.0, 1.0, 0.0), (1.0, 0.0, 1.0)], hard=NODE_REASON),
    ]


FAMILIES = {
    "plane-plane": plane_plane,
    "saddle-plane": saddle_plane,
    "bump-plane": bump_plane,
    "cylinder-plane": cylinder_plane,
    "sphere-plane": sphere_plane,
    "cylinder-cylinder": cylinder_cylinder,
}


def rotated(placement, v):
    x, y, z, _ = placement
    return tuple(v[0] * x[i] + v[1] * y[i] + v[2] * z[i] for i in range(3))


def placed_point(placement, scale, p):
    # Rounded far below the tolerance, so that files do not depend on the
    # platform's last digits.
    q = rotated(placement, tuple(scale * c for c in p))
    return [round(q[i] + placement[3][i], 10) for i in range(3)]


def placed_surface(s, name, placement, scale):
    out = {"name": name}
    for key, value in s.items():
        if key == "poles":
            out[key] = [placed_point(placement, scale, p) for p in value]
        elif key == "weights":
            out[key] = [round(w, 15) for w in value]
        else:
            out[key] = value
    return out


def case_for(family, r, placement_name, scale):
    placement = PLACEMENTS[placement_name]
    expected = {
        "overlap": r["overlap"],
        "curves": r["curves"],
        "points": r["points"],
        "contact": r["contact"],
        "witnesses": [placed_point(placement, scale, w) for w in r["witnesses"]],
    }
    if r["hard"]:
        expected["may_fail_closed"] = r["hard"]
    a, b = family.split("-")
    return {
        "format": "stereon-corpus-case/1",
        "id": f"nurbs-ssi/{family}/{r['name']}-{placement_name}-x{int(scale)}",
        "kind": "hard" if r["hard"] else "reference",
        "summary": f"A B-spline {a} and a B-spline {b}, {r['name']}, placed at {placement_name}, "
                   f"scaled by {int(scale)}.",
        "added": "2026-11-10",
        "source": "Generated by tools/generate_nurbs_corpus.py",
        "licence": "MPL-2.0",
        "references": ["ADR-0010", "ADR-0029"],
        "tags": ["generated", "nurbs", r["name"]],
        "inputs": {
            "surfaces": [placed_surface(r["a"], "A", placement, scale),
                         placed_surface(r["b"], "B", placement, scale)],
            "tolerance": TOLERANCE * scale,
            "scale": scale,
        },
        "operation": {"name": "intersect", "operands": ["A", "B"]},
        "expected": expected,
    }


def all_cases():
    for family, regimes in FAMILIES.items():
        for r in regimes():
            for placement_name in PLACEMENTS:
                for scale in SCALES:
                    yield case_for(family, r, placement_name, scale)


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
    suite = args.root / "corpus" / "nurbs-ssi"
    if suite.exists():
        for path in sorted(suite.glob("*/*.json")):
            if path.resolve() not in expected_paths:
                if args.check:
                    stale.append(str(path.relative_to(args.root)))
                else:
                    path.unlink()
    if stale:
        print(f"generate_nurbs_corpus: {len(stale)} stale case(s), first {stale[0]}", file=sys.stderr)
        return 1
    print(f"generate_nurbs_corpus: {count} cases {'match' if args.check else 'written'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
