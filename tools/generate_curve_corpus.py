#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""Generate corpus/curve-curve/: pairs of 3D curves with a known closest approach.

For every unordered pair of curve types, the second curve is translated so
that its point at a chosen parameter lies a chosen gap from the first
curve's point at its chosen parameter, along their common normal. That
makes the two parameters a critical point of the distance, at exactly that
gap, whichever curves they are. The gaps straddle the resolution
(ADR-0001), so a case states whether the curves meet there.

Usage: generate_curve_corpus.py [--root <repository root>] [--check]
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

RESOLUTION = 1e-5
GAPS = {"touching": 0.0, "within": 0.4, "beyond": 1.5, "apart": 10.0}
TYPES = ("line", "circle", "ellipse", "parabola", "hyperbola", "bspline")
PARAMETERS = {
    "line": (1.1, -0.7),
    "circle": (0.7, 2.5),
    "ellipse": (1.1, 4.0),
    "parabola": (0.6, -1.3),
    "hyperbola": (0.5, -1.0),
    "bspline": (0.3, 0.65),
}
BSPLINE_LOCAL = ((-2.0, 0.0), (-1.0, 1.5), (1.0, -1.5), (2.0, 0.0))
FRAMES = (
    ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
    ((0.0, 0.6, 0.8), (1.0, 0.0, 0.0)),
)
ORIGIN = (1.0, 2.0, 3.0)


def add(p, q, s=1.0):
    return tuple(a + s * b for a, b in zip(p, q))


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def length(v):
    return math.sqrt(sum(c * c for c in v))


def on_frame(frame, x, y):
    return add(tuple(x * c for c in frame[0]), tuple(y * c for c in frame[1]))


def bezier(poles, t):
    work = list(poles)
    while len(work) > 1:
        work = [(a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])) for a, b in zip(work, work[1:])]
    return work[0]


def local(kind, t):
    """A curve's point and tangent in its own frame's plane."""
    if kind == "line":
        return (t, 0.0), (1.0, 0.0)
    if kind == "circle":
        return (2.0 * math.cos(t), 2.0 * math.sin(t)), (-2.0 * math.sin(t), 2.0 * math.cos(t))
    if kind == "ellipse":
        return (3.0 * math.cos(t), 1.5 * math.sin(t)), (-3.0 * math.sin(t), 1.5 * math.cos(t))
    if kind == "parabola":
        return (t * t / 2.0, t), (t, 1.0)
    if kind == "hyperbola":
        return (2.0 * math.cosh(t), math.sinh(t)), (2.0 * math.sinh(t), math.cosh(t))
    point = bezier(BSPLINE_LOCAL, t)
    hodograph = [(3.0 * (b[0] - a[0]), 3.0 * (b[1] - a[1])) for a, b in zip(BSPLINE_LOCAL, BSPLINE_LOCAL[1:])]
    return point, bezier(hodograph, t)


def description(kind, name, origin, frame):
    """The curve as the case file states it."""
    x, y = frame
    curve = {"name": name, "type": kind}
    if kind == "line":
        curve.update(origin=origin, direction=x, domain=[-4.0, 4.0])
    elif kind == "circle":
        curve.update(center=origin, radius=2.0, x_axis=x, y_axis=y, domain=[0.0, 2.0 * math.pi])
    elif kind == "ellipse":
        curve.update(center=origin, major=3.0, minor=1.5, x_axis=x, y_axis=y, domain=[0.0, 2.0 * math.pi])
    elif kind == "parabola":
        curve.update(vertex=origin, focal=0.5, x_axis=x, y_axis=y, domain=[-3.0, 3.0])
    elif kind == "hyperbola":
        curve.update(center=origin, major=2.0, minor=1.0, x_axis=x, y_axis=y, domain=[-2.0, 2.0])
    else:
        poles = [add(origin, on_frame(frame, p[0], p[1])) for p in BSPLINE_LOCAL]
        curve.update(degree=3, knots=[0.0] * 4 + [1.0] * 4, poles=[list(p) for p in poles], domain=[0.0, 1.0])
    for key, value in curve.items():
        if isinstance(value, tuple):
            curve[key] = list(value)
    return curve


def case_for(a, b, placement, gap_name, ratio):
    ta = PARAMETERS[a][placement]
    tb = PARAMETERS[b][placement]
    pa, da = local(a, ta)
    pb, db = local(b, tb)
    point_a = add(ORIGIN, on_frame(FRAMES[0], *pa))
    tangent_a = on_frame(FRAMES[0], *da)
    tangent_b = on_frame(FRAMES[1], *db)
    normal = cross(tangent_a, tangent_b)
    sine = length(normal) / (length(tangent_a) * length(tangent_b))
    if sine < 0.2:
        raise SystemExit(f"{a}-{b} placement {placement}: tangents too close to parallel ({sine})")
    unit = tuple(c / length(normal) for c in normal)
    gap = ratio * RESOLUTION
    target = add(point_a, unit, gap)
    # Rounded to a grid far below the gaps: the last digits of a platform's
    # cos and cosh would otherwise change the files between machines.
    origin_b = tuple(round(c, 12) for c in add(target, on_frame(FRAMES[1], *pb), -1.0))
    case_id = f"curve-curve/{a}-{b}/{gap_name}-{placement}"
    return {
        "format": "stereon-corpus-case/1",
        "id": case_id,
        "kind": "reference",
        "summary": f"A {a} and a {b} whose closest approach near the stated parameters is {ratio} resolutions.",
        "added": "2026-11-09",
        "source": "Generated by tools/generate_curve_corpus.py",
        "licence": "MPL-2.0",
        "references": ["ADR-0001", "ADR-0029"],
        "tags": ["generated", gap_name],
        "inputs": {
            "curves": [
                description(a, "A", ORIGIN, FRAMES[0]),
                description(b, "B", origin_b, FRAMES[1]),
            ]
        },
        "operation": {"name": "intersect", "operands": ["A", "B"]},
        "expected": {"meets": ratio <= 1.0, "parameters": [ta, tb], "gap": gap},
    }


def all_cases():
    for i, a in enumerate(TYPES):
        for b in TYPES[i:]:
            for placement in range(2):
                for gap_name, ratio in GAPS.items():
                    yield case_for(a, b, placement, gap_name, ratio)


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
        print(f"generate_curve_corpus: {len(stale)} stale case(s), first {stale[0]}", file=sys.stderr)
        return 1
    print(f"generate_curve_corpus: {count} cases {'match' if args.check else 'written'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
