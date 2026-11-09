#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""Generate corpus/curved-ssi/: pairs of analytic surfaces with a known answer.

Every unordered pair of analytic surfaces with a cone, a torus or an ellipsoid
in it, in configurations whose topology follows from a closed-form argument:
coaxial profiles meeting, a section of known type, a surface whose distance
from an axis or a centre bounds where it can meet the other. The expected
result never comes from the kernel. Each configuration is placed three ways
and at two scales.

Usage: generate_surface_corpus.py [--root <repository root>] [--check]
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

EX, EY, EZ = (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)
TAU = 2.0 * math.pi
PLACEMENTS = {
    "origin": (EX, EY, EZ, (0.0, 0.0, 0.0)),
    # A rational rotation, so that the frame is orthonormal to rounding.
    "tilted": ((0.6, 0.8, 0.0), (-0.48, 0.36, 0.8), (0.64, -0.48, 0.6), (12.0, -7.0, 5.0)),
    # Near the edge of the size box of ADR-0001, as the analytic-ssi sweeps.
    "far": (EX, EY, EZ, (480000.0, -300000.0, 210000.0)),
}
SCALES = (1.0, 100.0)
TOLERANCE = 1e-5
TANGENT_REASON = "a tangency or node away from the closed forms; the trace may fail closed"
COINCIDENT_REASON = "coincidence over part of a chart, or at a pole, is not resolved yet"
POLE_REASON = "the section passes through a pole of the chart, where its pcurve has no longitude"
BUDGET_REASON = "a traced branch of the larger scale may exceed the sample budget before it certifies"
EXACT_ONLY = ("origin", "far")


def plane(origin, x=EX, y=EY):
    return {"type": "plane", "origin": origin, "x_axis": x, "y_axis": y}


def plane_x(d):
    """The plane x = d, parallel to the z axis."""
    return plane((d, 0.0, 0.0), EY, EZ)


def plane_y(d):
    return plane((0.0, d, 0.0), EZ, EX)


def sphere(center, radius):
    return {"type": "sphere", "center": center, "radius": radius, "x_axis": EX, "y_axis": EY}


def cylinder(origin, radius, v=(-8.0, 8.0)):
    """A cylinder about the z axis direction, bounded along it."""
    return {"type": "cylinder", "origin": origin, "radius": radius, "x_axis": EX, "y_axis": EY, "v": v}


def cone(origin, radius, angle, v=(-8.0, 8.0)):
    return {"type": "cone", "origin": origin, "radius": radius, "semi_angle": angle, "x_axis": EX, "y_axis": EY, "v": v}


def torus(center, major, minor):
    return {"type": "torus", "center": center, "major": major, "minor": minor, "x_axis": EX, "y_axis": EY}


def ellipsoid(center, a, b, c):
    return {"type": "ellipsoid", "center": center, "a": a, "b": b, "c": c, "x_axis": EX, "y_axis": EY}


O = (0.0, 0.0, 0.0)


def z_at(h):
    return (0.0, 0.0, h)


def x_at(d):
    return (d, 0.0, 0.0)


def regime(
    name,
    a,
    b,
    curves=0,
    points=0,
    contact="crossing",
    witnesses=(),
    overlap=False,
    hard=None,
    placements=tuple(PLACEMENTS),
    hard_large=None,
):
    """A configuration and its answer. `hard` allows a refusal everywhere,
    `hard_large` only at the larger scale."""
    return {
        "placements": placements,
        "hard_large": hard_large,
        "name": name,
        "a": a,
        "b": b,
        "curves": curves,
        "points": points,
        "contact": contact,
        "witnesses": list(witnesses),
        "overlap": overlap,
        "hard": hard,
    }


def on_ellipsoid_line(p, d, axes):
    """The two points of the line p + s d on the origin-centred ellipsoid."""
    qa = sum((d[i] / axes[i]) ** 2 for i in range(3))
    qb = 2.0 * sum(p[i] * d[i] / axes[i] ** 2 for i in range(3))
    qc = sum((p[i] / axes[i]) ** 2 for i in range(3)) - 1.0
    root = math.sqrt(qb * qb - 4.0 * qa * qc)
    return [tuple(p[i] + s * d[i] for i in range(3)) for s in ((-qb - root) / (2 * qa), (-qb + root) / (2 * qa))]


def circle_circle(p, r, q, s):
    """The two meets of two circles of a plane."""
    dx, dy = q[0] - p[0], q[1] - p[1]
    d = math.hypot(dx, dy)
    along = (d * d + r * r - s * s) / (2.0 * d)
    across = math.sqrt(r * r - along * along)
    mid = (p[0] + along * dx / d, p[1] + along * dy / d)
    return [(mid[0] - across * dy / d, mid[1] + across * dx / d), (mid[0] + across * dy / d, mid[1] - across * dx / d)]


def plane_cone():
    a = math.pi / 6.0
    t = math.tan(a)
    c = cone(O, 0.0, a)
    yield regime("perpendicular", plane(z_at(2.0)), c, curves=1, witnesses=[(2 * t, 0.0, 2.0), (-2 * t, 0.0, 2.0)])
    yield regime("lower-nappe", plane(z_at(-3.0)), c, curves=1, witnesses=[(3 * t, 0.0, -3.0)])
    for degrees in (20, 40):
        th = math.radians(degrees)
        tilted = plane(z_at(2.0), (math.cos(th), 0.0, -math.sin(th)), EY)
        up, down = 2.0 / (1.0 + t * math.tan(th)), 2.0 / (1.0 - t * math.tan(th))
        yield regime(f"tilted-{degrees}", tilted, c, curves=1, witnesses=[(up * t, 0.0, up), (-down * t, 0.0, down)])
    yield regime("beyond-box", plane(z_at(20.0)), c)
    yield regime("apex", plane(O), c, points=1, witnesses=[O])


def plane_torus():
    T = torus(O, 3.0, 1.0)
    for h in (0.0, 0.5):
        z = math.sqrt(1.0 - h * h)
        yield regime(f"across-{h}".replace(".", "p"), plane(z_at(h)), T, curves=2, witnesses=[(3 + z, 0.0, h), (3 - z, 0.0, h)])
    yield regime("across-top", plane(z_at(1.0)), T, curves=1, contact="tangential", witnesses=[(3.0, 0.0, 1.0)])
    yield regime("across-above", plane(z_at(1.5)), T)
    yield regime("meridian", plane_x(0.0), T, curves=2, witnesses=[(0.0, 4.0, 0.0), (0.0, -2.0, 0.0)])
    y = math.sqrt(9.0 - 1.0)
    yield regime("two-ovals", plane_x(1.0), T, curves=2, witnesses=[(1.0, y, 1.0), (1.0, -y, 1.0)])
    yield regime("node", plane_x(2.0), T, curves=1, hard=TANGENT_REASON, witnesses=[(2.0, 0.0, 0.0)])
    y, z = math.sqrt(9.0 - 6.25), math.sqrt(0.75)
    yield regime("one-oval", plane_x(2.5), T, curves=1, witnesses=[(2.5, y, 1.0), (2.5, 0.0, z), (2.5, 0.0, -z)])
    yield regime("touching", plane_x(4.0), T, points=1, witnesses=[(4.0, 0.0, 0.0)], hard=TANGENT_REASON)
    yield regime("apart", plane_x(4.5), T)


def plane_ellipsoid():
    axes = (4.0, 3.0, 2.0)
    E = ellipsoid(O, *axes)
    for d in (0.0, 2.0, 3.9):
        hard = POLE_REASON if d == 0.0 else None
        yield regime(f"x-{d}".replace(".", "p"), plane_x(d), E, curves=1, witnesses=on_ellipsoid_line((d, 0.0, 0.0), EY, axes), hard=hard)
    # The section is decided exactly on the doubles, so a rotation's rounding
    # would make this tangency a small ellipse or a miss.
    yield regime("x-touching", plane_x(4.0), E, points=1, witnesses=[(4.0, 0.0, 0.0)], hard=TANGENT_REASON, placements=EXACT_ONLY)
    yield regime("x-apart", plane_x(4.5), E)
    yield regime("y-1", plane_y(1.0), E, curves=1, witnesses=on_ellipsoid_line((0.0, 1.0, 0.0), EX, axes))
    yield regime("y-apart", plane_y(3.5), E)
    n = (1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0)
    support = math.sqrt(sum((axes[i] * n[i]) ** 2 for i in range(3)))
    root5 = math.sqrt(5.0)
    x = (2.0 / root5, -1.0 / root5, 0.0)
    y = (n[1] * x[2] - n[2] * x[1], n[2] * x[0] - n[0] * x[2], n[0] * x[1] - n[1] * x[0])
    for d in (0.0, 1.0, 2.5):
        p = tuple(d * axes[i] ** 2 * n[i] / support**2 for i in range(3))
        yield regime(f"oblique-{d}".replace(".", "p"), plane(p, x, y), E, curves=1, witnesses=on_ellipsoid_line(p, x, axes))
    p = tuple(1.1 * support * axes[i] ** 2 * n[i] / support**2 for i in range(3))
    yield regime("oblique-apart", plane(p, x, y), E)


def sphere_cone():
    a = math.pi / 6.0
    s, co = math.sin(a), math.cos(a)
    c = cone(O, 0.0, a)

    def along(px, pz, radius):
        """The two points of the generator (s, c) at `radius` from (px, pz)."""
        mid = px * s + pz * co
        root = math.sqrt(mid * mid - (px * px + pz * pz - radius * radius))
        return [((mid + e * root) * s, 0.0, (mid + e * root) * co) for e in (-1.0, 1.0)]

    yield regime("coaxial-apart", sphere(z_at(4.0), 1.5), c)
    yield regime(
        "coaxial-touching",
        sphere(z_at(4.0), 2.0),
        c,
        curves=1,
        contact="tangential",
        witnesses=[(2 * co * s * 2, 0.0, 4 * co * co)],
        hard=TANGENT_REASON,
    )
    yield regime("coaxial-two", sphere(z_at(4.0), 3.0), c, curves=2, witnesses=along(0.0, 4.0, 3.0))
    # The upper circle is at v = 8.04: past the default box.
    yield regime("coaxial-around-apex", sphere(z_at(4.0), 5.0), cone(O, 0.0, a, (-9.0, 9.0)), curves=2, witnesses=along(0.0, 4.0, 5.0))
    px, pz = 3.0, 2.0
    d = abs(px * co - pz * s)
    yield regime("off-axis-apart", sphere((px, 0.0, pz), 1.0), c)
    foot = px * s + pz * co
    yield regime("off-axis-touching", sphere((px, 0.0, pz), d), c, points=1, witnesses=[(foot * s, 0.0, foot * co)])
    yield regime("off-axis-loop", sphere((px, 0.0, pz), 2.0), c, curves=1, witnesses=along(px, pz, 2.0))


def sphere_torus():
    T = torus(O, 3.0, 1.0)
    yield regime("concentric-inside-hole", sphere(O, 1.5), T)
    yield regime("concentric-inner-circle", sphere(O, 2.0), T, curves=1, contact="tangential", witnesses=[(2.0, 0.0, 0.0)])
    x = 17.0 / 6.0
    z = math.sqrt(9.0 - x * x)
    yield regime("concentric-two", sphere(O, 3.0), T, curves=2, witnesses=[(x, 0.0, z), (x, 0.0, -z)])
    yield regime("concentric-outer-circle", sphere(O, 4.0), T, curves=1, contact="tangential", witnesses=[(4.0, 0.0, 0.0)])
    yield regime("concentric-around", sphere(O, 5.0), T)
    yield regime("in-the-tube", sphere(x_at(3.0), 0.5), T)
    outer = circle_circle((0.0, 0.0), 4.0, (3.0, 0.0), 2.0)
    inner = circle_circle((0.0, 0.0), 2.0, (3.0, 0.0), 2.0)
    yield regime(
        "across-the-tube",
        sphere(x_at(3.0), 2.0),
        T,
        curves=2,
        witnesses=[(p[0], p[1], 0.0) for p in outer + inner],
    )
    yield regime("outside-apart", sphere(x_at(6.0), 1.5), T)
    yield regime("outside-touching", sphere(x_at(6.0), 2.0), T, points=1, witnesses=[(4.0, 0.0, 0.0)])
    near = circle_circle((0.0, 0.0), 4.0, (6.0, 0.0), 3.0)
    yield regime("outside-loop", sphere(x_at(6.0), 3.0), T, curves=1, witnesses=[(p[0], p[1], 0.0) for p in near])


def sphere_ellipsoid():
    axes = (4.0, 3.0, 2.0)
    E = ellipsoid(O, *axes)

    def in_plane(i, j, k, radius):
        """A point with coordinate k zero where |p| = radius on the ellipsoid."""
        square_j = (1.0 - radius * radius / axes[i] ** 2) / (1.0 / axes[j] ** 2 - 1.0 / axes[i] ** 2)
        p = [0.0, 0.0, 0.0]
        p[i], p[j] = math.sqrt(radius * radius - square_j), math.sqrt(square_j)
        return p

    yield regime("inside", sphere(O, 1.5), E)
    yield regime("touching-short-axis", sphere(O, 2.0), E, points=2, witnesses=[z_at(2.0), z_at(-2.0)], hard=TANGENT_REASON)
    p = in_plane(1, 2, 0, 2.5)
    yield regime("around-short-axis", sphere(O, 2.5), E, curves=2, witnesses=[tuple(p), (p[0], p[1], -p[2])])
    yield regime("node", sphere(O, 3.0), E, curves=1, witnesses=[(0.0, 3.0, 0.0)], hard=TANGENT_REASON)
    p = in_plane(1, 0, 2, 3.5)
    yield regime("around-long-axis", sphere(O, 3.5), E, curves=2, witnesses=[tuple(p), (-p[0], p[1], p[2])])
    yield regime("touching-long-axis", sphere(O, 4.0), E, points=2, witnesses=[x_at(4.0), x_at(-4.0)], hard=TANGENT_REASON)
    yield regime("outside", sphere(O, 5.0), E)


def cylinder_cone():
    a = math.pi / 6.0
    t = math.tan(a)
    both = cone(O, 0.0, a)
    upper = cone(O, 0.0, a, (0.1, 8.0))
    yield regime("coaxial-both-nappes", cylinder(O, 1.5, (-9.0, 9.0)), both, curves=2, witnesses=[(1.5, 0.0, 1.5 / t), (1.5, 0.0, -1.5 / t)])
    yield regime("coaxial-one-nappe", cylinder(O, 1.5, (-9.0, 9.0)), upper, curves=1, witnesses=[(1.5, 0.0, 1.5 / t)])
    yield regime("coaxial-beyond-box", cylinder(O, 1.5, (-9.0, 9.0)), cone(O, 0.0, a, (5.0, 8.0)))
    yield regime("axis-outside-both-nappes", cylinder(x_at(1.0), 0.4, (-9.0, 9.0)), both, curves=2, witnesses=[(1.4, 0.0, 1.4 / t), (1.4, 0.0, -1.4 / t)], hard_large=BUDGET_REASON)
    yield regime("axis-outside-one-nappe", cylinder(x_at(1.0), 0.4, (-9.0, 9.0)), upper, curves=1, witnesses=[(1.4, 0.0, 1.4 / t), (0.6, 0.0, 0.6 / t)])
    yield regime("cone-axis-inside", cylinder(x_at(0.4), 1.0, (-9.0, 9.0)), both, curves=2, witnesses=[(1.4, 0.0, 1.4 / t), (1.4, 0.0, -1.4 / t)])


def cylinder_torus():
    T = torus(O, 3.0, 1.0)
    z = math.sqrt(0.75)
    yield regime("coaxial-two", cylinder(O, 3.5), T, curves=2, witnesses=[(3.5, 0.0, z), (3.5, 0.0, -z)])
    yield regime("coaxial-crest", cylinder(O, 3.0), T, curves=2, witnesses=[(3.0, 0.0, 1.0), (3.0, 0.0, -1.0)])
    yield regime("coaxial-outer", cylinder(O, 4.0), T, curves=1, contact="tangential", witnesses=[(4.0, 0.0, 0.0)])
    yield regime("coaxial-inner", cylinder(O, 2.0), T, curves=1, contact="tangential", witnesses=[(2.0, 0.0, 0.0)])
    yield regime("coaxial-apart", cylinder(O, 5.0), T)
    yield regime("over-the-tube", cylinder(x_at(3.0), 0.5), T, curves=2, witnesses=[(3.5, 0.0, z), (3.5, 0.0, -z)])
    yield regime("past-the-outside", cylinder(x_at(3.5), 1.0), T, curves=1, witnesses=[(2.5, 0.0, z), (2.5, 0.0, -z)])
    x = 14.0 / 6.0
    y = math.sqrt(9.0 - x * x)
    yield regime("past-both-sides", cylinder(x_at(3.0), 2.0), T, curves=2, witnesses=[(x, y, 1.0), (x, -y, 1.0)])
    yield regime("touching-outside", cylinder(x_at(4.5), 0.5), T, points=1, witnesses=[(4.0, 0.0, 0.0)], hard=TANGENT_REASON)
    yield regime("apart", cylinder(x_at(6.0), 0.5), T)


def cylinder_ellipsoid():
    axes = (4.0, 3.0, 2.0)
    E = ellipsoid(O, *axes)

    def height(x):
        return 2.0 * math.sqrt(1.0 - x * x / 16.0)

    z = height(2.0)
    yield regime("inside-short-sides", cylinder(O, 2.0), E, curves=2, witnesses=[(2.0, 0.0, z), (2.0, 0.0, -z)], hard_large=BUDGET_REASON)
    yield regime("node", cylinder(O, 3.0), E, curves=1, witnesses=[(0.0, 3.0, 0.0)], hard=TANGENT_REASON)
    z = height(3.5)
    yield regime("between-sides", cylinder(O, 3.5), E, curves=2, witnesses=[(3.5, 0.0, z), (-3.5, 0.0, z)])
    yield regime("touching-long-axis", cylinder(O, 4.0), E, points=2, witnesses=[x_at(4.0), x_at(-4.0)], hard=TANGENT_REASON)
    yield regime("outside", cylinder(O, 5.0), E)
    S = ellipsoid(O, 3.0, 3.0, 2.0)
    z = 2.0 * math.sqrt(1.0 - 4.0 / 9.0)
    yield regime("spheroid-two", cylinder(O, 2.0), S, curves=2, witnesses=[(2.0, 0.0, z), (2.0, 0.0, -z)])
    yield regime("spheroid-equator", cylinder(O, 3.0), S, curves=1, contact="tangential", witnesses=[x_at(3.0)], hard=TANGENT_REASON)
    yield regime("spheroid-outside", cylinder(O, 4.0), S)


def cone_cone():
    a, b = math.pi / 6.0, math.pi / 4.0
    ta, tb = math.tan(a), math.tan(b)
    wide = cone(z_at(2.0), 0.0, b, (-9.0, 9.0))
    narrow = cone(O, 0.0, a, (-9.0, 9.0))
    z1, z2 = 2.0 * tb / (tb - ta), 2.0 * tb / (ta + tb)
    yield regime("coaxial-unequal", narrow, wide, curves=2, witnesses=[(z1 * ta, 0.0, z1), (z2 * ta, 0.0, z2)])
    yield regime("coaxial-equal", narrow, cone(z_at(2.0), 0.0, a, (-9.0, 9.0)), curves=1, witnesses=[(ta, 0.0, 1.0)])
    yield regime("coaxial-nappes-apart", cone(O, 0.0, a, (0.1, 9.0)), cone(z_at(-1.0), 0.0, a, (-9.0, -0.1)))
    yield regime("coincident", narrow, narrow, overlap=True, hard=COINCIDENT_REASON)


def cone_torus():
    a = math.pi / 4.0
    c = cone(O, 0.0, a)
    root = math.sqrt(0.5)
    yield regime("coaxial-two", c, torus(z_at(3.0), 3.0, 1.0), curves=2, witnesses=[(3 - root, 0.0, 3 - root), (3 + root, 0.0, 3 + root)])
    h = 3.0 + math.sqrt(2.0)
    foot = 3.0 + math.sqrt(0.5)
    yield regime("coaxial-touching", c, torus(z_at(h), 3.0, 1.0), curves=1, contact="tangential", witnesses=[(foot, 0.0, foot)], hard=TANGENT_REASON)
    yield regime("coaxial-apart", c, torus(z_at(6.0), 3.0, 1.0))


def cone_ellipsoid():
    a = math.pi / 4.0
    E = ellipsoid(O, 4.0, 3.0, 2.0)
    z = 1.0 / math.sqrt(1.0 / 16.0 + 1.0 / 4.0)
    yield regime("both-nappes", cone(O, 0.0, a), E, curves=2, witnesses=[(z, 0.0, z), (z, 0.0, -z)])
    yield regime("one-nappe", cone(O, 0.0, a, (0.1, 8.0)), E, curves=1, witnesses=[(z, 0.0, z)])
    yield regime("beyond-box", cone(O, 0.0, a, (6.0, 8.0)), E)


def torus_torus():
    A = torus(O, 3.0, 1.0)
    for h, name in ((1.0, "shifted-1"), (1.5, "shifted-1p5")):
        p = circle_circle((3.0, 0.0), 1.0, (3.0, h), 1.0)
        yield regime(name, A, torus(z_at(h), 3.0, 1.0), curves=2, witnesses=[(q[0], 0.0, q[1]) for q in p])
    yield regime("shifted-touching", A, torus(z_at(2.0), 3.0, 1.0), curves=1, contact="tangential", witnesses=[(3.0, 0.0, 1.0)], hard=TANGENT_REASON)
    yield regime("shifted-apart", A, torus(z_at(3.0), 3.0, 1.0))
    p = circle_circle((3.0, 0.0), 1.0, (4.5, 0.0), 1.0)
    yield regime("concentric-two", A, torus(O, 4.5, 1.0), curves=2, witnesses=[(q[0], 0.0, q[1]) for q in p])
    yield regime("concentric-touching", A, torus(O, 5.0, 1.0), curves=1, contact="tangential", witnesses=[(4.0, 0.0, 0.0)])
    yield regime("concentric-apart", A, torus(O, 6.0, 1.0))
    yield regime("coincident", A, A, overlap=True)


def torus_ellipsoid():
    T = torus(O, 3.0, 1.0)

    def crossings(a, c):
        """Where the tube's circle (3 + cos v, sin v) meets the ellipse, upper half."""
        f = lambda v: ((3.0 + math.cos(v)) / a) ** 2 + (math.sin(v) / c) ** 2 - 1.0
        lo, hi = 0.0, math.pi
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if f(mid) > 0.0 else (lo, mid)
        v = 0.5 * (lo + hi)
        return 3.0 + math.cos(v), math.sin(v)

    for a, c, name in ((3.0, 1.5, "through-the-tube"), (2.5, 2.0, "wide-through-the-tube")):
        x, z = crossings(a, c)
        yield regime(name, T, ellipsoid(O, a, a, c), curves=2, witnesses=[(x, 0.0, z), (x, 0.0, -z)])
    yield regime("in-the-hole", T, ellipsoid(O, 1.5, 1.5, 1.2))


def ellipsoid_ellipsoid():
    E = ellipsoid(O, 4.0, 3.0, 2.0)

    def along_x(d):
        s = math.sqrt(1.0 - d * d / 64.0)
        return [(d / 2.0, 3.0 * s, 0.0), (d / 2.0, -3.0 * s, 0.0), (d / 2.0, 0.0, 2.0 * s)]

    for d in (2.0, 4.0, 7.0):
        yield regime(f"shifted-x-{int(d)}", E, ellipsoid(x_at(d), 4.0, 3.0, 2.0), curves=1, witnesses=along_x(d))
    yield regime("shifted-x-touching", E, ellipsoid(x_at(8.0), 4.0, 3.0, 2.0), points=1, witnesses=[x_at(4.0)], hard=TANGENT_REASON)
    yield regime("shifted-x-apart", E, ellipsoid(x_at(9.0), 4.0, 3.0, 2.0))
    s = math.sqrt(0.75)
    yield regime("shifted-z-2", E, ellipsoid(z_at(2.0), 4.0, 3.0, 2.0), curves=1, witnesses=[(4 * s, 0.0, 1.0), (0.0, 3 * s, 1.0)])
    yield regime("nested", E, ellipsoid(O, 2.0, 1.5, 1.0))
    yield regime("coincident", E, E, overlap=True, hard=COINCIDENT_REASON)


FAMILIES = {
    "plane-cone": plane_cone,
    "plane-torus": plane_torus,
    "plane-ellipsoid": plane_ellipsoid,
    "sphere-cone": sphere_cone,
    "sphere-torus": sphere_torus,
    "sphere-ellipsoid": sphere_ellipsoid,
    "cylinder-cone": cylinder_cone,
    "cylinder-torus": cylinder_torus,
    "cylinder-ellipsoid": cylinder_ellipsoid,
    "cone-cone": cone_cone,
    "cone-torus": cone_torus,
    "cone-ellipsoid": cone_ellipsoid,
    "torus-torus": torus_torus,
    "torus-ellipsoid": torus_ellipsoid,
    "ellipsoid-ellipsoid": ellipsoid_ellipsoid,
}


def rotated(placement, v):
    x, y, z, _ = placement
    return tuple(v[0] * x[i] + v[1] * y[i] + v[2] * z[i] for i in range(3))


def placed_point(placement, scale, p):
    # Rounded far below the tolerance: the last digits of a platform's
    # trigonometry would otherwise change the files between machines.
    q = rotated(placement, tuple(scale * c for c in p))
    return [round(q[i] + placement[3][i], 10) for i in range(3)]


def placed_axis(placement, v):
    return [round(c, 14) for c in rotated(placement, v)]


LENGTHS = ("radius", "major", "minor", "a", "b", "c")
POINTS = ("origin", "center")


def placed_surface(surface, name, placement, scale):
    out = {"name": name}
    for key, value in surface.items():
        if key in POINTS:
            out[key] = placed_point(placement, scale, value)
        elif key in ("x_axis", "y_axis"):
            out[key] = placed_axis(placement, value)
        elif key in LENGTHS:
            out[key] = value * scale
        elif key == "v":
            out["box"] = [[0.0, TAU], [value[0] * scale, value[1] * scale]]
        else:
            out[key] = value
    return out


def case_for(family, r, placement_name, scale):
    placement = PLACEMENTS[placement_name]
    scale_name = f"x{int(scale)}"
    case_id = f"curved-ssi/{family}/{r['name']}-{placement_name}-{scale_name}"
    expected = {
        "overlap": r["overlap"],
        "curves": r["curves"],
        "points": r["points"],
        "contact": r["contact"],
        "witnesses": [placed_point(placement, scale, w) for w in r["witnesses"]],
    }
    reason = r["hard"]
    if scale == SCALES[-1] and r["hard_large"]:
        reason = r["hard_large"]
    if reason:
        expected["may_fail_closed"] = reason
    a, b = family.split("-")
    return {
        "format": "stereon-corpus-case/1",
        "id": case_id,
        "kind": "hard" if reason else "reference",
        "summary": f"A {a} and a {b}, {r['name']}, placed at {placement_name}, scaled by {int(scale)}.",
        "added": "2026-11-09",
        "source": "Generated by tools/generate_surface_corpus.py",
        "licence": "MPL-2.0",
        "references": ["ADR-0010", "ADR-0029"],
        "tags": ["generated", r["name"]],
        "inputs": {
            "surfaces": [
                placed_surface(r["a"], "A", placement, scale),
                placed_surface(r["b"], "B", placement, scale),
            ],
            "tolerance": TOLERANCE,
            "scale": scale,
        },
        "operation": {"name": "intersect", "operands": ["A", "B"]},
        "expected": expected,
    }


def all_cases():
    for family, regimes in FAMILIES.items():
        for r in regimes():
            for placement_name in r["placements"]:
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
    suite = args.root / "corpus" / "curved-ssi"
    if suite.exists():
        for path in sorted(suite.glob("*/*.json")):
            if path.resolve() not in expected_paths:
                if args.check:
                    stale.append(str(path.relative_to(args.root)))
                else:
                    path.unlink()
    if stale:
        print(f"generate_surface_corpus: {len(stale)} stale case(s), first {stale[0]}", file=sys.stderr)
        return 1
    print(f"generate_surface_corpus: {count} cases {'match' if args.check else 'written'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
