#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""Generate the parametric sweeps of the analytic-ssi corpus (RS-0001).

Every Tier A degeneracy of RS-0001's table is swept over the offsets
delta / tau in {0, 0.1, 0.5, 0.9, 0.99, 1.01, 1.1, 2, 10} from the degenerate
configuration, on both sides of it where the sides differ, over the
tolerances tau in {1e-5, 1e-4, 1e-3, 1e-2} mm, and over two placements: at
the origin, and translated to (480, -300, 210) m, near the edge of the
+-500 m size box (ADR-0001). Frames are aligned with the coordinate axes, so an offset of
zero is an exact degeneracy of the doubles at either placement.

Each case is written to corpus/analytic-ssi/sweep-<degeneracy>/<case>.json.
Its nominal configuration is computed by check_corpus.py's exact rational
oracle from the doubles written, never from the intended offset, so a
rounded offset far from the origin is judged as the doubles it became. Its
tolerant expectation follows RS-0001's contract: within the tolerance of
the degeneracy, the degenerate configuration, snapped unless the doubles
are exactly degenerate; beyond it, the nominal one.

Usage: generate_corpus.py [--root <repository root>] [--check]

With --check, nothing is written: the script fails if a generated case is
missing from the corpus or differs from what it would write, or a file in a
sweep directory is not one it generates.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_corpus  # noqa: E402

ADDED = "2026-11-08"
RATIOS = ("0", "0.1", "0.5", "0.9", "0.99", "1.01", "1.1", "2", "10")
TOLERANCES = ("1e-5", "1e-4", "1e-3", "1e-2")
PLACEMENTS = {
    "origin": (0.0, 0.0, 0.0),
    # Near the edge of the +-500 m size box of ADR-0001, along no axis, so
    # every coordinate carries an offset.
    "far": (480000.0, -300000.0, 210000.0),
}
NEAR_THRESHOLD = {"0.9", "0.99", "1.01", "1.1"}
Z = (0.0, 0.0, 1.0)
X = (1.0, 0.0, 0.0)


def moved(base: tuple, offset: tuple) -> list[float]:
    return [base[0] + offset[0], base[1] + offset[1], base[2] + offset[2]]


def plane(point, normal):
    return {"type": "plane", "point": point, "normal": list(normal)}


def sphere(centre, radius):
    return {"type": "sphere", "centre": centre, "radius": radius}


def cylinder(point, axis, radius):
    return {"type": "cylinder", "point": point, "axis": list(axis), "radius": radius}


# Each degeneracy: (description, degenerate configuration, the sides it
# has, a function of (placement, signed offset) giving the two surfaces).
# The offset grows the gap: a positive one separates the surfaces or moves
# a tangency outward, a negative one makes them cross.
DEGENERACIES = {
    "plane-plane-parallel": (
        "parallel planes", "coincident", (1,),
        lambda o, d: (plane(moved(o, (0, 0, 0)), Z), plane(moved(o, (0, 0, d)), Z)),
    ),
    "plane-sphere-tangent": (
        "a plane and a sphere of radius 2", "tangent_point", (1, -1),
        lambda o, d: (plane(moved(o, (0, 0, 2 + d)), Z), sphere(moved(o, (0, 0, 0)), 2.0)),
    ),
    "plane-cylinder-tangent": (
        "a plane parallel to the axis of a cylinder of radius 2", "tangent_line", (1, -1),
        lambda o, d: (plane(moved(o, (2 + d, 0, 0)), X), cylinder(moved(o, (0, 0, 0)), Z, 2.0)),
    ),
    "sphere-sphere-external": (
        "spheres of radii 2 and 1.5 side by side", "tangent_point_external", (1, -1),
        lambda o, d: (sphere(moved(o, (0, 0, 0)), 2.0), sphere(moved(o, (3.5 + d, 0, 0)), 1.5)),
    ),
    "sphere-sphere-internal": (
        "a sphere of radius 1.5 inside one of radius 2", "tangent_point_internal", (1, -1),
        lambda o, d: (sphere(moved(o, (0, 0, 0)), 2.0), sphere(moved(o, (0.5 - d, 0, 0)), 1.5)),
    ),
    "sphere-cylinder-coaxial": (
        "a sphere of radius 2 centred on the axis of a cylinder", "tangent_circle", (1, -1),
        lambda o, d: (sphere(moved(o, (0, 0, 0)), 2.0), cylinder(moved(o, (0, 0, 0)), Z, 2.0 + d)),
    ),
    "cylinder-cylinder-external": (
        "parallel cylinders of radii 2 and 1.5 side by side", "tangent_line_external", (1, -1),
        lambda o, d: (cylinder(moved(o, (0, 0, 0)), Z, 2.0), cylinder(moved(o, (3.5 + d, 0, 0)), Z, 1.5)),
    ),
    "cylinder-cylinder-internal": (
        "a cylinder of radius 1.5 inside a parallel one of radius 2", "tangent_line_internal", (1, -1),
        lambda o, d: (cylinder(moved(o, (0, 0, 0)), Z, 2.0), cylinder(moved(o, (0.5 - d, 0, 0)), Z, 1.5)),
    ),
}


def name_of(text: str) -> str:
    return text.replace(".", "p")


def case_for(family: str, side: int, ratio: str, tolerance: str, placement: str) -> tuple[str, dict]:
    description, degenerate, _, make = DEGENERACIES[family]
    tau = float(tolerance)
    delta = side * float(ratio) * tau
    first, second = make(PLACEMENTS[placement], delta)
    first = {"name": "A", **first, "tolerance": tau}
    second = {"name": "B", **second, "tolerance": tau}
    nominal = check_corpus.classify_pair(first, second)
    assert nominal["tier"] == "A", family
    within = Fraction(ratio) < 1
    configuration = degenerate if within else nominal["configuration"]
    snapped = within and nominal["configuration"] != degenerate
    where = "at the origin" if placement == "origin" else "near the edge of the size box"
    if ratio == "0":
        stem = f"touch-tol-{tolerance}-{placement}"
        offset = "exactly degenerate"
    else:
        stem = f"{'out' if side > 0 else 'in'}-{name_of(ratio)}tol-tol-{tolerance}-{placement}"
        offset = f"{ratio} tolerances {'apart' if side > 0 else 'across'}"
    case_id = f"analytic-ssi/sweep-{family}/{stem}"
    tags = ["sweep", "coincidence" if degenerate == "coincident" else "tangency"]
    tags.append("snap" if snapped else ("threshold" if ratio in NEAR_THRESHOLD else "robust"))
    if placement == "far":
        tags.append("size-box")
    case = {
        "format": check_corpus.CASE_FORMAT,
        "id": case_id,
        "kind": "hard" if ratio in NEAR_THRESHOLD or ratio == "0" else "reference",
        "summary": f"Sweep of {description} {where}: {offset} from {degenerate}, tolerance {tolerance} mm.",
        "added": ADDED,
        "source": "Generated by tools/generate_corpus.py (the parametric sweeps of RS-0001)",
        "licence": "MPL-2.0",
        "references": ["RS-0001", "ADR-0010"],
        "tags": tags,
        "inputs": {"surfaces": [first, second]},
        "operation": {"name": "intersect", "operands": ["A", "B"]},
        "expected": {
            "nominal": nominal,
            "tolerant": {"status": "proposed", "configuration": configuration, "snapped": snapped},
        },
    }
    return case_id, case


def generate() -> dict[str, dict]:
    cases = {}
    for family, (_, _, sides, _) in DEGENERACIES.items():
        for placement in PLACEMENTS:
            for tolerance in TOLERANCES:
                for ratio in RATIOS:
                    for side in sides if ratio != "0" else (1,):
                        case_id, case = case_for(family, side, ratio, tolerance, placement)
                        assert case_id not in cases, case_id
                        cases[case_id] = case
    return cases


SHORT_ARRAY = re.compile(r"\[\s*([^\[\]{}]*?)\s*\]")


def text_of(case: dict) -> str:
    """Indented JSON with arrays of numbers or strings on one line, as the
    seed cases are written."""
    text = json.dumps(case, indent=2)
    text = SHORT_ARRAY.sub(lambda m: "[" + re.sub(r",\s*", ", ", m.group(1)) + "]", text)
    return text + "\n"


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)
    corpus = args.root / "corpus"
    cases = generate()
    wanted = {corpus / (case_id + ".json"): text_of(case) for case_id, case in cases.items()}
    if args.check:
        problems = []
        for path, text in wanted.items():
            if not path.exists() or path.read_text(encoding="utf-8") != text:
                problems.append(f"{path.relative_to(args.root)}: missing or not as generated")
        for family in DEGENERACIES:
            for path in sorted((corpus / "analytic-ssi" / f"sweep-{family}").glob("*")):
                if path not in wanted:
                    problems.append(f"{path.relative_to(args.root)}: not a generated case")
        for problem in problems:
            print(problem)
        print(f"generate_corpus: {len(wanted)} generated cases, {len(problems)} problem(s)")
        return 1 if problems else 0
    for path, text in wanted.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8", newline="\n")
    print(f"generate_corpus: wrote {len(wanted)} cases")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
