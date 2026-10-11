#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""Validate the geometry corpus (corpus/README.md) and its ground truth.

Every case file corpus/<suite>/<family>/<case>.json is checked against the
case format stereon-corpus-case/1: required fields, an id equal to its path,
and a licence row in corpus/LICENSES.md for anything not MPL-2.0.

For the analytic-ssi suite (RS-0001) the check also recomputes the expected
nominal configuration of every surface pair. Each input double is converted
to an exact fraction and every decision is the sign of a polynomial in those
fractions, so the result is the configuration of the doubles themselves,
not of the decimal text and not of a floating-point evaluation. The
tolerant expectation is the claim under test in RS-0001 and is only
checked for shape.

Usage: check_corpus.py [--root <repository root>] | --self-test
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from fractions import Fraction
from pathlib import Path

CASE_FORMAT = "stereon-corpus-case/1"
KINDS = {"regression", "hard", "reference"}
REQUIRED = ("format", "id", "kind", "summary", "added", "source", "licence", "inputs", "operation", "expected")
DATE = re.compile(r"\d{4}-\d{2}-\d{2}")
REFERENCE = re.compile(r"(?:ADR|RS)-\d{4}")

SURFACE_FIELDS = {
    "plane": ("point", "normal"),
    "sphere": ("centre", "radius"),
    "cylinder": ("point", "axis", "radius"),
}
TOLERANT_STATUS = {"proposed", "open"}
RESOLUTION = 1e-5  # ADR-0001 linear resolution, mm

Vec = tuple[Fraction, Fraction, Fraction]


# --- Exact vector arithmetic -------------------------------------------------


def vec(values: list[float]) -> Vec:
    x, y, z = (Fraction(v) for v in values)
    return (x, y, z)


def sub(a: Vec, b: Vec) -> Vec:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def dot(a: Vec, b: Vec) -> Fraction:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a: Vec, b: Vec) -> Vec:
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def is_zero(a: Vec) -> bool:
    return a == (0, 0, 0)


def compare(a: Fraction, b: Fraction) -> int:
    return (a > b) - (a < b)


# --- Nominal oracle for surface pairs ----------------------------------------


def classify_pair(first: dict, second: dict) -> dict:
    """Return {"tier": "A", "configuration": ...} or {"tier": "B"} for the exact input doubles."""
    order = ("plane", "sphere", "cylinder")
    a, b = sorted((first, second), key=lambda s: order.index(s["type"]))
    kinds = (a["type"], b["type"])
    if kinds == ("plane", "plane"):
        return {"tier": "A", "configuration": plane_plane(a, b)}
    if kinds == ("plane", "sphere"):
        return {"tier": "A", "configuration": plane_sphere(a, b)}
    if kinds == ("plane", "cylinder"):
        return {"tier": "A", "configuration": plane_cylinder(a, b)}
    if kinds == ("sphere", "sphere"):
        return {"tier": "A", "configuration": sphere_sphere(a, b)}
    if kinds == ("sphere", "cylinder"):
        return sphere_cylinder(a, b)
    return cylinder_cylinder(a, b)


def plane_plane(a: dict, b: dict) -> str:
    n1, n2 = vec(a["normal"]), vec(b["normal"])
    if not is_zero(cross(n1, n2)):
        return "line"
    return "coincident" if dot(n1, sub(vec(b["point"]), vec(a["point"]))) == 0 else "disjoint"


def plane_sphere(plane: dict, sphere: dict) -> str:
    n = vec(plane["normal"])
    s = dot(n, sub(vec(sphere["centre"]), vec(plane["point"])))
    r = Fraction(sphere["radius"])
    return {1: "disjoint", 0: "tangent_point", -1: "circle"}[compare(s * s, r * r * dot(n, n))]


def plane_cylinder(plane: dict, cylinder: dict) -> str:
    n, axis = vec(plane["normal"]), vec(cylinder["axis"])
    if dot(n, axis) == 0:
        s = dot(n, sub(vec(cylinder["point"]), vec(plane["point"])))
        r = Fraction(cylinder["radius"])
        return {1: "disjoint", 0: "tangent_line", -1: "two_lines"}[compare(s * s, r * r * dot(n, n))]
    return "circle" if is_zero(cross(n, axis)) else "ellipse"


def two_round(distance_squared: Fraction, r1: Fraction, r2: Fraction, contact: str, crossing: str) -> str:
    """Configuration of two circles of radii r1, r2 whose centres are sqrt(distance_squared) apart."""
    if distance_squared == 0 and r1 == r2:
        return "coincident"
    outer = compare(distance_squared, (r1 + r2) ** 2)
    if outer > 0:
        return "disjoint"
    if outer == 0:
        return f"{contact}_external"
    inner = compare(distance_squared, (r1 - r2) ** 2)
    if inner < 0:
        return "disjoint"
    return f"{contact}_internal" if inner == 0 else crossing


def sphere_sphere(a: dict, b: dict) -> str:
    d = sub(vec(b["centre"]), vec(a["centre"]))
    return two_round(dot(d, d), Fraction(a["radius"]), Fraction(b["radius"]), "tangent_point", "circle")


def sphere_cylinder(sphere: dict, cylinder: dict) -> dict:
    axis = vec(cylinder["axis"])
    if not is_zero(cross(axis, sub(vec(sphere["centre"]), vec(cylinder["point"])))):
        return {"tier": "B"}
    rs, rc = Fraction(sphere["radius"]), Fraction(cylinder["radius"])
    return {"tier": "A", "configuration": {1: "two_circles", 0: "tangent_circle", -1: "disjoint"}[compare(rs, rc)]}


def cylinder_cylinder(a: dict, b: dict) -> dict:
    a1, a2 = vec(a["axis"]), vec(b["axis"])
    if not is_zero(cross(a1, a2)):
        return {"tier": "B"}
    w = cross(a1, sub(vec(b["point"]), vec(a["point"])))
    distance_squared = dot(w, w) / dot(a1, a1)
    configuration = two_round(distance_squared, Fraction(a["radius"]), Fraction(b["radius"]), "tangent_line",
                              "two_lines")
    return {"tier": "A", "configuration": configuration}


# --- Case validation ---------------------------------------------------------


def check_vector(value, label: str, nonzero: bool) -> list[str]:
    if not (isinstance(value, list) and len(value) == 3 and all(isinstance(v, (int, float)) for v in value)):
        return [f"{label} must be a list of three numbers"]
    if not all(math.isfinite(v) for v in value):
        return [f"{label} must be finite"]
    if nonzero and all(v == 0 for v in value):
        return [f"{label} must not be the zero vector"]
    return []


def check_surfaces(surfaces) -> tuple[dict, list[str]]:
    if not isinstance(surfaces, list) or not surfaces:
        return {}, ["inputs.surfaces must be a non-empty list"]
    by_name: dict[str, dict] = {}
    problems = []
    for index, surface in enumerate(surfaces):
        label = f"inputs.surfaces[{index}]"
        name, kind = surface.get("name"), surface.get("type")
        if not isinstance(name, str) or name in by_name:
            problems.append(f"{label}: missing or duplicate name")
            continue
        if kind not in SURFACE_FIELDS:
            problems.append(f"{label}: type must be one of {', '.join(SURFACE_FIELDS)}")
            continue
        for field in SURFACE_FIELDS[kind]:
            if field == "radius":
                radius = surface.get("radius")
                if not isinstance(radius, (int, float)) or not math.isfinite(radius) or radius <= 0:
                    problems.append(f"{label}: radius must be a positive number")
            else:
                problems += check_vector(surface.get(field), f"{label}.{field}", field in ("normal", "axis"))
        tolerance = surface.get("tolerance")
        if not isinstance(tolerance, (int, float)) or not tolerance >= RESOLUTION:
            problems.append(f"{label}: tolerance must be a number >= the resolution {RESOLUTION} mm")
        by_name[name] = surface
    return by_name, problems


def check_tolerant(tolerant, is_pair: bool) -> list[str]:
    if not isinstance(tolerant, dict) or tolerant.get("status") not in TOLERANT_STATUS:
        return [f"expected.tolerant.status must be one of {', '.join(sorted(TOLERANT_STATUS))}"]
    if tolerant["status"] == "open":
        return [] if isinstance(tolerant.get("question"), str) else ["an open tolerant expectation names its question"]
    if "outcome" in tolerant:
        return [] if isinstance(tolerant["outcome"], str) else ["expected.tolerant.outcome must be a string"]
    if is_pair:
        if not isinstance(tolerant.get("configuration"), str) or not isinstance(tolerant.get("snapped"), bool):
            return ["a proposed pair expectation has a configuration and snapped: true|false"]
        return []
    if not isinstance(tolerant.get("pairs"), dict) or not isinstance(tolerant.get("snaps"), list):
        return ["a proposed arrangement expectation has pairs and snaps"]
    return []


def check_analytic_ssi(case: dict) -> list[str]:
    surfaces, problems = check_surfaces(case["inputs"].get("surfaces"))
    if problems:
        return problems
    operation = case["operation"]
    operands = operation.get("operands")
    if operation.get("name") not in ("intersect", "arrange"):
        return ["operation.name must be intersect or arrange"]
    if not isinstance(operands, list) or any(o not in surfaces for o in operands):
        return ["operation.operands must name surfaces from inputs.surfaces"]
    is_pair = operation["name"] == "intersect"
    if is_pair and len(operands) != 2:
        return ["intersect takes exactly two operands"]
    if not is_pair and len(operands) < 3:
        return ["arrange takes three or more operands"]

    nominal = case["expected"].get("nominal")
    problems += check_tolerant(case["expected"].get("tolerant"), is_pair)
    if is_pair:
        computed = classify_pair(surfaces[operands[0]], surfaces[operands[1]])
        recorded = {k: v for k, v in (nominal or {}).items() if k in ("tier", "configuration")}
        if computed["tier"] == "B":
            recorded.pop("configuration", None)
        if recorded != computed:
            problems.append(f"expected.nominal is {recorded}, but the exact input doubles give {computed}")
        return problems

    pairs = (nominal or {}).get("pairs")
    if not isinstance(pairs, dict):
        return problems + ["expected.nominal.pairs must map 'A-B' to a configuration"]
    for i, first in enumerate(operands):
        for second in operands[i + 1:]:
            key = f"{first}-{second}"
            computed = classify_pair(surfaces[first], surfaces[second])
            expected = computed.get("configuration", "tier_b")
            if pairs.get(key) != expected:
                problems.append(f"expected.nominal.pairs[{key}] is {pairs.get(key)!r}, "
                                f"but the exact input doubles give {expected!r}")
    return problems


SUITE_CHECKS = {"analytic-ssi": check_analytic_ssi}


def check_case(case: dict, case_id: str, licensed: str) -> list[str]:
    missing = [field for field in REQUIRED if field not in case]
    if missing:
        return [f"missing field(s): {', '.join(missing)}"]
    problems = []
    if case["format"] != CASE_FORMAT:
        problems.append(f"format must be {CASE_FORMAT}")
    if case["id"] != case_id:
        problems.append(f"id {case['id']!r} does not match the file path {case_id!r}")
    if case["kind"] not in KINDS:
        problems.append(f"kind must be one of {', '.join(sorted(KINDS))}")
    if not isinstance(case["added"], str) or not DATE.fullmatch(case["added"]):
        problems.append("added must be a date YYYY-MM-DD")
    if case["licence"] != "MPL-2.0" and f"`{case_id}`" not in licensed:
        problems.append(f"licence {case['licence']} needs a row for `{case_id}` in corpus/LICENSES.md")
    for reference in case.get("references", []):
        if not REFERENCE.fullmatch(reference):
            problems.append(f"reference {reference!r} is not ADR-NNNN or RS-NNNN")
    suite_check = SUITE_CHECKS.get(case_id.split("/")[0])
    if suite_check and not problems:
        problems += suite_check(case)
    return problems


def check_tree(root: Path) -> int:
    corpus = root / "corpus"
    licences = corpus / "LICENSES.md"
    if not licences.is_file():
        print("check_corpus: corpus/LICENSES.md is missing", file=sys.stderr)
        return 1
    licensed = licences.read_text(encoding="utf-8")
    count = 0
    files = 0
    for path in sorted(corpus.rglob("*.json")):
        files += 1
        relative = path.relative_to(corpus).with_suffix("").as_posix()
        if relative.count("/") != 2:
            problems = ["case files live at corpus/<suite>/<family>/<case>.json"]
        else:
            try:
                case = json.loads(path.read_text(encoding="utf-8"))
            except json.JSONDecodeError as error:
                problems = [f"not valid JSON: {error}"]
            else:
                problems = check_case(case, relative, licensed)
        for message in problems:
            print(f"{path.relative_to(root).as_posix()}: {message}", file=sys.stderr)
            count += 1
    if count:
        print(f"check_corpus: {count} problem(s)", file=sys.stderr)
        return 1
    print(f"check_corpus: {files} case(s), no problems")
    return 0


def plane(point, normal):
    return {"type": "plane", "point": point, "normal": normal}


def sphere(centre, radius):
    return {"type": "sphere", "centre": centre, "radius": radius}


def cylinder(point, axis, radius):
    return {"type": "cylinder", "point": point, "axis": axis, "radius": radius}


SELF_TEST = [
    (plane([0, 0, 0], [0, 0, 1]), plane([0, 0, 0], [1, 0, 0]), "line"),
    (plane([0, 0, 5], [0, 0, 1]), plane([3, 4, 5], [0, 0, 2]), "coincident"),
    (plane([0, 0, 0], [0, 0, 1]), plane([0, 0, 1e-9], [0, 0, -1]), "disjoint"),
    (plane([0, 0, 10], [0, 0, 1]), sphere([0, 0, 0], 10), "tangent_point"),
    (plane([0, 0, 0.3], [0, 0, 1]), sphere([0, 0, 0.2], 0.1), "circle"),  # 0.3 - 0.2 < 0.1 in doubles
    (plane([0, 0, 0.1], [0, 0, 1]), sphere([0, 0, 0], 0.1), "tangent_point"),
    (plane([0, 0, 5], [0, 0, 1]), cylinder([0, 0, 0], [1, 0, 0], 5), "tangent_line"),
    (plane([0, 0, 2], [0, 0, 1]), cylinder([0, 0, 0], [1, 0, 0], 5), "two_lines"),
    (plane([0, 0, 3], [0, 0, 1]), cylinder([0, 0, 0], [0, 0, 7], 5), "circle"),
    (plane([0, 0, 0], [0, 1, 1]), cylinder([0, 0, 0], [0, 0, 1], 5), "ellipse"),
    (sphere([0, 0, 0], 3), sphere([5, 0, 0], 2), "tangent_point_external"),
    (sphere([0, 0, 0], 5), sphere([3, 0, 0], 2), "tangent_point_internal"),
    (sphere([0, 0, 0], 5), sphere([1, 0, 0], 2), "disjoint"),
    (sphere([0, 0, 0], 5), sphere([0, 0, 0], 5), "coincident"),
    (sphere([0, 0, 1], 5), cylinder([0, 0, 0], [0, 0, 1], 3), "two_circles"),
    (sphere([0, 0, 1], 3), cylinder([0, 0, 0], [0, 0, 1], 3), "tangent_circle"),
    (sphere([1e-9, 0, 0], 5), cylinder([0, 0, 0], [0, 0, 1], 3), None),
    (cylinder([0, 0, 0], [0, 0, 1], 2), cylinder([5, 0, 0], [0, 0, -3], 3), "tangent_line_external"),
    (cylinder([0, 0, 0], [0, 0, 1], 3), cylinder([0, 0, 9], [0, 0, 2], 3), "coincident"),
    (cylinder([0, 0, 0], [0, 0, 1], 3), cylinder([1, 0, 0], [0, 0, 1], 2), "tangent_line_internal"),
    (cylinder([0, 0, 0], [0, 0, 1], 3), cylinder([5, 0, 0], [1e-12, 0, 1], 3), None),
]


def self_test() -> int:
    failures = []
    for first, second, expected in SELF_TEST:
        for a, b in ((first, second), (second, first)):
            result = classify_pair(a, b)
            got = result.get("configuration") if result["tier"] == "A" else None
            if got != expected:
                failures.append(f"{a['type']}-{b['type']}: expected {expected}, got {got}")
    for failure in failures:
        print("self-test: " + failure, file=sys.stderr)
    if failures:
        return 1
    print(f"check_corpus self-test: {len(SELF_TEST)} surface pairs pass in both orders")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument("--self-test", action="store_true", help="run the built-in cases and exit")
    args = parser.parse_args()
    if args.self_test:
        return self_test()
    return check_tree(args.root.resolve())


if __name__ == "__main__":
    sys.exit(main())
