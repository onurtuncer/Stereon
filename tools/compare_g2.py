#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""Compare the G2 records of two platforms (ADR-0007, G2).

`stereon_determinism --report <file>` writes, for every Boolean of
tests/support/BooleanCases.hpp, the result's entity counts and every
vertex's position and tolerance in index order, with a digest of its bits.
Two platforms' files agree when every case is present in both, the entity
counts are equal, and every vertex lies within the larger of its two
tolerances of the other's. How many results are bit-identical is reported
and never fails the comparison: G1 promises identical bits on one
executable and one machine, and two compilers' libms need not round alike.
Volumes and areas join the comparison when the kernel computes them (the
Boolean prototype's M11).

Usage: compare_g2.py <first.json> <second.json> [--summary <markdown>]
       compare_g2.py --self-test
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path


def finite(value: object) -> bool:
    """A real number that is neither infinite nor NaN (booleans excluded)."""
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def validate(report: object, label: str) -> list[str]:
    """What makes a report unfit to compare: empty, not a list of cases, a
    case named twice, or a case whose record is incomplete. A complete case
    has a name, counts of every entity kind as non-negative integers, as
    many vertices as its vertex count, each three finite coordinates and a
    finite non-negative tolerance, and its digest. Comparing an incomplete
    report would pass whatever it leaves out."""
    if not isinstance(report, list) or not report:
        return [f"the {label} report has no cases"]
    problems: list[str] = []
    names: set[str] = set()
    for position, case in enumerate(report):
        name = case.get("name") if isinstance(case, dict) else None
        where = f"the {label} report's case {name if isinstance(name, str) else position}"
        if not isinstance(name, str) or not name:
            problems.append(f"{where}: no name")
            continue
        if name in names:
            problems.append(f"{where}: named twice")
        names.add(name)
        counts = case.get("counts")
        if not isinstance(counts, dict) or not counts or not all(
                isinstance(n, int) and not isinstance(n, bool) and n >= 0 for n in counts.values()):
            problems.append(f"{where}: no counts, or a count that is not a non-negative integer")
            continue
        vertices = case.get("vertices")
        if not isinstance(vertices, list) or len(vertices) != counts.get("vertex"):
            problems.append(f"{where}: the vertices do not match the vertex count")
            continue
        for index, vertex in enumerate(vertices):
            if not (isinstance(vertex, list) and len(vertex) == 4 and all(finite(c) for c in vertex)
                    and vertex[3] >= 0.0):
                problems.append(f"{where}: vertex {index} is not three finite coordinates and a "
                                "finite non-negative tolerance")
        if not isinstance(case.get("bits"), str):
            problems.append(f"{where}: no digest")
    return problems


def compare(first: list[dict], second: list[dict]) -> tuple[list[str], int]:
    """The problems found, and how many cases are bit-identical. Both reports
    are validated first; an invalid one is not compared."""
    problems = validate(first, "first") + validate(second, "second")
    if problems:
        return problems, 0
    identical = 0
    theirs = {case["name"]: case for case in second}
    for case in first:
        name = case["name"]
        other = theirs.pop(name, None)
        if other is None:
            problems.append(f"{name}: missing from the second report")
            continue
        if case["counts"] != other["counts"]:
            problems.append(f"{name}: entity counts differ, {case['counts']} against {other['counts']}")
            continue
        # Equal counts and validated reports: the vertex lists are equally long.
        for index, (mine, yours) in enumerate(zip(case["vertices"], other["vertices"], strict=True)):
            tolerance = max(mine[3], yours[3])
            distance = math.dist(mine[:3], yours[:3])
            if not distance <= tolerance:
                problems.append(f"{name}: vertex {index} is {distance:.3g} mm from the other "
                                f"platform's, beyond the tolerance {tolerance:.3g} mm")
        identical += case["bits"] == other["bits"]
    for name in theirs:
        problems.append(f"{name}: missing from the first report")
    return problems, identical


def summary(first: str, second: str, cases: int, problems: list[str], identical: int) -> str:
    lines = [
        "## G2: the Booleans across platforms (ADR-0007)",
        "",
        f"{first} against {second}: {cases} results, entity counts and every vertex within the "
        "larger of its two tolerances.",
        "",
        f"- Agreeing: {'yes' if not problems else 'no'}",
        f"- Bit-identical: {identical} of {cases} (reported, not required)",
    ]
    lines += [f"- {problem}" for problem in problems]
    return "\n".join(lines) + "\n"


def self_test() -> int:
    vertex = [1.0, 2.0, 3.0, 1e-5]
    base = [{"name": "a", "counts": {"vertex": 1}, "vertices": [vertex], "bits": "x"}]
    nudged = [{"name": "a", "counts": {"vertex": 1}, "vertices": [[1.0, 2.0, 3.0 + 5e-6, 1e-5]],
               "bits": "y"}]
    far = [{"name": "a", "counts": {"vertex": 1}, "vertices": [[1.0, 2.0, 3.1, 1e-5]], "bits": "y"}]
    two = [{"name": "a", "counts": {"vertex": 2}, "vertices": [vertex, vertex], "bits": "x"}]
    # Each report claims what it holds; these leave something out or hold
    # what no tolerance can compare, and must fail even against themselves.
    missing = [{"name": "a", "counts": {"vertex": 2}, "vertices": [vertex], "bits": "x"}]
    infinite = [{"name": "a", "counts": {"vertex": 1}, "vertices": [[1.0, 2.0, 3.0, math.inf]],
                 "bits": "x"}]
    negative = [{"name": "a", "counts": {"vertex": 1}, "vertices": [[1.0, 2.0, 3.0, -1.0]],
                 "bits": "x"}]
    nan = [{"name": "a", "counts": {"vertex": 1}, "vertices": [[math.nan, 2.0, 3.0, 1e-5]],
            "bits": "x"}]
    twice = base + base
    cases = [
        ("identical", base, base, 0, 1),
        ("within tolerance", base, nudged, 0, 0),
        ("beyond tolerance", base, far, 1, 0),
        ("counts differ", base, two, 1, 0),
        ("case missing", base, [base[0] | {"name": "b"}], 2, 0),
        ("second report empty", base, [], 1, 0),
        ("both reports empty", [], [], 2, 0),
        ("a vertex missing from the second", two, missing, 1, 0),
        ("a vertex missing from both", missing, missing, 2, 0),
        ("an infinite tolerance", infinite, infinite, 2, 0),
        ("a negative tolerance", negative, base, 1, 0),
        ("a coordinate that is NaN", nan, base, 1, 0),
        ("a case named twice", twice, twice, 2, 0),
    ]
    failures = 0
    for label, first, second, problems, identical in cases:
        found, same = compare(first, second)
        if len(found) != problems or same != identical:
            print(f"self-test: {label}: {found}, {same} identical", file=sys.stderr)
            failures += 1
    if failures:
        return 1
    print(f"compare_g2 self-test: {len(cases)} cases pass")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("reports", nargs="*", type=Path)
    parser.add_argument("--summary", type=Path, help="append a Markdown summary to this file")
    parser.add_argument("--self-test", action="store_true", help="run the built-in cases and exit")
    args = parser.parse_args()
    if args.self_test:
        return self_test()
    if len(args.reports) != 2:
        parser.error("two reports are needed")
    first, second = (json.loads(path.read_text(encoding="utf-8")) for path in args.reports)
    problems, identical = compare(first, second)
    text = summary(args.reports[0].name, args.reports[1].name, len(first), problems, identical)
    print(text, end="")
    if args.summary:
        with args.summary.open("a", encoding="utf-8") as out:
            out.write(text)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
