#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""Reject tolerance writes outside ShapeBuilder::SetTolerance (ADR-0050).

A vertex's or edge's tolerance changes only through SetTolerance(), which
records the change with its cause. The records live in the ShapeBuilder's
tables, which only src/topo/ can reach: other libraries see a Shape's
records read-only and change tolerances through the builder's public API.
So this check reads src/topo/ and rejects every assignment to a member
named Tolerance (`x.Tolerance = ...`, `p->Tolerance += ...`), except:

- on an EdgeSpec declared in the same file, which is a request to create
  an edge, not a write to a record;
- on a line marked as a reviewed exception: SetTolerance's own write, and
  the records AddVertex and AddEdge create:
      record.Tolerance = tolerance; // STEREON-LINT-ALLOW(tolerance-write): a new vertex

It works on the token level, so it cannot see a write through a reference
to the field (`double& t = r.Tolerance; t = 1.0;`). Review still has to look
for those. Comments and string literals are ignored.

Usage: check_tolerance_writes.py [--root <repository root>] | --self-test
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

CHECKED_DIRECTORY = Path("src") / "topo"
SOURCE_SUFFIXES = {".cppm", ".cpp", ".hpp", ".h", ".ixx"}
ALLOW_MARKER = "STEREON-LINT-ALLOW(tolerance-write)"

# `name.Tolerance =`, `name->Tolerance +=`, ...; not `==`.
WRITE = re.compile(r"(?:\b(\w+)\s*(?:\.|->)|(?:\.|->))\s*Tolerance\s*(?:[-+*/]?=)(?!=)")
EDGE_SPEC = re.compile(r"\bEdgeSpec\s+(\w+)\b")


def blank_comments_and_strings(text: str) -> str:
    """Replace comments, string and character literals with spaces, keeping line breaks."""

    def blank(match: re.Match[str]) -> str:
        return re.sub(r"[^\n]", " ", match.group(0))

    pattern = re.compile(
        r'//[^\n]*|/\*.*?\*/|R"(?P<delim>[^(\s]*)\(.*?\)(?P=delim)"|"(?:\\.|[^"\\\n])*"|\'(?:\\.|[^\'\\\n])*\'',
        re.DOTALL,
    )
    return pattern.sub(blank, text)


def check_text(text: str) -> list[tuple[int, str]]:
    """Return (line number, message) for each violation in a source text."""
    original_lines = text.splitlines()
    code = blank_comments_and_strings(text)
    specs = set(EDGE_SPEC.findall(code))
    problems = []
    for number, (line, original) in enumerate(zip(code.splitlines(), original_lines), start=1):
        if ALLOW_MARKER in original:
            continue
        for match in WRITE.finditer(line):
            if match.group(1) in specs:
                continue
            problems.append((number, "a tolerance written outside ShapeBuilder::SetTolerance(); "
                                     "change it through SetTolerance() with its cause (ADR-0050)"))
    return problems


def check_tree(root: Path) -> int:
    count = 0
    files = 0
    directory = root / CHECKED_DIRECTORY
    for path in sorted(directory.rglob("*")):
        if not path.is_file() or path.suffix not in SOURCE_SUFFIXES:
            continue
        files += 1
        for line, message in check_text(path.read_text(encoding="utf-8")):
            print(f"{path.relative_to(root).as_posix()}:{line}: {message}", file=sys.stderr)
            count += 1
    if count:
        print(f"check_tolerance_writes: {count} violation(s)", file=sys.stderr)
        return 1
    print(f"check_tolerance_writes: {files} file(s) in {CHECKED_DIRECTORY.as_posix()}, no violations")
    return 0


SELF_TEST_BAD = [
    "record.Tolerance = tolerance;",
    "r.Tolerance = 2.0 * r.Tolerance;",
    "m_Tables.Vertices[v].Tolerance = t;",
    "edge->Tolerance += gap;",
    "e.Tolerance *= 2.0;",
    "Modify(i, [&](EdgeRecord& e) { e.Tolerance = t; });",
]
SELF_TEST_GOOD = [
    "EdgeSpec spec;\nspec.Tolerance = original.Tolerance;",
    "const double t = record.Tolerance;",
    "if (r.Tolerance == t) {}",
    "if (r.Tolerance <= t) {}",
    "if (r.Tolerance >= t) {}",
    "if (r.Tolerance != t) {}",
    "double Tolerance = 0.0;",
    "// r.Tolerance = t; in a comment",
    'Fail("r.Tolerance = t");',
    "r.Tolerance = t; // STEREON-LINT-ALLOW(tolerance-write): SetTolerance itself",
    "STEREON_TRY_VOID(SetTolerance(kept, t, ToleranceCause::MergeVertices, \"merge\"));",
]


def self_test() -> int:
    failures = []
    for text in SELF_TEST_BAD:
        if not check_text(text):
            failures.append(f"not flagged: {text}")
    for text in SELF_TEST_GOOD:
        problems = check_text(text)
        if problems:
            failures.append(f"wrongly flagged: {text} ({problems[0][1][:40]}...)")
    for failure in failures:
        print("self-test: " + failure, file=sys.stderr)
    if failures:
        return 1
    print(f"check_tolerance_writes self-test: {len(SELF_TEST_BAD)} bad and {len(SELF_TEST_GOOD)} good cases pass")
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
