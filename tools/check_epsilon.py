#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""Ban hidden epsilons and raw sign tests in topology code (ADR-0001, ADR-0002).

In src/topo/, src/boolean/ and src/mesh/, tolerances come from
ToleranceOf(entity) or ctx.GetResolution(), and topological decisions go
through predicates in Stereon.Robust. This check rejects:

- comparisons with a floating-point literal, e.g. `d < 1e-9`, `x > 0.0`,
  `0.5 <= t`, `a == 0.0`
- machine-epsilon constants: std::numeric_limits<T>::epsilon(),
  DBL_EPSILON, FLT_EPSILON, LDBL_EPSILON

It works on the token level, so it cannot see types: a comparison of a
double with an integer literal (`x > 0`) is not caught. Review still has to
look for those. Comments and string literals are ignored.

A justified exception is marked on the same line:
    if (t < 0.5) { // STEREON-LINT-ALLOW(epsilon): splits the span in half, not a tolerance

This is a script rather than a clang-tidy check because clang-tidy cannot
parse the GCC module and reflection flags the build uses (ADR-0016).

Usage: check_epsilon.py [--root <repository root>] | --self-test
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

CHECKED_LIBRARIES = ("topo", "boolean", "mesh")
SOURCE_SUFFIXES = {".cppm", ".cpp", ".hpp", ".h", ".ixx"}
ALLOW_MARKER = "STEREON-LINT-ALLOW(epsilon)"

# A floating literal: has a decimal point or an exponent (integers are not matched).
FLOAT = r"(?:\d[\d']*\.[\d']*(?:[eE][-+]?\d+)?|\.\d[\d']*(?:[eE][-+]?\d+)?|\d[\d']*[eE][-+]?\d+)[fFlL]?"
# <, >, <=, >=, ==, != but not <<, >>, ->, <=>, or template angle brackets next to other operators.
COMPARE = r"(?:(?<![<>=!\-])(?:<=|>=|==|!=)(?![=>])|(?<![<>\-])<(?![<=])|(?<![<>\-=])>(?![>=]))"
LITERAL_RIGHT = re.compile(COMPARE + r"\s*[-+]?\s*" + FLOAT + r"(?![\w.])")
LITERAL_LEFT = re.compile(r"(?<![\w.])" + FLOAT + r"\s*" + COMPARE)
EPSILON_CONSTANT = re.compile(r"\bnumeric_limits\s*<[^;]*?>\s*::\s*epsilon\b|\b(?:L?DBL|FLT|LDBL)_EPSILON\b")


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
    code_lines = blank_comments_and_strings(text).splitlines()
    problems = []
    for number, (code, original) in enumerate(zip(code_lines, original_lines), start=1):
        if ALLOW_MARKER in original:
            continue
        if LITERAL_RIGHT.search(code) or LITERAL_LEFT.search(code):
            problems.append((number, "comparison with a floating-point literal; use a Stereon.Robust "
                                     "predicate or ToleranceOf()/GetResolution() (ADR-0001, ADR-0002)"))
        if EPSILON_CONSTANT.search(code):
            problems.append((number, "machine epsilon is not a geometric tolerance; use ToleranceOf() or "
                                     "ctx.GetResolution() (ADR-0001)"))
    return problems


def check_tree(root: Path) -> int:
    count = 0
    files = 0
    for library in CHECKED_LIBRARIES:
        directory = root / "src" / library
        if not directory.is_dir():
            continue
        for path in sorted(directory.rglob("*")):
            if not path.is_file() or path.suffix not in SOURCE_SUFFIXES:
                continue
            files += 1
            for line, message in check_text(path.read_text(encoding="utf-8")):
                print(f"{path.relative_to(root).as_posix()}:{line}: {message}", file=sys.stderr)
                count += 1
    if count:
        print(f"check_epsilon: {count} violation(s)", file=sys.stderr)
        return 1
    print(f"check_epsilon: {files} file(s) in {', '.join(CHECKED_LIBRARIES)}, no violations")
    return 0


SELF_TEST_BAD = [
    "if (distance < 1e-9) {",
    "if (d <= 1.0e-6) return;",
    "return x > 0.0;",
    "bool zero = (area == 0.);",
    "if (0.5 <= t) {}",
    "if (std::abs(a - b) < .001f) {}",
    "auto r = s >= -1e-12;",
    "if (x != 0.0L) {}",
    "const double e = std::numeric_limits<double>::epsilon();",
    "if (d < DBL_EPSILON) {}",
    "if (d > 1'000.5) {}",
]
SELF_TEST_GOOD = [
    "if (distance < ToleranceOf(edge)) {",
    "if (Orient3D(a, b, c, d) > 0) {}",
    "std::cout << 0.5 << '\\n';",
    "auto shifted = value >> 2;",
    "auto x = p->Weight * 0.5;",
    "auto order = a <=> b;",
    "std::array<double, 3> v{0.5, 1.0, 2.0};",
    "double half = 0.5;",
    "// if (d < 1e-9) in a comment",
    'Log("tolerance < 1e-9");',
    "if (t < 0.5) { // STEREON-LINT-ALLOW(epsilon): span midpoint, not a tolerance",
    "for (int i = 0; i < 3; ++i) {}",
    "const auto eps = ctx.GetResolution();",
    "x = y >= 3 ? 1.5 : 2.5;",
]


def self_test() -> int:
    failures = []
    for line in SELF_TEST_BAD:
        if not check_text(line):
            failures.append(f"not flagged: {line}")
    for line in SELF_TEST_GOOD:
        problems = check_text(line)
        if problems:
            failures.append(f"wrongly flagged: {line} ({problems[0][1][:40]}...)")
    for failure in failures:
        print("self-test: " + failure, file=sys.stderr)
    if failures:
        return 1
    print(f"check_epsilon self-test: {len(SELF_TEST_BAD)} bad and {len(SELF_TEST_GOOD)} good cases pass")
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
