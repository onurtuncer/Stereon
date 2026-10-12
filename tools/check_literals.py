#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""Ban magic numbers (ADR-0024).

Two patterns are rejected in the tracked C++ sources:

- a floating literal that spells a <numbers> constant, one of its common
  multiples (2 pi, pi/2, pi/4) or a degree conversion factor, anywhere:
  write std::numbers::pi (or pi_v<T>) instead;
- in src/ and bench/, a decimal floating literal with ten or more
  significant digits: a value that precise is a named constant that should
  carry its formula, e.g. (3.0 + 16.0 * Epsilon) * Epsilon rather than
  3.3306690738754706e-16. Tests are exempt, because expected values are
  literals by nature. Hexadecimal floats (0x1p-53) are always fine.

Comments and string literals are ignored. A justified exception is marked on
the same line:
    const double k = 2.7182818284; // STEREON-LINT-ALLOW(literal): reference value from the paper

This is a script rather than a clang-tidy check because clang-tidy cannot
parse the GCC module and reflection flags the build uses (ADR-0016).

Usage: check_literals.py [--root <repository root>] | --self-test
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

SOURCE_SUFFIXES = (".cppm", ".cpp", ".hpp", ".h", ".ixx")
CHECKED_DIRECTORIES = ("src", "tests", "bench", "python")
PRECISION_DIRECTORIES = ("src", "bench")
PRECISION_DIGITS = 10
ALLOW_MARKER = "STEREON-LINT-ALLOW(literal)"

# The leading significant digits of each constant, as they appear in a decimal
# literal. Five digits are enough to recognise the constant and few enough that
# a value such as 3.14 (which may be anything) is left alone.
KNOWN_CONSTANTS = {
    "pi": "3.1415",
    "2 pi": "6.2831",
    "pi / 2": "1.5707",
    "pi / 4": "0.7853",
    "1 / pi": "0.3183",
    "e": "2.7182",
    "sqrt2": "1.4142",
    "1 / sqrt2": "0.7071",
    "sqrt3": "1.7320",
    "1 / sqrt3": "0.5773",
    "ln2": "0.6931",
    "ln10": "2.3025",
    "log2e": "1.4426",
    "log10e": "0.4342",
    "egamma": "0.5772",
    "phi": "1.6180",
    "pi / 180": "0.0174532",
    "180 / pi": "57.2957",
}

# A decimal floating literal: has a decimal point or an exponent. Hexadecimal
# floats start with 0x and are not matched (the lookbehind excludes an 'x').
FLOAT = re.compile(
    r"(?<![\w.])(?P<body>\d[\d']*\.[\d']*(?:[eE][-+]?\d+)?|\.\d[\d']*(?:[eE][-+]?\d+)?|\d[\d']*[eE][-+]?\d+)[fFlL]?(?![\w.])"
)


def blank_comments_and_strings(text: str) -> str:
    """Replace comments, string and character literals with spaces, keeping line breaks."""

    def blank(match: re.Match[str]) -> str:
        return re.sub(r"[^\n]", " ", match.group(0))

    pattern = re.compile(
        r'//[^\n]*|/\*.*?\*/|R"(?P<delim>[^(\s]*)\(.*?\)(?P=delim)"|"(?:\\.|[^"\\\n])*"|\'(?:\\.|[^\'\\\n])*\'',
        re.DOTALL,
    )
    return pattern.sub(blank, text)


def significant_digits(literal: str) -> int:
    """Significant decimal digits of a literal's mantissa (leading zeros do not count)."""
    mantissa = re.split(r"[eE]", literal, maxsplit=1)[0].replace("'", "").replace(".", "")
    return len(mantissa.lstrip("0"))


def known_constant(literal: str) -> str | None:
    """The <numbers> constant the literal spells, if any."""
    digits = literal.replace("'", "")
    # Normalise "3.1415e0" style mantissas only by their digit prefix.
    mantissa = re.split(r"[eE]", digits, maxsplit=1)[0]
    if mantissa.startswith("."):
        mantissa = "0" + mantissa
    for name, prefix in KNOWN_CONSTANTS.items():
        if mantissa.startswith(prefix):
            return name
    return None


def check_text(text: str, check_precision: bool) -> list[tuple[int, str]]:
    """Return (line, message) for every violation in one file's text."""
    problems = []
    raw_lines = text.split("\n")
    for number, line in enumerate(blank_comments_and_strings(text).split("\n"), start=1):
        if ALLOW_MARKER in raw_lines[number - 1]:
            continue
        for match in FLOAT.finditer(line):
            literal = match.group("body")
            constant = known_constant(literal)
            if constant is not None:
                problems.append((number, f"{literal} spells {constant}: use std::numbers (ADR-0024)"))
            elif check_precision and significant_digits(literal) >= PRECISION_DIGITS:
                problems.append(
                    (number, f"{literal} has {significant_digits(literal)} significant digits: name it and write its formula (ADR-0024)")
                )
    return problems


def walked_sources(root: Path) -> list[Path]:
    """Every source file under the checked directories, in a fixed order."""
    files = []
    for directory in CHECKED_DIRECTORIES:
        base = root / directory
        if base.is_dir():
            files.extend(f for f in sorted(base.rglob("*")) if f.is_file() and f.suffix in SOURCE_SUFFIXES)
    return files


def tracked_sources(root: Path) -> list[Path]:
    """The tracked source files under the checked directories.

    Falls back to walking the tree when git cannot answer: CI runs the
    check inside a container where the checkout belongs to another user and
    git refuses the repository ("dubious ownership"), and an exported tree
    has no git at all.
    """
    try:
        result = subprocess.run(
            ["git", "-C", str(root), "ls-files", "-z", "--", *[f"{d}/" for d in CHECKED_DIRECTORIES]],
            check=True,
            capture_output=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return walked_sources(root)
    files = [root / name for name in result.stdout.decode("utf-8").split("\0") if name]
    return [f for f in files if f.suffix in SOURCE_SUFFIXES]


def check_tree(root: Path) -> int:
    problems = 0
    files = 0
    for path in tracked_sources(root):
        files += 1
        relative = path.relative_to(root).as_posix()
        check_precision = relative.split("/", 1)[0] in PRECISION_DIRECTORIES
        for line, message in check_text(path.read_text(encoding="utf-8"), check_precision):
            print(f"{relative}:{line}: {message}")
            problems += 1
    if problems:
        print(f"check_literals: {problems} violation(s) in {files} file(s)", file=sys.stderr)
        return 1
    print(f"check_literals: {files} file(s), no violations")
    return 0


SELF_TEST_BAD = [
    "const double angle = t * 3.141592653589793;",
    "constexpr double TwoPi = 6.283185307179586;",
    "return x * 0.017453292519943295;",
    "const auto half = 1.5707963f;",
    "double e = 2.718281828459045;",
    "double r = 1.41421356237;",
    "x = .70710678;",
    "constexpr double Bound = 3.3306690738754706e-16;",  # precision rule
    "const double k = 1.2345678901;",
]

SELF_TEST_GOOD = [
    "const double angle = t * std::numbers::pi;",
    "constexpr double Epsilon = 0x1p-53;",
    "constexpr double Splitter = 0x1p27 + 1.0;",
    "constexpr double Bound = (3.0 + 16.0 * Epsilon) * Epsilon;",
    "const double half = 0.5 * width;",
    "if (count > 3) return 1e-5;",
    "// 3.141592653589793 in a comment is fine",
    'const char* text = "pi is 3.14159";',
    "double k = 2.7182818284; // STEREON-LINT-ALLOW(literal): reference value from the paper",
    "x = 3.14;",  # too short to be sure
    "constexpr std::size_t BodyFaces = 1'000'000;",
]

SELF_TEST_TESTS_ONLY_GOOD = [
    "REQUIRE(value == 0.333333333333);",  # expected values in tests may be long
]


def self_test() -> int:
    failures = 0
    for line in SELF_TEST_BAD:
        if not check_text(line, check_precision=True):
            print(f"self-test: expected a violation: {line}")
            failures += 1
    for line in SELF_TEST_GOOD:
        if check_text(line, check_precision=True):
            print(f"self-test: unexpected violation: {line}")
            failures += 1
    for line in SELF_TEST_TESTS_ONLY_GOOD:
        if check_text(line, check_precision=False):
            print(f"self-test: unexpected violation in test code: {line}")
            failures += 1
        if not check_text(line, check_precision=True):
            print(f"self-test: expected a precision violation in library code: {line}")
            failures += 1
    if failures:
        return 1
    print(
        f"check_literals self-test: {len(SELF_TEST_BAD)} bad, {len(SELF_TEST_GOOD)} good and "
        f"{len(SELF_TEST_TESTS_ONLY_GOOD)} test-only cases pass"
    )
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        return self_test()
    return check_tree(args.root.resolve())


if __name__ == "__main__":
    sys.exit(main())
