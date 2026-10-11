#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""Enforce the function complexity limits of ADR-0022.

Every function, member function and function template is measured with
lizard, which reads the sources without compiling them: cyclomatic complexity
(decision points plus one) and non-comment lines. The limits follow the
ADR-0015 tiers:

    Tier C   src/intersect, src/boolean, src/blend        complexity <= 10, lines <= 60
    Tier B   src/robust, src/topo, src/build, src/mesh    complexity <= 12, lines <= 80
    Tier A   everything else (src/core, src/io, python, bench, ...)
                                                          complexity <= 15, lines <= 100
    tests    tests/                                       complexity <= 15, lines <= 150

Lambdas count with the function that contains them. The contract specifiers
STEREON_PRE, STEREON_POST and their HOT forms are blanked before measuring;
lizard would otherwise take the macro for the function's name.

A function that legitimately exceeds a limit carries, on the line before it,

    // complexity-exempt: <reason>

The reason is required, and the number of exemptions is reported on every
run so that a rising count is visible at phase gates.

This is a script rather than a clang-tidy check because clang-tidy cannot
parse the GCC module and reflection flags the build uses (ADR-0016).

Usage: check_complexity.py [--root <repository root>] | --self-test

Exit codes: 0 when every function is within its limits, 1 on violations,
77 when lizard is not installed (CTest reads 77 as "skipped").
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

SKIP_EXIT_CODE = 77
SOURCE_SUFFIXES = (".cppm", ".cpp", ".hpp", ".h")
EXEMPT_MARKER = re.compile(r"//\s*complexity-exempt:(?P<reason>.*)$")
CONTRACT_SPECIFIER = re.compile(r"\bSTEREON_(?:HOT_)?(?:PRE|POST)\s*\(")
# Lines that may sit between an exemption comment and the function it exempts.
PREAMBLE = re.compile(r"^\s*(?:template\b|requires\b|\[\[|$)")


@dataclass(frozen=True)
class Limits:
    tier: str
    complexity: int
    lines: int


TIER_C = Limits("Tier C", 10, 60)
TIER_B = Limits("Tier B", 12, 80)
TIER_A = Limits("Tier A", 15, 100)
TESTS = Limits("tests", 15, 150)

# Directory prefix (POSIX, relative to the root) -> limits. Longest match wins.
DIRECTORY_LIMITS = {
    "src/intersect/": TIER_C,
    "src/boolean/": TIER_C,
    "src/blend/": TIER_C,
    "src/robust/": TIER_B,
    "src/topo/": TIER_B,
    "src/build/": TIER_B,
    "src/mesh/": TIER_B,
    "tests/": TESTS,
}


@dataclass
class Function:
    name: str
    line: int
    complexity: int
    lines: int
    exempt: bool


@dataclass
class FileResult:
    path: str
    functions: list[Function]
    violations: list[tuple[int, str]]
    exemptions: int


def limits_for(relative_path: str) -> Limits:
    best = None
    for prefix, limits in DIRECTORY_LIMITS.items():
        if relative_path.startswith(prefix) and (best is None or len(prefix) > len(best[0])):
            best = (prefix, limits)
    return best[1] if best else TIER_A


def blank_contract_specifiers(code: str) -> str:
    """Replace each STEREON_PRE(...)/STEREON_POST(...) with spaces, keeping line breaks."""
    out = []
    position = 0
    for match in CONTRACT_SPECIFIER.finditer(code):
        if match.start() < position:
            continue
        depth = 1
        end = match.end()
        while end < len(code) and depth:
            if code[end] == "(":
                depth += 1
            elif code[end] == ")":
                depth -= 1
            end += 1
        out.append(code[position : match.start()])
        out.append(re.sub(r"[^\n]", " ", code[match.start() : end]))
        position = end
    out.append(code[position:])
    return "".join(out)


def exemption_before(lines: list[str], start_line: int) -> tuple[bool, str | None]:
    """(exempt, problem): the marker on the line before the function, skipping template
    and attribute lines; `problem` names a marker without a reason."""
    index = start_line - 2  # 0-based line before the function
    while index >= 0 and PREAMBLE.match(lines[index]) and not lines[index].strip().startswith("//"):
        index -= 1
    if index < 0:
        return False, None
    match = EXEMPT_MARKER.search(lines[index])
    if not match:
        return False, None
    if not match.group("reason").strip():
        return True, "complexity-exempt marker without a reason"
    return True, None


def analyze(relative_path: str, code: str) -> FileResult:
    import lizard  # noqa: PLC0415  (checked in main())

    limits = limits_for(relative_path)
    info = lizard.analyze_file.analyze_source_code(relative_path, blank_contract_specifiers(code))
    lines = code.splitlines()
    functions = []
    violations = []
    exemptions = 0
    for f in info.function_list:
        exempt, problem = exemption_before(lines, f.start_line)
        functions.append(Function(f.name, f.start_line, f.cyclomatic_complexity, f.nloc, exempt))
        if problem:
            violations.append((f.start_line, f"{f.name}: {problem} (ADR-0022)"))
            continue
        over = []
        if f.cyclomatic_complexity > limits.complexity:
            over.append(f"cyclomatic complexity {f.cyclomatic_complexity} exceeds the {limits.tier} limit of {limits.complexity}")
        if f.nloc > limits.lines:
            over.append(f"{f.nloc} non-comment lines exceed the {limits.tier} limit of {limits.lines}")
        if not over:
            continue
        if exempt:
            exemptions += 1
            continue
        for message in over:
            violations.append((f.start_line, f"{f.name}: {message} (ADR-0022); split along the case analysis"))
    return FileResult(relative_path, functions, violations, exemptions)


def tracked_sources(root: Path) -> list[Path]:
    """The tracked C++ sources, from git; every file under src/, tests/, bench/ and python/
    when the tree is not a repository."""
    try:
        listing = subprocess.run(
            ["git", "-C", str(root), "ls-files", "-z", "--", *[f"*{suffix}" for suffix in SOURCE_SUFFIXES]],
            check=True, capture_output=True, text=True, encoding="utf-8",
        ).stdout
        return [root / name for name in listing.split("\0") if name]
    except (OSError, subprocess.CalledProcessError):
        files = []
        for directory in ("src", "tests", "bench", "python"):
            if (root / directory).is_dir():
                files.extend(p for p in (root / directory).rglob("*") if p.suffix in SOURCE_SUFFIXES)
        return sorted(files)


def check_tree(root: Path) -> int:
    results = [analyze(path.relative_to(root).as_posix(), path.read_text(encoding="utf-8"))
               for path in tracked_sources(root) if path.is_file()]
    violations = 0
    for result in results:
        for line, message in result.violations:
            print(f"{result.path}:{line}: {message}", file=sys.stderr)
            violations += 1
    functions = [(f, r.path) for r in results for f in r.functions]
    exemptions = sum(r.exemptions for r in results)
    if violations:
        print(f"check_complexity: {violations} violation(s)", file=sys.stderr)
        return 1
    if functions:
        worst, where = max(functions, key=lambda item: item[0].complexity)
        peak = f"max complexity {worst.complexity} ({worst.name} in {where})"
    else:
        peak = "no functions"
    print(f"check_complexity: {len(functions)} functions in {len(results)} files, {peak}, "
          f"{exemptions} exemption(s), no violations")
    return 0


# Self-test material. BRANCHY has cyclomatic complexity 11: ten ifs plus one.
BRANCHY = "int Branchy(int n)\n{\n" + "".join(f"    if (n == {i}) {{\n        return {i};\n    }}\n" for i in range(10)) + "    return 0;\n}\n"
LAMBDA = """int Outer(int n)
{
    auto inner = [&](int k) {
        if (k > 1) {
            return k;
        } else if (k < -1) {
            return -k;
        }
        return 0;
    };
    for (int i = 0; i < n; ++i) {
        if (i % 2 == 0 && i % 3 == 0) {
            n += inner(i);
        }
    }
    return n;
}
"""
CONTRACTS = """export module M;
export struct Interval {
    constexpr Interval(double lower, double upper) STEREON_HOT_PRE(lower <= upper) : m_Lower(lower), m_Upper(upper) {}
    double m_Lower;
    double m_Upper;
};
export template <typename T>
constexpr T Clamp(T v, T lo, T hi) STEREON_PRE(lo <= hi) STEREON_POST(r : r >= lo) { return v < lo ? lo : (v > hi ? hi : v); }
"""
LONG = "int Long(int n)\n{\n" + "".join(f"    n += {i};\n" for i in range(120)) + "    return n;\n}\n"


def self_test() -> int:
    failures = []

    def expect(condition: bool, text: str) -> None:
        if not condition:
            failures.append(text)

    flagged = analyze("src/boolean/Branchy.cppm", BRANCHY)
    expect(len(flagged.violations) == 1 and "complexity 11" in flagged.violations[0][1],
           "complexity 11 under src/boolean/ is not flagged")
    expect(analyze("src/topo/Branchy.cppm", BRANCHY).violations == [], "complexity 11 under src/topo/ (Tier B, limit 12) is flagged")
    expect(analyze("src/core/Branchy.cppm", BRANCHY).violations == [], "complexity 11 under src/core/ is flagged")

    exempt = analyze("src/boolean/Branchy.cppm", "// complexity-exempt: reference implementation kept literal\n" + BRANCHY)
    expect(exempt.violations == [] and exempt.exemptions == 1, "an exempted function is flagged or not counted")
    exempt_template = analyze("src/boolean/Branchy.cppm", "// complexity-exempt: dispatch table\ntemplate <int N>\n" + BRANCHY)
    expect(exempt_template.violations == [] and exempt_template.exemptions == 1, "an exemption above a template line is not honoured")
    no_reason = analyze("src/boolean/Branchy.cppm", "// complexity-exempt:\n" + BRANCHY)
    expect(len(no_reason.violations) == 1 and "without a reason" in no_reason.violations[0][1], "a marker without a reason passes")
    unused = analyze("src/core/Branchy.cppm", "// complexity-exempt: not needed here\n" + BRANCHY)
    expect(unused.exemptions == 0, "an exemption on a function within its limits is counted")

    lam = analyze("src/core/Lambda.cppm", LAMBDA)
    expect(len(lam.functions) == 1 and lam.functions[0].complexity == 6,
           f"lambda branches are not counted with the enclosing function: {[(f.name, f.complexity) for f in lam.functions]}")

    names = [f.name for f in analyze("src/core/Contracts.cppm", CONTRACTS).functions]
    expect(any(n.startswith("Interval::Interval") for n in names) and any(n.startswith("Clamp") for n in names)
           and not any("STEREON" in n for n in names), f"contract specifiers are not blanked: {names}")

    expect(analyze("tests/core/LongTests.cpp", LONG).violations == [], "120 lines in tests/ (limit 150) are flagged")
    long_core = analyze("src/core/Long.cppm", LONG)
    expect(len(long_core.violations) == 1 and "non-comment lines" in long_core.violations[0][1], "120 lines under src/core/ (limit 100) pass")

    expect(limits_for("src/blend/Fillet.cppm") is TIER_C and limits_for("bench/X.cpp") is TIER_A
           and limits_for("python/StereonModule.cpp") is TIER_A, "directory limits are misassigned")

    for failure in failures:
        print("self-test: " + failure, file=sys.stderr)
    if failures:
        return 1
    print("check_complexity self-test: 12 cases pass")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument("--self-test", action="store_true", help="run the built-in cases and exit")
    args = parser.parse_args()
    try:
        import lizard  # noqa: F401
    except ImportError:
        print("check_complexity: lizard is not installed (pip install lizard); skipping", file=sys.stderr)
        return SKIP_EXIT_CODE
    if args.self_test:
        return self_test()
    return check_tree(args.root.resolve())


if __name__ == "__main__":
    sys.exit(main())
