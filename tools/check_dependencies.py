#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""Check every third-party dependency's licence against ADR-0008.

ADR-0008 allows dependencies under MPL-2.0-compatible licences only (MIT,
BSD, Apache-2.0, MPL-2.0, Boost) and copyleft code only outside the
libraries. This script checks three things:

1. Every FetchContent_Declare in cmake/StereonDependencies.cmake has a
   `# Dependency licence: <SPDX expression>` line in the comments above it,
   and the expression is on the allow-list below.
2. With --deps-dir (the build's _deps directory, as CTest passes), every
   fetched dependency's licence file reads as the licence declared for it,
   so a wrong declaration is caught, not only a missing one.
3. OpenCASCADE (LGPL-2.1 with the OCCT exception) is found only by the test
   and benchmark CMake files, never by a library (ADR-0008: it is linked
   into those executables only).

Usage: check_dependencies.py [--root <repo>] [--deps-dir <build>/_deps] | --self-test
"""

from __future__ import annotations

import argparse
import re
import sys
import tempfile
from pathlib import Path

DEPENDENCY_FILE = Path("cmake/StereonDependencies.cmake")
MARKER = "Dependency licence:"

# MPL-2.0-compatible licences (ADR-0008), as SPDX expressions, with a phrase
# that the licence text itself contains.
ALLOWED = {
    "MIT": "Permission is hereby granted, free of charge",
    "BSD-2-Clause": "Redistributions in binary form must reproduce",
    "BSD-3-Clause": "Neither the name",
    "BSL-1.0": "Boost Software License",
    "Apache-2.0": "Apache License",
    "Apache-2.0 WITH LLVM-exception": "LLVM Exceptions to the Apache 2.0 License",
    "MPL-2.0": "Mozilla Public License",
    "Zlib": "This software is provided 'as-is'",
    "ISC": "Permission to use, copy, modify, and/or distribute",
}

# Licence files a dependency may ship, most specific first.
LICENCE_FILES = ("LICENSE", "LICENSE.txt", "LICENSE.md", "LICENCE", "COPYING", "LICENSE_1_0.txt")

# Where OpenCASCADE may be found: test and benchmark executables only.
OCCT_ALLOWED = ("tests/", "bench/")


def declarations(text: str) -> list[tuple[str, str | None, int]]:
    """(name, declared licence or None, line) for every FetchContent_Declare.

    The licence comes from the nearest `# Dependency licence:` line above the
    declaration and below the previous one."""
    lines = text.splitlines()
    found = []
    previous = 0
    for number, line in enumerate(lines):
        match = re.match(r"\s*FetchContent_Declare\(\s*(\w+)", line)
        if not match:
            continue
        licence = None
        for above in lines[previous:number]:
            comment = re.match(r"\s*#\s*" + re.escape(MARKER) + r"\s*(.+?)\s*$", above)
            if comment:
                licence = comment.group(1)
        found.append((match.group(1), licence, number + 1))
        previous = number + 1
    return found


def check_declarations(text: str) -> list[str]:
    problems = []
    for name, licence, line in declarations(text):
        if licence is None:
            problems.append(f"{DEPENDENCY_FILE.as_posix()}:{line}: {name} has no '# {MARKER} <SPDX>' line")
        elif licence not in ALLOWED:
            problems.append(f"{DEPENDENCY_FILE.as_posix()}:{line}: {name} is under '{licence}', which is "
                            "not on ADR-0008's allow-list")
    return problems


def licence_text(source: Path) -> str | None:
    for name in LICENCE_FILES:
        path = source / name
        if path.is_file():
            return path.read_text(encoding="utf-8", errors="replace")
    return None


def check_fetched(text: str, deps_dir: Path) -> list[str]:
    """Each fetched dependency's licence file matches its declaration."""
    problems = []
    for name, licence, _ in declarations(text):
        source = deps_dir / f"{name.lower()}-src"
        if licence not in ALLOWED or not source.is_dir():
            continue  # undeclared ones are reported already; unfetched ones are optional
        found = licence_text(source)
        if found is None:
            problems.append(f"{name}: no licence file in {source}")
        elif ALLOWED[licence].lower() not in " ".join(found.split()).lower():
            problems.append(f"{name}: the licence file in {source} does not read as {licence}")
    return problems


def check_occt(root: Path) -> list[str]:
    problems = []
    for path in sorted(root.rglob("CMakeLists.txt")) + sorted(root.rglob("*.cmake")):
        relative = path.relative_to(root).as_posix()
        if relative.startswith(("build", "_deps")) or "/_deps/" in relative:
            continue
        if re.search(r"find_package\(\s*OpenCASCADE", path.read_text(encoding="utf-8")):
            if not relative.startswith(OCCT_ALLOWED):
                problems.append(f"{relative}: finds OpenCASCADE (LGPL-2.1) outside the test and "
                                "benchmark executables")
    return problems


def self_test() -> int:
    good = ("# Dependency licence: MIT\nFetchContent_Declare(alpha\n    GIT_TAG v1)\n"
            "# Something else\n# Dependency licence: Apache-2.0 WITH LLVM-exception\n"
            "set(X ON)\nFetchContent_Declare(beta\n    GIT_TAG v2)\n")
    assert check_declarations(good) == [], check_declarations(good)
    missing = "# MIT\nFetchContent_Declare(gamma)\n"
    assert "no '# Dependency licence" in check_declarations(missing)[0]
    copyleft = "# Dependency licence: GPL-3.0-only\nFetchContent_Declare(delta)\n"
    assert "not on ADR-0008's allow-list" in check_declarations(copyleft)[0]
    # A declaration's licence line does not carry over to the next one.
    carried = "# Dependency licence: MIT\nFetchContent_Declare(a)\nFetchContent_Declare(b)\n"
    assert len(check_declarations(carried)) == 1
    # A fetched licence that is not the declared one is caught.
    with tempfile.TemporaryDirectory() as deps:
        source = Path(deps) / "alpha-src"
        source.mkdir()
        (source / "LICENSE").write_text("GNU GENERAL PUBLIC LICENSE Version 3", encoding="utf-8")
        mismatch = check_fetched("# Dependency licence: MIT\nFetchContent_Declare(alpha)\n",
                                 Path(deps))
        assert mismatch and "does not read as MIT" in mismatch[0], mismatch
        (source / "LICENSE").write_text("Permission is hereby granted, free of charge, to any",
                                        encoding="utf-8")
        assert check_fetched("# Dependency licence: MIT\nFetchContent_Declare(alpha)\n",
                             Path(deps)) == []
    print("check_dependencies: self-test passed")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument("--deps-dir", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        return self_test()
    text = (args.root / DEPENDENCY_FILE).read_text(encoding="utf-8")
    problems = check_declarations(text) + check_occt(args.root)
    if args.deps_dir:
        problems += check_fetched(text, args.deps_dir)
    for problem in problems:
        print(problem)
    count = len(declarations(text))
    checked = " and their fetched licence texts" if args.deps_dir else ""
    print(f"check_dependencies: {count} dependencies{checked}, "
          f"{'no violations' if not problems else f'{len(problems)} violation(s)'}")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
