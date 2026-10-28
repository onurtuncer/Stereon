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
4. The debug viewer's dependencies, GLFW and Dear ImGui, are linked only
   under viewer/ (ADR-0035), and nothing under src/ includes their headers,
   an OpenGL header or anything of viewer/ (ADR-0034: the kernel stays
   headless).
5. A copyleft dependency admitted as an optional adapter (ADR-0008), today
   the AD backend CppAD under EPL-2.0 (ADR-0045), is linked only by test
   executables, and nothing under src/ includes its headers.

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
# Copyleft licences ADR-0008 admits only for optional, separately built
# adapters, with the text their licence file must contain.
ADAPTER_ONLY = {
    "EPL-2.0 OR GPL-2.0-or-later": "Eclipse Public License",
}
# The adapter targets and where they may be named: test executables and the
# file that declares them.
ADAPTER_TARGETS = re.compile(r"(?<![\w-])(stereon_cppad|stereon::cppad)(?![\w-])")
ADAPTER_ALLOWED = ("tests/", DEPENDENCY_FILE.as_posix())
ADAPTER_INCLUDE = re.compile(r'^\s*#\s*include\s*[<"]cppad/', re.MULTILINE)

LICENCE_FILES = ("LICENSE", "LICENSE.txt", "LICENSE.md", "LICENCE", "COPYING", "LICENSE_1_0.txt")

# Where OpenCASCADE may be found: test and benchmark executables only.
OCCT_ALLOWED = ("tests/", "bench/")

# The viewer's dependency targets, and where they may be named: the viewer
# and the file that declares them.
VIEWER_TARGETS = re.compile(r"(?<![\w-])(glfw|stereon_imgui)(?![\w-])")
VIEWER_ALLOWED = ("viewer/", DEPENDENCY_FILE.as_posix())
# What no kernel source may include: a window system, a GUI toolkit, OpenGL,
# or the viewer.
VIEWER_INCLUDE = re.compile(
    r'^\s*#\s*include\s*[<"](?:GLFW/|imgui|GL/|OpenGL/|GLES\d*/|glad/|viewer/)', re.MULTILINE)
SOURCE_SUFFIXES = (".cppm", ".cpp", ".hpp", ".h", ".ixx")


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
        elif licence not in ALLOWED and licence not in ADAPTER_ONLY:
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
        expected = ALLOWED.get(licence) or ADAPTER_ONLY.get(licence)
        if expected is None or not source.is_dir():
            continue  # undeclared ones are reported already; unfetched ones are optional
        found = licence_text(source)
        if found is None:
            problems.append(f"{name}: no licence file in {source}")
        elif expected.lower() not in " ".join(found.split()).lower():
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


def cmake_files(root: Path) -> list[Path]:
    files = sorted(root.rglob("CMakeLists.txt")) + sorted(root.rglob("*.cmake"))
    return [p for p in files if not p.relative_to(root).as_posix().startswith(("build", "_deps"))
            and "/_deps/" not in p.relative_to(root).as_posix()]


def check_viewer_boundary(root: Path) -> list[str]:
    """The viewer's dependencies stay in viewer/, and src/ stays headless."""
    problems = []
    for path in cmake_files(root):
        relative = path.relative_to(root).as_posix()
        if relative.startswith(VIEWER_ALLOWED):
            continue
        for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            code = line.split("#", 1)[0]
            if VIEWER_TARGETS.search(code):
                problems.append(f"{relative}:{number}: names a debug viewer dependency outside "
                                "viewer/ (ADR-0035)")
    for path in sorted((root / "src").rglob("*")):
        if path.suffix in SOURCE_SUFFIXES and path.is_file():
            for match in VIEWER_INCLUDE.finditer(path.read_text(encoding="utf-8")):
                number = path.read_text(encoding="utf-8").count("\n", 0, match.start()) + 1
                problems.append(f"{path.relative_to(root).as_posix()}:{number}: includes a window "
                                "system, GUI or OpenGL header or the viewer (ADR-0034)")
    return problems


def check_adapter_boundary(root: Path) -> list[str]:
    """Copyleft adapters stay out of the libraries (ADR-0008, ADR-0045)."""
    problems = []
    for path in cmake_files(root):
        relative = path.relative_to(root).as_posix()
        if relative.startswith(ADAPTER_ALLOWED):
            continue
        for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            if ADAPTER_TARGETS.search(line.split("#", 1)[0]):
                problems.append(f"{relative}:{number}: links a copyleft adapter outside the test "
                                "executables (ADR-0008)")
    for path in sorted((root / "src").rglob("*")):
        if path.suffix in SOURCE_SUFFIXES and path.is_file():
            text = path.read_text(encoding="utf-8")
            for match in ADAPTER_INCLUDE.finditer(text):
                number = text.count("\n", 0, match.start()) + 1
                problems.append(f"{path.relative_to(root).as_posix()}:{number}: includes CppAD, "
                                "a copyleft adapter, in a library (ADR-0008, ADR-0045)")
    return problems


def self_test_adapter_boundary() -> None:
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        (root / "tests").mkdir()
        (root / "src/curve").mkdir(parents=True)
        (root / "tests/CMakeLists.txt").write_text(
            "target_link_libraries(t PRIVATE stereon::cppad)\n", encoding="utf-8")
        (root / "src/curve/CMakeLists.txt").write_text("# stereon::cppad\n", encoding="utf-8")
        (root / "src/curve/A.cpp").write_text("#include <vector>\n", encoding="utf-8")
        assert check_adapter_boundary(root) == [], check_adapter_boundary(root)
        (root / "src/curve/CMakeLists.txt").write_text(
            "target_link_libraries(c PRIVATE stereon_cppad)\n", encoding="utf-8")
        (root / "src/curve/A.cpp").write_text("#include <cppad/cppad.hpp>\n", encoding="utf-8")
        problems = check_adapter_boundary(root)
        assert len(problems) == 2, problems


def self_test_viewer_boundary() -> None:
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        (root / "viewer/app").mkdir(parents=True)
        (root / "src/core").mkdir(parents=True)
        (root / "viewer/app/CMakeLists.txt").write_text(
            "target_link_libraries(stn-view PRIVATE stereon_imgui glfw)\n", encoding="utf-8")
        (root / "src/core/CMakeLists.txt").write_text(
            "# glfw in a comment is fine\nadd_library(stereon_core)\n", encoding="utf-8")
        (root / "src/core/Core.cpp").write_text("#include <vector>\n", encoding="utf-8")
        assert check_viewer_boundary(root) == [], check_viewer_boundary(root)
        (root / "src/core/CMakeLists.txt").write_text(
            "target_link_libraries(stereon_core PRIVATE glfw)\n", encoding="utf-8")
        (root / "src/core/Core.cpp").write_text("#include <vector>\n#include <GLFW/glfw3.h>\n",
                                                encoding="utf-8")
        problems = check_viewer_boundary(root)
        assert len(problems) == 2, problems
        assert problems[0].startswith("src/core/CMakeLists.txt:1:"), problems
        assert problems[1].startswith("src/core/Core.cpp:2:"), problems


def self_test() -> int:
    self_test_viewer_boundary()
    self_test_adapter_boundary()
    adapter = "# Dependency licence: EPL-2.0 OR GPL-2.0-or-later\nFetchContent_Declare(cppad)\n"
    assert check_declarations(adapter) == [], check_declarations(adapter)
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
    problems = (check_declarations(text) + check_occt(args.root) + check_viewer_boundary(args.root)
                + check_adapter_boundary(args.root))
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
