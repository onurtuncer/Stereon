#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""Fail on an upward dependency between Stereon libraries.

The layering is the order of add_subdirectory() calls in src/CMakeLists.txt:
each library may depend only on libraries listed before it. A library
depends on another when any file under src/<library>/ does one of:

- imports its module:           import Stereon.Topo;   export import Stereon.Topo:Part;
- includes its public headers:  #include <Stereon/Topo/...>
- names its CMake target:       stereon::topo, stereon_topo

Usage: check_layering.py [--root <repository root>]
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

SOURCE_SUFFIXES = {".cppm", ".cpp", ".hpp", ".h", ".ixx"}

IMPORT_RE = re.compile(r"^\s*(?:export\s+)?import\s+Stereon\.(\w+)", re.MULTILINE)
INCLUDE_RE = re.compile(r'^\s*#\s*include\s*[<"]Stereon/(\w+)/', re.MULTILINE)
TARGET_RE = re.compile(r"\bstereon(?:::|_)(\w+)\b")
SUBDIR_RE = re.compile(r"^\s*add_subdirectory\(\s*(\w+)\s*\)", re.MULTILINE)


def strip_comments(text: str, cmake: bool) -> str:
    if cmake:
        return re.sub(r"#[^\n]*", "", text)
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.DOTALL)
    return re.sub(r"//[^\n]*", "", text)


def layering(root: Path) -> list[str]:
    text = strip_comments((root / "src" / "CMakeLists.txt").read_text(encoding="utf-8"), cmake=True)
    return SUBDIR_RE.findall(text)


def dependencies(path: Path, libraries: set[str]) -> list[tuple[int, str, str]]:
    """Return (line, kind, library) for every reference to a Stereon library."""
    is_cmake = path.name == "CMakeLists.txt" or path.suffix == ".cmake"
    raw = path.read_text(encoding="utf-8")
    text = strip_comments(raw, cmake=is_cmake)
    found = []
    patterns = [("target", TARGET_RE)] if is_cmake else [("import", IMPORT_RE), ("include", INCLUDE_RE)]
    for kind, pattern in patterns:
        for match in pattern.finditer(text):
            library = match.group(1).lower()
            if library in libraries:
                line = text.count("\n", 0, match.start()) + 1
                found.append((line, kind, library))
    return found


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent.parent)
    args = parser.parse_args()
    root: Path = args.root.resolve()

    order = layering(root)
    if not order:
        print("check_layering: no libraries found in src/CMakeLists.txt", file=sys.stderr)
        return 2
    rank = {library: index for index, library in enumerate(order)}

    errors = []
    for library in order:
        directory = root / "src" / library
        if not directory.is_dir():
            errors.append(f"src/CMakeLists.txt lists '{library}', but src/{library}/ does not exist")
            continue
        for path in sorted(directory.rglob("*")):
            if not path.is_file():
                continue
            if path.suffix not in SOURCE_SUFFIXES and path.name != "CMakeLists.txt" and path.suffix != ".cmake":
                continue
            for line, kind, used in dependencies(path, set(order)):
                if used != library and rank[used] > rank[library]:
                    relative = path.relative_to(root).as_posix()
                    errors.append(f"{relative}:{line}: {library} must not depend on {used} ({kind}); "
                                  f"{used} is above it in the layering")

    for directory in sorted((root / "src").iterdir()):
        if directory.is_dir() and directory.name not in rank:
            errors.append(f"src/{directory.name}/ is not listed in src/CMakeLists.txt")

    if errors:
        print("Library layering violations (order: " + " -> ".join(order) + "):", file=sys.stderr)
        for error in errors:
            print("  " + error, file=sys.stderr)
        return 1
    print(f"check_layering: {len(order)} libraries, no upward dependencies")
    return 0


if __name__ == "__main__":
    sys.exit(main())
