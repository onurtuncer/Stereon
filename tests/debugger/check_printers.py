#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""Check Stereon's debugger pretty-printers in a scripted session (ADR-0004).

Runs GDB or LLDB in batch mode on tests/debugger/DebuggerTarget.cpp, stops
in StereonDebuggerStop(), prints every variable of main() and compares each
line with the expected summary below. Both debuggers must produce the same
summaries, so the GDB and LLDB printers cannot drift apart.

Usage: check_printers.py --debugger gdb|lldb --target <exe> --printers <script> [--checked-handles]

Exit codes: 0 when every summary matches, 1 on a mismatch, 77 when the
debugger is not installed (CTest reads 77 as "skipped").
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

SKIP_EXIT_CODE = 77

# Variable -> the summary its printer must produce (a prefix where the rest
# depends on the debugger's own layout).
EXPECTED = {
    "vector": "(1, 2.5, -3)",
    "point": "(4, 5, 6)",
    "uv": "(0.5, -1.5)",
    "interval": "[1, 2]",
    "emptyInterval": "[]",
    "box": "[(0, 0, 0), (1, 2, 3)]",
    "emptyBox": "[]",
    "vertex": "Vertex#12",
    "invalidEdge": "Edge#-",
    "coedge": "Coedge#42 (reversed)",
    "face": "Face#0@0",
    "lineCurve": "Curve(Line) {Origin = (0, 0, 1), Direction = (1, 0, 0)}",
    "circleCurve": "Curve(Circle) {Center = (1, 2, 3), Radius = 2.5, XAxis = (1, 0, 0), YAxis = (0, 1, 0)}",
    "userCurve": "Curve((anonymous namespace)::UserCurve) {Base = {Origin = (0, 0, 1), Direction = (1, 0, 0)}}",
    "emptyCurve": "Curve (empty)",
    "planarCurve": "Curve2d(Line2d) {Origin = (0, 0), Direction = (0, 1)}",
    "cylinderSurface": "Surface(Cylinder) {Origin = (0, 0, 0), Radius = 1.5, XAxis = (1, 0, 0), YAxis = (0, 1, 0), RightHanded = false}",
    "ellipseCurve": "Curve(Ellipse) {Center = (0, 0, 0), MajorRadius = 3, MinorRadius = 1, XAxis = (1, 0, 0), YAxis = (0, 1, 0)}",
    "torusSurface": "Surface(Torus) {Center = (0, 0, 0), MajorRadius = 4, MinorRadius = 0.5, XAxis = (1, 0, 0), YAxis = (0, 1, 0), RightHanded = true}",
    "arena": "PagedArena of 6 entries in 2 pages, 1 shared",
    "snapshot": "PagedArena of 6 entries in 2 pages, 1 shared",
    "journal": "J#7",
    "sign": "Negative",
}

# Handles of a live shape. In a build with checked handles the printers ask
# the process to describe the entity (ADR-0004, safeguard 4); otherwise they
# show the plain handle.
SHAPE_HANDLES_CHECKED = {
    "wireVertex": "Vertex#1 [(1, 0, 0), tol 1e-07, 1 edges]",
    "wireEdge": "Edge#0 [Spatial Line, tol 1e-07, faces none]",
}
SHAPE_HANDLES_PLAIN = {"wireVertex": "Vertex#1", "wireEdge": "Edge#0"}


def native_path(path: Path) -> str:
    """The path as the debugger takes it. Under MSYS2 a native Windows
    debugger does not understand /d/... paths, so they become D:/...."""
    text = path.as_posix()
    match = re.match(r"^/([a-zA-Z])/(.*)$", text)
    if match and sys.platform in ("msys", "cygwin", "win32"):
        return f"{match.group(1).upper()}:/{match.group(2)}"
    return text


def gdb_session(target: Path, printers: Path, names: list[str]) -> list[str]:
    commands = [
        "set pagination off",
        "set print pretty off",
        f"source {native_path(printers)}",
        "break StereonDebuggerStop",
        "run",
        "up",
    ] + [f'printf "@@{name}@@ "\nprint {name}' for name in names]
    with tempfile.NamedTemporaryFile("w", suffix=".gdb", delete=False) as script:
        script.write("\n".join(commands) + "\n")
    result = subprocess.run(["gdb", "-batch", "-nx", "-x", script.name, str(target)],
                            capture_output=True, text=True, timeout=300)
    return (result.stdout + result.stderr).splitlines()


def lldb_session(target: Path, printers: Path, names: list[str]) -> list[str]:
    commands = [
        f"command script import {native_path(printers)}",
        "breakpoint set --name StereonDebuggerStop",
        # Disabling ASLR needs the personality syscall, which container
        # security profiles block ("personality set failed").
        "settings set target.disable-aslr false",
        "process launch",
        "up",
    ] + [f'script print("@@{name}@@ ", end="")\nframe variable {name}' for name in names]
    with tempfile.NamedTemporaryFile("w", suffix=".lldb", delete=False) as script:
        script.write("\n".join(commands) + "\nquit\n")
    result = subprocess.run(["lldb", "--batch", "--no-lldbinit", "-s", script.name, str(target)],
                            capture_output=True, text=True, timeout=300)
    return (result.stdout + result.stderr).splitlines()


def summaries(lines: list[str]) -> dict[str, str]:
    """name -> the printed value after its marker.

    GDB prints "$N = value" on the marker's line; LLDB echoes the command
    there and prints "(type) name = value" on a following line. The search
    runs from the marker to the next one.
    """
    marker = re.compile(r"@@(\w+)@@ ")
    starts = [(index, match) for index, line in enumerate(lines) if (match := marker.search(line))]
    found = {}
    for position, (index, match) in enumerate(starts):
        name = match.group(1)
        end = starts[position + 1][0] if position + 1 < len(starts) else len(lines)
        value = re.compile(r"(?:\$\d+|\b" + re.escape(name) + r") = (.*)$")
        candidates = [lines[index][match.end():]] + lines[index + 1:end]
        for text in candidates:
            hit = value.search(text)
            if hit:
                found[name] = hit.group(1).strip()
                break
    return found


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--debugger", choices=["gdb", "lldb"], required=True)
    parser.add_argument("--target", type=Path, required=True)
    parser.add_argument("--printers", type=Path, required=True)
    parser.add_argument("--checked-handles", action="store_true",
                        help="the target was built with STEREON_CHECKED_HANDLES")
    args = parser.parse_args()
    expected = EXPECTED | (SHAPE_HANDLES_CHECKED if args.checked_handles else SHAPE_HANDLES_PLAIN)

    if shutil.which(args.debugger) is None:
        print(f"check_printers: {args.debugger} not found; skipping")
        return SKIP_EXIT_CODE

    session = gdb_session if args.debugger == "gdb" else lldb_session
    lines = session(args.target.resolve(), args.printers.resolve(), list(expected))
    found = summaries(lines)

    failures = 0
    for name, summary in expected.items():
        actual = found.get(name)
        if actual is None or not actual.startswith(summary):
            print(f"{name}: expected {summary!r}, got {actual!r}")
            failures += 1
    if failures:
        print(f"check_printers ({args.debugger}): {failures} of {len(expected)} summaries wrong")
        print("--- session output ---")
        print("\n".join(lines[-80:]))
        return 1
    print(f"check_printers ({args.debugger}): {len(expected)} summaries match")
    return 0


if __name__ == "__main__":
    sys.exit(main())
