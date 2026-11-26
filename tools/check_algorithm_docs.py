#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
"""Check that each tracked first-party source belongs to a Sphinx method page."""
import argparse
import ast
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
SUFFIXES = set(".c .cc .cpp .cppm .cxx .h .hh .hpp .hxx .inl .ipp .ixx .py .sh .ps1 .cmake .js .ts .rs .cu .cuh".split())


def source_file(path):
    """Recognise implementation, test and build source tracked by this repository."""
    return path.name == "CMakeLists.txt" or path.suffix.lower() in SUFFIXES


def check(root, base=None):
    """Validate file ownership, citations and newly introduced routine descriptions."""
    raw = subprocess.check_output(["git", "ls-files", "-z"], cwd=root)
    sources = {p for item in raw.decode().split("\0") if item and source_file(p := Path(item))}
    seen = {}
    errors = []
    keys = set(re.findall(r"@\w+\s*\{([^,]+),", (root / "docs/references.bib").read_text()))
    for page in sorted((root / "docs/algorithms").glob("*.rst")):
        content = page.read_text()
        for group in re.findall(r":cite:p:`([^`]+)`", content):
            errors.extend(f"{page.name}: missing citation {key}" for key in group.split(",") if key not in keys)
        for block in re.findall(r"^\.\. algorithm-sources::\n((?:   [^\n]+\n)+)", content, re.M):
            for line in block.splitlines():
                path = Path(line.strip())
                if path in seen:
                    errors.append(f"{path}: duplicate source entry in {page.name} and {seen[path]}")
                seen[path] = page.name
    errors.extend(f"{p}: missing Sphinx algorithm entry" for p in sorted(sources - seen.keys()))
    errors.extend(f"{p}: stale Sphinx algorithm entry" for p in sorted(seen.keys() - sources))
    if base:
        errors.extend(check_new_routines(root, base, seen))
    return len(sources), errors


def routines(path, code):
    """Index routine signatures and whether they have renderable implementation notes."""
    if path.suffix == ".py":
        records = {}
        def visit(node, prefix):
            """Preserve class and enclosing function names in the inventory key."""
            for child in ast.iter_child_nodes(node):
                if isinstance(child, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
                    key = prefix + (child.name,)
                    if not isinstance(child, ast.ClassDef):
                        records[key] = (child.name, child.lineno, bool(ast.get_docstring(child)))
                    visit(child, key)
                else:
                    visit(child, prefix)
        visit(ast.parse(code), ())
        return records
    if path.suffix not in {".c", ".cpp", ".cc", ".cxx", ".cppm", ".ixx", ".h", ".hh", ".hpp", ".hxx", ".inl", ".ipp"}:
        return {}
    import lizard
    from check_complexity import blank_contract_specifiers
    lines = code.splitlines()
    bodies = lizard.analyze_file.analyze_source_code("source.cpp", blank_contract_specifiers(code))
    records = {}
    for fn in bodies.function_list:
        index = fn.start_line - 2
        while index >= 0:
            line = lines[index].strip()
            if line.startswith("///") or re.search(r"[;{}]|^#|^//", line):
                break
            index -= 1
        documented = index >= 0 and lines[index].strip().startswith("///")
        # A repeated signature can occur in separate local scopes or preprocessor branches.
        signature = fn.long_name
        key = (signature, 0)
        while key in records:
            key = (signature, key[1] + 1)
        records[key] = (fn.name.split("::")[-1], fn.start_line, documented)
    return records


def check_new_routines(root, base, owners):
    """Compare named definitions with the PR base, including private functions in old files."""
    raw = subprocess.check_output(["git", "diff", "--name-only", "-z", "--diff-filter=AM", "--no-renames", base, "--"], cwd=root)
    for filename in filter(None, raw.decode().split("\0")):
        path = Path(filename)
        if path not in owners:
            continue
        previous = subprocess.run(["git", "show", f"{base}:{filename}"], cwd=root,
                                  capture_output=True, text=True)
        before = routines(path, previous.stdout if previous.returncode == 0 else "")
        current = routines(path, (root / path).read_text())
        page = (root / "docs/algorithms" / owners[path]).read_text().split(".. algorithm-sources::")[0]
        for key in sorted(current.keys() - before.keys()):
            name, line, documented = current[key]
            if not documented and f"``{name}``" not in page:
                yield f"{filename}:{line}: document {name} in its source or in {owners[path]}"


def existing_review_candidates(root):
    """List kernel definitions without an adjacent doc comment or exact manual name."""
    manual = "\n".join(p.read_text() for p in (root / "docs").rglob("*.rst")
                       if "_build" not in p.parts)
    names_in_manual = set(re.findall(r"``([^`\n]+)``", manual))
    for path in sorted((root / "src").rglob("*")):
        if not path.is_file() or not source_file(path):
            continue
        relative = path.relative_to(root)
        for name, line, documented in routines(relative, path.read_text()).values():
            if not documented and name not in names_in_manual:
                yield f"{relative}:{line}: {name}"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", help="revision against which to check new named routines")
    parser.add_argument("--audit-existing", action="store_true",
                        help="print kernel routines to inspect for existing algorithm coverage")
    args = parser.parse_args()
    if args.audit_existing:
        candidates = list(existing_review_candidates(ROOT))
        for candidate in candidates:
            print(candidate)
        print(f"Existing kernel definitions suggested for manual review: {len(candidates)}")
    count, findings = check(ROOT, args.base)
    for finding in findings:
        print(finding)
    print(f"Algorithm documentation inventory: {count} source files, {len(findings)} findings")
    sys.exit(bool(findings))
