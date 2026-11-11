#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
"""Validate ADR-0053 inventory; --release refuses incomplete evidence.

Development CI validates records and retains base-branch obligations, without
claiming historical backfill is complete. Release mode additionally requires
reviewed inventory completeness, every defect closed with evidence, and all
MR-01/02/03/04 cells of every built-in ordered pair qualified.
Evidence references identify files or HTTPS records; their scientific validity
and external CI outcomes still require review, not a string-based assertion.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import subprocess
import sys

LEDGER = "docs/planning/intersection-metamorphic.json"
RELATIONS = {f"MR-{n:02}" for n in range(1, 5)}
STATES = {"untriaged", "fixture captured", "MR implemented", "verified", "exception reviewed", "not a defect"}


def handler_ids(root):
    def types(path, alias):
        text = (root / path).read_text()
        match = re.search(r"using\s+" + alias + r"\s*=\s*TypeList<([^>]+)>", text)
        if not match:
            raise ValueError(f"cannot read built-in catalogue {alias}")
        return [name.strip() for name in match[1].split(",")]
    curves = types("src/curve/Visit.cppm", "BuiltinCurveTypes")
    planar = types("src/curve/Visit.cppm", "BuiltinCurve2dTypes")
    surfaces = types("src/surface/Visit.cppm", "BuiltinSurfaceTypes")
    families = {"curve-curve": (curves, curves), "planar": (planar, planar),
                "curve-surface": (curves, surfaces), "surface-surface": (surfaces, surfaces)}
    return {f"{family}:{a}:{b}" for family, (left, right) in families.items()
            for a in left for b in right}


def reference(root, value):
    if not isinstance(value, str) or not value.strip():
        return False
    if value.startswith("https://"):
        return bool(re.fullmatch(r"https://[^/\s]+/\S+", value))
    path, _, anchor = value.partition("#")
    target = (root / path).resolve()
    if not target.is_relative_to(root.resolve()) or not target.is_file():
        return False
    return not anchor or anchor in target.read_text(encoding="utf-8")


def evidence(root, row, names, errors, label):
    for name in names:
        if not reference(root, row.get(name)):
            errors.append(f"{label}: missing/invalid evidence {name}")


def check(root, data, release=False, base=None, expected=None):
    errors = []
    if data.get("version") != 1:
        errors.append("unsupported ledger version")
    records = data["defects"]
    ids = [row["id"] for row in records]
    if len(set(ids)) != len(ids):
        errors.append("duplicate defect ID")
    if expected is None:
        narrative = (root / "docs/planning/intersection-metamorphic-backfill.md").read_text()
        leads = {"PR-" + number for number in re.findall(r"/pull/(\d+)", narrative)}
        errors.extend(f"historical lead missing from inventory: {name}"
                      for name in sorted(leads - set(ids)))
    if base:
        missing = {r["id"] for r in base["defects"]} - set(ids)
        errors.extend(f"removed obligation: {name}" for name in sorted(missing))
    for row in records:
        label = row["id"]
        if not isinstance(label, str) or not label or row.get("state") not in STATES:
            errors.append(f"{label}: invalid ID/state")
            continue
        evidence(root, row, ["source"], errors, label)
        state = row["state"]
        if state in {"verified", "exception reviewed"}:
            evidence(root, row, ["fixture", "direct", "oracle", "replay", "red", "green", "review"], errors, label)
        if state == "verified":
            evidence(root, row, ["metamorphic"], errors, label)
            if not row.get("relations") or not set(row["relations"]) <= {f"MR-{n:02}" for n in range(1, 10)}:
                errors.append(f"{label}: missing/invalid relations")
        if state in {"exception reviewed", "not a defect"}:
            evidence(root, row, ["rationale", "review"], errors, label)
        if release and state not in {"verified", "exception reviewed", "not a defect"}:
            errors.append(f"{label}: open obligation ({state})")
    if release:
        if data.get("inventory_complete") is not True:
            errors.append("historical inventory is not complete")
        evidence(root, data, ["inventory_review"], errors, "inventory")
    matrix = data["handlers"]
    expected = handler_ids(root) if expected is None else expected
    if set(matrix) != expected:
        errors.append("handler matrix differs from built-in ordered-pair catalogue")
    for name, cells in matrix.items():
        if set(cells) != RELATIONS:
            errors.append(f"{name}: MR-01/02/03/04 must all be declared")
        for relation, cell in cells.items():
            label = f"{name}/{relation}"
            state = cell.get("state")
            if state not in {"untriaged", "verified", "exception reviewed"}:
                errors.append(f"{label}: invalid state (unresolved is not a pass)")
            if state in {"verified", "exception reviewed"}:
                evidence(root, cell, ["domain", "review"], errors, label)
            if state == "verified":
                evidence(root, cell, ["test", "oracle", "replay", "green"], errors, label)
            if state == "exception reviewed":
                evidence(root, cell, ["rationale"], errors, label)
            if release and state == "untriaged":
                errors.append(f"{label}: untriaged relation")
    return errors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--release", action="store_true")
    parser.add_argument("--base", help="git revision whose defect IDs must be retained")
    args = parser.parse_args()
    root = args.root.resolve()
    try:
        data = json.loads((root / LEDGER).read_text())
        base = None
        if args.base:
            # Only the migration's base lacks this file; check the tree before
            # reading it so an invalid ref or git error is never ignored.
            tree = subprocess.check_output(["git", "ls-tree", args.base, "--", LEDGER], cwd=root, text=True)
            if tree.strip():
                base = json.loads(subprocess.check_output(["git", "show", f"{args.base}:{LEDGER}"], cwd=root, text=True))
        errors = check(root, data, args.release, base)
    except (OSError, ValueError, KeyError, TypeError, AttributeError, subprocess.CalledProcessError) as exc:
        print(f"check_metamorphic: invalid inventory: {exc}", file=sys.stderr)
        return 1
    for error in errors:
        print(error, file=sys.stderr)
    print(f"check_metamorphic: {len(data['defects'])} defects, {len(data['handlers'])} ordered pairs; "
          f"{'release' if args.release else 'inventory only'}: {'FAIL' if errors else 'PASS'}")
    return bool(errors)


if __name__ == "__main__":
    sys.exit(main())
