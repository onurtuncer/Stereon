#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
"""Check ADR numbering, index completeness, local links and numeric references."""

import argparse
from pathlib import Path
import re
import tempfile


def check(root):
    directory = root / "docs" / "adr"
    files = sorted(directory.glob("[0-9][0-9][0-9][0-9]-*.md"))
    errors = []
    numbers = {}
    for path in files:
        number = path.name[:4]
        if number in numbers:
            errors.append(f"duplicate ADR number: {number}")
        numbers[number] = path
    index_path = directory / "README.md"
    if not index_path.is_file():
        return ["missing ADR index"]
    index = index_path.read_text(encoding="utf-8")
    entries = re.findall(r"\[([0-9]{4})\]\(([^)]+)\)", index)
    for number, path in numbers.items():
        if number != "0000" and entries.count((number, path.name)) != 1:
            errors.append(f"missing or duplicate index entry: {path.name}")
    for number, target in entries:
        if number not in numbers or numbers[number].name != target:
            errors.append(f"incorrect index target: {number}: {target}")
    for path in [index_path, *files]:
        text = path.read_text(encoding="utf-8")
        for number in set(re.findall(r"\bADR-([0-9]{4})\b", text)):
            if number not in numbers:
                errors.append(f"{path.name}: missing ADR-{number}")
        for target in re.findall(r"\]\(([^\s)]+)\)", text):
            if re.match(r"[a-zA-Z][a-zA-Z0-9+.-]*:", target) or target.startswith("#"):
                continue
            target = target.split("#", 1)[0]
            if target and not (path.parent / target).exists():
                errors.append(f"{path.name}: broken local link: {target}")
    return errors


def self_test():
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        directory = root / "docs" / "adr"
        directory.mkdir(parents=True)
        index = directory / "README.md"
        adr = directory / "0001-example.md"
        index.write_text("[0001](0001-example.md)\n", encoding="utf-8")
        adr.write_text("# ADR-0001\n", encoding="utf-8")
        assert not check(root)
        adr.write_text("# ADR-0001\nADR-0002\n", encoding="utf-8")
        assert any("missing ADR-0002" in error for error in check(root))
        second = directory / "0002-second.md"
        second.write_text("# ADR-0002\n", encoding="utf-8")
        assert any("index entry" in error for error in check(root))
        index.write_text(index.read_text() + "[0002](0002-second.md)\n", encoding="utf-8")
        assert not check(root)
        second.write_text("[broken](missing.md)\n", encoding="utf-8")
        assert any("broken local link" in error for error in check(root))
        (directory / "0002-duplicate.md").write_text("duplicate\n", encoding="utf-8")
        assert any("duplicate ADR number" in error for error in check(root))
    print("ADR checker self-test passed")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return 0
    errors = check(args.root)
    for error in errors:
        print(error)
    if errors:
        return 1
    print("ADR references and index passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
