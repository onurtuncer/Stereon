#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
"""Exercise missing, stale and duplicate source records and unresolved citations."""
from pathlib import Path
import subprocess
import tempfile
import unittest

from check_algorithm_docs import check, routines


class InventoryTest(unittest.TestCase):
    def test_private_definitions(self):
        """New routines in existing files need documentation in Python and C++."""
        for filename, code in [("sample.py", 'def old():\n    """Return the baseline."""\n    return 0\ndef added():\n    return 1\n'),
                               ("sample.cppm", "/// Return the baseline.\nstatic int old() { return 0; }\nstatic int added() { return 1; }")]:
            result = {name: documented for name, _, documented in routines(Path(filename), code).values()}
            self.assertEqual(result, {"old": True, "added": False})

    def test_inventory_changes(self):
        """Reject missing and duplicated files, stale paths and unknown bibliography keys."""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            subprocess.run(["git", "init", "-q"], cwd=root, check=True)
            pages = root / "docs/algorithms"
            pages.mkdir(parents=True)
            (root / "docs/references.bib").write_text("@book{source, year={2026}}")
            (root / "solver.cpp").touch()
            subprocess.run(["git", "add", "solver.cpp"], cwd=root, check=True)
            page = pages / "solver.rst"
            self.assertEqual(len(check(root)[1]), 1)
            entry = ".. algorithm-sources::\n   solver.cpp\n"
            page.write_text(":cite:p:`source`\n\n" + entry)
            self.assertEqual(check(root), (1, []))
            page.write_text(":cite:p:`unknown`\n\n" + entry + entry + "   deleted.hpp\n")
            self.assertEqual(len(check(root)[1]), 3)


if __name__ == "__main__":
    unittest.main()
