#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
"""Negative and positive release-gate regressions, no C++ build required."""
import copy
import tempfile
import unittest
from pathlib import Path

from check_metamorphic import check, reference, RELATIONS


class GateTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root / "proof.md").write_text("fixture test oracle replay red green review domain rationale")
        fields = {k: "proof.md#" + k for k in
                  ("fixture", "direct", "oracle", "replay", "red", "green", "review", "metamorphic")}
        fields.update(direct="proof.md#test", metamorphic="proof.md#test")
        self.data = {"version": 1, "inventory_complete": True, "inventory_review": "proof.md#review",
                     "defects": [dict(id="bug-1", state="verified", source="https://github.com/a/b/pull/1",
                                      relations=["MR-01"], **fields)],
                     "handlers": {"pair": {mr: dict(state="verified", domain="proof.md#domain",
                                  review="proof.md#review", test="proof.md#test", oracle="proof.md#oracle",
                                  replay="proof.md#replay", green="proof.md#green") for mr in RELATIONS}}}

    def errors(self, **kwargs):
        return check(self.root, self.data, expected={"pair"}, **kwargs)

    def test_complete_release(self):
        self.assertEqual([], self.errors(release=True))

    def test_open_backfill_is_development_only(self):
        self.data["inventory_complete"] = False
        self.data["defects"][0]["state"] = "untriaged"
        self.assertEqual([], self.errors())
        self.assertTrue(self.errors(release=True))

    def test_every_evidence_field_is_required(self):
        for key in ("fixture", "direct", "oracle", "replay", "red", "green", "review", "metamorphic"):
            with self.subTest(key=key):
                value = self.data["defects"][0].pop(key)
                self.assertTrue(self.errors())
                self.data["defects"][0][key] = value

    def test_duplicate_and_removed_obligations(self):
        base = copy.deepcopy(self.data)
        self.data["defects"].append(copy.deepcopy(self.data["defects"][0]))
        self.assertTrue(self.errors())
        self.data["defects"] = []
        self.assertTrue(self.errors(base=base))

    def test_matrix_cannot_omit_pairs_or_relations(self):
        del self.data["handlers"]["pair"]["MR-04"]
        self.assertTrue(self.errors())
        self.data["handlers"] = {}
        self.assertTrue(self.errors())

    def test_unresolved_and_untriaged_not_release_passes(self):
        cell = self.data["handlers"]["pair"]["MR-01"]
        cell["state"] = "unresolved"
        self.assertTrue(self.errors())
        cell["state"] = "untriaged"
        self.assertEqual([], self.errors())
        self.assertTrue(self.errors(release=True))

    def test_exception_needs_independent_regression_and_review(self):
        row = self.data["defects"][0]
        row["state"] = "exception reviewed"
        self.assertTrue(self.errors())
        row["rationale"] = "proof.md#rationale"
        self.assertEqual([], self.errors(release=True))
        del row["direct"]
        self.assertTrue(self.errors())

    def test_dismissed_candidate_requires_review(self):
        row = self.data["defects"][0]
        row["state"] = "not a defect"
        self.assertTrue(self.errors())
        row["rationale"] = "proof.md#rationale"
        self.assertEqual([], self.errors(release=True))
        del row["review"]
        self.assertTrue(self.errors())

    def test_complete_inventory_needs_review(self):
        del self.data["inventory_review"]
        self.assertTrue(self.errors(release=True))

    def test_matrix_exception_needs_domain_and_review(self):
        cell = self.data["handlers"]["pair"]["MR-02"]
        cell["state"] = "exception reviewed"
        self.assertTrue(self.errors())
        cell["rationale"] = "proof.md#rationale"
        self.assertEqual([], self.errors(release=True))
        del cell["domain"]
        self.assertTrue(self.errors())

    def test_unknown_version_and_defect_state_fail(self):
        self.data["version"] = 2
        self.assertTrue(self.errors())
        self.data["version"] = 1
        self.data["defects"][0]["state"] = "done"
        self.assertTrue(self.errors())

    def test_refs_and_anchors_must_exist(self):
        for ref in ("missing.md", "proof.md#absent", "../outside", "", "https://", "http://a/b"):
            self.assertFalse(reference(self.root, ref), ref)


if __name__ == "__main__":
    unittest.main()
