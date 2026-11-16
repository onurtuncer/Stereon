#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
"""Negative and positive release-gate regressions, no C++ build required."""
import copy
import json
import shutil
import tempfile
import unittest
from pathlib import Path

from check_metamorphic import check, check_campaign, reference, RELATIONS


class GateTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root / "proof.md").write_text("fixture test oracle replay red green review domain rationale", encoding="utf-8")
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

    def test_implemented_relation_requires_test_but_is_not_qualified(self):
        cell = self.data["handlers"]["pair"]["MR-01"]
        cell["state"] = "MR implemented"
        self.assertEqual([], self.errors())
        self.assertTrue(self.errors(release=True))
        del cell["test"]
        self.assertTrue(self.errors())

    def test_implemented_defect_needs_replay_and_relation(self):
        row = self.data["defects"][0]
        row["state"] = "MR implemented"
        self.assertEqual([], self.errors())
        self.assertTrue(self.errors(release=True))
        del row["replay"]
        self.assertTrue(self.errors())

    def test_refs_and_anchors_must_exist(self):
        for ref in ("missing.md", "proof.md#absent", "../outside", "", "https://", "http://a/b"):
            self.assertFalse(reference(self.root, ref), ref)


class CampaignTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        repo = Path(__file__).resolve().parents[1]
        self.data = json.loads((repo / "docs/planning/intersection-metamorphic.json").read_text(encoding="utf-8"))
        self.path = self.root / "tools/mutation/intersection-backfill.json"
        catalogue = json.loads((repo / "tools/mutation/intersection-backfill.json").read_text(encoding="utf-8"))
        paths = {m["file"] for m in catalogue["mutants"]}
        paths.update(["tools/mutation/intersection-backfill.json", "tests/intersect/MetamorphicRegressionTests.cpp"])
        for relative in paths:
            target = self.root / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(repo / relative, target)
        self.catalogue = catalogue

    def errors(self):
        self.path.write_text(json.dumps(self.catalogue), encoding="utf-8")
        return check_campaign(self.root, self.data)

    def test_current_campaign_is_consistent(self):
        self.assertEqual([], self.errors())

    def test_stale_and_ambiguous_sites_are_rejected(self):
        row = self.catalogue["mutants"][0]
        path = self.root / row["file"]
        original = path.read_text(encoding="utf-8")
        path.write_text(original.replace(row["find"], ""), encoding="utf-8")
        self.assertTrue(self.errors())
        path.write_text(original + row["find"], encoding="utf-8")
        self.assertTrue(self.errors())

    def test_fault_must_change_production_source(self):
        row = self.catalogue["mutants"][0]
        row["replace"] = row["find"]
        self.assertTrue(self.errors())
        row["file"] = "tests/intersect/MetamorphicRegressionTests.cpp"
        self.assertTrue(self.errors())

    def test_selection_and_diagnostic_must_match_existing_test(self):
        row = self.catalogue["mutants"][0]
        row["tests"] = ".*"
        self.assertTrue(self.errors())
        row["diagnostic"] = "nonexistent"
        self.assertTrue(self.errors())

    def test_catalogue_cannot_be_empty_or_duplicate(self):
        row = self.catalogue["mutants"][0]
        self.catalogue["mutants"].append(row)
        self.assertTrue(self.errors())
        self.catalogue["mutants"] = []
        self.assertTrue(self.errors())

    def test_fault_cannot_lose_critical_or_repair_identity(self):
        row = self.catalogue["mutants"][0]
        row["critical"] = False
        self.assertTrue(self.errors())
        row["critical"] = True
        row["repair_commit"] = "0" * 40
        self.assertTrue(self.errors())


if __name__ == "__main__":
    unittest.main()
