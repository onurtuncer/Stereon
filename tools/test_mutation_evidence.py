#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
"""Evidence must be a repeated named assertion failure, not an empty run."""
import unittest
from unittest.mock import Mock
from mutate import Harness, Outcome, select_mutants


class EvidenceTests(unittest.TestCase):
    def classify(self, results):
        harness = Harness(Mock(), Mock())
        harness.apply = Mock(return_value=True)
        harness.keep_log = Mock()
        harness.build_and_test = Mock(side_effect=results)
        mutant = dict(id="PR-205", file="source.cpp", replace="fault", diagnostic="MR PR-205")
        outcome = Outcome("PR-205", True, "historical", "contract", "source.cpp", "selection")
        harness.classify(mutant, outcome)
        harness.tree.restore.assert_called_once_with("source.cpp")
        return outcome, harness

    def test_repeated_declared_failure_is_killed(self):
        result, h = self.classify([("fail", "1 - intersect: MR PR-205 (Failed)\n")] * 2)
        self.assertEqual("killed", result.outcome)
        self.assertTrue(result.diagnostic_found)
        self.assertIn("repeat run", h.keep_log.call_args.args[1])

    def test_different_failure_is_not_proof(self):
        for second in ("2 - unrelated (Failed)\n", "1 - intersect: MR PR-205 (SEGFAULT)\n"):
            result, _ = self.classify([("fail", "1 - intersect: MR PR-205 (Failed)\n"), ("fail", second)])
            self.assertEqual("inconclusive", result.outcome)

    def test_production_compile_error_is_invalid(self):
        result, _ = self.classify([("invalid", "source.cpp: error: bad code")])
        self.assertEqual("invalid", result.outcome)
        self.assertFalse(result.diagnostic_found)

    def test_timeout_and_empty_selection_are_inconclusive(self):
        for kind in ("timeout", "empty"):
            result, _ = self.classify([(kind, "")])
            self.assertEqual("inconclusive", result.outcome)

    def test_flaky_failure_is_inconclusive(self):
        result, _ = self.classify([("fail", "1 - intersect: MR PR-205 (Failed)\n"), ("pass", "")])
        self.assertEqual("inconclusive", result.outcome)

    def test_unknown_empty_and_duplicate_catalogues_are_errors(self):
        for rows, selection in (([], ""), ([{"id": "a"}], "typo"), ([{"id": "a"}, {"id": "a"}], "")):
            with self.assertRaises(ValueError):
                select_mutants(rows, selection)
        self.assertEqual([{"id": "a"}], select_mutants([{"id": "a"}, {"id": "b"}], "a"))


if __name__ == "__main__":
    unittest.main()
