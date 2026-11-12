#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""Regression cases for release evidence validation; no native build required."""

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import check_release_benchmarks as checker


def convergence(cases=10000, incorrect=0):
    """The schema emitted by ProjectionCorpus.cpp, with one optional surface failure."""
    def tally(failures):
        return {"cases": cases, "converged": cases - failures,
                "incorrect_success": failures, "worst_excess": 0.0, "failures": []}
    return {"corpus": "b13-v1", "seed": 20261101, "allowance": 1e-9,
            "cases": 2 * cases, "converged": 2 * cases - incorrect,
            "incorrect_success": incorrect, "rate": (2 * cases - incorrect) / (2 * cases),
            "lower_bound_95": 0.9998,
            "by_kind": {"curves": tally(0), "surfaces": tally(incorrect)}}


class ReleaseEvidenceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.evidence = self.root / "convergence.json"
        self.run = checker._passing_run()
        self.machine = {key: "x" for key in checker.REQUIRED_MACHINE} | {"commit": "abc"}
        self.write_evidence(convergence())

    def write_evidence(self, data):
        self.evidence.write_text(json.dumps(data), encoding="utf-8")

    def gates(self):
        return checker.check(self.run, self.machine, "abc", self.evidence)

    def test_valid_corpus_and_full_run_pass(self):
        self.assertTrue(checker.released(self.gates()))

    def test_real_incorrect_success_fails(self):
        self.write_evidence(convergence(incorrect=1))
        self.assertEqual(checker.convergence_gate(self.evidence).status, "FAIL")

    def test_invalid_counts_are_rejected_without_coercion(self):
        for field, value in (("converged", 20001), ("cases", -1),
                             ("incorrect_success", -1), ("converged", 20000.5),
                             ("converged", "20000"), ("incorrect_success", False),
                             ("converged", float("inf"))):
            with self.subTest(field=field, value=value):
                data = convergence()
                data[field] = value
                self.write_evidence(data)
                self.assertEqual(checker.convergence_gate(self.evidence).status, "INVALID")

    def test_missing_or_wrong_corpus_identity_is_rejected(self):
        for field, value in (("corpus", None), ("corpus", "other"), ("seed", 42),
                             ("allowance", 1.0)):
            with self.subTest(field=field, value=value):
                data = convergence()
                if value is None:
                    del data[field]
                else:
                    data[field] = value
                self.write_evidence(data)
                self.assertEqual(checker.convergence_gate(self.evidence).status, "INVALID")

    def test_short_corpus_is_not_release_evidence(self):
        self.write_evidence(convergence(cases=1))
        self.assertEqual(checker.convergence_gate(self.evidence).status, "INVALID")

    def test_subtotals_must_match(self):
        data = convergence()
        data["by_kind"]["surfaces"].update(converged=9999, incorrect_success=1)
        self.write_evidence(data)
        self.assertEqual(checker.convergence_gate(self.evidence).status, "INVALID")

    def test_malformed_structure_is_invalid_not_an_exception(self):
        for data in (None, [], {}, {**convergence(), "by_kind": []},
                     {**convergence(), "by_kind": {"curves": None}}):
            with self.subTest(data=data):
                self.write_evidence(data)
                self.assertEqual(checker.convergence_gate(self.evidence).status, "INVALID")

    def test_error_repetition_blocks_release_even_with_enough_good_samples(self):
        self.run["benchmarks"].append({"name": "B1.3/CurveProjection/Stereon",
                                       "run_type": "iteration", "error_occurred": True,
                                       "error_message": "projection failed"})
        self.assertFalse(checker.released(self.gates()))

    def test_false_error_flag_is_a_valid_sample(self):
        for entry in self.run["benchmarks"]:
            entry["error_occurred"] = False
        self.assertTrue(checker.released(self.gates()))

    def test_counter_gates_require_five_samples_on_each_needed_side(self):
        for key, name, _, _, check in checker.br.COUNTER_ROWS:
            if check is None:
                continue
            libraries = ("Stereon", "OCCT") if check[0] == "ratio below" else ("Stereon",)
            for library in libraries:
                with self.subTest(key=key, library=library):
                    self.run = checker._passing_run()
                    prefix = f"{key}/{library}"
                    matches = [e for e in self.run["benchmarks"] if e["name"] == prefix]
                    self.run["benchmarks"] = [e for e in self.run["benchmarks"]
                                               if e["name"] != prefix] + matches[:1]
                    gates = self.gates()
                    self.assertEqual(checker._status(gates, f"{key}:{name}"), "INVALID")
                    self.assertFalse(checker.released(gates))

    def test_counter_must_be_present_and_finite_in_every_repetition(self):
        key = "B0.4/HeapPerFace/1000000"
        name = "heap_bytes_per_face"
        for bad in (None, float("nan"), float("inf"), -1.0):
            with self.subTest(bad=bad):
                self.run = checker._passing_run()
                entry = next(e for e in self.run["benchmarks"] if e["name"] == key + "/Stereon")
                if bad is None:
                    del entry[name]
                else:
                    entry[name] = bad
                self.assertEqual(checker._status(self.gates(), f"{key}:{name}"), "INVALID")

    def test_cli_enforce_and_report_verdict_agree(self):
        run_path = self.root / "results.json"
        machine_path = self.root / "machine.txt"
        machine_path.write_text("\n".join(f"{k}={v}" for k, v in self.machine.items()), encoding="utf-8")
        output = self.root / "gates.json"
        command = [sys.executable, str(Path(checker.__file__)), str(run_path),
                   "--machine", str(machine_path), "--candidate", "abc",
                   "--convergence", str(self.evidence), "--json", str(output)]
        for fault in (None, "counts", "benchmark error", "short counter"):
            run = checker._passing_run()
            data = convergence()
            if fault == "counts":
                data["converged"] = 20001
            elif fault == "benchmark error":
                run["benchmarks"].append({"name": "B1.3/CurveProjection/Stereon",
                                          "error_occurred": True, "error_message": "failed"})
            elif fault == "short counter":
                prefix = "B0.4/HeapPerFace/1000000/Stereon"
                matches = [e for e in run["benchmarks"] if e["name"] == prefix]
                run["benchmarks"] = [e for e in run["benchmarks"]
                                     if e["name"] != prefix] + matches[:1]
            run_path.write_text(json.dumps(run), encoding="utf-8")
            self.write_evidence(data)
            for enforce in (False, True):
                with self.subTest(fault=fault, enforce=enforce):
                    result = subprocess.run(command + (["--enforce"] if enforce else []),
                                            capture_output=True, text=True)
                    self.assertEqual(result.returncode, int(fault is not None and enforce), result.stderr)
                    self.assertEqual(json.loads(output.read_text(encoding="utf-8"))["verdict"],
                                     "FAIL" if fault else "PASS")


if __name__ == "__main__":
    unittest.main()
