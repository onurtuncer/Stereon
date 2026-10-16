#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""Source-level mutation harness for ADR-0032 (mutation-based verification).

Each mutant in a catalogue (tools/mutation/*.json) is one conceptual fault:
an exact, unique piece of production source and its replacement, the
contract it challenges, the test targets to build (one name or a list) and
the CTest selection that should detect it. The harness

1. checks out the current commit into a disposable git worktree (never the
   developer's tree or a branch) and configures it with the reference GCC
   preset;
2. builds and runs every selected test set unmodified: the baseline must
   pass before any mutant counts;
3. applies each mutant alone, rebuilds incrementally, runs its tests and
   restores the file, classifying the result as ADR-0032 section 3 defines:

   Killed       the tests fail, reproducibly (run twice), or a compile-time
                test (a static_assert in a file under tests/) rejects the
                mutant;
   Survived     the tests pass, and an execution trap at the same site shows
                the tests reach it;
   Uncovered    the tests pass, and also pass with the trap: never reached;
   Equivalent   the catalogue records an expert rationale that the mutant
                cannot be observed, and it survived (a kill refutes the
                claim and is reported as such);
   Invalid      the production code does not build, or the mutant's source
                text is not unique;
   Inconclusive a timeout, an infrastructure failure or a flaky result.

The mutation score is killed / (killed + survived + uncovered); equivalent,
invalid and inconclusive mutants are excluded and listed. A critical mutant
must be killed with its declared diagnostic: the failing tests must include
one of the named ones ("diagnostic" is a string or a list). Each mutant's
build and test output is kept in <work>/logs as diagnostic evidence. The
exit status is non-zero if the baseline fails or any critical mutant is not
killed as declared.

Mutants never touch tests or oracles, and nothing is committed: the
worktree is removed at the end. Run from anywhere inside the repository:

    python tools/mutate.py --catalogue tools/mutation/pilot.json
"""

import argparse
import json
import platform
import re
import shutil
import subprocess
import sys
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path

SCORED = ("killed", "survived", "uncovered")
TRAP = "__builtin_trap()"


@dataclass
class Outcome:
    id: str
    critical: bool
    operator: str
    contract: str
    location: str
    tests: str
    outcome: str = "inconclusive"
    detail: str = ""
    failed_tests: list = field(default_factory=list)
    diagnostic_found: bool = False
    seconds: float = 0.0


def run(command, cwd, timeout):
    """Runs a command; returns (exit status or None on timeout, output)."""
    try:
        done = subprocess.run(command, cwd=cwd, capture_output=True, text=True,
                              errors="replace", timeout=timeout, check=False)
        return done.returncode, done.stdout + done.stderr
    except subprocess.TimeoutExpired as expired:
        output = (expired.stdout or b"") + (expired.stderr or b"")
        if isinstance(output, bytes):
            output = output.decode(errors="replace")
        return None, output


def targets_of(mutant):
    """The build targets of a mutant: "target" is one name or a list."""
    target = mutant["target"]
    return [target] if isinstance(target, str) else list(target)


def diagnostics_of(mutant):
    """The declared diagnostics: "diagnostic" is one string or a list."""
    declared = mutant.get("diagnostic", [])
    return [declared] if isinstance(declared, str) else list(declared)


def compile_time_failures(output):
    """The static assertions in test files that a build rejected: a
    compile-time test detecting the mutant, not a broken build."""
    found = []
    for line in output.splitlines():
        match = re.search(r"([^\s:]*tests[/\\][^\s:]+):(\d+):\d+: error: static assertion "
                          r"failed(.*)", line)
        if match:
            name = Path(match.group(1)).name
            found.append(f"compile-time: {name}:{match.group(2)}{match.group(3)}".strip())
    return found


def no_tests(output):
    return "No tests were found" in output


def failed_tests(output):
    """The names of the failed tests in CTest's summary."""
    names = []
    for line in output.splitlines():
        match = re.match(r"\s*\d+ - (.*?) \((Failed|Timeout|SEGFAULT|Subprocess aborted|"
                         r"Exception|Child aborted|Not Run|[A-Za-z ]+)\)\s*$", line)
        if match:
            names.append((match.group(1).strip(), match.group(2)))
    return names


class Worktree:
    """A disposable checkout of HEAD with its own build directory."""

    def __init__(self, repo, work, preset):
        self.repo = repo
        self.root = work / "tree"
        self.preset = preset
        self.build = self.root / "build" / preset

    def create(self):
        self.remove()
        self.root.parent.mkdir(parents=True, exist_ok=True)
        status, output = run(["git", "worktree", "add", "--detach", str(self.root), "HEAD"],
                             self.repo, 300)
        if status != 0:
            sys.exit("cannot create the worktree:\n" + output)

    def remove(self):
        run(["git", "worktree", "remove", "--force", str(self.root)], self.repo, 300)
        if self.root.exists():
            shutil.rmtree(self.root, ignore_errors=True)
        run(["git", "worktree", "prune"], self.repo, 60)

    def configure(self):
        status, output = run(["cmake", "--preset", self.preset], self.root, 1800)
        if status != 0:
            sys.exit("cannot configure the worktree:\n" + output[-4000:])

    def build_targets(self, targets, timeout):
        return run(["cmake", "--build", str(self.build), "--target", *targets], self.root,
                   timeout)

    def test(self, selection, test_timeout, timeout):
        return run(["ctest", "--test-dir", str(self.build), "-R", selection,
                    "--timeout", str(test_timeout), "--output-on-failure"], self.root,
                   timeout)

    def restore(self, relative):
        run(["git", "checkout", "--", relative], self.root, 60)


class Harness:
    def __init__(self, tree, args):
        self.tree = tree
        self.args = args

    def apply(self, mutant, replacement):
        """Replaces the mutant's unique source text; False if not unique."""
        path = self.tree.root / mutant["file"]
        text = path.read_text(encoding="utf-8")
        if text.count(mutant["find"]) != 1:
            return False
        path.write_text(text.replace(mutant["find"], replacement), encoding="utf-8",
                        newline="\n")
        return True

    def build_and_test(self, mutant):
        """('invalid'|'compile-killed'|'timeout'|'empty'|'pass'|'fail', output) for
        the current source."""
        status, output = self.tree.build_targets(targets_of(mutant), self.args.build_timeout)
        if status is None:
            return "timeout", output
        if status != 0:
            return ("compile-killed" if compile_time_failures(output) else "invalid"), output
        test_timeout = mutant.get("testTimeout", self.args.test_timeout)
        status, test_output = self.tree.test(mutant["tests"], test_timeout,
                                             self.args.run_timeout)
        output += test_output
        if no_tests(test_output):
            return "empty", output
        if status is None or any(kind == "Timeout" for _, kind in failed_tests(test_output)):
            return "timeout", output
        return ("pass" if status == 0 else "fail"), output

    def trap_reaches(self, mutant):
        """Whether the tests reach the mutant's site: an execution trap there
        makes them fail. None when the trap cannot be built."""
        find = mutant["find"]
        indent = find[: len(find) - len(find.lstrip())]
        probe = indent + TRAP + ";\n" + find if mutant["site"] == "statement" else \
            "(" + TRAP + ", " + find + ")"
        if not self.apply(mutant, probe):
            return None
        try:
            result, _ = self.build_and_test(mutant)
        finally:
            self.tree.restore(mutant["file"])
        if result in ("invalid", "timeout", "empty"):
            return None
        return result in ("fail", "compile-killed")

    def classify(self, mutant, outcome):
        if not self.apply(mutant, mutant["replace"]):
            outcome.outcome, outcome.detail = "invalid", "the source text is not unique"
            return
        try:
            result, output = self.build_and_test(mutant)
            if result == "fail":
                again, _ = self.build_and_test(mutant)
                result = "fail" if again == "fail" else "flaky"
        finally:
            self.tree.restore(mutant["file"])
        self.keep_log(mutant, output)
        self.record(mutant, outcome, result, output)

    def keep_log(self, mutant, output):
        logs = self.tree.root.parent / "logs"
        logs.mkdir(parents=True, exist_ok=True)
        (logs / f"{mutant['id']}.txt").write_text(output, encoding="utf-8")

    def record(self, mutant, outcome, result, output):
        if result == "invalid":
            first = next((l for l in output.splitlines() if "error" in l), "build failed")
            outcome.outcome, outcome.detail = "invalid", first.strip()[:300]
        elif result in ("timeout", "flaky", "empty"):
            outcome.outcome = "inconclusive"
            outcome.detail = {"timeout": "a test timed out",
                              "flaky": "failed once, passed once",
                              "empty": "the test selection matches no test"}[result]
        elif result in ("fail", "compile-killed"):
            outcome.failed_tests = compile_time_failures(output) + \
                [name for name, _ in failed_tests(output)]
            declared = diagnostics_of(mutant)
            outcome.diagnostic_found = not declared or any(
                d in failed for d in declared for failed in outcome.failed_tests)
            outcome.outcome = "killed"
            if mutant.get("equivalent"):
                outcome.detail = "the equivalence claim is refuted by this kill"
        else:
            self.record_pass(mutant, outcome)

    def record_pass(self, mutant, outcome):
        if mutant.get("equivalent"):
            outcome.outcome, outcome.detail = "equivalent", mutant["equivalent"]
            return
        reached = self.trap_reaches(mutant)
        if reached is None:
            outcome.outcome = "survived"
            outcome.detail = "coverage unknown: the execution trap could not be built here"
        elif reached:
            outcome.outcome, outcome.detail = "survived", "the tests reach the site"
        else:
            outcome.outcome, outcome.detail = "uncovered", "the tests never reach the site"


def baseline(tree, mutants, args):
    """Every target builds and every selection passes, unmodified."""
    targets = sorted({t for m in mutants for t in targets_of(m)})
    status, output = tree.build_targets(targets, args.build_timeout * 4)
    if status != 0:
        sys.exit("baseline: the targets do not build:\n" + output[-4000:])
    for selection in sorted({m["tests"] for m in mutants}):
        status, output = tree.test(selection, args.test_timeout, args.run_timeout * 2)
        if no_tests(output):
            sys.exit(f"baseline: the selection '{selection}' matches no test")
        if status != 0:
            sys.exit(f"baseline: tests '{selection}' fail unmodified:\n" + output[-4000:])


def tool_version(command):
    status, output = run(command, None, 60)
    return output.splitlines()[0].strip() if status == 0 and output else "unknown"


def summarise(outcomes, catalogue, revision, seconds):
    counts = {name: sum(o.outcome == name for o in outcomes)
              for name in ("killed", "survived", "uncovered", "equivalent", "invalid",
                           "inconclusive")}
    denominator = sum(counts[name] for name in SCORED)
    critical = [o for o in outcomes if o.critical]
    return {
        "catalogue": catalogue,
        "revision": revision,
        "compiler": tool_version(["g++", "--version"]),
        "cmake": tool_version(["cmake", "--version"]),
        "platform": platform.platform(),
        "seconds": round(seconds, 1),
        "counts": counts,
        "score": {"killed": counts["killed"], "denominator": denominator,
                  "value": counts["killed"] / denominator if denominator else None},
        "critical": {"total": len(critical),
                     "detected": sum(o.outcome == "killed" and o.diagnostic_found
                                     for o in critical)},
        "mutants": [asdict(o) for o in outcomes],
    }


def markdown(report):
    score = report["score"]
    value = f"{score['value']:.2f}" if score["value"] is not None else "n/a"
    lines = [
        f"# Mutation report: {Path(report['catalogue']).name}",
        "",
        f"- Revision: `{report['revision']}`",
        f"- Compiler: {report['compiler']}; {report['cmake']}; {report['platform']}",
        f"- Wall time: {report['seconds']} s",
        f"- Score: {score['killed']} killed / {score['denominator']} scored = {value} "
        "(killed / (killed + survived + uncovered))",
        f"- Critical faults detected with their diagnostic: {report['critical']['detected']}"
        f" of {report['critical']['total']}",
        f"- Counts: {', '.join(f'{k} {v}' for k, v in report['counts'].items())}",
        "",
        "| Mutant | Critical | Outcome | Contract | Detail |",
        "| --- | --- | --- | --- | --- |",
    ]
    for m in report["mutants"]:
        detail = m["detail"]
        if m["outcome"] == "killed":
            detail = "failed: " + "; ".join(m["failed_tests"][:3])
            if m["critical"]:
                detail += " (declared diagnostic " + \
                          ("found)" if m["diagnostic_found"] else "MISSING)")
        detail = detail.replace("|", "\\|")
        lines.append(f"| {m['id']} | {'yes' if m['critical'] else ''} | {m['outcome']} | "
                     f"{m['contract']} | {detail} |")
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--catalogue", required=True, type=Path)
    parser.add_argument("--work", type=Path, default=Path("build/mutation"))
    parser.add_argument("--only", default="", help="comma-separated mutant ids")
    parser.add_argument("--build-timeout", type=int, default=1800)
    parser.add_argument("--test-timeout", type=int, default=600)
    parser.add_argument("--run-timeout", type=int, default=3600)
    parser.add_argument("--keep", action="store_true", help="keep the worktree")
    args = parser.parse_args()

    repo = Path(run(["git", "rev-parse", "--show-toplevel"], None, 60)[1].strip())
    revision = run(["git", "rev-parse", "HEAD"], repo, 60)[1].strip()
    catalogue = json.loads((repo / args.catalogue).read_text(encoding="utf-8"))
    mutants = catalogue["mutants"]
    if args.only:
        wanted = set(args.only.split(","))
        mutants = [m for m in mutants if m["id"] in wanted]
    work = (repo / args.work).resolve()
    tree = Worktree(repo, work, catalogue.get("preset", "gcc16-debug"))
    start = time.monotonic()
    tree.create()
    try:
        tree.configure()
        baseline(tree, mutants, args)
        harness = Harness(tree, args)
        outcomes = []
        for mutant in mutants:
            began = time.monotonic()
            outcome = Outcome(mutant["id"], mutant.get("critical", False), mutant["operator"],
                              mutant["contract"], mutant["file"], mutant["tests"])
            harness.classify(mutant, outcome)
            outcome.seconds = round(time.monotonic() - began, 1)
            print(f"{outcome.id}: {outcome.outcome} ({outcome.seconds} s) {outcome.detail}",
                  flush=True)
            outcomes.append(outcome)
    finally:
        if not args.keep:
            tree.remove()
    report = summarise(outcomes, str(args.catalogue), revision, time.monotonic() - start)
    work.mkdir(parents=True, exist_ok=True)
    (work / "report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    (work / "report.md").write_text(markdown(report), encoding="utf-8")
    print(markdown(report))
    missed = [o.id for o in outcomes if o.critical
              and not (o.outcome == "killed" and o.diagnostic_found)]
    if missed:
        print("critical faults not detected as declared: " + ", ".join(missed))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
