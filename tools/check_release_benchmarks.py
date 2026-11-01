#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""The benchmark release gates of v0.1, checked rather than rendered.

tools/bench_report.py renders a reference run as a table and always succeeds;
a row below its target is a "no" in the table, not a failure. This script reads
the same Google Benchmark JSON and machine record, with the same parser and
catalogue, and decides every gate the release needs (docs/benchmarks.md,
ADR-0003, the Phase 0 exit gate of TODO.md), each as one of

    PASS            the measurement meets the target
    FAIL            it does not, or the two libraries did not compute the same thing
    INCONCLUSIVE    the row's spread exceeds MAX_SPREAD and its margin to the
                    target is within the spread, so a re-run could turn it
    MISSING         a row, a library's twin or a required evidence file is absent
    INVALID         a timing is not a positive finite number, or a row has fewer
                    repetitions than the protocol's
    REPORTED        a row without a target, shown and never gating
    NOT_COMPARABLE  a row whose libraries do not solve the same problem
                    (B1.4's approximation): shown, never a speed pass

Decisions use the full-precision medians, not the rounded figures of the report:
0.899 is below the parity band's 0.90 and 1.984 below 2.

    tools/check_release_benchmarks.py results.json results-topology.json \\
        --machine machine.txt --candidate <sha> --json gate.json --out gate.md --enforce

Without --enforce it writes the same verdicts and exits 0; with it, it exits 1
unless every required gate is PASS. --self-test checks the checker.

Re-runs: an INCONCLUSIVE gate is re-measured by dispatching the reference
workflow again on the same commit; every run is kept, and the release report
cites all of them (docs/benchmarks.md, *Repetitions and spread*).
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import bench_report as br  # noqa: E402  (the same parser and catalogue)

# The repetitions of the approved protocol (docs/benchmarks.md, revised and
# approved 2026-10-31): curve and surface rows 15, topology rows 5.
CURVE_REPETITIONS = 15
TOPOLOGY_REPETITIONS = 5
TOPOLOGY_PREFIXES = ("B0.3/", "B0.4/", "B0.5/", "ADR-0004/", "ADR-0005/")

# A gate is INCONCLUSIVE when its row's spread (the largest deviation of a
# repetition from the median, over both libraries) exceeds this and the
# relative margin to the target is no larger than the spread.
MAX_SPREAD = 0.10

# The parity band, inclusive: a ratio OCCT / Stereon of at least this meets
# a "≥ parity" target (docs/benchmarks.md; bench_report.PARITY_BAND).
PARITY_FLOOR = 1.0 - br.PARITY_BAND

# ADR-0003's virtual fallback: at most this many nanoseconds per call over the
# direct call on the concrete type.
FALLBACK_NS = 5.0

# Two libraries' results of the same case agree to this, relative, or the
# speed comparison is not of the same computation (the *_mm counters).
AGREEMENT = 1e-9

# B1.4's interpolation is the same task in both libraries (a cubic through the
# same points, chord-length parameters) only if both interpolated: each side's
# largest distance from the data points (fit_residual) within the tolerance the
# benchmark gives OCCT's GeomAPI_Interpolate, in millimetres.
FIT_CHECKED = {"B1.4/Interpolation/1000": 1e-7}

# B1.3's convergence target (docs/benchmarks.md): measured on a declared corpus
# and supplied as an evidence file; no timing proves it.
CONVERGENCE_TARGET = 0.9999

# The machine record's fields a release gate is tied to.
REQUIRED_MACHINE = ("commit", "compiler", "cpu", "occt_tag", "occt_commit", "occt_compiler",
                    "occt_build_type")

PASSING = {"PASS"}
NON_GATING = {"REPORTED", "NOT_COMPARABLE"}


@dataclass
class Gate:
    gate: str
    label: str
    target: str
    status: str
    stereon_ns: float | None = None
    occt_ns: float | None = None
    value: float | None = None
    threshold: float | None = None
    spread: float | None = None
    repetitions: dict[str, int] = field(default_factory=dict)
    note: str = ""


# --- Checks of one row --------------------------------------------------------


def expected_repetitions(key: str) -> int:
    return TOPOLOGY_REPETITIONS if key.startswith(TOPOLOGY_PREFIXES) else CURVE_REPETITIONS


def invalid(row: br.Row, libraries: tuple[str, ...]) -> str:
    """Why the row's timings cannot be used, or ''."""
    for lib in libraries:
        side = row.sides.get(lib)
        if side is None:
            continue
        if any(not (math.isfinite(t) and t > 0) for t in side.per_item_ns):
            return f"{lib}: a timing that is not a positive finite number"
        if len(side.per_item_ns) < expected_repetitions(row.key):
            return (f"{lib}: {len(side.per_item_ns)} repetitions, the protocol's "
                    f"{expected_repetitions(row.key)} required")
    return ""


def disagreement(row: br.Row) -> str:
    """The first *_mm counter the two libraries disagree on beyond AGREEMENT, or ''."""
    if any(lib not in row.sides for lib in br.LIBRARIES):
        return ""
    stereon, occt = (row.sides[lib] for lib in br.LIBRARIES)
    for name in sorted(set(stereon.counters) & set(occt.counters)):
        if not name.endswith("_mm"):
            continue
        a, b = stereon.counter(name), occt.counter(name)
        if a is None or b is None or not (math.isfinite(a) and math.isfinite(b)):
            return f"{name}: not a finite value on both sides"
        if abs(a - b) > AGREEMENT * max(abs(a), abs(b), 1.0):
            return f"{name}: Stereon {a!r} against OCCT {b!r}"
    return ""


def unfitted(row: br.Row) -> str:
    """For a fitting row, the side whose fit missed the data points, or ''."""
    tolerance = FIT_CHECKED.get(row.key)
    if tolerance is None:
        return ""
    for lib in br.LIBRARIES:
        residual = row.sides[lib].counter("fit_residual")
        if residual is None:
            return f"{lib} recorded no fit_residual"
        if not (math.isfinite(residual) and residual <= tolerance):
            return f"{lib}'s fit is {residual!r} mm from its points, beyond {tolerance}"
    return ""


def judged(margin: float, spread: float, met: bool) -> str:
    """PASS or FAIL, or INCONCLUSIVE when the spread is excessive and could turn it."""
    if spread > MAX_SPREAD and abs(margin) <= spread:
        return "INCONCLUSIVE"
    return "PASS" if met else "FAIL"


def time_gate(row: br.Row | None, key: str) -> Gate:
    """A row of the time table against its target."""
    label, target = br.ROWS.get(key, (key, "≥ parity"))
    gate = Gate(key, label, target, "MISSING")
    if row is None:
        gate.note = "the row did not run"
        return gate
    gate.repetitions = {lib: len(side.per_item_ns) for lib, side in row.sides.items()}
    gate.spread = row.spread
    for lib in br.LIBRARIES:
        if lib in row.sides:
            setattr(gate, f"{lib.lower()}_ns", row.sides[lib].median)
    reason = invalid(row, br.LIBRARIES)
    if reason:
        gate.status, gate.note = "INVALID", reason
        return gate
    if target == "report":
        gate.status, gate.value = "REPORTED", row.ratio
        return gate
    if key in br.NOT_LIKE_FOR_LIKE:
        gate.status, gate.value = "NOT_COMPARABLE", row.ratio
        gate.note = "the libraries do not solve the same problem; a speed ratio here is no pass"
        return gate
    missing = [lib for lib in br.LIBRARIES if lib not in row.sides]
    if missing:
        gate.note = f"no {', '.join(missing)} run"
        return gate
    reason = disagreement(row) or unfitted(row)
    if reason:
        gate.status, gate.note = "FAIL", f"not the same result: {reason}"
        return gate
    ratio = row.ratio
    gate.value = ratio
    if key in br.MUST_LEAD:
        gate.threshold = 1.0
        gate.status = judged(ratio - 1.0, row.spread, ratio > 1.0)
        gate.note = "must lead outright"
    else:
        gate.threshold = br.MIN_RATIO.get(key, PARITY_FLOOR)
        gate.status = judged(ratio / gate.threshold - 1.0, row.spread, ratio >= gate.threshold)
    return gate


def stereon_median(by_key: dict[str, br.Row], key: str) -> tuple[float | None, str]:
    """Stereon's median of a row, or None and why not."""
    row = by_key.get(key)
    if row is None or "Stereon" not in row.sides:
        return None, f"{key} did not run"
    reason = invalid(row, ("Stereon",))
    if reason:
        return None, f"{key}: {reason}"
    return row.sides["Stereon"].median, ""


def dispatch_gates(by_key: dict[str, br.Row]) -> list[Gate]:
    """ADR-0003 (revised 2026-11-01): VisitBuiltin over the same loop on a local copy, at most
    1.10, and the virtual fallback at most 5 ns per call over the direct call. The direct ratio
    is reported beside them."""
    gates: list[Gate] = []
    for curve, name in br.DISPATCH_CURVES.items():
        times: dict[str, float] = {}
        problems: list[str] = []
        for path in ("Direct", "DispatchCopy", "Dispatch", "Virtual"):
            value, why = stereon_median(by_key, f"B0.2/{path}/{curve}")
            if value is None:
                problems.append(why)
            else:
                times[path] = value
        spread = max((by_key[f"B0.2/{p}/{curve}"].sides["Stereon"].spread
                      for p in times if f"B0.2/{p}/{curve}" in by_key), default=0.0)
        overhead = Gate(f"B0.2/DispatchOverhead/{curve}", f"VisitBuiltin over the local copy, {name}",
                        f"≤ {br.DISPATCH_OVERHEAD:.2f}×", "MISSING", spread=spread,
                        threshold=br.DISPATCH_OVERHEAD)
        if "Dispatch" in times and "DispatchCopy" in times:
            overhead.value = times["Dispatch"] / times["DispatchCopy"]
            overhead.stereon_ns = times["Dispatch"]
            overhead.status = judged(br.DISPATCH_OVERHEAD / overhead.value - 1.0, spread,
                                     overhead.value <= br.DISPATCH_OVERHEAD)
            if "Direct" in times:
                overhead.note = f"VisitBuiltin / direct {times['Dispatch'] / times['Direct']:.4f} (reported)"
        else:
            overhead.note = "; ".join(problems)
        fallback = Gate(f"B0.2/FallbackOverhead/{curve}", f"Virtual call over the direct call, {name}",
                        f"≤ {FALLBACK_NS:.0f} ns per call", "MISSING", spread=spread, threshold=FALLBACK_NS)
        if "Virtual" in times and "Direct" in times:
            fallback.value = times["Virtual"] - times["Direct"]
            fallback.status = judged((FALLBACK_NS - fallback.value) / FALLBACK_NS, spread,
                                     fallback.value <= FALLBACK_NS)
        else:
            fallback.note = "; ".join(problems)
        gates += [overhead, fallback]
    return gates


def counter_gates(by_key: dict[str, br.Row]) -> list[Gate]:
    """The memory and copying gates of Phase 0 (B0.4, B0.5, ADR-0004, ADR-0005), at the gate size."""
    gates: list[Gate] = []
    for key, name, label, target, check in br.COUNTER_ROWS:
        if check is None:
            continue
        kind, bound = check
        gate = Gate(f"{key}:{name}", label, target, "MISSING", threshold=bound)
        row = by_key.get(key)
        stereon = row.sides["Stereon"].counter(name) if row and "Stereon" in row.sides else None
        occt = row.sides["OCCT"].counter(name) if row and "OCCT" in row.sides else None
        if stereon is None:
            gate.note = "no Stereon measurement"
        elif not math.isfinite(stereon) or stereon < 0:
            gate.status, gate.note = "INVALID", "not a finite, non-negative measurement"
        elif kind == "ratio below":
            if occt is None or not (math.isfinite(occt) and occt > 0):
                gate.note = "no OCCT measurement"
            else:
                gate.value = stereon / occt
                gate.status = "PASS" if gate.value < bound else "FAIL"
        else:
            gate.value = stereon
            gate.status = "PASS" if stereon < bound else "FAIL"
        gates.append(gate)
    return gates


def convergence_gate(path: Path | None) -> Gate:
    """B1.3's convergence on a declared corpus, from its evidence file."""
    gate = Gate("B1.3/Convergence", "Projection converged on the declared corpus",
                f"≥ {CONVERGENCE_TARGET:.2%}", "MISSING", threshold=CONVERGENCE_TARGET)
    if path is None:
        gate.note = "no convergence evidence supplied (--convergence)"
        return gate
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
        converged, total = int(data["converged"]), int(data["cases"])
        invalid_successes = int(data["incorrect_success"])
    except (OSError, ValueError, KeyError, TypeError) as error:
        gate.status, gate.note = "INVALID", f"unreadable evidence: {error}"
        return gate
    if total <= 0:
        gate.status, gate.note = "INVALID", "an empty corpus"
        return gate
    gate.value = converged / total
    gate.note = (f"{converged} of {total} converged, {invalid_successes} incorrect successes, "
                 f"corpus {data.get('corpus', 'unnamed')}")
    gate.status = "PASS" if gate.value >= CONVERGENCE_TARGET and invalid_successes == 0 else "FAIL"
    return gate


def metadata_gate(machine: dict[str, str], candidate: str | None) -> Gate:
    """The run is tied to a commit, a comparator build and a host."""
    gate = Gate("Run/Metadata", "Commit, comparator build and host recorded", "all present", "PASS")
    missing = [key for key in REQUIRED_MACHINE if not machine.get(key)]
    if missing:
        gate.status, gate.note = "MISSING", f"machine record lacks {', '.join(missing)}"
    elif candidate and machine["commit"] != candidate:
        gate.status, gate.note = "FAIL", f"the run measured {machine['commit']}, not the candidate {candidate}"
    else:
        gate.note = f"{machine['commit']} against OCCT {machine['occt_tag']} ({machine['occt_commit'][:9]})"
    return gate


# --- The whole check -----------------------------------------------------------


def required_time_keys() -> list[str]:
    """Every row of the catalogue with a target, plus the B0.2 rows that must lead."""
    keys = [key for key, (_, target) in br.ROWS.items()
            if target != "report" and not key.startswith(("B0.2/DispatchCopy", "B0.2/Virtual", "B0.2/Direct"))]
    return sorted(set(keys))


def check(results: dict, machine: dict[str, str], candidate: str | None,
          convergence: Path | None) -> list[Gate]:
    rows = br.read_rows(results)
    by_key = {row.key: row for row in rows}
    gates = [metadata_gate(machine, candidate)]
    gates += [time_gate(by_key.get(key), key) for key in required_time_keys()]
    # Report-only rows that ran are listed too, so the table shows them as not gating.
    gates += [time_gate(row, row.key) for row in rows
              if br.ROWS.get(row.key, ("", ""))[1] == "report"
              and not row.key.startswith(br.COUNTER_ONLY_PREFIXES)]
    gates += dispatch_gates(by_key)
    gates += counter_gates(by_key)
    gates.append(convergence_gate(convergence))
    return gates


def released(gates: list[Gate]) -> bool:
    return all(g.status in PASSING for g in gates if g.status not in NON_GATING)


def as_json(gates: list[Gate], machine: dict[str, str], candidate: str | None) -> str:
    return json.dumps({
        "release": "v0.1",
        "candidate": candidate or machine.get("commit"),
        "machine": machine,
        "protocol": {"curve_repetitions": CURVE_REPETITIONS, "topology_repetitions": TOPOLOGY_REPETITIONS,
                     "max_spread": MAX_SPREAD, "parity_floor": PARITY_FLOOR},
        "verdict": "PASS" if released(gates) else "FAIL",
        "gates": [asdict(g) for g in gates],
    }, indent=2, sort_keys=True)


def as_markdown(gates: list[Gate], machine: dict[str, str]) -> str:
    def number(v: float | None) -> str:
        return "—" if v is None else f"{v:.6g}"
    lines = [
        "# v0.1 benchmark release gates",
        "",
        f"Commit {machine.get('commit', '?')}, OCCT {machine.get('occt_tag', '?')}; "
        f"verdict **{'PASS' if released(gates) else 'FAIL'}**. Values are full-precision; ratios are OCCT / Stereon "
        "except the dispatch rows (Stereon over Stereon) and the counters.",
        "",
        "| Gate | Measure | Target | Value | Threshold | Spread | Status | Note |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for g in gates:
        spread = "—" if g.spread is None else f"{g.spread:.1%}"
        lines.append(f"| {g.gate} | {g.label} | {g.target} | {number(g.value)} | {number(g.threshold)} "
                     f"| {spread} | {g.status} | {g.note} |")
    return "\n".join(lines) + "\n"


# --- Self-test -----------------------------------------------------------------


def _series(name: str, per_item_ns: float, count: int, **counters: float) -> list[dict]:
    """`count` repetitions of one benchmark at a fixed time per item, given as the iteration's
    time so that the fixture's value is exact (items per second would round it on the way back)."""
    return [br._entry(name, i, per_item_ns, None, **counters) for i in range(count)]


def _passing_run() -> dict:
    """Every required row at a comfortable pass, with the protocol's repetitions."""
    entries: list[dict] = []
    for key in required_time_keys():
        reps = expected_repetitions(key)
        threshold = br.MIN_RATIO.get(key, 1.0 if key in br.MUST_LEAD else PARITY_FLOOR)
        fit = {"fit_residual": 1e-12} if key in FIT_CHECKED else {}
        entries += _series(f"{key}/Stereon", 10.0, reps, **fit)
        entries += _series(f"{key}/OCCT", 10.0 * threshold * 1.5, reps, **fit)
    for curve in br.DISPATCH_CURVES:
        for path, ns in (("Direct", 10.0), ("DispatchCopy", 12.0), ("Dispatch", 12.0), ("Virtual", 13.0)):
            entries += _series(f"B0.2/{path}/{curve}/Stereon", ns, CURVE_REPETITIONS)
    for key, name, _, _, check in br.COUNTER_ROWS:
        if check is None:
            continue
        kind, bound = check
        stereon = bound * 0.5 * (1000.0 if kind == "ratio below" else 1.0)
        entries += _series(f"{key}/Stereon", 1.0, TOPOLOGY_REPETITIONS, **{name: stereon})
        if kind == "ratio below":
            entries += _series(f"{key}/OCCT", 1.0, TOPOLOGY_REPETITIONS, **{name: 1000.0})
    return {"benchmarks": entries}


def _replace(results: dict, key: str, library: str, per_item_ns: float, count: int | None = None) -> dict:
    """The run with one side of one row re-measured."""
    kept = [e for e in results["benchmarks"] if not e["name"].startswith(f"{key}/{library}")]
    reps = count if count is not None else expected_repetitions(key)
    return {"benchmarks": kept + _series(f"{key}/{library}", per_item_ns, reps)}


def _status(gates: list[Gate], gate: str) -> str:
    return next(g.status for g in gates if g.gate == gate)


def self_test() -> int:
    machine = {key: "x" for key in REQUIRED_MACHINE} | {"commit": "abc"}
    evidence = Path(__file__).with_name(".self-test-convergence.json")
    evidence.write_text(json.dumps({"converged": 10000, "cases": 10000, "incorrect_success": 0,
                                    "corpus": "self-test"}), encoding="utf-8")
    try:
        run = _passing_run()
        gates = check(run, machine, "abc", evidence)
        assert released(gates), [g for g in gates if g.status not in PASSING | NON_GATING]

        # Exactly at the thresholds: inclusive for "at least", strictly above 1 for "must lead".
        at = _replace(run, "B1.6/SurfaceProperties/bicubic", "OCCT", 10.0 * PARITY_FLOOR)
        at = _replace(at, "B1.1/CurveD0BatchFast/cubic", "OCCT", 20.0)
        gates = check(at, machine, "abc", evidence)
        assert _status(gates, "B1.6/SurfaceProperties/bicubic") == "PASS"
        assert _status(gates, "B1.1/CurveD0BatchFast/cubic") == "PASS"
        lead = _replace(run, "B0.2/Dispatch/line", "OCCT", 10.0)
        assert _status(check(lead, machine, "abc", evidence), "B0.2/Dispatch/line") == "FAIL"

        # Just below, where the report would print 1.98 or 0.90 and look like a pass.
        below = _replace(run, "B1.1/CurveD0BatchFast/cubic-rational", "OCCT", 19.84)
        below = _replace(below, "B1.6/SurfaceProperties/bicubic", "OCCT", 8.99)
        gates = check(below, machine, "abc", evidence)
        assert _status(gates, "B1.1/CurveD0BatchFast/cubic-rational") == "FAIL"
        assert _status(gates, "B1.6/SurfaceProperties/bicubic") == "FAIL"
        assert not released(gates)

        # The interpolation row counts only when both libraries interpolated.
        fit = [e for e in run["benchmarks"] if not e["name"].startswith("B1.4/Interpolation/")]
        unfit = fit + _series("B1.4/Interpolation/1000/Stereon", 10.0, 15, fit_residual=1e-12)
        unfit += _series("B1.4/Interpolation/1000/OCCT", 20.0, 15, fit_residual=1e-3)
        assert _status(check({"benchmarks": unfit}, machine, "abc", evidence), "B1.4/Interpolation/1000") == "FAIL"
        silent = fit + _series("B1.4/Interpolation/1000/Stereon", 10.0, 15)
        silent += _series("B1.4/Interpolation/1000/OCCT", 20.0, 15, fit_residual=1e-12)
        assert _status(check({"benchmarks": silent}, machine, "abc", evidence), "B1.4/Interpolation/1000") == "FAIL"

        # A row or a twin omitted; a timing that is not a positive finite number; too few repetitions.
        omitted = {"benchmarks": [e for e in run["benchmarks"]
                                  if not e["name"].startswith("B1.3/CurveProjection/")]}
        assert _status(check(omitted, machine, "abc", evidence), "B1.3/CurveProjection") == "MISSING"
        twinless = {"benchmarks": [e for e in run["benchmarks"]
                                   if not e["name"].startswith("B1.3/CurveProjection/OCCT")]}
        assert _status(check(twinless, machine, "abc", evidence), "B1.3/CurveProjection") == "MISSING"
        for bad in (math.nan, math.inf, 0.0, -1.0):
            broken = _replace(run, "B1.5/DegreeElevation", "Stereon", bad)
            assert _status(check(broken, machine, "abc", evidence), "B1.5/DegreeElevation") == "INVALID", bad
        short = _replace(run, "B1.5/DegreeElevation", "Stereon", 10.0, count=5)
        assert _status(check(short, machine, "abc", evidence), "B1.5/DegreeElevation") == "INVALID"

        # Report-only rows never gate, however slow; the approximation row is never a speed pass.
        slow = _series("B1.1/CurveD0Batch/cubic/Stereon", 100.0, 15) + _series("B1.1/CurveD0Batch/cubic/OCCT", 10.0, 15)
        slow += _series("B1.4/Approximation/1000/Stereon", 1.0, 15) + _series("B1.4/Approximation/1000/OCCT", 100.0, 15)
        gates = check({"benchmarks": run["benchmarks"] + slow}, machine, "abc", evidence)
        assert _status(gates, "B1.1/CurveD0Batch/cubic") == "REPORTED"
        assert _status(gates, "B1.4/Approximation/1000") == "NOT_COMPARABLE"
        assert released(gates)

        # A run from before ADR-0003's revision has no local-copy rows: incomplete, not a pass.
        old = {"benchmarks": [e for e in run["benchmarks"] if "/DispatchCopy/" not in e["name"]]}
        gates = check(old, machine, "abc", evidence)
        assert _status(gates, "B0.2/DispatchOverhead/bspline") == "MISSING"
        assert not released(gates)
        assert _status(gates, "B0.2/FallbackOverhead/bspline") == "PASS"

        # Excessive spread near the target is inconclusive; far from it, it is not.
        noisy = [e for e in run["benchmarks"] if not e["name"].startswith("B1.5/BezierDecomposition/")]
        noisy += _series("B1.5/BezierDecomposition/Stereon", 10.0, 15)
        # Median ratio 0.92, 2 % above the floor, with repetitions 24 % below the median.
        noisy += [br._entry("B1.5/BezierDecomposition/OCCT", i, 7.0 if i % 2 else 9.2, None)
                  for i in range(15)]
        assert _status(check({"benchmarks": noisy}, machine, "abc", evidence),
                       "B1.5/BezierDecomposition") == "INCONCLUSIVE"
        # The same spread far from the target leaves the verdict standing.
        far = [e for e in noisy if not e["name"].startswith("B1.5/BezierDecomposition/OCCT")]
        far += [br._entry("B1.5/BezierDecomposition/OCCT", i, 21.0 if i % 2 else 27.6, None)
                for i in range(15)]
        assert _status(check({"benchmarks": far}, machine, "abc", evidence),
                       "B1.5/BezierDecomposition") == "PASS"

        # Different results from the two libraries: not a speed comparison.
        differ = [e for e in run["benchmarks"] if not e["name"].startswith("B1.6/ArcLength/")]
        differ += _series("B1.6/ArcLength/rational-cubic/Stereon", 10.0, 15, length_mm=111.8)
        differ += _series("B1.6/ArcLength/rational-cubic/OCCT", 12.0, 15, length_mm=111.9)
        assert _status(check({"benchmarks": differ}, machine, "abc", evidence),
                       "B1.6/ArcLength/rational-cubic") == "FAIL"

        # The run must be the candidate's, with its comparator recorded; convergence evidence is required.
        assert _status(check(run, machine, "def", evidence), "Run/Metadata") == "FAIL"
        assert _status(check(run, {"commit": "abc"}, "abc", evidence), "Run/Metadata") == "MISSING"
        assert _status(check(run, machine, "abc", None), "B1.3/Convergence") == "MISSING"
        evidence.write_text(json.dumps({"converged": 9998, "cases": 10000, "incorrect_success": 0}),
                            encoding="utf-8")
        assert _status(check(run, machine, "abc", evidence), "B1.3/Convergence") == "FAIL"

        # The outputs: a deterministic JSON verdict and a table.
        gates = check(run, machine, "abc", None)
        document = json.loads(as_json(gates, machine, "abc"))
        assert document["verdict"] == "FAIL" and document["candidate"] == "abc"
        assert as_json(gates, machine, "abc") == as_json(gates, machine, "abc")
        assert "| B1.3/Convergence |" in as_markdown(gates, machine)
    finally:
        evidence.unlink(missing_ok=True)
    print("check_release_benchmarks.py: self-test passed")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("results", nargs="*", type=Path, help="Google Benchmark JSON of stereon_occt_bench")
    parser.add_argument("--machine", type=Path, help="the run's machine record (key=value lines)")
    parser.add_argument("--candidate", help="the commit the release is qualified for; the run must have measured it")
    parser.add_argument("--convergence", type=Path, help="B1.3's convergence evidence (JSON)")
    parser.add_argument("--json", type=Path, help="write the machine-readable verdict here")
    parser.add_argument("--out", type=Path, help="write the Markdown table here instead of standard output")
    parser.add_argument("--enforce", action="store_true", help="exit 1 unless every required gate passes")
    parser.add_argument("--self-test", action="store_true", help="check the checker")
    args = parser.parse_args()
    if args.self_test:
        return self_test()
    if not args.results:
        parser.error("results are required")
    results = br.merge_results([json.loads(path.read_text(encoding="utf-8")) for path in args.results])
    machine = br.read_machine(args.machine.read_text(encoding="utf-8")) if args.machine else {}
    gates = check(results, machine, args.candidate, args.convergence)
    if args.json:
        args.json.write_text(as_json(gates, machine, args.candidate), encoding="utf-8")
    table = as_markdown(gates, machine)
    if args.out:
        args.out.write_text(table, encoding="utf-8")
    else:
        sys.stdout.reconfigure(encoding="utf-8")  # the table has superscripts; consoles vary
        sys.stdout.write(table)
    if args.enforce and not released(gates):
        failing = [f"{g.gate} {g.status}" for g in gates if g.status not in PASSING | NON_GATING]
        print(f"release gates not met: {', '.join(failing)}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
