#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""Turn a run of stereon_occt_bench into the benchmark report of docs/benchmarks.md.

The harness (bench/occt/OcctBench.cpp, bench/occt/OcctTopoBench.cpp and
bench/TopoBench.cpp) names every benchmark <ID>/<operation>/<variant>/<library>,
with a Stereon and an OCCT twin on the same geometry; the topology rows (B0.3
to B0.5) take the sheet's face count as an argument, which Google Benchmark
appends, <ID>/<operation>/<library>/<faces>, and which stays in the row's key.
Run it with Google Benchmark's JSON output and repetitions,

    stereon_occt_bench --benchmark_repetitions=15 \\
        --benchmark_enable_random_interleaving=true \\
        --benchmark_out=results.json --benchmark_out_format=json

and this script pairs the twins and reports, per row, the median over the
repetitions of the CPU time per item (docs/performance.rst, Method), the
speed ratio OCCT / Stereon, the largest relative deviation of any repetition
from its median (the noise of the row), the target of docs/benchmarks.md and
whether it is met. A ratio within 10 % of 1 is parity, unless the row has a
ratio target of its own (B0.3: at least 5x). The memory and copying gates
(B0.4, B0.5, ADR-0005) are counters, not times, and have a table of their
own: the median of each counter over the repetitions, and its target.

Several results files (one per invocation, say with different repetitions)
are read as one run. The machine file is key=value lines, one per line,
written by the workflow that ran the benchmarks
(.github/workflows/occt-reference.yml) and rendered as the report's Setup
table, so a report never travels without its machine, compiler, OCCT build
and commit.

Usage: bench_report.py <results.json>... [--machine <machine.txt>] [--out <report.md>]
                       [--title <title>] | --self-test
"""

from __future__ import annotations

import argparse
import json
import statistics
import sys
from dataclasses import dataclass
from pathlib import Path

PARITY_BAND = 0.10  # a ratio within this of 1 is reported as parity
LIBRARIES = ("Stereon", "OCCT")

# Row labels and targets of the Phase 1 catalogue (docs/benchmarks.md), keyed by
# the benchmark name without its library. Unknown names are reported as they are.
ROWS: dict[str, tuple[str, str]] = {
    "B1.1/CurveD0/cubic": ("Curve D0, cubic, 24 poles", "≥ parity"),
    "B1.1/CurveD0/cubic-rational": ("Curve D0, rational cubic", "≥ parity"),
    "B1.1/CurveD0/degree7": ("Curve D0, degree 7", "≥ parity"),
    "B1.1/CurveD2/cubic": ("Curve D0–D2, cubic", "≥ parity"),
    "B1.1/CurveD2/cubic-rational": ("Curve D0–D2, rational cubic", "≥ parity"),
    "B1.1/SurfaceD0/bicubic": ("Surface D0, bicubic 12 × 12", "≥ parity"),
    "B1.1/SurfaceD0/bicubic-rational": ("Surface D0, rational bicubic", "≥ parity"),
    "B1.1/SurfaceD1/bicubic": ("Surface D0–D1, bicubic", "≥ parity"),
    "B1.3/CurveProjection": ("Curve projection, rational cubic", "≥ parity"),
    "B1.3/SurfaceProjection": ("Surface projection, bicubic", "≥ parity"),
    "B1.4/Interpolation/1000": ("Interpolation, 1000 points, cubic", "≥ parity"),
    "B1.4/Approximation/1000": ("Approximation, 1000 points, 1e-3 mm", "≥ parity at equal tolerance"),
    "B1.5/KnotInsertion/100": ("100 knot insertions", "≥ parity"),
    "B1.5/BezierDecomposition": ("Bezier decomposition, 21 segments", "≥ parity"),
    "B1.6/CurveProperties/rational-cubic": ("Frenet frame and curvature, rational cubic", "≥ parity"),
    "B1.6/SurfaceProperties/bicubic": ("Principal curvatures and directions, bicubic", "≥ parity"),
    "B1.6/ArcLength/rational-cubic": ("Arc length of the whole curve, rational cubic", "≥ parity at equal accuracy"),
    "B1.6/PointAtLength/64": ("Parameter at a given length", "≥ parity"),
}

# The topology rows (B0.3 to B0.5) at both sheet sizes: 10^4 faces for a quick
# look, 10^6 for the body the gates of ADR-0004 and ADR-0005 name. A target of
# "report" is informational: the row has no gate of its own.
TOPOLOGY_SIZES = {"10000": "10⁴ faces", "1000000": "10⁶ faces"}
GATE_SIZE = "1000000"
for _size, _faces in TOPOLOGY_SIZES.items():
    _gate = _size == GATE_SIZE
    ROWS.update({
        f"B0.3/SheetBuild/{_size}": (f"Build the sheet, {_faces}", "report"),
        f"B0.3/FacesOfEdges/{_size}": (f"Faces of every edge, {_faces}", "≥ 5×" if _gate else "report"),
        f"B0.3/EdgesOfVertices/{_size}": (f"Edges of every vertex, {_faces}", "≥ 5×" if _gate else "report"),
        f"B0.3/AncestorMap/{_size}": (f"OCCT's ancestor map alone, {_faces}", "report"),
        f"B0.5/ReplaceSurface/{_size}": (f"Replace one face's surface, {_faces}", "report"),
        f"B0.5/DeleteFace/{_size}": (f"Delete one face, {_faces}", "report"),
        f"ADR-0005/EditSequence/{_size}": (f"100 edits without Compact(), {_faces}", "report"),
        f"ADR-0004/Compact/{_size}": (f"Compact() after 50 deletions, {_faces}", "report"),
    })

# Rows whose target is a ratio of their own rather than parity: OCCT / Stereon
# must reach it (docs/benchmarks.md, B0.3: adjacency at least 5x faster).
MIN_RATIO = {
    f"B0.3/FacesOfEdges/{GATE_SIZE}": 5.0,
    f"B0.3/EdgesOfVertices/{GATE_SIZE}": 5.0,
}

# Rows whose timed loop is not what they measure (B0.4 times the build; the
# measure is a counter), left out of the time table.
COUNTER_ONLY_PREFIXES = ("B0.4/",)

# Counters the report keeps, besides lengths (*_mm), which the notes show.
REPORTED_COUNTERS = {
    "heap_bytes_per_face", "heap_bytes_per_edge", "bytes_per_face", "bytes_per_edge",
    "copied_fraction", "pages_copied", "pages", "memory_ratio", "heap_delta_bytes",
}

# The memory and copying gates, one line each: (row key, counter, label, target,
# check). A check is ("ratio below", r): Stereon's value over OCCT's below r;
# ("below", x): Stereon's value below x; or None, reported without a gate.
COUNTER_ROWS: list[tuple[str, str, str, str, tuple[str, float] | None]] = []
for _size, _faces in TOPOLOGY_SIZES.items():
    _gate = _size == GATE_SIZE
    COUNTER_ROWS += [
        (f"B0.4/HeapPerFace/{_size}", "heap_bytes_per_face", f"Live heap per face, {_faces}",
         "< 50% of OCCT" if _gate else "report", ("ratio below", 0.5) if _gate else None),
        (f"B0.4/HeapPerFace/{_size}", "heap_bytes_per_edge", f"Live heap per edge, {_faces}", "report", None),
        (f"B0.4/BytesPerFace/{_size}", "bytes_per_face", f"Table pages per face, {_faces}", "report", None),
        (f"B0.5/ReplaceSurface/{_size}", "copied_fraction", f"Pages copied, one surface replaced, {_faces}",
         "< 1% of body" if _gate else "report", ("below", 0.01) if _gate else None),
        (f"B0.5/ReplaceSurface/{_size}", "heap_delta_bytes", f"Heap added, one surface replaced, {_faces}",
         "report", None),
        (f"B0.5/DeleteFace/{_size}", "copied_fraction", f"Pages copied, one face deleted, {_faces}",
         "< 1% of body" if _gate else "report", ("below", 0.01) if _gate else None),
        (f"B0.5/DeleteFace/{_size}", "heap_delta_bytes", f"Heap added, one face deleted, {_faces}",
         "report", None),
        (f"ADR-0005/EditSequence/{_size}", "memory_ratio",
         f"Body and its edited copy over the body, 100 edits, {_faces}",
         "< 1.2×" if _gate else "report", ("below", 1.2) if _gate else None),
    ]

# Rows whose two libraries do not compute the same thing, so the ratio is
# indicative only (bench/reports/phase1-2026-10-15.md, Notes).
NOT_LIKE_FOR_LIKE = {"B1.4/Approximation/1000"}

# Rows quoted per iteration (one fit, one batch of insertions) rather than per
# counted item, as docs/performance.rst and the earlier reports quote them. The
# harness counts the points of a fit as items, which would give a time per point.
PER_ITERATION = {"B1.4/Interpolation/1000", "B1.4/Approximation/1000", "B1.5/KnotInsertion/100"}

# Friendly labels for the machine file's keys, in the order the Setup table shows them.
SETUP_LABELS = (
    ("date", "Date"),
    ("commit", "Commit"),
    ("ref", "Ref"),
    ("workflow_run", "Workflow run"),
    ("runner", "Runner"),
    ("cpu", "CPU"),
    ("cpus", "CPUs"),
    ("memory", "Memory"),
    ("kernel", "Kernel"),
    ("compiler", "Compiler"),
    ("occt_tag", "OCCT tag"),
    ("occt_commit", "OCCT commit"),
    ("occt_compiler", "OCCT compiler"),
    ("occt_build_type", "OCCT build type"),
    ("repetitions", "Repetitions"),
)

TIME_UNITS_IN_NS = {"ns": 1.0, "us": 1e3, "ms": 1e6, "s": 1e9}


@dataclass
class Side:
    """One library's repetitions of one row."""

    per_item_ns: list[float]
    counters: dict[str, list[float]]

    @property
    def median(self) -> float:
        return statistics.median(self.per_item_ns)

    def counter(self, name: str) -> float | None:
        """The median of a counter over the repetitions, or None if the row has none."""
        values = self.counters.get(name)
        return statistics.median(values) if values else None

    @property
    def spread(self) -> float:
        """Largest relative deviation of a repetition from the median."""
        m = self.median
        return max(abs(t - m) / m for t in self.per_item_ns) if m > 0 else 0.0


@dataclass
class Row:
    key: str
    sides: dict[str, Side]

    @property
    def ratio(self) -> float | None:
        if any(lib not in self.sides for lib in LIBRARIES):
            return None
        return self.sides["OCCT"].median / self.sides["Stereon"].median

    @property
    def spread(self) -> float:
        return max(side.spread for side in self.sides.values())


# --- Reading Google Benchmark output -----------------------------------------


def split_name(name: str) -> tuple[str, str] | None:
    """'B1.1/CurveD0/cubic/Stereon' -> ('B1.1/CurveD0/cubic', 'Stereon'); None if not a twin.

    An argument Google Benchmark appends after the library stays in the key:
    'B0.3/FacesOfEdges/OCCT/1000000' -> ('B0.3/FacesOfEdges/1000000', 'OCCT')."""
    key, _, library = name.rpartition("/")
    argument = ""
    if library.isdigit() and key:
        argument = library
        key, _, library = key.rpartition("/")
    if library not in LIBRARIES or not key:
        return None
    return (f"{key}/{argument}" if argument else key), library


def per_item_ns(entry: dict, key: str) -> float:
    """CPU time per item of one repetition, in nanoseconds: per counted item, or per iteration where the
    row is quoted that way or nothing was counted."""
    items_per_second = entry.get("items_per_second")
    if items_per_second and key not in PER_ITERATION:
        return 1e9 / float(items_per_second)
    return float(entry["cpu_time"]) * TIME_UNITS_IN_NS[entry.get("time_unit", "ns")]


def read_rows(results: dict) -> list[Row]:
    """Pair the repetitions of every Stereon benchmark with its OCCT twin.

    Rows come back sorted by key, so the report reads by benchmark ID and two
    runs diff line by line, whatever order the interleaving ran them in."""
    rows: dict[str, Row] = {}
    for entry in results.get("benchmarks", []):
        # With repetitions the file holds every repetition and the aggregates;
        # the medians below are ours, so aggregates are skipped.
        if entry.get("run_type", "iteration") != "iteration":
            continue
        if "error_occurred" in entry:
            continue
        parts = split_name(entry["name"])
        if parts is None:
            continue
        key, library = parts
        row = rows.setdefault(key, Row(key, {}))
        side = row.sides.setdefault(library, Side([], {}))
        side.per_item_ns.append(per_item_ns(entry, key))
        for name, value in entry.items():
            if (name.endswith("_mm") or name in REPORTED_COUNTERS) and isinstance(value, (int, float)):
                side.counters.setdefault(name, []).append(float(value))
    return [rows[key] for key in sorted(rows)]


def merge_results(files: list[dict]) -> dict:
    """Several runs' results as one: their benchmark entries together."""
    return {"benchmarks": [entry for results in files for entry in results.get("benchmarks", [])]}


def read_machine(text: str) -> dict[str, str]:
    machine: dict[str, str] = {}
    for line in text.splitlines():
        key, sep, value = line.partition("=")
        if sep and key.strip():
            machine[key.strip()] = value.strip()
    return machine


# --- Formatting ----------------------------------------------------------------


def format_time(ns: float) -> str:
    """Three significant digits in the unit that keeps the number between 1 and 1000."""
    for unit, scale in (("ns", 1.0), ("µs", 1e3), ("ms", 1e6)):
        value = ns / scale
        text = f"{value:.3g}"
        if "e" not in text or unit == "ms":  # "e" means the rounded value reached 1000
            break
    if "e" in text:
        text = f"{value:.0f}"
    return f"{text} {unit}"


def format_ratio(ratio: float) -> str:
    """Two decimals near parity, where the verdict turns on the second one; fewer as the lead grows."""
    if ratio < 2:
        return f"{ratio:.2f}"
    return f"{ratio:.1f}" if ratio < 10 else f"{ratio:.0f}"


def verdict(row: Row) -> tuple[str, str]:
    """(ratio column, met column) for a row."""
    informational = ROWS.get(row.key, ("", ""))[1] == "report"
    ratio = row.ratio
    if ratio is None:
        if informational:
            return "—", "—"
        missing = [lib for lib in LIBRARIES if lib not in row.sides]
        return "—", f"no {', '.join(missing)} run"
    text = format_ratio(ratio)
    if row.key in NOT_LIKE_FOR_LIKE:
        return text, "see note (not like for like)"
    if row.key in MIN_RATIO:
        return text, "yes" if ratio >= MIN_RATIO[row.key] else "no"
    if informational:
        return text, "—"
    if abs(ratio - 1) <= PARITY_BAND:
        return text, "parity"
    return text, "yes" if ratio > 1 else "no"


def format_counter(name: str, value: float) -> str:
    """A counter as the report quotes it: a fraction in percent, bytes whole, a ratio with its sign."""
    if name == "copied_fraction":
        return f"{value:.3%}"
    if name == "memory_ratio":
        return f"{value:.3f}×"
    return f"{value:,.0f} B" if "bytes" in name else f"{value:.6g}"


def counter_rows(rows: list[Row]) -> list[str]:
    """The lines of the memory and copying table, for the gates whose rows ran."""
    by_key = {row.key: row for row in rows}
    out = []
    for key, name, label, target, check in COUNTER_ROWS:
        row = by_key.get(key)
        if row is None:
            continue
        values = {lib: row.sides[lib].counter(name) if lib in row.sides else None for lib in LIBRARIES}
        stereon, occt = values["Stereon"], values["OCCT"]
        if stereon is None:
            continue
        relative = stereon / occt if occt else None
        met = "—"
        if check is not None:
            kind, bound = check
            if kind == "ratio below":
                met = "no OCCT run" if relative is None else ("yes" if relative < bound else "no")
            else:
                met = "yes" if stereon < bound else "no"
        cells = [format_counter(name, v) if v is not None else "—" for v in (stereon, occt)]
        share = f"{relative:.0%}" if relative is not None else "—"
        out.append(f"| {key.split('/', 1)[0]} | {label} | {cells[0]} | {cells[1]} | {share} | {target} | {met} |")
    return out


def row_label(key: str) -> tuple[str, str, str]:
    """(ID, benchmark label, target) for a row key."""
    label, target = ROWS.get(key, (key.split("/", 1)[-1], "≥ parity"))
    return key.split("/", 1)[0], label, target


def render(rows: list[Row], machine: dict[str, str], title: str) -> str:
    out = [f"# {title}", ""]
    if machine:
        out += ["## Setup", "", "| | |", "| --- | --- |"]
        shown = set()
        for key, label in SETUP_LABELS:
            if key in machine:
                out.append(f"| {label} | {machine[key]} |")
                shown.add(key)
        for key in sorted(set(machine) - shown):
            out.append(f"| {key} | {machine[key]} |")
        out.append("")
    out += [
        "## Results",
        "",
        "Median CPU time over the repetitions, both libraries interleaved: per evaluated parameter, projected "
        "point, property or length query; per fit (B1.4), per 100 insertions and per decomposition (B1.5); "
        "per build, per query over the whole sheet or per edit (B0.3, B0.5, ADR-0004, ADR-0005). "
        "Ratio = OCCT / Stereon, so > 1 means Stereon is faster; within 10 % of 1 is parity, and a row with a "
        "ratio target of its own (B0.3: ≥ 5×) is met when it reaches it. A target of \"report\" is not a gate. "
        "Spread is the largest relative deviation of any repetition from its median, over both libraries.",
        "",
        "| ID | Benchmark | Stereon | OCCT | Ratio | Spread | Target | Met |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for row in rows:
        if row.key.startswith(COUNTER_ONLY_PREFIXES):
            continue
        ident, label, target = row_label(row.key)
        ratio, met = verdict(row)
        times = [format_time(row.sides[lib].median) if lib in row.sides else "—" for lib in LIBRARIES]
        out.append(f"| {ident} | {label} | {times[0]} | {times[1]} | {ratio} | {row.spread:.0%} | {target} | {met} |")
    memory = counter_rows(rows)
    if memory:
        out += [
            "",
            "## Memory and copying",
            "",
            "Median over the repetitions of each counter. Live heap is the bytes the sheet holds in live "
            "allocations with the shape alive, measured the same way for both libraries (bench/Memory.hpp); "
            "pages copied are the pages of the edited snapshot not shared with the body (ADR-0004, ADR-0005). "
            "Stereon / OCCT below 100 % means Stereon uses less.",
            "",
            "| ID | Measure | Stereon | OCCT | Stereon / OCCT | Target | Met |",
            "| --- | --- | --- | --- | --- | --- | --- |",
            *memory,
        ]
    notes = []
    for row in rows:
        counters = {lib: side.counters for lib, side in row.sides.items() if side.counters}
        for name in sorted({n for c in counters.values() for n in c if n.endswith("_mm")}):
            values = ", ".join(
                f"{statistics.median(c[name]):.17g} ({lib})" for lib, c in counters.items() if name in c
            )
            notes.append(f"- **{row.key}**, `{name}`: {values}.")
    if NOT_LIKE_FOR_LIKE & {row.key for row in rows}:
        notes.append(
            "- **Approximation is not like for like.** Stereon returns the fewest-pole least-squares cubic within "
            "the budget; OCCT's `GeomAPI_PointsToBSpline` also optimises parameters and continuity. "
            "Read the ratio as indicative."
        )
    if notes:
        out += ["", "## Notes", "", *notes]
    out.append("")
    return "\n".join(out)


# --- Self-test -----------------------------------------------------------------


def _entry(name: str, index: int, cpu_ns: float, items: int | None, **counters: float) -> dict:
    entry = {
        "name": name,
        "run_name": name,
        "run_type": "iteration",
        "repetitions": 3,
        "repetition_index": index,
        "threads": 1,
        "iterations": 1000,
        "real_time": cpu_ns,
        "cpu_time": cpu_ns,
        "time_unit": "ns",
    }
    if items is not None:
        entry["items_per_second"] = items / (cpu_ns * 1e-9)
    entry.update(counters)
    return entry


def _aggregate(name: str, cpu_ns: float) -> dict:
    return {"name": f"{name}_median", "run_name": name, "run_type": "aggregate", "aggregate_name": "median",
            "cpu_time": cpu_ns, "time_unit": "ns"}


def self_test() -> int:
    n = 4096
    benchmarks = [
        # A 2x lead, median taken from three repetitions, the aggregate ignored.
        _entry("B1.1/CurveD0/cubic/Stereon", 0, 100 * n, n),
        _entry("B1.1/CurveD0/cubic/Stereon", 1, 90 * n, n),
        _entry("B1.1/CurveD0/cubic/Stereon", 2, 110 * n, n),
        _aggregate("B1.1/CurveD0/cubic/Stereon", 100 * n),
        _entry("B1.1/CurveD0/cubic/OCCT", 0, 200 * n, n),
        _entry("B1.1/CurveD0/cubic/OCCT", 1, 190 * n, n),
        _entry("B1.1/CurveD0/cubic/OCCT", 2, 500 * n, n),  # a stray slow repetition; the median ignores it
        # Parity, one item per iteration, no items counter.
        _entry("B1.5/BezierDecomposition/Stereon", 0, 22_900, None),
        _entry("B1.5/BezierDecomposition/OCCT", 0, 24_000, None),
        # OCCT ahead, with a counter that the notes report.
        _entry("B1.6/ArcLength/rational-cubic/Stereon", 0, 200_000, 1, length_mm=111.80161013322777),
        _entry("B1.6/ArcLength/rational-cubic/OCCT", 0, 150_000, 1, length_mm=111.80161013322775),
        # Not like for like, in milliseconds, quoted per fit although the harness counts 1000 points.
        _entry("B1.4/Approximation/1000/Stereon", 0, 3.8e6, 1000),
        _entry("B1.4/Approximation/1000/OCCT", 0, 141e6, 1000),
        # An unknown twin and a benchmark without a twin.
        _entry("B9.9/Mystery/Stereon", 0, 10, None),
        _entry("B9.9/Mystery/OCCT", 0, 10, None),
        _entry("B1.1/Lonely/Stereon", 0, 10, None),
        _entry("Unpaired/thing", 0, 10, None),
    ]
    rows = read_rows({"benchmarks": benchmarks})
    assert [row.key for row in rows] == sorted(row.key for row in rows), "rows are sorted by key"
    by_key = {row.key: row for row in rows}

    cubic = by_key["B1.1/CurveD0/cubic"]
    assert len(cubic.sides["Stereon"].per_item_ns) == 3, "aggregates must be skipped"
    assert abs(cubic.sides["Stereon"].median - 100) < 1e-9
    assert abs(cubic.sides["OCCT"].median - 200) < 1e-9
    assert verdict(cubic) == ("2.0", "yes"), verdict(cubic)
    assert abs(cubic.spread - 1.5) < 1e-9, cubic.spread  # the 500 ns repetition against a 200 ns median

    bezier = by_key["B1.5/BezierDecomposition"]
    assert format_time(bezier.sides["Stereon"].median) == "22.9 µs"
    assert verdict(bezier) == ("1.05", "parity"), verdict(bezier)

    length = by_key["B1.6/ArcLength/rational-cubic"]
    assert verdict(length) == ("0.75", "no"), verdict(length)
    assert length.sides["OCCT"].counters == {"length_mm": [111.80161013322775]}

    approx = by_key["B1.4/Approximation/1000"]
    assert verdict(approx) == ("37", "see note (not like for like)"), verdict(approx)
    assert format_time(approx.sides["Stereon"].median) == "3.8 ms"
    assert format_time(approx.sides["OCCT"].median) == "141 ms"

    assert verdict(by_key["B1.1/Lonely"]) == ("—", "no OCCT run")
    assert "Unpaired/thing" not in by_key and "Unpaired" not in by_key
    assert row_label("B9.9/Mystery") == ("B9.9", "Mystery", "≥ parity")

    assert format_time(94.4) == "94.4 ns" and format_time(999.6) == "1 µs" and format_time(0.34e6) == "340 µs"
    assert format_ratio(0.88) == "0.88" and format_ratio(1.15) == "1.15" and format_ratio(1.96) == "1.96"
    assert format_ratio(2.04) == "2.0" and format_ratio(3.44) == "3.4" and format_ratio(36.6) == "37"

    machine = read_machine("cpu=Intel Core i5-10210U\ncpus=4\nno equals sign here\nocct_tag=V8_0_1\n=empty key\n")
    assert machine == {"cpu": "Intel Core i5-10210U", "cpus": "4", "occt_tag": "V8_0_1"}

    report = render(rows, machine, "Test report")
    assert "| CPU | Intel Core i5-10210U |" in report
    assert "| B1.1 | Curve D0, cubic, 24 poles | 100 ns | 200 ns | 2.0 | 150% | ≥ parity | yes |" in report
    assert "`length_mm`: 111.80161013322777 (Stereon), 111.80161013322775 (OCCT)" in report
    assert "not like for like" in report
    _self_test_topology()
    print("bench_report.py: self-test passed")
    return 0


def _self_test_topology() -> None:
    """The B0 rows: the size Google Benchmark appends, the 5x target, the counter gates, several files."""
    assert split_name("B0.3/FacesOfEdges/OCCT/1000000") == ("B0.3/FacesOfEdges/1000000", "OCCT")
    assert split_name("B1.4/Interpolation/1000/Stereon") == ("B1.4/Interpolation/1000", "Stereon")
    assert split_name("B0.3/FacesOfEdges/1000000") is None
    ms = 1e6
    first = {"benchmarks": [
        # 6x on one query (met), 3x on the other (not met).
        _entry("B0.3/FacesOfEdges/Stereon/1000000", 0, 400 * ms, None),
        _entry("B0.3/FacesOfEdges/OCCT/1000000", 0, 2400 * ms, None),
        _entry("B0.3/EdgesOfVertices/Stereon/1000000", 0, 1000 * ms, None),
        _entry("B0.3/EdgesOfVertices/OCCT/1000000", 0, 3000 * ms, None),
        # OCCT only, informational.
        _entry("B0.3/AncestorMap/OCCT/1000000", 0, 2300 * ms, None),
    ]}
    second = {"benchmarks": [
        # Memory: the counters' medians, 25 % of OCCT; the build time is not reported.
        _entry("B0.4/HeapPerFace/Stereon/1000000", 0, 900 * ms, None, heap_bytes_per_face=100.0),
        _entry("B0.4/HeapPerFace/Stereon/1000000", 1, 900 * ms, None, heap_bytes_per_face=120.0),
        _entry("B0.4/HeapPerFace/Stereon/1000000", 2, 900 * ms, None, heap_bytes_per_face=90.0),
        _entry("B0.4/HeapPerFace/OCCT/1000000", 0, 15000 * ms, None, heap_bytes_per_face=400.0),
        # Copying: 0.2 % of the body (met); the 100 edits at 1.3x (not met).
        _entry("B0.5/DeleteFace/Stereon/1000000", 0, 50 * ms, None, copied_fraction=0.002),
        _entry("ADR-0005/EditSequence/Stereon/1000000", 0, 9000 * ms, None, memory_ratio=1.3),
    ]}
    rows = read_rows(merge_results([first, second]))
    by_key = {row.key: row for row in rows}
    assert verdict(by_key["B0.3/FacesOfEdges/1000000"]) == ("6.0", "yes")
    assert verdict(by_key["B0.3/EdgesOfVertices/1000000"]) == ("3.0", "no")
    assert verdict(by_key["B0.3/AncestorMap/1000000"]) == ("—", "—")
    assert verdict(by_key["ADR-0005/EditSequence/1000000"]) == ("—", "—")
    assert by_key["B0.4/HeapPerFace/1000000"].sides["Stereon"].counter("heap_bytes_per_face") == 100.0
    memory = counter_rows(rows)
    assert "| B0.4 | Live heap per face, 10⁶ faces | 100 B | 400 B | 25% | < 50% of OCCT | yes |" in memory
    assert "| B0.5 | Pages copied, one face deleted, 10⁶ faces | 0.200% | — | — | < 1% of body | yes |" in memory
    assert any("1.300×" in line and line.endswith("| < 1.2× | no |") for line in memory), memory
    report = render(rows, {}, "Topology")
    assert "## Memory and copying" in report
    assert "Faces of every edge, 10⁶ faces" in report
    assert "Build the sheet" not in report  # no SheetBuild entries in this run
    assert "| B0.4 | " not in report.split("## Memory and copying")[0], "B0.4 build times stay out"


# --- Main ------------------------------------------------------------------------


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("results", nargs="*", type=Path,
                        help="Google Benchmark JSON output of stereon_occt_bench, one file per invocation")
    parser.add_argument("--machine", type=Path, help="key=value lines describing the machine and the builds")
    parser.add_argument("--out", type=Path, help="write the Markdown report here instead of standard output")
    parser.add_argument("--title", default="Phase 1 benchmark report — curves and surfaces (B1.1–B1.6)")
    parser.add_argument("--self-test", action="store_true", help="run the built-in checks and exit")
    args = parser.parse_args()

    if args.self_test:
        return self_test()
    if not args.results:
        parser.error("a results file is required (or --self-test)")

    results = merge_results([json.loads(path.read_text(encoding="utf-8")) for path in args.results])
    machine = read_machine(args.machine.read_text(encoding="utf-8")) if args.machine else {}
    rows = read_rows(results)
    if not rows:
        names = ", ".join(str(path) for path in args.results)
        print(f"bench_report.py: no Stereon/OCCT twins in {names}", file=sys.stderr)
        return 1
    repetitions = {len(side.per_item_ns) for row in rows for side in row.sides.values()}
    machine.setdefault("repetitions", ", ".join(str(r) for r in sorted(repetitions)))
    report = render(rows, machine, args.title)
    if args.out:
        args.out.write_text(report, encoding="utf-8")
    else:
        sys.stdout.buffer.write(report.encode("utf-8"))  # the report is UTF-8 whatever the console is
    return 0


if __name__ == "__main__":
    sys.exit(main())
