#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""Keep the trend of Stereon's own benchmarks across the commits on main.

docs/benchmarks.md asks for results stored per commit and plotted as trends,
and for a regression of more than 10 % in any Stereon benchmark between
commits to be explained. The Benchmarks workflow (.github/workflows/bench.yml)
runs stereon_bench on every push to main and keeps the runs on the
`bench-data` branch, one small file per run, through three commands:

    record   one run of Google Benchmark JSON (with repetitions) to a run file:
             per benchmark, the median CPU time per item over the repetitions
             (per iteration where nothing is counted) and its spread, the
             largest relative deviation of a repetition from that median;
             with the commit, date and machine from a key=value machine file,
             as tools/bench_report.py reads it.
    compare  the latest run file against the one before: every benchmark
             slower by more than the threshold (10 %) is a regression, listed
             in Markdown and, with --annotate, as a GitHub warning. A change
             of runner CPU is called out, since hosted runners differ and a
             comparison across two CPUs is only indicative. Regressions never
             fail the job: timings on a shared runner cannot gate a merge.
    render   every run file to one static HTML page, a small chart per
             benchmark over the last runs, which the Docs workflow publishes
             with the manual.

Run files are named so that they sort in the order the runs were made
(<UTC time>-<commit>.json), and nothing else orders them.

Usage: bench_trend.py record <results.json> --machine <machine.txt> --out <run.json>
       bench_trend.py compare <runs dir> [--threshold 0.10] [--out <summary.md>] [--annotate]
       bench_trend.py render <runs dir> --out <index.html> [--last 60] [--repo-url <url>]
       bench_trend.py --self-test
"""

from __future__ import annotations

import argparse
import html
import json
import statistics
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

SCHEMA = 1
THRESHOLD = 0.10  # docs/benchmarks.md: a regression of more than 10 % is explained
LAST_RUNS = 60  # runs per chart
TIME_UNITS_IN_NS = {"ns": 1.0, "us": 1e3, "ms": 1e6, "s": 1e9}


# --- record ------------------------------------------------------------------------


def per_item_ns(entry: dict) -> float:
    """CPU time per counted item of one repetition, or per iteration where nothing is counted."""
    items_per_second = entry.get("items_per_second")
    if items_per_second:
        return 1e9 / float(items_per_second)
    return float(entry["cpu_time"]) * TIME_UNITS_IN_NS[entry.get("time_unit", "ns")]


def read_machine(text: str) -> dict[str, str]:
    machine: dict[str, str] = {}
    for line in text.splitlines():
        key, sep, value = line.partition("=")
        if sep and key.strip():
            machine[key.strip()] = value.strip()
    return machine


def summarise(results: dict) -> dict[str, dict]:
    """Per benchmark: the median time per item over its repetitions, the spread and the count."""
    times: dict[str, list[float]] = {}
    for entry in results.get("benchmarks", []):
        # The medians are ours; aggregates and failed runs are skipped.
        if entry.get("run_type", "iteration") != "iteration" or "error_occurred" in entry:
            continue
        times.setdefault(entry["name"], []).append(per_item_ns(entry))
    summary = {}
    for name in sorted(times):
        values = times[name]
        median = statistics.median(values)
        spread = max(abs(v - median) / median for v in values) if median > 0 else 0.0
        summary[name] = {"ns": median, "spread": spread, "repetitions": len(values)}
    return summary


def record(results: dict, machine: dict[str, str]) -> dict:
    return {
        "schema": SCHEMA,
        "commit": machine.get("commit", ""),
        "date": machine.get("date", ""),
        "cpu": machine.get("cpu", ""),
        "machine": machine,
        "benchmarks": summarise(results),
    }


# --- compare -----------------------------------------------------------------------


@dataclass
class Change:
    name: str
    before: float
    after: float
    spread: float  # the larger of the two runs' spreads

    @property
    def relative(self) -> float:
        return self.after / self.before - 1


def load_runs(directory: Path) -> list[dict]:
    """Every run file in the directory, oldest first."""
    runs = []
    for path in sorted(directory.glob("*.json")):
        run = json.loads(path.read_text(encoding="utf-8"))
        if run.get("schema") == SCHEMA:
            runs.append(run)
    return runs


def changes(before: dict, after: dict) -> list[Change]:
    """The benchmarks the two runs share, slowest change first."""
    out = []
    for name, now in after["benchmarks"].items():
        then = before["benchmarks"].get(name)
        if then is None or then["ns"] <= 0:
            continue
        out.append(Change(name, then["ns"], now["ns"], max(then["spread"], now["spread"])))
    return sorted(out, key=lambda c: c.relative, reverse=True)


def format_time(ns: float) -> str:
    """Three significant digits in the unit that keeps the number between 1 and 1000."""
    for unit, scale in (("ns", 1.0), ("µs", 1e3), ("ms", 1e6), ("s", 1e9)):
        if ns / scale < 999.5 or unit == "s":
            return f"{ns / scale:.3g} {unit}"
    raise AssertionError("unreachable")


def short(commit: str) -> str:
    return commit[:7] if commit else "?"


def compare_report(runs: list[dict], threshold: float) -> tuple[str, list[Change]]:
    """(Markdown summary, regressions) of the latest run against the one before."""
    if len(runs) < 2:
        return "No earlier run on record; nothing to compare yet.\n", []
    before, after = runs[-2], runs[-1]
    all_changes = changes(before, after)
    regressions = [c for c in all_changes if c.relative > threshold]
    improvements = [c for c in all_changes if c.relative < -threshold]
    out = [f"Run of `{short(after['commit'])}` against `{short(before['commit'])}`: "
           f"{len(all_changes)} benchmarks compared, {len(regressions)} slower and "
           f"{len(improvements)} faster by more than {threshold:.0%}.", ""]
    if before.get("cpu") != after.get("cpu"):
        out += [f"**The runner's CPU changed** ({before.get('cpu') or '?'} → {after.get('cpu') or '?'}), "
                "so the comparison is indicative only.", ""]
    for title, rows in (("Slower", regressions), ("Faster", improvements)):
        if rows:
            out += [f"### {title}", "", "| Benchmark | Before | After | Change | Spread |",
                    "| --- | --- | --- | --- | --- |"]
            out += [f"| `{c.name}` | {format_time(c.before)} | {format_time(c.after)} | {c.relative:+.1%} "
                    f"| {c.spread:.0%} |" for c in rows]
            out.append("")
    return "\n".join(out), regressions


# --- render ------------------------------------------------------------------------

CHART_WIDTH, CHART_HEIGHT, PAD = 300, 90, 6

PAGE_STYLE = """
:root { --bg: #ffffff; --fg: #1f2328; --muted: #656d76; --line: #0969da; --grid: #d0d7de;
        --bad: #cf222e; --card: #f6f8fa; }
@media (prefers-color-scheme: dark) {
  :root { --bg: #0d1117; --fg: #e6edf3; --muted: #8d96a0; --line: #4493f8; --grid: #30363d;
          --bad: #f85149; --card: #161b22; }
}
body { background: var(--bg); color: var(--fg); margin: 0 auto; padding: 16px; max-width: 1280px;
       font: 14px/1.45 system-ui, -apple-system, "Segoe UI", sans-serif; }
h1 { font-size: 1.5rem; margin: 0 0 4px; } h2 { font-size: 1.15rem; margin: 28px 0 8px; }
p { color: var(--muted); margin: 0 0 8px; }
.grid { display: grid; gap: 12px; grid-template-columns: repeat(auto-fill, minmax(300px, 1fr)); }
.card { background: var(--card); border: 1px solid var(--grid); border-radius: 6px; padding: 8px 10px; }
.name { font: 12px ui-monospace, SFMono-Regular, Consolas, monospace; overflow-wrap: anywhere; }
.now { font-variant-numeric: tabular-nums; } .worse { color: var(--bad); } .muted { color: var(--muted); }
svg { display: block; width: 100%; height: auto; margin-top: 4px; }
svg polyline { fill: none; stroke: var(--line); stroke-width: 1.5; }
svg circle { fill: var(--line); } svg line { stroke: var(--grid); }
a { color: var(--line); }
"""


def group_of(name: str) -> str:
    """The section a benchmark is shown in: its benchmark ID or function, before the first slash."""
    return name.split("/", 1)[0].split("<", 1)[0]


def chart(points: list[tuple[dict, float]], repo_url: str) -> str:
    """An inline SVG of one benchmark's time over the runs, one titled point per run."""
    values = [v for _, v in points]
    low, high = min(values), max(values)
    span = high - low or high or 1.0
    step = (CHART_WIDTH - 2 * PAD) / max(len(points) - 1, 1)
    xy = [(PAD + i * step, CHART_HEIGHT - PAD - (v - low) / span * (CHART_HEIGHT - 2 * PAD))
          for i, v in enumerate(values)]
    dots = []
    for (x, y), (run, value) in zip(xy, points):
        label = html.escape(f"{short(run['commit'])} {run.get('date', '')}: {format_time(value)}"
                            f" ({run.get('cpu') or 'unknown CPU'})")
        dot = f'<circle cx="{x:.1f}" cy="{y:.1f}" r="2.5"><title>{label}</title></circle>'
        if repo_url and run.get("commit"):
            dot = f'<a href="{html.escape(repo_url)}/commit/{html.escape(run["commit"])}">{dot}</a>'
        dots.append(dot)
    line = " ".join(f"{x:.1f},{y:.1f}" for x, y in xy)
    return (f'<svg viewBox="0 0 {CHART_WIDTH} {CHART_HEIGHT}" role="img" '
            f'aria-label="{html.escape(format_time(low))} to {html.escape(format_time(high))}">'
            f'<line x1="{PAD}" y1="{CHART_HEIGHT - PAD}" x2="{CHART_WIDTH - PAD}" y2="{CHART_HEIGHT - PAD}"/>'
            f'<polyline points="{line}"/>{"".join(dots)}</svg>')


def card(name: str, points: list[tuple[dict, float]], repo_url: str, threshold: float) -> str:
    now = points[-1][1]
    trend = ""
    if len(points) > 1:
        relative = now / points[-2][1] - 1
        cls = "worse" if relative > threshold else "muted"
        trend = f' <span class="{cls}">{relative:+.1%}</span>'
    low, high = min(v for _, v in points), max(v for _, v in points)
    return (f'<div class="card"><div class="name">{html.escape(name)}</div>'
            f'<div class="now">{format_time(now)}{trend} '
            f'<span class="muted">range {format_time(low)}–{format_time(high)}, {len(points)} runs</span></div>'
            f'{chart(points, repo_url)}</div>')


def render(runs: list[dict], last: int, repo_url: str, threshold: float) -> str:
    runs = runs[-last:]
    names = sorted({name for run in runs for name in run["benchmarks"]})
    groups: dict[str, list[str]] = {}
    for name in names:
        groups.setdefault(group_of(name), []).append(name)
    latest = runs[-1] if runs else {"commit": "", "date": "", "cpu": ""}
    out = ["<!doctype html>", '<html lang="en"><head><meta charset="utf-8">',
           '<meta name="viewport" content="width=device-width, initial-scale=1">',
           "<title>Stereon benchmark trend</title>", f"<style>{PAGE_STYLE}</style></head><body>",
           "<h1>Stereon benchmark trend</h1>",
           f"<p>stereon_bench on every push to main: the median CPU time per evaluated item (per iteration "
           f"where nothing is counted) over the repetitions of each run, the last {len(runs)} runs. "
           f"Latest: {html.escape(short(latest['commit']))} on {html.escape(latest.get('date', ''))}, "
           f"{html.escape(latest.get('cpu') or 'unknown CPU')}. Hosted runners differ from run to run, so a "
           f"step in a chart is read against the CPU in its points' titles. The change shown is against the "
           f"run before; above {threshold:.0%} it is marked.</p>"]
    for group, members in groups.items():
        out += [f"<h2>{html.escape(group)}</h2>", '<div class="grid">']
        for name in members:
            points = [(run, run["benchmarks"][name]["ns"]) for run in runs if name in run["benchmarks"]]
            out.append(card(name, points, repo_url, threshold))
        out.append("</div>")
    out += ["</body></html>", ""]
    return "\n".join(out)


# --- Self-test -----------------------------------------------------------------------


def _entry(name: str, index: int, cpu_ns: float, items: int | None = None, run_type: str = "iteration") -> dict:
    entry = {"name": name, "run_type": run_type, "repetition_index": index, "iterations": 100,
             "cpu_time": cpu_ns, "time_unit": "ns"}
    if items:
        entry["items_per_second"] = items / (cpu_ns * 1e-9)
    return entry


def self_test() -> int:
    results = {"benchmarks": [
        _entry("Direct/Line", 0, 200.0, items=100), _entry("Direct/Line", 1, 220.0, items=100),
        _entry("Direct/Line", 2, 180.0, items=100), _entry("Direct/Line_median", 0, 200.0, run_type="aggregate"),
        _entry("B0.3/FacesOfEdges/Stereon/10000", 0, 1.5e6),
        {"name": "Broken", "run_type": "iteration", "error_occurred": True, "cpu_time": 1, "time_unit": "ns"},
    ]}
    machine = read_machine("commit=abcdef0123456789\ndate=2026-10-25\ncpu=AMD EPYC 7763\nnot a pair\n")
    first = record(results, machine)
    assert first["commit"] == "abcdef0123456789" and first["cpu"] == "AMD EPYC 7763"
    line = first["benchmarks"]["Direct/Line"]
    assert abs(line["ns"] - 2.0) < 1e-12, line  # per item: 200 ns over 100 items
    assert line["repetitions"] == 3 and abs(line["spread"] - 0.1) < 1e-12, line
    assert set(first["benchmarks"]) == {"Direct/Line", "B0.3/FacesOfEdges/Stereon/10000"}
    assert first["benchmarks"]["B0.3/FacesOfEdges/Stereon/10000"]["ns"] == 1.5e6

    second = json.loads(json.dumps(first))
    second["commit"], second["cpu"] = "1234567890abcdef", "Intel Xeon 8370C"
    second["benchmarks"]["Direct/Line"]["ns"] = 2.5  # 25 % slower
    second["benchmarks"]["B0.3/FacesOfEdges/Stereon/10000"]["ns"] = 1.2e6  # 20 % faster
    second["benchmarks"]["New/Thing"] = {"ns": 5.0, "spread": 0.0, "repetitions": 1}

    with tempfile.TemporaryDirectory() as directory:
        runs_dir = Path(directory)
        # Written out of order: the file names, not the writes, set the order.
        (runs_dir / "20261025T120000Z-1234567890ab.json").write_text(json.dumps(second), encoding="utf-8")
        (runs_dir / "20261024T120000Z-abcdef012345.json").write_text(json.dumps(first), encoding="utf-8")
        (runs_dir / "notes.json").write_text(json.dumps({"schema": 99}), encoding="utf-8")
        runs = load_runs(runs_dir)
    assert [run["commit"] for run in runs] == ["abcdef0123456789", "1234567890abcdef"]

    summary, regressions = compare_report(runs, THRESHOLD)
    assert [c.name for c in regressions] == ["Direct/Line"], regressions
    assert abs(regressions[0].relative - 0.25) < 1e-12
    assert "| `Direct/Line` | 2 ns | 2.5 ns | +25.0% | 10% |" in summary, summary
    assert "| `B0.3/FacesOfEdges/Stereon/10000` | 1.5 ms | 1.2 ms | -20.0% |" in summary, summary
    assert "2 benchmarks compared, 1 slower and 1 faster" in summary, summary
    assert "CPU changed" in summary
    assert compare_report(runs[:1], THRESHOLD)[1] == []

    page = render(runs, LAST_RUNS, "https://github.com/onurtuncer/Stereon", THRESHOLD)
    assert page.count('<div class="card">') == 3, "one chart per benchmark"
    assert "<h2>Direct</h2>" in page and "<h2>B0.3</h2>" in page and "<h2>New</h2>" in page
    assert '<span class="worse">+25.0%</span>' in page
    assert "/commit/1234567890abcdef" in page
    assert render(runs, 1, "", THRESHOLD).count("<circle") == 3, "only the last run, one point per chart"
    escaped = render([{"commit": "", "date": "", "cpu": "", "benchmarks": {"A<int>/x&y": {"ns": 1.0}}}],
                     LAST_RUNS, "", THRESHOLD)
    assert "A&lt;int&gt;/x&amp;y" in escaped and "<h2>A</h2>" in escaped

    assert format_time(2.0) == "2 ns" and format_time(999.4) == "999 ns" and format_time(1500.0) == "1.5 µs"
    assert format_time(3.2e9) == "3.2 s"
    print("bench_trend.py: self-test passed")
    return 0


# --- Main --------------------------------------------------------------------------


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--self-test", action="store_true", help="run the built-in checks and exit")
    commands = parser.add_subparsers(dest="command")

    rec = commands.add_parser("record", help="one run's Google Benchmark JSON to a run file")
    rec.add_argument("results", type=Path)
    rec.add_argument("--machine", type=Path, required=True, help="key=value lines: commit, date, cpu, ...")
    rec.add_argument("--out", type=Path, required=True)

    cmp_ = commands.add_parser("compare", help="the latest run against the one before")
    cmp_.add_argument("runs", type=Path)
    cmp_.add_argument("--threshold", type=float, default=THRESHOLD)
    cmp_.add_argument("--out", type=Path, help="write the Markdown summary here instead of standard output")
    cmp_.add_argument("--annotate", action="store_true", help="print a GitHub warning per regression")

    ren = commands.add_parser("render", help="every run to one HTML page of charts")
    ren.add_argument("runs", type=Path)
    ren.add_argument("--out", type=Path, required=True)
    ren.add_argument("--last", type=int, default=LAST_RUNS)
    ren.add_argument("--repo-url", default="")
    ren.add_argument("--threshold", type=float, default=THRESHOLD)

    args = parser.parse_args()
    if args.self_test:
        return self_test()
    if args.command == "record":
        results = json.loads(args.results.read_text(encoding="utf-8"))
        run = record(results, read_machine(args.machine.read_text(encoding="utf-8")))
        if not run["benchmarks"]:
            print(f"bench_trend.py: no benchmark results in {args.results}", file=sys.stderr)
            return 1
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(run, indent=1, sort_keys=True) + "\n", encoding="utf-8")
        return 0
    if args.command == "compare":
        summary, regressions = compare_report(load_runs(args.runs), args.threshold)
        if args.out:
            args.out.write_text(summary, encoding="utf-8")
        else:
            sys.stdout.buffer.write(summary.encode("utf-8"))
        if args.annotate:
            for c in regressions:
                print(f"::warning title=Benchmark slower by {c.relative:.0%}::{c.name}: "
                      f"{format_time(c.before)} -> {format_time(c.after)} (spread {c.spread:.0%})")
        return 0
    if args.command == "render":
        runs = load_runs(args.runs)
        if not runs:
            print(f"bench_trend.py: no run files in {args.runs}", file=sys.stderr)
            return 1
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(render(runs, args.last, args.repo_url, args.threshold), encoding="utf-8")
        return 0
    parser.error("a command (record, compare, render) or --self-test is required")
    return 2


if __name__ == "__main__":
    sys.exit(main())
