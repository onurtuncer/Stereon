# ADR-0053 owner review package (v0.2)

- **Date:** 2026-11-26
- **Revision reviewed:** `main` at [`09c0505`](https://github.com/onurtuncer/Stereon/commit/09c05055c2e2737dfe7810307337bb2b7ec8e22b) (#342)
- **Prepared by:** the v0.2 verification agent (Agent E). This file only collects evidence. It records no review: the owner's review and the `inventory_review` field of [the inventory](intersection-metamorphic.json) are the owner's to write.
- **Check:** `python tools/check_metamorphic.py --release`

Every `verified` cell and defect record in the inventory rests today on Agent E's documented review, as each relations document says. The owner's review of those cells, defects and exceptions, and the owner's `inventory_review`, are still to be recorded. This package lists, for each family, what a cell cites and where the cited test ran green at the revision above, so that the review can be done from one page.

## CI runs at 09c0505

| Run | Job | What it covers |
| --- | --- | --- |
| Linux, push | [`gcc16-debug`](https://github.com/onurtuncer/Stereon/actions/runs/45829925674/job/139943711808), [`gcc16-release`](https://github.com/onurtuncer/Stereon/actions/runs/45829925674/job/139943711584), [`gcc16-asan`](https://github.com/onurtuncer/Stereon/actions/runs/45829925674/job/139943711702) | every per-PR test cited below |
| Windows, push | [`gcc16-debug`](https://github.com/onurtuncer/Stereon/actions/runs/45829925840/job/139943736324), [`gcc16-windows`](https://github.com/onurtuncer/Stereon/actions/runs/45829925840/job/139943737037) | the same |
| macOS, push | [`gcc16-debug`](https://github.com/onurtuncer/Stereon/actions/runs/45829925624/job/139943711047) | the same |
| Robustness, dispatched | [exhaustive](https://github.com/onurtuncer/Stereon/actions/runs/45834398857/job/139959311443) | the hidden `[exhaustive]` planar tests the planar cells cite |
| Intersection metamorphic evidence, dispatched | [fault replay](https://github.com/onurtuncer/Stereon/actions/runs/45834385843/job/139959267243) | every fault of `tools/mutation/intersection-backfill*.json`: 30 killed of 30, none survived |
| Intersection release evidence, dispatched | [ADR-0053 release qualification](https://github.com/onurtuncer/Stereon/actions/runs/45834378374/job/139959242404) | fails: the open items under *What still blocks the release check* |

Each cited test was checked by name in these logs and is listed as `Passed`. None of the per-PR tests is `[!mayfail]`.

## Relation cells by family

States are those on `main` at 09c0505. A cell's domain, oracle and replay are in its family's relations document.

### Curve–curve: 36 ordered pairs, 144 cells ([relations](intersection-metamorphic-relations-curve-curve.md))

| Cells | State | Test (`tests/intersect/`) | Green at 09c0505 |
| --- | --- | --- | --- |
| MR-01, all 36 pairs | verified | `CurveCorpusRelationTests.cpp`: "MR-01: curve-curve corpus claims hold under exact rigid motions" | Linux, Windows, macOS |
| MR-02, all 36 pairs | verified | `CurveCorpusTests.cpp`: "every curve-curve corpus case meets exactly when its gap is within the resolution" | Linux, Windows, macOS |
| MR-03, all 36 pairs | verified | `CurveCorpusRelationTests.cpp`: "MR-03: curve-curve corpus claims hold under power-of-two scales" | Linux, Windows, macOS |
| MR-04, all 36 pairs | verified | `CurveCorpusRelationTests.cpp`: "MR-04: curve-curve corpus claims hold in equivalent charts" | Linux, Windows, macOS |

Review points: the claim is local, a constructed closest approach at gaps of 0, 0.4, 1.5 and 10 resolutions, with no global completeness. Rational 3D B-splines are outside this corpus.

### Planar: 36 ordered pairs, 144 cells ([relations](intersection-metamorphic-relations-planar.md))

| Cells | State | Test (`tests/intersect/PlanarCorpusTests.cpp`) | Green at 09c0505 |
| --- | --- | --- | --- |
| MR-01, MR-02, MR-03, all 36 pairs | MR implemented | "no planar-curves corpus case is answered wrongly, both ways" (`[exhaustive]`) | Robustness |
| MR-04, all 36 pairs | MR implemented | "MR-04: every planar-curves corpus case answers alike in an equivalent chart, both ways" (`[exhaustive]`) | Robustness |

These cells stay `MR implemented` under the audit correction of 2026-11-12 until completeness is strict. In the Robustness run, corpus completeness has 0 undocumented refusals. MR-04 completeness has 48, all in reversed identical arcs of two B-splines (`BSplineCurve2d` with itself). #351 makes both strict where the evidence allows, and will move 143 cells; `planar:BSplineCurve2d:BSplineCurve2d/MR-04` stays open.

### Curve–surface: 42 ordered pairs, 168 cells ([relations](intersection-metamorphic-relations-curve-surface.md))

| Cells | State | Test (`tests/intersect/CurveSurfaceRelationTests.cpp`) | Green at 09c0505 |
| --- | --- | --- | --- |
| MR-01, line, circle, ellipse and B-spline curve against all 7 surfaces (28 pairs) | verified | "MR-01: curve-surface crossings move with exact rigid motions" | Linux, Windows, macOS |
| MR-03, the same 28 pairs | verified | "MR-03: curve-surface crossings scale with the geometry and its resolution" | Linux, Windows, macOS |
| MR-04, the same 28 pairs | verified | "MR-04: curve-surface crossings keep their place in an equivalent curve chart" | Linux, Windows, macOS |
| MR-01, MR-03, MR-04, hyperbola and parabola against all 7 surfaces (42 cells) | MR implemented; `verified` in #346 | the same three tests | Linux, Windows, macOS |
| MR-02, all 42 pairs | exception reviewed | none: no exposed surface-curve operation (ADR-0053 section 1) | n/a |

Review points:
- The fixture is one transverse crossing at a known parameter per pair. Tangencies and overlaps are left to the handler tests.
- The MR-02 exception is a statement about the API, and needs the owner's explicit review (ADR-0053 section 3).

### Surface–surface: 49 ordered pairs, 196 cells ([relations](intersection-metamorphic-relations-surface-surface.md))

| Cells | State | Test (`tests/intersect/SurfaceRelationTests.cpp`) | Green at 09c0505 |
| --- | --- | --- | --- |
| MR-01, all 49 pairs | verified | "MR-01: surface-surface answers move with exact rigid motions" | Linux, Windows, macOS |
| MR-02, all 49 pairs | verified | "MR-02: surface-surface answers hold with the operands exchanged" | Linux, Windows, macOS |
| MR-03, all 49 pairs | verified | "MR-03: surface-surface answers scale with the geometry and its tolerances" | Linux, Windows, macOS |
| MR-04, all 49 pairs | verified | "MR-04: surface-surface answers keep their place in equivalent surface charts" | Linux, Windows, macOS |

Review points:
- There is one transverse closed-form configuration per unordered pair.
- No fault spot-check was run for these cells when they were qualified (#329).
- Plane:Cylinder and Cylinder:Plane MR-04 were closed by AUDIT-20261122-plane-cylinder-periodic-box (#342).

## Defect records

| Records | State | Evidence | Red at 09c0505 |
| --- | --- | --- | --- |
| PR-112, PR-112-budget, PR-120, PR-159, PR-194, PR-200, PR-205, PR-206, PR-207, PR-218, PR-230, PR-230-pole, PR-236, PR-236-tangent, PR-238, PR-240, PR-241, PR-244, PR-247, PR-263-translated-boundary, SSI-translated-boundary, AUDIT-20261112-cone-polish, AUDIT-20261113-narrow-sphere-seam, CI-272-ellipse-plane, AUDIT-20261122-plane-cylinder-periodic-box (25) | verified | [triage](intersection-metamorphic-triage.md), batches 1–3 | each fault killed twice in the fault replay above. PR-200's own fault is masked on current code; its two-site replay is manual (batch 1). PR-240 shares PR-241's fault. |
| PR-75, PR-90, PR-106, PR-126, PR-135, PR-267, PR-281, PR-295 (8) | not a defect | [triage](intersection-metamorphic-triage.md): no intersection result changed, or a documented limitation lifted | n/a |
| PR-239 (1) | not a defect | the owner's decision of 2026-11-12: a contract change | n/a |
| PR-271, PR-289, PR-304, PR-312, PR-315, PR-316, PR-317, PR-318 (8) | untriaged on `main`; `MR implemented` in #348, `verified` when its replay is green | triage batch 4 (#348) | #348's fault replay |
| PR-301, PR-322, PR-323 (3) | untriaged | triage, *Leads still open after batch 4* | none yet |

## What still blocks the release check

At 09c0505 `--release` reports 199 errors:

| Errors | What | Who closes it |
| --- | --- | --- |
| 144 | planar cells `MR implemented` | #351 (143 cells); `planar:BSplineCurve2d:BSplineCurve2d/MR-04` needs a kernel answer for reversed identical B-spline arcs, or the owner's decision to document the refusal |
| 42 | curve–surface hyperbola and parabola cells `MR implemented` | #346 |
| 11 | untriaged defect records | #348 (8); PR-301, PR-322 and PR-323 remain |
| 1 | `inventory_complete` is false | see below |
| 1 | `inventory_review` missing | the owner only |

`inventory_complete` is false because the historical audit read merged pull requests only. The audit after #263 did not read closed issues, failing corpus entries or regression-test comments since #263. The owner decides whether the merged-PR audit, with batches 1–4, is the complete historical inventory. If so, the owner records `inventory_complete: true` together with `inventory_review`. If not, the scope still to be read is those three sources, from #263 (`cbeddfb`) to the release revision.

## What the owner records

1. A review of the cells and defect records above, or of a sample with the reason for it. The `review` fields stay as they are, and the owner's review is referenced from `inventory_review`.
2. The MR-02 exception of the 42 curve–surface pairs.
3. `inventory_complete` and `inventory_review` in `intersection-metamorphic.json`.
4. The decision on `planar:BSplineCurve2d:BSplineCurve2d/MR-04` (see #351).
