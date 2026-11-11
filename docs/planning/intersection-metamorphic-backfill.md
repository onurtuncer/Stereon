# Intersection defect-to-metamorphic regression ledger

- **Status:** Initial backfill inventory (incomplete pending full PR, issue and corpus audit)
- **Date:** 2026-11-11
- **Policy:** [ADR-0053](../adr/0053-metamorphic-testing.md)
- **Scope:** Confirmed repairs and candidate defects identified from merged Stereon PR titles,
  recent regression tests and the intersection adversarial suite.

## Rules

Every **confirmed** intersection defect is a permanent obligation under
ADR-0053: preserve a failing direct reproducer, add at least one applicable
metamorphic regression with independent geometric/semantic checks, and show
that the original fault is detected. If no sound nontrivial relation applies,
record a reviewed mathematical exception and an independent oracle.

**This ledger is not a claim of historical completeness.** PR titles
identify leads, not necessarily unique or confirmed defects. The backfill
is not complete until all relevant merged repairs, regression comments,
issues and corpus failures are examined. Keep entries visible until their
verification evidence is linked; do not silently change status to complete.

## Known repair entries

| PR / defect | Failure mode | Candidate MR | Current evidence | Backfill state |
| --- | --- | --- | --- | --- |
| [#200](https://github.com/onurtuncer/Stereon/pull/200) | Quadric root placement from far line origin | MR-01, MR-04 | Direct tests exist in original PR; compare with #205 failure | untriaged |
| [#205](https://github.com/onurtuncer/Stereon/pull/205) | Far-origin line collapses two crossings | MR-01, MR-04 | Deterministic direct + transformed checks in `tests/intersect/AdversarialTests.cpp`; randomized origin-shift suite already present | MR implemented (verify faulty revision and CI) |
| [#206](https://github.com/onurtuncer/Stereon/pull/206) | Singleton-domain coplanar curve lost | MR-04, MR-07 | Direct regression in `tests/intersect/ConicSurfaceTests.cpp` or corresponding conic test; link exact location after inspection | fixture captured |
| [#207](https://github.com/onurtuncer/Stereon/pull/207) | Tiny traced branch/sliver discarded | MR-01, MR-03 | Direct fixture from repair PR needs audit | untriaged |
| [#218](https://github.com/onurtuncer/Stereon/pull/218) | Line–cone result depends on line origin | MR-04, MR-01 | Inspect repair regression and add transform-specific test | untriaged |
| [#230](https://github.com/onurtuncer/Stereon/pull/230) | Plane–ellipsoid tangent and longitude-box errors | MR-02, MR-04 | Original PR's tangent/pole tests already swap operands; test stability and link replay evidence | MR implemented (verify faulty revision and CI) |
| [#236](https://github.com/onurtuncer/Stereon/pull/236) | Two corpus-declared SSI configurations incorrect | MR-01, MR-02, MR-03 | Inspect corpus case identities and direct fixtures | untriaged |
| [#238](https://github.com/onurtuncer/Stereon/pull/238) | Property counterexamples from CI seeds | Relation depends on exact counterexample | Recover seeds and fixes from PR | untriaged |
| [#239](https://github.com/onurtuncer/Stereon/pull/239) | 3D line resolution semantics differed from spline | MR-02, MR-03 | Inspect direct test and relation applicability | untriaged |
| [#240](https://github.com/onurtuncer/Stereon/pull/240), [#241](https://github.com/onurtuncer/Stereon/pull/241) | B-spline curve meeting on domain boundary | MR-04, MR-02 | Re-landed repair; check one root-cause entry vs two attempts | untriaged |
| [#244](https://github.com/onurtuncer/Stereon/pull/244) | Near-axis sphere–cone contact misclassified as points | MR-02, MR-01 | Direct regression in `tests/intersect/ElementaryPairsTests.cpp` | fixture captured |
| [#247](https://github.com/onurtuncer/Stereon/pull/247) | Unequal near-parallel cones falsely coincident | MR-02 | Existing direct regression extended to swap operands in `tests/intersect/ElementaryPairsTests.cpp` | MR implemented (verify faulty revision and CI) |

## Candidate repairs needing defect classification

The following merged PRs also touch intersections, robustness or regressions,
but their defects and required MRs must be established from the diffs before
counting or closing obligations:

- [#75](https://github.com/onurtuncer/Stereon/pull/75),
  [#90](https://github.com/onurtuncer/Stereon/pull/90),
  [#106](https://github.com/onurtuncer/Stereon/pull/106),
  [#112](https://github.com/onurtuncer/Stereon/pull/112),
  [#120](https://github.com/onurtuncer/Stereon/pull/120),
  [#126](https://github.com/onurtuncer/Stereon/pull/126),
  [#135](https://github.com/onurtuncer/Stereon/pull/135),
  [#159](https://github.com/onurtuncer/Stereon/pull/159),
  [#194](https://github.com/onurtuncer/Stereon/pull/194).

## Exit criteria

- All confirmed intersection defects in the audited history are identified.
- Each entry links its minimal direct reproducer, the metamorphic test and
  MR, seeded replay, independent semantic checks, and CI evidence.
- Every verified entry has demonstrated red-on-defect/green-on-repair
  behavior, or a documented reviewed exception with independent oracle.
- No unidentified or untriaged repair is hidden behind an aggregate pass
  percentage.
- The ledger remains updated for each future intersection repair PR.

## Review note

Tests added alongside this ledger have **not** been built or run by this
documentation update. Their CI outcome and red-on-defect proof are still
needed before changing them to `verified`.

## Executable inventory and release gate

[`intersection-metamorphic.json`](intersection-metamorphic.json) is the
machine-readable obligation and handler/MR inventory. The table above explains
the historical leads; the JSON owns their checked states. A candidate that
proves not to be a defect stays in the JSON as `not a defect`, with `rationale`
and `review` evidence. Do not delete it. PR #240 and #241 remain separate leads
until reviewed; shared tests may discharge both without losing either ID.

Development CI and CTest run `python tools/check_metamorphic.py`; PR CI also
uses `--base` to refuse deleted obligations. This validates records, not release
readiness. All 163 ordered built-in operand pairs are inventoried automatically
against the curve and surface catalogues, with MR-01/02/03/04 cells. Initial
`untriaged` cells are deliberately not claims of coverage or applicability.

Before qualifying v0.2, run `python tools/check_metamorphic.py --release`.
The **Intersection release evidence** workflow runs this on v0.2 tags and on
manual dispatch. Its failure blocks qualification; it does not roll back a tag
or replace the other release gates. Do not publish a qualified release without
this check passing. Existing open obligations intentionally make it fail.

Each `verified` defect must supply `fixture`, `direct`, `metamorphic`, `oracle`,
`replay`, `red`, `green`, `review`, and applicable `relations`. An
`exception reviewed` defect instead supplies `rationale` (considered relations
and why none apply); its direct regression, independent oracle, replay and
red/green evidence remain mandatory. Handler cells require `domain` and
`review`; verified cells also require `test`, `oracle`, `replay`, and `green`;
reviewed exceptions require `rationale`, including unsupported domains.

Evidence is an existing repository path, optionally followed by `#` and an
exact test name/text anchor, or an HTTPS permalink to a review/CI artifact.
Use immutable revisions/run IDs for external evidence. The checker validates
references and mandatory evidence presence; a reviewer must still inspect the
scientific claims and CI outcomes. `inventory_complete` and `inventory_review`
may be set only after the historical audit, not merely after tests pass.
