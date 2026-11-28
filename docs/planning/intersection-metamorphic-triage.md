# Intersection metamorphic triage, batch 1

- **Date:** 2026-11-12
- **Base revision:** `main` at `227d79ab24c423b67aba155e5eb7cc5125f681f6`
- **Policy:** [ADR-0053](../adr/0053-metamorphic-testing.md), [backfill ledger](intersection-metamorphic-backfill.md)
- **Toolchain:** GCC 16 (MSYS2 UCRT64), preset `gcc16-debug`, Windows 11
- **Reviewer:** the v0.2 verification agent (Agent E), independent of the repairs' authors. This is the triage review the JSON's `review` fields cite. It is not the owner's `inventory_review`, which this batch does not write.

Each lead below was triaged from its merged diff, its PR description and the regression tests it added. Leads that change no intersection code are `not a defect` for this ledger: they stay recorded with the rationale. Confirmed intersection defects get the following:

- a direct regression (the repair PR's own test);
- a metamorphic regression in `tests/intersect/MetamorphicRegressionTests.cpp` with an independent oracle;
- a fault reintroduction in `tools/mutation/intersection-backfill.json`;
- the red/green record below.

## Leads that are not intersection defects

### PR-75

The audit fixes of 2026-10-20 (`0bbd9d51e4`) touch `src/core`, `src/curve` (B-spline conversion of trimmed curves), `src/io`, `src/surface` (sweeps, trimmed surfaces) and `src/topo`. They touch nothing under `src/intersect/`, and no intersection result changes. Their regressions live in the curve, io, surface, topo and probe suites. Not an intersection defect. The curve conversion and sweep findings belong to those libraries' own regression suites (ADR-0032 mutation catalogue, finding F9).

### PR-90

Test expectations only (`e9504da675`): `tests/topo/GeometricValidityTests.cpp` and `tests/topo/OperationTests.cpp`, with no kernel change. Not a defect of any kind in the intersection handlers.

### PR-106

Transaction identities, non-finite validation samples and display sampling (`05eb8f0699`): `src/topo`, `src/mesh`. No intersection code. Not an intersection defect.

### PR-126

The Boolean classification findings (`5fdde67bc5`) are in `src/boolean` (winding sampling, seeds, certified rays). The commit's intersection hunks (`src/intersect/Exact.hpp`, and the call sites in `CylinderCylinder.cpp`, `CylinderCylinderSkew.cpp`, `PlaneCylinder.cpp`, `PlanePlane.cpp`) defer the exact `Expansion` behind a `std::function` and call `Value()` where the value was read before. That is a behaviour-preserving evaluation-order change: the filtered sign reads the same interval first and builds the same exact expansion when it straddles zero. No intersection answer was wrong before it. Not an intersection defect. The ray decisions it repaired are Boolean classification (`Detail::RayMeets`), outside the `Intersect` handlers.

### PR-135

Boolean cavities, `Transformed` on trimmed curves and surfaces, and concurrent archive saves (`1236aaa786`): `src/boolean`, `src/curve/Trimmed.cpp`, `src/surface/Trimmed.cpp`, `src/io`. No intersection handler changes. The cavity placement uses classifier ray crossings (`Detail::FaceCrossingsOf`), which is Boolean classification. Not an intersection defect.

## Confirmed intersection defects

For PR-207 and PR-238 the replay is the campaign's: `python tools/mutate.py --catalogue tools/mutation/intersection-backfill.json --only <ID>` on committed HEAD. It builds a disposable worktree, requires the named test to pass on the repaired code, then restores the fault and requires the named test to fail twice. Every fixture is deterministic; there is no random seed. PR-200's replay is manual, at two sites; see the red and green record.

### PR-200

- **Defect:** `LineSurface.cpp` placed a line's crossings of a sphere or cylinder from coefficients formed at the line's own origin. Those carry eps |w|² for an origin |w| away, so from 24 units the crossings lay 1.4e-13 off the cylinder, about 600 units in the last place. Repair `e8db99b53bf4396bf1f48d92b195f7b06452b85c` polishes each root from the line re-centred at its rounded point.
- **Direct regression:** `tests/intersect/LineSurfaceTests.cpp`, "a crossing from a far line origin lies on the quadric to rounding".
- **Metamorphic regression:** "MR PR-200 far line origin crossings on the quadric".
  - MR-04: the origin moved back by m ∈ {0, 8, 64} directions, exact on these dyadic inputs.
  - MR-01: the three cyclic axis permutations, exact.
  - Both the sphere and the cylinder.
- **Oracle:** the closed-form distance of each crossing from the sphere's centre, or from the cylinder's axis, minus the radius. It must be within 32 eps |p| plus 4 eps |t| |d|, the line's own parameter grid. The parameter must move by exactly m (to 1e-12 m). The kernel's projection is not used.
- **Fault:** the polish removed (`t = rounded`) together with PR-205's rounded coefficients. The single-site fault is masked on current `main`; see the red and green record.

### PR-207

- **Defect:** `ImplicitTrace.cpp` dropped a traced chain too short to fit a guide, a sliver of a branch at a corner of the search box. The result was empty, claiming a separation without proof (ADR-0029). Repair `f2138cf159a179e31d8adf56397a04e82e5814f3` reports it as a point with a certified bound.
- **Direct regression:** `tests/intersect/ImplicitTraceTests.cpp`, "a sliver of a branch at a corner of the box is a point, not a separation".
- **Metamorphic regression:** "MR PR-207 traced sliver at a box corner is a point".
  - MR-01: the three cyclic axis permutations, exact. A translation or scale is not used, because it could round away the 1e-13 premise.
  - MR-02: the operands exchanged.
- **Oracle:** the closed-form sliver point (2 cos u0, 2 sin u0, sqrt(4 cos u0 - 1)) at u0 = pi/6, and the closed-form distances from the sphere's centre and the cylinder's axis, each within the point's bound + 1e-12. The result must not be empty.
- **Fault:** the sliver point is computed but not pushed into the result.

### PR-238

- **Defect:** `PlaneEllipsoid.cpp` let a sample with a non-finite longitude rate reach the interpolation. A plane made through a pole and rounded off it by a unit in the last place escapes the exact pole test, so the error escaped as `InvalidInput` instead of failing closed with `IntersectionFailed`. Repair `721c41b86ccb1cbf1103e265e568dd8eb77e889b` refuses such a sample. The companion PeriodicBSpline change in that PR was a property defect, not a kernel one, per its description.
- **Direct regression:** `tests/intersect/PlaneEllipsoidTests.cpp`, "a plane made through a pole and rounded off it fails closed, as on the sphere".
- **Metamorphic regression:** "MR PR-238 rounded pole plane fails closed on ellipsoid and sphere".
  - MR-02: the operands exchanged.
  - MR-04: the round ellipsoid replaced by the sphere, the same point set in an equivalent chart.
  - No rigid motion or scale: either would re-round the plane and change the premise.
- **Oracle:** the premise, checked independently (the pole's signed distance from the rounded plane is at most 1e-15), and the fail-closed contract. Every order and representation must fail with `IntersectionFailed`. This is a fail-closed contract regression, not resolved geometry.
- **Fault:** the finite-rate check is disabled.

## Leads left open

After batch 2, PR-194 and PR-239 remained; both are now resolved below (owner decisions of 2026-11-12). None is closed. Closing needs a metamorphic test and a fault that is actually killed, and this batch has not produced either for these leads.

## Red and green record

### Fault replay of PR-207 and PR-238 (tools/mutate.py)

- **Command:** `PYTHONUTF8=1 python tools/mutate.py --catalogue tools/mutation/intersection-backfill.json --only PR-200,PR-207,PR-238 --work build/metamorphic-triage`, run on committed HEAD `71d7b5e671eb8fc6950d4ce029933cea7253c36a` (this branch's first commit, on `main` 227d79a).
- **Toolchain:** g++ (Rev4, Built by MSYS2 project) 16.2.0; cmake 4.2.3-msvc3; Windows-11-10.0.26200-SP0; preset `gcc16-debug`; wall time 199.8 s.
- **Baseline (green on repair):** the disposable worktree built `stereon_intersect_tests` unmodified, and every selected test passed before any fault was applied.
- **PR-207:** killed. "intersect: MR PR-207 traced sliver at a box corner is a point" failed in both runs; the declared diagnostic was found.
- **PR-238:** killed. "intersect: MR PR-238 rounded pole plane fails closed on ellipsoid and sphere" failed in both runs; the declared diagnostic was found.
- **PR-200:** survived; the tests reach the site. See below.

### Two-site replay of PR-200 (manual)

Removing the polish alone (the single-site fault) is not observable on current `main`. Neither the MR test nor the PR's own direct test ("a crossing from a far line origin lies on the quadric to rounding") fails. The later repair of PR-205 (`3a147886`) forms the quadratic with `QuadraticOf` instead of rounded doubles, which already places these crossings to rounding. The defect PR-200 repaired existed with both rounded coefficients and no polish, so it was reintroduced at both sites: the catalogue's PR-205 `find`/`replace` and the PR-200 one (`t = rounded`), applied together to `src/intersect/LineSurface.cpp` in the working tree, rebuilt with `cmake --build --preset gcc16-debug --target stereon_intersect_tests`.

- **Red:** "MR PR-200 far line origin crossings on the quadric" failed on each of two runs, with the same assertions. The closed-form distance from the cylinder was 1.4155e-13 against the bound 3.441e-14 and 1.4178e-13 against 4.233e-14, the defect's 1.41e-13 of PR #200. The direct test failed as well.
- **Controls:** with the PR-205 fault alone, MR PR-200 passes (168 assertions) while MR PR-205 fails (162 of 900); with the PR-200 fault alone, both pass. So the test detects PR-200's fault in the code it was repaired in, and on current `main` the polish and `QuadraticOf` are redundant protections. That redundancy is why the single-site PR-200 fault is not in the CI catalogue: `mutate.py` takes one site per fault, and that fault would survive by design. The workflow cannot replay this evidence; this record is it.
- **Green:** the source restored (`git status` clean under `src/`), rebuilt, and MR PR-200 passed (168 assertions).

### Green on main

At 227d79a with this branch's tests, built with `gcc16-debug`, `stereon_intersect_tests "MR PR-*"` passed: 11 test cases, 6861 assertions, the eight earlier campaign tests included.

## Batch 2

- **Date:** 2026-11-12
- **Branch:** `verify/metamorphic-triage-2`, stacked on batch 1 (b78bb2d, on `main` 227d79a)
- **Reviewer, toolchain and meaning of "verified":** as for batch 1. "Verified" means this agent's documented review, with the owner's review still pending.

This batch closes the earlier campaign's records (PR-205, 206, 218, 230, 230-pole, 240, 241, 244, 247, PR-263-translated-boundary) and the fixtures captured in batch 1 (PR-112, 120, 159, 236), one ID per root cause. PR-194 (sub-IDs) and PR-239 (defect or contract change) stay open for the owner.

### Campaign replay of the eight earlier faults

- **Command:** `PYTHONUTF8=1 python tools/mutate.py --catalogue tools/mutation/intersection-backfill.json --only PR-205,PR-206,PR-218,PR-244,PR-247,PR-230,PR-230-pole,PR-241 --work build/metamorphic-campaign`, run on committed HEAD `b78bb2d76ac6d3af90baf7abaf83990cf0c3df1b`.
- **Toolchain:** g++ (Rev4, Built by MSYS2 project) 16.2.0; cmake 4.2.3-msvc3; Windows-11-10.0.26200-SP0; preset `gcc16-debug`; wall time 2131.2 s.
- **Baseline:** built unmodified, every selected test passed.
- **Result:** score 8 / 8. Each fault was killed by its own named test, failing in both runs, with the declared diagnostic found: PR-205, PR-206, PR-218, PR-244, PR-247, PR-230, PR-230-pole, PR-241. PR-240 and PR-241 are one root cause (the shared boundary search, as the backfill ledger records), so the PR-241 kill discharges both IDs. PR-244 and PR-247 are fail-closed contract regressions: their red shows the contract breaks without the repair, not resolved geometry.
- **Green on main:** see "Green on main" in the batch 1 record. The same MR tests pass, together with the batch 2 tests below (7 test cases, 3852 assertions).

### PR-112

- **Defect:** coincident cylinders returned the whole requested boxes as their overlap, whatever the other's box, even boxes apart (`CylinderCylinder.cpp`, repair `93fd35313fd0047a2fc1d6536bd32137a27f120f`: `ClipOverlap`).
- **Direct regression:** "coincident cylinders intersect their requested height domains" and "coincident cylinder arcs split at the periodic seam".
- **Metamorphic regression:** "MR PR-112 coincident cylinders clipped to the requested domains". Every placement of the campaign: MR-01, and MR-03 with heights and the resolution scaled but not angles. Both orders (MR-02).
- **Oracle:** the closed-form intersection of the requested intervals: heights [0, s] against [2s, 3s] are apart, against [0.5s, 2s] they give [0.5s, s]; arcs [0, 1] and [2, 3] are apart.
- **Fault:** the overlap is returned with the unclipped boxes.

### PR-112-budget

- **Defect:** the plane-cylinder ellipse pcurve converted and allocated a sample count that a tiny resolution made unbounded (`PlaneCylinder.cpp`, same repair: `SampleBudget`). This is a new ID, split from PR-112 as a separate root cause.
- **Direct regression:** "ellipse pcurve work is bounded before converting or allocating its samples".
- **Metamorphic regression:** "MR PR-112-budget plane cylinder pcurve work is bounded". Every placement with the resolution 1e-20 times the scale: the spacing scales as the radius, so the sample count is invariant (MR-01, MR-03). Both orders (MR-02).
- **Oracle:** the count is about 2e5, beyond the budget 2^16, so the result must be `CertificationFailed` naming the sample budget. This is a fail-closed contract regression.
- **Fault:** the budget is raised by 2^40. The overflow guard is kept, so the fault reaches no undefined conversion.

### PR-120

- **Defect:** skew cylinders did not clip B's angular sector. This is PR #120's description; the repair is in its parent `3a2ba895be415efe0198af878a2630738fd1a097` (#118, `SectorFit` in `CylinderCylinderSkew.cpp`), since #120's own diff is tests only.
- **Direct regression:** "skew cylinders exclude a wide angular sector containing no intersection", "skew cylinder singleton sectors reject the opposite ray", "skew cylinder singleton sectors retain both generator intersections" and "skew cylinders clip wide sectors at either angular boundary".
- **Metamorphic regression:** "MR PR-120 skew cylinder angular sectors of B". Every placement (MR-01, MR-03), both orders (MR-02).
- **Oracle:** the closed form. A (radius 3 about z) and B (radius 2 along x through (0, 4, 0)) meet only where cos u_B <= -1/2. So the sector [-2, 2], also shifted by a turn, and the singleton {0} must be empty; the sector [2.5, 3.8], inside [2 pi / 3, 4 pi / 3], must not be.
- **Fault:** B's sector fit forced to Inside.

### PR-159

- **Defect:** skew cylinders internally tangent on both seams were refused. A sliver piece under the parameter merge gap 1e-9 of the span could not be fitted (`CylinderCylinderSkew.cpp`, repair `35f3d2370958b445985c2bb55c24cb9ab5d074d9`).
- **Direct regression:** "internally tangent cylinders touching on both seams meet at a node".
- **Metamorphic regression:** "MR PR-159 skew cylinders tangent on both seams meet at a node". Untranslated placements (MR-01 by exact permutation, MR-03 by powers of two) and both orders (MR-02). A translation could round the tangency's premise.
- **Oracle:** one node at the closed-form point of tangency, four branch ends there, and every branch sample on both cylinders by the closed-form distance from each axis.
- **Fault:** `CutsOf` merges cuts by the old parameter gap, 1e-9 of the span, instead of the piece's extent in space. The repair's other hunk, `BisectSeam` returning the seam end, is not observable alone: with it reverted and the merge kept, the direct test passes (8078 assertions, twice). The merge is the root cause.

### PR-236

- **Defect:** concentric spheres within the pair tolerance were an Exact overlap at distance 0 (`ElementaryPairs.cpp`, repair `252c1350cfd13429e7f406c70850dd8e288ee205`).
- **Direct regression:** "coaxial surfaces of revolution within the tolerance are declared, not exact".
- **Metamorphic regression:** "MR PR-236 concentric spheres within the tolerance are declared". Every placement with the tolerance scaled (MR-01, MR-03), both orders (MR-02).
- **Oracle:** Tolerant evidence with distance |R_b - R_a| (the closed form), and apart beyond the tolerance.
- **Fault:** the gap is ignored when declaring, which gives Exact evidence.

### PR-236-tangent

- **Defect:** a sphere on a cylinder's axis, with radii within the pair tolerance, failed the whole intersection, because the declared tangent circle was certified against the resolution (same repair). This is a new ID, split from PR-236 as a separate root cause.
- **Direct regression:** "a tangent circle declared within the tolerance is certified within it".
- **Metamorphic regression:** "MR PR-236-tangent sphere on a cylinder's axis within the tolerance". Every placement with lengths and the tolerance scaled (MR-01, MR-03), both orders (MR-02).
- **Oracle:** one tangential, certified curve, with its bound between 0.8 of the closed-form gap and the tolerance.
- **Fault:** the declared tangency is certified against the resolution only.

### PR-263-translated-boundary

- **Defect:** the PR-241 fixture failed natively at scale 1/4 translated by (8, -4, 2) (Actions run 43052745077, job 130633998895). That run is the original red on the defective revision. The pair search inverted an ill-conditioned Krawczyk matrix where no stationary point could exist (repair `cbeddfb8c0dffa25e5086bab19cc46fc6e7f5f53`, `NoStationaryPoint`).
- **Metamorphic regression:** "MR PR-263-translated-boundary spline end meeting away from the origin": the translated placements at every scale (MR-01, MR-03), reversal and knot maps (MR-04), both orders (MR-02).
- **Oracle:** the constrained closest pair of the two segments in closed form, as in PR-241's record.
- **Fault:** the stationary-point exclusion is disabled.

### Batch 2 fault replay

- **Command:** `PYTHONUTF8=1 python tools/mutate.py --catalogue tools/mutation/intersection-backfill.json --only PR-112,PR-112-budget,PR-120,PR-159,PR-236,PR-236-tangent,PR-263-translated-boundary --work build/metamorphic-batch2`, run on committed HEAD `ec0acc3f4b2313d2b8b788e5888c51ff648b583c`.
- **Toolchain:** g++ (Rev4, Built by MSYS2 project) 16.2.0; cmake 4.2.3-msvc3; Windows-11-10.0.26200-SP0; preset `gcc16-debug`; wall time 277.2 s.
- **Baseline:** built unmodified, every selected test passed.
- **Result:** score 7 / 7, every critical fault detected with its declared diagnostic. Each named MR test failed in both runs with its fault in: PR-112, PR-112-budget, PR-120, PR-159 (the merge-gap fault), PR-236, PR-236-tangent, PR-263-translated-boundary.
- **Pre-check:** before the catalogue entries were written, each fault was applied in the working tree and the original direct tests were run twice. All failed, except the PR-159 `BisectSeam` candidate, which passed and was therefore not used.

### Green on main, batch 2

On this branch (main 227d79a plus test-only changes), built with `gcc16-debug`, `stereon_intersect_tests "MR PR-*"`: all tests passed, 18 test cases, 10713 assertions. That covers every MR regression of the ledger, batch 1, batch 2 and the earlier campaign. The two `mutate.py` baselines (b78bb2d, ec0acc3) are green runs of the same tests in a clean checkout.

## PR-194 and PR-239

- **Date:** 2026-11-12
- **Branch:** `verify/pr-194`, stacked on batch 2 (#272)
- **Owner decisions (2026-11-12):** PR-194 stays one ledger entry and is not split. PR-239 is a contract change, not a defect.

### PR-194

One obligation covers six root causes in `b74d6e23ff00b2994222bdc611ab09161e256282`. Each has its direct regression from the PR, a metamorphic regression named "MR PR-194 ..." in `tests/intersect/MetamorphicRegressionTests.cpp`, and its own fault reintroduction:

| Root cause | Direct regression | MR test, relations | Fault |
| --- | --- | --- | --- |
| PlanePlane: a rotated, bounded coincidence encoded as the AABB of its image | "rotated bounded coincidence cannot be encoded as an AABB overlap" | "MR PR-194 rotated bounded plane coincidence", MR-01, MR-02 | the refusal of a non-box-preserving bounded overlap disabled |
| PlanePlane: the squared norm of a nonzero cross product underflowed | "near-parallel planes with an underflowed squared cross norm" | "MR PR-194 near-parallel planes underflowed cross norm", MR-01, MR-02 | `Normalized(cross)` restored for `cross / hypot` |
| PlaneCylinder: an uncut conic kept in a bounded plane box | "bounded plane patches cannot receive uncut conic branches" | "MR PR-194 bounded plane box gets no uncut conic branch", MR-01, MR-02 | the partial-ellipse guard disabled |
| SplineSurface: an undecided root reported as a contact | "unresolved spline tangency and near miss never become successful contacts", "an undecided spline root at a knot fails closed" | "MR PR-194 undecided spline root fails closed", MR-01, MR-04 | the `Decided` check disabled |
| Linearization: a result for other geometry differentiated | "JVP and VJP reject intersection results from stale geometry" | "MR PR-194 stale geometry is refused by JVP and VJP", MR-01, MR-03 | the consistent-contact check disabled |
| GenericSSI: the hull test bypassed the factory's knot and weight bounds | "direct NURBS construction cannot bypass factory weight bounds" | "MR PR-194 NURBS hull applies the factory bounds", MR-01, MR-02 | the factory checks skipped |

The oracles are closed forms, stated in each test's comment:
- the diamond |x - c| + |y - c| against the unit square;
- the x axis as the line of z = 0 and the tilted plane;
- the ellipse of the plane with normal (0, 1, 1) on the cylinder of radius 5;
- the double root of (1 - 2s)^2 and a root exactly at a knot;
- a stored point 0.01 off the moved plane;
- the factory's own verdict on the same data.

The PlaneCylinder, SplineSurface and Linearization rows are fail-closed contracts; the others are resolved geometry.

The campaign checker requires every mutant ID in `tools/mutation/intersection-backfill.json` to be a unique ledger record. So that catalogue holds one fault as `PR-194` (the PlanePlane AABB fault), and `tools/mutation/intersection-backfill-pr194.json` holds all six. The CI evidence workflow replays the first catalogue only; the six are replayed with `mutate.py --catalogue` on the second.

#### PR-194 fault replay

- **Commands:** `PYTHONUTF8=1 python tools/mutate.py --catalogue tools/mutation/intersection-backfill-pr194.json --work build/metamorphic-pr194`, then `... --catalogue tools/mutation/intersection-backfill.json --only PR-194`. Both ran on committed HEAD `97763d769e8bacd644ad63895c87d0c63123ab72`.
- **Toolchain:** g++ 16.2.0 (MSYS2), cmake 4.2.3-msvc3, Windows-11-10.0.26200, preset `gcc16-debug`; wall times 274.8 s and 172.2 s.
- **Baseline:** green in both runs before any fault was applied.
- **Red, all six root causes:** score 6 / 6, each fault killed in both of its runs by its own named MR test with the declared diagnostic found: plane-plane-aabb, plane-plane-underflow, plane-cylinder-uncut, spline-undecided, stale-linearization, nurbs-hull-bounds. The campaign entry `PR-194` (the plane-plane-aabb fault) is killed the same way, 1 / 1.
- **Green:** on this branch, `stereon_intersect_tests "MR PR-194*"` passes: 6 test cases, 282 assertions.

### PR-239

**Not a defect: a contract change** (owner decision, 2026-11-12). Before #239, `LineLine` decided 3D lines meeting exactly. That was the documented rule then, not a fault. The owner's decision of 2026-11-08 made a `Line` and the same segment as a degree-1 B-spline meet alike: within the Context's resolution. #239 (`0338cb3e44`) implemented that rule. There is no faulty behaviour to reintroduce, so red-on-fault does not apply: restoring the old code would test the old contract, not detect a defect.

The evidence that the new contract holds:
- **#239's direct tests** in `tests/intersect/LineLineTests.cpp` (`[line3d]`): "3D skew lines meet within the resolution, decided exactly" (the exact boundary, a gap of the resolution meets and the next double does not), "3D parallel lines within the resolution overlap, beyond it are apart", "3D segments a rounding apart meet at their ends", and the property "a 3D line meets as the same segment as a degree-1 B-spline does". All `[line3d]` tests pass on this branch: 15 test cases, 150 assertions.
- **MR-02 on Line:Line:** the curve-curve corpus runner "every curve-curve corpus case meets exactly when its gap is within the resolution" runs the line-line family in both orders.
- **MR-01, MR-03 and MR-04 on Line:Line:** the curve-curve relation tests of #275 (`tests/intersect/CurveCorpusRelationTests.cpp`, https://github.com/onurtuncer/Stereon/pull/275). They run the line-line cases under exact rotations and translation, power-of-two scales with the resolution scaled, and a line chart change.
- **#239's documented residual is closed.** The residual: nearly parallel segments touching end to end, whose supports' closest approach lies past both ends, met as lines but were missed as splines. PR-241 (`51a2795`) closed it with the spline boundary-minimum search. Its tests "nearly parallel segments touching end to end meet as lines and as splines" (`LineLineTests.cpp`) and "nearly parallel splines touching end to end meet at the ends" (`SplineCurvePairTests.cpp`) pass on this branch. With PR-241's campaign fault applied (`SearchEnds` disabled for non-empty searches), both fail in each of two runs. So that repair is what closes the residual, and no open defect remains from it.

## Batch 3: the open records of #305

- **Date:** 2026-11-22
- **Branch:** `verify/metamorphic-defects`, on `main` 898ed69
- **Reviewer and meaning of "verified":** as for batches 1 and 2: this agent's documented review, with the owner's review still pending.
- **Toolchain:** the shared disk was full during this work, so no local build was made. The red and green evidence is the pull request's native CI: the *Intersection metamorphic evidence* workflow (`tools/mutate.py` over every `intersection-backfill*.json`, GCC 16.2 on Linux, `gcc16-debug`), and the Linux, Windows and macOS CTest runs.

At #305 the ledger held four open defect records: `SSI-translated-boundary` (`untriaged`, although #279 repaired it), `AUDIT-20261112-cone-polish` and `AUDIT-20261113-narrow-sphere-seam` (`MR implemented` by #286, with no replayable red), and `CI-272-ellipse-plane` (`fixture captured`). Each now has a metamorphic regression named "MR <ID> ..." in `tests/intersect/MetamorphicRegressionTests.cpp` and a fault reintroduction that CI replays. Every fixture is deterministic; there is no random seed.

#### SSI-translated-boundary

- **Defect:** far from the origin the generic SSI refused (`IntersectionFailed`) a pair of B-spline surfaces that it answers at the origin, near a branch's end on the domain's boundary. There were two root causes, both repaired in #279 (`0965d8aa4335210199404d6d01f2412bac218e5b`):
  - *search-shift:* the patches' interval poles widen with the coordinates, while Krawczyk's acceptance of a root is measured in parameters. The search now runs on the surfaces moved exactly about the origin (`ShiftOf`).
  - *edge-halving:* a root just outside the other patch's cell left the operator undecided on a whole edge. An undecided edge is now halved, up to `EdgeHalvings` = 8 times.
- **Original red:** #279's description: on the base revision, refused at the shift (64, 0, 0) and beyond.
- **Direct regression:** "SSI-translated-boundary: far translations answer as at the origin" (`GenericSSITranslationTests.cpp`).
- **Metamorphic regression:** "MR SSI-translated-boundary far hyperbola branches end on the domain's sides". The same closed-form pair: 2x'^2 - 3y'^2 against the bowl x'^2 + y'^2 + 1/100, meeting where x'^2 - 4y'^2 = 1/100.
  - MR-01: the shifts (64, 0, 0) and (480000, -300000, 210000), with the identity and the cyclic turns of the axes.
  - MR-02: both orders.
  - MR-03: scale 4 with the tolerance and context.
  - MR-04: the first surface's u reversed.
- **Oracle:** the closed form above, checked on every branch sample in the surfaces' local coordinates. Also checked independently: each branch has both ends on one side of the unit square, x = 0 or x = 1, at |y'| = sqrt(0.06). That makes the boundary roots part of the claim.
- **Faults:** `tools/mutation/intersection-backfill-ssi-translated-boundary.json`, one per root cause: the shift forced to zero, and `EdgeHalvings` set to 0.

#### AUDIT-20261112-cone-polish

- **Defect:** a line's crossing of a cone was polished at origin + t direction, rounded twice. When a far origin cancels t direction, that loses the small offset, and near a generator the root's correction magnifies the loss. Found by CI seed 17298607298449989786 (Actions run 43251352135). Repaired in #286 (`76eb218f62b8828fc4dfe4ce50ed335491d91ac5`, `std::fma` in `PolishedOnCone`).
- **Direct regression:** "MR audit cone crossings retain their positions after an exact origin shift" (`LineRevolutionTests.cpp`).
- **Metamorphic regression:** "MR AUDIT-20261112-cone-polish far line origin keeps cone crossings". The shrunk seed's dyadic line and right cone.
  - MR-04: the line's origin moved back 2^19 directions.
  - MR-01: the three exact turns of the axes (identity and two cyclic), each with and without the translation (8, -4, 2).
- **Oracle:** the ideal cone's roots -191 ± 60 sqrt(10), computed independently of the kernel. Each point must also lie on the actual rounded cone's support in its own frame, and the parameters must differ by exactly 2^19.
- **Fault:** the `std::fma` point replaced by `Eval(line, t)`, the pre-repair code.

#### AUDIT-20261113-narrow-sphere-seam

- **Defect:** a plane's section of a sphere crossed the sphere's seam twice between adjacent samples of the trace's first scan. The fit clamped the excursion to the seam and never certified, so the trace failed with the sample budget. Found by coverage seed 13593231413355092879 (Actions run 43300991832, job 131471814362). Repaired in #286 (`SampleCuts` in `FitCertified`).
- **Direct regression:** "MR audit narrow sphere seam excursions remain complete certified circles" (`ImplicitTraceTests.cpp`). The audit record observed that this test fails without the repair, with the original sample-budget error. The observation was made before a workspace reset and was not replayed.
- **Metamorphic regression:** "MR AUDIT-20261113-narrow-sphere-seam plane section stays one whole circle". The shrunk seed's sphere and plane.
  - MR-01: a cyclic turn of the axes and the translation (8, -4, 2).
  - MR-02: both orders.
  - MR-03: scale 2 with the context.
  - MR-04: the sphere's chart turned a quarter about its axis.
- **Oracle:** the section's closed form. The circle lies in the plane, about the foot of the centre, with radius sqrt(1 - (5/8)^2). Every sample must be on it, the branches must sweep exactly one turn, and they must make one connected curve. Each branch must also be certified, and its pcurves must lie on both surfaces.
- **Fault:** the seam cuts found by the fitting samples are ignored (`if (false && !cuts.empty())`).

#### CI-272-ellipse-plane

- **Defect:** a conic is searched in quarter-turn arcs. A crossing exactly on a cut or on the seam was seen undecided from both arcs and refused: "unresolved root existence or multiplicity". Found by RapidCheck seed 17041817369036985825 (macOS Actions run 43273260370; [replay](ellipse-plane-ci-replay.md)). Repaired in #300 (`e8080c8aeb5490f6622f11ab52c68b885916fafd`): `ConicRootsOf` searches again with the cuts moved.
- **Direct regression:** "MR macOS ellipse plane shrunk linearization fixture resolves" (`EllipsePlaneReplayTests.cpp`).
- **Metamorphic regression:** "MR CI-272-ellipse-plane crossings on the seam and a cut are decided". The shrunk ellipse and the plane through its centre.
  - MR-01: the three exact turns of the axes, each with and without the translation (8, -4, 2).
  - MR-03: scales 1/4, 1 and 4, with the context.
  - MR-04: the ellipse reversed, and its frame turned half (t' = t - pi).
- **Oracle:** the plane's equation on the ellipse is b sin(phi) sin(t) = 0. So there are two transverse crossings, at centre ± a x, independently of the kernel. Every reported point must be decided, crossing, inside the curve's domain, on the curve and on the plane, at one of those two points, and both points must be reported.
- **Fault:** the second search disabled (`if (false && !AllDecided(found))`).

#### Batch 3 fault replay

- **Command (CI):** `python3 tools/mutate.py --catalogue tools/mutation/intersection-backfill*.json --work build/metamorphic`. This is the *Intersection metamorphic evidence* workflow on this pull request.
- **Run:** [#330's *Eight historical faults (GCC 16)* job](https://github.com/onurtuncer/Stereon/actions/runs/45064572675/job/137337734813), at `07780ba` (main 898ed69 with test-only changes), `gcc:16.2` container, `gcc16-debug`.
- **Baseline:** green; every selected test passed before any fault was applied.
- **Red:** killed in both runs by the named MR test, with the declared diagnostic found:
  - `AUDIT-20261112-cone-polish`;
  - `AUDIT-20261113-narrow-sphere-seam`;
  - `CI-272-ellipse-plane`;
  - `SSI-translated-boundary-search-shift`.

  The 24 earlier faults were killed in the same run too (score 28 / 29).
- **Survived:** `SSI-translated-boundary-edge-halving` (`EdgeHalvings` = 0). The MR test still passes without the halving, as long as the shift is kept. Like PR-159's `BisectSeam` candidate, this root cause is not observable alone with these fixtures. It was removed from the catalogue, so that the workflow tests only faults some test detects, and it is recorded here. No test was loosened. At that revision the record stayed `MR implemented`. The one-edge parabola below supplies the missing evidence with the search shift enabled.
- **Green:** the same tests passed in every CTest job of #330: [Linux `gcc16-debug`](https://github.com/onurtuncer/Stereon/actions/runs/45064572641/job/137337734757), [Linux `gcc16-asan`](https://github.com/onurtuncer/Stereon/actions/runs/45064572641/job/137337735042), [Windows `gcc16-debug`](https://github.com/onurtuncer/Stereon/actions/runs/45064572891/job/137337735813), [Windows `gcc16-windows`](https://github.com/onurtuncer/Stereon/actions/runs/45064572891/job/137337735659), [macOS `gcc16-debug`](https://github.com/onurtuncer/Stereon/actions/runs/45064572663/job/137337734763) and [coverage](https://github.com/onurtuncer/Stereon/actions/runs/45064572648/job/137337734570).
- **Review:** this agent's review of the fixtures, oracles and runs above. On it, `AUDIT-20261112-cone-polish`, `AUDIT-20261113-narrow-sphere-seam` and `CI-272-ellipse-plane` are `verified`. The owner's review is pending.

#### SSI-translated-boundary edge-halving replay (2026-11-25)

- **Historical replay.** On Linux with GCC 16.2, the original direct fixture from
  `75722ba82` on its parent `cbece058f` failed at (64, 0, 0), with one failure
  in 1067 assertions. On main `1fc0086` the direct fixture and the existing
  far-hyperbola MR passed all 7036 assertions.
- **Edge oracle.** In local coordinates the patch is
  `(s*u, s*v, s*(u + (v-1/2)^2 - 1/16))` on the unit square, with
  `s = 2^-12`. Its section by the wider plane z = 0 is one arc,
  `u = 1/16 - (v-1/2)^2`, with both ends on u = 0 at v = 1/4 and 3/4.
  The regression requires one certified branch, both distinct endpoints
  in either traversal order, and sampled residuals against both equations.
  It repeats at x offsets 0 and 1 and in both operand orders. The pair
  tolerance is 5e-5 and the world-coordinate oracle tolerance is 1e-8.
- **Edge red.** On [fixture revision 93e1411](https://github.com/onurtuncer/Stereon/commit/93e1411c2b3b55ea7338b1f29ae4931285884743),
  replacing only `EdgeHalvings = 8` with `0` made the public intersection
  refuse the arc. Two independent runs failed at `result.has_value()` with
  the search shift enabled. The named test is “MR SSI-translated-boundary
  resolves both roots on one patch edge”.
- **Edge green.** With halving enabled, the same regression passed all 308
  assertions in GCC 16.2 debug on Linux. All three targeted tests passed
  7344 assertions. The clean-worktree campaign on `93e1411` with CMake
  4.4.5 scored 2/2 critical faults killed. Both failed twice with their
  declared test diagnostics; no fault survived or was inconclusive.
- **Edge review.** The two endpoint ordinates follow from the patch
  equation at u = 0. Checking both prevents a duplicate endpoint from
  passing. The translated and exchanged cases use the same geometric oracle.
  The far-hyperbola regression and the one-edge regression detect the
  search-shift and edge-halving faults separately. Replay from the fixture
  revision with
  `python tools/mutate.py --catalogue tools/mutation/intersection-backfill-ssi-translated-boundary.json`.

## Inventory audit after #263

- **Date:** 2026-11-22
- **Branch:** `verify/metamorphic-inventory-audit`, on `main` 898ed69
- **Reviewer:** the v0.2 verification agent (Agent E). This is the agent's triage, not the owner's `inventory_review`.
- **Scope:** every merge after #263 (`cbeddfb`) that touches `src/intersect/` or the certification it relies on (`src/surface/Certify.cpp`), up to #323. Each was read from its merged description and diff summary. Those already in the inventory under another ID are SSI-translated-boundary (#279), the two AUDIT records (#286) and CI-272-ellipse-plane (#300).

A lead is a **defect** when the kernel, inside its documented supported domain, answered wrongly, or refused what the ledger said it answers. That matches the precedent of PR-159 and CI-272, refusals that were counted. A lead is **not a defect** when the earlier outcome was a refusal or an undecided report that the ledger itself listed as a limitation, or when the PR changes no result.

#### Confirmed defect leads, recorded untriaged

Each needs a fixture, a metamorphic regression and a replayed fault before it can close. None has them yet.

| Record | Merge | What was wrong before the repair |
| --- | --- | --- |
| PR-271 | `ca74138` | A B-spline's contact of order four with a line or conic was refused, although the ledger promised touches "reported bracketed, never dropped". Found by the planar corpus. |
| PR-289 | `51906ee` | Skew cylinders were answered in one operand order and refused in the other: an MR-02 failure on `analytic-ssi/cylinder-cylinder/skew-tilted-non-unit`. |
| PR-301 | `0b8f782` | A conic's crossing brought back across a full turn's seam could be placed at about -1e-17, outside the conic's domain: a wrong answer, about 2 % of property runs. |
| PR-304 | `814391b` | `ClosestApproach` placed its free parameter off the rounding it should reach. |
| PR-312 | `9438a7e` | `ClosestApproach` returned A 1.5e-7 short of its end when B landed on its own, a distance 1.4e-7 above an independent oracle (Nightly seed). |
| PR-315 | `eeb3367` | Eight far-placed `nurbs-ssi` reference fixtures were refused by `CertifyDeviation`. |
| PR-316 | `8a1a760` | Two B-spline curves: a proved crossing bracketed so that the JVP refused the primal point as not on both curves (seed 15971369584168069721). |
| PR-317 | `776459d` | Six `sphere-plane/near-pole` reference fixtures were refused. |
| PR-318 | `f0da86e` | A B-spline surface request ending on an interior knot answered branches outside the requested domain: a wrong answer. |
| PR-322 | `898ed69` | Skew-cylinder B-seam isolation by a fixed 256-step scan could miss seam crossings (#195). Replaced by fail-closed interval isolation. |
| PR-323 | `d9f9207` | `CertifyDeviation`'s B-spline partial enclosure had a rounding margin below the rounding it covers over short knot spans: an unproved certificate, though no wrong enclosure was observed. |

#### Leads that are not intersection defects

- **PR-267** (`3c0ee49`): two B-splines' crossings within rounding of a knot were reported undecided, never dropped. That was the ledger's listed limitation ("Remaining: ... undecided"). The PR decides them, a capability extension.
- **PR-281** (`48b3fd9`): a tangency of two B-splines and exactly refined copies were corpus cases marked `may_fail_closed`, refused by contract. The PR answers them.
- **PR-295** (`1b7b528`): the final diff has no kernel change. Its repair was dropped in favour of #300, whose defect is CI-272-ellipse-plane. It adds a replay fixture only.

No record is needed for the features and refactors in the same range. They add capabilities, products or speed, and the outcome they replace was a documented `NotImplemented` or no result change: #269, #273, #277, #282, #287, #288, #292, #296, #305, and before #263 #249, #253, #255, #256, #257.

#### What the audit does not establish

It reads merged pull requests only. Closed issues, failing corpus entries and regression-test comments since #263 were not audited, so `inventory_complete` stays false. The owner decides whether the scope above, together with the earlier batches, is the complete historical audit.

## Batch 4: the defect leads of the inventory audit after #263

- **Date:** 2026-11-26
- **Branch:** `verify/metamorphic-defects-after-263`, on `main` 09c0505 (#342)
- **Reviewer and meaning of "verified":** as for batches 1 to 3: this agent's documented review, with the owner's review still pending.
- **Toolchain:** the shared disk was full, so the red and green evidence is native CI: the *Intersection metamorphic evidence* workflow (`tools/mutate.py` over every `intersection-backfill*.json`, GCC 16.2 on Linux, `gcc16-debug`) and the Linux, Windows and macOS CTest jobs of the pull request.

The audit above recorded eleven confirmed leads as `untriaged`. This batch takes eight of them, whose repairs are single root causes (PR-317: two) with a fault that can be put back at one site. For each: the repair's own reproducer as the direct regression, a metamorphic regression named "MR PR-<n> ..." in `tests/intersect/MetamorphicRegressionTests.cpp` with an independent oracle, and a fault reintroduction in `tools/mutation/intersection-backfill.json` (PR-317: `intersection-backfill-pr317.json`, one fault per root cause). Every fixture is deterministic; there is no random seed. A fault that preserves today's interfaces is an *adapted* reintroduction of the original fault, as for the earlier batches; each says below what it restores.

PR-301, PR-322 and PR-323 follow the eight. Their faults cannot all be put back at one site on current code, so their red is a replay on the revision before the repair, in native CI (*Historical replays of PR-301, PR-322 and PR-323*). PR-301 also has a fault in the catalogue. The batch ends with the leads merged after #323 (*Leads after #323*).

#### PR-271

- **Defect:** a B-spline's contact of order four with a line or a conic was refused ("the subdivision budget is spent"), although the ledger promised touches "reported bracketed, never dropped". Along the spline the other curve's implicit equation vanishes there with its first three derivatives, so no piece about the touch can be proved apart, monotone or convex in doubles. Found by the planar corpus (`osculating` cases). Repaired in #271 (`ca74138`): `IsolateContacts` proves a derivative of order three or more of one sign and reports the contact as an undecided, bracketed touch.
- **Direct regression:** "a B-spline's contact of order four with a parabola is reported as a touch" (`PlanarSplineOsculationTests.cpp`).
- **Metamorphic regression:** "MR PR-271 a B-spline contact of order four with a parabola is a touch".
  - MR-01: a quarter turn and the translation (1024, -512).
  - MR-02: both orders.
  - MR-03: scales 2^-10 and 2^10.
  - MR-04: the spline reversed, which keeps the touch at s = 1/2.
- **Oracle:** the quartic Bezier (3(2s - 1), 9(2s - 1)^2 + (2s - 1)^4) is y = x^2 + (x/3)^4 with exact poles, touching y = x^2 to order four at s = 1/2. One undecided tangential point there, within the fourth root of the coefficients' rounding.
- **Fault:** the search for a derivative of order three or more of one sign is skipped (`LowestSigned(piece, 3, ...)` replaced by no order), which is the pre-repair `Isolate`.

#### PR-289

- **Defect:** the skew cylinders of `analytic-ssi/cylinder-cylinder/skew-tilted-non-unit` were answered in one order and refused in the other ("a branch of the skew cylinders exceeds the sample budget"): an MR-02 failure. The curve was parametrised by the first operand's angle. Repaired in #289 (`51906ee`): the parametrising cylinder is chosen from the geometry, and the other one is tried after a refusal.
- **Direct regression:** "two skew cylinders are answered in either order, the parametrisation chosen canonically" (`SkewOrderTests.cpp`).
- **Metamorphic regression:** "MR PR-289 skew cylinders are answered in either order".
  - MR-02: both orders.
  - MR-01: the cyclic turns of the axes and the translation (8, -4, 2).
  - MR-03: scales 1/4 and 4, with the radii, the pair tolerance and the context.
- **Oracle:** the axes are 4 apart, between |3 - 2| and 3 + 2, so the section is one loop. Every sample of every certified branch must lie on both cylinders, by the closed-form distance from each axis.
- **Fault:** the curve parametrised by the first operand, with no second attempt (`return run(false)`). It restores both halves of the defect at once: one block replaces the canonical choice and the fallback, either of which alone would hide it.

#### PR-304

- **Defect:** `ClosestApproach` stopped a free parameter about 1e-8 short when the other rested on an end (shrunk macOS seed 9376910471991398838). `Refine` took a step only if phi = |A - B|^2/2 did not rise, and near the minimiser a correct step lowers phi below its rounding. Repaired in #304 (`814391b`): a step that shrinks the free gradient and raises phi by no more than its rounding is taken too.
- **Direct regression:** "closest approach places the free parameter to rounding when the other rests on an end" (`ClosestApproachTests.cpp`).
- **Metamorphic regression:** "MR PR-304 closest approach places the free parameter to rounding". The shrunk pair, B moved by 17 x 17 small steps.
  - MR-01: the cyclic turns of the axes and the translation (8, -4, 2), exact on the dyadic centres.
  - MR-02: both orders, the ends exchanged.
  - MR-03: scales 1/4 and 4, with the context; the angles do not scale.
- **Oracle:** B's parameter is the angle of A's end point projected into B's plane, atan2 in B's frame, to 1e-13.
- **Fault:** the descent test alone (`m.Phi <= l.Phi`).

#### PR-312

- **Defect:** `ClosestApproach` returned A 1.5e-7 short of its end when B landed on its own, a distance 1.4e-7 above the long double oracle (Nightly, 2026-11-17). The step that landed B on its end moved A by a rounding, and that ended the loop. Repaired in #312 (`9438a7e`): a step that changes the active set is not convergence.
- **Direct regression:** "closest approach moves the other parameter on once one lands on an end" (`ClosestApproachTests.cpp`).
- **Metamorphic regression:** "MR PR-312 closest approach steps on after one parameter lands on an end".
  - MR-01: the cyclic turns of the axes and the translation (8, -4, 2).
  - MR-02: both orders.
  - MR-03: scales 1/4 and 4, with the context.
- **Oracle:** the direct regression's long double oracle puts the minimum at the corner where both arcs end. The MR test requires that corner: both ends `Upper`, and the distance equal to the corner's closed-form distance to 1e-12 relative. It also checks independently that the corner is a first-order minimum on the box: phi does not decrease into either domain.
- **Fault:** the active-set clause removed from the convergence test.

#### PR-315

- **Defect:** eight far `nurbs-ssi` reference fixtures, in both orders, were refused by `CertifyDeviation`. Far from the origin the branches were certified against the input surfaces, whose rational enclosure widens every homogeneous component by a margin scaled to the largest pole: a floor of about |P| 1.5e-8, 0.007 mm at |P| = 6e5, at any span width. Repaired in #315 (`eeb3367`): each branch is certified on the surfaces moved exactly about the origin, then moved back with its rounding.
- **Direct regression:** "far-placed rational SSI branches certify against the input surfaces" (`SurfaceCorpusTests.cpp`).
- **Metamorphic regression:** "MR PR-315 far rational SSI branches certify against the input surfaces". The exact rational unit sphere and the plane z = 3/5.
  - MR-01: at the origin and at the corpus's far origin (480000, -300000, 210000), each also with the axes turned cyclically.
  - MR-02: both orders.
- **Oracle:** the circle of radius 4/5 about the axis in the plane, in closed form. Every sample of every certified branch must lie on it within the pair tolerance, with each bound within the resolution, and the branches together must make one curve going once round it.
- **Fault:** the branch moved back first and certified against the input surfaces, rebuilt exactly from the moved ones. This is the pre-repair certification at the current interfaces.

#### PR-316

- **Defect:** a decided crossing of two cubic Beziers was placed at the middle of brackets 8e-13 wide, the width at which Krawczyk's refinement on the subdivided pieces stalls. The curves' points at its parameters were 1.5e-13 apart, and `IntersectJvp` refused the primal's own point as stale (seed 15971369584168069721). Repaired in #316 (`8a1a760`): each decided crossing is narrowed on the splines' own poles.
- **Direct regression:** "two cubic Beziers' proved crossings lie on both curves to rounding" (`PlanarLinearizationTests.cpp`).
- **Metamorphic regression:** "MR PR-316 proved crossings of two cubic Beziers lie on both curves".
  - MR-01: a quarter turn and the translation (8, -4), exact on the dyadic poles.
  - MR-02: both orders.
  - MR-03: scales 1/4 and 4.
  - MR-04: either cubic reversed, t' = 1 - t.
- **Oracle:** each operand is evaluated at its reported parameter by de Casteljau in long double, independently of the kernel. The two points, and the reported position, must be within 64 units in the last place of the coordinates' size. There must be two decided crossings, at the identity placement's positions, placed.
- **Fault:** the tightening skipped.

#### PR-317

- **Defect:** the six `sphere-plane/near-pole` reference fixtures were refused: the plane z = 0.999 against the exact rational unit sphere, a small circle 2.56 degrees from tangency. Two root causes, both repaired in #317 (`776459d`):
  - *partials-about-patch:* a rational patch's partials were enclosed by the quotient rule about the origin, widened by the patch's distance from it times the width of the weight's derivative. No advance direction or boundary root was proved near the pole, and the search refused the pair as a tangential contact. Now they are enclosed about the patch's middle.
  - *roots-in-frame:* off the world axes (the corpus's tilt), the boundary roots were not isolated, and the search ran out of its budget. Now they are isolated in a frame about either patch's normal.

  #317 also proves an advance direction in such a frame. Its own ablation shows no fixture that needs that third change once the other two are in, so it has no fault here.
- **Direct regression:** "a plane a thousandth below a rational sphere's pole cuts its small circle" (`GenericSSITests.cpp`).
- **Metamorphic regression:** "MR PR-317 a plane near the pole of a rational sphere cuts its small circle".
  - MR-01: the identity, a cyclic turn of the axes, the corpus's tilt, and its far translation.
  - MR-02: both orders.
- **Oracle:** the circle of radius sqrt(1 - 0.999^2), in closed form, checked as for PR-315.
- **Faults:** `tools/mutation/intersection-backfill-pr317.json`, one per root cause:
  - `PR-317-partials-about-patch`: the partials enclosed about the origin;
  - `PR-317-roots-in-frame`: no isolation in a frame.

#### PR-318

- **Defect:** a request on a B-spline surface ending exactly on an interior knot kept the neighbouring Bezier patch whole. That patch meets the request only in the knot edge, and the cut there could not be split off. So intersections beyond the request were answered, outside the requested domain: a wrong answer. Repaired in #318 (`f0da86e`): such a patch is skipped, and a cut that cannot be enclosed fails.
- **Direct regression:** "a requested B-spline surface knot excludes curve crossings beyond it" (`CurveSplineSurfaceTests.cpp`).
- **Metamorphic regression:** "MR PR-318 a request ending on a surface knot excludes crossings beyond it". The degree-1 spline plane S(u, v) = (u, v, 0) with a u knot at 1/2, and the line along z through (x, 1/2, -1), for x = 1/4 and 3/4.
  - MR-01: the cyclic turns of the axes and the translation (8, -4, 2).
  - MR-03: scales 1/4 and 4.
  - MR-04: the surface's u reversed (u' = 1 - u, which exchanges the halves) and the line run downwards (t' = 2 - t).
  - MR-02 does not apply to that half: there is no exposed surface-curve operation.
  - The generic SSI half is #318's other regression: the exact rational unit sphere requested over u in [0, 1/2] or [1/2, 1], ending on its knot, against the plane z = 1/2. It runs at the identity, under a cyclic turn and under the corpus's tilt (MR-01), and in both orders (MR-02).
  - The [first fault replay](https://github.com/onurtuncer/Stereon/actions/runs/45896367150/job/140171146116) (`572e8c7`, before this half was added) showed the fault surviving the curve-surface half alone. The curve-surface path does not reach the restored patch. The sphere is the fixture on which the defect was found (#317's finding).
- **Oracle:** for the curve-surface half, the crossing S(x, 1/2) at the line's t = 1, in closed form. The half of the request that contains u = x must have exactly that point; the other half, which ends on the knot, must have none. For the sphere, half of the circle of radius sqrt(3)/2: every branch must stay inside the requested u and on the circle, and the branches together must turn pi.
- **Fault:** a cut that meets a patch only in its edge keeps the patch whole, and such a patch counts as overlapping the request. This is the pre-repair `Clipped` at the current interfaces, applied as one block over the two adjacent functions.

#### PR-301

- **Defect:** #300 decides a crossing on the seam of a full turn by searching again with the cuts turned, then bringing the roots past the upper end back a period. Each bracket brought back is widened by a few units for the rounding of 2 pi. It could therefore reach outside the domain, and the root (or its polish inside the bracket) could be placed just outside `Domain(curve)`. That is a wrong answer, in about 2 % of the translation property's runs (seed 152652381). Repaired in #301 (`0b8f782`): every bracket and root of the turned search is kept in the domain.
- **Direct and metamorphic regression:** "MR PR-301 seam crossings of a full turn stay in the conic domain". #301's own tests were property seeds, so this deterministic fixture is new. An ellipse (semi-axes 3/2 and 3/5) against the plane through its centre and its major axis: it crosses on the seam (t = 0) and on a cut (t = pi).
  - MR-01: 16 turns about an oblique axis, each with 8 translations, so that the rounding about 2 pi varies.
  - MR-03: scales 1/4 and 4 with the context.
  - MR-04: the ellipse reversed.
  - MR-02 does not apply: there is no exposed surface-curve operation.
- **Oracle:** the plane's equation on the ellipse is b sin(t) = 0, so the crossings are the centre -+ a x, in closed form. Each must be decided, transverse and on the curve, with its parameter in `Domain(curve)`.
- **Fault:** the clamp of `BroughtBack` removed. The fault replay shows whether it is observable on current code; the historical replay does not depend on it.

#### PR-322

- **Defect:** the skew-cylinder handler looked for the curve's crossings with B's seam by a fixed scan of 256 steps of A's angle. Two crossings within one step were missed together, and the branch was not cut at B's seam, so the answer lost both seam crossings (issue #195). Repaired in #322 (`898ed69`): fail-closed interval isolation.
- **Direct regression:** "two B seam crossings inside a former scan cell are both isolated" (`CylinderCylinderSkewTests.cpp`).
- **Metamorphic regression:** "MR PR-322 two B seam crossings in one former scan cell are both isolated". #322's fixture:
  - A, of radius 3 about z;
  - B, of radius 2 along x through (0, 4, 0), with frame (0, -q, s), (0, -s, -q), where q = (1 + delta)/2 and delta = 1e-4 or 1e-7;
  - A's angle requested on an asymmetric window about the fold.

  The relations:
  - MR-01: the cyclic turns of the axes and the translation (8, -4, 2).
  - MR-02: both orders, with the boxes exchanged.
  - MR-03: scales 1/4 and 4 with the context.
  - MR-04: B's handedness reversed, which reverses its angle.
- **Oracle:** B's seam generator meets A at (+-sqrt(9 - y^2), y, 2s), with y = 4 - 2q, in closed form. Both points must be ends of certified branches. Every sample must lie on both cylinders, by the distance from each axis.
- **Fault:** none in the catalogue. #322 replaced the scan with a new isolation, and no single site of the new code restores the scan. The red is the historical replay.

#### PR-326

- **Defect (issue #336):** two coaxial profile circles touching within the pair tolerance were declared tangent at the radical axis's point. Near an internal tangency that point lies gap (R_p + R_q)/2d from the contact, outside both spheres by more than the tolerance. Example: a sphere of radius 1.5 inside one of radius 2, 0.005 from touching, at tolerance 0.01, put the contact at (2.0152, 0, 0). The OCCT corpus comparison of #326 found it, and #326 repaired it (`8f955e0`, `TouchAlong`). The audit after #263 read #326 as a comparison only, so this record is new.
- **Direct regression:** "spheres touching within the tolerance touch between their surfaces" (`ElementaryPairsTests.cpp`).
- **Metamorphic regression:** "MR PR-326 spheres touching within the tolerance touch between their surfaces". Gaps of 1/2 and 9/10 of the tolerance, inside and outside, internal and external.
  - MR-01: the cyclic turns of the axes and the translation (8, -4, 2).
  - MR-02: both orders.
  - MR-03: scales 1/4 and 4, with the tolerance and the context.
- **Oracle:** the contact lies within half the gap of each sphere, by the closed-form distance from each centre.
- **Fault:** the radical axis's point restored for the tolerant tangency.

#### Historical replays of PR-301, PR-322 and PR-323

Each replay is a scratch branch that is never merged. For PR-301 and PR-322 it is the revision before the repair with only the regression under test added. For PR-323 it is current `main` with the pre-repair code put back. The *Linux* workflow was dispatched on each branch.

- **PR-301:** `verify/replay-pr301-red` at `bfcbe56`: `8ba9d6d`, the parent of #301, plus the MR test.
  - In [`gcc16-debug`](https://github.com/onurtuncer/Stereon/actions/runs/45858817643/job/140044149022), "MR PR-301 seam crossings of a full turn stay in the conic domain" failed 72 of its 12288 assertions.
  - Every failure is `Contains(Domain(curve), point.Parameter)`, with the parameter reported as 6.28319, past the domain's end. They occur at scales 1/4, 1 and 4 and in both charts.
  - The run's other failure, "conic derivative fixtures retain unresolved domain ends and tangencies", was already failing on `8ba9d6d`; #304 repaired it.
- **PR-322:** `verify/replay-pr322-red` at `44fa41b`: `d9f9207`, the parent of #322, plus the MR test.
  - In [`gcc16-debug`](https://github.com/onurtuncer/Stereon/actions/runs/45858785888/job/140044036294), "MR PR-322 ..." was the only failure, 144 of 57222 assertions.
  - Every failure is `atPlus` or `atMinus` at delta = 1e-7: both seam crossings missed together, in all 18 placements, both handednesses and both orders.
  - Every delta = 1e-4 case passes, and so do the checks of the samples on both cylinders.
- **PR-323:** `verify/replay-pr323-red` at `11ff162`: `main` 09c0505 with `ValueAt`'s margin again taken from the surface's largest homogeneous pole for every spline, partials included, as before #323.
  - In [`gcc16-debug`](https://github.com/onurtuncer/Stereon/actions/runs/45858979933/job/140044715911), the three oracle tests of #323 failed, and nothing else did: "a B-spline surface's jet covers its rounding on every span, polynomial or rational", "a partial that vanishes by cancellation is enclosed to its poles' rounding" and "a partial's enclosure covers the worst rounding of its own poles".
- **Green:**
  - The PR-301, PR-322 and PR-326 MR tests: this pull request's CTest jobs (see *Batch 4 fault replay*).
  - The three PR-323 tests passed on `main` 09c0505 in [Linux `gcc16-debug`](https://github.com/onurtuncer/Stereon/actions/runs/45829925674/job/139943711808).

#### PR-323: metamorphic exception

- **Defect:** `CertifyDeviation`'s B-spline enclosure (`src/surface/Certify.cpp`, `ValueAt`) widened a partial's value by a margin sized to the surface's largest pole. A partial's poles are a degree over a knot span times differences of the surface's, so over short spans the margin fell below the rounding it covers. The certificate was unproved there, though no wrong enclosure was observed. Repaired in #323 (`d9f9207`): the margin comes from the poles each sum combines.
- **Direct regression and oracle:** "a B-spline surface's jet covers its rounding on every span, polynomial or rational" (`tests/surface/CertifyTests.cpp`).
  - Its oracle is exact in doubles: Cox–de Boor from the definition, on integer poles, dyadic knots and dyadic parameters.
  - It covers 518400 components: degrees 1 to 4, rational and polynomial, spans from 2^10 to 2^-45.
  - Two tests accompany it: "a partial's enclosure covers the worst rounding of its own poles" and "a partial that vanishes by cancellation is enclosed to its poles' rounding".
- **Relations considered:**
  - MR-01 (rigid motion), MR-02 (exchange) and MR-03 (scale) act on an intersection's operands. This defect is in a point enclosure, and it made no intersection answer wrong in any run measured. Two intersection runs would agree and expose nothing.
  - MR-04 (equivalent parameterization) is what the oracle already applies, by scaling the knot vectors by powers of two from 2^10 to 2^-45. That is the rechart under which the defect appears, checked against the exact value rather than against another run.
  - A metamorphic test of intersections would be weaker than the exact oracle, so none is invented (ADR-0053 section 3).
- **Red:** the historical replay above.
- **Green:** the three tests on `main` 09c0505.
- **Review:** this agent's review of the margin argument in `ValueAt`'s comment, the oracle and the runs above. On it the record is `exception reviewed`. ADR-0053 asks for an exception to be reviewed explicitly; the owner's review is pending.

#### Leads after #323

The inventory audit after #263 stopped at #323. These merges since then, up to `main` 09c0505, touch `src/intersect/` or `src/surface/Certify.cpp`. Each was read from its description and diff.

- **PR-326** (`8f955e0`): a defect, recorded above. Its kernel commit is inside an OCCT comparison PR.
- **PR-334** (`431d9e9`): not a defect. A crossing of two B-splines at a spline's own end, or at a corner, was reported undecided. That was the ledger's listed limitation (2D curve–curve, *Remaining*), and the PR decides such crossings. The point it replaced was labelled `Crossing`, but undecided points carry no claim about their contact.
- **PR-337** (`7cb4564`): not a defect. Transverse nodes on knot lines and seams were refused with `IntersectionFailed`, the generic SSI's documented limitation (step 3 of the critical path). The PR answers them.
- **PR-339** (`3b4afef`): verified chart-continuity defect in the skew-cylinder fitter. See the historical contract and replay below.
- **PR-342** (`09c0505`): already recorded, as AUDIT-20261122-plane-cylinder-periodic-box.
- #310 (which closed issue #309) changed tests only. #333, #338 and #341 add a corpus, Booleans and evidence, and change no intersection result. None of them needs a record.

#### Batch 4 fault replay

- **Command (CI):** `python3 tools/mutate.py --catalogue tools/mutation/intersection-backfill*.json --work build/metamorphic`, the *Intersection metamorphic evidence* workflow on this pull request.
- **Run:** [#348's *Eight historical faults (GCC 16)* job](https://github.com/onurtuncer/Stereon/actions/runs/45974949086/job/140424249217) at `95755ef` (`main` with test-only changes), `gcc:16.2` container, `gcc16-debug`.
- **Baseline:** green. Every selected test passed before any fault was applied.
- **Batch 4 red:** each fault was killed in both runs by its named MR test, with the declared diagnostic found: PR-271, PR-289, PR-301, PR-304, PR-312, PR-315, PR-316, PR-317-partials-about-patch, PR-317-roots-in-frame, PR-318 and PR-326. The 30 earlier faults were killed in the same run: score 41 / 41, none survived. PR-322's red is its historical replay above.
- **Batch 4 green:** at the same head, every MR test of this batch passed in each native CTest job of #348: [Linux (gcc16-debug)](https://github.com/onurtuncer/Stereon/actions/runs/45974949273/job/140424249958), [Linux (gcc16-release)](https://github.com/onurtuncer/Stereon/actions/runs/45974949273/job/140424250089), [Linux (gcc16-asan)](https://github.com/onurtuncer/Stereon/actions/runs/45974949273/job/140424249883), [Windows (gcc16-debug)](https://github.com/onurtuncer/Stereon/actions/runs/45974949160/job/140424250795), [Windows (gcc16-windows)](https://github.com/onurtuncer/Stereon/actions/runs/45974949160/job/140424250506), [macOS (gcc16-debug)](https://github.com/onurtuncer/Stereon/actions/runs/45974949208/job/140424250365).
- **Batch 4 review:** this is this agent's review of the fixtures, oracles and runs above. On it, PR-271, PR-289, PR-301, PR-304, PR-312, PR-315, PR-316, PR-317, PR-318, PR-322 and PR-326 are `verified`. The owner's review is pending.

#### PR-339 skew-cylinder chart continuity

**Historical contract.** The [manual at the defective revision](https://github.com/onurtuncer/Stereon/blob/e84edffd785590fd01f393303b5e3721dc4f9c2d/docs/intersection.rst#L141-L195) includes skew cylinders, caps fitting at 2^14 samples, and permits failure at the cap. It separately permits unresolved or inseparable seam cuts and a 2^16 seam-isolation budget. The two reported fixtures triggered the fitting cap. Neither needed a cut at the offending seam. An inconsistent periodic representative introduced a jump into a smooth pcurve. The refusal was safe; the chart discontinuity that caused it was a defect inside the supported geometry.

**Fixture and oracle.** In `CylinderCylinderSkewSeamTests.cpp`, the parametrising cylinder's seam misses the far fold by 2^-53 radians. In the second fixture B's seam lies above the loop by the 2^-52 input offset with A's seam just outside the fold. Rotating the parameter frames preserves the world-space intersection (MR-04). The tests require an answer at these offsets. The common Check routine compares independently bracketed crossings against the algebraic count and checks the returned pieces against the fixed reference curve. The separate concern about the FoldTouch threshold had no reproduced wrong branch and carries no independent red claim.

**Replay and review.** [PR #339's red-on-main table](https://github.com/onurtuncer/Stereon/pull/339) records the two named failures with the solver from e84edff, GCC 16.2 debug on Windows. This is the PR author's local replay. The repaired PR head 0d6d7b1 passed [Windows debug](https://github.com/onurtuncer/Stereon/actions/runs/45570577102/job/139055930299) and [Linux debug](https://github.com/onurtuncer/Stereon/actions/runs/45570577088/job/139055929913). Review covers the two reproduced chart-dependent refusals. The documented fitting and isolation limits continue to apply.
