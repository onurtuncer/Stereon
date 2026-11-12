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

After batch 2 these two remain. Each is a confirmed or likely intersection defect with a direct regression already on `main`. PR-194 is `fixture captured`, pending its split into sub-IDs; PR-239 stays `untriaged`. None is closed. Closing needs a metamorphic test and a fault that is actually killed, and this batch has not produced either for these leads.

### PR-194

At least five root causes in one PR (`b74d6e23ff`), each with its own direct test. PlanePlane reported a rotated bounded coincidence as an axis-aligned overlap, and the squared cross norm of near-parallel planes underflowed. PlaneCylinder gave bounded plane patches uncut conic branches. SplineSurface reported unproved tangencies, near misses and knot roots as contacts. Linearization (JVP/VJP) accepted stale geometry. GenericSSI let direct NURBS construction bypass the factory weight bounds. Open: ADR-0053 asks for one traceable entry per root cause, so splitting this lead into sub-IDs (as PR-230-pole) is a ledger-structure decision for the coordinator.

### PR-239

An owner semantic decision (2026-11-08, `0338cb3e44`): 3D lines now meet within the resolution, as degree-1 B-splines do. The former exact-meeting rule was documented behaviour, not an implementation fault, so whether this lead is a defect or a contract change is the owner's classification. A fault reintroduction would be the whole former `LineLine.cpp`, not one statement. Its property test already compares each line with its degree-1 B-spline (a representation relation, MR-04). Open pending that classification.

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

Filled in from the replay of the seven batch 2 faults.
