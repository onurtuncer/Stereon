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

Each of these is a confirmed or likely intersection defect with a direct regression already on `main`. Each is recorded as `fixture captured`, except PR-239, which stays `untriaged`. None is closed. Closing needs a metamorphic test and a fault that is actually killed, and this batch has not produced either for these leads.

### PR-112

Two intersection root causes in one PR (`93fd35313f`). First, coincident cylinders ignored the requested height domains, and their arcs were not split at the periodic seam (`CylinderCylinder.cpp`). Second, the plane-cylinder ellipse pcurve work was not bounded before allocation (`PlaneCylinder.cpp`). Direct tests: "coincident cylinders intersect their requested height domains", "coincident cylinder arcs split at the periodic seam", "coincident cylinder overlap coverage agrees with an angular membership oracle" and "ellipse pcurve work is bounded before converting or allocating its samples". Open: the domain defect needs an MR test (MR-01, and MR-04 with a seam-shifted chart) and a fault reintroduction. The work bound is a resource contract whose MR applicability must be reviewed. The topo rollback-identity part is not an intersection defect.

### PR-120

The merged diff (`470cbf7be4`) is tests only, and its description says no production change was needed for its four regression cases. The sector clipping repair it describes is in its parent `3a2ba895be415efe0198af878a2630738fd1a097` (#118), which changed `CylinderCylinderSkew.cpp`. Open: the defect is confirmed, but the repair revision must be pinned to the #118 hunk before a fault can be reintroduced. The four direct tests already run in both operand orders (MR-02).

### PR-159

Skew cylinders internally tangent on both seams were refused (`35f3d23709`): a sliver piece under the merge gap could not be fitted. This was a refusal of a supported case, not a wrong answer. Open: it needs an MR test (MR-02, and MR-04 with an exact seam rotation) and a fault restoring the old merge gap in `CutsOf`.

### PR-194

At least five root causes in one PR (`b74d6e23ff`), each with its own direct test. PlanePlane reported a rotated bounded coincidence as an axis-aligned overlap, and the squared cross norm of near-parallel planes underflowed. PlaneCylinder gave bounded plane patches uncut conic branches. SplineSurface reported unproved tangencies, near misses and knot roots as contacts. Linearization (JVP/VJP) accepted stale geometry. GenericSSI let direct NURBS construction bypass the factory weight bounds. Open: ADR-0053 asks for one traceable entry per root cause, so splitting this lead into sub-IDs (as PR-230-pole) is a ledger-structure decision for the coordinator.

### PR-236

Two ElementaryPairs defects (`252c1350cf`). Concentric spheres within the tolerance were an Exact overlap at distance 0. A sphere on a cylinder's axis, with radii within the tolerance, failed the whole intersection. Open: each needs an MR test (MR-02, MR-03 with the tolerance scaled, MR-01 by exact permutation) and a fault reintroduction. The first is a contract (evidence kind) defect, the second a refusal.

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
