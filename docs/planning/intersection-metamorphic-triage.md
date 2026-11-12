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

The replay for all three is the campaign's: `python tools/mutate.py --catalogue tools/mutation/intersection-backfill.json --only <ID>` on committed HEAD. It builds a disposable worktree, requires the named test to pass on the repaired code, then restores the fault and requires the named test to fail twice. Every fixture is deterministic; there is no random seed.

### PR-200

- **Defect:** `LineSurface.cpp` placed a line's crossings of a sphere or cylinder from coefficients formed at the line's own origin. Those carry eps |w|² for an origin |w| away, so from 24 units the crossings lay 1.4e-13 off the cylinder, about 600 units in the last place. Repair `e8db99b53bf4396bf1f48d92b195f7b06452b85c` polishes each root from the line re-centred at its rounded point.
- **Direct regression:** `tests/intersect/LineSurfaceTests.cpp`, "a crossing from a far line origin lies on the quadric to rounding".
- **Metamorphic regression:** "MR PR-200 far line origin crossings on the quadric".
  - MR-04: the origin moved back by m ∈ {0, 8, 64} directions, exact on these dyadic inputs.
  - MR-01: the three cyclic axis permutations, exact.
  - Both the sphere and the cylinder.
- **Oracle:** the closed-form distance of each crossing from the sphere's centre, or from the cylinder's axis, minus the radius. It must be within 32 eps |p| plus 4 eps |t| |d|, the line's own parameter grid. The parameter must move by exactly m (to 1e-12 m). The kernel's projection is not used.
- **Fault:** the polish is removed (`t = rounded`).

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

## Red and green record

Filled in from the replay run below.
