# Priority intersection metamorphic campaign

This campaign implements tests and fault reintroductions for eight identified
faults. It does **not** close the historical inventory or claim native CI has
passed. Records stay `MR implemented` until the baseline, fault-detection logs,
and independent review are linked. The matrix covers only the domains below,
not every behavior of each handler.

## Replay and effective context

All fixtures are deterministic; no random generator or seed is used. The
repository revision fixes `Context::Default()` in `src/core/Context.cppm`.
`Placement::Policy()` selects serial execution and scales **only** linear
resolution and absolute tolerance ceiling by the geometric scale. Defaults:
linear resolution 1e-5 mm, angular resolution 1e-11 radians, absolute ceiling
1e-2 mm, tolerance growth allowed, geometric validity, no sinks, default pool
configuration (unused by serial execution). Other policies remain those of
that revision's default Context. PR-244 uses the serial default Context and
passes its resolution explicitly as pair tolerance.

Placements: scales {1/4, 1, 4}; three proper cyclic Cartesian frames;
translations {(0,0,0), (8,-4,2)} mm. Thus 18 declared placements, with no
rejection sampling. Far origins shift line parameter by {2^20, 2^27}; circle
charts rotate by pi/2 and subtract pi/2 from the parameter. PR-244 deliberately
uses operand permutation only: translating a 1e-15 displacement could erase
it by rounding and violate the transformation's premise. Failure output names
fixture, placement, scale, origin shift, chart and operand order as applicable.

```sh
cmake --preset gcc16-debug
cmake --build --preset gcc16-debug --target stereon_intersect_tests
ctest --preset gcc16-debug -R '^intersect: MR PR-' --output-on-failure
python tools/mutate.py --catalogue tools/mutation/intersection-backfill.json --work build/metamorphic
```

The final command uses a disposable checkout of **committed HEAD**. It first
builds and passes each selected test on the repaired code, then restores one
historical fault at a time. A named assertion must fail twice. Empty test
selections, production compile errors, timeouts, crashes, differing rerun
failures and unknown fault IDs do not qualify as evidence. The logs preserve
baseline build/tests and both fault runs. The JSON report records revision,
compiler, CMake and platform; upload the complete logs with that report.

## PR-205

Unit sphere and x-directed line starting at (-2,0,0), then shifted far back
along its own direction. Exact independent answers are (-1,0,0) and (1,0,0),
transformed by each placement. Both runs must return two decided crossings,
no overlaps, and points evaluating on both supports. Parameters differ by
the origin shift, independent of geometric scale. MR-01/03/04.

Fault: restore rounded quadratic coefficients and discriminant instead of
`QuadraticOf`, the relevant source change of repair
`3a147886ddc7ad652dacc98d2272cf6f98a0a3a5`. This previously collapsed roots.
The catalogue deliberately targets the sphere path; it does not claim to
cover every planar/conic handler repaired in that PR.

## PR-206

Coplanar radius-1.5 circle clipped by [-2,2] x [-1,1]. Singleton angle 0 is
inside: exactly one singleton overlap, no isolated points; angle pi/2 is
outside: an empty result. Rotating the circle chart remaps the angle without
changing either point. All placements scale the plane's UV box because its
axes remain unit length. Evaluate the returned pcurve on the plane and check
against the independently known point. MR-01/03/04.

Fault: remove the singleton-domain branch of `InsideBox`, restoring the
relevant omission fixed by `806217bae3c15b93b5fcdb2255534a5fec11b727`.

## PR-218

Cone of radius 1 at z=0 and half-angle pi/4; horizontal line at z=1/2.
The independent crossings are x=+-1.5, y=0, z=1/2. Near and far origins must
both resolve these crossings under all placements. Same independent-support
and parameter checks as PR-205. MR-01/03/04.

Fault: set the recentering parameter to zero, so the enclosed quadratic is
again formed at the far origin. This is a controlled reintroduction of the
cancellation fixed by `b265d641e61cc6b8a43fe8127e46394d1cbec108`, adapted to
the current helper interfaces, not a claim that the entire old revision runs.

## PR-244

Sphere radius 2 at (offset,0,4), offsets {1e-15,1e-7,5e-6}, against the cone
at the origin of half-angle pi/6, restricted to [0,2pi] x [-8,8]. With pair
tolerance 1e-5 the contract permits curves or an explicit unresolved failure,
never isolated meridian points masquerading as the near-axis contact.
Operand permutation must preserve success/failure; on success independently
check pcurve consistency on both supports in both orders. MR-02.

Fault: replace `rho > tol` with `rho > 0`, as before repair
`aae5992bfa7d14150624bfd3e8a3189ee7ab267f`. This is **fail-closed contract
evidence**, not proof that unresolved geometry was solved. No general
sphere/cone robustness pass percentage is inferred from it.

## PR-247

Coaxial cones with equal apex and unequal half-angles pi/4 and pi/4+2^-42.
The independent premise is unequal support sets. In this documented narrow
angle band the current contract must fail with `IntersectionFailed`, not
publish tolerant coincidence. Check both operand orders and all placements.
MR-01/02/03. This is **fail-closed contract evidence**, not resolved geometry.

Fault: bypass `UnresolvedConeAngles`' rejection, preserving the helper call
so an unused-function warning cannot masquerade as detection. This reverses
the relevant behavior of `9514e8d7fbe0517b193a319b61bf962b266e56ee`.

## PR-230

The exact-predicate repair uses the evaluator's **rounded** third axis,
not an exact cross product of the stored X/Y axes. Start with sqrt(1/2)
components in the XY frame and semi-axes (3,2,1). For both handednesses and
both poles, require one tangent point; a plane one representable step beyond
that pole must be separate. Exercise operand exchange, three cyclic axis
permutations and dyadic scaling (9 placements). Check returned parameters
on both surfaces independently against the known evaluator pole.
MR-01/02/03. No translation: rounding a translated tangent plane can change
whether the input doubles actually describe tangency. This restriction is a
precondition of this exact-contact fixture, not a general handler exemption.

Fault: restore the exact-field cross product in `AlongFrame`, reversing the
axis part of `124b5211cad558f837ebf378df13bcd2de6f659c`.

## PR-230-pole

An axis-aligned ellipsoid of semi-axes (3,2,1), in all 18 placements, touches
a plane at each pole. Longitude boxes [1,2] and [1+2pi,2+2pi] must both retain
the point: longitude is arbitrary at a pole. Check one isolated point, no
branches, membership of returned UV in both boxes, and independent evaluation
on both supports, in both operand orders. Restricting latitude to [-pi/4,pi/4]
must instead prove separation. MR-01/02/03/04; plane UV boxes scale with length,
ellipsoid angular boxes do not.

Fault: remove pole-specific longitude selection and always call
`InPeriodicDomain`, reversing the second root cause in PR-230. The ledger
retains PR-230 and adds PR-230-pole so both independent faults need evidence.

## Ellipsoid combined applicability

The Plane/Ellipsoid and Ellipsoid/Plane MR-01/02/03 cells cover the union of
the two fixture families above, with their distinct preconditions. The
rounded-axis case deliberately excludes translations; the pole-box case
includes them. MR-04 covers the pole-box longitude change only. These entries
do not assert general coverage of oblique sections, seams or near-tangencies.

## PR-241

PR-240 and its re-land PR-241 share the missing domain-boundary search.
Keep both ledger entries; they may share this test and eventual replay proof.
Two degree-one splines extend away from their nearest endpoints (1,0,0) and
(1+g,g,0), with g=2^-20 and slope 2^-23 on the second segment. Their endpoint
distance sqrt(2)*g is below resolution; exactly one meeting must be found.
With g=2^-15 the same construction is separated. All lengths and the Context
scale together across all 18 placements. Reverse both curves, remap their
parameter domain from [0,1] to [2,6], and exchange operands. MR-01/02/03/04.

Check the independently known endpoint parameters and points on each curve.
The representative position lies on operand A: exchanging operands changes
it to the other endpoint, so byte equality would be wrong. This covers the
end/end boundary defect; end/interior and general higher-degree extrema
remain outside this fixture's coverage.

Fault: disable `SearchEnds` for nonempty inputs while retaining a syntactic
reference to its helper (avoiding unused-function build failures). This
reintroduces the relevant omission addressed by
`51a2795f39031179e3bcd586505e14be57e676f5`, not an entire historical build.

## Remaining qualification

Every campaign record stays `MR implemented`, not `verified`, until native
GCC 16 baseline/fault logs and review are linked. Other historical leads,
unreviewed relations and the complete-inventory audit remain open.
