# Numerical audit replay, 2026-11-12

Base: main `3c0ee49c3d09b09897c23e7a1b6ce0de556f7174`.
Native validation uses GCC 16.2, the gcc16-debug preset and its default floating-point flags.

## Cone polishing

CI seed `17298607298449989786` exposed origin-shift sensitivity near a cone generator.
Separate multiplication and addition in the polishing point lose the small offset
when the distant origin cancels the line parameter. Evaluate those coordinates
with `std::fma` before forming the correction; no tolerance is increased.

`MR audit cone crossings retain their positions after an exact origin shift`
in LineRevolutionTests.cpp uses dyadic geometry, an exact parameter-origin shift
of 2^19 and a translation (8,-4,2). Its independent ideal cone roots are
-191 +/- 60 sqrt(10); actual rounded-cone support residuals are also checked.
MR-04 compares shifted parameters and positions, and MR-01 repeats the case
under translation. The 40 assertions and original CI seed pass after repair.

## Test oracle corrections

Seed `8908293537530491705` compares exact ellipsoid separation with a sphere's
tolerance contact. Cardinalities need not agree in that boundary band. Check the
sphere's contact against both supports and the boundary distance instead. Keep
ordinary branch geometry checks. A deterministic exact-separation/tolerant-contact
fixture records this distinction. The seed replay passes 218 assertions in five
matching tests.

Seeds `17041817369036985825` and `16155406228833616885` encounter an unresolved
ellipse domain endpoint and cylinder tangency before the derivative property's
conditioning precondition. Only IntersectionFailed may be discarded by that
precondition; answered cases still undergo all derivative and transpose checks.
Both seeds pass. Deterministic fixtures retain both unresolved outcomes and a
regular shifted-chart control requiring two crossings. These are outstanding
kernel capability limitations, not successful intersection qualification.
Replaying the #283 boundary change alone did not resolve either seed.

Replay using RC_PARAMS=seed=<seed> and the test names above; the derivative test
is `circle and ellipse crossings: second-order Taylor decay and the transpose`.
The three new deterministic fixtures pass 55 assertions. This evidence is
implemented, not independently reviewed release qualification.

## Follow-up: coverage branch certification failure (2026-11-13)

The initial #286 head passed Linux, Windows and macOS workflows, but coverage
failed one of 1071 CTest cases (four unrelated tool tests were skipped).
[Coverage job](https://github.com/onurtuncer/Stereon/actions/runs/43300991832/job/131471814362)
failed `branch points: second-order Taylor decay and the transpose` with
RapidCheck seed `13593231413355092879` (Catch seed `3487449154`).

Both sphere properties shrink to the same generator integers:
`[-40, -40, -40, 8, -20, -39, -24, -5]`. In RandomPlaneSphere these mean
center (-5,-5,-5), radius 1, normal (-2.5,-4.875,-3), and signed plane offset
-5/8. Construct the plane with PlaneWithNormal(center + offset *
(normal / Length(normal)), normal), exactly as the existing generator does.

The failure occurs in Meet at BranchLinearizationTests.cpp:151, before either
derivative identity is evaluated: "a traced branch exceeds the sample budget
before it certifies within the resolution 1e-05 mm". This is an unresolved
primal certification failure, not evidence that a derivative is incorrect.
Its relationship to this change has not yet been established by a base/head
replay. Do not suppress the assertion or increase tolerances to make CI green.

Replay with the gcc16-coverage preset and:
```sh
RC_PARAMS=seed=13593231413355092879 build/gcc16-coverage/tests/stereon_intersect_tests \
  'branch points: second-order Taylor decay and the transpose'
```
Next: reproduce on base and head, isolate the pcurve certification budget
failure, then add direct and metamorphic regression coverage for the repair.


## Narrow seam excursion repair (2026-11-13)

The shrunk sphere-plane case was reproduced directly. Its circle crosses the
sphere longitude seam twice between adjacent points of the initial 64-step
scan. Chart then clamps the missed excursion to the seam: increasing fitting
samples leaves a deviation plateau near 0.0002500838 mm, above the unchanged
0.00001 mm resolution. This is not a derivative failure.

FitCertified now inspects its unwrapped fitting samples for additional seam
crossings before fitting/clamping. It returns bracketed, bisected cuts to
AddBranch, which retries the resulting pieces independently. The initial scan,
resolution, sample budget and mandatory deviation certificates are unchanged.
The number of inserted cuts is bounded; exhaustion fails closed. Endpoint
rounding slivers use the existing geometric reach criterion. Sampling discovers
cuts but does not replace certification or prove arbitrary seam completeness.

The regression `MR audit narrow sphere seam excursions remain complete certified circles`
in ImplicitTraceTests.cpp covers 24 configurations: translations, exchange,
dyadic scaling and a quarter-turn sphere chart (MR-01 through MR-04). Its
independent plane-section circle oracle checks radial/plane residuals, total
angular coverage and connectivity; both surface pcurves must certify.

Observed GCC 16.2 debug results before the workspace reset:
- Regression: 19,080 assertions passed.
- Entire [trace] group: 20,888 assertions in 17 tests passed.
- Restoring the original ImplicitTrace.cpp makes the final regression fail at
  its first required successful result with the original sample-budget error.
- The full RapidCheck seed replay was started, but its completion was not
  observed before the workspace reset. It is not counted as passing evidence.

The same source/test edits were reconstructed after the reset. Formatting,
complexity, inventory and literal checks are rerun; CI must validate the
published commit, including the coverage configuration. Qualification remains
MR implemented, not reviewed/verified release evidence.
