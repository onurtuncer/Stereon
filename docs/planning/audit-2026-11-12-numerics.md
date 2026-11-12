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
