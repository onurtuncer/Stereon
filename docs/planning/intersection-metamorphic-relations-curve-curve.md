# Curve–curve relation cells: MR-01 to MR-04 for every ordered pair of 3D curves

- **Date:** 2026-11-12
- **Branch:** `verify/relations-curve-curve`, stacked on the planar relations branch (#274)
- **Policy:** [ADR-0053](../adr/0053-metamorphic-testing.md), [inventory](intersection-metamorphic.json)
- **Toolchain:** GCC 16.2 (MSYS2 UCRT64), preset `gcc16-debug`, Windows 11
- **Reviewer:** the v0.2 verification agent (Agent E). "Verified" means this documented review, pending the owner's review, as for the [planar cells](intersection-metamorphic-relations-planar.md).

The 36 ordered pairs of `Line`, `Circle`, `Ellipse`, `BSplineCurve`, `Hyperbola` and `Parabola` each have four cells, 144 in all. The oracle is the existing 3D corpus [`corpus/curve-curve/`](../../corpus/curve-curve/README.md). `tools/generate_curve_corpus.py` writes 168 cases, 8 for each of the 21 unordered families. It translates the second curve so that its point at a stated parameter lies a stated gap from the first curve's point at its stated parameter, along their common normal. Those parameters are then a critical point of the distance at exactly that gap. The claim of a case is local: the curves meet near the stated parameters exactly when the gap is within the resolution. The gaps are 0, 0.4, 1.5 and 10 resolutions, so every family has meeting and non-meeting cases. The kernel is never consulted for the answer.

#### Domain

The domain is the corpus's: the 21 families, each with local claims around one constructed closest approach at gaps of 0, 0.4, 1.5 and 10 resolutions, with the cubic B-splines polynomial. The claim says nothing about meetings elsewhere on the pair, and the relations inherit that locality. Out of scope: rational B-splines in 3D (not in this corpus), global completeness, and user types.

#### Relations as the tests run them

- **MR-02, operand exchange.** This is the existing test "every curve-curve corpus case meets exactly when its gap is within the resolution" (`tests/intersect/CurveCorpusTests.cpp`). It runs every case in both orders and checks the claim with the parameters exchanged.
- **MR-01, rigid motion.** New: "MR-01: curve-curve corpus claims hold under exact rigid motions" (`tests/intersect/CurveCorpusRelationTests.cpp`). Every case is run under the two cyclic axis permutations and a quarter turn about z, all exact rotations, each with and without the translation (8, -4, 2). The translation alone is also run, so there are seven motions per case, in both orders. Points move, unit frames turn, and parameters keep their meaning. Translation rounding (about 1e-15) is far inside the claim's margin (at least 0.4 resolution, about 4e-6).
- **MR-03, uniform scale.** New: "MR-03: curve-curve corpus claims hold under power-of-two scales". Scales 1/4 and 4, with the Context's resolution and tolerance ceiling scaled with them: ADR-0053 transforms the policy separately from the geometry. A parabola's parameter is a length and scales; the others do not.
- **MR-04, equivalent reparameterization.** New: "MR-04: curve-curve corpus claims hold in equivalent charts". Each operand in turn is moved to an equivalent chart, built from the case data:
  - a line's origin moved two directions on with its direction doubled, t' = (t - 2) / 2;
  - a circle's frame turned a quarter, t' = t - pi/2;
  - an ellipse's frame turned half, t' = t - pi;
  - a hyperbola's or a parabola's y axis negated, t' = -t;
  - a B-spline reversed (knots reflected, poles reversed), t' = a + b - t.

  The domains are mapped in the same way. Each case's claim must hold at the closed-form mapped parameters, in both orders.

The 4 BSplineCurve:BSplineCurve cells previously marked `MR implemented` by the PR-241 campaign test are qualified here by the corpus relations. The campaign test still runs, and its evidence is unchanged.

#### Replay

```sh
cmake --build --preset gcc16-debug --target stereon_intersect_tests
build/gcc16-debug/tests/stereon_intersect_tests "every curve-curve corpus case meets exactly when its gap is within the resolution"
build/gcc16-debug/tests/stereon_intersect_tests "MR-01: curve-curve*"
build/gcc16-debug/tests/stereon_intersect_tests "MR-03: curve-curve*"
build/gcc16-debug/tests/stereon_intersect_tests "MR-04: curve-curve*"
```

The fixtures are the committed corpus files; there is no seed. The tests run in every pull request (not hidden).

#### Green

On this branch (`main` 227d79a with test-only changes), all four tests pass:

| Test | Assertions |
| --- | --- |
| MR-02 (existing) | 673 |
| MR-01 | 4705 |
| MR-03 | 1345 |
| MR-04 | 1345 |

Every run returned a result (the tests `REQUIRE` it), and every claim held.

#### Spot-checks: faults the cited tests kill

Each fault was applied alone to a working tree with no committed change to `src/`, rebuilt, and run twice; the source was then restored.

- **MR-02, operand order.** The fault: `Intersect` routes an ordered pair whose first type comes later in the catalogue to the mirrored handler without exchanging the parameters. The existing both-orders corpus test kills it in both runs. It fails in the families whose handler is not symmetric (circle, ellipse, hyperbola, parabola and B-spline pairs), the only ones where the fault changes anything.
- **MR-03, scale.** The fault: the curve pair search uses the default resolution instead of the Context's. Only the new MR-03 test kills it, in both runs, in every family the pair search handles. The existing corpus test, MR-01 and MR-04 all pass with it, because they run at scale 1. The relation catches a fault the plain corpus misses.

No MR-04-specific fault site was found. `LineLine` is exact and the conic pieces read the frames as given, so no sampled fault there would affect only a changed chart. MR-04's evidence is its green run and the kills above, which it shares: it runs both orders, so it also kills the MR-02 fault.
