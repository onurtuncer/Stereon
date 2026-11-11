# ADR-0053: Contract-driven metamorphic testing and intersection-bug promotion

- **Status:** Proposed
- **Date:** 2026-11-11
- **Owner:** Onur Tuncer
- **Phase:** 2 (v0.2 intersections), extending to Boolean, topology and machining queries
- **Related:** ADR-0001, ADR-0002, ADR-0007, ADR-0018, ADR-0026, ADR-0028, ADR-0029, ADR-0032, ADR-0041, ADR-0049, ADR-0050, ADR-0052

## Context

A CAD operation can produce an apparently plausible answer even when its geometry,
classification or topology is wrong. An independent ground-truth answer is often
expensive or unavailable for near-tangent and coincident intersections, trimmed
surfaces, spline intersections, Boolean reconstruction and machining queries.

A *metamorphic relation* (MR) specifies how the result must behave under a
transformation with known mathematical meaning. For operation `F`, input
transformation `T`, and appropriate output transformation `R_T`, a typical
relation is `F(T(x)) ~= R_T(F(x))`; discrete and topological outputs need
separate semantic comparison. Such a relation is a second oracle, not proof
that either run was absolutely correct.

Stereon already exercises some of these ideas in
`tests/intersect/AdversarialTests.cpp` (including exact translations,
powers-of-two scaling and parameter/origin changes). They must become an
explicit, contract-driven, auditable verification policy. In particular, a
confirmed intersection defect must never be fixed with only a point-example
test when a valid metamorphic relation can expose its failure mode.

## Options

1. **Handpicked regression fixtures only.** Cheap but misses equivalent
   configurations and correlated assumptions.
2. **Unconstrained randomized property tests.** Wider inputs but no
   necessary relationship to a geometric contract; difficult to reproduce.
3. **Contract-driven metamorphic tests with mandatory promotion of
   intersection defects.** More engineering effort, but reproducible and
   sensitive to coordinate, tolerance and classification faults.

## Decision

**Adopt option 3.** Stereon shall use reproducible metamorphic relations as a
first-class layer of validation, without replacing analytic oracles,
differential checks, independent structural validators, fuzzing or mutation
testing (ADR-0032/0041).

Each handler shall declare its applicable relations, domain restrictions,
expected semantic mapping, comparison rule and evidence level. A relation
cannot be claimed merely because a transformed invocation returned success.

### 1. Initial relation catalogue

| ID | Relation | What must agree | Preconditions and caveats |
| --- | --- | --- | --- |
| MR-01 | Rigid-motion equivariance | Transform intersection positions and geometric branches; preserve invariant classification and incidence | Initially proper rotations plus translations; quantify floating-point error from large offsets. Reflections are a separate orientation-aware case. |
| MR-02 | Operand permutation | Geometric result sets and appropriate contact semantics | Swap operand-specific parameters, normals, sidedness, provenance and curve orientation; not byte equality. |
| MR-03 | Positive uniform-scale covariance | Geometric result after scaling, dimensional policy and justified error bounds | Scale all relevant length quantities and context; never silently scale angles, dimensionless values or parametric thresholds. |
| MR-04 | Equivalent reparameterization | Geometric sets, appropriately remapped parameters, signs and derivatives | Supported charts and representations only; use the chain rule and handle singular/periodic charts explicitly. |
| MR-05 | Boolean set identities | Set-theoretic region, not identical B-rep cell decomposition | Supported valid solids, defined regularized-set semantics and resolved numerical classifications. |
| MR-06 | Topological invariance | Applicable connectivity, incidence, orientability, manifoldness, Euler characteristic and genus | Homeomorphic transformations; compare only properties defined for the model category. Euler characteristic alone is insufficient. |
| MR-07 | Composition / inverse transform | Equivalence after forward and inverse transformations | Numerical round-trip error accounted for; no bitwise floating-point assumption. |
| MR-08 | Partition and reconstruction | Covered region and consistent boundary incidence after imprint/split/reconstruction | Valid construction and explicit overlap/gap policy. |
| MR-09 | Machining query equivariance | Rigid-motion-invariant collision class and clearance, transformed contact and trajectory information | Supported tool/workpiece and bounded-motion contracts per ADR-0049. |

An intersecting pair of unbounded primitives may admit MR-01/02/03 but no
meaningful MR-04; conversely, a periodic curve can admit MR-04 without a
numerically stable equivalence near a seam. **Applicability is a testable
precondition, not a mechanism for excluding difficult supported inputs.**

For line/curve, curve/surface and surface/surface intersection, the initial
mandatory applicability review covers MR-01, MR-02, MR-03 and MR-04.
For a particular operand pair, MR-02 applies only if reversing the operands is
an exposed and semantically compatible operation.

### 2. Numerical and semantic comparison

No new global epsilon is introduced. The authority for model meaning,
pair tolerance and growth is ADR-0050 where it supersedes ADR-0001; distinct
numeric purposes and dimensions remain governed by ADR-0028.

- Preserve the distinction between exact predicate decisions, model
  coincidence, numerical uncertainty and healing.
- Transform **policy** separately from **geometry**; for MR-03 run both
  a scale-covariant test (dimensional context transformed) and, where useful,
  a fixed-context stress test (different question, separately labeled).
- Budget input transformation error, base-operation error, transformed-run
  error and comparison arithmetic. A computed residual alone is not a proof
  of geometric accuracy; never widen an acceptance threshold until it passes.
- Compare mathematical sets and typed classifications, not output ordering,
  entity IDs, or raw serialization. Map parameter domains, normals,
  transitions, orientation, provenance and branch traversal direction as
  required by the operation contract.
- Classifications are expected to agree only when the declared relation and
  stability conditions entail agreement. Distinguish `valid`,
  `violated`, `unresolved`, `inapplicable` and `unsupported`.
  An unresolved result cannot count as a pass.
- Independently check primal geometry on both sides (for example, every
  reported point lies on both supports and in the requested domains; a
  branch's certificate has the promised coverage). An entirely empty result
  must not pass only because the transformed result is also empty.
- For topological outputs, run an independent traversal/validator rather
  than trusting producer-maintained counts or success flags.
- Degenerate or discontinuous sensitivities follow ADR-0052's explicit
  `Nondifferentiable` outcomes; do not manufacture derivatives to satisfy
  an otherwise inapplicable transformation test.

### 3. Mandatory promotion of *every confirmed intersection defect*

**Every confirmed intersection bug, past or future, SHALL be inventoried and
promoted into a permanent metamorphic regression obligation.** This is a
merge and release requirement, not a recommendation.

For each defect, the repair PR shall carry:

1. A unique bug/PR identifier and the smallest reproducing input, preserved
   in a deterministic test or versioned corpus with replay information.
2. A previously failing *direct* regression proving the actual defect and
   its contract (not simply a crash/no-crash assertion).
3. At least one applicable **metamorphic regression**, chosen to exercise
   that same failure mode and neighboring configurations, naming the MR ID.
   It must test base and transformed results with semantic or independent
   geometric checks, not just compare two potentially wrong outputs.
4. Seeds/transformation parameters, complete effective Context, fixture
   identifier, and sufficient provenance to replay the failure.
5. Evidence that the case fails on the defective revision or a controlled
   reintroduction of the original fault, and passes on the repaired revision.

If no nontrivial mathematically sound relation applies to the defect's
supported domain, **do not invent one**. Record a *metamorphic exception*
in the defect ledger, identifying each considered MR and why it is
inapplicable, plus a stronger independent oracle or invariant-based
regression. Such an exception must be reviewed explicitly and may not be
used to omit the direct regression. A mathematically unsupported relation
is less valuable than an honest, independent oracle.

For previously fixed intersection bugs, construct a *backfill ledger* by
reviewing closed bug issues, merged repair PRs, existing regression test
comments and failing corpus entries. Do not mark the migration complete
until each confirmed defect is linked to a passing metamorphic test or a
reviewed exception. Duplicate symptoms of one root cause may share a
parameterized test, but each defect retains a traceable ledger entry.

Suggested columns for the ledger:

`Defect ID | Root cause | Handler | Minimal fixture | Direct regression |
MR ID/transformation | Oracle/semantic checks | Replay | Backfill state | Exception rationale`

Permitted backfill states are `untriaged`, `fixture captured`,
`MR implemented`, `exception reviewed`, and `verified`; only the
last two can close the obligation, with a reviewed exception requiring
the independent regression to be passing.

Do not rewrite an old bug test only to rename it "metamorphic": the second
execution must involve a meaningful transformation and an asserted,
mathematically justified relationship.

### 4. Execution architecture and evidence

Use the existing Catch2/RapidCheck and CMake/CTest infrastructure; extend
`tests/intersect/` and the existing adversarial/property suites instead
of introducing a second framework solely for this ADR.

A test case consists of:

1. A valid baseline fixture and independently checked baseline result.
2. A named, deterministic or seeded transformation with preconditions.
3. The transformed invocation under a recorded effective Context.
4. A relation-specific semantic comparison and transformed-result checks.
5. A replay artifact for any violation.

Failure reports shall record the repository revision, handler, MR, bug ID
when relevant, input and transformation parameters, seed/generator version,
numerical policy and units, compiler/build preset, and expected/actual
outcome. Failure reduction must preserve input validity and applicability.

Test families:

- **PR:** bounded deterministic fixtures, regression cases and small
  seeded property samples; mandatory for changed supported handlers.
- **Scheduled:** broader random and adversarial configurations, scale
  ranges, reparameterizations and degeneracies, with replay artifacts.
- **Release/research:** reproducible campaigns with independent reference
  solutions, defined datasets and comparative baselines where applicable.

The release dashboard shall report handler/MR coverage, applicable cases,
verified backfill obligations, unresolved outcomes, exceptions, failed
cases, and CI evidence. Raw pass percentage without exposure counts and
applicability data is not a robustness guarantee.

### 5. Research and performance discipline

Metamorphic results may support Stereon's reproducible CAD research only
when the relation definitions, baselines, uncertainty policies, hardware,
test seeds and denominators are published with the results.

A transformed pair with matching wrong answers is not evidence of absolute
correctness. Report independently established correctness separately from
metamorphic consistency. Keep correctness campaigns separate from
uninstrumented timing benchmarks (ADR-0043).

## Consequences

- Bug fixes now carry an explicit second-oracle obligation and traceable
  failure replay; CI and review cost increase.
- Existing adversarial tests become a foundation for standardized
  relations, not a reason to duplicate coverage or to claim an inventory
  that has not been performed.
- Some mathematically singular or tolerance-boundary cases need explicit
  applicability analysis and independent-oracle exceptions.
- The bug-to-MR ledger becomes a release artifact and can support
  numerical robustness publications without overstating proof.

## Verification

This ADR is implemented for v0.2 intersections only when:

1. All applicable v0.2 handlers have an explicit MR-01/02/03/04 matrix.
2. Representative direct and metamorphic regression tests exercise point,
   overlap, tangent, separated and unresolved result categories.
3. A recorded defect's original failure is detected and its transformed
   regression independently checks the returned geometry/classification.
4. Test cases reproduce from recorded seeds, contexts and fixtures.
5. Classification and numeric comparators obey ADR-0028/0029/0050.
6. Intersection PR CI executes the required tests, with no silent skips.
7. The backfill ledger enumerates every *identified* historical confirmed
   intersection defect; unreviewed historical PRs remain visible as
   `untriaged`, not silently considered complete.
8. Release checks refuse to treat open mandatory defect obligations or
   unresolved required outcomes as passes.

MR-05/06/08 follow the Boolean and topology milestones; MR-09 follows the
machining-query milestone. MR-07 is reusable from the first implementation.

**The ADR defines the policy. Its status remains Proposed until the
verification implementation and historical backfill evidence are reviewed.**
