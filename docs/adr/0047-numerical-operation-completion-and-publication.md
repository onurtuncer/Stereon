# ADR-0047: Numerical operation completion and publication

- **Status:** Proposed
- **Date:** 2026-10-28
- **Owner:** Onur Tuncer
- **Phase:** 0–1, applying to later numerical operations
- **Related:** [ADR-0002](0002-numeric-robustness.md), [ADR-0006](0006-error-handling.md), [ADR-0007](0007-threading-and-determinism.md), [ADR-0016](0016-cpp26-adoption-and-contracts.md), [ADR-0026](0026-modelling-operations-by-contracts.md), [ADR-0028](0028-purpose-specific-numerical-tolerances.md), [ADR-0029](0029-intersection-and-classification-outcomes.md), [ADR-0040](0040-atomic-topology-edit-publication.md), [ADR-0041](0041-layered-validation-evidence.md), [ADR-0042](0042-cpp26-contract-validation-responsibilities.md)

## Context

Valid input does not guarantee a valid or accurate computed result. Internal construction can bypass factory checks: parameter translation can collapse a knot interval, weight scaling can leave the supported range, and periodic extension can overflow. Conversely, degree elevation can return structurally valid data representing the wrong curve when ill-conditioned inverse steps amplify rounding.

A successful return must therefore mean more than completion of a loop or construction of an object. Structural validity, geometric fidelity and completion of the requested task are distinct obligations. A mathematically exact algorithm is not automatically stable in floating point; a small residual is not automatically a small forward error; sampling is not a global certificate.

The motivating findings were reproduced in an isolated source harness at revision `052e76f`. They motivate this policy, not a claim about native verification or the disposition of every finding at subsequent revisions. This draft is based on the repository at `1776a52`.

## Options

1. **Trust valid inputs and the mathematical algorithm.** Minimal overhead, but permits invalid-success and numerically corrupted-success results.
2. **Run an expensive universal geometric validator after every calculation.** Unclear guarantees, excessive cost in evaluators, and still insufficient for operation-specific fidelity or completeness.
3. **Require operation-specific completion evidence before publication.** Reuse structural checks, establish the promised numerical properties, and fail explicitly when the required evidence cannot be obtained.

## Decision

Adopt option 3. Extend ADR-0006's result boundary and ADR-0029's completion semantics to numerical construction, transformation, fitting, integration, projection and representation-changing operations. Specialized intersection/topology assurance requirements remain binding.

### 1. Declare the success contract

Each operation documents:

- supported input domain, including degree/order, parameter ranges, weights, degeneracies and representation limits;
- the requested invariant or approximation, its metric, units, domain and numerical policy;
- whether accuracy evidence is a rigorous bound, a justified floating-point error bound, an estimator, or a sample-based measurement;
- completeness requirements and finite iteration/subdivision/work budgets, when relevant;
- permitted fallbacks and whether they preserve parameterisation, structure, topology and accuracy;
- what success establishes, what failure means and what limitations remain.

Do not silently upgrade an estimate to a certificate or downgrade a certified request to a sampled approximation. Existing explicitly limited APIs may continue to promise limited assurance: this ADR does not declare all projection leaves certified or fitting deviations global. Such limitations must remain visible in documentation, result semantics and release evidence. A stronger consuming operation must obtain the additional evidence it needs.

### 2. Establish three distinct obligations

| Obligation | Required check or argument | What it does not prove |
| --- | --- | --- |
| Representation validity | Counts, finite required fields, ordering, multiplicities, nonempty required domains, supported weights, index ranges and type-specific invariants | That the result represents the intended geometry |
| Numerical/geometric fidelity | Preservation or approximation of the specified quantity over the declared domain, at the promised assurance level | That all requested solutions or components were found |
| Completion | Required domain/tasks processed; convergence and accuracy obligations met within budgets; unresolved work accounted for | Stronger assurance than the API promises |

A successful result must meet all applicable obligations. A query that does not construct geometry need not repeat irrelevant construction checks. A deliberate empty result is permitted only when its operation-specific semantics establish emptiness; an empty work list caused by failure is not evidence.

### 3. Validate candidates before publication

Construct into local candidate state and publish only after the operation's obligations hold. Reuse the invariant implementation shared with checked factories; do not maintain weaker duplicate checks for generated data. Avoid an unnecessary public-factory round trip or full copy when an internal validator can check the candidate directly.

Internal unchecked aggregate construction is permitted only where the function establishes the invariants and makes that obligation explicit. Any optimized path that skips an invariant check must document why its construction preserves the invariant and have tests covering that argument. Data-dependent arithmetic must not inherit an unchecked-path exemption merely because its inputs passed a factory.

Check results of arithmetic that can invalidate a representation: finite differences and total knot width, representability and ordering after translation, positive supported weights after scaling, and checked size/index arithmetic before allocation or access. Not every type requires every value finite: intentionally unbounded domains retain their declared representation rules.

These checks belong at operation/publication boundaries. Do not add full knot-vector scans, heap allocation or model-wide validation to scalar evaluator hot loops. Low-level evaluators retain their documented preconditions and numerical behavior under ADR-0016; a fallible operation consuming them must detect outcomes inconsistent with its own success contract.

### 4. Account for conditioning and roundoff

Do not infer numerical fidelity from real-arithmetic identities alone. Inversion, division by small differences, cancellation, repeated reconstruction and transformations of poorly scaled parameters require an operation-specific stability argument or a posteriori evidence appropriate to the promised guarantee.

A reconstruction residual bounds geometric error only through a justified relationship. Include conditioning or an equivalent bound when needed; do not introduce one arbitrary global condition-number cutoff. Rescaling or normalization is allowed only if it preserves the declared semantics and supported representation.

For degree elevation, matching degree/counts/knots is insufficient. Establish shape preservation to its documented floating-point allowance. If inverse knot removal cannot establish it, use a stable alternative or fail. Keeping elevated Bézier pieces with extra knots is a fallback only if the API explicitly permits a nonminimal representation; it must not silently break a minimal-pole promise.

For concatenation, normalization must preserve the curve and supported weight range. If shifting the second parameter domain collapses a span, report inability to construct the requested representation; do not silently change parameterisation. For periodic conversion, validate both seam semantics and the extended representation; an approximate periodic-recognition shortcut cannot bypass the requested closure contract.

### 5. Consume numerical evidence explicitly

Residuals, convergence flags, error estimates, bounds and removal counts are part of the operation's decision. Do not discard them when they determine whether a result meets its contract. Numerical helpers should return named result structures with `[[nodiscard]]` where ignoring the result would lose necessary evidence.

An intentional discard requires a local explanation and a verified alternative argument. Casting a result to void is not evidence. Prefer separate helper interfaces when one caller genuinely needs a different guarantee rather than repeatedly bypassing checks.

An iteration limit, stalled step, precision limit or exhausted work budget is not convergence. A zero or non-finite residual must also be interpreted under the method's assumptions; underflow and overflow do not constitute success evidence.

### 6. Preserve failure distinctions

Use the existing error vocabulary; no enum renumbering or new code is introduced by this documentation change.

| Situation | Meaning and usual code |
| --- | --- |
| Invalid caller data | `InvalidInput`, or `ParseError` at a parser boundary |
| Valid input outside the explicitly supported representation/algorithm domain | `NotSupported`; `NotImplemented` only for a planned, absent implementation |
| Supported computation cannot settle its required conclusion within its work or precision limits | `Unresolved` |
| Required accuracy certificate cannot be established | `CertificationFailed` |
| Candidate violates representation/post-validation requirements | `InvalidResult`; this remains a defect to investigate, not normal success |
| Operation-specific tolerance-growth ceiling exceeded | `ToleranceExceeded` under the existing tolerance policy |
| Degeneracy with an existing defined failure contract | `DegenerateResult`; do not label mere numerical uncertainty as proven degeneracy |

These are semantic choices, not interchangeable fallback labels. A specialized API may retain an established failure code, but its diagnostic must distinguish the reason. Unsupported representability must be declared; it cannot be invented after a faulty implementation fails an otherwise supported case.

Diagnostics identify the operation/stage, relevant domain or entity, effective tolerance/budget, observed residual or bound and termination reason where applicable. Do not put dimensionless quantities into millimetre-labelled tolerance fields. Use existing diagnostic messages until typed fields are implemented; future structured payloads require an implementation/schema review.

No failed operation publishes a partial ordinary result. Exploratory partial-result APIs, if added, require a distinct contract and cannot feed complete-result consumers without additional validation. `std::bad_alloc` may propagate under ADR-0006, but existing inputs and published snapshots must remain unchanged. Publication and associated metadata/history remain atomic under ADR-0040.

### 7. Test the advertised domain and algebraic relationships

Each numerical operation needs focused examples plus properties drawn from its full declared domain. Generators must deliberately cover supported maximum degree/order, rational/nonrational cases, repeated and strongly uneven knots, translated and rescaled parameter domains, permitted weight limits and near-degenerate configurations. Ordinary random sampling alone is not a boundary strategy. Reject invalid generated inputs explicitly so a campaign does not silently test only an easy subset.

Use operation-specific metamorphic tests, including where applicable:

- evaluation before/after degree elevation or knot insertion, with parameter correspondence preserved;
- knot-removal deviation against its reported bound;
- reversal twice and reparameterisation with mapped parameters;
- uniform weight scaling within supported limits;
- translation/scaling covariance and periodic seam consistency;
- archive round trips of successfully generated geometry, where that geometry is in the stored subset.

Expected numerical differences must use justified tolerances tied to scale and the operation's guarantee; not every transformation is bit-identical. Preserve exact bit comparisons where an API or ADR-0007 explicitly promises them. Apply transformations only when the transformed input remains representable and supported.

Use independent analytic, higher-precision or comparator oracles where suitable. Dense samples help expose errors but do not certify unsampled behavior. Never loosen a test tolerance merely to absorb a discovered blow-up. Keep minimized failures as permanent regressions, with the numerical policy needed to replay them.

### 8. Separate optimization from changes in numerical semantics

Allocation removal, SIMD and workspace reuse must preserve the established completion/failure contract and required arithmetic ordering. Compare optimized and reference paths on boundary cases, not only representative throughput fixtures. A performance gain cannot justify dropped evidence, a weaker guarantee or successful publication of an invalid candidate.

Prefer separate reviewable changes for algorithm repair and performance refactoring. If they must be coupled, state why and provide evidence that isolates both effects. Measure operation-level performance and allocations after correctness checks are included.

## Migration and release evidence

| Work item | Required disposition |
| --- | --- |
| B-spline rearrangements | Common candidate-invariant validation; elevation fidelity evidence; periodic and concatenation representability checks |
| Iterative geometry algorithms | Inventory termination conditions, ignored evidence and budget-exhaustion behavior |
| Fitting, projection and integration | Reconcile estimator/sample/certificate terminology with actual guarantees; preserve declared limitations |
| Property tests | Expand domain coverage beyond low-degree, moderately scaled fixtures; retain minimized audit cases |
| Publication/failure paths | Verify unchanged inputs on numerical failure and injected allocation failure |
| Release documents | Separate implemented, verified and accepted; link exact commit, corpus, native-toolchain and comparator results |

This inventory is not a statement that all rows are currently missing or all existing algorithms are defective. Record evidence and close each obligation individually. Optional capabilities may be explicitly deferred or withdrawn from release exposure; advertising them as supported requires their success contract to hold.

## Consequences

Success becomes a meaningful, operation-specific claim. Some formerly successful unstable cases will return explicit failure until a sound algorithm supports them. Candidate validation and fidelity checks may add cost; benchmark that cost at the boundary instead of weakening the guarantee.

The policy avoids a universal expensive validator and does not require a generic certificate object on every small result. New public diagnostics or changes in guarantees follow the compatibility and ADR review process. Documentation alone does not complete implementation or release qualification.

## Verification

| Criterion | Pass evidence |
| --- | --- |
| AC-0047-01: Declared contracts | Reviewed inventory of supported domains, metrics, assurance, work limits and fallback semantics for targeted operations |
| AC-0047-02: Candidate invariants | Native regressions show collapsed spans, overflowing periodic extensions and unsupported generated weights cannot be published as valid results |
| AC-0047-03: Fidelity | Elevation and related operations pass independent/property checks across their full supported range, including the uneven-knot audit fixtures; residual use has a justified error relationship |
| AC-0047-04: Completion | Forced iteration/subdivision limits and stalled/precision-limited cases produce the documented failures, never false convergence or false emptiness |
| AC-0047-05: Evidence consumption | Review finds no unexplained discard of success-critical bounds, residuals or convergence flags; regression tests exercise their rejection paths |
| AC-0047-06: Metamorphic coverage | Deterministic campaigns exercise the listed relationships and boundary families, with explicit valid-input filtering and justified tolerances |
| AC-0047-07: Failure atomicity | Numerical failures and injected allocation failures leave inputs and published snapshots unchanged; no partial ordinary result escapes |
| AC-0047-08: Integration qualification | Native GCC 16 tests and applicable sanitizers, archive/comparator checks and determinism runs pass; performance reports include required checks and exact revision evidence |

**Acceptance:** Proposed until owner review, migration disposition and linked verification evidence. This ADR extends rather than replaces the specialized accepted contracts in ADR-0006 and ADR-0029. It does not itself mark v0.1 ready or assert that the motivating defects are fixed.
