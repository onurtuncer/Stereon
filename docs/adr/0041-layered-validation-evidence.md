# ADR-0041: Use layered validation with revision-bound evidence

- **Status:** Proposed
- **Date:** 2026-10-20
- **Owner:** Onur Tuncer
- **Phase:** 0
- **Related:** Validity checking; topology categories; topology/geometry separation; operation contracts; independent validation at Freeze.

## Context

One `is_valid` flag cannot distinguish intact references from a valid solid or a correct modelling result. A topologically sound face may have an invalid pcurve; a valid solid may still be the wrong Boolean result. Rechecking all geometry at every Euler edit is expensive, but reusing stale validation evidence is unsound.

ADR-0038 requires independent structural validation at Freeze, and accepted ADR-0027 requires the independent validator at Freeze and Commit. This ADR preserves that requirement and defines complementary evidence for geometry, model category and operation semantics.

## Decision

Preserve ADR-0018's cumulative Structural, Geometric and Solid check levels and its build/Context scheduling. Add property-specific evidence identifying checks, scope, numerical assurance and the input state. These evidence families supplement the existing levels; they do not replace them with a new incompatible API. Public success requires all properties promised by the selected operation contract, established by validation, reviewed construction or current scoped evidence as permitted by ADR-0026. Context-selected checker levels remain unchanged; disabling an optional check does not produce evidence that the check ran.

### Validation responsibilities

| Family | Examples | Does not establish |
| --- | --- | --- |
| Handle and structural | Liveness, ownership, reciprocal Next/Prev, radial membership, endpoint incidence, root reachability. | Geometric embedding or a material region. |
| Association and geometry | Parameter domains, endpoint agreement, pcurve/support agreement, trim consistency under the selected tolerance context. | Complete solid semantics or the requested operation's meaning. |
| Category | Manifold vertex links, required closure/orientation, connectedness, permitted exterior/cavity relationships. | That an operation produced the intended region. |
| Operation semantics | Partition/union/intersection relationship, frame conditions, lineage and attribute guarantees. | Unrelated properties not included in the contract. |

Topology and category properties can overlap; the report names actual properties rather than inferring all of one family from a label. A lower-dimensional or mixed model is not invalid merely because it is not a solid. Category assertions follow the operation's declared supported domain and ADR-0018. This proposal does not import unintegrated category/ownership drafts or change the existing entity model.

### Required scheduling

- **On entry:** check handle/parameter safety and validate or reuse the required input evidence.
- **After each primitive edit in debug and sanitizer builds:** perform the required local postchecks over the complete affected neighborhood, including all radial/coedge uses and dependent associations.
- **At Freeze and Commit:** run the independent structural validator required by ADR-0038 and ADR-0027. An incremental replacement for this full check requires a later explicit decision and demonstrated coverage of preserved invariants.
- **Before public operation publication:** establish the selected output category, geometric and operation-specific guarantees. Unaffected evidence may be reused only under the dependency rules below.
- **During testing and diagnostics:** run broader cross-checks and independent semantic oracles where applicable. Optional extra checks may change cost, not the meaning of success.

The independent structural validator reads authoritative model records rather than trusting only an operator's dirty list, cached counts or success flag. It bounds traversal by available records/visited sets and checks references before dereference, so corrupt cycles and handles produce diagnostics rather than unbounded loops or unsafe access.

### Evidence identity and invalidation

Evidence records at least model/snapshot or builder identity, revision, property/profile, entity/dependency scope, numerical context, method/assurance, and contract/validator version as applicable. It cannot be forged by ordinary API callers merely setting a flag.

A relevant change to topology, geometry, associations, ownership or tolerance context invalidates dependent evidence. Initial implementations may invalidate all evidence on any builder mutation. Selective reuse is an optimization requiring explicit dependency closure and tests. Evidence for an immutable snapshot remains valid for that snapshot, not automatically for a derived candidate.

Geometric evidence is not silently reinterpreted after a tolerance change. Revalidation or an explicitly reviewed transfer rule is required. Deserialized evidence is untrusted by default; rerun required checks unless a separate validated persistence scheme is adopted. Do not reuse evidence merely because entity IDs or storage addresses happen to match.

### Assurance and unresolved outcomes

Reports distinguish established, violated, unresolved and not checked for relevant properties. A timeout or numerically ambiguous result cannot become a pass. A required unresolved property blocks publication under that contract; explicit weaker result categories or approximation policies require their own contracts.

Sampling may provide diagnostic evidence but must not be labelled a complete guarantee over a continuous curve or surface without an applicable bound. A certificate here is an inspectable record of a stated method and assurance, not automatically an exact mathematical proof. Construction evidence and numerical validation must describe their limitations honestly.

### Efficient local checks

The affected neighborhood includes indirect dependencies: splitting an edge requires visiting every incident coedge, including repeated seam uses and supported non-manifold radial uses. Shared geometry changes require either isolation through immutable/copy-on-write supports or invalidation of every affected consumer.

Tests compare local checks with an independent full validator on generated and adversarial fixtures. They deliberately corrupt data omitted from a claimed dirty set to ensure publication-level validation detects that omission. Optimization shall not replace checks with the Euler equation alone.

## Options

- **One global valid bit:** rejected because scope, assurance and staleness are hidden.
- **Only local operator checks:** rejected as the publication gate; correlated implementation defects can escape.
- **Full geometric validation after every splice:** rejected as the universal schedule; valid construction evidence and unaffected-property reuse are permitted.
- **Accept unresolved as probably valid:** rejected unless an explicit different contract provides a weaker, accurately named result.

## Consequences

Callers can distinguish structural validity, model category and semantic correctness. Evidence enables measured optimization but adds metadata and dependency-management obligations. Conservative invalidation is the initial correctness baseline. The cost of independent structural Freeze validation remains visible and must be benchmarked rather than silently removed.

## Verification

- A structurally sound face with an invalid pcurve fails the relevant geometric property.
- A closed non-manifold complex does not receive a manifold-shell or Solid claim.
- A geometrically valid but incorrect split result fails the operation-semantic check or cannot establish the required evidence.
- Stale revision, changed geometry, changed tolerance, wrong model and deserialized untrusted evidence are rejected or revalidated.
- Corrupt loops, invalid handles and missing radial membership are detected with bounded validator execution.
- Unresolved and not-checked properties cannot satisfy mandatory publication requirements.
- Local/full-validator differential tests and targeted mutations detect incorrect splices, omitted incidence updates and missed invalidation.

**Acceptance owner:** independent validator maintainer with geometry and modelling-operation review.

## Related decisions

- [ADR-0006](0006-error-handling.md)
- [ADR-0016](0016-cpp26-adoption-and-contracts.md)
- [ADR-0018](0018-validity-checking.md)
- [ADR-0026](0026-modelling-operations-by-contracts.md)
- [ADR-0027](0027-controlled-topology-mutation-and-euler-operators.md)
- [ADR-0038](0038-topology-traversal-lifetime-and-validation.md)
- [ADR-0039](0039-executable-topology-contracts.md)
- [ADR-0040](0040-atomic-topology-edit-publication.md)
- [ADR-0042](0042-cpp26-contract-validation-responsibilities.md)
