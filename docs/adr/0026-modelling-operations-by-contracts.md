# ADR-0026: Specify Modelling Operations by Contracts

- **Status:** Proposed
- **Date:** 2026-10-13
- **Owner:** Onur Tuncer
- **Phase:** 0
- **Scope:** Public kernel modelling operations and their observable behavior

## Context

A CAD kernel must define what an operation means, which inputs it accepts, what its output guarantees, and what callers can rely on when it fails. Function signatures and algorithm descriptions alone do not establish these properties.

An operation may produce structurally consistent topology while violating geometric constraints. Numerical uncertainty, unsupported degeneracies, or resource exhaustion may prevent an otherwise meaningful mathematical operation from completing. An empty geometric result may also be a legitimate success rather than a failure.

Without explicit contracts, callers must infer these distinctions from implementation behavior. This makes algorithm replacement, regression testing, application integration, and recovery unreliable.

Ian Stroud's *Boundary Representation Modelling Techniques*, particularly its treatment of stepwise modelling algorithms and modelling operator definition, provides conceptual motivation. The contract format and failure guarantees below are Stereon design decisions, not prescriptions attributed to the book.

## Options

### Document only signatures and nominal examples

This is inexpensive initially but leaves degenerate cases, failure recovery, and result guarantees undefined. Rejected for supported public operations.

### Permit best-effort partial mutation on failure

This can simplify implementations but transfers recovery responsibility to every caller and risks corrupting subsequent operations. Rejected as the default public-operation contract. Explicitly staged user workflows remain possible through separately committed operations.

### Use one universal validity guarantee for all operations

This obscures distinctions among solids, sheets, wires, approximation, and intermediate geometric results. Rejected. Each operation declares its result category and guarantees within shared terminology.

### Require exact arithmetic and complete geometric certification for every operation

This would impose an implementation strategy and assurance level beyond what this behavioral ADR establishes. Rejected as a universal requirement. Individual operations may provide stronger guarantees and must state them accurately.

## Decision

**Every public modelling operation shall have a documented, versioned behavioral contract specifying accepted inputs, result semantics, success guarantees, and failure behavior.**

**A failed operation shall leave its input models and previously published model state unchanged.** A successful operation may publish only a result satisfying the guarantees declared by its contract.

These requirements apply to operations that construct or change geometry, topology, or their associations, including primitives, transformations, splitting, sewing, Boolean operations, sweeping, and healing as they become available. Internal algorithms may use narrower contracts, but they must collectively uphold the public operation's guarantees.

This decision specifies behavioral contracts, which are distinct from C++ language contracts. It supplements ADR-0005 (immutable results), ADR-0006 (`std::expected` and `OperationResult`), ADR-0016 (programmer obligations), and ADR-0018 (validation levels); it does not reopen those choices.

## Required contract fields

Each public operation shall document the following fields before being treated as a supported API:

| Field | Required content |
| --- | --- |
| Purpose and semantics | The mathematical or geometric transformation, including the meaning of the result. |
| Input categories | Accepted solids, shells, sheets, wires, geometric objects, collections, and their required validity conditions. |
| Parameters and context | Parameter domains, units, orientation conventions, tolerance context, and relevant execution options. |
| Input validity and programmer obligations | Checked data requirements and their validation evidence, separately from low-level programmer obligations governed by ADR-0016. |
| Supported cases | Treatment of disconnected inputs, empty inputs, coincidence, tangency, degeneracy, and relevant representation restrictions. |
| Success postconditions | Output category, structural and geometric guarantees, and the relation between the result and the requested operation. |
| Mutation and identity | Immutable input/output and sharing behavior under ADR-0005, plus history guarantees under ADR-0011; shape-local handles do not provide cross-operation identity. |
| Failure outcomes | Distinguishable failure categories, required diagnostics, and the unchanged-state guarantee. |
| Numerical assurance | Tolerance policy used, validation performed, and the limits of any exactness or geometric-validity claim. |
| Reproducibility | ADR-0007 G1 applies; identify any additional guarantees across algorithm versions or outside the G1 envelope. |
| Sensitivities | Derivative capability, parameter dependencies, correspondence, regularity assumptions, and statuses under ADR-0021; unsupported products are explicit. |

A contract may explicitly state that a behavior is unspecified. It shall not silently imply guarantees that the implementation does not provide.

## Contract rules

### 1. Preconditions define the supported domain

Contracts shall distinguish invalid input from valid input that lies outside the implementation's supported domain.

Examples include an invalid reference, a non-finite parameter, an invalid solid boundary, a valid but unsupported non-manifold input, or a valid tangent configuration that a particular implementation cannot resolve.

Public operations validate external data and return `KernelError` for invalid data, as required by ADR-0006. A user-provided shape, numerical parameter, or imported geometry must never trigger a C++ contract violation merely because it is invalid.

Handle/store misuse, low-level evaluator domain obligations, and internal programming invariants remain governed by ADR-0016. Their contract violations are not converted into recoverable `KernelError` values. The word "precondition" in the behavioral specification must not blur this boundary.

Validation evidence may be reused only when it applies to the current immutable model and the required checks. This proposal creates no public unchecked-input escape hatch and does not remove any ADR-0018 mandatory checks.

### 2. Success has explicit meaning

An operation shall report success only when its declared postconditions are established by construction, validation, or a documented combination of both.

The contract shall identify which guarantees are structural, which are geometric, and which depend on numerical tolerances or approximation. A structural validator passing shall not be described as proof of complete solid validity.

Legitimate empty results shall be represented as successful empty results where the operation's semantics permit them. For example, the regularized intersection of disjoint solids is empty; it is not an algorithm failure.

Representational freedom across algorithm versions shall be explicit. Within ADR-0007's G1 envelope, the same executable, hardware, inputs, and Context settings must produce bit-identical results across thread counts and schedulers; decomposition and ordering are not free to vary there. Cross-version stable naming is not introduced by this ADR.

### 3. Failure is distinguishable from a valid result

The public API shall distinguish at least the following semantic outcomes when relevant:

- **Invalid input:** a required input or parameter condition is violated.
- **Unsupported case:** the input is meaningful but outside the documented implementation scope.
- **Numerical or convergence failure:** the implementation cannot establish the required result within its numerical or iteration limits.
- **Resource limit or cancellation:** execution ends before a result can be committed.
- **Invalid result detected by validation:** a post-check fails and produces `KernelError::InvalidResult` with the validity report, as specified by ADR-0018.

Internal programmer-invariant violations remain C++ contract violations under ADR-0016, not a recoverable failure category.

Use ADR-0006's `ErrorCode` and structured diagnostics for these outcomes; any additional cancellation/resource-limit codes are defined when the corresponding feature is introduced. Shape-producing operations return `std::expected<OperationResult, KernelError>`. `std::bad_alloc` may propagate under ADR-0006, while preserving immutable inputs. The API must not encode every failure as an empty shape, a null reference without explanation, or a nominally successful but unvalidated result.

Diagnostics shall identify the operation and failed stage and, where available, the relevant input entities and failed condition. Numerically unresolved cases shall be reported as unresolved rather than asserted to be geometrically impossible.

### 4. Public operations provide the strong failure guarantee

Before success is committed, candidate changes shall remain private to the operation or otherwise reversible. Validation required by the contract shall occur before publication.

On any failure delivered through the supported API, inputs and previously published model state shall remain observably unchanged. This includes topology, geometry, associations, externally visible identities, and model metadata. Shared objects and aliases shall not expose partial edits.

Implementation follows ADR-0004 and ADR-0005: a private mutable builder and copy-on-write pages produce a new immutable shape. There is no public in-place shape mutation or kernel-level undo log. The publication step must not expose a partial shape, history, or diagnostics if allocation or validation fails.

Diagnostic logs, execution statistics, allocation activity, and semantically transparent caches may change. Such caches must remain correct or be invalidated and must not change subsequent model semantics. The guarantee does not claim recovery from process termination, hardware failure, or other events outside the supported failure model.

If an operation receives invalid input, failure preserves that input unchanged; it does not promise to repair it.

### 5. Approximation and healing must be explicit

An operation shall not silently expand tolerances, perturb geometry, remove small features, or invoke healing outside its documented contract.

Approximate or healing behavior may be supported through an explicit operation or policy. Its contract shall describe the permitted changes, applicable bounds or acceptance criteria, and the information returned to the caller. The result must still satisfy the postconditions of that selected policy.

A numerical fallback that changes implementation technique without changing the promised semantics is permitted without redefining the operation. A fallback that weakens the result guarantee is not.

ADR-0001 defines the tolerance model and healing budget; ADR-0006 requires tolerance changes and warnings in `OperationResult`, with `Context` warning promotion. ADR-0009 defines edge certificates. This proposal does not relax those requirements.

### 6. Composition preserves public guarantees

An operation composed of multiple modelling steps shall specify a contract for the whole operation. Successful intermediate steps shall not become visible if a later step fails, unless the API explicitly defines a sequence of separately committed operations.

Private intermediate structures may temporarily violate final-model invariants when permitted by their internal contracts. They shall not be exposed as valid public models, and algorithms accessing them must respect the internal representation's assumptions.

Batch operations shall state whether the batch is atomic or whether each item commits separately. A partially completed batch must not be presented as atomic success.

### 7. Contract changes are API changes

Changes to accepted inputs, output semantics, failure guarantees, or promised identity behavior shall be reviewed as API compatibility changes even when function signatures remain unchanged.

Algorithm substitutions and performance improvements are permitted when the observable contract remains satisfied. Regression tests shall therefore target documented semantics rather than undocumented details of one algorithm's output.

## Illustrative contract: regularized solid intersection

This example illustrates the required specificity. It does not select a Boolean implementation or claim existing Stereon support.

| Field | Example specification |
| --- | --- |
| Purpose | Compute the regularized intersection of the occupied regions of two inputs: the closure of the interior of their set intersection. |
| Inputs | Two bounded valid solid regions within the supported representation and numerical domain; a region may contain disconnected components. |
| Parameters | An explicit or documented default numerical context, including applicable tolerances and execution limits. |
| Output | An empty region or a valid collection of solid components representing the regularized intersection under the declared numerical guarantees. |
| Boundary-only contact | Isolated point, curve, or surface contact contributes no solid volume to the regularized result. If the implementation cannot resolve the contact reliably, it reports failure. |
| Degenerate cases | Coincident and tangent boundaries are either supported under specified rules or produce a documented unsupported/unresolved outcome. |
| Mutation | Inputs remain unchanged on both success and failure. A successful result is published separately; any shared storage must obey the API's aliasing policy. |
| Failure | No result is published as valid. A diagnostic identifies the failure category and the affected stage or entities when available. |
| Representation | ADR-0007 G1 reproducibility applies. Cross-version face decomposition or identifiers are not promised; history follows ADR-0011. |

The regularized semantics in this example deliberately differ from a general set intersection that preserves lower-dimensional contact. A future API offering the latter must state that distinction.

## Consequences

### Benefits

- Applications can recover from failures without inspecting or rebuilding possibly corrupted inputs.
- Empty, unsupported, and numerically unresolved results have distinct meanings.
- Tests can assess semantic correctness across algorithm changes.
- Numerical limitations and representation restrictions become visible engineering constraints.
- Composite operations can preserve the same guarantees as primitive operations.

### Costs and obligations

- Contracts require maintenance alongside implementation changes.
- Staging and validation introduce runtime and memory costs that must be measured.
- Diagnostics require deliberate design rather than incidental error strings.
- Shared geometry, caches, aliases, and identity updates must participate in failure-safe publication.
- Supported-domain limits must remain honest even when a best-effort output looks plausible.

## Verification

As each operation becomes supported, its implementation and tests shall demonstrate:

1. The required contract fields are present and agree with the public API.
2. Representative successful results satisfy the declared structural and geometric postconditions.
3. Invalid, unsupported, and numerically unresolved cases produce distinguishable outcomes where applicable.
4. Valid empty results are distinguishable from failure.
5. Controlled failures at meaningful stages, including before publication, leave inputs and published state unchanged through all supported aliases.
6. Cancellation and resource-limit outcomes satisfy the same failure guarantee when those mechanisms are exposed.
7. Approximation or healing occurs only under the selected documented policy.
8. Multi-step operations do not leak intermediate changes when a later step fails.
9. G1 tests compare results bit for bit across schedulers/thread counts. Cross-version semantic tests avoid assuming stable decomposition or identifiers unless promised.

These are requirements for future implementation, not assertions that the repository currently meets them. Validation should use representative adversarial cases and targeted failure injection; this ADR does not require exhaustive geometric certification of every possible input.

## Relationship to existing decisions

- [ADR-0005](0005-immutability-and-sharing.md): immutable inputs, builders, and structural sharing.
- [ADR-0006](0006-error-handling.md): result types, error categories, diagnostics, and allocation failure.
- [ADR-0007](0007-threading-and-determinism.md): G1 reproducibility; not optional per operation.
- [ADR-0016](0016-cpp26-adoption-and-contracts.md): C++ contracts for programmer obligations only.
- [ADR-0018](0018-validity-checking.md): validation levels and build/Context policies. Each operation names applicable levels and guarantees; no new universal full-certification claim is introduced.
- [ADR-0021](0021-shape-derivatives.md): per-operation sensitivity declarations.
- [ADR-0025](0025-separate-topology-from-geometry.md): separation of structural and geometric responsibilities.

This proposal supplements these decisions and supersedes none. Specific operation domains, cancellation/resource-budget mechanisms, validation-evidence APIs, and performance targets remain implementation work. Already specified result types, tolerances, immutability, and determinism are not deferred or reopened.

## References

1. Ian Stroud, *Boundary Representation Modelling Techniques*, Springer, 2006. DOI: [10.1007/978-1-84628-616-2](https://doi.org/10.1007/978-1-84628-616-2). Relevant reading: “Stepwise modelling algorithms” and “Modelling operator definition.” The publisher's contents and scope have been verified; detailed chapter-level verification remains pending. This ADR does not attribute its failure guarantees or API design to a specific passage.
2. Companion proposal: [ADR-0025: Separate Topology from Geometry](0025-separate-topology-from-geometry.md).

