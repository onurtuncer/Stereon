# ADR-0054: Topological event detection and classification

- **Status:** Proposed
- **Date:** 2026-11-11
- **Owner:** Onur Tuncer
- **Phase:** 2 (intersection-event foundations); 3–4 (B-rep integration); subsequent parametric tracking and optimization stages
- **Related:** [ADR-0002](0002-numeric-robustness.md), [ADR-0006](0006-error-handling.md), [ADR-0007](0007-threading-and-determinism.md), [ADR-0010](0010-ssi-strategy.md), [ADR-0013](0013-boolean-classification.md), [ADR-0021](0021-shape-derivatives.md), [ADR-0025](0025-separate-topology-from-geometry.md), [ADR-0026](0026-modelling-operations-by-contracts.md), [ADR-0027](0027-controlled-topology-mutation-and-euler-operators.md), [ADR-0028](0028-purpose-specific-numerical-tolerances.md), [ADR-0029](0029-intersection-and-classification-outcomes.md), [ADR-0030](0030-entity-lineage-in-operation-results.md), [ADR-0032](0032-mutation-based-kernel-verification.md), [ADR-0040](0040-atomic-topology-edit-publication.md), [ADR-0041](0041-layered-validation-evidence.md), [ADR-0047](0047-numerical-operation-completion-and-publication.md), [ADR-0052](0052-space-curve-meeting-derivatives.md), [ADR-0053](0053-metamorphic-testing.md)

## Context

A continuously varying geometric model can undergo a discrete change in its intersection arrangement or B-rep topology. Examples include the birth or disappearance of an intersection loop, tangential contact, an edge or face collapsing, the splitting or merging of components, and creation or closure of a tunnel. Merely comparing endpoint B-reps can miss intermediate events. Numerical difficulty, near-zero residuals and Jacobian singularities can indicate candidate events but do not themselves establish a topological change.

Stereon already requires explicit intersection outcomes (ADR-0029), separate geometry and topology (ADR-0025), operation-scoped lineage (ADR-0030), purpose-specific tolerances (ADR-0028), and derivative-validity boundaries (ADR-0021). Without a common event contract, intersection handlers, Boolean topology builders, and future parametric optimization workflows can independently interpret the same critical configuration in inconsistent ways.

Two distinctions are essential:

1. An **intersection-arrangement event** changes the structure of a relation between two inputs; it need not alter either input object's own topology. Two tangent circles, for example, do not individually change topology when the number of their intersection points changes.
2. A **geometric degeneracy** (parallel normals, singular Jacobian, near-zero edge length) is evidence to investigate, not by itself proof that the intersection arrangement or B-rep topology changed.

The kernel should make supported events observable, attributable, testable and, where possible, certifiable without promising a universal solver for arbitrary parametric surfaces.

## Options

1. **Leave event handling implicit in intersections, Booleans and healing.** Minimal new API, but critical configurations are difficult to diagnose, compare, trace and reuse for optimization.
2. **Build an independent, universal topological-event engine.** Conceptually uniform, but duplicates intersection/classification machinery, introduces inconsistent predicates and has unrealistic near-term completeness requirements.
3. **Add a staged, cross-cutting event contract implemented by existing operation families.** Reuses established predicates and witnesses; defines classification, evidence, parameter localization and lineage without assuming universal coverage.
4. **Compare only pre-change and post-change topology.** Useful as an additional consistency check, but transient or simultaneous events can be missed, and the critical parameter is not isolated.

## Decision

**Adopt option 3.** Stereon shall support explicit detection and classification of topological and intersection-arrangement events, initially for documented analytic intersection cases, later for B-rep reconstruction and one-parameter geometry families. Detection is observational and does not mutate geometry or topology. No complete general-purpose event detector is implied by this ADR.

The rules below apply to each event class as it becomes supported. A feature that cannot certify an event must report limited evidence or an unresolved outcome, never invent certainty or silently publish changed topology on the strength of a heuristic candidate.

### 1. Distinguish event domains and classification certainty

Every event belongs to a declared domain:

- **Intersection arrangement:** creation, disappearance, split, merge or contact of intersection components.
- **B-rep topology:** changes in entity incidence/connectivity, component membership, or topological holes/tunnels arising from a completed operation.
- **Geometric degeneracy:** a noteworthy critical geometric condition without a demonstrated topology change.

The first two domains describe different objects and must not be conflated. Candidate detection shall not automatically cause B-rep mutation.

Event records shall distinguish:

- **Candidate:** a condition suggests that an event may occur;
- **Certified:** the specified event is established within the declared geometry, parameter domain and numerical policy;
- **Unresolved:** evidence is insufficient to determine whether or which event occurred.

A successful no-event result must establish the absence of the **supported** event classes over the declared parameter/domain coverage. Not finding roots or seeds is not a completeness proof. Unsupported geometries, exhausted budgets and ambiguous cases shall use existing explicit failure/unresolved conventions (ADRs 0006, 0029 and 0047); a partial observation must not masquerade as complete success.

### 2. Extensible event taxonomy

An initial semantic vocabulary shall cover:

| Event kind | Meaning |
| --- | --- |
| `Contact` | Previously disjoint entities attain contact; no implied input-topology change. |
| `IntersectionBirth` / `IntersectionDeath` | An intersection component appears / disappears along an oriented parameter sweep. |
| `IntersectionSplit` / `IntersectionMerge` | Intersection components branch or coalesce along an oriented sweep. |
| `EdgeCollapse` / `FaceCollapse` | A B-rep edge or face loses its valid geometric extent, requiring a separately verified topology response. |
| `ComponentSplit` / `ComponentMerge` | A B-rep connected component separates / joins in the resulting model. |
| `TopologyChange` | A verified B-rep topology change not yet covered by a specialized kind. |
| `DegenerateConfiguration` | A geometric degeneracy with no demonstrated topology transition. |
| `Unclassified` | A detected critical configuration whose event class is unresolved. |

This is a **conceptual taxonomy**, not a mandate to introduce these exact C++ identifiers before the corresponding operation exists. Birth/death and split/merge are directional labels; parameter reversal must transform them appropriately. Coincidence/overlap may persist over an interval and must not be reduced to a single isolated event.

Global Betti numbers or counts of faces/edges alone are not sufficient to classify all local changes. Local adjacency, incidence, curve graph structure, orientation and supported topological invariants may be required.

### 3. Result and witness contract

A supported detector shall produce a structured event result describing, where applicable:

| Field | Required semantics |
| --- | --- |
| Kind and domain | Event class and whether it refers to geometry, intersection arrangement or B-rep topology. |
| Assurance / completion | Candidate, certified or unresolved; declared coverage and outstanding work. |
| Critical parameter | A bounded interval enclosing an event parameter for a tracked family; omitted when not applicable. |
| Affected entities | Qualified references to participating geometric/topological entities and owning snapshots. |
| Witness | Contact point, critical curve, parameter tuple, singularity test, adjacency change, or other justified evidence. |
| Lineage | Operation-scoped entity correspondence where an operation can establish it under ADR-0030. |
| Diagnostics | Predicate/solver outcome, numerical policy, error bounds or limitations, and relevant operation context. |

Critical values shall be expressed as intervals when exact localization is unavailable; approximate roots shall not be represented as certified exact values. Event witnesses must respect ADR-0029's distinction between geometric components, accuracy and completeness.

Entity references and retained snapshots shall follow ADR-0030: no global identifiers inferred from equal indices, pointer addresses or geometric similarity. A detector may provide prospective affected entities before publication, but a finalized lineage manifest belongs to the operation result that published topology.

### 4. Detection is separate from mutation and publication

Detectors are read-only with respect to input geometry and topology. They may be invoked during an operation, but do not own Euler operators, healing, Boolean assembly or feature regeneration. Existing modeling operations remain responsible for decisions and for publishing complete, valid topology (ADRs 0027 and 0040).

A candidate event cannot authorize topology mutation by itself. Operations needing a definite classification must obtain evidence at their required assurance level or return an explicit unresolved/failure outcome, consistent with ADRs 0029 and 0047.

### 5. Reuse robust geometric predicates and isolate critical parameters

Candidate triggers may include tangency, a loss of Jacobian rank, degenerate trims, near-zero feature extent or intersection-graph changes. A rank deficiency is a **candidate**, not sufficient proof of a topology transition. Event classification must inspect relevant local neighborhoods or before/after sides under a documented method.

For parameterized geometry, techniques may include analytic critical equations, root isolation, adaptive subdivision, interval arithmetic and adaptive-precision predicates. Sign-change sampling alone is insufficient: tangencies and even-multiplicity events may not change sign. Supported tracking must specify its parameter domain, coverage guarantee, work limits and numerical policy.

Events with overlapping critical-parameter intervals cannot be arbitrarily ordered. Simultaneous or unresolved event groups must be representable. An endpoint comparison is an independent check, not a replacement for tracking.

### 6. Numerical and architecture policy

Detection shall use Stereon's purpose-specific tolerance policy (ADR-0028) and numeric robustness strategy (ADR-0002). Distinguish root localization, approximation, classifier resolution and bounded topology repair. Do not introduce a new global epsilon or silently upgrade tolerant coincidence to exact equality (ADR-0050).

Keep geometry detection separate from topology classification (ADR-0025); reuse the existing intersection and classification outcome model (ADR-0029), context/determinism rules (ADR-0007), and operation result/lineage model (ADR-0030). Prefer value-oriented results, bounded computation and explicit outcomes in C++26 modules; module names are an implementation decision, not fixed here. Event output is deterministic for the same declared input, context and execution guarantees.

### 7. Relationship to shape derivatives

ADR-0021 already defines derivatives as operation contracts. This ADR **extends, but does not replace**, that decision: topology events identify boundaries across which ordinary smooth-branch derivatives may no longer be valid. An event indication may invalidate or condition an `IntersectJvp`/`IntersectVjp` product, but it does not supply a derivative through a topology change.

Respect ADR-0052's explicit choice that meetings of space curves do not acquire a derivative from `Intersect`; any future `ClosestApproach` derivative is a separately named product. Generalized/nonsmooth derivatives and automatic topology-changing optimization require an additional decision if proposed.

## Verification

Each supported event kind shall have independently specified reference configurations and negative controls, with a documented domain and expected assurance. Minimum evidence:

1. **Analytic reference cases:** externally tangent and separated circles/spheres, plane–sphere tangency, controlled cylinder intersection cases, and an independently constructed intersection-component birth/death example; do not assert a case is certified without an algorithm that proves it.
2. **Metamorphic invariants:** rigid-motion equivariance of witnesses and invariance of event parameters/class; operand symmetry where applicable; monotone invertible reparameterization of critical intervals; correct birth/death reversal under sweep reversal; subdivision consistency for certified events.
3. **Numerical stress:** near-tangencies, coincidence/overlap, very small features, translated and scaled coordinates, multiple nearby or simultaneous events, trimming-boundary degeneracies and exhausted search budgets.
4. **Topology/lineage evidence:** where a B-rep update is involved, validate resulting incidence, component/connectivity changes and affected-entity correspondence against ADR-0030 and topology validation contracts.
5. **Regression, metamorphic and mutation testing:** every confirmed intersection-event defect inherits ADR-0053's mandatory defect-promotion contract: a retained defect-ledger entry, minimal direct reproducer, an applicable metamorphic regression with an independent oracle, complete context and replay data, and red-on-defect/green-on-repair evidence. Where no sound relation applies, retain a reviewed mathematical exception explaining each considered relation and a stronger independent regression. Historical defects remain open obligations until that evidence is reviewed. Sweep reversal, monotone reparameterization and subdivision relations supplement ADR-0053's catalogue; they do not replace its evidence or release requirements. Include classification-logic mutation tests under ADR-0032 where practicable.
6. **Performance:** benchmark candidate-trigger overhead, localization/classification cost, and unresolved-event rate against supported geometry families. Certification need not be invoked on every uncomplicated query.

An event detector must not report a certified no-event result for an interval merely because no candidate was sampled. Metamorphic checks complement independent oracles and cannot substitute for geometric certification.

## Delivery and release gates

### Stage A — Intersection events (phase 2 / v0.2 as coverage permits)

Introduce internal candidate/witness semantics using existing intersection handlers. Support specific analytically verifiable tangencies and birth/death cases; add explicit candidate/certified/unresolved status and metamorphic tests. Do not make all of ADR-0054 a blanket v0.2 gate.

**Gate:** supported analytic cases have independent tests; certified and unresolved outcomes are not conflated; no unsupported completeness assertion. Intersection-event defects and applicable handler relations participate in ADR-0053's machine-readable inventory and strict release-evidence check; an open obligation cannot count as a pass.

### Stage B — B-rep events (phases 3–4)

Integrate event observations with intersection graph reconstruction, edge/face collapse, topology operators, Boolean publication and operation lineage.

**Gate:** affected entities and published topology are consistent with validated B-rep results, with explicit failures where resolution is insufficient.

### Stage C — Parameterized event tracking (later)

Locate critical parameter intervals, track multiple/simultaneous events, and support parameter sweeps without silently skipping declared supported event classes.

**Gate:** reference parameter families and reversed/reparameterized sweeps have consistent bounded outcomes.

### Stage D — CAD–CAE and derivative integration (later)

Consume event data through ADR-0021's derivative-status contract and future analysis or optimization workflows.

**Gate:** derivative validity is correctly restricted to declared regular branches; event detection never fabricates generalized derivatives.

The capability begins as internal infrastructure. Promote it to a public API only after semantics, lifetime/identity rules, certification levels and supported-domain documentation stabilize. Unsupported event classes do not become release blockers unless a specific dependent modeling operation requires them for correctness.

## Consequences

**Benefits:** better diagnosis of difficult intersection and Boolean cases; explicit numerical and topological assurance; reusable provenance for topology-changing operations; a foundation for robust parameter sweeps, manufacturing analysis and topology-aware differentiation; stronger automated regression evidence.

**Costs:** classification logic, witness storage, root localization and numerical stress testing add work and possible runtime overhead. Generic surfaces may remain unresolved. The main failure risks are treating a heuristic degeneracy as a certified event and duplicating existing intersection logic; both are prohibited by the contracts above.

**Non-goals:** a universal certified event solver; automatic healing or mutation performed by the detector; globally persistent topology naming; differentiability through arbitrary topology changes; exact ordering of unresolved simultaneous events.

## Acceptance criteria (Stage A)

- [ ] Supported event domains, status and witness contracts are specified in the chosen code interfaces.
- [ ] Initial analytic intersection handlers report supported candidates and classifications, including unresolved outcomes.
- [ ] Critical parameter localization returns a justified interval or explicit limitation where applicable.
- [ ] Independently known analytic positive and negative cases have regression tests.
- [ ] Applicable rigid-motion, operand-symmetry and parameter-sweep metamorphic relations are tested.
- [ ] Numerical near-degeneracy and resource-limit behavior is tested.
- [ ] The existing intersection result, tolerance, lineage and publication contracts are preserved.
- [ ] Benchmark evidence measures impact on the primal intersection path.

**Decision summary:** Adopt explicit but staged topological-event detection; classify and substantiate supported events without mutating inputs, begin in intersections, then integrate topology, parameter tracking and derivative-validity consumers.
