# ADR-0055: Geometric assurance profiles and representation-independent query contracts

- **Status:** Proposed
- **Date:** 2026-11-11
- **Owner:** Onur Tuncer
- **Phase:** 2 (v0.2 inventory and intersection evidence), 3–6 (topology, Booleans, machining, exchange and sensitivities)
- **Related:** [ADR-0002](0002-numeric-robustness.md), [ADR-0006](0006-error-handling.md), [ADR-0007](0007-threading-and-determinism.md), [ADR-0009](0009-parameter-space-representation.md), [ADR-0010](0010-ssi-strategy.md), [ADR-0013](0013-boolean-classification.md), [ADR-0021](0021-shape-derivatives.md), [ADR-0025](0025-separate-topology-from-geometry.md), [ADR-0026](0026-modelling-operations-by-contracts.md), [ADR-0028](0028-purpose-specific-numerical-tolerances.md), [ADR-0029](0029-intersection-and-classification-outcomes.md), [ADR-0030](0030-entity-lineage-in-operation-results.md), [ADR-0041](0041-layered-validation-evidence.md), [ADR-0044](0044-brep-to-tri-dexel-conversion.md), [ADR-0047](0047-numerical-operation-completion-and-publication.md), [ADR-0049](0049-machining-collision-and-clearance-queries.md), [ADR-0050](0050-tolerant-topology.md), [ADR-0053](0053-metamorphic-testing.md), [ADR-0054](0054-topological-event-detection.md)
- **Planning:** [Geometric computing convergence](../planning/geometric-computing-roadmap.md); [assurance inventory](../planning/geometric-assurance-inventory.md)

## Context

Stereon's planned progression from intersections through solids, Booleans, machining and differentiable geometry increasingly exposes geometric queries to consumers outside a conventional modeling operation. A decision that a pair of shapes intersects, a computed intersection curve, a validated Boolean, a clearance query and a shape sensitivity are different promises. For instance, an exact orientation decision does not imply an error-bounded constructed curve; a small geometric residual does not prove that every intersection branch was found; a watertight mesh does not, by itself, prove that it represents the intended B-rep solid.

The existing ADRs already define the authoritative, operation-specific requirements: robust arithmetic (0002), exact/error outcomes (0006), surface intersection completeness and accuracy (0010), Boolean classification (0013), derivatives (0021), intersection results (0029), lineage (0030), numerical completion and publication (0047), tolerant topology (0050), metamorphic testing (0053) and topological events (0054). This ADR does **not** replace those decisions. It establishes a common language for their promises and requires a reviewable, release-scoped inventory, so applications do not silently treat an estimated result as certified or a narrow-domain algorithm as general.

## Options

1. **Leave assurance entirely local to each API.** Preserves specialization but makes cross-family consumption and release qualification difficult to compare.
2. **Create one universal `CertifiedGeometryResult` and require proofs for every operation.** Attractive uniformity but false equivalence between decision, construction, completeness and topology; invasive ABI change and high runtime cost.
3. **Introduce a multidimensional geometric-assurance profile for documented operation families, using existing result types and evidence contracts.** Makes guarantees comparable without demanding a universal solver, arithmetic model or result wrapper.

## Decision

**Adopt option 3.** For each supported geometry-query or geometry-producing operation family, Stereon shall publish a **geometric-assurance profile** describing the domain and independently justified guarantees. Profiles are a documentation and verification obligation first. Their implementation may reuse existing result fields, error codes, validators and diagnostics; they do not require a universal public C++ type or a new global certification layer.

### 1. Assurance has independent dimensions

A profile MUST address the following, even if the answer for a dimension is `NotApplicable`, `NotSupported` or `Unverified`:

| Dimension | Required meaning |
| --- | --- |
| Supported domain | Geometry types, trims, parameter ranges, scale/numeric limits, dimensional assumptions, and excluded degeneracies. |
| Decision validity | What predicates/classifications are established; exact, rigorously bounded, estimated, or unresolved. |
| Construction accuracy | What metric and domain a geometric error bound applies to, including units and how the bound is obtained. |
| Completeness | Whether all components/roots/events within the requested domain are accounted for, and the reason empty output is valid. |
| Topological validity | Which incidence, connectivity, orientation and manifold/solid invariants are established, when a B-rep or related complex is produced. |
| Regularity and events | Conditions permitting smooth tracking/sensitivities; known critical configurations and unresolved event classes where applicable. |
| Provenance and correspondence | Which input entities, snapshots, representations and numeric policies justify the result; lineage when established by ADR-0030. |
| Computation and evidence | Budget/failure modes, test/corpus references, implementation revision, CI and comparison evidence, including unverified gaps. |

These dimensions MUST NOT be collapsed into a single `certified` Boolean, a confidence percentage or an overall score. Distinguish an exact predicate on nominal floating-point inputs from exact real-world measurement, and a rigorous enclosure from a sample-based error estimate.

### 2. Existing operation semantics remain authoritative

- ADR-0029 continues to define dimension-aware intersection results, completion and geometric contact evidence. This ADR does not declare `Decided` or `Certified` for cases those handlers leave unresolved.
- ADR-0047 remains authoritative for numerical construction and publication; input validity, geometric fidelity and task completion remain distinct.
- ADR-0013 and topology-validity contracts govern Boolean results. A successful geometric curve without complete fragment classification does not license a successful Boolean.
- ADR-0054 handles event detection as an observational service. A singular Jacobian alone is not proof that topology changed.
- ADR-0021 governs derivative validity; successful primal computation does not establish a valid derivative.
- ADR-0050 governs tolerant equality and growth. No new global epsilon, hidden snapping, or automatic escalation from tolerance-based coincidence to exact coincidence is introduced.
- ADR-0006 governs operation errors, including unsupported or unresolved cases. When requested assurance cannot be established, the operation MUST fail or report the explicitly limited result permitted by its established API contract; it MUST NOT silently return a stronger success.

### 3. Representation-independent *meaning*, representation-specific guarantees

Geometric queries SHOULD declare geometric semantics independently of how they are implemented (analytic equation, NURBS, trim, interval subdivision, controlled tessellation or later implicit representation). This does not require identical algorithms, guarantees or performance for all representations.

At any representation conversion or hybrid algorithm boundary:

1. Identify authoritative input geometry and the transformed/derived representation.
2. Declare the approximation metric, units, full domain, orientation and known error or coverage bounds.
3. Preserve or explicitly invalidate entity/face correspondence and topology claims.
4. Do not infer exact intersection, separation, watertightness or Boolean topology from a discrete proxy without a documented soundness argument.
5. If the required proof is unavailable, return the existing unresolved/unsupported outcome rather than silently relaxing a query.

Use of hybrid B-rep/mesh or implicit methods in the **production** kernel is a later scoped decision, supported by a research specification and comparison evidence. This ADR is not blanket approval to substitute approximate proxies for exact CAD geometry.

### 4. Predicate and construction separation

Operations that make a discrete branch or topology decision MUST identify the evidence sufficient for that decision. Approximate construction of points or curves is independently checked against the promised geometric error metric; neither predicate robustness nor a small residual alone establishes completeness. Fast filtered evaluation and robust fallbacks are permitted by ADR-0002, and expensive certification MAY be conditional on the requested assurance and the difficult-case path.

A geometric witness (e.g. contact point, intersection branch, orientation predicate, separation bound) MAY support multiple consumers, but each consumer must check that the witness's domain, assurance and context satisfy its own contract.

### 5. Evidence and release discipline

Maintain [the assurance inventory](../planning/geometric-assurance-inventory.md) for each declared operation family, initially intersections. An entry records:

- the published promise and excluded cases;
- the implementing file(s) and relevant API / ADR;
- direct analytic or independent oracle cases, metamorphic relations, known defects, and comparator fixtures where applicable;
- coverage/completeness evidence and unresolved outcomes;
- the test revision, CI run and reviewer approval that establish the gate.

Inventory entries MUST distinguish **documented intent**, **implemented behavior**, and **verified release evidence**. A planning row is never a certification claim. Every release report links the inventory snapshot applicable to its supported domain; no test run or measured advantage is inferred from an ADR alone.

Requirements and status are maintained per supported family. Do not block an entire release on unimplemented broad-domain methods that the release explicitly excludes. Conversely, a documented limitation is not an excuse for silently passing invalid results within the promised domain.

### 6. Composable queries without architecture erosion

Stereon SHOULD expose reusable predicates, spatial queries and validation witnesses through appropriate existing core/robust/intersection/topology modules. Reuse existing classification methods and provenance structures rather than creating a second intersection engine or a parallel B-rep store.

This ADR does not mandate one monolithic geometric-computing module, a single public result wrapper, runtime certificate storage for every hot-path query, new entity identities, or tight coupling between geometry and topology. New public interfaces require separate API review under ADRs 0003/0019/0020/0025/0026.

### 7. Research scope

The assurance inventory supplies reproducible benchmarks for three research questions without predetermining the answers:

- **Certified curved intersection:** how to establish both branch completeness and bounded geometric error at tractable cost.
- **Computational topology:** how to localize critical configurations and validate event/topology consistency with qualified lineage.
- **Representation-aware geometric queries:** when a derived mesh/implicit representation can accelerate a query while preserving the requested geometric and topological contract.

Novel algorithms that are not yet justified belong in `docs/research/`; adoption requires a separate decision. Differentiable CAD integration consumes the assurance/event evidence but remains governed by ADR-0021.

## Consequences

**Benefits:** geometric-query contracts can be compared across CAD/CAM/CAE uses; release gates become auditable; certificates and uncertainties remain honest; algorithms can evolve without conflating representation, accuracy and topology; a strong foundation exists for future computational-geometry and differentiable-CAD research.

**Costs:** maintaining profiles and a verified evidence inventory requires reviews, traceable tests and CI links. Some operations will legitimately remain `Unverified` or `Unresolved`, and public coverage claims may become narrower but more accurate.

**Risk:** generic assurance terminology could encourage superficial checkbox compliance. Mitigation: require concrete witness types, metrics, domains, independent checks and revision-bound evidence; no summary flag or claimed completeness from sampling.

## Verification

1. **Contract review:** choose one analytic intersection, one spline intersection, one Boolean/topology-producing operation and one manufacturing query; map each into distinct assurance dimensions without claiming that all dimensions are available.
2. **Counterexample test:** prove that an unresolved tangency, missed-branch fixture, or incomplete classification cannot be treated as complete success by a consuming operation in its supported domain.
3. **Evidence traceability:** each supported inventory entry links its implementation, supported-domain documentation, executable fixtures and CI evidence (or a visible missing-evidence status).
4. **Representation boundary checks:** for any new proxy-based query, construct independent cases that expose missing features, poor approximation, altered orientation or topology, and verify bounded result or explicit unresolved failure.
5. **Metamorphic and differential checks:** apply ADR-0053 where mathematically valid, with independent analytic or geometric oracles; transformed outputs are compared semantically, not by raw IDs or byte serialization.
6. **Benchmark reproducibility:** record comparable input domains, platform, build flags, accuracy targets, memory, timings and invalid-success counts; performance claims never displace correctness or completeness evidence.

## Staging and release policy

- **v0.2:** publish the first intersection assurance inventory and audit the currently supported handler families. Close or explicitly carry the v0.2 cross-cutting evidence gaps. ADR-0053/0054 Stage A provide metamorphic and event evidence where supported.
- **v0.3:** extend profiles to topology queries, lineage, analytic ray/face and tri-dexel conversion where those operations are delivered.
- **v0.4–v0.6:** attach profiles to planar, curved and broader Boolean release subsets, including topology validation and unresolved classification.
- **v0.7–v0.8:** extend profiles to mesh associations, machining queries, conversion and import/export; investigate hybrid methods in research before adopting new production semantics.
- **v0.9–v1.0:** require published supported-domain assurance for finishing and sensitivity workflows, including explicit event/derivative validity, before stable qualification.

These gates supplement, not replace, existing version milestones and ADR acceptance policies. **Acceptance of this ADR does not certify any current geometric implementation.**

## Stage-one acceptance checklist

- [ ] Inventory defines the supported-domain and assurance dimensions for each currently advertised v0.2 intersection family.
- [ ] At least one reviewed family per result category has direct oracle, metamorphic, unresolved and complete-empty tests.
- [ ] Evidence references an exact implementation revision, test fixtures and CI execution.
- [ ] A consumer test demonstrates that an incomplete/unresolved result does not produce an unjustified successful higher-level operation.
- [ ] No existing operation result, tolerance, topology, lineage or derivative contract is weakened.
- [ ] Release documentation distinguishes intended, implemented and verified guarantees.

**Decision summary:** Adopt release-scoped, multidimensional geometric-assurance profiles and representation-independent query semantics, while retaining operation-specific results, explicit uncertainty, numerical discipline and the current phase/release architecture.
