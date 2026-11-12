# Computational geometry convergence roadmap

**Status:** Proposed planning overlay, 2026-11-11. **No new release numbers or acceptance claims.**  
**Owner:** Onur Tuncer  
**Governance:** [Release roadmap](roadmap-and-versioning.md), [ADR-0055](../adr/0055-geometric-assurance-profiles.md), [assurance inventory](geometric-assurance-inventory.md)

## Objective

Evolve Stereon from a precise B-rep CAD kernel into a **computational-geometry-grade engineering kernel**: supported operations state what they have decided, constructed, covered, validated and left unresolved. Reuse geometric witnesses across modeling, topology, CAM and CAD–CAE without confusing exact predicates with exact geometry, sampled proxies with complete representations, or B-rep validity with functional correctness.

This **overlays** the existing M1–M7 milestones and v0.2–v1.0 release mapping. It does not create an eighth milestone, expand v1.0 into OCCT/CGAL parity, prescribe a single universal result type, or approve a new production geometry representation.

## Architectural workstreams

| Track | Research/engineering objective | Existing authorities | First evidence |
| --- | --- | --- | --- |
| A. Robust geometry | Trustworthy geometric decisions, complete branch discovery and bounded constructions | ADR-0002, 0010, 0028, 0029, 0047, 0050, 0055 | v0.2 handler assurance inventory and difficult-case corpus |
| B. Computational topology | Explicit incidence, event classification, lineage and validated publication | ADR-0025, 0027, 0030, 0040, 0051, 0054 | v0.3 topology queries and event witnesses; v0.4 Boolean reconstruction |
| C. Reusable geometric queries | Representation-independent meaning with representation-specific guarantees | ADR-0009, 0044, 0049, 0055 | v0.3 analytic rays; v0.5 trimmed queries; v0.7 conversions/meshing |
| D. Manufacturing geometry | Bounded contact, clearance, tri-dexel and offset queries | ADR-0044, 0049, 0050 | v0.3 analytic subset, v0.7 machining benchmark |
| E. Differentiable CAD | Valid regular-branch sensitivities and explicit critical events | ADR-0021, 0045, 0052, 0054, 0055 | v0.2 intersection linearisation; v0.9 workflow gate |

### Dependency order

1. Trust predicates and intersection completeness before allowing their results to drive discrete topology decisions.
2. Prove topological consistency and operation-scoped lineage before advertising broad Boolean success.
3. Export reusable geometric queries only with a clearly stated domain, error model and representation boundary.
4. Integrate derivatives only on declared regular branches, with event detection marking unsupported transitions.
5. Research alternative representations *after* independent oracles and invariant comparisons exist; a mesh is not an exact B-rep by assertion.

These are logical dependencies. Independent work may proceed in parallel (e.g. analytic ray/face queries need not wait for generic SSI).

## Release-by-release work packages and exit evidence

| Release | Incremental computational-geometry capability | Required evidence / release boundary |
| --- | --- | --- |
| **v0.2 — Intersections** | Publish assurance profiles for existing 2D curve–curve, 3D curve–curve, curve–surface, elementary and generic SSI families. Audit filtered predicates, root isolation, contact handling, bounds and proved empty results. Backfill metamorphic regressions (ADR-0053). Begin supported analytic contact/event witnesses (ADR-0054). | Reviewed [inventory](geometric-assurance-inventory.md), linked corpus and CI, exact vs bounded vs unresolved distinctions; existing [v0.2 ledger](v0.2-acceptance-ledger.md) remains binding. No claim that every spline tangency is decided. |
| **v0.3 — Solids** | Structured adjacency/incidence and `ShapeQueries`, verified topology/geometry validation, operation-scoped lineage, B-rep event consequences for supported edits. Analytic ray-face and bounded tri-dexel conversion. | Validity/correspondence and replay fixtures; no event inferred from entity count alone. Preserve 0.3-T3 and 0.3-F1/F2 gates. |
| **v0.4 — Planar Booleans** | Evidence-bearing splitting, classification, region reconstruction and lineage for planar/polyhedral operations. | Complete supported-domain examples including contact/overlap, independent set/volume checks, topology checks, failure-on-unresolved, metamorphic identities. |
| **v0.5 — Curved Booleans** | Selected analytic/trimmed surface intersection and Boolean cases; trim-aware closest points and bounded conversion. | Seam/pole/tangent/coplanar fixtures, explicit unresolved cases, certified or otherwise declared error/completeness limits. Do not infer broad freeform success. |
| **v0.6 — Robust Boolean gate** | Broaden verified engineering-part coverage; aggregate assurance and tolerance/lineage evidence across pipeline stages. | Existing Phase 4 corpus/benchmark gates, reproducible OCCT comparison with pinned environment, invalid-success count, declared unsupported domain. |
| **v0.7 — Meshing/CAM** | Error-bounded mesh-to-B-rep associations, tri-dexel fields, tool-contact and clearance queries with reproducible performance testing. | Bound conversion/discretization errors separately from sampling/grid resolution; benchmark B5.7 measures timing, memory, accuracy and failures. |
| **v0.8 — Exchange** | Imported/exported geometry preserves or explicitly changes assurance and provenance; representation conversions report lost guarantees. | STEP AP242/healing fixtures, round-trip geometric/topological checks and appropriately bounded import tolerances. |
| **v0.9 — Finishing/sensitivities** | Fillet/offset/shell topology-event reporting; differentiable regular-branch operations and engineering integration. | Analytic/finite-difference derivative comparisons in regular regions, event-driven invalidation at critical configurations, bounded finishing domain. |
| **v1.0 — Stable supported kernel** | Publish an audit-backed supported-domain catalog for geometry and topology operations; stable API and reproducible qualification. | Exact release commit, linked inventory, CI/corpus, comparison builds and release report; no implied universal certified geometry or OCCT/CGAL feature parity. |

## First execution queue (highest priority)

| Order | Task | Input / dependency | Deliverable | Evidence required |
| --- | --- | --- | --- | --- |
| 1 | Audit the v0.2 assurance matrix | ADR-0029, ADR-0047, [v0.2 ledger](v0.2-acceptance-ledger.md) | Supported domains and `decided/bounded/unresolved` claims by handler | Line-linked implementation, fixtures, and current CI per family |
| 2 | Finish historical metamorphic promotion | ADR-0053, [backfill](intersection-metamorphic-backfill.md) | Direct + metamorphic + independent oracle for each confirmed defect | Red-on-defect/green-on-repair, seed and revision links |
| 3 | Stage analytical event witnesses | ADR-0054, existing intersection outcomes | Candidate/certified/unresolved tests for circle/circle, sphere/sphere, plane/sphere | Verified critical parameter interval; no false topology assertion |
| 4 | Close risky intersections | v0.2 ledger's spline tangency, overlap, knot and periodic gaps | Supported resolution or explicit correctly-scoped failures | Adversarial/regression corpus, error bounds and completeness explanation |
| 5 | Integrate topology event consequences | ADR-0030, 0051 and v0.3 topology builders | Affected entities, lineage and B-rep validation | Independent topology invariants and event tests |
| 6 | Audit Boolean evidence | v0.4–v0.6 and ADR-0013 | Per-stage assurances and pinned OCCT comparisons | Same domains, tolerances, versions, invalid-success counts |
| 7 | Explore hybrid acceleration in research | Validated B-rep baseline, representative cases | Research specification with bounded conversion and failure analysis | Go/no-go decision before any new production representation |

Tasks 1–4 are the near-term queue; other rows are dependencies for later release gates. Marking a task in a planning document does not imply any code was implemented or verified.

## Assurance verification rules

- Every supported family separates *decision validity*, *construction accuracy*, *completeness*, *topological validity*, *event/derivative regularity* and *provenance* per ADR-0055. Not all are applicable to every query.
- Direct analytic references and independent geometric checks are the primary oracles. Metamorphic relations and OCCT/CGAL comparisons are supplemental, not proofs of the result.
- Tests transform coordinate frames, operand order, scale (with correctly transformed dimensional policy), parameter charts and scheduling only where those transformations preserve the specified problem.
- Test suites must distinguish incorrect success, correct decided results, correctly unresolved cases and unsupported requests. Review failure rates with the supported-domain denominator.
- Comparison runs pin code version, compiler, platform, numeric policy, input corpus, output equivalence test and budgets; report time, memory, completeness/accuracy and invalid success separately.
- A release report identifies the exact inventory revision and row status. A result not independently assessed is **not verified**, even if an earlier ledger describes an implementation.

## Post-v1.0 research opportunities (not release commitments)

**Certified curved surface intersections:** target topology completeness and bounded deviation simultaneously on difficult bicubic/NURBS cases, with explicit scalability limits.

**Parameterized topology-event localization:** isolate contacts, loop births/deaths and collapses over a parameter interval; preserve direction-aware semantics and operation-scoped lineage.

**Hybrid B-rep / mesh / implicit reasoning:** investigate conservative detection and error-bounded recovery, with mathematical arguments for any topological guarantees. Open a research specification (per [research rules](../research/README.md)) before selecting a production architecture.

**CAD–CAE and differentiated operations:** extend ADR-0021's sensitivities only to supported regular operations; event-aware optimization is a consumer of geometry guarantees, not an excuse to force derivatives across nonsmooth topology changes.

**Parallel query acceleration:** investigate broad-phase spatial structures, batched predicates and deterministic parallelism after correctness and reproducible baselines exist. GPU support is a research option, not an implicit new v1.0 requirement.

## Maintenance

- The [primary roadmap](roadmap-and-versioning.md) owns release numbering and M1–M7 phase gates.
- [ADR-0055](../adr/0055-geometric-assurance-profiles.md) owns the proposed cross-family assurance decision; preexisting accepted ADRs still govern their specialized operations.
- The [inventory](geometric-assurance-inventory.md) owns assurance coverage and evidence status; the [v0.2 ledger](v0.2-acceptance-ledger.md) owns v0.2 release acceptance.
- The [metamorphic backfill](intersection-metamorphic-backfill.md) owns the historical intersection-defect promotion queue.
- Review this overlay on each release gate and whenever an ADR changes relevant numerical or topology semantics.
