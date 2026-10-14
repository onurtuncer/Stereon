# Stereon roadmap and versioning

**Planning revision:** 0.2 • 14 October 2026

**Owner:** Onur Tuncer  
**Status:** Proposed release plan; acceptance evidence is required before tagging a release.

This plan embeds six modelling milestones derived from Ian Stroud's *Boundary Representation Modelling Techniques* into the existing release sequence. It retains Phase 0 and the existing release anchors, adding explicit v0.2/v0.4/v0.5/v0.7/v0.9 gates between them. Milestone IDs M1–M6 below describe dependencies; they are not phase numbers, package versions, or the capability matrix's Mxx row IDs. R0–R5 remain long-term capability stages, and v1.0 does not imply comprehensive OCCT parity.

## Release mapping

| Release / implementation phase | Milestone contribution | Required outcome |
| --- | --- | --- |
| Phase 0 / 0.0.x | M1: executable modelling foundations, initial subset | Debug viewer, controlled builder prototype, separate validity checks, numerical policies, diagnostic dumps, replayable fixtures, and a minimal persistence prototype. Existing Phase 0 gates remain binding. |
| v0.1 Geometry / Phase 1 | M1 completion for available geometry; M4 geometry groundwork | Analytic and NURBS evaluation, fitting, projection and parameter-domain handling. Record tolerances and approximation limits; support save/load of the declared geometry subset. |
| v0.2 Intersections / Phase 2 | M2/M3 planar intersection services; M4 curved intersection groundwork | Inspectable intersection/classification results, including tangency, coincidence, overlap and unresolved/unsupported outcomes. Satisfy the applicable SSI gates before consuming these services in broad Booleans. |
| v0.3 Solids / Phase 3 | M2: planar construction and interrogation; M4: curved and trimmed construction | Planar faces with holes, polyhedral solids, extrusion, sectioning, imprinting and classification; then the advertised primitives, revolve, sweep and loft, trimming and sewing. Deliver operation history and native shape persistence for this supported subset. |
| v0.4 Planar Booleans / Phase 4a | M3: planar Boolean pipeline | Union, intersection and difference for the declared planar/polyhedral domain, using inspectable splitting/imprinting, classification and assembly; validate lineage and attribute propagation. |
| v0.5 Curved Booleans / Phase 4b | M4: bounded curved Boolean integration | Extend the validated planar pipeline to a declared analytic and trimmed-spline subset; check seams, poles, tangencies, overlap and intersection certificates. No broad real-part robustness claim yet. |
| v0.6 Booleans / Phase 4c | M3/M4 robust engineering Boolean gate | Broaden the declared real-part domain and pass the existing Phase 4 corpus, success-rate, comparison, tolerance-growth and determinism gates. |
| v0.7 Meshing / Phase 5a | M4/M6 engineering workflow integration | Non-destructive display/analysis tessellation, mesh-to-B-rep associations, tagged duct workflow and usable Python construction/Boolean/mesh calls; validate native migrations. Public STEP/healing coverage remains the v0.8 gate. |
| v0.8 Exchange / Phase 5b | M4 import/repair integration; M6 initial product/application services | STEP AP242 read/write with explicit healing, derived display/analysis meshes, boundary-condition tags, native migrations, and usable Python workflows. |
| v0.9 Finishing / Phase 6a | M5: local and finishing operations; M6 bounded application support | Bounded constant-radius fillets, chamfers, shelling/thickening and shape sensitivities; history/tag integration and Phase 6 algorithm gates. Public API remains pre-1.0. |
| v1.0 Stable kernel / Phase 6b | M1–M6 supported-domain release qualification | Stable documented public API, compatibility/migration notes, supported-platform qualification and end-to-end engineering acceptance. Re-run applicable release gates after final changes; do not expand to full OCCT parity by implication. |
| Post-v1.0 / R3–R5 workstreams as applicable | M5/M6 extended scope | Draft-angle operations, variable-radius blends, feature verification/recognition, full parametric naming/history, richer assemblies/documents, drawing services, unfolding, medial-axis methods and manufacturing queries, with separate users and acceptance datasets. |

Geometry and topology work may overlap. Curved geometry evaluation need not wait for planar Booleans. However, broad Boolean release claims depend on a validated planar pipeline, curved intersection/trim support, and a declared input domain. A bounded Phase 0 Boolean prototype is research evidence, not completion of v0.6.

## Milestones and exit evidence

### M1 — Executable modelling foundations

Deliver controlled edits, topology/geometry validation, purpose-specific numerical policies, diagnostic dumps, replay journals and minimal versioned save/load for the available representation. Euler operators establish topology invariants; they do not establish geometric validity (Stroud §4.2).

Exit evidence: seeded structural and geometric defects are distinguished; failed edits leave the input unchanged; successful edits pass the applicable validity level; supported values survive round trips with their units and numerical context. Phase 0 may use diagnostic/prototype persistence; production `.stn` persistence and compatibility gates remain due by v0.3 under ADR-0012. No prototype silently establishes a public archival format.

### M2 — Planar construction and interrogation

Deliver planar faces with holes, polyhedral solids, extrusion, sectioning, imprinting and point classification. Treat geometric coincidence separately from topological identification; touching assembly occurrences are not automatically sewn together.

Exit evidence: declared-domain cases cover cavities, disconnected components, contact configurations and features near resolution limits; section and imprint boundaries are independently checked; applicable Structural/Geometric/Solid checks pass. Classification reports boundary and ambiguous cases explicitly. This subgate precedes M3 broad Boolean work.

### M3 — Planar Boolean pipeline

Expose intersection discovery, splitting/imprinting, classification and result assembly as separately diagnosable stages. Build union, intersection and difference on the planar supported subset for v0.4 before expanding to curved solids in v0.5 and the broad v0.6 robustness gate.

Exit evidence: disjoint, contained, crossing, coincident, edge-contact and vertex-contact cases have explicit expected semantics; outputs pass topology/geometry checks and analytic mass-property checks; splits, merges and deletions have correct lineage; attribute conflicts are reported. Unsupported or unresolved cases return structured errors. The planar subgate does not replace the Phase 4 corpus, success-rate, comparison or determinism gates in ADR-0013 and benchmarks B4.1–B4.6.

### M4 — Curved and trimmed modelling

Stage work across v0.1 geometry, Phase 2 intersections, v0.3 construction, v0.5/v0.6 curved Booleans and v0.7/v0.8 mesh/exchange integration. Deliver analytic/spline geometry, trimming, periodic seams, degenerate boundaries, curve/surface intersection services and sewing with explicit identification policies.

Exit evidence: topology and geometric boundaries agree across seams, poles, tangencies and gaps; intersection branches and coincidence regions are checked within the supported domain; certificates and cached representations obey ADR-0009/0010; sewing reports free edges, rejected identifications and tolerance changes. Boolean integration requires the applicable M3 gate and curved intersection/trim checks; successful evaluator tests alone do not establish curved Boolean support.

### M5 — Local and finishing operations

Deliver bounded constant-radius fillets, chamfers, offsets and shelling/thickening for v0.9, qualified for stable supported-domain use at v1.0. Draft-angle operations and variable-radius blends remain post-v1.0 scope.

Exit evidence: supported configurations and limits are documented; collapse, self-intersection and topology-changing events are detected; results validate or return structured errors; lineage and tags survive edits. Shape sensitivities distinguish fixed-topology derivatives from nondifferentiable events per ADR-0021. Retain the Phase 6 benchmark gates B6.1–B6.4.

### M6 — Product and application services

Begin lineage and attribute infrastructure with v0.3; exercise it in v0.4–v0.6 and publish bounded mesh/Python workflows in v0.7, followed by exchange/healing in v0.8. Validate finishing-operation integration and supported interface stability for v1.0. Full feature recognition, assemblies, document services and manufacturing reasoning remain later capabilities.

Exit evidence: split/merge lineage composes through chained operations and compaction; metadata conflicts and invalidated references are visible; the tagged duct workflow preserves inlet/outlet/wall semantics through cut, validation, watertight mesh and external volume meshing (ADR-0011/0014). Persistent naming across regeneration remains a separate post-v1.0 problem; operation history is not a claim that it is solved.

## Architectural dependencies

The book motivates these acceptance requirements, not an automatic change of ADR status. Check and review the relevant records before implementing their gates.

| Concern | Existing record / planning action |
| --- | --- |
| Numerical policy and purpose-specific tolerances | ADR-0001/0002 and [ADR-0028](../adr/0028-purpose-specific-numerical-tolerances.md): purpose-specific controls, dimensional conversions, observable changes and persisted numerical context; policy in M1/v0.1, intersection use in v0.2, native shape context by v0.3. |
| Intersection/classification outcomes | ADR-0010/0013 and [ADR-0029](../adr/0029-intersection-and-classification-outcomes.md): dimension-aware results and independent side-transition evidence; v0.2 intersection gate, v0.3 queries and v0.4–v0.6 Booleans. RS-0001 remains research. |
| Entity lineage | ADR-0011 and [ADR-0030](../adr/0030-entity-lineage-in-operation-results.md): qualified many-to-many relations, explicit split/merge events and composition by v0.3, extended through v0.4–v0.6. |
| Attribute propagation | [ADR-0031](../adr/0031-attribute-propagation-and-conflicts.md): per-key split/merge/revalidation rules and conflict diagnostics in v0.3–v0.6; tagged engineering mesh gate in v0.7. |
| Coincidence versus topology identification | ADR-0004/0009/0018 plus construction/sewing contracts; record identification rules before expanding sewing or simplification. |
| Native schema and migrations | ADR-0012; prototype early, public shape persistence by v0.3. |
| Derived tessellation and tags | ADR-0014; mesh generation must not alter authoritative B-rep geometry/topology. |
| Reproducible failures | ADR-0006/0023; journals, dumps, numerical context and stage diagnostics accompany corpus failures. |
| Controlled edits and operation contracts | ADR-0025/0026/0027 alongside ADR-0004/0005/0018 define association, operation and controlled-edit contracts; review before phase acceptance. |

## Versioning and release gates

1. **Product version.** `CMakeLists.txt` is the source for the numeric kernel version (`LibraryVersion`, documentation and bindings follow it). The current 0.0.1 remains unchanged by this planning update. Phase 0 uses 0.0.x; the advertised releases target every minor from 0.1.0 through 0.9.0, followed by 1.0.0. Each minor has the bounded scope and acceptance gate in the release table. Phase 4 is split across 0.4–0.6, Phase 5 across 0.7–0.8 and Phase 6 across 0.9–1.0; implementation phases need not match minor numbers.
2. **Pre-1.0 changes.** New capabilities or incompatible public API changes after a tagged milestone require a minor-version increment with migration notes. Patch releases contain compatible fixes and documentation updates; 0.0.x remains explicitly experimental. Candidate tags use `vX.Y.Z-rc.N` and do not assert that the final gate has passed.
3. **Release evidence.** Tag a final milestone only after all mandatory supported-domain acceptance fixtures, applicable ADR gates, supported-platform CI, determinism checks and required benchmark gates pass. Publish the exact implementation commit, corpus revision, comparator build, domain limits, structured failures and observed invalid-success count in a release report. A benchmark target is not an achieved measurement. M1–M6 evidence complements existing gates; it does not weaken them.
4. **Post-1.0 API policy.** Within the declared stable public API, additive compatible capabilities increment the minor version, compatible fixes increment the patch, and incompatible changes increment the major. C++ source compatibility and binary compatibility are distinct: compiler/runtime/configuration requirements and rebuild rules follow ADR-0019. Stable API does not promise cross-compiler C++ ABI compatibility.
5. **Native format version.** `.stn` schema versions are independent of product versions. Follow ADR-0012: every schema change increments the format version; explicit migrations read older public formats; unknown required chunks fail and unknown optional chunks warn. Keep fixtures from each tagged public-format release. Prototype dump versions must be labelled separately from public `.stn` schema versions.
6. **C ABI version.** Version/capability queries and ABI-major compatibility follow ADR-0020 independently of product versions. Do not advertise a stable C ABI until its ownership, error and old-consumer compatibility gates pass. C ABI implementation is not a prerequisite for early geometry releases.
7. **Planning versions.** This document's revision and the capability matrix's document version identify planning edits, not delivered software. Phase numbers, M1–M6 milestone IDs and R0–R5 capability stages must never be presented as package versions or completion percentages.

Use an acceptance ledger per release: milestone/subgate, supported domain, owner, required ADR decision, fixture IDs, implementation commit, result/report link and unresolved limitations. All gates above are planned requirements; no new capability completion or audit grade is asserted by this document.

## Source and related plans

Stroud (2006), *Boundary Representation Modelling Techniques*: chapters 3–7 for representation and stepwise operations; §8.3 for information preservation; §10.3 for faceting; §§11.1–11.3 for persistence; §13.5 for sewing; §§16.1–16.3 for tolerances and diagnostics. Modern staging, immutable mesh derivation and versioning details here are Stereon planning choices, not quotations or claims that the book specifies them.

- [Capability matrix and historical audit](stereon-capability-matrix.md)
- [Phase 0 implementation checklist](../../TODO.md)
- [Benchmark specifications and existing phase gates](../benchmarks.md)
- [Architecture decision records](../adr/README.md)
