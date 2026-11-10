# Stereon roadmap and versioning

**Planning revision:** 0.5 • 10 November 2026 (tolerant topology, ADR-0050, as interim steps of v0.2 and v0.3; face adjacency and canonical surface queries, ADR-0051, in v0.3 before the planar Booleans)  
**Previous revision:** 0.4 • 3 November 2026 (intersection linearisation for sensitivities moved to v0.2, ADR-0021 §9)

**Owner:** Onur Tuncer  
**Status:** Proposed release plan; acceptance evidence is required before tagging a release.

This plan embeds six modelling milestones derived from Ian Stroud's *Boundary Representation Modelling Techniques* into the existing release sequence, and a seventh, M7 machining services, that is Stereon's own requirement. It retains Phase 0 and the existing release anchors, adding explicit v0.2/v0.4/v0.5/v0.7/v0.9 gates between them. Milestone IDs M1–M7 below describe dependencies; they are not phase numbers, package versions, or the capability matrix's Mxx row IDs. R0–R5 remain long-term capability stages, and v1.0 does not imply comprehensive OCCT parity.

## Release mapping

| Release / implementation phase | Milestone contribution | Required outcome |
| --- | --- | --- |
| Phase 0 / 0.0.x | M1: executable modelling foundations, initial subset | Debug viewer, controlled builder prototype, separate validity checks, numerical policies, diagnostic dumps, replayable fixtures, and a minimal persistence prototype. Existing Phase 0 gates remain binding. |
| v0.1 Geometry / Phase 1 | M1 completion for available geometry; M4 geometry groundwork | Analytic and NURBS evaluation, fitting, projection and parameter-domain handling. Record tolerances and approximation limits; support save/load of the declared geometry subset. |
| v0.2 Intersections / Phase 2 | M2/M3 planar intersection services; M4 curved intersection groundwork | Inspectable intersection/classification results, including tangency, coincidence, overlap and unresolved/unsupported outcomes. Satisfy the applicable SSI gates before consuming these services in broad Booleans. Every intersection handler exposes the residual linearisation of its results and a regularity classification for sensitivities, with a first intersection JVP (ADR-0021 §9, stage 4). Tolerances originate in intersections, so their discipline is fixed here (ADR-0050, interim steps 0.2-T1 and 0.2-T2): typed growth causes, a ceiling relative to the model, the sum rule for pair tests, and intersection results that report the fit error and gaps that become tolerances. |
| v0.3 Solids / Phase 3 | M2: planar construction and interrogation; M4: curved and trimmed construction; M7: analytic machining queries | Planar faces with holes, polyhedral solids, extrusion, sectioning, imprinting and classification; then the advertised primitives, revolve, sweep and loft, trimming and sewing. Deliver operation history and native shape persistence for this supported subset. Certified curve discretisation and batched ray–face intersection for analytic faces; tri-dexel conversion of analytic solids (ADR-0044). Pcurves on every face coedge, `tighten`/`heal` and certified containment checks (ADR-0050, 0.3-T3); the face adjacency graph with tolerance-aware convexity and canonical surface classification, in an owned query workspace (ADR-0051, 0.3-F1 and 0.3-F2). |
| v0.4 Planar Booleans / Phase 4a | M3: planar Boolean pipeline | Union, intersection and difference for the declared planar/polyhedral domain, using inspectable splitting/imprinting, classification and assembly; validate lineage and attribute propagation. Face splitting reads ADR-0051's convexity together with ADR-0013's side evidence. |
| v0.5 Curved Booleans / Phase 4b | M4: bounded curved Boolean integration; M7: trimmed-spline machining queries | Extend the validated planar pipeline to a declared analytic and trimmed-spline subset; check seams, poles, tangencies, overlap and intersection certificates. No broad real-part robustness claim yet. Tri-dexel conversion and trim-aware closest-point queries extended to trimmed B-spline faces. |
| v0.6 Booleans / Phase 4c | M3/M4 robust engineering Boolean gate | Broaden the declared real-part domain and pass the existing Phase 4 corpus, success-rate, comparison, tolerance-growth and determinism gates. The tolerance-growth gate is ADR-0050's: output tolerances no worse than the largest input tolerance plus the SSI fit error, every growth with a typed cause. |
| v0.7 Meshing / Phase 5a | M4/M6 engineering workflow integration; M7: machining performance gate | Non-destructive display/analysis tessellation, mesh-to-B-rep associations, tagged duct workflow and usable Python construction/Boolean/mesh calls; validate native migrations. Public STEP/healing coverage remains the v0.8 gate. Tri-dexel fields with per-ray interval operations and tile-local updates, tool-contact and clearance queries, and benchmark B5.7 against OCCT. |
| v0.8 Exchange / Phase 5b | M4 import/repair integration; M6 initial product/application services | STEP AP242 read/write with explicit healing, derived display/analysis meshes, boundary-condition tags, native migrations, and usable Python workflows. Imported tolerances map onto ADR-0050 (cause `Import`, rejected or healed above `res_max`); an imported NURBS cylinder is recognised as a cylinder (ADR-0051). |
| v0.9 Finishing / Phase 6a | M5: local and finishing operations; M6 bounded application support; M7: tool offsets | Bounded constant-radius fillets, chamfers, shelling/thickening and shape sensitivities; history/tag integration and Phase 6 algorithm gates. Public API remains pre-1.0. Offset surfaces for tool-centre geometry with cusp, collapse and self-intersection detection, from M5's offsets. |
| v1.0 Stable kernel / Phase 6b | M1–M7 supported-domain release qualification | Stable documented public API, compatibility/migration notes, supported-platform qualification and end-to-end engineering acceptance. Re-run applicable release gates after final changes; do not expand to full OCCT parity by implication. |
| Post-v1.0 / R3–R5 workstreams as applicable | M5/M6 extended scope | Draft-angle operations, variable-radius blends, feature recognisers (separate libraries over ADR-0051's kernel queries) and feature verification, full parametric naming/history, richer assemblies/documents, drawing services, unfolding, medial-axis methods, swept-volume stock removal and manufacturing reasoning beyond M7, with separate users and acceptance datasets. |

Geometry and topology work may overlap. Curved geometry evaluation need not wait for planar Booleans. However, broad Boolean release claims depend on a validated planar pipeline, curved intersection/trim support, and a declared input domain. A bounded Phase 0 Boolean prototype is research evidence, not completion of v0.6.

## Interim steps

Interim steps break a release gate into reviewable stages with their own evidence. They are not releases and do not change the product version; each is recorded in the release's acceptance ledger.

| Step | Release | Content | Depends on | Exit evidence |
| --- | --- | --- | --- | --- |
| 0.2-T1 | v0.2 | ADR-0050: typed tolerance-growth causes on `SetTolerance`; `res_max = min(1e-3 × bounding diagonal, Context ceiling)` read by rule G04; a CI check that nothing else writes a tolerance; `coincident(a, b)` with the sum rule, and the Boolean prototype and Euler code moved to it | ADR-0001, ADR-0018 (existing fields, G01–G04) | Boolean prototype fixtures and ADR-0001's merge-policy test pass under the sum rule; the CI check rejects a seeded direct write |
| 0.2-T2 | v0.2 | ADR-0050: intersection results that create edges report the certified fit error and the gaps that become tolerances, with their causes | 0.2-T1; ADR-0010, ADR-0029 §6 | Every edge-creating intersection result in the corpus carries its tolerance evidence |
| 0.3-T3 | v0.3 | ADR-0050: pcurves on every face coedge; `tighten` and `heal`; certified G01/G02 bounds; the tolerant regression suite (gaps at 1e-6, 1e-4 and 1e-2 of the bounding diagonal) | 0.2-T2; ADR-0009 | The suite validates in CI; tightening lowers tolerances to the certified error without moving authoritative geometry |
| 0.3-F1 | v0.3 | ADR-0051: the `ShapeQueries` workspace keyed by snapshot and numerical context; tolerance-widened convexity; the face adjacency graph; loops and vertex convexity | 0.2-T1; ADR-0036, ADR-0037 | Synthetic-corpus edges classified as constructed, or `indeterminate` where an injected gap makes them so; no stale result in the edit fuzz test |
| 0.3-F2 | v0.3 | ADR-0051: canonical classification (analytic type check; freeform plane/sphere/cylinder/cone/torus fits with boundary verification); face roles; UV sampling; a reference recogniser in a separate library | 0.3-F1 | An in-kernel NURBS cylinder reports `Cylinder` within `res_abs`; the reference recogniser uses only the public API |

The v0.4 planar Boolean gate requires 0.3-F1; the v0.6 tolerance-growth gate and the v0.8 STEP gate close ADR-0050's verification (accepted 2026-11-10) and complete ADR-0051's acceptance.

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

Tolerant topology (ADR-0050) is M4's discipline for gaps: tolerances originate in v0.2's intersections (0.2-T1, 0.2-T2) and are tightened and certified in v0.3 (0.3-T3).

Exit evidence: topology and geometric boundaries agree across seams, poles, tangencies and gaps; intersection branches and coincidence regions are checked within the supported domain; certificates and cached representations obey ADR-0009/0010; sewing reports free edges, rejected identifications and tolerance changes. Boolean integration requires the applicable M3 gate and curved intersection/trim checks; successful evaluator tests alone do not establish curved Boolean support.

### M5 — Local and finishing operations

Deliver bounded constant-radius fillets, chamfers, offsets and shelling/thickening for v0.9, qualified for stable supported-domain use at v1.0. Draft-angle operations and variable-radius blends remain post-v1.0 scope.

Exit evidence: supported configurations and limits are documented; collapse, self-intersection and topology-changing events are detected; results validate or return structured errors; lineage and tags survive edits. Shape sensitivities distinguish fixed-topology derivatives from nondifferentiable events per ADR-0021. Retain the Phase 6 benchmark gates B6.1–B6.4.

### M6 — Product and application services

Begin lineage and attribute infrastructure with v0.3; exercise it in v0.4–v0.6 and publish bounded mesh/Python workflows in v0.7, followed by exchange/healing in v0.8. Validate finishing-operation integration and supported interface stability for v1.0. The kernel provides the queries that recognisers need (ADR-0051: face adjacency, convexity, canonical surfaces, face roles, UV samples) from v0.3; recognisers themselves, assemblies, document services and manufacturing reasoning remain later capabilities above the kernel.

Exit evidence: split/merge lineage composes through chained operations and compaction; metadata conflicts and invalidated references are visible; the tagged duct workflow preserves inlet/outlet/wall semantics through cut, validation, watertight mesh and external volume meshing (ADR-0011/0014). Persistent naming across regeneration remains a separate post-v1.0 problem; operation history is not a claim that it is solved.

### M7 — Machining services

Stereon's own requirement, not derived from Stroud. Provide the kernel queries a machining application needs: certified discretisation, ray and distance queries against trimmed faces, B-rep to tri-dexel conversion for in-process stock (ADR-0044), tool-contact and clearance queries, and tool offsets. Toolpath strategy, cutter libraries, stock-update policy and postprocessing stay in applications above the kernel. Display sampling is never machining discretisation: `SampleCurve` and tessellation bounds are for display only.

Stage the work across releases. v0.3: certified curve subdivision with conservative deviation bounds, batched ray–trimmed-face intersection for analytic faces, and tri-dexel conversion of analytic solids. v0.5: the same for trimmed B-spline faces, and trim-aware closest-point queries. v0.7: the tri-dexel field type with per-ray interval operations and tile-local updates, tool-contact and clearance queries, and the B5.7 performance gate. v0.9: offset surfaces for tool-centre geometry, built on M5's offsets. v1.0: qualification of the supported domain. Ray–face intersection needs no general SSI and may run ahead of it; offsets and contact queries on freeform faces depend on M4's curve/surface intersection services.

Exit evidence: every discretisation and every tri-dexel endpoint carries a certified error bound, kept apart from grid resolution; features finer than the grid are reported, and unresolved rays or queries return structured errors rather than missing material; fields are bit-identical across schedulers and thread counts (ADR-0007); endpoint derivatives on fixed face branches agree with finite differences and events are reported (ADR-0021); benchmark B5.7 is measured on identical rays, trims and accuracy requirements against OCCT, with memory and failure rates. A speed claim against OCCT is made per benchmark, never as a general one.

## Architectural dependencies

The book motivates these acceptance requirements, not an automatic change of ADR status. Check and review the relevant records before implementing their gates.

| Concern | Existing record / planning action |
| --- | --- |
| Numerical policy and purpose-specific tolerances | ADR-0001/0002 and [ADR-0028](../adr/0028-purpose-specific-numerical-tolerances.md): purpose-specific controls, dimensional conversions, observable changes and persisted numerical context; policy in M1/v0.1, intersection use in v0.2, native shape context by v0.3. |
| Intersection/classification outcomes | ADR-0010/0013 and [ADR-0029](../adr/0029-intersection-and-classification-outcomes.md): dimension-aware results and independent side-transition evidence; v0.2 intersection gate, v0.3 queries and v0.4–v0.6 Booleans. RS-0001 remains research. |
| Intersection sensitivities | ADR-0021 §6 and §9: residuals, tracking constraints, `F_z`/`F_p` and a regularity classification in every intersection handler from v0.2, before the generic SSI is built; projections and moving trims after; shape sensitivities by v0.9. |
| Entity lineage | ADR-0011 and [ADR-0030](../adr/0030-entity-lineage-in-operation-results.md): qualified many-to-many relations, explicit split/merge events and composition by v0.3, extended through v0.4–v0.6. |
| Attribute propagation | [ADR-0031](../adr/0031-attribute-propagation-and-conflicts.md): per-key split/merge/revalidation rules and conflict diagnostics in v0.3–v0.6; tagged engineering mesh gate in v0.7. |
| Coincidence versus topology identification | ADR-0004/0009/0018 plus construction/sewing contracts; record identification rules before expanding sewing or simplification. |
| Tolerant topology | [ADR-0050](../adr/0050-tolerant-topology.md) with ADR-0001/0009/0018/0029: typed growth causes, model-relative ceiling, sum rule for pair tests, tightening; 0.2-T1/T2 in v0.2, 0.3-T3 in v0.3, growth gate in v0.6, import in v0.8. On acceptance it supersedes ADR-0001's pair rule and ceiling. |
| Face adjacency and canonical surfaces | [ADR-0051](../adr/0051-face-adjacency-and-canonical-surfaces.md) with ADR-0013/0036/0037/0050: kernel queries in an owned workspace; convexity is never classification evidence; recognisers stay above the kernel; 0.3-F1/F2 in v0.3, consumed by the v0.4 Booleans. |
| Native schema and migrations | ADR-0012; prototype early, public shape persistence by v0.3. |
| Derived tessellation and tags | ADR-0014; mesh generation must not alter authoritative B-rep geometry/topology. |
| Machining queries and tri-dexel fields | [ADR-0044](../adr/0044-brep-to-tri-dexel-conversion.md) with ADR-0021/0028/0029: certified endpoints, explicit unresolved rays, deterministic tiles; analytic subset by v0.3, trimmed splines by v0.5, B5.7 gate by v0.7. |
| Reproducible failures | ADR-0006/0023; journals, dumps, numerical context and stage diagnostics accompany corpus failures. |
| Controlled edits and operation contracts | ADR-0025/0026/0027 alongside ADR-0004/0005/0018 define association, operation and controlled-edit contracts; review before phase acceptance. |

## Versioning and release gates

1. **Product version.** `CMakeLists.txt` is the source for the numeric kernel version (`LibraryVersion`, documentation and bindings follow it). The current version is 0.1.0, the v0.1 "Geometry" milestone: Phase 0's exit gate met and every gate of the [v0.1 acceptance ledger](v0.1-acceptance-ledger.md) done, with the evidence, limits and migration notes in the [release report](v0.1-release-report.md). Phase 0 uses 0.0.x; the advertised releases target every minor from 0.1.0 through 0.9.0, followed by 1.0.0. Each minor has the bounded scope and acceptance gate in the release table. Phase 4 is split across 0.4–0.6, Phase 5 across 0.7–0.8 and Phase 6 across 0.9–1.0; implementation phases need not match minor numbers.
2. **Pre-1.0 changes.** New capabilities or incompatible public API changes after a tagged milestone require a minor-version increment with migration notes. Patch releases contain compatible fixes and documentation updates; 0.0.x remains explicitly experimental. Candidate tags use `vX.Y.Z-rc.N` and do not assert that the final gate has passed.
3. **Release evidence.** Tag a final milestone only after all mandatory supported-domain acceptance fixtures, applicable ADR gates, supported-platform CI, determinism checks and required benchmark gates pass. Publish the exact implementation commit, corpus revision, comparator build, domain limits, structured failures and observed invalid-success count in a release report. A benchmark target is not an achieved measurement. M1–M7 evidence complements existing gates; it does not weaken them.
4. **Post-1.0 API policy.** Within the declared stable public API, additive compatible capabilities increment the minor version, compatible fixes increment the patch, and incompatible changes increment the major. C++ source compatibility and binary compatibility are distinct: compiler/runtime/configuration requirements and rebuild rules follow ADR-0019. Stable API does not promise cross-compiler C++ ABI compatibility.
5. **Native format version.** `.stn` schema versions are independent of product versions. Follow ADR-0012: every schema change increments the format version; explicit migrations read older public formats; unknown required chunks fail and unknown optional chunks warn. Keep fixtures from each tagged public-format release. Prototype dump versions must be labelled separately from public `.stn` schema versions.
6. **C ABI version.** Version/capability queries and ABI-major compatibility follow ADR-0020 independently of product versions. Do not advertise a stable C ABI until its ownership, error and old-consumer compatibility gates pass. C ABI implementation is not a prerequisite for early geometry releases.
7. **Planning versions.** This document's revision and the capability matrix's document version identify planning edits, not delivered software. Phase numbers, M1–M7 milestone IDs and R0–R5 capability stages must never be presented as package versions or completion percentages.

Use an acceptance ledger per release: milestone/subgate, supported domain, owner, required ADR decision, fixture IDs, implementation commit, result/report link and unresolved limitations. All gates above are planned requirements; no new capability completion or audit grade is asserted by this document.

## Source and related plans

Stroud (2006), *Boundary Representation Modelling Techniques*: chapters 3–7 for representation and stepwise operations; §8.3 for information preservation; §10.3 for faceting; §§11.1–11.3 for persistence; §13.5 for sewing; §§16.1–16.3 for tolerances and diagnostics. Modern staging, immutable mesh derivation and versioning details here are Stereon planning choices, not quotations or claims that the book specifies them.

- [Capability matrix and historical audit](stereon-capability-matrix.md)
- [Phase 0 implementation checklist](../../TODO.md)
- [Benchmark specifications and existing phase gates](../benchmarks.md)
- [Architecture decision records](../adr/README.md)
