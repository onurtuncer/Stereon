# Stereon capability matrix
Version 0.2 • 12 October 2026 • Source and CI audit added

**Status:** Proposed planning baseline (not an ADR).
**Owner:** Onur Tuncer.
**Review date:** 2026-10-12.
**Review trigger:** topology prototype completion or next phase gate.

See the [planning index](README.md) for scope and follow-up work.

## Publication update — 13 October 2026

The detailed audit below describes commit `cf1dbaa925fb9c60fc642ca6962f329af7d52047` on 12 October; it is historical evidence, not a claim about today's main branch.

This document is being added on top of [25effbdcf26d1f35e3e6ba3314bb54512b43b5bd](https://github.com/onurtuncer/Stereon/commit/25effbdcf26d1f35e3e6ba3314bb54512b43b5bd). Since the audit, [PR #27](https://github.com/onurtuncer/Stereon/pull/27) changed arena mutation to scoped callbacks and included the surface-trigonometry module-loading fix; [PR #28](https://github.com/onurtuncer/Stereon/pull/28) added predicate numerical-domain handling and tests. These changes deserve a future targeted reassessment of F01/F02 and their evidence. The 54 row grades have not been re-audited here.

All nine workflow runs returned success for the publication base, including [Linux](https://github.com/onurtuncer/Stereon/actions/runs/38391298134). Accordingly, the failed Release/TSan builds and recommendation to restore CI in the historical appendix are no longer current blockers at that base. No new OCCT comparison evidence is claimed.

## Objective and scope
Build an independent industrial CAD platform that can rival the open-source OCCT platform across its functional domains and demonstrated reliability. Comprehensive coverage is the destination; releases establish useful subsets along the way.

This matrix now includes a source and CI audit of [onurtuncer/Stereon at cf1dbaa925fb](https://github.com/onurtuncer/Stereon/tree/cf1dbaa925fb9c60fc642ca6962f329af7d52047). Grades describe only the explicitly identified implemented subset, never completion of an entire capability family. Proposed acceptance criteria and stages remain planning recommendations, not delivery commitments. Detailed evidence and gaps appear in the audit appendix.

Reference: OCCT documentation identifies version 8.0.1 at the time of preparation. Before comparative testing, pin its exact source tag/commit, dependencies, build flags, platform, and enabled modules. Documentation coverage is not evidence that OCCT succeeds on every case.

Scope classes:
- **B — baseline:** a functional family represented in the open-source OCCT platform; detailed subfeature coverage must be inventoried before a parity claim.
- **D — differentiation:** an additional Stereon objective; this label does not assert that competitors lack it.
- **E — ecosystem extension:** separately tracked beyond the open-source baseline. Commercial SDKs and third-party applications must not be silently counted as OCCT core.

Functional equivalence does not require the same classes, internal representations, or algorithms. External components may satisfy a capability if their integration, redistribution terms, maintenance, and operational dependencies are documented. An OCCT-backed adapter may aid migration, but does not count as independent Stereon implementation of the delegated operation.

## Maturity and evidence
Track two dimensions per row: capability maturity and comparative evidence.

| Code | Capability maturity |
|---|---|
| U | Unassessed: no repository evidence reviewed |
| 0 | Assessed and absent |
| 1 | Experimental implementation with stated limitations |
| 2 | Defined subset supported with automated acceptance cases and documented errors |
| 3 | Production candidate: public contract, regression coverage, independent checks, supported-platform CI, and known limitations |
| 4 | Industrial evidence: external application use and sustained maintenance on representative workloads |

| Code | Comparative evidence |
|---|---|
| — | No comparison performed |
| C | Comparison measured and published; parity not established |
| P | Parity gate passed for an explicitly bounded subfeature and corpus |
| A | Demonstrated advantage for a named metric, without hiding regressions |

All 54 rows are source-assessed: **6 at maturity 2, 10 at maturity 1, and 38 at maturity 0**. No row is graded 3 or 4. Every row remains **— for comparative evidence**, because no executed Stereon-versus-OCCT comparison was found. CODEOWNERS names @onurtuncer as current repository owner; future delivery owners remain to be assigned. These counts are an inventory, not a completion percentage. A grade 2 is restricted to the subset named in the appendix. Split broad rows before scheduling work.

## Release stages
Stages below indicate first useful support, not completion of the entire row.

| Stage | Intended usable outcome | Exit evidence |
|---|---|---|
| R0 — contracts and measurement | Foundations on which later modules can depend | Agreed validity, tolerance, ownership, error and history contracts; versioned fixtures and comparator harness |
| R1 — geometry and representation | Construct, interrogate, validate, save and inspect defined geometry subsets | All declared R1 acceptance fixtures pass; supported-domain limits documented |
| R2 — useful solid modeling | Independently build and exchange simple engineering parts | Repeatable construction → modification → export → independent read workflow |
| R3 — engineering breadth | Freeform parts, advanced modifications and repair | Representative difficult cases plus application pilot; no silent invalid successes on gate corpus |
| R4 — application platform | Build an engineering application around Stereon | Assembly/document/viewer integration, lifecycle tests, stable supported interfaces |
| R5 — industrial parity | Demonstrated breadth against the frozen OCCT baseline | All baseline subfeatures mapped and assessed; comparative gates met and limitations disclosed |

Workstreams overlap. For example, exchange and healing prototypes should expose representation problems before R3. Visualization needs an early diagnostic viewer even though its full platform maturity comes later.

## A. Foundations and geometric representation
Source family mapping: [S1]. Acceptance requirements are proposed Stereon contracts.

| ID | Capability / class | First stage | Required acceptance evidence | Maturity / comparison |
|---|---|---|---|---|
| F01 | Numerical policy, units and tolerances / B | R0 | Explicit linear, angular and parametric tolerances; equivalent results after unit conversion; scale and conditioning limits documented | 2 / — |
| F02 | Ownership, allocation and object lifetime / B | R0 | Copy/share/mutation semantics tested; invalid handles detected; allocation failures preserve valid state | 2 / — |
| F03 | Errors, progress and cancellation / B | R0 | Typed outcomes identify failing operation and entities; cancellation leaves no partially committed model | 1 / — |
| F04 | Transforms, placements and spatial bounds / B | R1 | Composition/inverse identities; orientation behavior under reflection; bounds conservatively contain tested geometry | 2 / — |
| G01 | Analytic 2D/3D curves and surfaces / B | R1 | Values, derivatives, periods and domains checked against independent formulas, including seams and poles | 2 / — |
| G02 | Bézier, B-spline and rational geometry / B | R1 | Knot insertion and degree elevation preserve geometry within budget; repeated knots and extreme weights covered | 0 / — |
| G03 | Trimming, periodicity and curve-on-surface representation / B | R1 | 3D edge and surface parameter curve agree within declared tolerance; seam crossings and degenerate boundaries validated | 1 / — |
| G04 | Interpolation, approximation and conversion / B | R2 | Residuals and requested continuity measured; unmet fitting budgets return explicit failure or qualified output | 0 / — |

## B. Topology and geometric algorithms
Source family mapping: [S1, S2].

| ID | Capability / class | First stage | Required acceptance evidence | Maturity / comparison |
|---|---|---|---|---|
| T01 | Boundary representation and incidence / B | R1 | Vertices through solids and mixed collections; orientation, shared entities and permitted non-manifold configurations have explicit invariants | 0 / — |
| T02 | Validation and point classification / B | R1 | Seeded defects detected; inside/outside/on-boundary classification checked on analytic cases and ambiguous neighborhoods | 0 / — |
| T03 | Construction, copying and traversal / B | R1 | Shared topology preserved when intended; deep copies independent; traversals terminate and respect occurrence orientation | 0 / — |
| T04 | Modification history / B | R2 | Generated, modified and deleted entity relations correct through chained operations, including splits and merges | 0 / — |
| Q01 | Projection and closest-point queries / B | R1 | Residual and distance checks; multiple solutions and singular points handled; global versus local guarantees stated | 2 / — |
| Q02 | Intersections / B | R2 | Curve/curve, curve/surface and surface/surface coverage enumerated; tangency, coincidence and small components included; returned geometry validated | 0 / — |
| Q03 | Extrema, distance and spatial acceleration / B | R2 | Analytic minima and independent references; acceleration produces no false exclusions relative to exhaustive tests | 1 / — |
| Q04 | Length, area, volume, centroid and inertia / B | R1 | Analytic references, rigid-transform identities and numerical convergence; invalid domains reported explicitly | 0 / — |

## C. Construction and modification
Source family mapping: [S2]. Each operation must declare allowed input types and output semantics.

| ID | Capability / class | First stage | Required acceptance evidence | Maturity / comparison |
|---|---|---|---|---|
| M01 | Primitives, extrusions and revolutions / B | R2 | Closed valid results where promised; analytic mass properties; partial angles, reversed profiles and degenerate inputs covered | 0 / — |
| M02 | General sweeps, pipes and lofts / B | R3 | Frame/twist policy documented; section correspondence verified; continuity and self-intersection outcomes measured | 0 / — |
| M03 | Union, intersection, difference, section and splitting / B | R2 | Set membership and volume identities, valid boundaries, history and perturbation tests; touching and coincident cases explicit | 0 / — |
| M04 | Sewing and shell/solid assembly / B | R2 | Gaps assessed against budget; orientations consistent; open shells never silently promoted to closed solids | 0 / — |
| M05 | Offsets, thickening and hollowing / B | R3 | Offset distance checked away from singularities; feature collapse and self-intersection produce defined outcomes | 0 / — |
| M06 | Fillets and chamfers / B | R3 | Requested radius/chord policy and continuity checked; corner interactions, short edges and variable laws represented in corpus | 0 / — |
| M07 | Drafts and local feature construction / B | R3 | Geometric angle/depth constraints verified; adjacent boundaries reconstructed; failed edits roll back | 0 / — |
| M08 | Surface construction and filling / B | R3 | Boundary fit and requested continuity verified; multi-patch and incompatible constraints tested | 0 / — |
| M09 | Simplification and local feature removal / B | R3 | Removed feature scope declared; retained regions checked against displacement budget; topology and history validated | 0 / — |

## D. Diagnosis and repair
Source family mapping: [S3]. Repair changes must be observable.

| ID | Capability / class | First stage | Required acceptance evidence | Maturity / comparison |
|---|---|---|---|---|
| H01 | Shape diagnosis / B | R1 | Localized reports for gaps, missing parameter curves, orientation errors and degeneracies; clean fixtures avoid spurious errors | 0 / — |
| H02 | Wire, edge, face and shell repair / B | R2 | Before/after defect report; admissible displacement and tolerance inflation enforced; unrepairable cases remain explicit | 0 / — |
| H03 | Representation regularization / B | R3 | Splitting, merging and continuity adjustments preserve intended geometry within budget; all entity replacements recorded | 0 / — |
| H04 | Controlled repair pipelines / B | R3 | Configurations replayable; each repair step logged; no silent feature deletion or unit rescaling | 0 / — |

## E. Exchange and persistence
Source family mapping: [S1, S4]. Format names are not sufficient coverage claims: track read/write direction, schema, entity class and attributes.

| ID | Capability / class | First stage | Required acceptance evidence | Maturity / comparison |
|---|---|---|---|---|
| X01 | Native model persistence / B | R1 | Geometry, topology and sharing round-trip; corrupt input rejected; version migration policy tested | 0 / — |
| X02 | STEP geometry read/write / B | R2 | Exact supported schema/entity inventory; cross-reader checks of dimensions, units, geometry and solid validity | 0 / — |
| X03 | STEP assemblies and engineering metadata / B | R4 | Placements, shared instances, names, colors and supported PMI semantics compared field by field; losses reported | 0 / — |
| X04 | IGES read/write / B | R3 | Supported entity inventory; trimmed surfaces and units checked; unsupported constructs diagnosed | 0 / — |
| X05 | Mesh and visualization exchange / B | R2 | Per-format coverage inventory, initially STL; orientation, scale and attributes checked where representable; additional baseline formats completed by R5 | 0 / — |
| X06 | OCCT migration bridge / D | R2 | Bounded model conversion plus application porting example; provenance, conversion loss and runtime dependency disclosed | 0 / — |

## F. Meshing, visualization and application services
Source family mapping: [S1, S5]. Broad parity requires subfeature inventories within these rows.

| ID | Capability / class | First stage | Required acceptance evidence | Maturity / comparison |
|---|---|---|---|---|
| V01 | Boundary tessellation and mesh association / B | R1 | Stated chordal/angular targets measured; shared edges consistent; orientation and face association retained | 0 / — |
| V02 | CAD inspection and interactive selection / B | R1 | Wireframe/shaded inspection and subshape picking first; later instancing, clipping, transparency, materials and large-scene budgets | 0 / — |
| V03 | Drawing and annotation services / B | R4 | Hidden-line, section and dimension fixtures; visibility and placement checked on analytic scenes | 0 / — |
| A01 | Documents, attributes and transactions / B | R4 | Undo/redo, copy/paste and rollback recover geometry and metadata; transaction isolation contract documented | 0 / — |
| A02 | Persistent references and regeneration support / B | R3 | References survive defined edits; ambiguous or deleted targets reported; dependency update and failure propagation tested | 0 / — |
| A03 | Assemblies and document exchange integration / B | R4 | Shared prototypes, nested placements and metadata survive editing and persistence; large-assembly limits measured | 0 / — |
| A04 | Scripting and diagnostic harness / B | R1 | Reproducible headless operation scripts; fixtures emit machine-readable results; all advertised examples run in CI | 1 / — |

## G. Product quality and deliberate advantages
These are Stereon product requirements rather than claims about missing OCCT functionality.

| ID | Capability / class | First stage | Required acceptance evidence | Maturity / comparison |
|---|---|---|---|---|
| P01 | SDK packaging, interfaces and bindings / D | R1 | Clean consumer builds; stable C boundary and Python binding considered against existing ADRs; ownership and exception behavior tested | 1 / — |
| P02 | Concurrency and reproducibility / D | R0 | Shared-read/concurrent-edit rules explicit; race checks; repeated runs semantically equivalent; any bitwise guarantee separately bounded | 1 / — |
| P03 | Performance and memory discipline / D | R1 | Paired workloads report median/p95 time and peak memory; accuracy held constant; initialization and I/O separated | 1 / — |
| P04 | CI, fuzzing and regression maintenance / D | R0 | Sanitizers, malformed inputs, minimized reproducers, platform matrix and tracked failures; corpus licensing documented | 2 / — |
| P05 | Documentation, releases and support / D | R1 | Executable tutorials, migration notes, compatibility policy, issue triage and named maintainers before production claims | 1 / — |
| D01 | Shape-derivative contracts / D | R0 | Scalar/AD compatibility decided early; later gradient checks on smooth fixed-topology cases; topology events and undefined derivatives explicitly signaled | 1 / — |
| D02 | Operation provenance and automation / D | R0 | Inputs, configuration and entity lineage sufficient for replay; machine-readable diagnostics; bounded recovery avoids silent geometry changes | 1 / — |
| D03 | Simulation integration / D | R3 | CAD parameters → geometry → mesh association → solver sensitivities demonstrated; full-loop gradients validated on a named engineering problem | 0 / — |

## H. Separately tracked ecosystem scope
These remain expansion objectives. Their separate classification does not reduce the full-platform ambition.

| ID | Capability / class | Initial planning position | Evidence required before commitment | Maturity / comparison |
|---|---|---|---|---|
| E01 | Proprietary CAD translators / E | Evaluate after neutral exchange | Actual customer formats, access/licensing route, read/write coverage and external round-trip tests | 0 / — |
| E02 | Advanced volume meshing and specialist mesh processing / E | Integration study alongside R3 | Solver requirements, element quality and boundary fidelity; ownership of external dependencies | 0 / — |
| E03 | Full sketch constraints, feature authoring and assembly solving / E | Application-layer roadmap | Defined user workflows and solver integration; distinguish these from low-level geometric construction constraints | 0 / — |
| E04 | Specialized CAM, recognition and reconstruction tools / E | Demand-led modules | Named use cases and standalone acceptance corpora; separate product support commitments | 0 / — |

## Comparative acceptance protocol
The following are proposed planning gates, not measured results or established industry standards.

1. **Freeze the scope before running.** Record supported subfeatures, model families, input validity rules, tolerances, unit conventions, platform and resource limits. Unsupported cases stay visible in coverage reports.
2. **Use several independent corpora.** Include analytic constructions, difficult synthetic degeneracies, licensed real engineering models, and an external holdout set. Partition by model family, not merely by near-duplicate file.
3. **Judge semantics independently.** Use analytic references, validated high-precision calculations where suitable, and independent geometric checks. Treat OCCT as a comparator, not the correctness oracle. Entity numbering and face counts need not match.
4. **Separate outcomes.** Count verified success, explicit rejection, unsupported input, invalid-success, timeout and crash. Rejection is safer than false success, but is still a failure to fulfill an in-scope operation.
5. **Measure quality.** Check boundary validity, geometry discrepancy, topology requirements, conservation identities and attribute preservation. Approximate sampling is useful evidence, not a proof of global error bounds.
6. **Measure cost on equivalent work.** Hold required accuracy and input preparation constant; disclose healing and caching. Report latency distributions and memory on the same hardware.
7. **Publish adverse cases.** Include Stereon-only failures, OCCT-only failures, both-fail cases and known ambiguity. Retain regressions in the corpus.
8. **Have another party reproduce the result.** A maintainer-written benchmark alone is insufficient evidence of industrial maturity.

Provisional per-subfeature parity gate for R5:
- 100% of deterministic contract/conformance fixtures pass.
- Zero crashes and zero silent invalid-success outcomes observed on the declared release corpus. This is a finite-test result, not a universal guarantee.
- Every baseline subfeature is mapped; omissions preclude an unqualified full-coverage claim.
- Verified-success rate is no lower than the pinned baseline on each predeclared critical family; publish counts and paired statistical uncertainty.
- Critical customer workflows pass individually; aggregate rates cannot conceal their failure.
- Runtime and memory fit declared workload budgets. As an initial investigation threshold, review median or p95 runtime above 2× the comparator and peak memory above 2×; these are proposed review triggers, not automatic acceptance or claims of performance parity.
- External holdout results and limitations accompany the release report.

Set corpus sizes and statistical margins after a pilot establishes case diversity and baseline failure rates. A large number of near-duplicate boxes is not broad evidence.

## Dependency order and early decisions
The main modeling dependency is numerical policy → geometry/parameterization → topology/validation → intersections/classification → Booleans → advanced modification and repair. Exchange and meshing should exercise each layer as soon as useful subsets exist. Document identity and operation history must be designed before dependent applications become difficult to migrate.

Before R1, resolve or confirm existing ADRs for:
- Error and tolerance propagation; approximate versus exact/adaptive methods.
- Object identity, ownership, mutability and concurrency.
- Parameter-curve consistency, degeneracy and non-manifold representation.
- Modification history and persistent-reference ambiguity.
- Serialization evolution, public interfaces and dependency boundaries.
- Differentiation scope and behavior at nonsmooth/topology-changing events.

Do not replace existing ADR decisions on the strength of this planning document. Reconcile this matrix against the repository before turning recommendations into requirements.

## First execution package
1. Source/CI audit completed in this version; use the evidence appendix to create subfeature issues. Refresh against a new pinned commit after meaningful implementation changes.
2. Pin the OCCT comparison build and create a subfeature inventory for each B row.
3. Establish the harness and a small diagnostic seed corpus spanning analytic, periodic, trimmed, degenerate and imported geometry.
4. Build the first vertical workflow: construct a part → modify → validate → compute properties → tessellate → STEP export → independent read.
5. Assign owners and identify the numerical/representation risks that block R2. Estimate dates only after this audit.

For each row, maintain: ID; subfeature; class; owner; maturity; implementation commit; acceptance fixtures; valid-input domain; tolerance contract; benchmark baseline; comparative outcome; limitations; dependencies; target release. Changes to test scope and accepted limitations require a visible review trail.

## Sources and interpretation
The following official documentation establishes the reference families. The proposed architecture, stages, acceptance tests and thresholds are original planning recommendations, not statements made by these sources. Links may track newer documentation later; archive the chosen baseline when benchmarking.

- **S1:** OCCT introduction and platform overview: https://occt3d.com/dev/doc/overview/html/
- **S2:** Modeling Algorithms user guide: https://occt3d.com/dev/doc/overview/html/occt_user_guides__modeling_algos.html
- **S3:** Shape Healing user guide: https://occt3d.com/dev/doc/overview/html/occt_user_guides__shape_healing.html
- **S4:** STEP Translator user guide: https://occt3d.com/dev/doc/overview/html/occt_user_guides__step.html
- **S5:** OCAF user guide: https://occt3d.com/dev/doc/overview/html/occt_user_guides__ocaf.html

This matrix includes Stereon source/CI evidence. Full API/subfeature completeness against OCCT source has not yet been audited, and no OCCT parity is claimed.


## Repository audit — evidence and limits

**Snapshot:** [cf1dbaa925fb9c60fc642ca6962f329af7d52047](https://github.com/onurtuncer/Stereon/commit/cf1dbaa925fb9c60fc642ca6962f329af7d52047), the observed main-branch head when this audit began. The repository was not modified. Evidence links below are commit-pinned and do not silently follow later changes.

**Method:** inventoried the tracked tree; inspected implementations, test registrations and relevant test cases, module stubs, roadmap/ADRs, CI workflows and the exact commit's workflow/job results. This is a capability inventory and targeted source review, not an exhaustive numerical correctness proof, security audit or coverage certification.

**Local verification:** layering passed for 11 libraries; epsilon check passed on 3 files (currently mostly stubs); literal checks passed on 71 files; corpus validation passed on 28 cases; the corpus self-test passed 21 surface pairs in both operand orders. These are source/reference-data checks. No local C++ build was performed: this environment has GCC 13.3, no GCC 16 command, and no CMake. C++ execution evidence comes from the repository's CI.

**Inventory:** 18 test-source files contain 107 TEST_CASE declarations. This is not the number of executed CTest tests or a coverage percentage. There are 28 corpus JSON fixtures and 24 numbered ADRs excluding the template. All 24 ADRs are marked Proposed, including some whose ideas are already partly implemented.

### CI at the audited commit

| Workflow / job | Observed result | Interpretation |
| --- | --- | --- |
| [Linux Debug and ASan](https://github.com/onurtuncer/Stereon/actions/runs/38261612048) | Build and Test steps succeeded | Positive execution evidence for existing test suites, not platform parity |
| [Linux Release](https://github.com/onurtuncer/Stereon/actions/runs/38261612048/job/114853245787) | Build failed; Test skipped | SurfaceTests.cpp compilation fails while loading a compiled module |
| [Linux TSan](https://github.com/onurtuncer/Stereon/actions/runs/38261612048/job/114853246107) | Build failed before tests | No passing TSan execution evidence for this commit |
| [Windows](https://github.com/onurtuncer/Stereon/actions/runs/38261612112) and [macOS](https://github.com/onurtuncer/Stereon/actions/runs/38261611898) | Workflows succeeded | GCC platform coverage recorded; not a guarantee for every supported configuration |
| [Clang and MSVC](https://github.com/onurtuncer/Stereon/actions/runs/38261611866) | Both jobs' Configure, Build and Test steps succeeded | Verified individual steps, not merely the non-blocking workflow badge; policy still treats these as non-reference compilers |
| [Python](https://github.com/onurtuncer/Stereon/actions/runs/38261611876) | Succeeded | Current binding is only a version attribute, not geometry functionality |
| [Benchmarks](https://github.com/onurtuncer/Stereon/actions/runs/38261612029) | Succeeded | Workflow builds and smoke-runs microbenchmarks; does not run an OCCT comparator |
| [Lint](https://github.com/onurtuncer/Stereon/actions/runs/38261611957), [Coverage](https://github.com/onurtuncer/Stereon/actions/runs/38261612059), [Docs](https://github.com/onurtuncer/Stereon/actions/runs/38261611930) | Succeeded | Infrastructure evidence; coverage percentage was not retrieved |

Release and TSan logs both report `Stereon.Surface: error: failed to read compiled module cluster 79: Bad file data`, followed by a failure to load `::sin@Stereon.Surface` at `tests/support/Dual.hpp:65`. This identifies a build/module-loading blocker. It is not evidence of a failed geometric assertion, and the compiler/root cause has not been reproduced or proven in this audit.

The contract probes in [tests/CMakeLists.txt](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/tests/CMakeLists.txt) explicitly encode a known failure for imported module-template contracts. Passing CI must therefore not be interpreted as universal contract enforcement. See the [toolchain probe record](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/docs/toolchain/gcc16-feature-probe.md).

### Capability evidence register

Maturity 0 means the named functional capability is absent from the inspected production source. A proposed ADR, empty target, planned example or supporting primitive does not raise that grade. Existing supporting pieces are identified where useful. Grade 1 means partial groundwork; grade 2 means a defined, tested subset. No grade asserts production readiness. Current repository ownership follows CODEOWNERS; no new assignees have been appointed.

#### F01 — maturity 2; comparison —

Context tolerance controls, outward-rounded intervals, and adaptive predicates are implemented and tested.

- Implementation/design evidence: [src/core/Context.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/core/Context.cppm); [src/core/Interval.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/core/Interval.cppm); [src/robust/Predicates.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/robust/Predicates.cppm).
- Test evidence: [tests/core/ContextTests.cpp](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/tests/core/ContextTests.cpp); [tests/core/IntervalTests.cpp](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/tests/core/IntervalTests.cpp); [tests/robust/PredicateTests.cpp](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/tests/robust/PredicateTests.cpp).
- Remaining gap: Unit conversion, model-scale validation and per-entity tolerance propagation through actual operations remain unproved.

#### F02 — maturity 2; comparison —

Typed/generational handles, checked access and paged copy-on-write storage provide a tested foundation.

- Implementation/design evidence: [src/core/Handle.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/core/Handle.cppm); [src/core/Arena.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/core/Arena.cppm); [src/core/SmallBuffer.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/core/SmallBuffer.cppm).
- Test evidence: [tests/core/HandleTests.cpp](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/tests/core/HandleTests.cpp); [tests/core/ArenaTests.cpp](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/tests/core/ArenaTests.cpp); [tests/curve/CurveTests.cpp](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/tests/curve/CurveTests.cpp).
- Remaining gap: No Shape/ShapeBuilder ownership lifecycle yet; allocation-failure injection and application-level lifetime gates remain.

#### F03 — maturity 1; comparison —

KernelError, expected-based propagation, contracts and journal identifiers exist.

- Implementation/design evidence: [src/core/Error.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/core/Error.cppm); [src/core/include/Stereon/Core/Try.hpp](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/core/include/Stereon/Core/Try.hpp); [src/core/Journal.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/core/Journal.cppm).
- Test evidence: [tests/core/ErrorTests.cpp](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/tests/core/ErrorTests.cpp); [tests/core/ContractTests.cpp](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/tests/core/ContractTests.cpp); [tests/probes/ContractProbe.cpp](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/tests/probes/ContractProbe.cpp).
- Remaining gap: No model-operation cancellation/progress/rollback implementation; module contract enforcement has documented compiler caveats.

#### F04 — maturity 2; comparison —

Scalar-generic matrices/transforms and interval-based boxes have unit and property tests.

- Implementation/design evidence: [src/core/Matrix.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/core/Matrix.cppm); [src/core/Box.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/core/Box.cppm).
- Test evidence: [tests/core/MatrixTests.cpp](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/tests/core/MatrixTests.cpp); [tests/core/BoxTests.cpp](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/tests/core/BoxTests.cpp).
- Remaining gap: No B-rep placements; built-in curved geometry transforms assume similarity transforms, not general affine closure.

#### G01 — maturity 2; comparison —

Line/Circle in 2D and 3D and Plane/Cylinder surfaces support evaluation, derivatives, projection and reversal.

- Implementation/design evidence: [src/curve/Line.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/curve/Line.cppm); [src/curve/Circle.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/curve/Circle.cppm); [src/surface/Plane.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/surface/Plane.cppm); [src/surface/Cylinder.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/surface/Cylinder.cppm).
- Test evidence: [tests/curve/CurveTests.cpp](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/tests/curve/CurveTests.cpp); [tests/surface/SurfaceTests.cpp](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/tests/surface/SurfaceTests.cpp).
- Remaining gap: A bounded analytic subset only: no implemented sphere, cone, torus or general conics; evaluator input contracts still incomplete.

#### G02 — maturity 0; comparison —

No Bézier/B-spline/NURBS production types found; roadmap only.

- Implementation/design evidence: [src/curve/Curve.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/curve/Curve.cppm); [src/surface/Surface.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/surface/Surface.cppm); [TODO.md](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/TODO.md).
- Test evidence: No implementation acceptance tests found.
- Remaining gap: Implement scalar-generic rational evaluation, knot operations, continuity and independent reference tests.

#### G03 — maturity 1; comparison —

Curve2d wrappers and periodic circle/cylinder domains are useful prerequisites.

- Implementation/design evidence: [src/curve/Erased.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/curve/Erased.cppm); [src/curve/Circle.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/curve/Circle.cppm); [src/surface/Cylinder.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/surface/Cylinder.cppm); [docs/adr/0009-parameter-space-representation.md](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/docs/adr/0009-parameter-space-representation.md).
- Test evidence: [tests/curve/CurveTests.cpp](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/tests/curve/CurveTests.cpp); [tests/surface/SurfaceTests.cpp](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/tests/surface/SurfaceTests.cpp).
- Remaining gap: No trimmed face/edge representation, coedge parameter curves or 3D/2D consistency certificates.

#### G04 — maturity 0; comparison —

No fitting, interpolation or geometry conversion implementation found.

- Implementation/design evidence: [src/curve/Curve.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/curve/Curve.cppm); [src/surface/Surface.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/surface/Surface.cppm); [README.md](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/README.md).
- Test evidence: No implementation acceptance tests found.
- Remaining gap: Build fitting/continuity contracts and residual-based acceptance cases.

#### T01 — maturity 0; comparison —

Topo module explicitly empty; handles and arenas are prerequisites, not a B-rep.

- Implementation/design evidence: [src/topo/Topo.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/topo/Topo.cppm); [docs/adr/0004-topology-storage-and-handles.md](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/docs/adr/0004-topology-storage-and-handles.md).
- Test evidence: No implementation acceptance tests found.
- Remaining gap: Implement incidence stores, coedges/radial rings, ShapeBuilder/Freeze and representation invariants.

#### T02 — maturity 0; comparison —

No shape validity checker or point-in-solid classifier; policy enums do not implement checks.

- Implementation/design evidence: [src/topo/Topo.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/topo/Topo.cppm); [docs/adr/0018-validity-checking.md](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/docs/adr/0018-validity-checking.md).
- Test evidence: No implementation acceptance tests found.
- Remaining gap: Implement structural validation with seeded defects before geometric/solid classification.

#### T03 — maturity 0; comparison —

No topological construction, copying or adjacency traversal.

- Implementation/design evidence: [src/topo/Topo.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/topo/Topo.cppm); [docs/adr/0004-topology-storage-and-handles.md](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/docs/adr/0004-topology-storage-and-handles.md).
- Test evidence: No implementation acceptance tests found.
- Remaining gap: Prove shared incidence, immutable snapshots, orientation and traversal termination on actual shapes.

#### T04 — maturity 0; comparison —

History relations and attribute propagation are designed in a proposed ADR only.

- Implementation/design evidence: [docs/adr/0011-history-and-naming.md](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/docs/adr/0011-history-and-naming.md); [src/topo/Topo.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/topo/Topo.cppm).
- Test evidence: No implementation acceptance tests found.
- Remaining gap: Implement OperationResult/history composition and split/merge/deletion tests.

#### Q01 — maturity 2; comparison —

Projection is implemented for the four analytic families and dispatched through erased wrappers.

- Implementation/design evidence: [src/curve/Line.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/curve/Line.cppm); [src/curve/Circle.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/curve/Circle.cppm); [src/surface/Plane.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/surface/Plane.cppm); [src/surface/Cylinder.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/surface/Cylinder.cppm).
- Test evidence: [tests/curve/CurveTests.cpp](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/tests/curve/CurveTests.cpp); [tests/surface/SurfaceTests.cpp](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/tests/surface/SurfaceTests.cpp).
- Remaining gap: No general trimmed/NURBS projection or solution enumeration; circle/cylinder axial ambiguity chooses parameter zero.

#### Q02 — maturity 0; comparison —

Intersect module is empty. Analytic SSI corpus and rational reference classification are test tooling only.

- Implementation/design evidence: [src/intersect/Intersect.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/intersect/Intersect.cppm); [tools/check_corpus.py](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/tools/check_corpus.py); [corpus/analytic-ssi/README.md](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/corpus/analytic-ssi/README.md).
- Test evidence: [tools/check_corpus.py](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/tools/check_corpus.py).
- Remaining gap: Implement analytic intersection results and run kernel outputs against fixtures; then extend to general curves/surfaces.

#### Q03 — maturity 1; comparison —

Point distances, interval boxes and built-in bounding queries exist.

- Implementation/design evidence: [src/core/Vector.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/core/Vector.cppm); [src/core/Box.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/core/Box.cppm); [src/curve/Circle.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/curve/Circle.cppm).
- Test evidence: [tests/core/VectorTests.cpp](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/tests/core/VectorTests.cpp); [tests/core/BoxTests.cpp](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/tests/core/BoxTests.cpp).
- Remaining gap: No general geometric extrema, shape distance or spatial acceleration index.

#### Q04 — maturity 0; comparison —

No geometric integration or mass-property engine found; vector length is not this capability.

- Implementation/design evidence: [src/build/Build.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/build/Build.cppm); [docs/benchmarks.md](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/docs/benchmarks.md).
- Test evidence: No implementation acceptance tests found.
- Remaining gap: Implement length/area/volume/centroid/inertia with analytic references and convergence tests.

#### M01 — maturity 0; comparison —

Build module empty; Cylinder is an unbounded surface evaluator, not a solid constructor.

- Implementation/design evidence: [src/build/Build.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/build/Build.cppm); [src/surface/Cylinder.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/surface/Cylinder.cppm).
- Test evidence: No implementation acceptance tests found.
- Remaining gap: Construct and validate a closed box and cylinder, then extrusions and revolutions.

#### M02 — maturity 0; comparison —

No sweeps, pipes or loft algorithms.

- Implementation/design evidence: [src/build/Build.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/build/Build.cppm); [README.md](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/README.md).
- Test evidence: No implementation acceptance tests found.
- Remaining gap: Requires topology plus curve/surface construction, frame laws and self-intersection policy.

#### M03 — maturity 0; comparison —

Boolean module explicitly empty.

- Implementation/design evidence: [src/boolean/Boolean.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/boolean/Boolean.cppm); [docs/adr/0013-boolean-classification.md](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/docs/adr/0013-boolean-classification.md).
- Test evidence: No implementation acceptance tests found.
- Remaining gap: Depends on validated B-rep, intersections and classification; no Boolean performance claim is currently supported.

#### M04 — maturity 0; comparison —

No sewing or shell/solid assembly implementation.

- Implementation/design evidence: [src/build/Build.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/build/Build.cppm); [src/topo/Topo.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/topo/Topo.cppm).
- Test evidence: No implementation acceptance tests found.
- Remaining gap: Implement bounded gap handling and closure/orientation checks.

#### M05 — maturity 0; comparison —

Blend module empty; offset/thickening/hollowing not implemented.

- Implementation/design evidence: [src/blend/Blend.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/blend/Blend.cppm); [README.md](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/README.md).
- Test evidence: No implementation acceptance tests found.
- Remaining gap: Define singularity and feature-collapse outcomes before advancing to general shapes.

#### M06 — maturity 0; comparison —

No fillet or chamfer implementation; variable-radius blends are outside current repository v1.0 scope.

- Implementation/design evidence: [src/blend/Blend.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/blend/Blend.cppm); [README.md](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/README.md).
- Test evidence: No implementation acceptance tests found.
- Remaining gap: Constant-radius/chamfer first; retain variable-radius support in the broader parity roadmap.

#### M07 — maturity 0; comparison —

No draft or local feature construction implementation.

- Implementation/design evidence: [src/build/Build.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/build/Build.cppm); [src/blend/Blend.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/blend/Blend.cppm).
- Test evidence: No implementation acceptance tests found.
- Remaining gap: Add explicit coverage and geometric acceptance cases to long-term roadmap.

#### M08 — maturity 0; comparison —

No surface filling/skinning construction engine.

- Implementation/design evidence: [src/surface/Surface.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/surface/Surface.cppm); [src/build/Build.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/build/Build.cppm).
- Test evidence: No implementation acceptance tests found.
- Remaining gap: Requires spline foundation and continuity-constrained construction.

#### M09 — maturity 0; comparison —

No simplification or feature-removal implementation.

- Implementation/design evidence: [src/build/Build.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/build/Build.cppm); [src/blend/Blend.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/blend/Blend.cppm).
- Test evidence: No implementation acceptance tests found.
- Remaining gap: Add a dedicated workflow and retained-region deviation budget.

#### H01 — maturity 0; comparison —

No shape diagnosis implementation; validity-level configuration is preparatory.

- Implementation/design evidence: [src/topo/Topo.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/topo/Topo.cppm); [docs/adr/0018-validity-checking.md](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/docs/adr/0018-validity-checking.md).
- Test evidence: No implementation acceptance tests found.
- Remaining gap: Start with deterministic structural defect reports and offending-entity IDs.

#### H02 — maturity 0; comparison —

No geometry/topology healing implementation.

- Implementation/design evidence: [src/io/Io.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/io/Io.cppm); [docs/benchmarks.md](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/docs/benchmarks.md).
- Test evidence: No implementation acceptance tests found.
- Remaining gap: Keep repair separate from diagnosis; enforce displacement and tolerance budgets.

#### H03 — maturity 0; comparison —

No splitting/merging/regularization operations.

- Implementation/design evidence: [src/topo/Topo.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/topo/Topo.cppm); [src/build/Build.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/build/Build.cppm).
- Test evidence: No implementation acceptance tests found.
- Remaining gap: Require history and geometry-deviation evidence for every representation change.

#### H04 — maturity 0; comparison —

No operational repair pipeline.

- Implementation/design evidence: [src/io/Io.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/io/Io.cppm); [src/core/Journal.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/core/Journal.cppm).
- Test evidence: No implementation acceptance tests found.
- Remaining gap: Journaling alone does not supply repair replay or rollback.

#### X01 — maturity 0; comparison —

Io empty; native format and reflection serialization remain proposed.

- Implementation/design evidence: [src/io/Io.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/io/Io.cppm); [docs/adr/0012-native-file-format.md](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/docs/adr/0012-native-file-format.md).
- Test evidence: No implementation acceptance tests found.
- Remaining gap: Implement versioned shape persistence and corrupt-input/round-trip tests.

#### X02 — maturity 0; comparison —

No STEP geometry reader or writer.

- Implementation/design evidence: [src/io/Io.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/io/Io.cppm); [README.md](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/README.md).
- Test evidence: No implementation acceptance tests found.
- Remaining gap: Inventory AP242 subset and implement independent round-trip checks.

#### X03 — maturity 0; comparison —

No STEP assembly/metadata implementation.

- Implementation/design evidence: [src/io/Io.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/io/Io.cppm); [README.md](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/README.md).
- Test evidence: No implementation acceptance tests found.
- Remaining gap: Separate placement/instancing, styling and PMI coverage from plain geometry.

#### X04 — maturity 0; comparison —

No IGES implementation; not explicitly listed in current Io module scope.

- Implementation/design evidence: [src/io/Io.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/io/Io.cppm); [README.md](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/README.md).
- Test evidence: No implementation acceptance tests found.
- Remaining gap: Add IGES read/write to the full-parity backlog.

#### X05 — maturity 0; comparison —

STL/OBJ/glTF named as planned Io scope, with no code beyond an empty module.

- Implementation/design evidence: [src/io/Io.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/io/Io.cppm).
- Test evidence: No implementation acceptance tests found.
- Remaining gap: Implement one format with orientation/scale tests, then enumerate remaining baseline formats.

#### X06 — maturity 0; comparison —

OCCT comparison/bridge described in benchmark plans; no bridge implementation found.

- Implementation/design evidence: [docs/benchmarks.md](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/docs/benchmarks.md); [TODO.md](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/TODO.md).
- Test evidence: No implementation acceptance tests found.
- Remaining gap: Pin OCCT and provide a minimal independently checked conversion path.

#### V01 — maturity 0; comparison —

Mesh module explicitly empty.

- Implementation/design evidence: [src/mesh/Mesh.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/mesh/Mesh.cppm); [docs/adr/0014-mesh-generation-boundary.md](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/docs/adr/0014-mesh-generation-boundary.md).
- Test evidence: No implementation acceptance tests found.
- Remaining gap: Add boundary tessellation, shared-edge consistency and mesh/face association.

#### V02 — maturity 0; comparison —

Debug viewer listed as TODO; no tracked viewer directory at this commit.

- Implementation/design evidence: [TODO.md](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/TODO.md); [README.md](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/README.md).
- Test evidence: No implementation acceptance tests found.
- Remaining gap: Start with diagnostic points/curves/bounds before full CAD scene interaction.

#### V03 — maturity 0; comparison —

No drawing/annotation engine; drafting explicitly excluded from repository v1.0.

- Implementation/design evidence: [README.md](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/README.md).
- Test evidence: No implementation acceptance tests found.
- Remaining gap: Retain as post-v1.0/full-platform work, not as completed kernel scope.

#### A01 — maturity 0; comparison —

No document/attribute/transaction system; arena snapshots are only a building block.

- Implementation/design evidence: [src/core/Arena.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/core/Arena.cppm); [docs/adr/0005-immutability-and-sharing.md](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/docs/adr/0005-immutability-and-sharing.md).
- Test evidence: No implementation acceptance tests found.
- Remaining gap: Build document identity, attribute persistence and undo/redo semantics.

#### A02 — maturity 0; comparison —

History ADR reserves future persistent IDs; full naming/regeneration is not implemented.

- Implementation/design evidence: [docs/adr/0011-history-and-naming.md](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/docs/adr/0011-history-and-naming.md).
- Test evidence: No implementation acceptance tests found.
- Remaining gap: Add application-level reference resolution with explicit ambiguous/deleted outcomes.

#### A03 — maturity 0; comparison —

No assembly or document exchange integration.

- Implementation/design evidence: [src/topo/Topo.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/topo/Topo.cppm); [src/io/Io.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/io/Io.cppm).
- Test evidence: No implementation acceptance tests found.
- Remaining gap: Needs instances/placements, document services and metadata exchange.

#### A04 — maturity 1; comparison —

CTest and Python source/corpus tools provide a diagnostic harness subset.

- Implementation/design evidence: [tests/CMakeLists.txt](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/tests/CMakeLists.txt); [tools/check_corpus.py](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/tools/check_corpus.py); [python/StereonModule.cpp](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/python/StereonModule.cpp).
- Test evidence: [tests/CMakeLists.txt](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/tests/CMakeLists.txt).
- Remaining gap: No headless geometry operation runner; Python module exposes only __version__.

#### P01 — maturity 1; comparison —

CMake component targets and a nanobind version module exist; ABI/package policy documented.

- Implementation/design evidence: [CMakeLists.txt](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/CMakeLists.txt); [src/CMakeLists.txt](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/CMakeLists.txt); [python/StereonModule.cpp](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/python/StereonModule.cpp); [docs/adr/0019-library-linking-packaging-and-abi.md](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/docs/adr/0019-library-linking-packaging-and-abi.md); [docs/adr/0020-cpp-c-abi-and-python-bindings.md](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/docs/adr/0020-cpp-c-abi-and-python-bindings.md).
- Test evidence: [.github/workflows/python.yml](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/.github/workflows/python.yml).
- Remaining gap: No implemented stable C ABI or installed consumer SDK; reconcile matrix recommendation with proposed C++-first/Python-direct policy.

#### P02 — maturity 1; comparison —

Thread-safe journal tests, shared-arena read/copy-edit tests and execution wrappers exist.

- Implementation/design evidence: [src/core/Arena.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/core/Arena.cppm); [src/core/Journal.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/core/Journal.cppm); [src/core/include/Stereon/Core/Exec.hpp](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/core/include/Stereon/Core/Exec.hpp).
- Test evidence: [tests/core/ArenaTests.cpp](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/tests/core/ArenaTests.cpp); [tests/core/ContextTests.cpp](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/tests/core/ContextTests.cpp); [tests/core/ExecTests.cpp](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/tests/core/ExecTests.cpp).
- Remaining gap: No geometry-operation determinism across thread counts; pinned TSan build fails before tests.

#### P03 — maturity 1; comparison —

Robust, arena and SIMD microbenchmarks compile and smoke-run in CI.

- Implementation/design evidence: [bench/RobustBench.cpp](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/bench/RobustBench.cpp); [bench/ArenaBench.cpp](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/bench/ArenaBench.cpp); [bench/SimdBench.cpp](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/bench/SimdBench.cpp); [.github/workflows/bench.yml](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/.github/workflows/bench.yml).
- Test evidence: [.github/workflows/bench.yml](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/.github/workflows/bench.yml).
- Remaining gap: No implemented OCCT comparison, industrial shape workload or stored comparative trend report found.

#### P04 — maturity 2; comparison —

Catch2/RapidCheck, sanitizer jobs, lint and licensed corpus validation are operational subsets.

- Implementation/design evidence: [tests/CMakeLists.txt](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/tests/CMakeLists.txt); [.github/workflows/linux.yml](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/.github/workflows/linux.yml); [.github/workflows/lint.yml](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/.github/workflows/lint.yml); [corpus/LICENSES.md](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/corpus/LICENSES.md).
- Test evidence: [tests/robust/PredicateTests.cpp](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/tests/robust/PredicateTests.cpp); [tools/check_corpus.py](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/tools/check_corpus.py).
- Remaining gap: Fuzzing and OCCT harness still TODO; Release/TSan CI failures prevent a clean gate; corpus tests are not kernel SSI tests.

#### P05 — maturity 1; comparison —

Contribution/governance files, 24 proposed ADRs, documentation and CODEOWNERS exist.

- Implementation/design evidence: [CONTRIBUTING.md](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/CONTRIBUTING.md); [.github/CODEOWNERS](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/.github/CODEOWNERS); [docs/adr/README.md](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/docs/adr/README.md); [docs/index.rst](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/docs/index.rst); [README.md](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/README.md).
- Test evidence: [.github/workflows/deploy-docs.yml](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/.github/workflows/deploy-docs.yml).
- Remaining gap: Pre-alpha API explicitly unstable; maturity and support evidence limited; do not treat proposed ADRs as accepted.

#### D01 — maturity 1; comparison —

Scalar-generic primitives/evaluators and forward dual-number checks make AD readiness concrete.

- Implementation/design evidence: [src/core/Scalar.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/core/Scalar.cppm); [src/core/Vector.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/core/Vector.cppm); [src/core/Matrix.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/core/Matrix.cppm); [docs/adr/0021-shape-derivatives.md](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/docs/adr/0021-shape-derivatives.md).
- Test evidence: [tests/support/Dual.hpp](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/tests/support/Dual.hpp); [tests/core/VectorTests.cpp](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/tests/core/VectorTests.cpp); [tests/core/MatrixTests.cpp](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/tests/core/MatrixTests.cpp); [tests/curve/CurveTests.cpp](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/tests/curve/CurveTests.cpp); [tests/surface/SurfaceTests.cpp](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/tests/surface/SurfaceTests.cpp).
- Remaining gap: Erased geometry interfaces use double; no full sensitivity API, correspondence/status implementation, Boolean derivatives or topology-event tests.

#### D02 — maturity 1; comparison —

Log/journal sinks, MemoryJournal and error journal IDs implemented with tests.

- Implementation/design evidence: [src/core/Log.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/core/Log.cppm); [src/core/Journal.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/core/Journal.cppm); [src/core/Error.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/core/Error.cppm).
- Test evidence: [tests/core/ContextTests.cpp](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/tests/core/ContextTests.cpp).
- Remaining gap: Journal source explicitly defers replay to first operations; entity lineage and geometry-operation recovery absent.

#### D03 — maturity 0; comparison —

A tagged duct-to-volume-mesher example is specified in the meshing ADR, not implemented.

- Implementation/design evidence: [docs/adr/0014-mesh-generation-boundary.md](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/docs/adr/0014-mesh-generation-boundary.md); [docs/adr/0021-shape-derivatives.md](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/docs/adr/0021-shape-derivatives.md).
- Test evidence: No implementation acceptance tests found.
- Remaining gap: Build tagged CAD → surface mesh → external volume mesh workflow before full solver sensitivity coupling.

#### E01 — maturity 0; comparison —

No proprietary translator modules found.

- Implementation/design evidence: [src/io/Io.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/io/Io.cppm).
- Test evidence: No implementation acceptance tests found.
- Remaining gap: Customer demand and licensing route must precede translator commitments.

#### E02 — maturity 0; comparison —

External volume mesher integration planned; internal volume meshing is deliberately outside kernel scope.

- Implementation/design evidence: [docs/adr/0014-mesh-generation-boundary.md](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/docs/adr/0014-mesh-generation-boundary.md); [src/mesh/Mesh.cppm](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/mesh/Mesh.cppm).
- Test evidence: No implementation acceptance tests found.
- Remaining gap: Implement and validate Gmsh/snappyHexMesh adapters after watertight tagged surface output.

#### E03 — maturity 0; comparison —

No full sketch/assembly solver or feature-authoring engine; parametric history outside v1.0.

- Implementation/design evidence: [README.md](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/README.md); [docs/adr/0011-history-and-naming.md](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/docs/adr/0011-history-and-naming.md).
- Test evidence: No implementation acceptance tests found.
- Remaining gap: Plan as application services and integration; operation history is not a parametric feature tree.

#### E04 — maturity 0; comparison —

No CAM, recognition or reconstruction product modules found.

- Implementation/design evidence: [src/CMakeLists.txt](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/src/CMakeLists.txt); [README.md](https://github.com/onurtuncer/Stereon/blob/cf1dbaa925fb9c60fc642ca6962f329af7d52047/README.md).
- Test evidence: No implementation acceptance tests found.
- Remaining gap: Keep separate demand-led subprojects with named users and acceptance datasets.

### Roadmap reconciliation

The broader objective agreed in this conversation exceeds the repository's current v1.0 scope. Keep both visible rather than treating v1.0 as comprehensive OCCT parity.

| Topic | Repository position | Matrix treatment |
| --- | --- | --- |
| Variable-radius blends | Explicitly outside v1.0 | M06 remains a later full-parity requirement |
| Drafting | Explicitly outside v1.0 | V03 retained for platform scope |
| Full persistent naming / parametric history | History maps first; parametric history outside v1.0 | T04 and A02/E03 are separate deliverables |
| C ABI and packaging | Proposed ADR-0019/0020; C++ first, initial direct Python bindings | P01 follows this sequencing; C ABI is not an early geometry prerequisite |
| Volume meshing | External mesher integration boundary | E02 accepts a supported adapter; do not require an internal volume mesher |
| IGES, full documents/assemblies, drawing services | Not established by current implementation roadmap | Explicit coverage additions required for comprehensive rivalry |
| Release numbering | Repository Phase 0, v0.1/v0.3/v0.6/v0.8/v1.0 | Matrix R0–R5 are capability stages, not replacements for repository versions |

### Recommended next implementation sequence

1. **Restore the current verification gate.** Resolve the Surface/Dual module-loading failure, then obtain passing Release and TSan builds/tests at one pinned revision. Retain explicit known contract limitations; do not hide them with green badges.
2. **Implement the missing model representation.** T01–T03 and G03: stores, coedges, shared parameterization, ShapeBuilder/Freeze, structural validation and manual box/cylinder fixtures. Handles and arenas are ready foundations to exercise here.
3. **Connect the existing SSI corpus to real algorithms.** Q02: implement a bounded analytic surface-pair subset, return explicit unsupported outcomes elsewhere, and compare actual outputs against the 28 seed fixtures. Tolerant expectations are provisional and require their own verification.
4. **Deliver the first solid-modeling vertical slice.** M01/M03/T04/Q04: construct box and cylinder solids, implement a bounded Boolean case, preserve operation history, validate the result and check analytic volume.
5. **Develop spline geometry in parallel with this dependency sequence as staffing permits.** G02/G04/Q01: NURBS evaluation, knot operations, fitting and projection. This is indispensable for broad parity, not an optional optimization.
6. **Add a usable inspection/exchange path.** V01/V02/X01/X02: diagnostic viewer, tessellation and native persistence, followed by a bounded STEP round-trip. Keep error/tolerance evidence visible at every step.

These are prioritized deliverables, not calendar estimates. Full OCCT comparison infrastructure is still needed before any row can receive P or A comparative status. The existing docs/benchmarks.md provides a useful starting specification, but its targets are not measurements.

### Overall assessment

Stereon is a functioning and meaningfully tested foundation for a future kernel, with a small analytic geometry layer. Its present production source does not yet implement a B-rep solid modeler. The most consequential next advance is to connect storage, representation, validation and analytic intersections into a valid solid-producing workflow. Keep the full-platform destination, but judge progress by those executable workflows.
