# Architecture Decision Records

This directory records Stereon's significant design decisions. Each ADR states the context, the options considered, the decision, its consequences, and how the decision is verified.

Open questions that have to be researched before they can be decided live in [`docs/research/`](../research/) as research specifications. A go decision there leads to a new or revised ADR here.

## Rules

- One file per decision, numbered sequentially: `NNNN-short-title.md`. Start from [`0000-template.md`](0000-template.md).
- Status moves from **Proposed** to **Accepted** after review and, where the ADR names one, after its verification prototype passes.
- Proposed ADRs may be revised in place; a substantive revision adds a **Revised** line under the date naming what changed. Accepted ADRs are never edited except to fix typos. A change of mind is a new ADR that supersedes the old one; the old one's status becomes "Superseded by ADR-MMMM".
- A pull request that contradicts an accepted ADR is rejected unless it comes with a new ADR.
- All ADRs are reviewed at every release gate.

## Index

| ADR | Title | Phase | Status |
| --- | --- | --- | --- |
| [0001](0001-tolerance-model.md) | Tolerance model | 0 | Accepted |
| [0002](0002-numeric-robustness.md) | Numeric robustness strategy | 0 | Accepted |
| [0003](0003-geometry-polymorphism.md) | Geometry polymorphism via external polymorphism | 0 | Accepted |
| [0004](0004-topology-storage-and-handles.md) | Topology storage and handles | 0 | Accepted |
| [0005](0005-immutability-and-sharing.md) | Immutability and structural sharing | 0 | Accepted |
| [0006](0006-error-handling.md) | Error handling | 0 | Accepted |
| [0007](0007-threading-and-determinism.md) | Threading and determinism | 0 | Proposed |
| [0008](0008-licence.md) | Licence | 0 | Accepted |
| [0009](0009-parameter-space-representation.md) | Parameter-space curves (pcurves) | 1 | Proposed |
| [0010](0010-ssi-strategy.md) | Surface–surface intersection strategy | 2 | Accepted in part |
| [0011](0011-history-and-naming.md) | Operation history and topological naming | 3 | Proposed |
| [0012](0012-native-file-format.md) | Native file format and versioning | 3 | Proposed |
| [0013](0013-boolean-classification.md) | Boolean classification | 4 | Accepted |
| [0014](0014-mesh-generation-boundary.md) | Mesh generation boundary | 5 | Proposed |
| [0015](0015-ai-contribution-policy.md) | AI contribution policy | 0 | Accepted |
| [0016](0016-cpp26-adoption-and-contracts.md) | C++26 adoption and contracts policy | 0 | Accepted |
| [0017](0017-naming-conventions.md) | Naming conventions (Hazel style) | 0 | Accepted |
| [0018](0018-validity-checking.md) | Validity checking | 0 | Proposed |
| [0019](0019-library-linking-packaging-and-abi.md) | Library linking, packaging and ABI policy | 0 | Proposed |
| [0020](0020-cpp-c-abi-and-python-bindings.md) | C++ API, C ABI and Python bindings | 0 | Proposed |
| [0021](0021-shape-derivatives.md) | Shape derivatives as an operation contract | 0 | Proposed |
| [0022](0022-function-complexity.md) | Function complexity limits | 0 | Accepted |
| [0023](0023-logging-and-journal-sinks.md) | Logging and journal sinks | 0 | Proposed |
| [0024](0024-named-constants.md) | Named constants | 0 | Accepted |
| [0025](0025-separate-topology-from-geometry.md) | Separate topology from geometry | 0 | Accepted |
| [0026](0026-modelling-operations-by-contracts.md) | Specify modelling operations by contracts | 0 | Proposed |
| [0027](0027-controlled-topology-mutation-and-euler-operators.md) | Controlled topology mutation and Euler operators | 0 | Accepted |
| [0028](0028-purpose-specific-numerical-tolerances.md) | Separate numerical tolerances by purpose | 0 | Proposed |
| [0029](0029-intersection-and-classification-outcomes.md) | Define intersection and classification outcomes explicitly | 2 | Accepted |
| [0030](0030-entity-lineage-in-operation-results.md) | Make entity lineage part of operation results | 3 | Accepted |
| [0031](0031-attribute-propagation-and-conflicts.md) | Specify attribute propagation and conflict policies | 3 | Proposed |
| [0032](0032-mutation-based-kernel-verification.md) | Mutation-based verification of kernel contracts | 0 | Accepted |
| [0033](0033-modular-sdk-product-scope.md) | Modular CAD SDK with an independent kernel | 0 | Proposed |
| [0034](0034-kernel-presentation-boundary.md) | Separate the modelling kernel from presentation and interaction | 0 | Proposed |
| [0035](0035-debug-viewer.md) | Debug viewer over dump files with GLFW, OpenGL and Dear ImGui | 0 | Proposed |
| [0036](0036-typed-topology-ranges.md) | Typed read-only topology ranges | 0 | Proposed |
| [0037](0037-topology-traversal-semantics.md) | Entity, use and orientation semantics for topology traversal | 0 | Proposed |
| [0038](0038-topology-traversal-lifetime-and-validation.md) | Topology traversal lifetime, invalidation and validation | 0 | Proposed |
| [0039](0039-executable-topology-contracts.md) | Make topology-operator contracts executable and auditable | 0 | Proposed |
| [0040](0040-atomic-topology-edit-publication.md) | Stage topology edits and publish complete results atomically | 0 | Proposed |
| [0041](0041-layered-validation-evidence.md) | Use layered validation with revision-bound evidence | 0 | Proposed |
| [0042](0042-cpp26-contract-validation-responsibilities.md) | Assign C++26 contracts and kernel validation distinct responsibilities | 0 | Proposed |
| [0043](0043-optional-performance-instrumentation.md) | Optional performance instrumentation and reproducible benchmarks | 0 | Proposed |
| [0044](0044-brep-to-tri-dexel-conversion.md) | B-rep to tri-dexel conversion with bounded endpoints and explicit unresolved rays | 3–5 | Proposed |
| [0045](0045-cppad-ad-backend.md) | CppAD as the automatic-differentiation backend | 0 | Proposed |
| [0046](0046-ownership-borrowing-and-pointer-lifetimes.md) | Explicit ownership, borrowing and pointer lifetimes | 0 | Proposed |
