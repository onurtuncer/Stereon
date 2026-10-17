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
| [0001](0001-tolerance-model.md) | Tolerance model | 0 | Proposed |
| [0002](0002-numeric-robustness.md) | Numeric robustness strategy | 0 | Proposed |
| [0003](0003-geometry-polymorphism.md) | Geometry polymorphism via external polymorphism | 0 | Proposed |
| [0004](0004-topology-storage-and-handles.md) | Topology storage and handles | 0 | Proposed |
| [0005](0005-immutability-and-sharing.md) | Immutability and structural sharing | 0 | Proposed |
| [0006](0006-error-handling.md) | Error handling | 0 | Proposed |
| [0007](0007-threading-and-determinism.md) | Threading and determinism | 0 | Proposed |
| [0008](0008-licence.md) | Licence | 0 | Proposed |
| [0009](0009-parameter-space-representation.md) | Parameter-space curves (pcurves) | 1 | Proposed |
| [0010](0010-ssi-strategy.md) | Surface–surface intersection strategy | 2 | Proposed |
| [0011](0011-history-and-naming.md) | Operation history and topological naming | 3 | Proposed |
| [0012](0012-native-file-format.md) | Native file format and versioning | 3 | Proposed |
| [0013](0013-boolean-classification.md) | Boolean classification | 4 | Proposed |
| [0014](0014-mesh-generation-boundary.md) | Mesh generation boundary | 5 | Proposed |
| [0015](0015-ai-contribution-policy.md) | AI contribution policy | 0 | Accepted |
| [0016](0016-cpp26-adoption-and-contracts.md) | C++26 adoption and contracts policy | 0 | Proposed |
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
| [0027](0027-controlled-topology-mutation-and-euler-operators.md) | Controlled topology mutation and Euler operators | 0 | Proposed |
| [0028](0028-purpose-specific-numerical-tolerances.md) | Separate numerical tolerances by purpose | 0 | Proposed |
| [0029](0029-intersection-and-classification-outcomes.md) | Define intersection and classification outcomes explicitly | 2 | Proposed |
| [0030](0030-entity-lineage-in-operation-results.md) | Make entity lineage part of operation results | 3 | Proposed |
| [0031](0031-attribute-propagation-and-conflicts.md) | Specify attribute propagation and conflict policies | 3 | Proposed |
| [0032](0032-mutation-based-kernel-verification.md) | Mutation-based verification of kernel contracts | 0 | Accepted |
