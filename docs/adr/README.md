# Architecture Decision Records

This directory records Stereon's significant design decisions. Each ADR states the context, the options considered, the decision, its consequences, and how the decision is verified.

## Rules

- One file per decision, numbered sequentially: `NNNN-short-title.md`. Start from [`0000-template.md`](0000-template.md).
- Status moves from **Proposed** to **Accepted** after review and, where the ADR names one, after its verification prototype passes.
- Accepted ADRs are never edited except to fix typos. A change of mind is a new ADR that supersedes the old one; the old one's status becomes "Superseded by ADR-MMMM".
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
| [0015](0015-ai-contribution-policy.md) | AI contribution policy | 0 | Proposed |
| [0016](0016-cpp26-adoption-and-contracts.md) | C++26 adoption and contracts policy | 0 | Proposed |
