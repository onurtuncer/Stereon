# ADR-0025: Separate Topology from Geometry

- **Status:** Accepted
- **Date:** 2026-10-13
- **Accepted:** 2026-10-17 by Onur Tuncer, as direction: verification item 2 holds (the curve and surface libraries evaluate with no B-rep); items 1 and 3–8 are acceptance gates for the topology work as its entity types land
- **Owner:** Onur Tuncer
- **Phase:** 0
- **Scope:** The conceptual B-rep model and the boundary between topology, geometry, and their associations

## Context

A boundary representation must describe both how a boundary is connected and where that boundary lies in space. These are related but distinct responsibilities.

Topology describes entities, incidence, boundary ordering, and orientation. Geometry describes points, curves, surfaces, and their evaluation. A topological edge, for example, is a bounded boundary entity; a supporting curve is a mathematical locus that may extend beyond that edge. A face is a bounded, oriented region associated with a supporting surface, not the surface itself.

Conflating these concepts couples connectivity operations to geometric representations and encourages incorrect identity rules. Coincident geometry does not necessarily represent one topological entity, and shared supporting geometry does not imply shared topology.

Ian Stroud's *Boundary Representation Modelling Techniques* is a principal conceptual reference for this design. This ADR records Stereon's own architectural decision, not a transcription of the book or a commitment to its implementation layout.

## Options

### Embed geometry and connectivity in one inseparable entity model

This can simplify an initial implementation, but makes geometric representation changes affect topology APIs and encourages implicit coupling between geometric equality and topological identity. Rejected as the architectural contract. Physical co-location for performance remains permissible if semantic separation is preserved.

### Infer connectivity from coordinates or geometric proximity

This reduces explicit connectivity storage but cannot reliably distinguish coincident, intentionally separate boundaries from shared boundaries. It also makes traversal depend on geometric computations and tolerance choices. Rejected for the authoritative B-rep model.

### Keep topology and geometry unrelated until export

This maximizes independence but leaves core modelling operations without authoritative support, trimming, and parameter mappings. Rejected: explicit associations are part of the B-rep model itself.

## Decision

**Stereon shall represent topology and geometry as separate conceptual models, connected through explicit, validated associations.**

This separation is mandatory at the semantic and API boundaries. It does not require a particular allocation strategy, class hierarchy, module layout, or physical storage arrangement.

### 1. Topology owns connectivity and orientation

The topology model shall describe topological identity, incidence, ordered boundaries, and orientation relative to a containing entity or use, where applicable.

Topological traversal and structural inspection shall not require evaluation of a curve or surface. Geometry types such as lines, circles, and B-spline surfaces shall not determine the identity or connectivity of topological entities.

ADR-0004 already selects typed index handles, coedges with radial rings, separate geometry pools, and non-manifold intermediate support. This ADR states the semantic separation underlying that representation; it does not replace those choices.

### 2. Geometry owns mathematical representation and evaluation

The geometry model shall describe geometric supports and their intrinsic parameter domains, including point values, curves, and surfaces. Geometric evaluation and any supported derivative evaluation shall be available without requiring a complete B-rep model.

Geometric objects shall not own B-rep connectivity or require authoritative backlinks to incident topological entities. Reverse lookup indexes may be maintained separately as derived data.

### 3. Associations describe how topology uses geometry

The B-rep association layer shall explicitly record the bindings and restrictions needed to interpret topology geometrically. Conceptually, these include:

| Association | Responsibility |
| --- | --- |
| Vertex to point | Associate a topological vertex with its geometric position. |
| Edge to geometric representation | Follow ADR-0009's edge kind and shared parameter interval. For `Exact` and `Spatial` edges the 3D curve is authoritative; surface-defined kinds use their pcurves and supporting surfaces, with a certified 3D cache where specified. |
| Face to supporting surface | Identify the surface supporting the face; the face's boundary topology and associated trimming data delimit its region. |
| Edge occurrence in a face to parameter-space curve | Describe that boundary occurrence in the supporting surface's parameter space, when required by the representation. |
| Topological direction to geometric parameter direction | Record their relationship explicitly rather than infer it from coordinate values. |

The cardinality of these associations shall accommodate periodic seams and representation-specific degeneracies. In particular, a single edge can have distinct parameter-space boundary occurrences on the same face. A degenerate edge shall not be forced to carry a fictitious nondegenerate 3D curve.

These requirements define association semantics; they do not prescribe a dedicated object type or separate allocation for every association. Under ADR-0004, pcurves belong to coedges and surfaces are reached through faces; an edge does not duplicate surface references. ADR-0009 remains authoritative for edge kinds, primary-side evaluation, shared parameterization, degeneracies, and certificates.

### 4. Identity is independent of geometric coincidence

Distinct topological entities may occupy the same geometric location. Distinct edges or faces may also refer to the same supporting curve or surface while representing different bounded regions.

Geometric equality, proximity, or shared geometry references shall not implicitly merge topological entities. Sewing, merging, or other changes to connectivity must be explicit modelling operations subject to their own contracts.

Conversely, sharing a geometric support shall not be interpreted as a topological adjacency relation.

### 5. Geometric edits do not silently redefine topology

A geometric support or association may be replaced while retaining topological identity when connectivity remains unchanged and the resulting B-rep satisfies its applicable validity requirements.

If an edit requires splitting, merging, deleting, or reconnecting entities, that change shall be expressed as a topology operation. Retaining connectivity alone does not establish geometric validity: an edit may still introduce self-intersection, collapse, or inconsistent trimming.

Updating shared geometry must have explicit effects on every dependent association. The implementation shall either preserve the old support for unaffected consumers or deliberately update and revalidate all affected consumers. Silent mutation that leaves dependent B-rep data stale is not permitted.

ADR-0005 supplies the implementation policy: immutable frozen shapes and geometry, with changes staged through a builder and copy-on-write storage. Existing input snapshots remain unchanged. Topological identity across operations is communicated by ADR-0011 history, not by treating shape-local handles as globally persistent identifiers.

### 6. Consistency is checked across the boundary

The separation shall preserve a clear distinction between:

- **Structural checks:** reference integrity, incidence, boundary ordering, and orientation constraints applicable to the supported topology.
- **Association and geometric checks:** valid parameter restrictions, endpoint agreement, and agreement between 3D boundary geometry and its surface-space representation where present.

Passing structural checks shall not be reported as proof of a valid geometric solid. Likewise, geometric coincidence shall not be treated as proof of topological connectivity.

Numerical acceptance criteria follow ADR-0001, with robustness and certificates governed by ADR-0002 and ADR-0009. Validation uses ADR-0018's Structural, Geometric, and Solid levels. This ADR introduces no global epsilon and does not weaken any existing certificate requirement.

## Rationale

This decision allows topology algorithms to operate on connectivity without depending on a specific curve or surface representation. It also permits geometric evaluators to be reused independently of the B-rep graph.

Explicit associations make boundedness, trimming, parameter direction, and seam behavior visible responsibilities. They avoid embedding assumptions such as “one curve equals one edge” or “one surface equals one face” in the core architecture.

The separation supports future geometric representation changes and derivative evaluation without requiring a redesign of topology. It does not, by itself, make topology-changing operations differentiable or guarantee persistent naming through arbitrary edits.

## Consequences

### Benefits

- Connectivity can be inspected and tested independently of geometric evaluation.
- Geometric supports can be reused without forcing topological sharing.
- Supporting geometry can be replaced while preserving valid, unchanged topology.
- Association inconsistencies become identifiable errors rather than hidden assumptions.
- Geometry representation choices remain independent of the connectivity representation.

### Costs and obligations

- Associations introduce bookkeeping and potentially additional indirection.
- Geometry changes require dependency tracking, cache invalidation, and appropriate revalidation.
- APIs must distinguish topological entities from geometric supports and their bounded uses.
- Serialization must preserve both models and their associations without inferring identity from coordinates.

## Verification

The initial implementation shall demonstrate the following behaviors as the relevant entity types become available:

1. A topology-only traversal performs no geometric evaluation.
2. A curve or surface evaluator can be used without constructing a B-rep graph.
3. Two topologically distinct edges can share a supporting curve with independent parameter intervals.
4. Coincident vertices or edges remain distinct until an explicit operation changes connectivity.
5. Replacing geometry while preserving valid connectivity does not unnecessarily replace topological identities.
6. An inconsistent geometric association can be diagnosed even when the topology is structurally valid.
7. Seam associations distinguish multiple parameter-space occurrences of an edge on a periodic face.
8. A geometric edit cannot silently leave other consumers of shared geometry with stale associations or caches.

These are implementation acceptance requirements, not claims that Stereon currently satisfies them. Seam and degeneracy support may be staged, but the association design must not rule them out.

## Relationship to existing decisions

- [ADR-0004](0004-topology-storage-and-handles.md) defines storage, handles, radial rings, and geometry pools.
- [ADR-0005](0005-immutability-and-sharing.md) defines immutable snapshots and sharing.
- [ADR-0009](0009-parameter-space-representation.md) defines authoritative edge geometry and associations.
- [ADR-0018](0018-validity-checking.md) defines validation levels.
- [ADR-0021](0021-shape-derivatives.md) defines derivative semantics; separation alone does not establish differentiability.

This proposal makes the semantic boundary explicit and supplements those decisions. It supersedes none of them. Further storage optimization and association APIs must preserve both this boundary and their existing guarantees.

## References

1. Ian Stroud, *Boundary Representation Modelling Techniques*, Springer, 2006. DOI: [10.1007/978-1-84628-616-2](https://doi.org/10.1007/978-1-84628-616-2). Relevant reading: “Datastructures and tools,” “Special representations,” and “Free-form geometry.” The publisher's contents and scope have been verified; detailed chapter-level verification remains pending. No specific passage or implementation prescription is attributed to the book here.
2. Open CASCADE Technology, [ShapeAnalysis_Edge documentation](https://occt3d.com/dev/doc/refman/html/class_shape_analysis___edge.html), illustrating explicit checks for agreement between edge geometry and parameter-space representations. This is supporting implementation evidence, not an API compatibility requirement.

