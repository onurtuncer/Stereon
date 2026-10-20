# ADR-0029: Define Intersection and Classification Outcomes Explicitly

- **Status:** Accepted
- **Date:** 2026-10-14
- **Revised:** 2026-10-20 (the decisions on acceptance below: the prototype's scope, analytic tangency evidence, declared tolerant coincidence, `ErrorCode::Unresolved`)
- **Accepted:** 2026-10-20 by Onur Tuncer, before the Boolean prototype that builds on it (`docs/planning/boolean-prototype-plan.md`, decision D7)
- **Owner:** Onur Tuncer
- **Phase:** 2 (v0.2 intersection contract); classification integrated in v0.3–v0.6
- **Extends:** ADR-0006, ADR-0010, ADR-0013, ADR-0026

## Context

An empty intersection can be a proved answer or the result of failing to find a branch. A small residual can establish proximity without proving intersection completeness. Parallel surface normals establish a first-order tangency condition but do not by themselves prove whether nearby surface fragments change sides. Treating these as one Boolean result can produce plausible but incorrect topology.

Stroud §§6.1 and 7.3 motivate explicit handling of intersection special cases, coincidence and tangency. Stereon's outcome algebra, side-transition evidence and failure contract below extend that guidance and the existing SSI/Boolean proposals.

## Options

1. **Return geometry or an empty collection.** Compact but conflates disjointness, unsupported cases and failed computation.
2. **One enum mixing contact geometry and solver status.** Better diagnostics but cannot represent several simultaneous components or independent accuracy/completeness evidence.
3. **Separate geometric components, assurance and completion status, with explicit classification evidence.** More detailed results; downstream decisions can check exactly what is established.

## Decision

Adopt option 3. Preserve ADR-0006's `std::expected` boundary, ADR-0010's separate completeness/accuracy obligations and ADR-0013's classification strategy. Names below specify semantics; concrete structs and any new error codes are reviewed with implementation, not declared already available.

### 1. Successful intersections describe their dimension and domain

A complete successful result contains zero or more components; different contact kinds may coexist. Each component identifies both inputs, its spatial and parameter domains, orientation where meaningful, bounds/certificates and numerical policy (ADR-0028).

| Geometric outcome | Required meaning |
| --- | --- |
| Empty / disjoint | Separation established over the full requested domain; not merely no seeds found. |
| Isolated point contact | A zero-dimensional intersection, with location/parameter bounds and contact evidence. Tangent spheres can produce a point rather than a curve. |
| Crossing branch | A curve of transverse intersection on the declared regular domain. |
| Tangential branch | A curve with established tangential contact; side transition is established separately. |
| Overlap region | A region where the requested trimmed domains coincide under the stated exact/tolerant semantics; return its parameter correspondence, extent and boundary components. |
| Full-domain coincidence | The requested domains coincide completely; record their orientation/correspondence. It is a particular overlap extent, not identity of topology or C++ objects. |

Coincident supporting surfaces do not imply overlapping trimmed faces. An overlap can have no boundary in its requested domain. A generic result must not force all contacts into B-spline curves; exact analytic geometry remains exact where available under ADR-0009.

Distinguish nominal geometry from any explicitly authorized tolerant placement/healing. A policy that treats close geometry as coincident must state that interpretation, its evidence and any actual perturbation. This ADR does not approve RS-0001's research placement strategy or invent a nominal intersection when one does not exist.

### 2. Completion and assurance cannot be inferred from contact type

Return complete success only after the operation's declared completeness and accuracy obligations are satisfied. Record the supported domain and applicable evidence; a collection of certified curves does not prove no other curves exist.

Invalid data, unsupported cases, exhausted subdivision/iteration, and failed accuracy certification remain distinct diagnostics. Use existing `InvalidInput`, `NotSupported`, `NotImplemented`, `IntersectionFailed` and `CertificationFailed` codes where their meanings fit. Unresolved classification needs a distinguishable diagnostic reason and stage; it must not masquerade as `Outside`, `OnBoundary` or a geometric impossibility. New enum values, if needed, are appended when implemented without renumbering existing codes.

Private working results may contain unresolved parameter boxes and partial components. On a public complete-result API, these can accompany failure diagnostics, but are not a successful intersection or a valid published shape. Any future exploratory partial-result API must have a distinct contract and cannot feed ordinary Boolean classification as complete evidence.

### 3. Contact geometry and side transition are separate

For each relevant pair of fragments adjacent across an intersection branch, establish one of:

- **SameSide:** the adjacent fragments have the same inside/outside state relative to the other solid.
- **OppositeSide:** their inside/outside states differ.
- **Unresolved:** the relation is not established; refine, classify independently with adequate evidence, or fail before it is needed.

A regular transverse crossing implies an opposite-side relation when the local solid/domain assumptions hold. Tangential normal alignment alone implies neither relation. For example, in dimensionless local coordinates, surfaces `z = 0` and `z = x^2` touch along `x = 0`, while `z = 0` and `z = x^3` cross with tangential contact along the same line. Both have parallel normals at contact. Singularity, trim boundaries and higher-order contact need their own analysis or an unresolved result.

Split branches where contact type or side-transition evidence changes. An unresolved transition interval cannot be crossed by blindly propagating a classification. For overlap boundaries use the overlap/orientation evidence and independently justified neighboring states, not a generic toggle rule. ADR-0013's revised propagation and cross-check rules use this evidence.

### 4. Classification has explicit boundary and uncertainty semantics

Point/region queries distinguish `Inside`, `Outside` and `OnBoundary` when established. Fragment classification additionally distinguishes `OnSame` and `OnOpposite` for oriented overlap; these are not classifications of an arbitrary point without face orientation.

`OnBoundary` means the query satisfies the declared boundary relation under the model's numerical policy. An interval straddling an undecided test is **Unresolved**, not evidence of boundary membership. When a query cannot establish a supported state within its work limits, the complete public query returns failure with the relevant entities/domains and ambiguity evidence.

Seeding, deterministic ray retries, subdivision and winding-number fallbacks must preserve the requested assurance. A fine tessellation is not by itself a proof about the B-rep. A tessellation fallback must bound approximation/numerical error and establish that the classification is insensitive to it; otherwise refine or report unresolved. Dense sampling is a regression oracle, not completeness evidence.

Boolean assembly consumes only established states and valid contact relations. Lower-dimensional contact follows the public Boolean's declared regularized or general-set semantics (ADR-0026); do not silently add point/curve contacts to a regularized solid result or report failure for a legitimate empty result.

### 5. Diagnostics and publication are reproducible

Record the failed stage, input-qualified entities, remaining parameter domains, effective numerical policy, bounds and exhausted budgets where applicable. G1 deterministic ordering applies to components, unresolved work items and diagnostic collections; thread scheduling must not select a different result. Logs and journal IDs may be diagnostic-only as in ADR-0023.

A failure publishes no partial shape, history or metadata. Input shapes remain unchanged under ADR-0026/0027. Bindings preserve the same distinctions rather than mapping all failures to an empty list.

### 6. Decisions on acceptance (2026-10-20)

- **The full model, the prototype's cases.** The outcome type carries every kind of §1 (isolated points and full coincidence included), so that later surfaces need no redesign. The Boolean prototype produces what planes and cylinders give: proved empty, crossing lines and conics, tangential lines, coplanar and coaxial overlaps.
- **Tangency between quadrics is proved, not assumed.** For a plane tangent to a cylinder and for two parallel cylinders in contact (externally or internally), the side relation follows from the closed form: the surfaces touch without crossing, so the relation is SameSide, decided with exact predicates and recorded as proved evidence. Any tangential contact without such a proof is Unresolved. The x³ case of §3, a tangential crossing, does not arise between quadrics.
- **Tolerant coincidence is declared, and moves nothing.** Under nominal geometry (the prototype's decision D1), two supports that are exactly parallel (decided exactly on the input) and within the pair tolerance of each other are an Overlap, with evidence that says it is tolerant and gives the distance. No surface is moved: the result keeps the representative's surface (the lower input index, as in ADR-0001's merge policy), and the edge tolerances on the other side grow to cover the distance, logged and reported. Supports that are not exactly parallel are never coincident; they give a nominal crossing.
- **`ErrorCode::Unresolved`.** An unresolved classification, or any other established-state question that a public operation needs answered and cannot settle, fails with a new code, appended when first used (the next value, nothing renumbered), and a message naming the stage and the input-qualified entities.

## Consequences

Algorithms need result containers for points and regions as well as curves, and evidence for side transitions. Uncertain cases may fail more often initially; those failures become explicit corpus work rather than silent topology errors. The SSI construction strategy and exact-predicate boundary remain unchanged.

## Verification

1. Analytic fixtures distinguish disjoint planes, a transverse line, sphere/plane tangent point, tangent cylinders, coincident supporting surfaces with disjoint trims, partial overlap and full-domain coincidence.
2. The `x^2` versus `x^3` contact fixtures establish SameSide versus OppositeSide, or explicitly reject the case as unresolved/unsupported; neither may be classified from normal angle alone.
3. Mixed components and changes in contact/transition along a branch retain all components and their domains. A known small loop cannot disappear into a successful empty result when subdivision is exhausted.
4. Force work exhaustion and failed certificates independently: return failure, preserving diagnostic partial evidence and unchanged inputs; never publish nominal success with an incomplete set.
5. An ambiguous near-boundary point stays unresolved until membership is established. Tessellation fallback cases near the deviation bound must refine or fail, not guess a side.
6. For regularized Booleans, disjoint and boundary-only intersections give valid empty results when established; coincident same/opposite orientations follow the operation's keep/remove rules.
7. Independent analytic values and adversarial corpus cases supplement differential/sampling checks. G1 repeats compare components, classifications and semantic diagnostics across thread counts.

## References and scope

- Ian Stroud, *Boundary Representation Modelling Techniques*, Springer, 2006, §6.1 (printed pp. 126–141) and §7.3 (pp. 229–230); special-case discussion reviewed in the supplied book.
- [ADR-0006](0006-error-handling.md), [ADR-0009](0009-parameter-space-representation.md), [ADR-0010](0010-ssi-strategy.md), [ADR-0013](0013-boolean-classification.md), [ADR-0026](0026-modelling-operations-by-contracts.md), [ADR-0028](0028-purpose-specific-numerical-tolerances.md).

This proposal supplements the listed decisions and reconciles their proposed contact/output rules; it supersedes none and accepts no research result. The verification cases become implementation gates as the corresponding domains are supported.

### Status of the verification

| Case | Gate |
| --- | --- |
| 1: disjoint planes, a transverse line, tangent cylinders, coincident supports with disjoint trims, partial overlap, full coincidence | The Boolean prototype's analytic handlers (M2, M3) |
| 1: the sphere–plane tangent point | The first handler for spheres |
| 2: x² against x³ | Generic SSI (ADR-0010), where non-quadric contact arises |
| 3, 4: mixed components; forced exhaustion and failed certificates | M3 and M4 |
| 5: an ambiguous near-boundary point stays unresolved | M7, classification |
| 6: empty and boundary-only regularized results; same and opposite overlaps | M8, assembly |
| 7: independent analytic values, adversarial corpus, G1 repeats | M3 (corpus), M10 (G1) |
