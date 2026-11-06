# ADR-0049: Machining collision and clearance queries with bounded motion

- **Status:** Proposed
- **Date:** 2026-11-06
- **Owner:** Onur Tuncer
- **Phase:** M7, staged after the required solid and query capabilities
- **Related:** ADR-0002, ADR-0006, ADR-0007, ADR-0018, ADR-0021, ADR-0028, ADR-0029, ADR-0033, ADR-0041, ADR-0043, ADR-0044, ADR-0047

## Context

Machining verification must distinguish holder or shank interference, insufficient
clearance, intended cutter contact, and gouging of protected geometry. A useful
kernel must answer these questions against a declared geometric domain, with
numerical evidence and identifiable entities. Surface intersection alone is not
a solid collision test: one solid may contain another without their boundaries
meeting. Positive separation may still violate the required clearance.

Checking poses independently can miss a collision between samples, particularly
when a long tool assembly rotates. Constructing a Boolean common at every pose
also performs unnecessary topology construction and does not check the motion
between poses. The required capability is geometric proximity and continuous
verification, with explicit completion evidence.

ADR-0044 proposes stock representation and ray queries. It does not define solid
proximity, machining contact semantics or continuous collision contracts. This
proposal complements it and ADR-0033; it supersedes neither. It records a target
architecture, not an assertion that these services exist or are production-ready.

## Options

1. **Boolean common at sampled poses.** Reuses modelling operations, but builds
   unnecessary topology, needs separate distance queries, and misses events
   between poses unless additional motion bounds are established.
2. **Triangle or dexel collision only.** Supports fast approximate checks, but
   unbounded tessellation error or material between rays prevents a proof about
   the original CAD solid.
3. **Dedicated proximity queries with conservative acceleration, geometric
   refinement and bounded continuous motion.** Reuses Stereon's geometry and
   classification contracts while allowing early exits and representation-specific
   acceleration. Requires distance, containment and motion verification work.

## Decision

Adopt option 3, delivered incrementally for explicitly supported domains.

### 1. Responsibility and dependency boundaries

- Stereon owns reusable static solid overlap, containment, contact and bounded
  minimum-distance queries, plus verification over a supplied bounded rigid motion.
  These are headless capabilities above their geometry/topology dependencies,
  with no dependency from foundational modules back into machining services.
- The machining application owns cutter definitions, cutting versus non-cutting
  assembly regions, allowed-contact rules, protected material, machine kinematics,
  controller interpolation, toolpath strategy, postprocessing and stock-update
  policy. It supplies tool envelopes as supported geometry with explicit units,
  frames, instance identities and transforms. Swept-solid construction remains an
  application responsibility as in ADR-0044; motion verification need not construct
  a swept B-rep.
- The application converts machine-axis motion into rigid-body trajectories for
  every relevant component, including moving fixtures and tables. Endpoint poses
  alone do not specify a five-axis machine trajectory.
- Kernel geometry outcomes and application permission decisions remain separate.
  A cutting tool is not globally exempt from collision: its cutting region may
  remove designated stock but may not penetrate protected finished-part material.
- Exact library/module names are deferred to implementation under ADR-0033.
  C++26, module boundaries, value semantics and existing numerical policies apply.
  No new dependency, empty package or full machine simulator is introduced here.

### 2. Static query contract

Solid queries require closed, oriented, valid solids in the declared supported
subset (ADR-0018). Surface-only queries, if provided, have separate semantics and
cannot establish solid disjointness. Invalid input and unsupported geometry
return structured errors, naming the affected entities.

Requests specify the query (overlap, clearance threshold, or minimum distance),
finite work budget, numerical policy and requested evidence. Let c >= 0 be the
required clearance in model units, and d the minimum distance between the closed
solid sets. Interior overlap is classified separately: unsigned d is zero for
both contact and overlap and is not penetration depth.

The following are semantic outcomes, not a commitment to a particular enum ABI:

| Outcome | Required evidence |
| --- | --- |
| Clear | Solids are disjoint and a certified lower bound establishes d > c. With c = 0, this proves separation. |
| ClearanceViolation | Solids are disjoint, d > 0 is established, and a certified upper bound establishes d <= c. No physical contact is implied. |
| Contact | Boundaries meet and interiors do not overlap, established by supported predicates/classification. A merely small distance is insufficient. |
| Interference | Solid interiors overlap, including full containment; a boundary crossing is not required. |
| Unresolved | Available evidence cannot establish the requested conclusion within the supported numerical/work limits. |

For positive clearance, Contact and Interference also fail the clearance
requirement. Equality d = c violates this contract; uncertainty spanning that
threshold is Unresolved unless another valid proof settles it. A zero-clearance
request must never infer Contact solely from a distance enclosure containing zero.

Distance results include an enclosure [dLower, dUpper], its metric and units,
validated witness points where available, and entity/instance references tied to
input shape revisions. The upper bound must correspond to feasible geometry;
a local stationary point does not establish a global minimum or lower bound.
Minimum-distance requests must meet their requested enclosure width, even if a
threshold-only query could already terminate. Contact normals at non-smooth
features may be non-unique; report that explicitly. Penetration depth is a
separate optional capability, not a fabricated negative unsigned distance.

Unresolved includes its reason, implicated entities and any valid partial bounds.
The default complete-query interface returns a structured failure for unresolved
work, consistent with ADR-0006 and ADR-0029. A diagnostic partial-result interface
must identify incompleteness explicitly. No timeout, unsupported leaf, exhausted
subdivision budget or non-convergence may become Clear. A proved interference
witness can finish an existence query early; claiming absence, all contacts or a
global minimum requires the corresponding complete evidence.

### 3. Conservative acceleration and geometric refinement

- Use bounding-volume hierarchies over objects and trimmed faces/patches. Bounds
  enclose geometry, numerical uncertainty and transforms; clearance-aware culling
  discards a pair only when separation beyond the requested threshold is proved.
- Prefer analytic envelopes and certified analytic queries where supported.
  Refine remaining candidates against the actual trimmed geometry, including
  edges, vertices, seams, poles and containment. Reuse shared classification and
  intersection services rather than creating inconsistent predicates.
- Meshes may accelerate queries when approximation bounds are certified and
  included in the distance/occupancy bounds. Ordinary display tessellation cannot
  certify clearance to the CAD model. Conservative overlap is only a candidate,
  not proof of interference.
- Cache acceleration data in an explicitly owned workspace keyed by shape
  revision and configuration. Reuse local-space geometry for rigid instances;
  invalidate affected bounds on geometry or stock changes. Follow ADR-0007 for
  parallel reductions, deterministic results and tie-breaking of witnesses.
- General Boolean construction is not a prerequisite for every query. Existing
  collision libraries may serve as optional adapters or differential references;
  no backend is selected here, and an adapter may not silently strengthen its
  backend's assurance or supported-motion claims.

### 4. Continuous verification

A continuous request specifies a closed time/parameter interval, both objects'
rigid motions and a documented interpolation law. The motion provider supplies
conservative swept bounds or certified displacement/speed bounds over subintervals,
including rotational displacement of all geometry about the actual pivot. A pose
callback without such bounds supports sampled diagnostics only.

Use conservative advancement or interval subdivision with valid distance and
motion bounds. Clear requires coverage of the entire requested interval: every
subinterval must be proved clear. Fixed-rate samples, adaptive samples without
bounds, and an assumed straight interpolation between machine poses cannot
establish this result. Discontinuities are explicit segment boundaries and must
not be silently interpolated across.

A witnessed threshold violation or interference may terminate an existence query.
A request for the first event additionally requires a certified event-time bracket
and proof that the preceding interval is clear. Initial contact/overlap is handled
explicitly. Report the supplied motion model, time units, verified interval and
remaining unresolved subintervals. Clearance crossing, first contact and first
interior overlap are distinct events. Budget exhaustion is Unresolved, never an
implicit proof that the move is safe.

### 5. Stock, tolerances and derivatives

ADR-0044 remains the owner of tri-dexel fields and their interval operations.
The application maintains the evolving stock; Stereon verifies queries against a
specified stock revision. Endpoint accuracy, transverse ray spacing and feature
preservation are separate error sources. A feature-detection heuristic or accurate
endpoints alone do not certify material between rays. Certified collision use
requires conservative material bounds over cells, adequate proved refinement, or
a geometric fallback; otherwise the query remains approximate or Unresolved.

Keep required machining clearance, allowed stock allowance, model tolerance,
algorithmic convergence limits, geometric approximation error and numerical error
bounds separate (ADR-0028). Account for each relevant uncertainty explicitly;
do not add an arbitrary epsilon or enlarge model tolerance to force a decision.
Physical machine accuracy, deflection and setup uncertainty belong to the
application's clearance/envelope policy. A geometric certificate is conditional
on its supplied geometry and motion, not a machine-safety certification.

Under ADR-0021, distance or clearance derivatives may be exposed only on a fixed,
regular closest-feature branch. Feature switches, non-unique minimizers and
contact events are reported as Nondifferentiable/NearEvent or NotSupported, rather
than differentiated through. Derivative support is not required for the first stage.

### 6. Delivery stages

1. **Static analytic subset:** holder/shank envelopes against static fixtures and
   protected part solids; solid overlap, containment and threshold clearance with
   documented bounded-distance support. Define the supported pairs before coding.
2. **Three-axis offline verification:** fixed-orientation straight translations,
   certified swept bounds and continuous clearance checks against static geometry.
3. **Broader geometry and rotation:** trimmed freeform faces, tool assemblies and
   five-axis rigid trajectories supplied by the application, with verified motion
   bounds and moving-fixture support.
4. **Evolving stock integration:** tri-dexel or other bounded stock queries,
   revision invalidation and application-defined cutting/contact policy.

Each stage must pass its verification gate before claiming support. These are
capability stages, not new version promises: the existing roadmap is unchanged.
The first prototype need not wait for unrestricted Boolean or general SSI support.

## Consequences

- Machining applications receive explainable geometric results without paying to
  construct intersection solids at every pose. Shared distance, containment and
  motion services also benefit assembly clearance checks.
- Conservative bounds and explicit unresolved cases may require refinement and
  may reject difficult queries. This is preferable to an unsupported clear result;
  usable completion rates must still be measured rather than hidden by failures.
- Implementation requires global-distance bounds, containment handling, acceleration
  structures, motion contracts and regression corpora beyond surface intersections.
- This proposal establishes no real-time guarantee, production qualification or
  speed advantage over OCCT/FCL. Those require scoped measurements at equal
  geometry, motion and assurance. No code or runtime behavior changes with this ADR.

## Verification

Before accepting this ADR, demonstrate stages 1 and 2 on a declared analytic
prototype domain and publish revision-bound evidence under ADR-0041. Later stages
have their own gates; prototype acceptance does not qualify unsupported domains.

- **Static analytic oracle:** separated, tangent, overlapping and fully contained
  solids, cavities and disconnected components. Clearance cases below, exactly at
  and above c. Known distances lie inside returned bounds; witnesses lie on the
  reported geometry within their bounds. No false Clear or false definitive event.
- **Adversarial geometry:** nearly tangent features, edges/vertices, seams, thin
  obstacles, large translations, mixed scales and invalid/non-manifold inputs.
  Every outcome is correct or explicitly unresolved. Each repaired defect adds a
  permanent regression; unresolved cases are counted, not dropped from evidence.
- **Motion oracle:** an obstacle crossed only between endpoint poses, a grazing
  move, initial overlap, threshold crossing and safe near misses. At stage 3 add
  rotation-only collisions at the holder tip, off-axis pivots and simultaneous
  fixture motion. A first-event bracket contains the oracle event and has a proved
  clear prefix. Include split-interval consistency and reversed-motion checks.
- **Completion and failure:** force work-budget exhaustion, unsupported geometry,
  missing motion bounds and stale workspace revisions. None becomes Clear. Every
  ordinary fixture in a stage's promised corpus must resolve within its declared
  budget; adversarial unresolved rates and reasons are reported separately.
- **Stock at stage 4:** features between all three ray families, local removal and
  stale stock bounds. No sampled-field result is upgraded to a CAD certificate
  without conservative occupancy evidence or a fallback.
- **Independent evidence:** analytic or high-precision bounded references are the
  correctness oracle. Differential checks against OCCT or FCL supplement them on
  matching supported domains; disagreement is investigated, not settled by voting.
- **Performance and determinism:** report cold construction and warm queries,
  throughput, latency distribution, peak memory, refinement work and unresolved
  rate at fixed clearance and error bounds. Record hardware, compiler, input and
  thread count under ADR-0043. Compare like-for-like static and continuous queries;
  repeat at multiple thread counts with identical declared results under ADR-0007.
