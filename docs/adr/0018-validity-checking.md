# ADR-0018: Validity checking

- **Status:** Proposed
- **Date:** 2026-10-10
- **Owner:** Onur Tuncer
- **Phase:** 0

## Context

Several decisions already lean on "the validity checker" without defining it. ADR-0001 relies on it for the containment invariant, ADR-0004 runs it after every operation in debug builds, ADR-0006 promises that an operation returns either a shape that passes it or an error, and ADR-0009 uses it to verify edge certificates. The checker is therefore one of the most important components in the kernel: it decides what may leave an operation. It needs a definition of what it checks, at what cost, and how it stays independent of the algorithms it judges, so that a shared wrong assumption does not approve a wrong result.

## Options

1. **Assertions inside each algorithm.** No separate component; each algorithm checks what it happens to know about. Coverage is uneven and the checks share the algorithm's assumptions.
2. **One monolithic `IsValid()`.** Simple to call; too slow to run after every operation, and one boolean answer does not help debugging.
3. **A levelled checker with its own traversal and predicates.** Cheap structural checks run always; geometric and solid-level checks run per policy. Findings name entities and rules.

## Decision

Adopt option 3.

### Levels

`Check(shape, level, ctx)` returns a `ValidityReport`. Levels are cumulative.

| Level | Checks | Cost |
| --- | --- | --- |
| `Structural` | Every index is in range and refers to a live slot (ADR-0004). `Next`/`Prev` of coedges are mutual inverses; every loop closes. Every radial ring closes and all its coedges reference the same edge. Every coedge belongs to exactly one loop; every loop to one face; every face to at most one shell. Edge kind and coedge count match the table in ADR-0009; the primary coedge is in the edge's ring. Shell faces are edge-connected. | O(n), no geometry evaluation |
| `Geometric` | Containment invariant (ADR-0001): vertex spheres contain edge ends; edge tubes contain the certified 3D curve and the images of all pcurves. Every `Intersection` and `Trim` edge has certificates and they are not stale. Pcurves lie within their surface's domain; each loop is closed in parameter space (seams handled). Tolerances are ≥ resolution and ≤ ceiling. Adjacent faces on a `Trim` edge reference the same surface; on an `Intersection` edge, different surfaces. | O(n) surface evaluations |
| `Solid` | Each closed shell has exactly two coedges per edge with opposite sense. Face orientation is consistent across every edge; the outer shell's normals point outward and void shells' inward (signed volume). Shells are nested correctly. No face–face self-intersection within a shell except along shared edges (uses ADR-0010). Volume is finite and non-zero. | up to O(n log n) with SSI |

### Rules

- **Independence.** The checker lives in `topo/` and uses only `robust/` predicates, direct geometry evaluation and its own traversal code. It must not use the adjacency caches, BVHs, classification results or intermediate data of the operation under check. The layering check (`tools/check_layering.py`) already prevents it from depending on `build/`, `boolean/` or `blend/`.
- **Findings, not booleans.** A `ValidityReport` lists findings, each with a stable rule ID (e.g. `S03 RadialRingOpen`, `G02 EdgeTubeViolated`), the entities involved via `Describe()`, the measured value and the limit. `IsValid()` is a convenience that checks for the absence of error-severity findings.
- **When it runs.** Debug, ASan and fuzzing presets run `Structural` and `Geometric` after every public operation and `Solid` after every operation that returns a solid (ADR-0004 safeguard 5). Release builds run `Structural` after every operation; `Geometric` and `Solid` run according to the `Context` policy, default `Geometric` on and `Solid` off. The ADR-0006 guarantee "returns a valid shape or an error" is verified by the CI corpus run with all levels on; in release it is a promise backed by that run plus the checks the policy enables.
- **Error mapping.** An input shape that fails a required level gives `KernelError::InvalidInput` with the report attached. A result that fails a post-check gives `KernelError::InvalidResult` with the report attached, records a journal entry for replay (ADR-0004), and is counted as a kernel bug in CI.
- **Checker bugs count double.** A false negative (approving an invalid shape) or a false positive (rejecting a valid one) is fixed before any other bug of the same severity, and every such case joins the checker's corpus.

## Consequences

- Every other ADR can now say precisely what "valid" means by naming a level and a rule.
- Operations become slower in debug and CI builds; acceptable, and release pays only for the `Structural` level by default.
- The `Solid` level depends on surface–surface intersection, so a self-intersection check is only as good as ADR-0010; the geometric checks (certificates, containment) do not depend on it.
- A second implementation of loop and ring traversal must be maintained for independence; it doubles as the reference for the production traversal in property tests.

## Verification

- A seeded-fault corpus: flipped face, open shell, missing coedge in a ring, coedge in two loops, pcurve outside its domain, `Trim` edge between faces on different surfaces, `Intersection` edge with a stale certificate, edge tolerance below the measured gap, self-intersecting shell, void shell outside the outer shell. Pass: each fault is found by the intended rule at the intended level, 100% of the time.
- Zero findings on the valid corpus (no false positives).
- `Structural` costs ≤ 5% of a boolean's run time on the 10⁶-face body; `Geometric` ≤ 50%.
- The layering check confirms the checker imports nothing from `build/`, `boolean/` or `blend/`.
- The end-to-end example (ADR-0014): a duct with inlet, outlet and wall tags, cut by a cylinder, produces a `Solid`-level clean report before meshing.
