# Geometric assurance inventory

**Status:** Initial review queue; no current-release certification asserted.  
**Established:** 2026-11-11. **Basis:** the [v0.2 acceptance ledger](v0.2-acceptance-ledger.md), last marked updated 2026-11-10 to main at #249.  
**Authority:** [ADR-0055](../adr/0055-geometric-assurance-profiles.md), [ADR-0029](../adr/0029-intersection-and-classification-outcomes.md), [ADR-0047](../adr/0047-numerical-operation-completion-and-publication.md).

## How to read this inventory

This is an **audit worklist**, not an algorithm capability matrix or release report. The descriptions below are paraphrases of implementation claims in the existing v0.2 ledger, not a fresh line-by-line source/CI assessment. Every row therefore starts **Pending review** until the domain, specific witness, test result and exact revision have been independently checked.

Do not mark a row Verified merely because a handler exists or a regression test was checked in. The audit distinguishes:

- **Planned:** an ADR, design or roadmap requests behavior, but implementation is not established.
- **Implemented / review pending:** implementation is documented or observed; the full assurance claim is not independently qualified.
- **Verified (revision-bound):** supported domain, direct oracle, numerical/semantic evidence, test run and reviewer are recorded.
- **Unresolved / limited:** a case or part of the declared domain cannot satisfy the requested guarantee; document explicit failure semantics.
- **Not applicable:** the assurance dimension is not meaningful for that operation.

Each assurance dimension must be recorded independently: decision validity, geometric accuracy, completeness, topology, event/derivative regularity and provenance.

## v0.2 baseline: families to audit

| Family | Existing reported scope (not a new certification claim) | Unresolved assurance question / next evidence | Current inventory status |
| --- | --- | --- | --- |
| 2D curve–curve | Exact or interval-based line, conic and B-spline handlers; tangency and overlap behavior varies by pair | Close spline/spline tangencies, nonidentical overlaps, knot-near crossings; reviewed proved-empty/complete-domain cases and comparator tests | Implemented / review pending |
| 3D curve–curve | Line/line resolution semantics; spline and conic closest-meeting searches, some with Krawczyk checks | Verify boundary minima, tangency fail-closed contracts, exact-vs-resolution semantics and independent OCCT/analytic corpus | Implemented / review pending |
| Curve–surface | Analytic quadrics/ellipsoid and spline/analytic handlers; spline-surface candidate isolation for regular crossings | Spline-surface tangencies, periodic cases, supported bounds and complete-empty cases, pinned differential results | Implemented / review pending |
| Analytic surface–surface | Elementary pairs and analytic special cases, including overlapping and tangential configurations | Case-by-case decision status, independent full-domain evidence and contact graph consistency | Implemented / review pending |
| Generic B-spline SSI | Interval-enclosed patch subdivision, normal-cone pruning, tracing and unresolved failures described in ledger | Verify global branch completeness, bounded trace error, singular/tangential behavior, coverage budgets and small-loop fixtures | Implemented / review pending |
| Intersection-to-Boolean integration | Analytic small-solid prototype and Boolean intersection consumption | Classifier/assembly proof, split topology validation, lineage, invalid-success count over a declared corpus | Implemented / review pending |
| Tolerant intersection evidence | Steps 0.2-T1/T2 tied to ADR-0050; typed tolerance causes and reported fit/gap evidence | Verify every relevant result and direct-tolerance-write prohibition at current head; link CI | Implemented / review pending |
| Intersection sensitivities | ADR-0021 stage-4 linearisation is in the cross-cutting v0.2 ledger | Verify per-handler regularity, finite-difference checks, nondifferentiable outcomes and current coverage | Implemented / review pending |
| Topological event witnesses | ADR-0054 Stage A describes analytic contact and critical-parameter cases | Implement and qualify candidate/certified/unresolved classifications with reference event sweeps | Planned |
| Metamorphic regression promotion | ADR-0053 and [backfill ledger](intersection-metamorphic-backfill.md) identify historical repairs | Close each confirmed defect with red-on-defect/green-on-repair and independent semantic oracle | In progress / not independently qualified |

These rows are intentionally broad **inventory seeds**, not a claim that each family has a uniform assurance level. The responsible implementation reviewer must split a row into supported geometry pairs, domains or algorithm paths when guarantees differ.

## Required per-family record

Copy the following fields into a reviewed family section or machine-readable evidence artifact, without changing the published API contract:

| Field | Entry required |
| --- | --- |
| Family / handler / version | Operation name, geometry pair, supported parameter/trim domain; product/release stage |
| Operation and architecture authority | Source/API path and relevant ADRs; exact implementation revision |
| Decision assurance | Exact predicate, rigorous enclosure, justified numerical decision, heuristic, unresolved, or not applicable |
| Construction assurance | Metric, units, domain, rigorous bound vs estimate/sample; source of error |
| Completeness | Search domain, exclusions, proved-pruning/separation rules, work budget and unresolved cases |
| Topology/lineage | What is validated (if anything), owner-qualified entity refs, topology validation and publication evidence |
| Regularity/events | Smooth-branch assumptions, event status, derivative-status contract if supported |
| Representation boundary | Authoritative B-rep/analytic/spline input, any derived discretization, conversion error and validity of mapping |
| Tests and independent oracle | Corpus IDs, analytic/reference cases, metamorphic relation IDs, repaired defects, differential comparator |
| CI/reviewer evidence | Exact commit, CI run URL/result, environment, review decision and date |
| Performance evidence | Comparable error/coverage policy, workload, build flags and measured time/memory/failure rate |
| Publication status | Planned, implemented/review pending, verified, unresolved/limited or not applicable **per dimension** |

A profile may be incomplete during development. Missing evidence MUST remain visible rather than becoming a silent successful gate.

## Cross-family audit questions

1. **Correctly empty?** What proves an empty result over the complete requested domain, rather than only at sampled seeds?
2. **Correctly constructed?** Is the geometric error bounded over the entire returned object or only sampled?
3. **Correct topology?** For a Boolean or B-rep construction, does a separate validator establish incidence/orientation and the intended solid region?
4. **Correct status?** Is an inconclusive tangency, near-coincidence or parameter singularity reported as such?
5. **Correct provenance?** Are input snapshots, context/tolerances, affected entities and representation changes traceable?
6. **Correct transformation?** Do appropriate rigid motion, operand swap, scale, parameter reversal and partition metamorphic relations preserve semantic results with justified error budgets?
7. **Correct comparison?** Are any OCCT/CGAL comparisons for actually equivalent operations and numerical policies, using pinned builds and independent correctness checks?

## Gate evidence placeholders

**v0.2 / 0.2-CG1:** Audit records created per supported intersection family; no verified row without a CI run, reviewer, implementation revision and executable oracles.  
**v0.2 / 0.2-CG2:** ADR-0053 backfill and ADR-0054 initial analytical event evidence linked where supported.  
**v0.3 / 0.3-CG3:** Topology query and event evidence added with lineage snapshots.  
**v0.4–v0.6:** Boolean pipeline assurances and independent topological validation added to the supported-domain release reports.  
**v0.7–v0.9:** Machining/meshing/exchange and sensitivities added to the same assurance scheme.

The gates and accepted release criteria live in the [primary roadmap](roadmap-and-versioning.md); this inventory is the evidence index, not a replacement for acceptance ledgers.
