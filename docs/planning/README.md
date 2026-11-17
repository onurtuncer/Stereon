# Stereon planning

**Owner:** Onur Tuncer  
**Audit date:** 2026-10-12  
**Publication review:** 2026-10-13  
**Next review:** after the topology prototype lands, or before the next phase gate.

- [Roadmap, modelling milestones and versioning](roadmap-and-versioning.md)
- [Computational geometry convergence roadmap](geometric-computing-roadmap.md): workstreams, dependency order, staged release evidence and research scope
- [Geometric assurance inventory](geometric-assurance-inventory.md): operation-family assurance profiles and revision-bound audit worklist
- [Capability matrix and repository audit](stereon-capability-matrix.md)
- [v0.1 acceptance ledger](v0.1-acceptance-ledger.md): gates, evidence and the OCCT operations still to add
- [v0.1 release report](v0.1-release-report.md): the release-wide record the roadmap requires (commit, corpus, comparator, limits, structured failures, invalid successes), final since the `v0.1.0` tag of 2026-11-03
- [v0.1 release audit of 2026-10-30](v0.1-audit-2026-10-30.md): the read-only audit of `fb7366a` that held the tag, its two numerical defects and their disposition
- [Boolean prototype plan](boolean-prototype-plan.md): the last Phase 0 gate (ADR-0001, ADR-0007), its milestones and the decisions of 2026-10-20
- [Curved Booleans](curved-booleans.md): design of cones, spheres, tori and ellipsoids in Union, Common and Cut, and of coincident curved faces, in slices (v0.2)
- [Placement and follow-up checklist](capability-matrix-checklist.md)
- [Existing benchmark specification](../benchmarks.md)
- [Architecture decision records](../adr/README.md)

The capability matrix is a **proposed planning baseline**, not an ADR or a delivery commitment. Its audit is pinned to commit `cf1dbaa925fb9c60fc642ca6962f329af7d52047`; maturity grades apply only to the subsets supported by that evidence.

The computational-geometry overlay and proposed ADR-0055 do not alter release numbering, automatically accept other proposed ADRs, or establish previously unverified numerical guarantees.

The R0–R5 stages describe the long-term objective of comprehensive OCCT competition. They do not rename the repository's existing release phases or silently accept proposed ADRs.

The checklist is retained as a follow-up planning document; unchecked boxes do not mean the associated capability is absent. Implementation assessments are in the matrix. OCCT comparison harnesses, dedicated acceptance-fixture directories and published comparison reports mentioned in the checklist are **planned additions**, not existing artifacts. The implemented microbenchmarks remain under `bench/`, and the current corpus remains under `corpus/`.

An exact OCCT benchmark commit and build configuration have not yet been established. The matrix's OCCT documentation version is a reference for scope, not a reproducible benchmark pin. Record those details with the comparison harness when it is implemented.


The matrix includes a publication update for main at `25effbdcf26d1f35e3e6ba3314bb54512b43b5bd`, distinguishing the now-passing CI from the original audit snapshot.

The [14 October roadmap integration](roadmap-and-versioning.md) maps Stroud-derived M1–M6 checkpoints to the existing release sequence, adds exit evidence and distinguishes product, native-format and C ABI versions. It does not refresh the historical source audit or accept proposed ADRs.
