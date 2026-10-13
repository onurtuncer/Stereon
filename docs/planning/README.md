# Stereon planning

**Owner:** Onur Tuncer  
**Audit date:** 2026-10-12  
**Publication review:** 2026-10-13  
**Next review:** after the topology prototype lands, or before the next phase gate.

- [Capability matrix and repository audit](stereon-capability-matrix.md)
- [Placement and follow-up checklist](capability-matrix-checklist.md)
- [Existing benchmark specification](../benchmarks.md)
- [Architecture decision records](../adr/README.md)

The capability matrix is a **proposed planning baseline**, not an ADR or a delivery commitment. Its audit is pinned to commit `cf1dbaa925fb9c60fc642ca6962f329af7d52047`; maturity grades apply only to the subsets supported by that evidence.

The R0–R5 stages describe the long-term objective of comprehensive OCCT competition. They do not rename the repository's existing release phases or silently accept proposed ADRs.

The checklist is retained as a follow-up planning document; unchecked boxes do not mean the associated capability is absent. Implementation assessments are in the matrix. OCCT comparison harnesses, dedicated acceptance-fixture directories and published comparison reports mentioned in the checklist are **planned additions**, not existing artifacts. The implemented microbenchmarks remain under `bench/`, and the current corpus remains under `corpus/`.

An exact OCCT benchmark commit and build configuration have not yet been established. The matrix's OCCT documentation version is a reference for scope, not a reproducible benchmark pin. Record those details with the comparison harness when it is implemented.


The matrix includes a publication update for main at `25effbdcf26d1f35e3e6ba3314bb54512b43b5bd`, distinguishing the now-passing CI from the original audit snapshot.
