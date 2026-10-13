# Stereon Capability Matrix Checklist

## Placement

- [ ] Place the capability matrix in the repository’s documentation or planning area, not in source code\.
- [ ] Use the recommended path:
  - [ ] `docs/planning/stereon-capability-matrix.md`
- [ ] If the repository has an established documentation structure, use the closest equivalent:
  - [ ] `docs/architecture/stereon-capability-matrix.md`
  - [ ] `docs/roadmap/stereon-capability-matrix.md`
  - [ ] `planning/capability-matrix.md`
- [ ] Confirm that the document is placed outside source code and follows the repository’s established documentation structure\.
- [ ] Treat the document as a proposed planning baseline, not an ADR\.
- [ ] Avoid placing it in `docs/adr/` unless the repository treats planning baselines as ADRs\.
- [ ] If the repository has a documentation index, add a link under planning, architecture, or project roadmap documentation\.

## Documentation and Repository Additions

### Capability Matrix and Planning Index

- [ ] Add the capability matrix:
  - [ ] `docs/planning/stereon-capability-matrix.md`
- [ ] Add or update the planning index:
  - [ ] `docs/planning/README.md`
  - [ ] Link to the capability matrix\.
  - [ ] Explain that the matrix is a proposed planning baseline\.
- [ ] Explicitly mark the document as a proposed planning baseline\.
- [ ] Add an owner to repository metadata or the planning index\.
- [ ] Add a review date to repository metadata or the planning index\.

### ADRs and Referenced Decisions

- [ ] Add separate ADRs for decisions referenced by the matrix:
  - [ ] `docs/adr/`
- [ ] Link existing ADRs from the matrix instead of duplicating their content\.
- [ ] Link to existing ADRs rather than duplicating their content\.

### Benchmark Documentation

- [ ] Add published benchmark comparison reports:
  - [ ] `docs/benchmarks/`
- [ ] Record the benchmark configuration separately once it has been established\.

## Infrastructure

### OCCT Benchmark Infrastructure

- [ ] Add OCCT benchmark infrastructure:
  - [ ] `benchmarks/occt/`
- [ ] Record the pinned OCCT version or commit\.
- [ ] Record the pinned OCCT commit separately once it has been established\.
- [ ] Record the build configuration\.
- [ ] Add the comparison harness\.

### Acceptance Fixtures and Corpus Metadata

- [ ] Add acceptance fixtures and corpus metadata:
  - [ ] `tests/fixtures/capability-matrix/`

## Validation and Readiness

### Capability Matrix Status

- [ ] Keep every capability\-matrix row at `U / —` until repository evidence has been reviewed\.

### Repository References and Planned Additions

- [ ] Verify that all referenced paths, ADRs, fixtures, benchmarks, and reports exist or are clearly marked as planned additions\.
- [ ] Confirm that links to existing ADRs resolve correctly\.
- [ ] Confirm that links from the planning index to the capability matrix resolve correctly\.
- [ ] Confirm that benchmark reports, fixtures, corpus metadata, and OCCT infrastructure are either present or explicitly identified as planned additions\.

### Final Review

- [ ] Confirm that the document is explicitly identified as a proposed planning baseline\.
- [ ] Confirm that ownership and review\-date information is recorded in repository metadata or the planning index\.
- [ ] Confirm that the pinned OCCT version or commit and benchmark configuration are recorded separately once established\.
- [ ] Confirm that all capability\-matrix rows remain at `U / —` until repository evidence has been reviewed\.
