## Summary

<!-- What does this change do, and why? -->

Fixes #<!-- issue number, if any -->

## Libraries touched

<!-- Tick every library the change touches. The highest tier sets the review rules (ADR-0015). -->

- [ ] Tier A: `core`, `curve`/`surface` evaluation, `io`, `python`, `viewer`, `bench`, tests, docs, build/CI
- [ ] Tier B: `robust`, `topo`, `build`, `mesh` (the reviewer must be able to explain the algorithm; property tests required)
- [ ] Tier C: `intersect`, `boolean`, `blend` (line-by-line review and module-owner approval required)

## Tests

- [ ] This PR adds a test that **fails without the change and passes with it**: <!-- name the test(s) -->
- [ ] Property-based tests are included (required for Tier B and Tier C)
- [ ] Correctness is checked by an independent referee, not only by tests written in this PR: differential tests (OpenCASCADE/Manifold), analytic values or the corpus <!-- required for AI-assisted Tier B/C changes -->
- [ ] Any failure I could not fix is added to `corpus/` as a regression case with a replay journal, not skipped

## Checklist

- [ ] The change follows all accepted ADRs, or this PR adds a new ADR that supersedes the one it contradicts
- [ ] No literal epsilons or raw sign tests in `topo/`, `boolean/` or `mesh/` (ADR-0001, ADR-0002)
- [ ] No global or `thread_local` mutable state; parallel code has a determinism test (ADR-0007)
- [ ] Every new source file has the MPL-2.0 header (`docs/license-header.txt`)
- [ ] Any new dependency has an MPL-compatible licence, named here: <!-- dependency + licence -->
- [ ] Public behaviour changes are documented
- [ ] Every commit is signed off (`git commit -s`)

## AI disclosure (ADR-0015)

- [ ] AI tools drafted a significant part of this change

<!-- If ticked, say which parts and how you verified them. You are responsible for every line and must be able to explain it in review. -->
