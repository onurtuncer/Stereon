# ADR-0015: AI contribution policy

- **Status:** Proposed
- **Date:** 2026-09-28
- **Owner:** Onur Tuncer
- **Phase:** 0

## Context

Much of Stereon's code will be drafted by AI coding agents. They are effective on well-specified modules (NURBS algorithms, parsers, bindings, tests) and unreliable on correctness-critical, degenerate-case-heavy code (surface–surface intersection, booleans, blends). Without a policy, the project risks a large codebase that passes obvious tests, fails on hard cases, and that no human fully understands.

## Options

1. **No AI-written code.** Safest; gives up most of the productivity gain.
2. **AI code treated like any other contribution.** Fast; review depth varies and critical code may merge without deep understanding.
3. **Tiered policy by library, with tests as the specification.**

## Decision

Adopt option 3.

- **Tests first.** Every change, AI-written or not, includes a test that fails before the change and passes after it. Phase exit gates are written as executable tests before implementation starts.
- **Tiers.**
  - *Tier A* (`core`, `curve`, `surface` evaluation, `io`, `python`, `viewer`, `bench`, tests): AI-written code allowed; standard review.
  - *Tier B* (`robust`, `topo`, `build`, `mesh`): AI drafts allowed; the reviewer must be able to explain the algorithm, and property tests are required.
  - *Tier C* (`intersect`, `boolean`, `blend`): AI drafts allowed only as proposals; line-by-line human review, and the named module owner approves.
- **ADRs are human-owned.** Agents receive the relevant ADRs as context; a change that contradicts an accepted ADR is rejected.
- **Independent referees.** AI-written changes are judged by differential tests against OpenCASCADE and Manifold, analytic values and the corpus — never by tests written in the same change alone.
- **Disclosure.** Pull requests state whether AI tools drafted significant parts of the change (a checkbox in the PR template).
- **Unfixed failures** become minimal regression cases with replay journals, not skipped tests.

## Consequences

- The project gains most of the AI productivity benefit where it is safe.
- Human attention is concentrated on Tier C, where it matters most.
- Module owners become a bottleneck for Tier C; pair ownership mitigates this.

## Verification

Quarterly audit: sample merged Tier B/C changes and confirm tests, reviewer sign-off and owner approval. Track the corpus success rate per release; a drop blocks the release.
