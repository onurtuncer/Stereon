# Contributing to Stereon

Thank you for your interest in Stereon. This document explains how to propose changes and what a pull request needs before it can be merged.

> **Status: pre-alpha.** Stereon is in Phase 0 (foundations). The build system, CI and core libraries are being set up, and APIs change without notice. Small fixes and discussion are welcome now; larger contributions are easiest once the Phase 0 foundations settle. Please open an issue before starting anything substantial.

## Ground rules

1. **Every change includes a test that fails before it and passes after it.** No exceptions for small changes, refactors that fix behaviour, or AI-assisted work.
2. **Accepted ADRs are binding.** Read the relevant records in [`docs/adr/`](docs/adr/) before you start. A pull request that contradicts an accepted ADR is rejected unless it comes with a new ADR that supersedes it.
3. **Failures are never skipped.** A case you cannot fix becomes a minimal regression case in `corpus/` with a replay journal. It is not marked as skipped or deleted.
4. **Be respectful.** Everyone taking part follows the [Code of Conduct](CODE_OF_CONDUCT.md).

## Before you start

- **Bugs:** open an issue using the bug template. Attach the smallest input that reproduces the problem (a native file or an operation journal if you have one) and the exact error or invalid result.
- **Features and design changes:** open an issue first. Anything that touches tolerances, topology storage, geometry dispatch, error handling, threading or the file format probably needs an ADR. Start from [`docs/adr/0000-template.md`](docs/adr/0000-template.md).
- **Security issues:** do not open a public issue. Follow [SECURITY.md](SECURITY.md).

## Building and testing

You need **GCC 16 or newer** (the reference compiler), CMake 3.30+ and Ninja. Clang and MSVC are not yet supported for development. They run in CI as non-blocking jobs until they support C++26 reflection and contracts ([ADR-0016](docs/adr/0016-cpp26-adoption-and-contracts.md)).

```bash
cmake --preset gcc16-debug
cmake --build --preset gcc16-debug
ctest --preset gcc16-debug
```

Use `gcc16-debug` while developing. It enforces contracts, enables checked handles (`STEREON_CHECKED_HANDLES`) and runs the validity checker after every public operation. Before opening a pull request, also run `gcc16-asan`. Run `gcc16-tsan` too if your change touches anything parallel.

On Windows, build in an MSYS2 UCRT64 shell with `mingw-w64-ucrt-x86_64-gcc`.

## What a pull request needs

- [ ] A test that fails without the change and passes with it.
- [ ] Property-based tests for new algorithms in `robust/`, `topo/`, `build/`, `mesh/`, `intersect/`, `boolean/` and `blend/`.
- [ ] All blocking CI jobs green (Linux, Windows and macOS on GCC 16, plus the sanitizer jobs).
- [ ] No new contract violations in the observe-mode corpus run.
- [ ] The MPL-2.0 header on every new source file (see [Licence headers](#licence-headers)).
- [ ] Documentation updated if public behaviour changes.
- [ ] Every commit signed off (see [Developer Certificate of Origin](#developer-certificate-of-origin)).
- [ ] The AI-disclosure checkbox in the pull request template filled in honestly.

Keep pull requests focused. One logical change per pull request is much easier to review than a mixed bag, and Tier C changes in particular get line-by-line review.

## Coding standards

Stereon is written in C++26 with C++20 modules. These rules come from the ADRs and are checked in review and, where possible, by tooling.

**Language and style**

- Format with the repository's `.clang-format` and keep `.clang-tidy` clean.
- Naming follows the [Hazel engine](https://github.com/TheCherno/Hazel) style ([ADR-0017](docs/adr/0017-naming-conventions.md)):

  | Element | Style | Example |
  | --- | --- | --- |
  | Namespaces, modules | PascalCase | `Stereon`, `Stereon.Core` |
  | Types, concepts, enums, enum values | PascalCase | `KernelError`, `CurveGeometry`, `ErrorCode::ToleranceExceeded` |
  | Functions and methods | PascalCase; accessors start with `Get` | `MakeBox`, `Orient3D`, `ctx.GetResolution()` |
  | Public struct fields, constants | PascalCase | `Version::Major`, `LibraryVersion` |
  | Private and protected members | `m_` + PascalCase | `m_Pages` |
  | Static members and file-scope statics | `s_` + PascalCase | `s_NextStoreId` |
  | Local variables, parameters | camelCase | `faceCount`, `tolerance` |
  | Macros | `STEREON_` + UPPER_SNAKE_CASE | `STEREON_TRY` |
  | Source files | PascalCase | `Version.cppm`, `VersionTests.cpp` |

  Library folders and CMake targets stay lowercase (`src/core`, `stereon::core`).
- Geometry types are plain aggregates with no base class. Their behaviour lives in free functions found by ADL: `Eval`, `Deriv`, `Domain`, `BoundingBox`, and so on ([ADR-0003](docs/adr/0003-geometry-polymorphism.md)).
- Prefer value semantics. Avoid owning raw pointers and `shared_ptr` graphs between topological entities.

**Correctness**

- **No exceptions.** Operations that can fail return `std::expected<T, KernelError>` marked `[[nodiscard]]`. `std::bad_alloc` is the only exception that may propagate ([ADR-0006](docs/adr/0006-error-handling.md)).
- **Contracts for programming errors.** Preconditions (parameter in domain, non-degenerate input) and key postconditions are written as contracts. A contract violation never becomes a `KernelError` ([ADR-0016](docs/adr/0016-cpp26-adoption-and-contracts.md)).
- **No invalid results.** An operation returns a shape that passes the validity checker or an error. Returning an invalid shape is a bug.
- **No hidden epsilons.** Code in `topo/`, `boolean/` and `mesh/` never compares against literal constants. It uses `ToleranceOf(entity)` or `ctx.GetResolution()` ([ADR-0001](docs/adr/0001-tolerance-model.md)).
- **Predicates go through `robust/`.** Every topological decision in `topo/`, `boolean/` and `mesh/` uses a predicate from `stereon::robust`. Direct sign tests on computed doubles are not allowed there ([ADR-0002](docs/adr/0002-numeric-robustness.md)).
- **Tolerances grow only explicitly.** Increasing an edge or vertex tolerance happens through a named operation that logs the entity, the old and new values and the reason.

**State and threading**

- **No global mutable state** and no `thread_local` state that can affect results. Caches, schedulers, tolerances and logging live in the `Context` passed to every operation ([ADR-0007](docs/adr/0007-threading-and-determinism.md)).
- **Determinism.** Results must be bit-identical for any thread count. Merge parallel results in a fixed order keyed by input entity index, never in completion order, and use fixed-order or compensated floating-point reductions. Every parallel algorithm needs a determinism test.
- **Handles stay typed.** Use `Id<Tag>` handles and their named accessors (`IsReversed(h)`, `BaseOf(h)`). Never convert handles to raw integers or mask bits by hand ([ADR-0004](docs/adr/0004-topology-storage-and-handles.md)).

**Library layering**

Libraries depend only on libraries below them: `core` → `robust` → `curve` / `surface` → `intersect` → `topo` → `build` → `boolean` → `blend`, with `mesh` and `io` on top. CI fails on an upward dependency.

## Dependencies

New third-party dependencies must use a licence compatible with MPL-2.0: MIT, BSD, Apache-2.0, MPL-2.0 or Boost. Copyleft dependencies are only allowed as optional, separately built adapters ([ADR-0008](docs/adr/0008-licence.md)). Name the dependency and its licence in the pull request description.

Test corpus files keep their original licences. Record the source and licence of every file you add in `corpus/LICENSES.md`, and only add files you have the right to redistribute.

## AI-assisted contributions

AI-assisted contributions are welcome under the same rules as any other change, plus the tiered review in [ADR-0015](docs/adr/0015-ai-contribution-policy.md):

| Tier | Libraries | What is required |
| --- | --- | --- |
| A | `core`, `curve` and `surface` evaluation, `io`, `python`, `viewer`, `bench`, tests | Standard review |
| B | `robust`, `topo`, `build`, `mesh` | The reviewer must be able to explain the algorithm; property tests are required |
| C | `intersect`, `boolean`, `blend` | AI drafts are proposals only; line-by-line human review and approval from the module owner listed in `CODEOWNERS` |

Further rules:

- Say in the pull request whether AI tools drafted a significant part of the change.
- Tests written in the same change are not enough to prove an AI-written change correct. It must also pass the independent referees: differential tests against OpenCASCADE and Manifold, analytic values and the corpus.
- You are responsible for every line you submit, however it was produced. You must understand it well enough to answer review questions.
- ADRs are written and owned by humans.

## Commit conventions

- Subject line: `<library>: <imperative summary>`, at most 72 characters, no trailing period. Examples: `core: add Interval with outward rounding`, `robust: fix Orient3D filter bound`, `docs: clarify tolerance ceiling in ADR-0001`.
- Use `build`, `ci`, `docs` or `tests` as the prefix for changes outside `src/`.
- After a blank line, the body explains *why* the change is made. Reference issues and ADRs (`Fixes #42`, `See ADR-0004`).
- Keep each commit buildable. Branch protection on `main` requires a linear history, so rebase rather than merge.

## Developer Certificate of Origin

Stereon uses the [Developer Certificate of Origin 1.1](https://developercertificate.org/) (DCO) instead of a contributor licence agreement. By signing off a commit, you certify that you wrote the change or otherwise have the right to submit it under the project's licence (MPL-2.0).

Sign off every commit with `git commit -s`. This adds a line like the following, using your real name and an email you can be reached at:

```
Signed-off-by: Jane Doe <jane@example.com>
```

Pull requests with unsigned commits cannot be merged. To sign off commits you have already made, run `git rebase --signoff main` and force-push your branch.

## Licence headers

Stereon is licensed under the [Mozilla Public License 2.0](LICENSE). Every source file starts with this header, using the comment syntax of its language. [`docs/license-header.txt`](docs/license-header.txt) has copy-ready versions for C++, CMake, Python and other languages:

```cpp
// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
```

Files that cannot hold a comment (Markdown, JSON, images) are covered by [`REUSE.toml`](REUSE.toml). `reuse lint` checks that every file is covered, and will run in CI.

## Review and merging

- `main` is protected. Every pull request needs one approving review and green blocking CI. Tier C changes also need approval from the module owner.
- Reviewers may ask for smaller pull requests, more tests or an ADR. This is normal, not a rejection.
- Maintainers squash or rebase on merge to keep the history linear.

## Questions

For questions about the design or where to start, open a GitHub Discussion or an issue labelled `question`.
