# ADR-0017: Naming conventions

- **Status:** Proposed
- **Date:** 2026-10-08
- **Owner:** Onur Tuncer
- **Phase:** 0

## Context

Stereon's public API is still almost empty, so a naming convention can be chosen now at almost no cost. Within a few weeks `core`, `robust` and the geometry scaffolding will add hundreds of public names, and the ADRs already name functions (`Orient3D`, `FacesOf`, `VisitBuiltin`). Changing the convention later would mean renaming the whole API and breaking every user. Naming also needs to be checkable by tooling rather than argued in review.

## Options

1. **Standard-library style** (snake_case everywhere, as in `std::` and Boost). Matches the standard library that Stereon code calls constantly; types, functions and variables look alike.
2. **Google C++ style** (PascalCase types and functions, snake_case variables, `kConstant`, trailing `_` on members). Widely documented; the `k` prefix and trailing underscores are less readable in long geometric expressions.
3. **Hazel engine style** ([TheCherno/Hazel](https://github.com/TheCherno/Hazel)): PascalCase for namespaces, types, functions and public fields; `m_`/`s_` prefixes for private and static members; camelCase for locals and parameters. Types, functions and members are distinct at a glance, and the convention is familiar to many game and graphics programmers.

## Decision

Adopt option 3, the Hazel style, for all Stereon C++ code.

| Element | Style | Example |
| --- | --- | --- |
| Namespaces | PascalCase | `Stereon`, `Stereon::Exec` |
| Module names | PascalCase, matching the namespace | `Stereon.Core`, partition `Stereon.Core:Version` |
| Types, concepts, enums, enum values | PascalCase | `KernelError`, `CurveGeometry`, `ErrorCode::ToleranceExceeded` |
| Free functions and methods | PascalCase | `MakeBox`, `Orient3D`, `Describe` |
| Accessors | `Get` prefix | `ctx.GetResolution()`, `error.GetMessage()` |
| Public struct fields | PascalCase | `Version::Major`, coedge `Next` |
| Constants (`constexpr`, `const` globals) | PascalCase | `LibraryVersion`, `PageSize` |
| Private and protected members | `m_` + PascalCase | `m_Pages` |
| Static members and file-scope statics | `s_` + PascalCase | `s_NextStoreId` |
| Local variables, function parameters | camelCase | `faceCount`, `tolerance` |
| Macros | `STEREON_` + UPPER_SNAKE_CASE | `STEREON_TRY`, `STEREON_CHECKED_HANDLES` |
| Source files | PascalCase | `Version.cppm`, `VersionTests.cpp` |

Further rules:

- **Name clashes between a function and a type are avoided by the function name.** Free functions that return a modified copy use past participles (`Transformed(curve, t)`, `Reversed(curve)`), so they do not collide with types such as `Transform`. Predicates on handles use `Is` (`IsReversed(h)`). Relationship queries use `Of` (`FacesOf(edge)`, `ToleranceOf(entity)`, `BaseOf(h)`).
- **Acronyms and dimensions** are written as words: `Orient2D`, `InSphere`, `BSplineCurve`, `PCurve`, `Id`.
- **Not covered by this ADR:** standard-library names and customisation points keep their own spelling (`begin`, `end`, `size`, `operator==`, `std::expected::transform`). Library folders and CMake targets stay lowercase (`src/core`, `stereon::core`), since they are not C++ identifiers.
- **Formatting is separate.** This ADR fixes names only; layout is set by `.clang-format`.

## Consequences

- Types, functions and members are distinguishable without tooling; `m_` and `s_` make object and global state visible in algorithm code.
- Stereon code reads differently from the standard library it calls. Accepted: the boundary is visible, which is useful in a kernel that wraps `std::` heavily.
- Names in ADRs 0001–0016 and the README were renamed to this style while they were still Proposed.
- Python bindings (nanobind) will need a rule for exposing PascalCase names: keep them, or map to snake_case for Python users. To be decided with the bindings.
- Range-based `for` and other standard customisation points must keep lowercase member names (`begin`, `end`), an exception reviewers must allow.

## Verification

`.clang-tidy` enforces the table above with `readability-identifier-naming`. Pass: a CI clang-tidy job reports zero naming violations on `src/` and `tests/`. Until that job exists, review checks names against this ADR.
