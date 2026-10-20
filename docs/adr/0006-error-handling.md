# ADR-0006: Error handling

- **Status:** Accepted
- **Date:** 2026-09-28
- **Revised:** 2026-10-10 (`OperationResult`; boundary between errors and contract violations)
- **Revised:** 2026-10-14 (outcome, lineage and metadata extension references)
- **Revised:** 2026-10-20 (two tiers for a handle from another shape; `Context` modifiers; the entity and report fields of `KernelError` are pending)
- **Accepted:** 2026-10-20 by Onur Tuncer: `Result` and `OperationResult` are built and used by every fallible operation, the kernel does not throw, and publication refuses an invalid shape with `InvalidResult` (see *Status of the verification*)
- **Owner:** Onur Tuncer
- **Phase:** 0

## Context

Kernel operations fail for legitimate reasons: a fillet radius too large, a boolean with an unresolvable degeneracy, a corrupt STEP file. Callers need to know what failed and where, without the kernel ever returning a silently invalid shape. Bindings (C, Python, C#) also need a uniform error model.

## Options

1. **Exceptions.** Idiomatic C++; invisible in signatures, costly on hot failure paths, awkward across C and language boundaries.
2. **Error codes with out-parameters.** Portable; loses context and is easy to ignore.
3. **`std::expected<T, KernelError>` with structured diagnostics.** Failure is explicit in the signature and carries rich context.

## Decision

Adopt option 3.

- Every public operation that can fail returns `std::expected<T, KernelError>`.
- **Shape-producing operations return `OperationResult`.** Booleans, sweeps, sewing, healing, import and every other operation that creates or changes topology return `std::expected<OperationResult, KernelError>`. `OperationResult` holds the `Shape`, the `History` (ADR-0011), the list of tolerance changes (entity, old value, new value, reason; ADR-0001) and structured warnings (`Diagnostic` with a code, message and entities). Success with a grown tolerance or a dropped attribute is thus visible to the caller, not only in a log, and the application can show it or reject the result. The `Context` can promote chosen warning codes to errors (for example "any tolerance growth fails"). Pure queries (mass properties, distance, classification of a point) return `std::expected<T, KernelError>` directly.
- `KernelError` holds: an `ErrorCode` enum (e.g. `InvalidInput`, `InvalidResult`, `ToleranceExceeded`, `IntersectionFailed`, `CertificationFailed`, `DegenerateResult`, `NotSupported`, `NotImplemented`, `ParseError`); a human-readable message; the entities involved (handles into the input shapes); the tolerance values at stake; an attached `ValidityReport` where one exists (ADR-0018); and an optional journal ID for replay.
- The kernel does not throw. `std::bad_alloc` may still propagate.
- **Errors versus contract violations.** Everything that enters the kernel from outside — user shapes, numeric parameters of public operations, files, geometry built from imported data — is *data* and is validated; invalid data yields `KernelError` (`InvalidInput`, `ParseError`, `NotSupported`), never a contract violation. Contracts (ADR-0016) guard what the kernel's own code is responsible for: handle validity and store membership, documented parameter domains of low-level evaluators, and internal invariants. A corrupt STEP file is a `ParseError`; a boolean on a shape that fails the validity checker is `InvalidInput`. A handle from another shape, or a stale or dead one, is checked in two tiers: an operation or a query that returns a `Result` (a traversal, a builder edit) refuses it with `InvalidInput`, which costs nothing next to its work, while an O(1) accessor that returns a value (`PositionOf`, `CurveOf`, ...) treats it as a contract violation, checked in the checked-handle builds, because a branch and a `Result` on every field read would cost more than the read. A public operation never has a precondition that user data can violate. The `Context` modifiers (`WithResolution`, ...) take values the caller chose rather than data, and guard them with contracts.
- An operation either returns a shape that passes the validity checker (ADR-0018) or returns an error. Returning an invalid shape is a bug; where a post-check catches it, the error is `InvalidResult` with the report attached.
- The C API maps `ErrorCode` to integer codes and exposes the message and entity list through accessor functions.

## Related result contracts

- [ADR-0028](0028-purpose-specific-numerical-tolerances.md): numerical diagnostics identify purpose and units; additional quantities are not unlabelled millimetre values.
- [ADR-0029](0029-intersection-and-classification-outcomes.md): public complete-result queries distinguish established empty/boundary outcomes from unresolved computation; partial evidence accompanies failure.
- [ADR-0030](0030-entity-lineage-in-operation-results.md): history references include owning input/output manifests and publish atomically with the shape.
- [ADR-0031](0031-attribute-propagation-and-conflicts.md): conflicts, losses and strict required-key checks complete before publication; geometry success alone does not override metadata-policy failure.

These proposals define future result requirements. Any additional enum codes are appended when implemented; current numeric error-code values are not changed here.

## Consequences

- Callers cannot ignore failures accidentally (`[[nodiscard]]` on all results).
- Diagnostics make corpus failures triageable automatically.
- Code becomes more verbose; monadic helpers (`and_then`, `transform`) and a `STEREON_TRY` macro keep it readable.
- Every shape-producing operation must assemble an `OperationResult`; a small amount of boilerplate per operation, which the builder's `Freeze()` and the history recorder supply automatically.
- Input validation at public entry points costs a `Structural` validity check per call (ADR-0018); accepted, because the alternative is terminating the host application on bad user data.

## Verification

Every corpus failure case must return a `KernelError` naming at least one involved entity. A CI check fails if any operation returns a shape that does not pass validation. A fuzzing run that feeds corrupt STEP files and invalid shapes to every public operation must produce only `KernelError`s, never a contract violation or a crash. Every tolerance change logged by ADR-0001 must also appear in the `OperationResult` of the operation that made it.

### Status of the verification

| Criterion | Status |
| --- | --- |
| A shape that fails validation is never published | Met: `Freeze` and `Commit` validate and return `InvalidResult` with the report's text (`src/topo/Builder.cpp`) |
| Tolerance changes appear in the `OperationResult` | Met: `Commit` reports every `SetTolerance` with its reason and warns `ToleranceGrown` (`tests/topo/OperationTests.cpp`) |
| Corpus failures name an entity | Pending: `KernelError` has no entity list or attached `ValidityReport` yet; the report's text is in the message. Both fields are added with the first operation that fails on a corpus case |
| Fuzzing yields only `KernelError`s | Pending: needs the fuzzing harness (TODO.md, section 9) |
