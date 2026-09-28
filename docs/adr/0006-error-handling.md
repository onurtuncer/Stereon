# ADR-0006: Error handling

- **Status:** Proposed
- **Date:** 2026-09-28
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
- `KernelError` holds: an `ErrorCode` enum (e.g. `InvalidInput`, `ToleranceExceeded`, `IntersectionFailed`, `DegenerateResult`, `NotImplemented`, `ParseError`); a human-readable message; the entities involved (handles into the input shapes); the tolerance values at stake; and an optional journal ID for replay.
- The kernel does not throw. `std::bad_alloc` may still propagate. Programming errors (violated preconditions) are contract violations (ADR-0016), not `KernelError`s.
- An operation either returns a shape that passes the validity checker or returns an error. Returning an invalid shape is a bug.
- The C API maps `ErrorCode` to integer codes and exposes the message and entity list through accessor functions.

## Consequences

- Callers cannot ignore failures accidentally (`[[nodiscard]]` on all results).
- Diagnostics make corpus failures triageable automatically.
- Code becomes more verbose; monadic helpers (`and_then`, `transform`) and a `STEREON_TRY` macro keep it readable.

## Verification

Every corpus failure case must return a `KernelError` naming at least one involved entity. A CI check fails if any operation returns a shape that does not pass validation.
