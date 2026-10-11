# ADR-0020: C++ API, C ABI and Python bindings

- **Status:** Proposed
- **Date:** 2026-10-11
- **Owner:** Onur Tuncer
- **Phase:** 0 (policy); implementation staged below
- **Related:** ADR-0003, ADR-0006, ADR-0007, ADR-0011, ADR-0016, ADR-0017, ADR-0018; companion ADR-0019 draft (library linking, packaging and ABI policy)

## Context

Stereon needs a full native interface for C++ applications, convenient scripting for geometry exploration and engineering workflows, and a binary integration boundary usable from other languages.

These are different requirements. C++ concepts, value types and type erasure are useful native abstractions, but their layouts and calling conventions are not a cross-compiler ABI. A C-compatible interface offers a narrower interoperability boundary, but still needs explicit versioning, ownership and platform conventions. Python provides the interactive environment needed for notebooks, regression cases, optimisation and CFD/FEA preparation.

ADR-0006 now defines shape-producing operations in terms of OperationResult: Shape, History, tolerance changes and structured warnings. Bindings must preserve this information rather than reduce every successful operation to a shape.

## Options

1. **C++ and direct Python bindings only.** Fast to develop; leaves other languages without a deliberately supported binary boundary.
2. **C ABI first; all language bindings use it.** Uniform boundary, but every early API experiment requires a C representation and a binding update.
3. **C++ native API, early direct Python bindings, staged C ABI for established operations.** Supports interactive development while allowing the binary interface to mature deliberately.
4. **C ABI only as the public interface.** Small binary boundary, but loses the native extension model and C++ ergonomics.

## Decision

Adopt option 3. Geometry and topology algorithms remain implemented once in the C++ kernel. Both the C adapter and Python bindings delegate to that implementation.

### 1. Native C++ interface

- The C++ API remains the most complete interface, including source-level custom geometry through ADR-0003.
- Preserve value semantics, immutable shapes, typed handles, structured results and explicit Context ownership.
- C++ consumers use compatible toolchains and configuration settings. No cross-compiler C++ ABI is promised.
- Experimental APIs may be exposed to Python without becoming part of the supported C ABI.

### 2. Early Python bindings

- Use nanobind to bind the C++ API directly. Python does not have to call through C.
- Build and distribute each extension with its matching kernel and runtime dependencies, following the companion linking policy.
- Python functions use snake_case and classes use PascalCase. Apply explicit, reviewed naming mappings rather than automatic renaming of every identifier.
- Shape-producing functions return a Python OperationResult exposing shape, history, tolerance changes and warnings. Do not silently discard diagnostics.
- Map KernelError to a Python exception carrying its code, message and structured diagnostics, including a validity report and replay identifier where available. Allocation failure becomes MemoryError.
- Expose batch operations and contiguous numerical buffers for evaluation and mesh data; avoid requiring a Python call per vertex or triangle.
- Default exported arrays to owned copies. A zero-copy view requires a pinned owner, documented lifetime and read-only access to immutable kernel storage.

### 3. C ABI scope

Introduce a small C-compatible interface for established operations before distributing a compiler-independent binary SDK.

Initial coverage should include:

- context creation and release;
- shape creation, import/export and Boolean operations as they become available;
- validity and mass-property queries;
- meshing with bulk extraction of coordinates, connectivity and tags;
- access to operation history, tolerance changes, warnings and errors.

Use opaque resource types and stn_ function names. Keep C++ classes, STL containers, templates, exceptions and internal topology handles out of public C declarations.

The first C ABI need not expose source-level geometry extension, schedulers or every low-level evaluator. Unsupported capabilities remain explicit rather than being approximated by a second implementation.

### 4. C representation, ownership and errors

- Public headers compile as C and as C++, with extern "C" guards and an explicit platform calling/export convention.
- Use explicitly sized integers for ABI status codes, flags and counts. Specify numeric representations, string encoding (UTF-8), buffer layout and alignment. Do not use C/C++ enum layout or bool layout as an unstated ABI assumption.
- Every allocation has one documented owner and a matching release function. Kernel-owned memory is released by the kernel; caller-owned buffers are never freed by it.
- Return operation-owned diagnostics, not a global or thread-local "last error". Reading diagnostics must not overwrite those of another operation.
- An operation returns a status plus an opaque operation-result resource on success and an opaque error resource on failure. Set unused output resource pointers to null. Specify borrowed versus retained access to shapes and history.
- Provide allocation-free status reporting for out-of-memory failures, where allocating an error object may itself fail.
- No C++ exception crosses the C boundary. Catch allocation failures and defensively contain other escaping C++ exceptions with an internal-error status; do not turn termination or contract violations into recoverable errors.
- Validate documented external data and buffer sizes before entering kernel algorithms. Arbitrary invalid pointers or use-after-release are caller errors; opaque pointer APIs cannot safely validate every such misuse.

### 5. Versioning and packaging

- Publish ABI version and capability queries independently of the kernel's product version.
- Once an ABI version is declared stable, keep its existing function signatures and status values compatible. Additive functions are capability/version gated; incompatible changes require a new ABI major version.
- Specify the supported operating system, architecture, calling convention and data layout. "C ABI" does not imply one binary works on every platform.
- Hosts and plugins exchanging live kernel resources use the same compatible kernel instance. Independent instances exchange serialized data, not live pointers or handles.
- Python may retain direct C++ bindings after the C ABI ships. Reconsider routing Python through C only if binary distribution requirements justify the additional adapter work.

### 6. Concurrency and callbacks

- Binding calls follow the kernel's documented Context and shape-sharing rules. No new global mutable state is introduced.
- Release Python's GIL around long-running native work only after inputs and lifetimes are secured, and only when that work does not access Python objects. Reacquire it before constructing Python results or invoking Python code.
- Initially exclude Python-defined geometry and callbacks from GIL-released worker execution. Their thread, exception and lifetime semantics require a separate design.
- The initial C ABI avoids callbacks. Any later callback API must define thread affinity, reentrancy, lifetime, cancellation and failure propagation.
- Retain input resources for the duration of a call. A binding must not let concurrent release invalidate in-flight work.

### 7. Implementation order

1. Implement the C++ operation and its correctness tests.
2. Expose useful operations through direct Python bindings for scripting and engineering validation.
3. Prototype the C boundary early with resource creation/release, one query and a structured failure.
4. Add established operations to the C ABI with explicit ownership and diagnostics.
5. Complete the binary SDK compatibility and packaging gates before advertising a stable C ABI.

This ADR defines a proposal and acceptance criteria; it does not claim that these interfaces or guarantees already exist. Its numbering follows the companion ADR-0019 draft, which must be reconciled with the repository index when submitted.

## Consequences

- C++ users retain the full native API; Python users get early access; other languages have a deliberate binary integration route.
- Algorithms are shared, but two adapters and their semantic parity tests must be maintained.
- The C ABI can stay smaller and more stable than the evolving C++ interface.
- Direct Python bindings remain tied to a compatible kernel build.
- Resource lifetimes, diagnostic fidelity and bulk data transfer are public API design responsibilities, not binding implementation details.

## Verification

Before declaring each interface supported:

1. Build a separate C translation unit against the installed C header and link/run it against the packaged kernel. For a compiler-independent SDK claim, also test a documented consumer compiler different from the kernel compiler.
2. Run the same operation fixtures through C++, Python and C where exposed. Compare validity, geometry within the operation's documented bounds, history, tolerance changes and diagnostic codes.
3. Test success and failure ownership paths, borrowed/retained lifetimes, null-output initialization and leak-free cleanup under applicable sanitizers.
4. Inject allocation failures: the C boundary returns an allocation-free failure status and Python raises MemoryError. Neither boundary returns an invalid success result.
5. Test buffer-size handling, count overflow, UTF-8 strings and mesh array layouts. Test that zero-copy views, if provided, retain their owners and cannot mutate immutable shapes.
6. Verify that long Python operations permit another Python thread to progress when the GIL is released, without permitting unsafe Context sharing or callback execution.
7. Once an ABI version is stable, run an old compiled C consumer against compatible newer kernel releases and test rejection of unsupported ABI versions.
8. Exercise both a self-contained Python package and the supported shared-kernel configuration without relying on the source tree.

Documentation review can approve the direction; stable SDK support requires these implementation gates to pass.
