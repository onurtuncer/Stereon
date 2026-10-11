# ADR-0019: Library linking, packaging and ABI policy

- **Status:** Proposed
- **Date:** 2026-10-10
- **Owner:** Onur Tuncer
- **Phase:** 0

## Context

Stereon has separate CMake targets for its architectural layers, but has not decided whether those boundaries should also be shared-library boundaries. Today, `stereon_add_library()` calls `add_library()` without an explicit type, so `BUILD_SHARED_LIBS` controls the individual targets. This does not establish a supported binary distribution or ABI policy.

The kernel is pre-alpha, uses C++26 modules and evolving toolchain features (ADR-0016), and exposes value types and configuration-dependent declarations, including checked handles (ADR-0004). C++ applications built from source, Python extensions and applications hosting multiple plugins have different packaging needs. A shared library alone does not make these interfaces binary-compatible.

This ADR proposes the target policy; it does not claim the packaging, export annotations or verification gates below are implemented.

## Options

1. **Shared library per subsystem.** Maps directly onto today's CMake targets; multiplies exported interfaces, loader dependencies and configuration compatibility boundaries.
2. **Static libraries only.** Simple source integration and deployment; unsuitable as the sole supported model when several plugins must exchange live kernel objects through one kernel instance.
3. **Modular internal targets, static by default, with one optional shared kernel.** Keeps layering independent of packaging and provides a deliberate boundary for binary distribution.
4. **Header-only kernel.** Avoids a separately linked kernel binary but exposes implementation, increases consumer compilation work and does not remove toolchain compatibility requirements.

## Decision

Adopt option 3.

### 1. Internal targets and delivered library

- Preserve the architectural layers and dependency checks. Internal components may be static or object libraries; they are not independently versioned shared-library interfaces.
- Provide one aggregate consumer target, `stereon::stereon`, backed by either a static kernel or a shared kernel. The shared artifact is named `Stereon` using the platform's suffix and naming conventions.
- Introduce a project-scoped `STEREON_BUILD_SHARED` option, default `OFF`. Do not change the parent project's `BUILD_SHARED_LIBS` or let it silently switch Stereon's internal components into separate shared libraries.
- Static and shared packages are built in separate build trees. Their consumer target name remains the same; package metadata identifies the selected linkage.
- The aggregate must include the required component objects and propagate the correct usage requirements. Merely linking component archives into an otherwise empty archive or shared target is not an acceptable implementation.
- Optional external mesher adapters remain separate targets (ADR-0014). This decision does not require third-party dependencies or the operating-system runtime to be statically linked.

### 2. Exports and compiler runtime

- Hide implementation symbols by default and explicitly export the supported public interface. C++ module exports and dynamic-linker symbol exports are separate requirements.
- Use a common `STEREON_API` annotation where required, with platform-appropriate import/export definitions and a static-build form.
- Compile static/object components with position-independent code when they feed shared libraries or Python extensions on platforms requiring it.
- Choose compiler-runtime linkage independently of kernel linkage. Record and validate the compiler, standard library, runtime and ABI-affecting build configuration.
- For the C++ interface, consumers and kernel must use compatible toolchains and matching declaration/layout-affecting settings, including `STEREON_CHECKED_HANDLES` and `STEREON_HOT_CONTRACTS`. Export these requirements through the package; diagnose known incompatible combinations.
- Do not treat a Windows `.lib` suffix as proof of static linkage: it may be an import library for a DLL.

### 3. C++ ABI and modules

- No stable C++ ABI is promised during pre-alpha. Rebuild the kernel and C++ consumers together when toolchain, kernel version or ABI-affecting settings change.
- Module interface sources and their usage requirements belong to the development package. Compiler-generated module artifacts are toolchain/configuration-specific and are not a portable binary interface.
- Binary packages identify the kernel revision/version, target architecture, compiler and standard-library versions, runtime choice and relevant configuration.
- Before v1.0, review the ABI policy explicitly; reaching v1.0 does not automatically promise cross-compiler C++ ABI compatibility.
- Version shared artifacts so incompatible builds cannot silently replace each other. The exact platform versioning scheme is an implementation deliverable.

### 4. C boundary and object ownership

- Before offering a compiler-independent binary SDK, provide a versioned C ABI with opaque handles and explicitly sized scalar types. Do not expose STL types, templates, C++ object layouts or exceptions through it.
- Kernel-owned resources are released through matching kernel functions. Caller-provided buffers have explicit size and lifetime rules; callbacks have documented calling convention and lifetime.
- Translate failures to C status/diagnostic results (ADR-0006), including allocation failures at the C boundary. Programming-error contract violations retain ADR-0016 semantics and are not recoverable status results.
- Define an ABI version query and reject unsupported ABI versions before exchanging objects.
- This C ABI is a packaging milestone, not a prerequisite for initial source-based C++ development.

### 5. Python and plugins

- A self-contained Python extension may link the static kernel privately; hide its internal kernel symbols and package the required runtime dependencies.
- When a host, Python extension and/or multiple plugins exchange live kernel objects, they must use one compatible shared kernel instance. Do not pass objects, handles or deallocation responsibilities between independent statically embedded kernel copies.
- Independent instances may exchange serialized data through supported formats, with each instance creating its own objects.
- Runtime-loaded geometry plugins are not implied by the source-level type-erasure extension mechanism in ADR-0003. A dynamically loaded extension needs a separately specified ABI and unload/lifetime contract.

### 6. Performance and implementation sequence

- Link-time optimization is optional. Benchmark static and shared builds with comparable compiler settings; static linking is not a geometric-performance guarantee.
- First implement the explicit static default and aggregate consumer target. Then implement the shared artifact, export control and installed-package tests. Add the C ABI before distributing a compiler-independent SDK.
- Advertise shared packaging as supported only after its verification gates pass. Review and acceptance of this ADR remain human-owned (ADR-0015).

## Consequences

- Architectural boundaries remain explicit without creating a DLL for every layer.
- Source-based users get a simple default; plugin hosts have a supported route to one kernel instance.
- Shared packaging requires export annotations, ownership rules, relocation tests and ABI configuration checks.
- Static consumers must relink when adopting kernel changes; independent static copies can duplicate code and state.
- Linking mode does not replace dependency and licensing review under ADR-0008.

## Verification

Before accepting the implemented linking policy:

1. Build and test static and shared variants on the reference toolchain for Linux, Windows and macOS. Confirm that shared packaging produces one kernel shared artifact, not one per internal layer.
2. Install each variant and build a separate consumer using only `find_package(Stereon CONFIG REQUIRED)` and `stereon::stereon`. Verify module consumption and that no source/build-tree paths are required.
3. Relocate the installed package and run its consumer using the documented loader setup. Inspect exported symbols and runtime dependencies; fail for unintended internal exports or missing required public symbols.
4. Test a self-contained Python extension importing and executing an available kernel operation. Extend the test to shape creation and release once those APIs exist.
5. Test a host and two plugins using the shared kernel: create, pass, query and release a resource across module boundaries using its documented ownership API. Add C ABI version mismatch and allocation-failure tests when that ABI is introduced.
6. Reject known incompatible consumer configurations, particularly checked-handle layout mismatches. Verify that package metadata records the configuration needed to reproduce a build.
7. Run representative correctness cases in both linking variants and the ADR-0007 thread-count determinism checks within each build. Report static/shared performance separately; do not impose cross-build bit identity.

These are implementation acceptance gates, not claims that a documentation-only proposal has passed them.
