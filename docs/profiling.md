# Profiling Stereon

The initial implementation of [ADR-0043](adr/0043-optional-performance-instrumentation.md)
uses Tracy **v0.12.2**, pinned to commit
`c556831ddc6fe26d2fce01c14c97205a9dad46d5` (BSD-3-Clause).
It is optional and disabled by default. An ordinary build does not fetch or
link Tracy. No profiling configuration is required of kernel consumers.

## Build and capture

Use the reference GCC 16 toolchain:

```bash
cmake --preset gcc16-profile
cmake --build --preset gcc16-profile
ctest --preset gcc16-profile
```

This uses optimisation with debug symbols and the release contract policy.
The backend builds as ordinary C++17, without the kernel's reflection,
contract, or module-scanning flags. The private adapter uses Tracy's C API
with a constexpr RAII guard; runtime calls are bypassed in constant evaluation.

Build or obtain the Tracy GUI or `tracy-capture` tool from the same pinned
version. Collector tools are developer tools and are not built by Stereon.
For a command-line capture, start the collector in one terminal:

```bash
tracy-capture -a 127.0.0.1 -o primitives.tracy
```

Then run the repeatable workload in another:

```bash
build/gcc16-profile/tests/stereon_profiling_workload 100000
```

The iteration count controls the opportunity to attach; increase it if the
process exits before connection. On-demand capture records only while connected.
The workload creates boxes and cylinders, checks their published validity level
and face counts, and exercises rejection of an invalid radius. It is a capture
driver and consumer-link check, not a Boolean benchmark or speedup claim.

Only localhost connections are enabled; discovery broadcasts, system tracing,
and Tracy crash handlers are disabled. No collector is required to execute or exit. A profiling build may
start Tracy's diagnostic threads; a default build does not. Capture buffers can
consume substantial memory, so choose a bounded workload.

## What is instrumented

| Scope | Meaning |
| --- | --- |
| `build.MakeBoxSolid` | Complete framed-box construction, including Commit |
| `build.MakeCylinderSolid` | Complete cylinder construction, including Commit |
| `topo.Commit` | Validation, history, and result publication |
| `topo.Freeze` | Structural validation and snapshot creation |
| `topo.ValidateStructure` | Independent structural validation |
| `topo.ValidateGeometry` | Independent geometric validation |

Structural and geometric validators sample face and edge **slot counts**,
including dead slots. These plots are observations at each call, not cumulative
totals or counts of every entity visited. Concurrent calls share each named plot;
do not interpret a plot as per-operation attribution. Scopes are nested, so their
inclusive durations must not be added together.

Boolean intersection/classification scopes, operation correlation across tasks,
and optional allocation tracking remain future work. Required validation is
neither skipped nor weakened by enabling profiling.

## Adding instrumentation

Include the private `Profiling.hpp` in the global module fragment of an
implementation unit. Never include it in an exported module interface or public
header. Libraries receive it through a PRIVATE CMake dependency.

```cpp
STEREON_PROFILE_SCOPE("component.Operation");
STEREON_PROFILE_COUNTER("component.candidate_pairs", candidateCount);
```

Use static string literals and side-effect-free counter expressions representable
as `int64_t`. Disabled macros do not evaluate their arguments. Counter arguments
are also skipped in constant evaluation. Scopes close on early returns, must stay
on one thread, and must not cross coroutine suspension or task migration. Use
braces for nested scopes and one scope declaration per source line.

Do not put timing or counters into decisions about geometry or error handling.
Start with coarse scopes; avoid a trace event for every arithmetic evaluation.

## Verification and benchmark separation

The adapter test checks constant evaluation, disabled argument elision, early
returns, and worker-thread use. The primitive workload checks public consumer
linkage and rejects leaked profiling definitions. The existing primitive and
topology tests supply broader correctness coverage. Linux CI runs the full test
suite with `gcc16-profile` as well as the default configurations.

For performance evidence, use a separate uninstrumented build and record the
revision, input, compiler flags, hardware, tolerances, validation/contract policy,
thread budget, warm-up and repetition counts. Use matching optimisation settings
when measuring profiling overhead. Do not compare the `bench` preset's
contract-ignored results with contract-enforced profiling timings as an overhead
measurement. Report cold/preparation costs separately from reuse-heavy execution.

This first slice does not establish a performance advantage or complete all
ADR-0043 acceptance criteria. Collector capture inspection, paired corpus runs,
and representative overhead measurements are still required before acceptance.
