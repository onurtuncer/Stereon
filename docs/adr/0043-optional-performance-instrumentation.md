# ADR-0043: Optional performance instrumentation and reproducible benchmarks

- **Status:** Proposed
- **Date:** 2026-10-20
- **Owner:** Onur Tuncer
- **Phase:** 0
- **Related:** ADR-0006 (errors), ADR-0007 (threading and determinism), ADR-0016 (C++26), ADR-0018 (validation), ADR-0019 (packaging), ADR-0023 (logging), ADR-0027 (topology mutation)

## Context

Boolean performance depends on candidate generation, numerical intersection,
classification, topology construction, and validation. Total elapsed time alone
does not explain which stage or input characteristic causes a regression.
Instrumentation must help diagnose these costs without becoming a public
dependency, changing modelling semantics, or preventing constant evaluation.

Profiling and benchmarking serve different purposes: a trace explains execution;
a reproducible benchmark establishes performance. Instrumentation can perturb
scheduling, allocation, and runtime, so trace timings alone cannot substantiate
speedup claims.

## Options

1. **Build a profiler.** Complete control, but substantial unrelated work on
   capture, storage, visualisation, and platform support.
2. **Use an existing profiler through a private adapter.** Reuses mature tools
   while keeping the kernel independent of their APIs.
3. **Use only external sampling and total operation timers.** Minimal integration,
   but limited visibility into algorithm stages and workload counts.

## Decision

Adopt option 2, with Tracy as the initial tracing backend. External sampling
profilers remain useful complementary tools. Do not build a profiler.

### Dependency and build boundary

- Provide an explicit `STEREON_ENABLE_PROFILING` build option, off by default.
  Disabled builds neither fetch nor link Tracy and do not evaluate arguments
  used solely for profiling.
- Keep named scopes and counters behind a small private adapter. No Tracy types,
  macros, headers, or configuration enter the exported API, C ABI, or consumer
  interface requirements. Any header bridge remains in implementation details.
- Pin the backend version and record its licence and build configuration.
  An enabled static build may need backend linkage; packaging must carry it
  privately without requiring applications to call Tracy.
- Backend replacement may use an ordinary reviewed change if these boundaries
  and verification criteria remain satisfied.

### Semantics and lifetime

- Profiling is observational. Timestamps, counters, connection state, and trace
  delivery never influence geometry, tolerances, classification, entity ordering,
  validation policy, errors, or transaction publication.
- Preserve constant evaluation. Runtime instrumentation is bypassed during
  constant evaluation; profiling-enabled builds must retain existing constexpr
  use cases. Verify the adapter with the supported compiler and module setup.
- Keep modelling caches and operation state in their existing owners.
  As a narrow qualification to the proposed ADR-0007 ban on global mutable
  state, an explicitly enabled diagnostic backend may own process-wide capture
  state and capture threads. Such state must not affect modelling results.
  No capture service starts in the default build.
- Scope lifetime must cover actual work, including asynchronous completion.
  Do not keep a thread-bound scope open across task migration. Use supported
  task correlation or worker scopes for parallel work.
- Missing collectors or dropped events do not fail a modelling operation.
  Capture completeness is reported separately. Runtime and memory overhead
  remain real costs and must be measured.

### Instrumentation scope

Start with coarse stages; add finer scopes only to investigate observed costs.

| Stage | Useful counters |
| --- | --- |
| Candidate generation | Entities considered, candidate face pairs |
| Intersection | Actual intersections, solver iterations, fallback calls |
| Trimming and classification | Fragments processed, unresolved outcomes |
| Topology and metadata | Entities created, modified, or reused |
| Validation and publication | Entities visited, validation failures |

Use bounded labels and operation correlation identifiers. Avoid per-evaluation
dynamic strings and shared atomic counters in arithmetic loops; aggregate
task-local counts where appropriate. Allocation tracking is separately opt-in.
Profiling traces do not replace logs, replay journals, or independent validators.
Retain mandatory checks and report their cost as part of the operation.

### Benchmark evidence

- Establish performance claims using optimised builds with profiling disabled.
  Keep diagnostic captures separate, using matching optimisation settings.
- Record source revision, compiler and flags, hardware, thread budget, input
  hashes or reproducible generators and seeds, tolerances, contract mode,
  validation policy, and requested history or metadata behaviour.
- Measure end-to-end latency through completed validation and publication.
  Declare whether import and preparation are included; report excluded costs.
- Report cold and reuse-heavy runs separately. Include preparation costs or
  explicitly state the amortisation schedule.
- Use repeated runs with a declared warm-up policy and repetition count.
  Report median, tail latency with adequate samples, variability, and memory.
- Compare against a pinned, appropriately configured baseline with equivalent
  accuracy, output requirements, and resource budgets. Report correctness
  failures, unresolved results, and timeouts; never silently drop them.
- Keep performance gates on controlled runners. General shared CI checks
  benchmark execution and correctness rather than enforcing noisy timing limits.
  This ADR establishes a measurement policy, not a promised speedup.

## Consequences

The kernel gains stage-level explanations without maintaining a profiler.
The adapter, benchmark corpus, and result metadata require maintenance.
Profiling builds have overhead and may perturb scheduling; their traces are
diagnostic evidence rather than the source of headline benchmark numbers.
Default consumers acquire no profiling dependency.

## Verification

Before acceptance, demonstrate the following on the supported reference toolchain:

1. Disabled builds work without Tracy available. Inspect a representative
   optimised object and link dependencies: no backend calls or symbols and no
   evaluation of profiling-only argument expressions remain.
2. Enabled builds capture named stages and counters for an available modelling
   operation. Add the Boolean stages as their implementations become available.
3. Existing constexpr checks compile with profiling both enabled and disabled.
4. Paired corpus runs retain equivalent geometry, topology, history, and error
   outcomes under identical modelling settings. Scheduler/thread-count runs
   continue to satisfy the applicable determinism policy.
5. Running without a collector completes normally; early-return and failed
   operation paths leave balanced scopes. Parallel captures correctly attribute
   worker execution rather than only submission time.
6. A reproducible benchmark report includes the metadata, correctness outcomes,
   cold/warm distinction, and timing statistics specified above.
7. Measure enabled instrumentation overhead against the disabled build on
   representative small and large cases; publish the result rather than assuming
   negligible cost.

Reference: [Tracy project and documentation](https://github.com/wolfpld/tracy).
