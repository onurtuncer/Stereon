# ADR-0028: Separate Numerical Tolerances by Purpose

- **Status:** Proposed
- **Date:** 2026-10-14
- **Owner:** Onur Tuncer
- **Phase:** 0 (policy); geometry in v0.1, intersections in v0.2, native shape persistence by v0.3
- **Extends:** ADR-0001, ADR-0002, ADR-0009, ADR-0012, ADR-0026

## Context

ADR-0001 separates modelling tolerance, numerical error bounds and healing budgets. It also defines linear and angular resolution. Algorithms still need distinct controls for relative comparisons, parameter-space steps, solver termination and approximation. Reusing model resolution for all these purposes mixes dimensions and allows a performance setting to change topology silently.

Stroud §16.1 distinguishes machine, length, angle, geometric, relative and polynomial tolerances and recommends preserving tolerances in saved models. The purpose types, assurance boundaries and persistence policy below are Stereon's interpretation of that guidance.

## Options

1. **One scalar epsilon for every comparison.** Easy to pass around; dimensionally unsafe and couples unrelated algorithms.
2. **Independent constants chosen inside each routine.** Locally convenient; impossible to audit, reproduce or coordinate through composed operations.
3. **Explicit purpose-specific policies with documented conversions and result bounds.** More interfaces and validation, but separates modelling semantics from numerical work controls.

## Decision

Adopt option 3. Preserve ADR-0001's model resolution, local tolerance, containment, growth ceiling and merge policies. This ADR adds purpose distinctions; it does not replace robust predicates with approximate equality or choose new global default numbers.

### 1. Every tolerance has a quantity, purpose and owner

| Quantity / purpose | Meaning and units | Owner / constraint |
| --- | --- | --- |
| Linear modelling resolution and entity tolerance | Spatial coincidence/containment allowance, mm | Model policy and entity records, ADR-0001; pair rules remain explicit and non-transitive. |
| Angular resolution | Angular comparison allowance, radians | Context policy; test normalized directions with valid norm bounds. A normal-angle test alone does not establish side transition (ADR-0029). |
| Relative comparison allowance | Dimensionless fractional allowance relative to a declared scale | Operation/query policy; name the scale and absolute floor in the compared quantity's units. Never silently increase entity tolerances. |
| Parametric allowance | Parameter interval or step allowance, in the curve/surface chart's coordinates | Evaluator/solver policy; parameter coordinates are not generally lengths and differ across charts. |
| Solver convergence targets | Residual, step and/or bracket criteria for the declared equation | Operation policy; each criterion declares its units or normalization and a bounded work budget. Convergence alone proves neither completeness nor geometric accuracy. |
| Approximation budget | Permitted deviation of a fitted or cached representation from its authority, usually mm, with separate angular/derivative budgets when requested | Operation policy; achieved whole-domain bounds/certificates follow ADR-0009/0010. Display tessellation has its own budgets. |
| Numerical error bound | Bound on arithmetic or construction error, with the computed quantity's units | Produced evidence, not a user-adjustable modelling tolerance. Machine epsilon is an arithmetic property, not a geometric threshold. |
| Healing budget / ceiling | Allowed explicit geometric/tolerance changes | Context and operation contract, ADR-0001/0026; never a solver stopping criterion. |

Use distinct named types or policy fields and purpose-specific accessors. APIs shall not accept an unlabeled `epsilon` for several rows. No implicit interchange of linear, angular, parametric and dimensionless controls is permitted. Public/imported settings validate finite values, permitted ranges and units before entering low-level programmer-contract APIs.

### 2. Conversions and error budgets are explicit

- Comparing squared distance with a linear allowance requires the corresponding squared limit and numerically safe arithmetic. Do not compare mm² with mm.
- A relative comparison can use `abs(a-b) <= absTol + relTol * scale`, only with a declared nonnegative scale in the same units as `a` and `b`. The absolute floor governs near-zero cases. Robust sign predicates retain their own exact/filtered semantics under ADR-0002; this formula is not a replacement for them.
- Mapping a parameter-space error into a spatial allowance requires an appropriate derivative/metric bound over the affected domain. Singular charts, poles and poorly conditioned maps require refinement, another chart or an explicit unresolved outcome; they do not justify assuming a universal UV-to-mm conversion.
- Distinguish a requested budget from an achieved error bound. Residual norms with different dimensions must be normalized or tested separately. A small iteration step can mean stagnation; it is not evidence of a small geometric residual.
- Composite operations must account for propagated errors using justified bounds. Shared or correlated contributions must not be hidden by treating each local budget as the entire output allowance. The final required certificate must cover the whole supported domain.
- Loosening a display or approximation setting must not silently redefine B-rep coincidence, remove features or enlarge a modelling tolerance. A construction exceeding its permitted output tolerance requires the explicit ADR-0001 growth path or failure.

### 3. Settings and changes are observable

Snapshot the effective numerical policy at operation entry, including caller overrides and algorithm policy/version identifiers needed for replay. Context settings remain immutable during a call; workers do not consult mutable global defaults.

Adaptive subdivision or iteration is permitted inside the declared policy. Changing solver technique while meeting the same promised bounds is not healing. Relaxing a requested acceptance target requires explicit caller authorization through the operation policy, records requested/effective values and the reason, and reports a diagnostic. Relaxing work criteria must never bypass a required output certificate. Actual entity tolerance growth remains in `OperationResult` even when logging is disabled, and obeys warning promotion and the ceiling.

Changing the Context for a later operation does not retroactively change stored entity tolerances or validate an old certificate under a different model policy. Compatibility between input model policies and the requested operation policy must be checked; any re-resolution/healing is explicit.

### 4. Persist model meaning and replay context

Under ADR-0012, native models store the units, model resolution, angular policy and other numerical policy values required to interpret their geometry, tolerances and certificates. Version this metadata with the schema; a reader must not silently substitute its current defaults for missing values. Migrations may supply a historical default only when the old schema defined it unambiguously; otherwise report the limitation or refuse interpretation that cannot be justified.

Solver work limits and algorithm-specific settings belong to the creating operation's journal/provenance when replay is promised. They are not all geometric invariants of the final shape. Persist effective overrides and policy identifiers in a serializable snapshot, not function pointers, sinks or a raw `Context` memory image. A journal needs enough input and build information for ADR-0007's replay envelope; storing a tolerance alone is insufficient.

Required model-semantic fields follow ADR-0012's required-chunk rules. If optional replay metadata is unavailable, report that exact replay is unavailable without claiming the model geometry is invalid.

## Consequences

Purpose distinctions make numerical changes reviewable and model files interpretable. They add validation and schema work. Existing `Context` resolution/angular/ceiling settings are a starting subset; the current error diagnostic's unlabelled millimetre values must gain purpose/unit information before it reports the additional quantities. This ADR does not claim those interfaces are implemented or prescribe a stable C++ layout.

## Verification

Implement these checks with the corresponding policy, algorithm and persistence features:

1. Unit-equivalent inputs within the supported size/resolution domain give equivalent outcomes after normalization; every length-valued allowance scales correctly and angular/dimensionless allowances do not.
2. Changing angular, convergence or display controls cannot silently modify entity tolerances. Every authorized relaxation/growth appears in structured results and replay context, including with no log sink.
3. Near-zero relative comparisons, squared-distance limits, rescaled parameterizations and a singular chart exercise the appropriate quantity and fail explicitly when a spatial bound cannot be established.
4. A stagnated solver with a small step and large residual does not report a certified intersection. Exhausted work and failed certification produce their documented errors.
5. Native round trips preserve model numerical semantics; old-schema fixtures use explicit migrations; unknown required policy data is rejected. Replay tests restore the recorded effective policy rather than current defaults.
6. A composed construction either establishes its final bound within the budget or takes an explicit authorized healing/failure path. G1 outcomes are identical across thread counts.

## References and scope

- Ian Stroud, *Boundary Representation Modelling Techniques*, Springer, 2006, §16.1, printed pp. 603–604 (reviewed in the supplied book).
- [ADR-0001](0001-tolerance-model.md), [ADR-0002](0002-numeric-robustness.md), [ADR-0009](0009-parameter-space-representation.md), [ADR-0012](0012-native-file-format.md), [ADR-0026](0026-modelling-operations-by-contracts.md), [ADR-0029](0029-intersection-and-classification-outcomes.md).

This proposal supplements the listed decisions, supersedes none, and changes no ADR acceptance status. Verification items are future acceptance requirements.
