# ADR-0021: Shape derivatives as an operation contract

- **Status:** Proposed
- **Date:** 2026-10-11
- **Owner:** Onur Tuncer
- **Phase:** 0 (contract); implementation staged below
- **Related:** ADR-0001, ADR-0002, ADR-0003, ADR-0004, ADR-0005, ADR-0006, ADR-0007, ADR-0008, ADR-0009, ADR-0010, ADR-0011, ADR-0013, ADR-0014, ADR-0019, ADR-0020

## Context

Stereon should support gradient-based engineering design, including CAD-driven CFD and FEM optimisation. Retrofitting this after geometric operations, topology management and public interfaces have stabilised would require invasive changes to their inputs and outputs.

Automatic differentiation (AD) of arithmetic alone is not enough. A useful derivative also needs a well-defined parameterisation, a correspondence between the entities evaluated at two parameter values, and a construction that is locally smooth. Boolean classification (ADR-0013), healing and the merge policy (ADR-0001), adaptive sampling and solver branch selection can all change the function being differentiated.

The existing decisions already provide most of the foundation: immutable shapes (ADR-0005), explicit tolerances (ADR-0001), structured results and errors (ADR-0006), operation history (ADR-0011), capability tiers for geometry (ADR-0003) and a deterministic `Context` (ADR-0007). This ADR extends them with sensitivity contracts. It does not claim that any derivative implementation exists yet.

The keywords MUST, SHOULD and MAY state the strength of each rule below.

## Options

1. **Add finite differences later.** Useful for checks and explicit fallbacks; does not establish dependency or correspondence semantics, and a perturbed evaluation can silently cross a branch event.
2. **Replace every scalar with one AD type.** Couples the whole kernel to one backend, slows the primal path, and does not resolve discrete topology or singular constructions.
3. **Differentiate only the final spline control points.** Useful as one parameterisation; insufficient for feature parameters, moving trims and retained design constraints.
4. **Differentiate every solver iteration by default.** Makes derivatives depend on stopping and branching behaviour; differentiating the converged equations expresses the intended construction better when it is regular.
5. **Guarantee differentiability through all topology changes.** Not a valid general contract for ordinary derivatives of B-rep operations; it needs nonsmooth or topology-optimisation methods.
6. **Differentiability as a declared per-operation contract.** Smooth numerical geometry is separated from discrete topology and classification; each operation declares what it can differentiate and under which assumptions, and reports when they fail.

## Decision

Adopt option 6. First-order derivatives are supported within declared regular branches and topology regimes, with explicit statuses and diagnostics when those assumptions do not hold. Options 1 and 3 remain available as a verification tool and as one parameterisation respectively.

Every new geometric operation MUST declare its parameter dependencies, correspondence rules, differentiability assumptions and derivative capabilities. An initial implementation MAY declare a product `NotSupported`; it MUST NOT silently return zero, detach a dependency, or imply that a derivative is valid outside its declared regime.

### 1. Derivative semantics

For a surface `x = S(u, v; p)` the API distinguishes:

| Quantity | Meaning |
| --- | --- |
| `S_u`, `S_v` | Geometric derivatives with respect to surface coordinates; what `Deriv` already returns (ADR-0003). |
| `S_p` | Design-parameter derivatives at fixed surface coordinates. |
| `dx/dp` | Total derivative of a tracked point under its declared correspondence rule. |
| `DJ(Ω)[V]` | Derivative of a shape functional under a boundary velocity field; generally assembled by an analysis integration, not by the kernel. |

For a moving surface coordinate the total derivative is

```text
dx/dp = S_p + S_u · du/dp + S_v · dv/dp.
```

The initial kernel contract covers parameter sensitivities, tracked boundary velocities and selected geometric functionals (area, volume). It does not claim that a CAD parameter gradient is automatically a parameterisation-independent shape gradient.

### 2. Scalar-generic numerical geometry

- Differentiable arithmetic is written so it runs on `double` and on supported active scalar types. Initial targets: vectors, transforms, analytic primitives, spline evaluation, and the normal and integration calculations that functionals need.
- Differentiable paths MUST NOT convert active values to `double` implicitly. Constants and mathematical functions preserve the scalar type. Any deliberate extraction of a primal value for a discrete decision is explicit and covered by the operation's branch contract (§5).
- **Fit with ADR-0003.** Scalar genericity applies to the built-in geometry aggregates and their free functions, reached through `VisitBuiltin`. The type-erased `Curve` and `Surface` wrappers and their virtual fallback stay `double`. Parameter derivatives of user geometry are a capability like the ADR-0003 tiers: a user type that provides parameter-derivative operations (an optional concept, named when the first slice lands) gains derivative support; one that does not reports `NotSupported` for those products. The kernel never differentiates user geometry by guessing.
- AD-library-specific behaviour is isolated behind internal adapters. This ADR selects no AD backend, does not make `Shape` scalar-templated, and does not require topology containers to hold active scalars. Any backend that is adopted must meet ADR-0008's licence rules and stay a private dependency (ADR-0019).
- Topology identifiers, connectivity, classifications and robust predicates stay discrete and `double`-based (ADR-0002). Robustness MUST NOT be weakened to make classification differentiable.

### 3. Evaluable dependencies and immutable results

- Parameter identities, units, scaling conventions and feature dependencies are explicit. Lengths follow ADR-0001 (millimetres); angles are radians.
- A sensitivity-capable construction retains, or gives access to, enough information to re-evaluate its parameter-dependent geometry. History (ADR-0011) alone is insufficient: knowing that one face came from another does not define its derivative. The dependency representation includes the operation's inputs and its evaluation rules for supported paths.
- The dependency graph is independent of any AD tape. Tapes, seeds, adjoints and solver workspaces live in a **sensitivity context** owned by the caller for one evaluation, never in the immutable `Shape` (ADR-0005). It is separate from the operation `Context` of ADR-0007, which keeps its role (tolerances, schedulers, caches); two sensitivity contexts share no mutable state.
- Caches account for the parameter state, the evaluation mode and the branch and tolerance policy. A primal-only cached value MUST NOT silently sever a derivative.
- Imported geometry with no feature history MAY expose control points or other explicit geometric parameters. Stereon does not invent design intent that the data does not carry.

### 4. Correspondence and trims

- A tracked sample identifies its source entity, its orientation where relevant, its local coordinates and its tracking convention. Supported conventions include fixed surface coordinates, fixed edge parameter, moving trim coordinates and a uniquely defined projection.
- Correspondence is defined within the declared regime by the operation's construction and its branch signature (§5), not by persistent naming, which ADR-0011 leaves out of scope for v1.0. Ambiguous or lost correspondence is reported. Entity indices are not correspondence: `Compact()` renumbering (ADR-0004) or any storage reordering MUST NOT redefine it.
- Trimmed geometry needs the trim curves, their parameter-space mappings and intersection endpoints. For ADR-0009 edges this means the pcurves on their shared parameter: an `Intersection` edge moves as the solution of `S₁(c₁(t)) = S₂(c₂(t))`, a `Trim` or `Boundary` edge with its pcurve, a `Spatial` edge with its independent 3D curve. Evaluating the untrimmed surface alone is never advertised as a total derivative of a trimmed boundary.
- Derivatives describe smooth geometry, not certificates. Moving a parameter does not carry an ADR-0009 certificate along; a re-evaluated primal is re-certified as usual.
- Spline degree and knot multiplicity changes are discrete structural changes. Variable knots, weights and other advanced parameters need their own regularity and admissibility contracts before they are supported.

### 5. Capability and validity are reported separately

An operation declares which derivative products it implements; the declaration is queryable before evaluation. Each evaluation also reports the validity of the requested derivative:

| Status | Meaning |
| --- | --- |
| `Valid` | Computed under the declared regularity assumptions and the applicable checks. Not a global mathematical certificate. |
| `NearEvent` | A local derivative is available, but diagnostics indicate proximity to a branch event or poor conditioning. |
| `Nondifferentiable` | The requested ordinary derivative is not defined under the selected semantics. |
| `NotSupported` | The requested product or operation has no implementation. |
| `Indeterminate` | Correspondence, convergence or regularity could not be established well enough to return a trustworthy derivative. |

- **Relation to ADR-0006.** A sensitivity evaluation returns `std::expected<SensitivityResult, KernelError>` (name provisional). Ordinary geometric and numerical failures of the evaluation, and invalid seeds or buffer sizes supplied by the caller, are `KernelError`s (`InvalidInput` for the latter). The statuses above describe a completed evaluation. `NotSupported` reuses the vocabulary of `ErrorCode::NotSupported`; it is a status here because an optimiser must be able to ask for many products and continue with those that exist.
- `NotSupported`, `Indeterminate` and `Nondifferentiable` results carry no usable derivative payload. A `NearEvent` payload is explicitly qualified for the caller.
- Results SHOULD include a branch/topology signature, the affected entities, the tolerance policy, solver residuals and conditioning or event diagnostics where available. Event distances or margins state their units and whether they are estimates or certified bounds. Absence of a detected event is not proof that a finite parameter step is safe.
- Fixed topology alone does not establish differentiability. Tangential intersections (the `Tangential` and `Overlap` contact types of ADR-0010), ambiguous closest points, coordinate singularities (`Degenerate` edges at poles and apexes, ADR-0009) and rank-deficient constraints need their own checks. The contract states branch and active-set assumptions where they apply.

### 6. Implicit constructions are differentiated through their equations

- For a supported construction defined by `F(z, p) = 0`, the directional derivative SHOULD be obtained from

  ```text
  F_z · ż = −F_p · ṗ.
  ```

- Implementations SHOULD expose the residual and a suitable linearisation, reuse factorisations where appropriate and provide transpose solves for reverse products. Explicit matrix inversion is not required.
- This applies to projections, intersections (ADR-0010) and constrained constructions. It needs a well-posed local system, adequate primal convergence and conditioning checks. Under-determined constructions define their parameterisation or supplementary constraints; a surface–surface intersection curve needs a tracking convention along the curve, which the shared edge parameter of ADR-0009 supplies.
- Differentiating a solver's iterations MAY be used when the derivative of that implemented algorithm is explicitly intended and justified. It is not the default substitute for differentiating a converged geometric construction.

### 7. JVP and VJP interfaces

For a sampled geometry map `x(p)` with Jacobian `J`, the sensitivity API accommodates both

```text
jvp(ṗ) = J · ṗ
vjp(x̄) = Jᵀ · x̄
```

- Dense Jacobian assembly MAY be provided for small problems but is never the only interface. Forward and reverse coverage MAY be staged and are discoverable independently.
- The API defines sample and parameter ordering, units, coordinate frames, seed conventions, and whether adjoints accumulate into or overwrite outputs. VJP inputs are covectors under a declared pairing. A continuous normal sensitivity density needs normals and integration weights before it can be treated as a discrete coordinate covector.
- **Determinism.** Reverse products accumulate in the fixed order of ADR-0007, keyed by sample and entity index, so a VJP is bit-identical for any thread count (G1).
- C++ is the primary interface. Python bindings (ADR-0020) expose the products on contiguous numeric buffers without exposing AD scalar types or tape internals. A later C ABI represents sensitivity contexts as opaque `stn_` resources with numeric buffers. Framework-specific autodiff adapters (PyTorch, JAX) stay outside the core contract.
- C++ names follow ADR-0017 when the first slice lands; the formulas above are notation, not API names.

### 8. Discrete operations and tolerance policy stay visible

- Tolerances (ADR-0001) are held fixed during a derivative evaluation unless a documented operation explicitly takes them as parameters. Healing, merging, tolerance growth and branch changes appear in the `OperationResult` (ADR-0006) and in the sensitivity diagnostics; a tolerance change between two evaluations is a branch event.
- Initial derivative validation uses a declared correspondence and a fixed sampling and connectivity regime. Adaptive remeshing, retessellation and topology-changing Boolean outcomes are outside that initial contract. The analysis-grade surface mesh of ADR-0014 is a natural source of tracked samples when each node records its face and (u, v).
- A supported Boolean derivative is the derivative of one regular local construction, not a guarantee that Booleans are globally differentiable. Crossing an event requires rebuilding the primal result and re-establishing correspondence and regularity before a new local derivative is computed.

### 9. Implementation order

1. **Foundation:** parameter identities, capability declarations, correspondence types, statuses and sensitivity-context ownership; scalar-generic numerical primitives without a permanent public AD backend.
2. **First vertical slice:** a parameterised primitive and a spline example, tracked boundary samples and JVPs, demonstrating the full parameter-to-sample dependency path.
3. **Optimisation interface:** VJPs, a documented discrete adjoint pairing, area and volume derivatives for supported regular shapes, and Python exposure.
4. **Implicit geometry:** projections and intersections with conditioning diagnostics and explicit tracking constraints, then moving trims and selected Boolean regimes.

Each stage declares its supported operations and products. Existing operations are classified `NotSupported` while they are brought into compliance; no claim of complete kernel differentiability is made.

## Consequences

- Geometry operations can take part in CAD-to-analysis optimisation without a later redesign of their inputs and outputs.
- Robust predicates and topology management keep their existing responsibilities (ADR-0002, ADR-0004).
- Directional products support large parameter sets without dense Jacobians.
- Unsupported and invalid derivative cases stay explicit to optimisers and analysis integrations.
- Immutable shapes stay reusable without sharing mutable differentiation state.
- Operation authors must define parameter dependencies, correspondence and regularity assumptions before derivative coverage is complete; review checks the declaration as it checks history recording (ADR-0011).
- Scalar-generic code, sensitivity contexts and implicit linearisations add implementation and maintenance work, and built-in geometry evaluators must compile for more than one scalar type. The ADR-0003 dispatch benchmark must show the `double` path unaffected.
- Forward and reverse coverage may differ during rollout.
- Event diagnostics can be incomplete or conservative; no certified topology-preserving step bounds are promised.
- Differentiating the CAD boundary does not differentiate volume meshing or an external CFD/FEM solver.

## Verification

The foundation is complete when a new geometric operation cannot be added without declaring its derivative capability, dependencies, correspondence and validity assumptions (enforced by a registration check in tests). The first derivative path is accepted when:

- A parameterised model produces an immutable primal shape and a separate sensitivity context.
- Analytic primitive checks agree with known coordinate and geometric-functional derivatives.
- For smooth test cases, the directional Taylor residual `‖x(p + h·v) − x(p) − h·J·v‖` shows second-order decay over a stated step range before round-off dominates.
- Implemented forward and reverse products satisfy `⟨J·v, w⟩ = ⟨v, Jᵀ·w⟩` to a stated tolerance.
- At least one moving-coordinate example exercises the full point derivative, including the coordinate-motion terms.
- Tangency, collapsing geometry, lost correspondence and a tolerance-triggered branch change each produce the documented status or `KernelError`, never a silently usable derivative.
- Sensitivity evaluation leaves the primal result and its immutability untouched; concurrent evaluations in separate sensitivity contexts give the same results as serial ones, bit for bit (ADR-0007 G1).
- The `double` evaluation benchmark of ADR-0003 stays within its pass criteria after the evaluators become scalar-generic.

Finite differences are verification tools, not the definition of a derivative. Tests document units, scales, tolerances and the regime in which perturbed samples keep their correspondence. Higher-order derivatives, generalised derivatives at events and certified step bounds are deferred.

## References

These inform the design; their libraries and algorithms are not mandated dependencies.

- [DLR adOCCT](https://github.com/dlr-sp/adOCCT): algorithmic differentiation applied to a CAD kernel (OCCT).
- Blondel et al., [Efficient and Modular Implicit Differentiation](https://arxiv.org/abs/2105.15183): differentiating solutions through their defining conditions.
- [SU2 software components](https://su2code.github.io/docs_v7/Software-Components/): separating analysis sensitivities from their projection onto geometric design variables.
