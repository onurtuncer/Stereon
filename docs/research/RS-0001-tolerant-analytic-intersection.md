# RS-0001: Tolerance-aware intersection and topology for analytic surfaces

- **Status:** Open (go/no-go pending)
- **Date:** 2026-10-11
- **Owner:** Onur Tuncer
- **Phase:** 0 side track. The spec, corpus and reference oracle need no kernel code; the kernel implementation waits for the Phase 0 topology prototype
- **Related:** ADR-0001, ADR-0002, ADR-0006, ADR-0007, ADR-0009, ADR-0010, ADR-0013, ADR-0018
- **Corpus:** [`corpus/analytic-ssi/`](../../corpus/analytic-ssi/)

## Question

Can intersection curves and the topology built from them be computed for a narrow set of analytic surfaces (planes, spheres and cylinders) when every input carries a modelling tolerance? And can the method guarantee that the whole output is consistent with **one** admissible placement of the inputs, not only pair by pair?

## Motivation

ADR-0010 assumes that the two surfaces are exact and that the contact type of a branch can be decided. ADR-0001 gives every vertex and edge a local tolerance and defines a merge policy for vertices and edges. Neither ADR says what intersection means when the *surfaces* are only known to within a tolerance. That is the normal state of imported data, of loosened edges, and of design intent ("this cylinder rests on that plane").

The difficulty is that the degenerate configurations that matter most in modelling, such as tangency and coincidence, never survive a perturbation:

- A sphere and a plane at distance `r + 0.5τ` are disjoint as numbers but tangent as intent. Under tolerance, either answer is defensible.
- Decisions made independently for each pair can contradict each other. If A is snapped to B and B to C, but A and C are more than a tolerance apart, no single placement of A, B and C realises all three answers. ADR-0001's merge policy makes this deterministic for vertices; nothing makes it consistent for surfaces.
- A topology built from inconsistent pair decisions can pass every local check and still be unrealisable, which is the silent-failure mode ADR-0006 forbids.

Kernels in production handle this with per-pair fuzzy values. OCCT's Boolean operations, for example, take a fuzzy value (`BOPAlgo_Options::SetFuzzyValue`). To our knowledge, these do not guarantee joint realisability. Establishing whether that is true is part of the novelty check below.

## Scope

**Surfaces.** Unbounded `Plane`, `Sphere` and `Cylinder`, given by points, directions and radii as IEEE-754 doubles. Direction vectors need not be unit length, because all decisions are made on the input doubles and nothing is normalised.

**Tier A: in scope for the go/no-go decision.** Pairs whose intersection is always a set of lines, conics or points, or a coincidence:

| Pair | Nominal configurations |
| --- | --- |
| plane–plane | disjoint (parallel), one line, coincident |
| plane–sphere | disjoint, tangent point, circle |
| plane–cylinder | disjoint, tangent line, two lines (axis ∥ plane); circle (axis ⊥ plane); ellipse (oblique) |
| sphere–sphere | disjoint, external or internal tangent point, circle, coincident |
| sphere–cylinder, centre on axis | disjoint, tangent circle, two circles |
| cylinder–cylinder, parallel axes | disjoint, external or internal tangent line, two lines, coincident |

**Tier B: kept out until the research question requires it.** Pairs whose intersection is a quartic space curve: sphere–cylinder off axis (Viviani's curve is a degenerate member) and cylinder–cylinder with skew or intersecting axes (two ellipses when the radii are equal and the axes intersect). The corpus contains Tier B cases so that the boundary of the method is tested. Inside the scope of this study, Tier B returns `NotSupported`.

Tier membership is itself a decision. Rotations are not admissible perturbations, so whether two axes are parallel is decided exactly on the input doubles: nearly parallel cylinders are Tier B. A sphere whose centre lies within tolerance of a cylinder's axis is snapped onto the axis (G3) and treated as Tier A.

**Arrangements.** Three or more Tier A surfaces whose pair decisions interact. These carry the joint-consistency claim, and they are the cases where a per-pair method is expected to fail.

**Out of scope:** cones, tori, NURBS, rotational perturbations (bounded by the angular resolution of ADR-0001), and bounded faces. Solid Booleans on boxes, spheres and cylinders belong to the ADR-0001 prototype, which uses this method once it exists.

## Novelty

The method contributes only if all three statements below survive the literature survey. If one does not, the study narrows to the statements that remain or stops (see the go/no-go criteria).

1. **Joint realisability for curved analytic surfaces.** The output (every pair configuration, every branch's contact type, every vertex incidence) is the exact intersection of one admissible placement of the inputs, and that placement is returned as the output geometry. Earlier work gives this for polyhedra (Segal 1990; Fang, Brüderlin and Zhu 1993) or reasons about imprecision without constructing a placement (epsilon geometry, Salesin, Stolfi and Guibas 1989). Production kernels decide pair by pair.
2. **Exact decisions on the input doubles.** Every configuration decision for a Tier A pair reduces to signs of polynomials in the input doubles. These are evaluated with interval filters and an exact fallback (ADR-0002), so the decision never depends on round-off, thread count or evaluation order.
3. **Explicit failure semantics.** The method separates three outcomes: decided robustly (every admissible placement agrees), decided by snapping (a degenerate configuration lies within tolerance and is chosen), and not realisable within the healing budget (`ToleranceExceeded`). Per-pair fuzzy methods merge the first two and have no third.

**Prior art to survey before the go/no-go decision.** Citations are from memory, and each must be checked against the source:

- Salesin, Stolfi, Guibas, *Epsilon geometry: building robust algorithms from imprecise computations*, SoCG 1989.
- Segal, *Using tolerances to guarantee valid polyhedral modeling results*, SIGGRAPH 1990.
- Fang, Brüderlin, Zhu, *Robustness in solid modelling: a tolerance-based intuitionistic approach*, Computer-Aided Design, 1993.
- Jackson, *Boundary representation modelling with local tolerances*, ACM Symposium on Solid Modeling, 1995.
- Hu, Patrikalakis, Ye, *Robust interval solid modelling*, parts I and II, Computer-Aided Design, 1996.
- Miller, Goldman, *Geometric algorithms for detecting and calculating all conic sections in the intersection of any two natural quadric surfaces*, Graphical Models and Image Processing, 1995.
- Sugihara, Iri, topology-oriented robust geometric algorithms.
- Patrikalakis, Maekawa, *Shape Interrogation for Computer Aided Design and Manufacturing*: the chapters on intersection and robustness.
- Current practice: OCCT fuzzy Booleans and `IntAna_QuadQuadGeo`, and published descriptions of tolerant modelling in Parasolid and ACIS.

The survey answers one question per statement: has this guarantee already been given for curved surfaces, and if so, by which method and with what limits?

## Contract

This section is the proposed behaviour of the method. If the decision is go, it becomes the basis of revisions to ADR-0001, ADR-0010 and ADR-0013.

### Admissible placements

Each input surface `S` with tolerance `τ` (resolution ≤ τ ≤ ceiling, ADR-0001) may be replaced by any surface of the same type obtained by these perturbations:

| Surface | Admissible perturbation |
| --- | --- |
| Plane | offset along its normal by at most `τ` |
| Sphere | centre displacement `Δc` and radius change `Δr` with `‖Δc‖ + |Δr| ≤ τ` |
| Cylinder | axis displacement `Δa` perpendicular to the axis and radius change `Δr` with `‖Δa‖ + |Δr| ≤ τ` |

Each perturbation moves every point of the surface by at most `τ`. Rotations are excluded in Tier A (see open question O2).

### Guarantees

- **G1 Realisability.** Every decided output is the exact intersection of one admissible placement of all inputs, and the output geometry *is* that placement. A snapped surface is stored as its perturbed form, not as the nominal one.
- **G2 Robust decisions where possible.** If every admissible placement of a pair gives the same configuration, that configuration is returned and nothing is snapped.
- **G3 Snapping toward degeneracy.** If a degenerate configuration (tangency, coincidence, centre on axis) is admissible for a pair, the degenerate configuration is returned. If several are admissible, the most degenerate one wins: coincidence over tangency. This encodes modelling intent: geometry placed within tolerance of touching is treated as touching.
- **G4 Joint consistency.** In an arrangement, snaps are made in the fixed order of the ADR-0001 merge policy (input index, then pair index). A snap is accepted only if the placement after it is still admissible for every surface. Otherwise the pair keeps its robust configuration if one exists. If none exists, the operation grows the tolerance within the healing budget (logged) or fails.
- **G5 Exact, deterministic decisions.** Every configuration decision is a sign of a polynomial in the input doubles and the tolerances, evaluated by `stereon::robust` (ADR-0002). There are no literal epsilons, and the result does not depend on thread count (ADR-0007).
- **G6 Typed, certified output.** Each branch carries its contact type (`Crossing`, `Tangential`, `Overlap`; ADR-0010, ADR-0013) and becomes an `Exact` edge (ADR-0009) whose certificate bounds its distance from both of the *placed* surfaces.

### When the method must report uncertainty or failure

| Situation | Outcome |
| --- | --- |
| A Tier B pair, or a surface type outside the scope | `KernelError::NotSupported`, naming the pair |
| The arrangement has no admissible placement realising a consistent configuration, and growth within the healing budget does not create one | `KernelError::ToleranceExceeded`, naming the conflicting surfaces and the growth that would be required |
| A snap or tolerance growth was applied | Success, with every snap and every growth listed in `OperationResult` (ADR-0006) |
| A certificate cannot be met | `KernelError::CertificationFailed` (ADR-0009). For Tier A this should never happen; if it does, it is a defect |

**Never:** a result that is not the intersection of an admissible placement, or a topology that the validity checker (ADR-0018) rejects.

An `Undecided` outcome is deliberately absent. For a single Tier A pair, the admissible placements form a connected set, and the configuration changes only across a degenerate locus. So if two configurations are admissible, a degenerate one is admissible too, and G2 or G3 always applies. Ambiguity appears only in arrangements, and there G4 resolves it or the operation fails. If the study finds a case that neither resolves, the contract needs an `AmbiguousConfiguration` outcome. That would require an ADR-0006 revision.

### Open questions

- **O1 Pair tolerance.** ADR-0001 says entity pairs use the larger of their two tolerances. Under the placement semantics above, both surfaces may move, which suggests the sum. The study must pick one and state the consequence for ADR-0001.
- **O2 Rotations.** Excluding rotations keeps every decision polynomial and simple. A plane tilted within the angular resolution moves points far from its origin by more than `τ` inside the size box, so rotations probably cannot simply be added to the admissible perturbations. Should they be admitted for Tier B?
- **O3 Snap priority.** G4 snaps in input order. It is not yet known whether an order exists that maximises the number of degenerate configurations kept (intent preserved) and can be computed in polynomial time.
- **O4 Snap or robust.** G3 prefers degeneracy whenever it is admissible. Some callers will want the opposite (nominal behaviour), so the preference may need to be a `Context` setting.

## Evidence

### Corpus

[`corpus/analytic-ssi/`](../../corpus/analytic-ssi/) holds the cases. Every case records two expected results:

- **`nominal`:** the exact configuration of the input doubles taken as exact numbers. This is mathematics with one right answer, and [`tools/check_corpus.py`](../../tools/check_corpus.py) recomputes it in rational arithmetic for every Tier A pair case.
- **`tolerant`:** the configuration this contract requires. It is the claim under test and is marked as proposed until the contract is accepted.

Case families, each swept over offsets `δ/τ ∈ {0, 0.1, 0.5, 0.9, 0.99, 1.01, 1.1, 2, 10}` from the degenerate configuration, over tolerances `τ ∈ {1e-5, 1e-4, 1e-3, 1e-2} mm`, and over placements near the origin and near the edge of the size box:

| Family | What it tests |
| --- | --- |
| Pair degeneracies (Tier A table) | Robust versus snapped decisions, G2 and G3 |
| Non-representable values (radii and offsets such as 0.1 or 1/3) | Decisions on the doubles, not the decimals, G5 |
| Chains A ≈ B ≈ C with A ≉ C (planes, spheres, coaxial cylinders) | Joint consistency, G4, and the non-transitive case of ADR-0001 |
| A surface tangent to two others at once (a sphere in a corner, a cylinder on two planes) | Snaps that compete for the same surface |
| Tier B near-degenerate pairs | The `NotSupported` boundary |

The seed cases committed with this spec cover each family once. The parametric sweeps are generated, not hand-written; the generator is a follow-up item.

### Baselines

1. **OCCT** `IntAna_QuadQuadGeo` for pairs, and `BOPAlgo` with and without a fuzzy value for arrangements, from the pinned OCCT build of [`docs/benchmarks.md`](../benchmarks.md).
2. **Nominal ablation.** The same exact predicates without snapping: intersect the input doubles exactly, then merge vertices per ADR-0001. This isolates what snapping and joint consistency contribute, independent of OCCT.
3. **Per-pair snapping ablation.** G3 without G4. If this matches the full method on every case, statement 1 of the novelty section has no evidence.

### Metrics

For every method and case: configuration agreement with `nominal` and with `tolerant`; whether the output is realisable (checked by the oracle, which verifies that the returned placement is admissible and that the returned configuration is its exact intersection); the robustness score of `docs/benchmarks.md` (valid, structured error, silently invalid); the total snap displacement and tolerance growth; and time per pair.

### Reference oracle

The go/no-go decision does not wait for the kernel. A reference implementation in Python with exact rational arithmetic (`fractions.Fraction`, so every double is converted exactly) implements the contract for Tier A pairs and arrangements. It is slow and that is acceptable: it produces the `tolerant` ground truth and checks realisability for every method. The kernel implementation is later measured against it.

## Go/no-go criteria

The decision is due at the end of Phase 0 (week 12) and is recorded in this file.

**Go** if all of the following hold:

1. **Novelty.** The survey finds no prior method that gives guarantee G1 for curved analytic surfaces with exact decisions. If one exists, a measurable improvement over it remains.
2. **Feasibility.** The reference oracle decides every Tier A case of the seed corpus and the generated sweeps, with every decision reduced to polynomial signs (no case needs a numeric tolerance of its own).
3. **Evidence.** On at least one arrangement family, a baseline (OCCT fuzzy or the per-pair ablation) produces unrealisable or invalid results on a reproducible fraction of the cases, while the oracle produces none. The cases must be ones a CAD user would recognise: resting contact, a part in a corner, stacked plates.

**No-go** if criterion 1 or 3 fails. The study then closes with:

- the corpus kept as the analytic seed of the ADR-0010 SSI corpus and of benchmark B2.3;
- the contract section turned into an ADR-0010 revision describing the best existing method, adopted as Stereon's analytic intersection behaviour;
- the oracle kept as the ground-truth tool for analytic SSI tests.

The work is useful to Stereon whichever way the decision goes.

## Deliverables

- [x] This specification
- [x] Corpus layout, case schema and `corpus/LICENSES.md` ([`corpus/README.md`](../../corpus/README.md))
- [x] Seed cases for every family in [`corpus/analytic-ssi/`](../../corpus/analytic-ssi/)
- [x] Nominal oracle for Tier A pairs in `tools/check_corpus.py`, run in CI
- [ ] Literature survey answering the three novelty statements
- [ ] Generator for the parametric sweeps
- [ ] Reference oracle for the `tolerant` contract, including arrangements
- [ ] Baseline runs (OCCT, nominal ablation, per-pair ablation)
- [ ] Go/no-go decision recorded here

## Decision

Pending.
