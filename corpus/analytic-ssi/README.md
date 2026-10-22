# Analytic surface–surface intersection

Pairs and arrangements of planes, spheres and cylinders, given with modelling tolerances. The suite is the evidence for [RS-0001](../../docs/research/RS-0001-tolerant-analytic-intersection.md). Whatever the go/no-go decision, it also seeds the analytic part of the ADR-0010 SSI corpus and of benchmark B2.3.

## Families

| Family | Contents |
| --- | --- |
| `plane-plane/`, `plane-sphere/`, `plane-cylinder/`, `sphere-sphere/`, `sphere-cylinder/`, `cylinder-cylinder/` | Tier A pairs: exact degeneracies, near-degeneracies on both sides of the tolerance threshold, non-unit directions, non-representable decimals, placements near the edge of the size box |
| `tier-b/` | Pairs whose intersection is a quartic, and near-Tier-A pairs that are Tier B as numbers |
| `arrangement/` | Three or more surfaces whose pair decisions interact: non-transitive chains and competing snaps |

The current cases are hand-written seeds, one or more per family. The parametric sweeps that RS-0001 describes will be generated.

## Inputs

`inputs.surfaces` is a list of named surfaces. All lengths are in millimetres. Direction vectors need not be unit length and are never normalised.

| `type` | Fields |
| --- | --- |
| `plane` | `point` (on the plane), `normal` |
| `sphere` | `centre`, `radius` |
| `cylinder` | `point` (on the axis), `axis`, `radius`; the cylinder is unbounded |

Every surface has a `tolerance` of at least the resolution, 1e-5 mm (ADR-0001). Within a seed case, all surfaces share one tolerance unless the case is about unequal tolerances.

## Operations

- `{"name": "intersect", "operands": ["A", "B"]}`: the configuration of one pair.
- `{"name": "arrange", "operands": ["A", "B", "C", ...]}`: the configuration of every pair, decided jointly.

## Expected results

`expected.nominal` is the exact configuration of the input doubles. `tools/check_corpus.py` recomputes it in rational arithmetic and fails if a case disagrees.

- Pair: `{"tier": "A", "configuration": <configuration>}`, or `{"tier": "B", "configuration": "quartic"}`. For Tier B, only the tier is checked.
- Arrangement: `{"pairs": {"A-B": <configuration>, ...}}`, with one key for every pair in operand order. A Tier B pair is `"tier_b"`.

`expected.tolerant` is what the RS-0001 contract requires. It is the claim under test, so it is checked for shape only.

- `{"status": "proposed", "configuration": ..., "snapped": true|false}` for a pair. `snapped` says whether the configuration comes from snapping to a degeneracy (G3) rather than holding for every admissible placement (G2).
- `{"status": "proposed", "pairs": {...}, "snaps": ["A-C", ...]}` for an arrangement. `snaps` lists the accepted snaps in the order they were made.
- `{"status": "proposed", "outcome": "NotSupported"}` when the contract requires an error.
- `{"status": "open", "question": "O1"}` when the answer depends on an open question of RS-0001.

The seed cases assume these provisional rules until RS-0001 settles them:

1. A pair's tolerance is the larger of the two (ADR-0001). This answers RS-0001's O1, decided for the Boolean prototype on 2026-10-20 (decision D2 of `docs/planning/boolean-prototype-plan.md`); `plane-plane/unequal-tolerances` is the case where the choice changes the answer.
2. In a snap, the surface with the higher operand index moves, by the smallest admissible displacement. A plane moves along its normal. A sphere or cylinder moves its centre or axis before changing its radius.
3. Pairs are decided in operand order (`A-B`, `A-C`, `B-C`, …), and a snap is accepted only if every surface's total displacement stays within its tolerance (G4).

## Configurations

| Pair | Configurations |
| --- | --- |
| plane–plane | `disjoint`, `line`, `coincident` |
| plane–sphere | `disjoint`, `tangent_point`, `circle` |
| plane–cylinder | `disjoint`, `tangent_line`, `two_lines`, `circle`, `ellipse` |
| sphere–sphere | `disjoint`, `tangent_point_external`, `tangent_point_internal`, `circle`, `coincident` |
| sphere–cylinder (centre on axis) | `disjoint`, `tangent_circle`, `two_circles` |
| cylinder–cylinder (parallel axes) | `disjoint`, `tangent_line_external`, `tangent_line_internal`, `two_lines`, `coincident` |

`disjoint` includes nested spheres and nested parallel cylinders. The contact type of each branch follows from the configuration: `tangent_*` is `Tangential`, `coincident` is `Overlap`, and every other branch is `Crossing` (ADR-0010, ADR-0013).
