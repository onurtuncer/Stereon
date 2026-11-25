# Boolean corpus

Booleans on solids made of planes and cylinders, with the result's topology known in closed form. `tests/boolean/CorpusTests.cpp` runs every case: the cases of 3 by 3 grids or smaller in the routine CI matrix, the larger ones in the *Boolean issue 160 large grids* workflow (`.github/workflows/issue160-large.yml`), on every pull request that touches the Boolean, on main and nightly.

No operation records a journal yet (ADR-0023), so the cases carry no `.journal`; their inputs are generated in full from the parameters below.

## Families

| Family | What it tests |
| --- | --- |
| `drill-grids/` | Issue #160: grids of parallel cylinders, disjoint or overlapping in one column, with caps level or offset |

## `inputs`

`inputs.operands` is a list of named operands. Each has a `type`:

- `drill-grid`: `n` by `n` cylinders of radius 2 and length 20, axis +z, the base of cylinder (i, j) at (10 i + 5, 10 j + 5, -5) + `offset`, for i, j in 0 .. n - 1; made with `MakeCylinderSolid` and united, on the serial scheduler, as a balanced tree of unions (`BooleanCases::Drills` in `tests/support/BooleanCases.hpp`). Neighbouring drills' axes are 10 apart and their radii 2, so the grid is n² disjoint solids of three faces each.

## `operation`

- `{"name": "build", "operands": ["A"]}`: the operand itself, whose making is the Boolean under test.
- `{"name": "union" | "common" | "cut", "operands": ["A", "B"]}`: the Boolean of the two operands on the serial scheduler.

## `expected`

- `solids`: the number of solids of the result.
- `faces` (optional): the number of faces of the result.
- `valid`: `true`; `ValidateTopology` finds the result valid, and at `ValidityLevel::Solid` reports nothing but `SelfIntersectionUnchecked`.

Two drills whose axes are 1 apart overlap and unite into one solid; drills whose axes are 10 or more apart stay separate. So two n by n grids overlapping in one column unite into 2n² - n solids, and two disjoint grids into 2n².
