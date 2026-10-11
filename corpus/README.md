# Geometry corpus

The corpus is how Stereon is judged: a collection of geometry cases, each with an operation and a way to tell whether the result is right. Every case runs in CI once the operation it exercises exists. Cases are added through pull requests or the *Corpus case* issue template.

## Rules

- **Regression cases are never deleted.** A case that once exposed a defect stays in the corpus permanently. If its format becomes obsolete, it is converted, not removed.
- **Every case states how correctness is judged**: analytic ground truth, agreement with another kernel, or "valid result" at minimum.
- **Only redistributable files.** Cases authored for Stereon are MPL-2.0 (covered by `REUSE.toml`). Anything else needs a row in [`LICENSES.md`](LICENSES.md) naming the case, its source and its licence (ADR-0008).
- **Inputs are doubles.** Numbers in case files are read as the IEEE-754 double nearest to the decimal text, as `strtod` and every JSON parser do. Ground truth refers to those doubles, not to the decimal text. For example, `0.3 - 0.2 < 0.1` holds for the doubles, so a decimal-looking tangency may be a small crossing; see `analytic-ssi/plane-sphere/decimal-tangent-double-circle`.
- `python tools/check_corpus.py` validates every case and recomputes ground truth where an oracle exists. It runs in CI.

## Layout

```
corpus/
  README.md              this file
  LICENSES.md            source and licence of every case not authored for Stereon
  <suite>/               one directory per purpose, e.g. analytic-ssi/
    README.md            what the suite tests, its families, its inputs and expected results
    <family>/            cases that test the same thing
      <case>.json        case metadata, and the inputs when they are small
      <case>.<ext>       attached input files (.step, .stn, journal), same stem as the case
```

The case ID is the path below `corpus/` without the extension, e.g. `analytic-ssi/plane-sphere/tangent-exact`. Use lower-case kebab-case names that describe the configuration, not a number.

| Suite | Purpose | Spec |
| --- | --- | --- |
| [`analytic-ssi/`](analytic-ssi/) | Plane, sphere and cylinder intersections under tolerance | [RS-0001](../docs/research/RS-0001-tolerant-analytic-intersection.md), ADR-0010 |

## Case format `stereon-corpus-case/1`

Each case is a JSON object with these fields:

| Field | Required | Meaning |
| --- | --- | --- |
| `format` | yes | Always `"stereon-corpus-case/1"` |
| `id` | yes | The case ID; must equal the file path below `corpus/` |
| `kind` | yes | `regression` (a defect was found with it; never deleted), `hard` (a degenerate or tricky situation), or `reference` (a typical case for coverage and benchmarks) |
| `summary` | yes | One sentence: what the case is and why it is interesting |
| `added` | yes | Date added, `YYYY-MM-DD` |
| `source` | yes | Who made it and where it came from |
| `licence` | yes | SPDX identifier; anything other than `MPL-2.0` needs a row in `LICENSES.md` |
| `references` | no | ADRs and research specs the case tests, e.g. `["RS-0001", "ADR-0010"]` |
| `tags` | no | Free-form labels for selecting cases, e.g. `tangency`, `snap` |
| `inputs` | yes | The inputs. The suite README defines the shape; attached files are listed in `inputs.files` |
| `operation` | yes | What to do with the inputs: `{"name": ..., "operands": [...]}` plus operation parameters |
| `expected` | yes | How to judge the result; the suite README defines the shape |
| `found` | for regressions | The defect: Stereon commit, issue number, and the wrong result or error that was observed |

A suite may add fields inside `inputs`, `operation` and `expected`, never at the top level. A change to the top-level fields is a new format version, and `check_corpus.py` must convert or accept every existing case before the old version is dropped.

## Workflow for a failing case

1. Reduce the failure to the smallest input that still fails, and save it with its operation journal if there is one.
2. Add it as a `regression` case with `found` filled in, in the suite and family that match.
3. The fix and the case land together; the case keeps running forever.
