# Research specifications

This directory holds Stereon's research questions: open problems where the right design is not known yet and has to be found out before anything can be decided. An ADR records a decision. A research specification (RS) records a question, the claim it tests, and the evidence that will settle it. It ends in a go/no-go decision.

## Rules

- One file per question, numbered sequentially: `RS-NNNN-short-title.md`.
- Each RS has these sections: **Question**, **Scope**, **Novelty** (what it adds beyond existing methods, with the prior art to check), **Contract** (what the method guarantees, and when it must report uncertainty or failure), **Evidence** (the corpus, baselines and metrics), **Go/no-go criteria**, **Deliverables** and **Decision**.
- Status moves from **Open** to **Go** or **No-go** when the decision is recorded in the file. An RS whose follow-up work is finished becomes **Closed**.
- An RS never changes the kernel's rules directly. A go decision leads to new ADRs or revisions of existing ones. A no-go decision records what was kept, such as the corpus or a reference implementation, and where it went.
- Evidence must be reproducible: cases live in [`corpus/`](../../corpus/), and every number in a decision names the commit and the cases it came from.

## Index

| RS | Title | Related ADRs | Status |
| --- | --- | --- | --- |
| [0001](RS-0001-tolerant-analytic-intersection.md) | Tolerance-aware intersection and topology for analytic surfaces | 0001, 0009, 0010, 0013 | Open |
