# Claude Code agents for the Stereon v0.2 release

Six reusable Markdown agent definitions live under [`.claude/agents/`](agents/). Claude Code recognizes each definition from YAML frontmatter. The repository's ordinary `main` branch and accepted ADRs remain authoritative; agent instructions do not change kernel contracts.

| Role | Agent name | Main ownership |
| --- | --- | --- |
| F — lead | `stereon-v02-coordinator` | Dependencies, contracts, PR integration and release decision |
| A | `stereon-v02-ci` | Linux/Windows/Lint/Docs CI |
| B | `stereon-v02-ssi` | Analytic and generic surface–surface intersection |
| C | `stereon-v02-curves` | 2D/3D curve pairs and curve–surface edge cases |
| D | `stereon-v02-sensitivities` | Intersection linearization, JVP/VJP and `ClosestApproach` |
| E | `stereon-v02-verification` | Metamorphic backfill, independent oracles, OCCT benchmarks and release evidence |

## Preferred mode: distinct worktree sessions

Start each command **in its own terminal**, from a trusted Stereon checkout with the committed agent definitions available. Use `claude --version` to confirm Claude Code is installed. If you are developing from a PR that adds these files, check out that branch locally or merge the PR first.

```bash
# Coordinator — read and manage task allocation / PR sequencing
claude --agent stereon-v02-coordinator

# One terminal/session per implementation specialist
claude --worktree v02-ci           --agent stereon-v02-ci
claude --worktree v02-ssi          --agent stereon-v02-ssi
claude --worktree v02-curves       --agent stereon-v02-curves
claude --worktree v02-sensitivities --agent stereon-v02-sensitivities
claude --worktree v02-verification --agent stereon-v02-verification
```

`--worktree` creates isolated sessions; ask each agent to confirm its real worktree, branch, base SHA and dirty state **before** editing. Claude Code's `isolation: worktree` frontmatter similarly isolates agents when invoked as ordinary subagents, but **do not rely on that field for a named Agent Team teammate**: in-process teammates can share the lead's checkout. Independently launched worktree sessions synchronize through PRs, not an automatically shared task list.

Only one agent may edit a shared source or test file at once. Changes to result contracts, common intersection helpers or public interfaces require prior coordinator allocation.

## Alternative: experimental Agent Teams for planning and reviews

Agent Teams can use the same agent definitions as teammate types. They have shared tasks and direct teammate messaging, but **they do not by themselves establish separate file-system worktrees**. Use them for read-only audits, planning, independent analyses or strictly nonoverlapping edits; keep concurrent implementation in distinct worktree sessions.

Enable this experimentally via a **local personal** `~/.claude/settings.json` (merge with existing content, never replace existing settings):

```json
{
  "env": {
    "CLAUDE_CODE_EXPERIMENTAL_AGENT_TEAMS": "1"
  }
}
```

Then launch the coordinator and request teammate types **by name**. The team lead is the main session; do not assume a special `.claude/teams/teams.json` file has any effect.

## Initial coordinator assignment

Paste the following into the coordinator session:

> Inspect Stereon main's current commit, supported-domain v0.2 ledger, roadmap, open PRs, actions and ADR status. Build an owner/PR/dependency/evidence board. Confirm that #262 is the prerequisite for stacked #263; check both live states. Keep PR #264 as documentation work unless merged. Assign A blocking CI; B difficult generic SSI; C curve intersection edge cases; D JVP/VJP and ClosestApproach; E independent verification. Give each implementation agent exclusive file ownership, a minimal independent oracle, direct and applicable metamorphic tests, a supported-domain contract, and native CI obligations. Do not edit the same checkout concurrently, alter accepted ADRs, mark an unresolved case successful, or merge unqualified PRs. Report the critical path and release gate status before integration.

## Non-negotiable release rules

1. **Domain and assurance first:** `empty` is only successful after declared-domain completeness is established; unresolved tangencies fail closed.
2. **Separate geometry from topology:** an intersection witness does not certify a B-rep reconstruction; a near tangency alone is not a topology event.
3. **No silent numerical weakening:** accepted tolerance policies, exact/bounded/estimated distinctions, budgets and error statuses must be honored.
4. **Preserve derivatives' mathematics:** differentiable regular branches are not permission to invent smooth derivatives through topology changes; ADR-0052 defines the space-curve meeting exception.
5. **A repaired bug needs proof:** direct reproducer, independent oracle, metamorphic relation when applicable and red-on-defect/green-on-repair evidence.
6. **Qualification is revision-bound:** native GCC16 platforms, B2.1–B2.3 and other applicable v0.2 release gates require actual CI/corpus/benchmark evidence. Do not tag from a planning document.

## Current dependency hints, not frozen status

At initial preparation the relevant PRs were:

- [#262](https://github.com/onurtuncer/Stereon/pull/262) — tolerance bounds and verification evidence infrastructure (draft at inspection).
- [#263](https://github.com/onurtuncer/Stereon/pull/263) — metamorphic fault backfill, stacked on #262 (draft at inspection).
- [#264](https://github.com/onurtuncer/Stereon/pull/264) — computational-geometry assurance roadmap and ADR-0055 (documentation PR at inspection).

**Always recheck GitHub and the current base SHA.** A PR, check failure, implementation claim or ADR status may change after these instructions were committed.

## Reviews and outputs

Implementation agents submit narrow PRs with exact tests, mathematical reasoning, verified CI links and explicit limitations. The coordinator decides merge order and writes the release-gate table. The verification agent independently reviews evidence rather than merely counting passing tests.

Useful authoritative files:

- [v0.2 acceptance ledger](../docs/planning/v0.2-acceptance-ledger.md)
- [Release roadmap](../docs/planning/roadmap-and-versioning.md)
- [Metamorphic backfill](../docs/planning/intersection-metamorphic-backfill.md)
- [Phase 2 benchmarks](../docs/benchmarks.md)
- [ADR index](../docs/adr/README.md)

No CI pass, PR merge or v0.2 release is implied by adding these agent definitions.
