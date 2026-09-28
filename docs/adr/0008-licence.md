# ADR-0008: Licence

- **Status:** Proposed
- **Date:** 2026-09-28
- **Owner:** Onur Tuncer
- **Phase:** 0

## Context

The licence decides who can use Stereon, how improvements flow back, and whether industrial users can embed it. It is very hard to change after outside contributions arrive, so it must be settled before the repository goes public.

## Options

1. **MIT / Apache-2.0 (permissive).** Maximum adoption; modified versions of the kernel may be kept closed.
2. **LGPL-2.1 (OCCT's licence).** Library-level copyleft; its relinking requirements are awkward for C++ templates, header code and static linking.
3. **GPL-3.0.** Strong copyleft; prevents use in closed applications.
4. **MPL-2.0.** File-level copyleft: modified Stereon files must stay open, while Stereon can be combined with code under other licences in a larger work. Includes an explicit patent grant and is GPL-compatible.

## Decision

Adopt MPL-2.0 for all kernel source.

- Every source file carries the MPL-2.0 header.
- Contributions are accepted under the contributor agreement described in `CONTRIBUTING.md`.
- Third-party dependencies must be compatible with MPL-2.0 (MIT, BSD, Apache-2.0, MPL-2.0, Boost). Copyleft dependencies are allowed only as optional, separately built adapters.
- Test corpus files carry their own source licences, recorded in `corpus/LICENSES.md`.

## Consequences

- Improvements to kernel files come back to the project; applications built on Stereon remain free to choose their own licence.
- Clear rules for templates and static linking, unlike LGPL.
- Dependency choices must be checked for licence compatibility in review.

## Verification

A CI licence check (e.g. REUSE or a header scanner) confirms every source file carries the header and every dependency is on the allowed list.
