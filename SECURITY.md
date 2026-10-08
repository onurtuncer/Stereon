# Security Policy

## Supported versions

Stereon is pre-alpha (Phase 0). No release is supported for production use yet. Security fixes are made on `main` only.

| Version | Supported |
| --- | --- |
| `main` | ✅ |
| Anything else | ❌ |

This table will be updated when the first release (v0.1) is published.

## Scope

Stereon reads geometry from files that may come from untrusted sources. Problems we treat as security issues include:

- Memory-safety bugs (out-of-bounds access, use-after-free, buffer overflows) that can be triggered by input data, especially through the STEP, STL, OBJ, glTF and native-format readers.
- Unbounded memory or CPU use caused by small malicious inputs (for example, decompression bombs or pathological knot vectors).
- Problems in the Python bindings that let crafted input escape the kernel's API contracts.
- Compromise of the build, CI or release pipeline, or of fetched dependencies.

These are ordinary bugs. Report them as public issues:

- Wrong geometry, invalid shapes or failed Boolean operations on well-formed input.
- Crashes that require misusing the C++ API, such as violating a documented precondition (contract).
- Problems found only in debug or sanitizer builds that do not affect release builds.

If you are unsure, report it privately. We would rather reclassify a report than have a real issue disclosed in public.

## Reporting a vulnerability

**Do not open a public issue, discussion or pull request for a security problem.**

Report it privately through GitHub: go to the repository's **Security** tab and choose **Report a vulnerability**. This opens a private advisory that only the maintainers can see.

If you cannot use GitHub's private reporting, email the maintainer, Prof. Dr. Onur Tuncer (Department of Aeronautical Engineering, Istanbul Technical University), at onur.tuncer@itu.edu.tr. Put `[Stereon security]` in the subject line.

Please include:

- A description of the issue and its impact.
- The affected component (for example, `io/step`) and the commit you tested.
- A minimal input file or code that reproduces it, and the compiler, platform and build preset you used.
- Any sanitizer output or stack trace.

## What to expect

- We acknowledge your report within **7 days**.
- We confirm or reject the issue and share a planned fix timeline within **30 days**.
- We develop the fix in a private fork and add a regression test, normally to the fuzzing corpus.
- We publish a GitHub security advisory when the fix is merged, crediting you unless you prefer to stay anonymous.

Stereon is maintained by volunteers. These are targets, not guarantees, but we will keep you informed if we need more time.

Please give us a reasonable chance to release a fix before you disclose the issue publicly. We aim to coordinate disclosure within 90 days of your report.
