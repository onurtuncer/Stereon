# ADR-0012: Native file format and versioning

- **Status:** Proposed (to be accepted before Phase 3)
- **Date:** 2026-09-28
- **Owner:** Onur Tuncer
- **Phase:** 3

## Context

STEP is the exchange format, but it is slow to parse, loses kernel-specific data (local tolerances, history, attributes) and cannot store regression cases exactly. The kernel needs a native format that round-trips a `Shape` bit-exactly, loads fast, and survives format evolution for years, because the regression corpus depends on it.

## Options

1. **Text format** (JSON/XML). Readable and diffable; large and slow.
2. **Third-party binary schema** (FlatBuffers, Cap'n Proto, Protobuf). Mature tooling; an external schema to keep in sync with the C++ types.
3. **Own chunked binary format with reflection-generated serialisers.** Direct mapping of the SoA arrays; C++26 reflection removes hand-written boilerplate.

## Decision

Adopt option 3, plus a JSON debug dump.

- File extension `.stn`. Little-endian, chunked layout: header (magic, format version, kernel version, resolution and units), then one chunk per array (geometry pools, each topology entity type, attributes, history), each with type tag, length and checksum.
- Chunks are optionally compressed (zstd).
- Serialisers are generated with static reflection over the geometry and topology structs (ADR-0016); hand-written code is allowed only for version migrations.
- **Versioning.** Readers accept all older format versions via explicit migration functions; unknown chunks are skipped with a warning. The format version increments on every schema change.
- A `stn-dump` tool writes a JSON view of any file for debugging and diffs.

## Consequences

- Exact round-trips make the regression corpus trustworthy.
- Loading is close to memory bandwidth for large models.
- Reflection support is required from the reference compiler; a hand-written fallback may be needed while Clang/MSVC lag.
- Migrations must be maintained forever once the format is public.

## Verification

Round-trip every corpus shape through write/read and compare arrays byte for byte. Load time for a 10⁶-face body < 1 s. Files written by each tagged release stay readable by all later releases (CI keeps one file per release).
