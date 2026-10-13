.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

B-rep topology
==============

The topology store, typed handles, ShapeBuilder and immutable Shape, Euler operators, and the validity checker.

Handles
-------

Entities are addressed by typed 32-bit handles (ADR-0004). ``VertexId``,
``EdgeId``, ``CoedgeId``, ``LoopId``, ``FaceId``, ``ShellId`` and ``SolidId``
are the instances of ``Id<EntityKind>``, one type per kind, so passing an
edge handle where a vertex is expected, or a raw integer anywhere, is a
compile error. A default-constructed handle is invalid; ``IsValid()`` says so.
Handles compare and hash by index, so they key ``std::map`` and
``std::unordered_map`` side tables, and they format as ``Vertex#12``.

A ``ShapeBuilder`` hands out ``BuilderId<Kind>``, an index with a
generation. ``GenerationTable<Kind>`` owns the generations, the alive bitset
and the free list: releasing a slot bumps its generation, so every handle to
the released entity is stale and ``IsAlive()`` says no. Slots released in a
session are not reused in that session, which keeps an index that is live
before and after an edit on the same entity; ``BeginSession()`` makes them
reusable. Builder handles format as ``Vertex#12@3``.

``Oriented<Kind>`` carries a sense with the handle, as a loop uses a coedge
forwards or reversed. The sense lives in the top bit of the index and is
reached only through ``BaseOf(h)``, ``IsReversed(h)`` and ``Reversed(h)``;
nothing masks it by hand.

With ``STEREON_CHECKED_HANDLES`` (the debug and sanitizer presets) every
handle also records its store and generation, and ``GenerationTable::Check``
verifies store, bounds and generation. A mismatch prints a report naming the
handle, the entity kind and the expected store, then raises a contract
violation. Release builds compile ``Id`` down to its 4-byte index, which a
``static_assert`` guards. Handles are local to one store; identity across
operations comes from history maps (ADR-0011).

Storage pages
-------------

Entity tables and geometry pools keep their entries in ``PagedArena<T, N>``
(ADR-0005): fixed-size pages of ``N`` entries, 256 unless a table chooses
otherwise, behind ``std::shared_ptr``. Copying an arena copies the page
pointers, so a snapshot of a million-entry table costs a few thousand pointer
copies and shares every page with its source. The first write to a page
copies the page when another arena still holds it and writes in place when
the arena holds it alone; memory grows with what an edit touched, not with
the size of the table. ``PushBack`` appends, ``Resize`` grows with
value-initialised entries or drops whole pages, and ``GetPageSpan(page)``
hands hot loops one page at a time, so the page lookup and the reference
counting stay out of inner loops.

Writes are scoped. ``Set(index, value)`` replaces an entry, and
``Modify(index, fn)`` and ``ModifyPage(page, fn)`` call ``fn`` with a
reference to the entry or a span of the page, made unique first and valid
only for the duration of the call; whatever ``fn`` returns, by value, is
returned. No function returns a mutable reference. A reference kept past the
copy-on-write step would otherwise write into a page that a later copy
shares, and change a snapshot that is supposed to be immutable. While ``fn``
runs, the arena may be read but not copied, moved, assigned, cleared,
resized or edited again. Every build checks this with a plain branch and
reports a violation from a module implementation unit, so the check does
not depend on the compiler honouring contracts in imported templates
(``docs/toolchain/gcc16-feature-probe.md``).

.. code-block:: cpp

   PagedArena<FaceRecord> faces = frozen;      // shares every page
   faces.Set(17, record);                      // copies one page
   faces.Modify(42, [](FaceRecord& f) { f.Tolerance *= 2.0; });
   faces.ModifyPage(0, [](std::span<FaceRecord> page) {   // one copy-on-write step
       for (FaceRecord& f : page) { f.Tolerance = 1e-5; }
   });

A ``const`` arena can be read, and copied, from any thread; the mutating
operations belong to one thread, the builder's (ADR-0007).
``SharedPageCount(a, b)`` and ``IsPageShared(page)`` report what a snapshot
and its edit still share, which the tests and ``bench/ArenaBench.cpp`` use
to measure the copy cost. The page size is a tuning constant per table, fixed
by the Phase 0 measurement of ADR-0005 over 128, 256 and 1 024 entries;
nothing in the public API of the stores depends on it.

.. todo::

   Document the SoA stores, ``ShapeBuilder``, ``Freeze()`` and ``Compact()``,
   the Euler operators and the validity checker as they land.

Related design decisions
------------------------

- `ADR-0004: Topology storage and handles <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0004-topology-storage-and-handles.md>`_
- `ADR-0005: Immutability and structural sharing <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0005-immutability-and-sharing.md>`_

