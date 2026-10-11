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

.. todo::

   Document the SoA stores, ``ShapeBuilder``, ``Freeze()`` and ``Compact()``,
   the Euler operators and the validity checker as they land.

Related design decisions
------------------------

- `ADR-0004: Topology storage and handles <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0004-topology-storage-and-handles.md>`_
- `ADR-0005: Immutability and structural sharing <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0005-immutability-and-sharing.md>`_

