.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Architecture overview
=====================

The library layers under src/, the dependency rules between them, the Context object, and how the design decisions fit together.

The Context
-----------

Every operation takes a ``const Context&`` besides its inputs, and a result
is a function of the two alone: the kernel has no global or thread-local
state (ADR-0007). The ``Context`` carries the tolerance settings of ADR-0001
(``GetResolution``, ``GetAngularResolution``, ``GetToleranceCeiling``, the
±500 m ``GetSizeBox`` and the policy for tolerance growth), the scheduler
kind and thread count of ADR-0007, the validity level checked after every
operation (ADR-0018), the compaction threshold of ADR-0004, and the two
diagnostic sinks of ADR-0023: a ``LogSink`` that receives formatted log
records above a level threshold, and a ``JournalSink`` that receives the
entry an operation records before it runs and answers with the ``JournalId``
a ``KernelError`` carries for replay.

A ``Context`` is a value. Its settings are fixed at construction and changed
by making a modified copy, ``ctx.WithResolution(1e-4)``, so one context can be
shared by any number of operations on any number of threads. The settings are
the programmer's: a value outside its documented domain is a contract
violation, not a ``KernelError``. Without sinks a context is silent;
``StderrLogSink()`` and ``MemoryJournal`` are the ready-made sinks for tests
and tools, and an application attaches its own logging framework behind a
``LogSink``.

.. todo::

   Write the rest of this chapter: the library layering diagram and the data
   flow. It is a placeholder for the documentation effort.

Related design decisions
------------------------

- `ADR-0003: Geometry polymorphism via external polymorphism <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0003-geometry-polymorphism.md>`_
- `ADR-0007: Threading and determinism <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0007-threading-and-determinism.md>`_
- `ADR-0021: Shape derivatives as an operation contract <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0021-shape-derivatives.md>`_
- `ADR-0023: Logging and journal sinks <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0023-logging-and-journal-sinks.md>`_

