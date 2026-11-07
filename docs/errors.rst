.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Errors and contracts
====================

A kernel meets two different kinds of failure. Data from outside it, such as
a user's shape, an imported file or an operation's parameters, can be wrong
in ways the kernel must report and survive. The kernel's own code can be
wrong too: a broken internal assumption is a bug, and continuing would make
it worse. Stereon handles the first with **results**, values that carry
either an answer or a structured error, and the second with **contracts**,
checks the compiler places on functions (ADR-0006, ADR-0016). This chapter
explains both and the boundary between them.

Results and errors
------------------

Every operation that can fail returns ``Stereon::Result<T>``, an alias for
``std::expected<T, KernelError>``. The kernel does not throw; only
``std::bad_alloc`` may still propagate. An operation returns either its
result, which for a shape-producing operation has passed the validity
checker (:doc:`topology`), or a ``KernelError``. A shape that is silently
broken is never returned.

A ``KernelError`` carries:

- an ``ErrorCode``, the reason a program can act on;
- a message for people, naming the operation, the entities and the values
  involved;
- the **tolerances at stake**, in millimetres, so that a caller can retry
  with a raised ceiling or report the gap (:doc:`tolerances`);
- the **journal entry** of the failed operation when the ``Context`` has a
  journal sink, from which the failure can be replayed (ADR-0023,
  :doc:`debugging`).

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - ``ErrorCode``
     - Meaning
   * - ``InvalidInput``
     - Input data failed validation: a non-finite coordinate, a knot vector
       out of order, a handle of another shape, a tolerance below the
       resolution.
   * - ``InvalidResult``
     - The result failed its validity post-check (ADR-0018). Within the
       kernel's supported domain this is a kernel bug.
   * - ``ToleranceExceeded``
     - A local tolerance would grow above the ceiling, or a fit cannot meet
       its budget (ADR-0001).
   * - ``IntersectionFailed``
     - An intersection could not be computed completely, for example a root
       that could not be isolated. It is never reported as an empty
       intersection (ADR-0029).
   * - ``CertificationFailed``
     - A deviation bound could not be certified (ADR-0009).
   * - ``DegenerateResult``
     - The result would be degenerate, such as a zero-length edge or a face
       of no area.
   * - ``NotSupported``
     - The input is valid but outside what the kernel handles: a user
       geometry type without the capabilities an algorithm needs, or a
       periodic curve beyond its domain.
   * - ``NotImplemented``
     - The case is planned but not yet implemented, such as an intersection
       pair without a handler. The message names it.
   * - ``ParseError``
     - A file or a stream is malformed (:doc:`data_exchange`).
   * - ``Unresolved``
     - A question about an established state, such as whether a point lies
       inside a solid, could not be settled within the work limits
       (ADR-0029).

The numbers of the codes are part of the C ABI and never change; new codes
are appended (ADR-0020).

Propagating errors
~~~~~~~~~~~~~~~~~~

``STEREON_TRY`` and ``STEREON_TRY_VOID``, from ``<Stereon/Core/Try.hpp>``,
return an error to the caller and otherwise unwrap the value. They are the
``?`` operator C++ lacks:

.. code-block:: cpp

   #include <Stereon/Core/Try.hpp>

   #include <numbers>

   import Stereon.Core;

   using namespace Stereon;

   Result<double> ParseRadius(double value)
   {
       if (!(value > 0.0)) {
           return Fail(ErrorCode::InvalidInput, "radius must be positive");
       }
       return value;
   }

   Result<double> Area(double value)
   {
       STEREON_TRY(const double radius, ParseRadius(value));
       return std::numbers::pi * radius * radius;
   }

``Fail(code, message)`` builds the error. ``WithTolerance(value)`` and
``WithJournalId(id)`` add to it.

Warnings
~~~~~~~~

An operation that succeeds may still have something to report: a tolerance
that grew, or a cached curve it dropped. These are not errors. They arrive
as structured diagnostics in the operation's ``OperationResult``, each with
a ``DiagnosticCode`` and the entities involved, and are reported even with
logging off (:doc:`topology`, :doc:`tolerances`).

Contracts
---------

A **contract** states what a function requires of its caller (a
*precondition*), what it guarantees on return (a *postcondition*), or what
must hold at a point inside it (an *assertion*). C++26 makes these part of
the language (P2900), and the compiler checks them according to the
build's **evaluation semantic**. Stereon writes them with the macros of
``<Stereon/Core/Contracts.hpp>``:

.. code-block:: cpp

   double Eval(const Line& line, double t)
       STEREON_PRE(IsFinite(t))
       STEREON_HOT_PRE(InDomain(line, t))
       STEREON_POST(r : IsFinite(r));

Contracts are for **programming errors only**: a documented precondition
broken by the caller, a handle used after its shape died, an internal
invariant that does not hold. A violated contract never becomes a
``KernelError``. Anything that can come from outside the kernel, such as a
user's shape, an operation's parameters or a file's contents, is validated
and rejected with a ``KernelError`` and never guarded by a contract
(ADR-0006, ADR-0016). The checked factories show the split. ``MakeCircle``
validates a radius from outside and returns ``InvalidInput``, while
``Eval`` on a circle assumes a valid one and checks only its own
arguments.

Two kinds of contract
~~~~~~~~~~~~~~~~~~~~~

**Boundary contracts**, ``STEREON_PRE``, ``STEREON_POST`` and ``STEREON_ASSERT``.
   Cheap checks on public entry points, compiled into every build.

**Hot-path contracts**, ``STEREON_HOT_PRE``, ``STEREON_HOT_POST`` and ``STEREON_HOT_ASSERT``.
   Checks inside inner loops and evaluators, compiled only when
   ``STEREON_HOT_CONTRACTS`` is on. Release builds turn them off.

The evaluation semantic applies to a whole translation unit, not to one
contract, so this split is how an optimised build can enforce its boundary
checks and drop its hot-path ones.

Evaluation semantics
~~~~~~~~~~~~~~~~~~~~

``STEREON_CONTRACT_SEMANTIC`` selects one of the four semantics of P2900:

.. list-table::
   :header-rows: 1
   :widths: 18 52 30

   * - Semantic
     - On a violation
     - Used by
   * - ``enforce``
     - Calls the violation handler, which reports the contract and its
       location, then terminates.
     - ``gcc16-debug``, ``gcc16-release`` and the sanitizer presets
   * - ``quick_enforce``
     - Terminates at once, without a report.
     - none (available)
   * - ``observe``
     - Calls the violation handler and continues. The run fails afterwards.
     - ``gcc16-observe``, the nightly corpus run
   * - ``ignore``
     - Does not evaluate the check.
     - ``bench``, so that timings measure the algorithms

The ``observe`` run is how the kernel learns about violations it would
otherwise stop at. It runs the whole test suite and corpus with every
violation reported, so one bug cannot hide the next (:doc:`testing`).

Checked handles
~~~~~~~~~~~~~~~

With ``STEREON_CHECKED_HANDLES``, as in ``gcc16-debug``, every handle to a
topological entity carries the identity of its store and the generation of
its slot, and every access checks both. Using a handle of another shape, or
one whose entity was deleted, then fails a contract instead of reading the
wrong entity (ADR-0004, :doc:`topology`). Release builds use plain 32-bit
indices.

Compiler support
~~~~~~~~~~~~~~~~

GCC 16 implements P2900 contracts. Under Clang and MSVC, which are not built
in CI, the macros expand to nothing, so code using them still
compiles. GCC 16.2 has bugs in contracts that cross module boundaries,
which the build works around (``stereon_anchor_contracts``, :doc:`building`).
``docs/toolchain/gcc16-feature-probe.md`` records them, and the workarounds
are removed as fixed compilers arrive.

Related design decisions
------------------------

- `ADR-0006: Error handling <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0006-error-handling.md>`_
- `ADR-0016: C++26 adoption and contracts <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0016-cpp26-adoption-and-contracts.md>`_
- `ADR-0020: C++, C ABI and Python bindings <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0020-cpp-c-abi-and-python-bindings.md>`_
- `ADR-0023: Logging and journal sinks <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0023-logging-and-journal-sinks.md>`_
- `ADR-0029: Intersection and classification outcomes <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0029-intersection-and-classification-outcomes.md>`_
