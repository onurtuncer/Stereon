.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Error handling and contracts
============================

KernelError and std::expected results, contract checks, and how build presets choose contract semantics.

Results and errors
------------------

Every operation that can fail returns ``Stereon::Result<T>``, an alias for
``std::expected<T, KernelError>``. The kernel does not throw; only
``std::bad_alloc`` may still propagate. A ``KernelError`` holds an
``ErrorCode``, a message and the tolerance values at stake.

``STEREON_TRY`` and ``STEREON_TRY_VOID``, from ``<Stereon/Core/Try.hpp>``,
return an error to the caller and otherwise unwrap the value:

.. code-block:: cpp

   #include <Stereon/Core/Try.hpp>

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
       return 3.141592653589793 * radius * radius;
   }

Invalid data from outside the kernel is a ``KernelError``. A broken
obligation of the kernel's own code, such as a handle from another shape,
is a contract violation instead (ADR-0016).

.. todo::

   Document contract checks and how build presets choose contract semantics.

Related design decisions
------------------------

- `ADR-0006: Error handling <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0006-error-handling.md>`_

