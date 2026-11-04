.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Python bindings
===============

Stereon's Python module, ``stereon``, binds the C++ kernel directly with
nanobind (ADR-0020). Python is for interactive work, scripting and tests
written in a notebook. Every algorithm is implemented once, in C++, and the
bindings delegate to it.

What exists today
-----------------

The module builds with the kernel, is tested in CI on every change, and
reports the kernel's version:

.. code-block:: python

   >>> import stereon
   >>> stereon.__version__
   '0.1.0'

That is all it exports so far. The bindings of the geometry, topology and
operation APIs are planned, following the rules below; until they exist,
Stereon is used from C++.

Building and importing
----------------------

Add ``-DSTEREON_BUILD_PYTHON=ON`` to any preset; Python 3.11 or newer and its
development headers are required (:doc:`building`):

.. code-block:: bash

   cmake --preset gcc16-debug -DSTEREON_BUILD_PYTHON=ON
   cmake --build --preset gcc16-debug
   PYTHONPATH=build/gcc16-debug/python python -c "import stereon; print(stereon.__version__)"

The extension is built in the build tree as ``stereon``. ``ctest`` runs the
test ``python: import stereon``, which imports it and checks its version
against the project's. There is no wheel or installed package yet. ADR-0019
sets out how the self-contained package will be built and distributed.

How the bindings will look
--------------------------

ADR-0020 fixes the rules the bindings follow as they grow:

- **Direct bindings.** Python calls the C++ API through nanobind, not
  through the C ABI. Experimental APIs may be bound in Python before they
  enter the stable C ABI.
- **Python names.** Functions are ``snake_case`` and classes ``PascalCase``,
  through an explicit, reviewed mapping rather than automatic renaming.
- **Results.** A shape-producing function returns an ``OperationResult``
  that exposes the shape, its history, its tolerance changes and its
  warnings. No diagnostic is dropped on the way to Python.
- **Errors.** A ``KernelError`` becomes a Python exception carrying its code,
  message and structured diagnostics, and allocation failure becomes
  ``MemoryError`` (:doc:`errors`).
- **Batches.** Evaluation and mesh data cross the boundary as contiguous
  numerical buffers, so that no Python call is made per point or per
  triangle.
- **Derivatives.** Sensitivity products take and return contiguous numeric
  buffers, without exposing AD scalar types or tapes (ADR-0021,
  :doc:`sensitivity`).

Implementation
--------------

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Code
     - Role
   * - ``python/StereonModule.cpp``
     - The nanobind module, ``NB_MODULE(stereon, module)``.
   * - ``python/CMakeLists.txt``
     - The extension target ``stereon_python``, output name ``stereon``,
       and its import test.
   * - ``.github/workflows/python.yml``
     - Builds and tests the module on every change.

Related design decisions
------------------------

- `ADR-0019: Library linking, packaging and ABI <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0019-library-linking-packaging-and-abi.md>`_
- `ADR-0020: C++, C ABI and Python bindings <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0020-cpp-c-abi-and-python-bindings.md>`_
- `ADR-0021: Shape derivatives <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0021-shape-derivatives.md>`_
