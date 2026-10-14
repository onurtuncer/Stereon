.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Testing and benchmarks
======================

Unit, property-based, fuzz and regression tests, the geometry corpus, and the OpenCASCADE benchmark plan in docs/benchmarks.md.

Differential tests against OpenCASCADE
--------------------------------------

Tests written alongside AI-assisted code do not prove it correct
(ADR-0015); an independent referee does. For geometry the referee is
OpenCASCADE. ``tools/occt/build_occt.sh`` builds a pinned OCCT, tag
``V8_0_1`` at commit ``b8f597c``, with the modelling modules only, and the
``OCCT`` workflow caches the installed tree under a key derived from that
script, so it is rebuilt only when the pin or an option changes.

Configuring with ``-DSTEREON_WITH_OCCT=ON`` and the OCCT install directory in
``CMAKE_PREFIX_PATH`` builds ``stereon_occt_tests``: the same line, circle,
plane and cylinder are built in both kernels and evaluated, with first and
second derivatives, over two periods of parameters, and must agree to a
relative 1e-13. The volume, area and face-count comparisons of
``docs/benchmarks.md`` join when Stereon can build bodies.

.. code-block:: bash

   tools/occt/build_occt.sh "$HOME/occt"            # once; most of an hour
   cmake --preset gcc16-release -DSTEREON_WITH_OCCT=ON -DCMAKE_PREFIX_PATH="$HOME/occt"
   cmake --build --preset gcc16-release --target stereon_occt_tests
   ctest --preset gcc16-release -R "^occt: "

OCCT is LGPL-2.1 with the OCCT exception. It is linked into the test and
benchmark executables only, which are built on request and never
distributed; no Stereon library depends on it (ADR-0008).

.. todo::

   Write the rest of this chapter: unit and property tests, fuzzing,
   regression cases and the corpus.

