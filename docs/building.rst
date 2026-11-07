.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Building and installing
=======================

This chapter builds the kernel, its tests and its optional components from
source, explains the build presets and options, and shows how to use the
libraries from another project and how to build this manual.

Requirements
------------

Stereon is written in C++26 and compiled as C++20 modules. It relies on
static reflection and contracts, which, as of 2026, only GCC 16 implements
completely, so **GCC 16 is the reference compiler** (ADR-0016). CI builds
every change with GCC 16.2 on Linux, Windows and macOS; those jobs block a
merge. Clang and MSVC are not built in CI until they support the same
features.

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Tool
     - Version and role
   * - GCC
     - 16 or newer, the only supported compiler. The CI reference is 16.2.0.
   * - CMake
     - 3.30 or newer, for C++20 module scanning.
   * - Ninja
     - The generator every preset uses; CMake builds modules only with
       Ninja here.
   * - Python
     - 3.11 or newer, optional: for the Python module, the source checks and
       this manual.
   * - Git
     - CMake fetches the dependencies from their repositories.

On **Windows**, build in an MSYS2 UCRT64 shell with the packages
``mingw-w64-ucrt-x86_64-gcc``, ``mingw-w64-ucrt-x86_64-cmake`` and
``mingw-w64-ucrt-x86_64-ninja``. On **macOS**, use GCC 16 from Homebrew, not
Apple Clang.

The first build
---------------

.. code-block:: bash

   git clone https://github.com/onurtuncer/Stereon.git
   cd Stereon
   cmake --preset gcc16-release      # configure into build/gcc16-release
   cmake --build --preset gcc16-release
   ctest --preset gcc16-release

The first configuration fetches the dependencies, and the first build takes
some ten minutes on a laptop. Each preset configures into
``build/<preset>``, so several can live side by side.

Presets
-------

A preset fixes the build type, how contracts are evaluated and which checks
are compiled in. For development use ``gcc16-debug``. Before opening a pull
request, run ``gcc16-release`` as well: it compiles with the optimiser's
warnings, such as ``-Wnull-dereference``, which a debug build does not see.

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Preset
     - Purpose
   * - ``gcc16-debug``
     - Development. Contracts enforced, checked handles on: every access to
       a topological entity is checked against its store and generation
       (ADR-0004, ADR-0016). Shape-producing operations also run the
       Geometric level of the validity checker.
   * - ``gcc16-release``
     - Optimised. Contracts at API boundaries are enforced; hot-path
       contracts are compiled out.
   * - ``gcc16-profile``
     - ``gcc16-release`` with symbols and Tracy instrumentation, for
       diagnostic timings (:doc:`performance`, ``docs/profiling.md``).
   * - ``gcc16-observe``
     - The CI corpus run. Contract violations are logged instead of
       terminating, and the run fails afterwards (ADR-0016).
   * - ``gcc16-coverage``
     - ``gcc16-debug`` instrumented for gcov; CI posts the report on pull
       requests.
   * - ``gcc16-asan``
     - AddressSanitizer and UndefinedBehaviorSanitizer. Linux and macOS
       only.
   * - ``gcc16-tsan``
     - ThreadSanitizer, with hot-path contracts off. Linux and macOS only.
   * - ``gcc16-windows``
     - ``gcc16-release`` with the GCC runtime linked statically, so the
       binaries run outside the MSYS2 shell.
   * - ``bench``
     - The benchmark suite: release, contracts ignored so that timings
       measure the algorithms, no tests (:doc:`performance`).
   * - ``clang-debug``, ``msvc-debug``
     - Local builds with other compilers, not run in CI.

Every preset except ``bench`` has a test preset of the same name, so
``ctest --preset <name>`` runs that build's tests.

Options
-------

The presets set these cache variables; any of them can be added to a
configure command, as in ``cmake --preset gcc16-debug -DSTEREON_BUILD_PYTHON=ON``.

.. list-table::
   :header-rows: 1
   :widths: 34 12 54

   * - Option
     - Default
     - Effect
   * - ``STEREON_BUILD_TESTS``
     - ``ON``
     - Unit, property and fuzz tests (Catch2, RapidCheck).
   * - ``STEREON_BUILD_BENCH``
     - ``OFF``
     - The benchmark suite (Google Benchmark).
   * - ``STEREON_BUILD_TOOLS``
     - ``ON``
     - The command-line tools, such as ``stn-dump`` (:doc:`debugging`).
   * - ``STEREON_BUILD_VIEWER``
     - ``OFF``
     - The debug viewer ``stn-view`` (GLFW, OpenGL, Dear ImGui).
   * - ``STEREON_BUILD_PYTHON``
     - ``OFF``
     - The Python module ``stereon`` (nanobind; :doc:`python`).
   * - ``STEREON_BUILD_DOCS``
     - ``OFF``
     - This manual and its API reference (Sphinx, Doxygen).
   * - ``STEREON_WITH_OCCT``
     - ``OFF``
     - Differential tests against a pinned OpenCASCADE (:doc:`testing`).
   * - ``STEREON_WITH_CPPAD``
     - ``OFF``
     - Integration tests of the CppAD backend (ADR-0045; :doc:`sensitivity`).
   * - ``STEREON_ENABLE_PROFILING``
     - ``OFF``
     - Private Tracy instrumentation (ADR-0043).
   * - ``STEREON_CONTRACT_SEMANTIC``
     - ``enforce``
     - How a contract is evaluated: ``ignore``, ``observe``, ``enforce`` or
       ``quick_enforce`` (:doc:`errors`).
   * - ``STEREON_HOT_CONTRACTS``
     - ``ON``
     - Compile the hot-path contracts (``STEREON_HOT_PRE`` and friends).
   * - ``STEREON_CHECKED_HANDLES``
     - ``OFF``
     - Handles carry their store and generation, and every access is
       checked (ADR-0004).
   * - ``STEREON_SANITIZER``
     - empty
     - ``address`` or ``thread``.
   * - ``STEREON_COVERAGE``
     - ``OFF``
     - gcov instrumentation.
   * - ``STEREON_WARNINGS_AS_ERRORS``
     - ``OFF``
     - Treat warnings as errors, as CI does.

Every Stereon target compiles with ``-Wall -Wextra -Wpedantic -Wconversion
-Wsign-conversion -Wshadow -Wnull-dereference`` and more, and CI adds
``-Werror``. Third-party code never sees these flags.

Dependencies
------------

CMake fetches every dependency at a pinned version
(``cmake/StereonDependencies.cmake``); nothing needs to be installed first.
Only the first two are linked into the libraries.

.. list-table::
   :header-rows: 1
   :widths: 22 18 60

   * - Dependency
     - Version
     - Use
   * - stdexec
     - ``nvhpc-26.05``
     - ``std::execution``, until the standard library ships it (ADR-0007).
   * - xsimd
     - 14.3.0
     - ``std::simd``, until the standard library ships it.
   * - Catch2
     - 3.16.1
     - The test framework.
   * - RapidCheck
     - pinned commit
     - Property-based tests.
   * - Google Benchmark
     - 1.9.5
     - The benchmark suite.
   * - nanobind
     - 3.1.0
     - The Python module.
   * - CppAD
     - 20260000.0
     - The AD backend's integration tests only. It is EPL-2.0 and is never
       linked into a library (ADR-0045).
   * - Tracy
     - 0.12.2
     - Profiling, with ``STEREON_ENABLE_PROFILING``.
   * - GLFW, Dear ImGui
     - 3.4, 1.91.8
     - The debug viewer.

Running the tests
-----------------

``ctest --preset <name>`` runs every test of a build, about eight hundred.
Each library has its own test executable, ``build/<preset>/tests/stereon_<library>_tests``,
which takes Catch2's filters:

.. code-block:: bash

   build/gcc16-debug/tests/stereon_intersect_tests "[linearization]"
   build/gcc16-debug/tests/stereon_topo_tests "G03*"

Property tests run a hundred random cases each by default. ``RC_PARAMS``
changes that, and replays a failure from its seed:

.. code-block:: bash

   RC_PARAMS="max_success=1000 seed=3" build/gcc16-debug/tests/stereon_intersect_tests "[property]"

On Windows, running an executable directly from outside the MSYS2 shell
needs ``C:\msys64\ucrt64\bin`` on ``PATH`` for the GCC runtime; ``ctest``
handles this itself. :doc:`testing` describes what the suites test.

Source checks
-------------

CI runs these checks on every change; ``ctest`` runs the Python ones
whenever Python is installed. ``pip install -r tools/requirements.txt``
installs the pinned ``clang-format`` and ``lizard`` they use.

.. code-block:: bash

   git ls-files -z -- '*.cpp' '*.cppm' '*.hpp' '*.h' | xargs -0 clang-format --dry-run --Werror
   python tools/check_layering.py      # no library imports one above it
   python tools/check_epsilon.py       # no literal tolerances in topo, boolean, mesh
   python tools/check_literals.py      # named constants (ADR-0024)
   python tools/check_complexity.py    # function size and complexity (ADR-0022)
   python tools/check_capabilities.py  # every operation declares its derivatives (ADR-0021)

Optional components
-------------------

**Python.**
   Add ``-DSTEREON_BUILD_PYTHON=ON``. The module is built in the build
   tree, so it can be imported with ``import stereon`` once that directory
   is on ``PYTHONPATH`` (:doc:`python`).

**Benchmarks.**
   ``cmake --preset bench && cmake --build --preset bench`` builds
   ``build/bench/bench/stereon_bench``. Its ``--benchmark_filter`` selects
   rows, as in ``--benchmark_filter=B1.1``. The comparison against
   OpenCASCADE needs ``STEREON_WITH_OCCT`` as well (:doc:`performance`).

**OpenCASCADE.**
   ``tools/occt/build_occt.sh`` builds the pinned OpenCASCADE 8.0 that the
   differential tests and the reference benchmarks compare against.
   Configure with ``-DSTEREON_WITH_OCCT=ON -DCMAKE_PREFIX_PATH=<its install>``.

**The viewer and the tools.**
   ``stn-dump`` prints a native file, and ``stn-view`` displays one
   (:doc:`debugging`).

Using Stereon from another project
----------------------------------

Stereon has no install rules or package configuration yet, and building it
from outside its own tree is not tested. A project builds it from source
alongside its own code, the way the kernel's own tests and tools consume it.
Each library is a CMake target ``stereon::<library>`` and a module
``Stereon.<Library>``, which re-exports the libraries below it.

.. code-block:: cmake

   include(FetchContent)
   FetchContent_Declare(stereon
       GIT_REPOSITORY https://github.com/onurtuncer/Stereon.git
       GIT_TAG        v0.1.0)
   set(STEREON_BUILD_TESTS OFF CACHE BOOL "" FORCE)
   FetchContent_MakeAvailable(stereon)

   add_executable(app main.cpp)
   # stereon_options carries the C++26 flags (-freflection, -fcontracts and
   # the contract semantic) the modules were compiled with; an importer must
   # match them.
   target_link_libraries(app PRIVATE stereon::surface stereon_options)
   stereon_anchor_contracts(app)

.. code-block:: cpp

   import Stereon.Surface;   // curves, surfaces and everything below them

The project needs GCC 16, CMake 3.30 and the Ninja generator.
``stereon_anchor_contracts`` works around a GCC 16.2 bug: a unit that
imports a Stereon module without parsing a contract of its own fails to link
the module's inline functions (``docs/toolchain/gcc16-feature-probe.md``).

Building this manual
--------------------

The manual is Sphinx with the Read the Docs theme. Its API reference comes
from Doxygen through Breathe, and its figures are TikZ.

.. code-block:: bash

   pip install -r docs/requirements.txt
   sphinx-build -b html docs build/html

That builds every chapter. The API pages stay empty unless Doxygen 1.17 runs
first, which ``cmake -DSTEREON_BUILD_DOCS=ON`` followed by the target
``Stereon_Docs`` does, or which ``STEREON_RUN_DOXYGEN=1`` asks ``conf.py`` to
do. The diagrams of :doc:`architecture` need Graphviz's ``dot``.

Each figure is a standalone LaTeX file ``docs/tikz/<name>.tex``, rendered to
``docs/_static/tikz/<name>.png`` with ``pdflatex`` and ``pdftoppm -r 200``:

.. code-block:: bash

   for tex in docs/tikz/*.tex; do
       name=$(basename "$tex" .tex)
       pdflatex -interaction=nonstopmode -halt-on-error -output-directory build/tikz "$tex"
       pdftoppm -r 200 -png -singlefile "build/tikz/$name.pdf" "docs/_static/tikz/$name"
   done

The documentation workflow does the same on every push to ``main`` and
publishes the result. Rendered figures are not committed.

Related design decisions
------------------------

- `ADR-0004: Topology storage and handles <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0004-topology-storage-and-handles.md>`_
- `ADR-0007: Threading and determinism <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0007-threading-and-determinism.md>`_
- `ADR-0016: C++26 adoption and contracts <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0016-cpp26-adoption-and-contracts.md>`_
- `ADR-0045: CppAD AD backend <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0045-cppad-ad-backend.md>`_
