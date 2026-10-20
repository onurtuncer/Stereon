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

Mutation testing
----------------

Coverage shows that a test executes code, not that it would notice the
code being wrong. ADR-0032 asks the opposite question: plant one plausible
fault in the implementation, run the unchanged tests, and see whether they
fail. ``tools/mutate.py`` does this on the reference compiler.

A catalogue, such as ``tools/mutation/pilot.json``, lists the mutants. Each
is one conceptual fault: an exact, unique piece of production source and
its replacement, the operator (a relational or arithmetic change, or a
domain fault such as accepting a stale handle), the ADR contract it
challenges, the test target to rebuild and the CTest selection that should
detect it. A **critical** mutant also names its diagnostic: the test whose
failure counts as detection.

The harness checks out the current commit into a disposable git worktree
under ``build/mutation``, configures it with the reference GCC preset and
builds and runs every selected test set unmodified; the baseline must pass
before any mutant counts. It then applies each mutant alone, rebuilds
incrementally, runs the selection and restores the file. Tests and oracles
are never mutated, and the worktree is removed at the end, so no mutant can
reach a branch or an artefact.

Each mutant is classified as ADR-0032 defines:

.. list-table::
   :header-rows: 1
   :widths: 18 82

   * - Outcome
     - Meaning
   * - Killed
     - The tests fail, and fail again when rerun.
   * - Survived
     - The tests pass, although they reach the mutated code: the harness
       plants an execution trap (``__builtin_trap()``) at the same site and
       sees the tests fail on it. A verification gap.
   * - Uncovered
     - The tests pass, and pass with the trap too: they never reach the
       site. An adequacy gap, not a kill.
   * - Equivalent
     - The catalogue records an expert rationale that no supported
       behaviour can tell the mutant apart, and it survived. Should it be
       killed, the claim is reported as refuted.
   * - Invalid
     - The mutant does not build, or its source text is not unique.
   * - Inconclusive
     - A test timed out, or the result was not reproducible.

The **mutation score** is killed / (killed + survived + uncovered); the
equivalent, invalid and inconclusive mutants are excluded and listed with
their reasons. The harness writes ``report.json`` and ``report.md`` and
exits non-zero if the baseline fails or a critical mutant is not killed
with its declared diagnostic.

.. code-block:: bash

   python3 tools/mutate.py --catalogue tools/mutation/pilot.json
   python3 tools/mutate.py --catalogue tools/mutation/pilot.json --only PRED-SIGN

The ``Mutation`` workflow runs the pilot catalogue weekly, as accepted
ADR-0032 requires, and on request, off the pull-request path, and publishes
the report. As ADR-0032 asks, it
reports rather than gates on a score; it fails only when the baseline
fails or a critical fault goes undetected. The pilot's results, and what
they changed in the tests, are in ``docs/planning/mutation-pilot.md``.

Fuzzing
-------

Bytes from outside the kernel must give a ``KernelError``, never a crash or a
contract violation (ADR-0006, ADR-0016). The fuzz targets in ``tests/fuzz/``
feed such bytes to the kernel's entry points:

- ``StnReadFuzz``: the ``.stn`` reader and its JSON view. An archive the
  reader accepts must also write and read back.
- ``BSplineFuzz``: the checked B-spline factory, then evaluation,
  derivatives, bounds and projection of whatever it accepts, through the
  erased wrapper whose parameter checks are compiled.

Each target exports ``LLVMFuzzerTestOneInput``, the libFuzzer entry point.
GCC has no libFuzzer, so ``tests/fuzz/StandaloneMain.cpp`` drives the targets
instead. It runs a target over its seed corpus, then over mutations of it:
flipped bits, boundary bytes, dropped, duplicated and inserted slices, and
overwritten eight-byte words. The mutations come from a seeded generator, so
a run can be repeated exactly:

.. code-block:: sh

   build/gcc16-asan/tests/stereon_fuzz_StnReadFuzz --runs 200000 --seed 20261020 tests/io/data

The driver counts contract violations and fails when there is one, writing
the input to ``fuzz-finding-<n>.bin``. That is what makes the ``observe``
build useful: violations are reported and the run continues, so one run
finds them all.

CTest runs every target on its corpus and 2000 mutations in every
configuration. The nightly workflow runs 200 000 mutations each under
AddressSanitizer and UBSan, with the day's date as the seed. It also runs the
whole suite in the ``gcc16-observe`` configuration.

A Clang build can link the targets with libFuzzer instead
(``-DSTEREON_LIBFUZZER=ON``), once Clang builds the kernel.
Determinism
-----------

ADR-0007 promises bit-identical results for any scheduler and thread count
on one executable and one machine (G1), and says every parallel algorithm
needs a determinism test. Parallel work runs through
``Detail::ParallelFor`` (``Stereon/Core/Parallel.hpp``). It cuts the index
range into chunks of fixed size whatever the thread count, runs them on the
calling thread for the serial scheduler, and otherwise on a pool of exactly
the ``Context``'s thread count (stdexec's ``static_thread_pool``). Each
chunk writes only its own indices' slots, so nothing depends on which
thread ran what or when. The operations built on it so far are
``ProjectAll`` onto B-spline curves and surfaces.

They are checked at three levels:

- Every pull request runs ``[determinism]`` in the test executables: the
  serial scheduler, 1, 2, 8 and 32 threads and the hardware's count, each
  twice, compared bit for bit with one another and with ``Project`` called
  point by point. It also runs the CTest case ``determinism: serial and 1,
  2, 8 and 32 threads give the same bytes``. The TSan job runs both, so a
  data race shows up there too.
- The ``Nightly`` workflow runs ``stereon_determinism`` once per
  configuration, each in its own process and twice, on 200 000 points. It
  compares the files byte for byte, puts their SHA-256 digests in the run's
  summary, and fails if any two differ. It also runs the full test suite on
  the debug and release presets.
- A new parallel operation adds its case to
  ``tests/support/DeterminismCases.hpp``, so that both checks cover it.

The cross-platform comparison (G2, Linux against Windows within tolerance)
is not run yet.

.. todo::

   Write the rest of this chapter: unit and property tests, fuzzing,
   regression cases and the corpus.

