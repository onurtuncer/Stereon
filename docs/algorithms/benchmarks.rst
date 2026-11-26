.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Benchmark pairing and reproducibility
=====================================

Pair measurements by operation, variant and implementation. bench_report.py uses median CPU time over repetitions and reports the largest relative deviation from the median alongside ratios. An unmatched pair cannot support a comparative speed claim. bench_trend.py orders stored runs and compares consecutive results against the configured threshold. Cross-platform geometry comparison checks topology and coordinate tolerances; bit identity is recorded separately. The benchmark manual defines the local workloads, thresholds and release criteria.

For experimental design and variability, see :cite:p:`kalibera2013benchmark`.

Detailed treatment: :doc:`/performance`.

Implementation notes
--------------------

.. algorithm-sources::
   bench/ArenaBench.cpp
   bench/BatchBench.cpp
   bench/BooleanBench.cpp
   bench/BroadPhaseBench.cpp
   bench/GeometryBench.cpp
   bench/HotPathBench.cpp
   bench/IntersectBench.cpp
   bench/IntersectCases.hpp
   bench/Memory.cpp
   bench/Memory.hpp
   bench/RobustBench.cpp
   bench/SimdBench.cpp
   bench/SurfaceDerivativeBench.cpp
   bench/TopoBench.cpp
   bench/convergence/ProjectionCorpus.cpp
   bench/occt/OcctBench.cpp
   bench/occt/OcctIntersectBench.cpp
   bench/occt/OcctTopoBench.cpp
   tools/bench_report.py
   tools/bench_trend.py
   tools/check_release_benchmarks.py
   tools/compare_g2.py
