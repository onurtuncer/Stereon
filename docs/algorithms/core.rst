.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Scalar geometry, spatial search and execution
=============================================

Frames, vectors and small matrices provide coordinates for the higher-level solvers. Interval arithmetic encloses their results by rounding the endpoints outwards. In the box-pair search, internal nodes enclose their descendants. The builder partitions entries around a median centre; disjoint node boxes eliminate every pair below them. Leaf comparisons supply the actual hits, which are emitted in deterministic input order. A box hit is a candidate for a narrow-phase test. See also :doc:`/performance` for the cost of producing candidates and :doc:`/tolerances` for the distinct model-resolution and entity-tolerance policies.

The storage and execution units below implement Stereon's own paging, scheduling, configuration and diagnostics. Their contracts and implementation notes are the primary sources for those engineering choices.

Method references: :cite:p:`ieee754,moore2009interval`.

For numerical conditioning and roundoff, see :cite:p:`higham2002accuracy`.

Detailed treatment: :doc:`/architecture`.

Implementation notes
--------------------

.. algorithm-sources::
   src/core/Allowances.cppm
   src/core/Arena.cpp
   src/core/Arena.cppm
   src/core/Box.cppm
   src/core/BoxPairs.cpp
   src/core/BoxPairs.cppm
   src/core/Context.cppm
   src/core/Core.cppm
   src/core/Enum.cppm
   src/core/Error.cppm
   src/core/EvaluatorChecks.cpp
   src/core/EvaluatorChecks.cppm
   src/core/Frame.cppm
   src/core/GeometryChecks.cpp
   src/core/GeometryChecks.cppm
   src/core/Handle.cpp
   src/core/Handle.cppm
   src/core/Journal.cppm
   src/core/Log.cpp
   src/core/Log.cppm
   src/core/Matrix.cppm
   src/core/Parallel.cpp
   src/core/Scalar.cppm
   src/core/Serialize.cpp
   src/core/Serialize.cppm
   src/core/Simd.cppm
   src/core/SmallBuffer.cppm
   src/core/Trigonometry.cpp
   src/core/Trigonometry.cppm
   src/core/Vector.cppm
   src/core/Version.cppm
   src/core/include/Stereon/Core/Contracts.hpp
   src/core/include/Stereon/Core/Exec.hpp
   src/core/include/Stereon/Core/Parallel.hpp
   src/core/include/Stereon/Core/Try.hpp
   src/core/private/Profiling.hpp
