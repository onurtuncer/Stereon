.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Build integration and worked examples
=====================================

The build scripts select dependencies, compiler features and optional capabilities. Standalone examples compose public kernel calls and expose their results. Compiler reproductions isolate a toolchain behaviour in a small program. These artefacts specify Stereon workflows; their verification and reproduction instructions are in :doc:`/building`, :doc:`/examples` and :doc:`/testing`.

Detailed treatment: :doc:`/building`.

Implementation notes
--------------------

.. algorithm-sources::
   CMakeLists.txt
   bench/CMakeLists.txt
   cmake/StereonDependencies.cmake
   cmake/StereonOptions.cmake
   cmake/StereonProfiling.cmake
   docs/CMakeLists.txt
   docs/_ext/stereon_algorithms.py
   docs/_ext/stereon_api.py
   docs/conf.py
   docs/tikz/CMakeLists.txt
   docs/toolchain/gcc-bugs/01-hidden-friend-eq/m.cppm
   docs/toolchain/gcc-bugs/01-hidden-friend-eq/main.cpp
   docs/toolchain/gcc-bugs/01-hidden-friend-eq/repro.sh
   docs/toolchain/gcc-bugs/02-stdexec-fmodules-ice/main.cpp
   docs/toolchain/gcc-bugs/02-stdexec-fmodules-ice/repro.sh
   docs/toolchain/gcc-bugs/03-source-location-tu-local/m.cppm
   docs/toolchain/gcc-bugs/03-source-location-tu-local/main.cpp
   docs/toolchain/gcc-bugs/03-source-location-tu-local/repro.sh
   docs/toolchain/gcc-bugs/04-contract-not-constant/main.cpp
   docs/toolchain/gcc-bugs/04-contract-not-constant/repro.sh
   docs/toolchain/gcc-bugs/06-bad-file-data-hidden-friend/a.cppm
   docs/toolchain/gcc-bugs/06-bad-file-data-hidden-friend/b.cppm
   docs/toolchain/gcc-bugs/06-bad-file-data-hidden-friend/c.cppm
   docs/toolchain/gcc-bugs/06-bad-file-data-hidden-friend/main.cpp
   docs/toolchain/gcc-bugs/06-bad-file-data-hidden-friend/repro.sh
   docs/toolchain/gcc-bugs/07-scrambled-aggregate/m.cppm
   docs/toolchain/gcc-bugs/07-scrambled-aggregate/main.cpp
   docs/toolchain/gcc-bugs/07-scrambled-aggregate/repro.sh
   docs/toolchain/gcc-bugs/08-recursive-lazy-load/a.cppm
   docs/toolchain/gcc-bugs/08-recursive-lazy-load/b.cppm
   docs/toolchain/gcc-bugs/08-recursive-lazy-load/c-part.cppm
   docs/toolchain/gcc-bugs/08-recursive-lazy-load/c.cppm
   docs/toolchain/gcc-bugs/08-recursive-lazy-load/main.cpp
   docs/toolchain/gcc-bugs/08-recursive-lazy-load/repro.sh
   docs/toolchain/gcc-bugs/09-var-template-explicit-spec/a.cppm
   docs/toolchain/gcc-bugs/09-var-template-explicit-spec/main.cpp
   docs/toolchain/gcc-bugs/09-var-template-explicit-spec/repro.sh
   docs/toolchain/gcc-bugs/10-mingw-format-dup/a.cpp
   docs/toolchain/gcc-bugs/10-mingw-format-dup/m.cppm
   docs/toolchain/gcc-bugs/10-mingw-format-dup/main.cpp
   docs/toolchain/gcc-bugs/10-mingw-format-dup/repro.sh
   docs/toolchain/gcc-bugs/11-mingw-make-shared-dup/a.cpp
   docs/toolchain/gcc-bugs/11-mingw-make-shared-dup/m.cppm
   docs/toolchain/gcc-bugs/11-mingw-make-shared-dup/main.cpp
   docs/toolchain/gcc-bugs/11-mingw-make-shared-dup/repro.sh
   docs/toolchain/gcc-bugs/12-contract-violation-conflict/m.cppm
   docs/toolchain/gcc-bugs/12-contract-violation-conflict/main.cpp
   docs/toolchain/gcc-bugs/12-contract-violation-conflict/repro.sh
   docs/toolchain/gcc-bugs/15-finish-member-declaration-ice/a.cppm
   docs/toolchain/gcc-bugs/15-finish-member-declaration-ice/b.cppm
   docs/toolchain/gcc-bugs/15-finish-member-declaration-ice/main.cpp
   docs/toolchain/gcc-bugs/15-finish-member-declaration-ice/repro.sh
   docs/toolchain/gcc-bugs/16-tu-has-violation/m.cppm
   docs/toolchain/gcc-bugs/16-tu-has-violation/main.cpp
   docs/toolchain/gcc-bugs/16-tu-has-violation/repro.sh
   docs/toolchain/gcc-bugs/17-imported-template-contracts-ignored/m.cppm
   docs/toolchain/gcc-bugs/17-imported-template-contracts-ignored/main.cpp
   docs/toolchain/gcc-bugs/17-imported-template-contracts-ignored/repro.sh
   docs/toolchain/gcc-bugs/ce-lib.ps1
   docs/toolchain/gcc-bugs/ce.ps1
   docs/toolchain/gcc-bugs/run-all.sh
   examples/CMakeLists.txt
   examples/boolean.cpp
   examples/curves.cpp
   examples/geometry_file.cpp
   examples/intersection.cpp
   examples/steinmetz.cpp
   python/CMakeLists.txt
   src/CMakeLists.txt
   src/blend/CMakeLists.txt
   src/boolean/CMakeLists.txt
   src/build/CMakeLists.txt
   src/core/CMakeLists.txt
   src/curve/CMakeLists.txt
   src/intersect/CMakeLists.txt
   src/io/CMakeLists.txt
   src/mesh/CMakeLists.txt
   src/robust/CMakeLists.txt
   src/surface/CMakeLists.txt
   src/topo/CMakeLists.txt
   tests/CMakeLists.txt
   tools/occt/build_occt.sh
   tools/stn-dump/CMakeLists.txt
   viewer/CMakeLists.txt
   viewer/app/CMakeLists.txt
   viewer/scene/CMakeLists.txt
