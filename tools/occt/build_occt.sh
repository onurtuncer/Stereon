#!/usr/bin/env bash
# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# Builds the pinned OpenCASCADE that Stereon's differential tests and OCCT
# benchmarks compare against (docs/benchmarks.md, TODO section 9).
#
#   tools/occt/build_occt.sh <install directory> [<work directory>]
#
# Everything that decides the build is in this file: the CI cache key is its
# hash, so changing the tag or an option rebuilds OCCT. Only the modules the
# comparisons need are built: FoundationClasses, ModelingData and
# ModelingAlgorithms (booleans, primitives, meshing, healing, properties).
# Visualization, DataExchange, ApplicationFramework and Draw are off, which
# also removes the Tcl/Tk, FreeType, X11 and OpenGL dependencies.
#
# OCCT is LGPL-2.1 with the OCCT exception. It is linked only into Stereon's
# OCCT test and benchmark executables, which CI builds and never
# distributes; no Stereon library depends on it (ADR-0008).

set -euo pipefail

OCCT_TAG="V8_0_1"
OCCT_COMMIT="b8f597c677811d1f9f4d8a97f5ae2825c0353a42"

install_dir="${1:?usage: build_occt.sh <install directory> [<work directory>]}"
work_dir="${2:-$(mktemp -d)}"

git clone --quiet --depth 1 --branch "${OCCT_TAG}" https://github.com/Open-Cascade-SAS/OCCT.git "${work_dir}/src"
actual="$(git -C "${work_dir}/src" rev-parse HEAD)"
if [ "${actual}" != "${OCCT_COMMIT}" ]; then
    echo "build_occt: tag ${OCCT_TAG} is at ${actual}, expected ${OCCT_COMMIT}" >&2
    exit 1
fi

cmake -S "${work_dir}/src" -B "${work_dir}/build" -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DINSTALL_DIR="${install_dir}" \
    -DBUILD_LIBRARY_TYPE=Shared \
    -DBUILD_CPP_STANDARD=C++17 \
    -DBUILD_MODULE_FoundationClasses=ON \
    -DBUILD_MODULE_ModelingData=ON \
    -DBUILD_MODULE_ModelingAlgorithms=ON \
    -DBUILD_MODULE_Visualization=OFF \
    -DBUILD_MODULE_ApplicationFramework=OFF \
    -DBUILD_MODULE_DataExchange=OFF \
    -DBUILD_MODULE_Draw=OFF \
    -DUSE_TK=OFF \
    -DUSE_FREETYPE=OFF \
    -DUSE_XLIB=OFF \
    -DUSE_OPENGL=OFF \
    -DUSE_TBB=OFF \
    -DBUILD_GTEST=OFF \
    -DBUILD_DOC_Overview=OFF
cmake --build "${work_dir}/build"
cmake --install "${work_dir}/build"

# Record what was built next to it, for benchmark reports.
{
    echo "tag=${OCCT_TAG}"
    echo "commit=${OCCT_COMMIT}"
    echo "compiler=$(c++ --version | head -1)"
    echo "build_type=Release"
} > "${install_dir}/stereon-occt-build.txt"
echo "build_occt: OCCT ${OCCT_TAG} installed in ${install_dir}"
