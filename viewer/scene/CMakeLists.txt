# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# The viewer's scene: dump files to display primitives, headless.
add_library(stereon_viewer_scene)
add_library(stereon::viewer_scene ALIAS stereon_viewer_scene)
target_sources(stereon_viewer_scene
    PUBLIC FILE_SET CXX_MODULES FILES Scene.cppm
    PRIVATE Primitives.cpp Scene.cpp)
target_link_libraries(stereon_viewer_scene PUBLIC stereon::io PRIVATE stereon_options)
target_compile_features(stereon_viewer_scene PUBLIC cxx_std_${STEREON_CXX_STANDARD})