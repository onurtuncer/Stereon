# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# The debug viewer (ADR-0035): an application over the kernel's public
# interfaces, never a part of it (ADR-0034). viewer/scene is headless and
# built with the tests, which exercise it without a window; viewer/app,
# stn-view, needs GLFW, OpenGL and Dear ImGui and is built only with
# STEREON_BUILD_VIEWER=ON. Nothing under src/ may refer to either
# (tools/check_dependencies.py).
add_subdirectory(scene)
if(STEREON_BUILD_VIEWER)
    add_subdirectory(app)
endif()