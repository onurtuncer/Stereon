# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# stn-view (ADR-0035): GLFW for the window, OpenGL 3.3 core through the
# loader in Gl.cpp, Dear ImGui for the panels. The GL code files import no
# Stereon module; Main.cpp does.
add_executable(stn-view Main.cpp Camera.cpp Gl.cpp Png.cpp Renderer.cpp)
target_link_libraries(stn-view PRIVATE stereon::viewer_scene stereon_imgui stereon_options)
stereon_anchor_contracts(stn-view)