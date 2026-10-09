# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# Libraries in layering order: each may depend only on those above it.
# tools/check_layering.py reads this list, so it is the single source of the
# layering rule (CONTRIBUTING.md, "Library layering").
add_subdirectory(core)
add_subdirectory(robust)
add_subdirectory(curve)
add_subdirectory(surface)
add_subdirectory(intersect)
add_subdirectory(topo)
add_subdirectory(build)
add_subdirectory(boolean)
add_subdirectory(blend)
add_subdirectory(mesh)
add_subdirectory(io)
