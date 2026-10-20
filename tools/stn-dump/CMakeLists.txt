# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# stn-dump (ADR-0012, ADR-0035): the JSON view of any .stn file, dumps
# included. Its tests are in tests/CMakeLists.txt.
add_executable(stn-dump StnDump.cpp)
target_link_libraries(stn-dump PRIVATE stereon::io stereon_options)
stereon_anchor_contracts(stn-dump)
