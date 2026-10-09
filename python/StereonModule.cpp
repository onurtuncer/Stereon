// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <nanobind/nanobind.h>

#include <format>

import Stereon.Core;

NB_MODULE(stereon, module)
{
    module.doc() = "Stereon B-rep geometric modeling kernel";
    const auto version = Stereon::LibraryVersion;
    const auto text = std::format("{}.{}.{}", version.Major, version.Minor, version.Patch);
    module.attr("__version__") = nanobind::str(text.c_str(), text.size());
}
