// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <cstdio>
#include <format>
#include <string>

module Stereon.Core;

import :Log;

namespace Stereon {

LogSink StderrLogSink()
{
    return [](const LogRecord& record) {
        // One fputs per record, so lines from several threads do not
        // interleave within a line.
        const std::string line = std::format("{}: {} ({}:{})\n", record.Level, record.Message,
                                             record.Location.FileName, record.Location.Line);
        std::fputs(line.c_str(), stderr);
    };
}

} // namespace Stereon
