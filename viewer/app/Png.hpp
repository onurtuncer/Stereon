// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// A PNG writer for stn-view's screenshots: RGBA, 8 bits per channel, the
// image data in uncompressed deflate blocks. Larger than a compressed PNG,
// but it needs no image library (ADR-0035 keeps the viewer's dependencies
// to GLFW and Dear ImGui), and every PNG reader reads it.

#pragma once

#include <cstdint>
#include <filesystem>
#include <span>
#include <vector>

namespace StereonView {

/// The PNG file's bytes for RGBA pixels, rows from the top.
[[nodiscard]] std::vector<std::uint8_t> EncodePng(int width, int height,
                                                  std::span<const std::uint8_t> rgba);

/// Writes the PNG; false if the file cannot be written.
[[nodiscard]] bool WritePng(const std::filesystem::path& path, int width, int height,
                            std::span<const std::uint8_t> rgba);

} // namespace StereonView
