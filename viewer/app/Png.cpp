// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Png.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <span>
#include <string_view>
#include <vector>

namespace StereonView {

namespace {

/// The most bytes one stored deflate block holds.
constexpr std::size_t StoredBlock = 65535;

[[nodiscard]] std::array<std::uint32_t, 256> CrcTable()
{
    std::array<std::uint32_t, 256> table{};
    for (std::uint32_t i = 0; i < 256; ++i) {
        std::uint32_t c = i;
        for (int bit = 0; bit < 8; ++bit) {
            c = (c & 1U) != 0 ? 0xEDB88320U ^ (c >> 1) : c >> 1;
        }
        table[i] = c;
    }
    return table;
}

[[nodiscard]] std::uint32_t Crc(const std::span<const std::uint8_t> bytes)
{
    static const std::array<std::uint32_t, 256> table = CrcTable();
    std::uint32_t crc = 0xFFFFFFFFU;
    for (const std::uint8_t b : bytes) {
        crc = table[(crc ^ b) & 0xFFU] ^ (crc >> 8);
    }
    return crc ^ 0xFFFFFFFFU;
}

void PutU32(std::vector<std::uint8_t>& out, const std::uint32_t value)
{
    for (int shift = 24; shift >= 0; shift -= 8) {
        out.push_back(static_cast<std::uint8_t>((value >> shift) & 0xFFU));
    }
}

void PutChunk(std::vector<std::uint8_t>& out, const std::string_view tag,
              const std::span<const std::uint8_t> data)
{
    PutU32(out, static_cast<std::uint32_t>(data.size()));
    const std::size_t start = out.size();
    out.insert(out.end(), tag.begin(), tag.end());
    out.insert(out.end(), data.begin(), data.end());
    PutU32(out, Crc(std::span<const std::uint8_t>(out).subspan(start)));
}

/// The zlib stream of `raw` in stored (uncompressed) blocks.
[[nodiscard]] std::vector<std::uint8_t> Stored(const std::vector<std::uint8_t>& raw)
{
    std::vector<std::uint8_t> z{0x78, 0x01};
    std::uint32_t a = 1;
    std::uint32_t b = 0;
    for (std::size_t at = 0; at < raw.size() || at == 0; at += StoredBlock) {
        const std::size_t length = std::min(StoredBlock, raw.size() - at);
        const bool last = at + length >= raw.size();
        z.push_back(static_cast<std::uint8_t>(last ? 1 : 0));
        z.push_back(static_cast<std::uint8_t>(length & 0xFFU));
        z.push_back(static_cast<std::uint8_t>(length >> 8));
        z.push_back(static_cast<std::uint8_t>(~length & 0xFFU));
        z.push_back(static_cast<std::uint8_t>((~length >> 8) & 0xFFU));
        for (std::size_t k = at; k < at + length; ++k) {
            z.push_back(raw[k]);
            a = (a + raw[k]) % 65521U;
            b = (b + a) % 65521U;
        }
        if (raw.empty()) {
            break;
        }
    }
    PutU32(z, (b << 16) | a);
    return z;
}

} // namespace

std::vector<std::uint8_t> EncodePng(const int width, const int height,
                                    const std::span<const std::uint8_t> rgba)
{
    const auto w = static_cast<std::size_t>(width);
    const auto h = static_cast<std::size_t>(height);
    std::vector<std::uint8_t> raw;
    raw.reserve(h * (1 + 4 * w));
    for (std::size_t row = 0; row < h; ++row) {
        raw.push_back(0); // filter: none
        const auto first = rgba.begin() + static_cast<std::ptrdiff_t>(row * w * 4);
        raw.insert(raw.end(), first, first + static_cast<std::ptrdiff_t>(w * 4));
    }
    std::vector<std::uint8_t> png{0x89, 'P', 'N', 'G', '\r', '\n', 0x1A, '\n'};
    std::vector<std::uint8_t> header;
    PutU32(header, static_cast<std::uint32_t>(width));
    PutU32(header, static_cast<std::uint32_t>(height));
    header.insert(header.end(), {8, 6, 0, 0, 0}); // 8 bits, RGBA, deflate, no filter, no interlace
    PutChunk(png, "IHDR", header);
    PutChunk(png, "IDAT", Stored(raw));
    PutChunk(png, "IEND", {});
    return png;
}

bool WritePng(const std::filesystem::path& path, const int width, const int height,
              const std::span<const std::uint8_t> rgba)
{
    const std::vector<std::uint8_t> png = EncodePng(width, height, rgba);
    std::ofstream file(path, std::ios::binary | std::ios::trunc);
    file.write(reinterpret_cast<const char*>(png.data()), // NOLINT: byte I/O
               static_cast<std::streamsize>(png.size()));
    file.close();
    return !file.fail();
}

} // namespace StereonView
