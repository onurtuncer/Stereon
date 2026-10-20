// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// stn-dump: the JSON view of a .stn file (ADR-0012, ADR-0035), dumps
// included: the header, the policy and every chunk with its records, the
// debug chunk's primitives among them. With --pretty, one member per line
// and flat arrays (a point, a knot vector) kept on one, so that two
// versions of a file diff line by line:
//
//     stn-dump --pretty before.stn > before.json
//     stn-dump --pretty after.stn > after.json
//     diff before.json after.json
//
// Exit status 0 on success, 1 for a file that cannot be read or parsed (the
// reason on stderr), 2 for bad arguments.

#include <cstddef>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <string_view>
#include <vector>

import Stereon.Io;

namespace {

constexpr std::string_view Usage = "usage: stn-dump [--pretty] <file.stn>\n";

/// The index just past the string literal that starts at `i`.
std::size_t SkipString(const std::string_view json, std::size_t i)
{
    for (++i; i < json.size(); ++i) {
        if (json[i] == '\\') {
            ++i;
        } else if (json[i] == '"') {
            return i + 1;
        }
    }
    return json.size();
}

/// The index just past the array that starts at `begin`, if it holds no
/// array or object; `begin` otherwise.
std::size_t FlatArrayEnd(const std::string_view json, const std::size_t begin)
{
    for (std::size_t i = begin + 1; i < json.size();) {
        const char c = json[i];
        if (c == '"') {
            i = SkipString(json, i);
        } else if (c == '[' || c == '{') {
            return begin;
        } else if (c == ']') {
            return i + 1;
        } else {
            ++i;
        }
    }
    return begin;
}

void NewLine(std::string& out, const int depth)
{
    out += '\n';
    out.append(static_cast<std::size_t>(2 * depth), ' ');
}

/// The JSON indented by two spaces per level, flat arrays on one line.
std::string Pretty(const std::string_view json)
{
    std::string out;
    int depth = 0;
    for (std::size_t i = 0; i < json.size();) {
        const char c = json[i];
        if (c == '"') {
            const std::size_t end = SkipString(json, i);
            out.append(json.substr(i, end - i));
            i = end;
            continue;
        }
        if (c == '[') {
            if (const std::size_t end = FlatArrayEnd(json, i); end != i) {
                out.append(json.substr(i, end - i));
                i = end;
                continue;
            }
        }
        if (c == '{' || c == '[') {
            out += c;
            NewLine(out, ++depth);
        } else if (c == '}' || c == ']') {
            NewLine(out, --depth);
            out += c;
        } else if (c == ',') {
            out += c;
            NewLine(out, depth);
        } else {
            out += c;
            if (c == ':') {
                out += ' ';
            }
        }
        ++i;
    }
    return out + '\n';
}

} // namespace

int main(const int argc, char** argv)
{
    const std::vector<std::string_view> args(argv + 1, argv + argc);
    bool pretty = false;
    std::string_view path;
    for (const std::string_view arg : args) {
        if (arg == "--pretty") {
            pretty = true;
        } else if (arg == "--help" || arg == "-h") {
            std::cout << Usage;
            return 0;
        } else if (path.empty() && !arg.starts_with("--")) {
            path = arg;
        } else {
            std::cerr << Usage;
            return 2;
        }
    }
    if (path.empty()) {
        std::cerr << Usage;
        return 2;
    }
    std::ifstream file{std::string(path), std::ios::binary};
    if (!file) {
        std::cerr << "stn-dump: cannot read \"" << path << "\"\n";
        return 1;
    }
    const std::vector<char> chars{std::istreambuf_iterator<char>(file),
                                  std::istreambuf_iterator<char>()};
    std::vector<std::byte> bytes;
    bytes.reserve(chars.size());
    for (const char c : chars) {
        bytes.push_back(static_cast<std::byte>(static_cast<unsigned char>(c)));
    }
    const Stereon::Result<std::string> json = Stereon::GeometryToJson(bytes);
    if (!json) {
        std::cerr << "stn-dump: " << json.error().GetMessage() << '\n';
        return 1;
    }
    std::cout << (pretty ? Pretty(*json) : *json + '\n');
    return 0;
}
