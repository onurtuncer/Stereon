// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The serialiser bridge of Records.cppm. This unit deliberately parses no
// contract and so does not include <Stereon/Core/Contracts.hpp>, unlike
// the library's other implementation units: see Records.cppm. It calls no
// inline module function that carries one.

module;

#include <cstddef>
#include <span>
#include <string>
#include <vector>

module Stereon.Io;

import :Records;

namespace Stereon::Detail {

template <typename T>
std::vector<std::byte> RecordBytes(const T& value)
{
    return Serialized(value);
}

template <typename T>
Result<T> RecordValue(const std::span<const std::byte> bytes)
{
    return Deserialized<T>(bytes);
}

template <typename T>
std::string RecordJson(const T& value)
{
    return ToJson(value);
}

// Every built-in type, in space and in the plane, every surface, and the
// tails of the composite records: the trims and the sweeps.
#define STEREON_ARCHIVE_RECORD(T)                                                                  \
    template std::vector<std::byte> RecordBytes<T>(const T&);                                      \
    template Result<T> RecordValue<T>(std::span<const std::byte>);                                 \
    template std::string RecordJson<T>(const T&);

STEREON_ARCHIVE_RECORD(Line)
STEREON_ARCHIVE_RECORD(Circle)
STEREON_ARCHIVE_RECORD(Ellipse)
STEREON_ARCHIVE_RECORD(Hyperbola)
STEREON_ARCHIVE_RECORD(Parabola)
STEREON_ARCHIVE_RECORD(BSplineCurve)
STEREON_ARCHIVE_RECORD(Line2d)
STEREON_ARCHIVE_RECORD(Circle2d)
STEREON_ARCHIVE_RECORD(Ellipse2d)
STEREON_ARCHIVE_RECORD(Hyperbola2d)
STEREON_ARCHIVE_RECORD(Parabola2d)
STEREON_ARCHIVE_RECORD(BSplineCurve2d)
STEREON_ARCHIVE_RECORD(Plane)
STEREON_ARCHIVE_RECORD(Cylinder)
STEREON_ARCHIVE_RECORD(Cone)
STEREON_ARCHIVE_RECORD(Sphere)
STEREON_ARCHIVE_RECORD(Torus)
STEREON_ARCHIVE_RECORD(BSplineSurface)
STEREON_ARCHIVE_RECORD(TrimRecord)
STEREON_ARCHIVE_RECORD(SurfaceTrimRecord)
STEREON_ARCHIVE_RECORD(ExtrusionRecord)
STEREON_ARCHIVE_RECORD(RevolutionRecord)

#undef STEREON_ARCHIVE_RECORD

} // namespace Stereon::Detail
