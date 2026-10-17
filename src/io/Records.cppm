// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <cstddef>
#include <span>
#include <string>
#include <vector>

export module Stereon.Io:Records;

import Stereon.Core;

/// The archive's bridge to the reflection serialiser: one value of a
/// built-in geometry type to and from its record bytes and JSON.
///
/// These are declared here and defined, with an explicit instantiation for
/// every built-in type, in Records.cpp, a unit that parses no contract. A
/// unit that parses one (as GeometryArchive.cpp does, for the contract
/// anchor) cannot instantiate the serialiser's templates: GCC 16.2's
/// built-in std::contracts::contract_violation then conflicts with the one
/// <meta> brings in (docs/toolchain/gcc16-feature-probe.md).

export namespace Stereon::Detail {

/// The trim of a trimmed-curve record (format version 2), after its basis
/// record: the range and the sense.
struct TrimRecord {
    double Lower = 0.0;
    double Upper = 0.0;
    bool Reversed = false;
};

template <typename T>
[[nodiscard]] std::vector<std::byte> RecordBytes(const T& value);

template <typename T>
[[nodiscard]] Result<T> RecordValue(std::span<const std::byte> bytes);

template <typename T>
[[nodiscard]] std::string RecordJson(const T& value);

} // namespace Stereon::Detail
