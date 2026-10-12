// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// Catch2 output for Stereon math types through their std::formatter, so a
/// failed REQUIRE prints "(1, 2, 3)" instead of "{?}". Include it after
/// `import Stereon.Core;`: the specialisations name the module's types.

#pragma once

#include <catch2/catch_tostring.hpp>

#include <cstddef>
#include <format>
#include <string>

namespace Catch {

template <Stereon::Scalar T>
    requires std::formattable<T, char>
struct StringMaker<Stereon::BasicVec2<T>> {
    static std::string convert(const Stereon::BasicVec2<T>& v) { return std::format("{}", v); }
};

template <Stereon::Scalar T>
    requires std::formattable<T, char>
struct StringMaker<Stereon::BasicVec3<T>> {
    static std::string convert(const Stereon::BasicVec3<T>& v) { return std::format("{}", v); }
};

template <Stereon::Scalar T>
    requires std::formattable<T, char>
struct StringMaker<Stereon::BasicPoint2<T>> {
    static std::string convert(const Stereon::BasicPoint2<T>& p) { return std::format("{}", p); }
};

template <Stereon::Scalar T>
    requires std::formattable<T, char>
struct StringMaker<Stereon::BasicPoint3<T>> {
    static std::string convert(const Stereon::BasicPoint3<T>& p) { return std::format("{}", p); }
};

template <Stereon::Scalar T, std::size_t N>
    requires std::formattable<T, char>
struct StringMaker<Stereon::BasicMat<T, N>> {
    static std::string convert(const Stereon::BasicMat<T, N>& m) { return std::format("{}", m); }
};

template <Stereon::Scalar T>
    requires std::formattable<T, char>
struct StringMaker<Stereon::BasicTransform<T>> {
    static std::string convert(const Stereon::BasicTransform<T>& t) { return std::format("{}", t); }
};

template <>
struct StringMaker<Stereon::Interval> {
    static std::string convert(const Stereon::Interval& i) { return std::format("{}", i); }
};

template <Stereon::AnyBox B>
struct StringMaker<B> {
    static std::string convert(const B& box) { return std::format("{}", box); }
};

template <Stereon::EntityKind K>
struct StringMaker<Stereon::Id<K>> {
    static std::string convert(const Stereon::Id<K>& h) { return std::format("{}", h); }
};

template <Stereon::EntityKind K>
struct StringMaker<Stereon::BuilderId<K>> {
    static std::string convert(const Stereon::BuilderId<K>& h) { return std::format("{}", h); }
};

template <Stereon::EntityKind K>
struct StringMaker<Stereon::Oriented<K>> {
    static std::string convert(const Stereon::Oriented<K>& h) { return std::format("{}", h); }
};

template <>
struct StringMaker<Stereon::LogLevel> {
    static std::string convert(const Stereon::LogLevel level) { return std::format("{}", level); }
};

template <>
struct StringMaker<Stereon::JournalId> {
    static std::string convert(const Stereon::JournalId id) { return std::format("{}", id); }
};

} // namespace Catch
