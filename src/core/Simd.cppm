// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <cstddef>
#include <version>

// Use std::simd once the standard library ships it; until then xsimd. The two
// APIs differ, so this partition defines a small common API instead of
// re-exporting either one.
#if defined(__cpp_lib_simd)
#include <simd>
#define STEREON_SIMD_USES_STD 1
#else
#include <xsimd/xsimd.hpp>
#define STEREON_SIMD_USES_STD 0
#endif

export module Stereon.Core:Simd;

/// Data-parallel arithmetic for evaluators (ADR-0007). Stereon code uses
/// Stereon::Simd::Vec and the functions below, never std::simd or xsimd
/// directly, so the switch to the standard library is a change in this file.
///
/// Vec<T> holds Width<T> lanes at the native width chosen at compile time.
/// Arithmetic (+, -, *, /) and comparisons use the operators of the underlying
/// type. All loads and stores are unaligned.
///
/// The std::simd branch follows the C++26 working draft but cannot be built
/// until a standard library defines __cpp_lib_simd (none does as of GCC 16).
export namespace Stereon::Simd {

/// True when std::simd is used; false while the xsimd polyfill is.
inline constexpr bool UsesStandardLibrary = STEREON_SIMD_USES_STD == 1;

#if STEREON_SIMD_USES_STD

template <class T>
using Vec = std::simd::vec<T>;

template <class T>
inline constexpr std::size_t Width = static_cast<std::size_t>(Vec<T>::size());

template <class T>
[[nodiscard]] inline Vec<T> Load(const T* source)
{
    return std::simd::unchecked_load<Vec<T>>(source, Width<T>);
}

template <class T>
inline void Store(const Vec<T>& value, T* destination)
{
    std::simd::unchecked_store(value, destination, Width<T>);
}

template <class T>
[[nodiscard]] inline T ReduceAdd(const Vec<T>& value)
{
    return std::simd::reduce(value);
}

template <class T>
[[nodiscard]] inline Vec<T> Min(const Vec<T>& a, const Vec<T>& b)
{
    return std::simd::min(a, b);
}

template <class T>
[[nodiscard]] inline Vec<T> Max(const Vec<T>& a, const Vec<T>& b)
{
    return std::simd::max(a, b);
}

/// Lane-wise absolute value.
template <class T>
[[nodiscard]] inline Vec<T> Abs(const Vec<T>& value)
{
    return std::simd::abs(value);
}

/// Lane-wise square root.
template <class T>
[[nodiscard]] inline Vec<T> Sqrt(const Vec<T>& value)
{
    return std::simd::sqrt(value);
}

/// Lane-wise fused multiply-add, a * b + c with one rounding.
template <class T>
[[nodiscard]] inline Vec<T> Fma(const Vec<T>& a, const Vec<T>& b, const Vec<T>& c)
{
    return std::simd::fma(a, b, c);
}

#else

template <class T>
using Vec = xsimd::batch<T>;

template <class T>
inline constexpr std::size_t Width = Vec<T>::size;

template <class T>
[[nodiscard]] inline Vec<T> Load(const T* source)
{
    return xsimd::load_unaligned(source);
}

template <class T>
inline void Store(const Vec<T>& value, T* destination)
{
    value.store_unaligned(destination);
}

template <class T>
[[nodiscard]] inline T ReduceAdd(const Vec<T>& value)
{
    return xsimd::reduce_add(value);
}

template <class T>
[[nodiscard]] inline Vec<T> Min(const Vec<T>& a, const Vec<T>& b)
{
    return xsimd::min(a, b);
}

template <class T>
[[nodiscard]] inline Vec<T> Max(const Vec<T>& a, const Vec<T>& b)
{
    return xsimd::max(a, b);
}

/// Lane-wise absolute value.
template <class T>
[[nodiscard]] inline Vec<T> Abs(const Vec<T>& value)
{
    return xsimd::abs(value);
}

/// Lane-wise square root.
template <class T>
[[nodiscard]] inline Vec<T> Sqrt(const Vec<T>& value)
{
    return xsimd::sqrt(value);
}

/// Lane-wise fused multiply-add, a * b + c with one rounding.
template <class T>
[[nodiscard]] inline Vec<T> Fma(const Vec<T>& a, const Vec<T>& b, const Vec<T>& c)
{
    return xsimd::fma(a, b, c);
}

#endif

/// Vec<T> with every lane set to value.
template <class T>
[[nodiscard]] inline Vec<T> Broadcast(T value)
{
    return Vec<T>(value);
}

} // namespace Stereon::Simd
