// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// Error propagation macros for Stereon::Result (ADR-0006). Modules cannot
/// export macros, so this is a header; include it in the global module
/// fragment of a module unit.
///
/// Both macros are used inside a function that returns a Result. On an error
/// they return it to the caller; the value types of the two Results may differ.
///
///     Result<Shape> MakeBlock(const Context& ctx, double size)
///     {
///         STEREON_TRY(const Profile profile, MakeSquare(ctx, size));
///         STEREON_TRY_VOID(CheckProfile(ctx, profile));
///         return Extrude(ctx, profile, size);
///     }
///
/// STEREON_TRY expands to several statements and declares a variable in the
/// enclosing scope, so it cannot be the body of an unbraced if or loop. The
/// declaration must not contain a top-level comma: write
/// STEREON_TRY(auto pair, ...) and unpack the pair on the next line rather
/// than STEREON_TRY(auto [a, b], ...).
///
/// The macros use only standard C++, not statement expressions, so they also
/// build with MSVC.

#pragma once

#include <expected>
#include <utility>

#define STEREON_TRY_CONCAT_IMPL(a, b) a##b
#define STEREON_TRY_CONCAT(a, b) STEREON_TRY_CONCAT_IMPL(a, b)

#define STEREON_TRY_IMPL(temp, declaration, ...)                                                   \
    auto temp = (__VA_ARGS__);                                                                     \
    if (!temp.has_value()) {                                                                       \
        return std::unexpected(std::move(temp).error());                                           \
    }                                                                                              \
    declaration = *std::move(temp)

/// Evaluates a Result; on success initialises `declaration` with its value,
/// on error returns the error from the enclosing function.
#define STEREON_TRY(declaration, ...)                                                              \
    STEREON_TRY_IMPL(STEREON_TRY_CONCAT(stereonTryResult_, __COUNTER__), declaration, __VA_ARGS__)

#define STEREON_TRY_VOID_IMPL(temp, ...)                                                           \
    do {                                                                                           \
        auto temp = (__VA_ARGS__);                                                                 \
        if (!temp.has_value()) {                                                                   \
            return std::unexpected(std::move(temp).error());                                       \
        }                                                                                          \
    } while (false)

/// Evaluates a Result and discards its value; on error returns the error from
/// the enclosing function.
#define STEREON_TRY_VOID(...)                                                                      \
    STEREON_TRY_VOID_IMPL(STEREON_TRY_CONCAT(stereonTryResult_, __COUNTER__), __VA_ARGS__)
