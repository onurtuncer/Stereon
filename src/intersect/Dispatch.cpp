// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

// The contract anchor (docs/toolchain/gcc16-feature-probe.md).
#include <Stereon/Core/Contracts.hpp>

#include <cmath>
#include <cstddef>
#include <expected>
#include <string>
#include <string_view>

module Stereon.Intersect;

import :Dispatch;

namespace Stereon {

namespace {

/// A requested domain: not empty and with no NaN bound. Unbounded is
/// allowed: an unbounded surface asked over its whole domain.
[[nodiscard]] bool IsDomain(const Interval& domain) noexcept
{
    return !std::isnan(domain.GetLower()) && !std::isnan(domain.GetUpper()) && !domain.IsEmpty();
}

[[nodiscard]] Result<void> CheckDomains(const bool valid, const std::string_view what)
{
    if (!valid) {
        return Fail(ErrorCode::InvalidInput,
                    std::string("Intersect: the requested ") + std::string(what) +
                        " domains must be non-empty and have no NaN bound");
    }
    return {};
}

} // namespace

Result<CurveCurveIntersection> Intersect(const Curve& a, const Interval& domainA, const Curve& b,
                                         const Interval& domainB, const Context& ctx)
{
    if (a.IsEmpty() || b.IsEmpty()) {
        return Fail(ErrorCode::InvalidInput, "cannot intersect an empty curve");
    }
    if (const Result<void> ok = CheckDomains(IsDomain(domainA) && IsDomain(domainB), "curve");
        !ok) {
        return std::unexpected(ok.error());
    }
    const std::size_t i = BuiltinIndexOf(a);
    const std::size_t j = BuiltinIndexOf(b);
    if (i == BuiltinCurveTypes::Size || j == BuiltinCurveTypes::Size) {
        return Fail(ErrorCode::NotSupported,
                    "curve intersection needs built-in curve types until user types can declare "
                    "the capabilities of tiers 2 and 3 (ADR-0003)");
    }
    return Detail::CurvePairTable[i][j](a, domainA, b, domainB, ctx);
}

Result<CurveCurveIntersection> Intersect(const Curve& a, const Curve& b, const Context& ctx)
{
    if (a.IsEmpty() || b.IsEmpty()) {
        return Fail(ErrorCode::InvalidInput, "cannot intersect an empty curve");
    }
    return Intersect(a, Domain(a), b, Domain(b), ctx);
}

Result<SurfaceSurfaceIntersection> Intersect(const Surface& a, const Box2& domainA,
                                             const Surface& b, const Box2& domainB,
                                             const Context& ctx)
{
    if (a.IsEmpty() || b.IsEmpty()) {
        return Fail(ErrorCode::InvalidInput, "cannot intersect an empty surface");
    }
    const bool valid =
        IsDomain(domainA.X) && IsDomain(domainA.Y) && IsDomain(domainB.X) && IsDomain(domainB.Y);
    if (const Result<void> ok = CheckDomains(valid, "surface"); !ok) {
        return std::unexpected(ok.error());
    }
    const std::size_t i = BuiltinIndexOf(a);
    const std::size_t j = BuiltinIndexOf(b);
    if (i == BuiltinSurfaceTypes::Size || j == BuiltinSurfaceTypes::Size) {
        return Fail(ErrorCode::NotSupported,
                    "surface intersection needs built-in surface types until user types can "
                    "declare the capabilities of tiers 2 and 3 (ADR-0003)");
    }
    return Detail::SurfacePairTable[i][j](a, domainA, b, domainB, ctx);
}

Result<SurfaceSurfaceIntersection> Intersect(const Surface& a, const Surface& b, const Context& ctx)
{
    if (a.IsEmpty() || b.IsEmpty()) {
        return Fail(ErrorCode::InvalidInput, "cannot intersect an empty surface");
    }
    return Intersect(a, Domain(a), b, Domain(b), ctx);
}

} // namespace Stereon
