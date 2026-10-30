// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// ADR-0046, section 3 (AC-0046-02): geometry inspection borrows from a named
// wrapper. Target<T>(), GetModel() and the borrowed accessors of the trimmed
// and swept types accept an lvalue receiver and reject a temporary at
// compile time, while owning results (a reversed or transformed copy) stay
// available from temporaries.

#include <catch2/catch_test_macros.hpp>

#include <numbers>
#include <utility>

import Stereon.Surface;

using Stereon::Curve;
using Stereon::Curve2d;
using Stereon::Cylinder;
using Stereon::ExtrusionSurface;
using Stereon::Line;
using Stereon::Plane;
using Stereon::Point3;
using Stereon::RevolutionSurface;
using Stereon::Surface;
using Stereon::TrimmedCurve;
using Stereon::TrimmedSurface;
using Stereon::Vec3;

namespace {

// Each borrow as a template over the wrapper, so that a rejected call is a
// false requirement rather than an error.
template <typename W, typename T>
constexpr bool TargetFromLvalue = requires(const W& w) { w.template Target<T>(); };
template <typename W, typename T>
constexpr bool TargetFromRvalue = requires(W&& w) { std::move(w).template Target<T>(); };
template <typename W>
constexpr bool ModelFromLvalue = requires(const W& w) { w.GetModel(); };
template <typename W>
constexpr bool ModelFromRvalue = requires(W&& w) { std::move(w).GetModel(); };
template <typename W>
constexpr bool BasisFromLvalue = requires(const W& w) { w.GetBasis(); };
template <typename W>
constexpr bool BasisFromRvalue = requires(W&& w) { std::move(w).GetBasis(); };
template <typename W>
constexpr bool ProfileFromRvalue = requires(W&& w) { std::move(w).GetProfile(); };
template <typename W>
constexpr bool AxisFromRvalue = requires(W&& w) { std::move(w).GetAxis(); };
template <typename W>
constexpr bool ReversedFromRvalue = requires(W&& w) { Reversed(std::move(w)); };

} // namespace

TEST_CASE("geometry inspection borrows from an lvalue wrapper only", "[surface][curve][borrowing]")
{
    STATIC_REQUIRE(TargetFromLvalue<Curve, Line>);
    STATIC_REQUIRE(ModelFromLvalue<Curve>);
    STATIC_REQUIRE_FALSE(TargetFromRvalue<Curve, Line>);
    STATIC_REQUIRE_FALSE(ModelFromRvalue<Curve>);
    STATIC_REQUIRE_FALSE(TargetFromRvalue<Curve2d, Stereon::Line2d>);

    STATIC_REQUIRE(TargetFromLvalue<Surface, Plane>);
    STATIC_REQUIRE(ModelFromLvalue<Surface>);
    STATIC_REQUIRE_FALSE(TargetFromRvalue<Surface, Plane>);
    STATIC_REQUIRE_FALSE(ModelFromRvalue<Surface>);

    STATIC_REQUIRE(BasisFromLvalue<TrimmedCurve>);
    STATIC_REQUIRE_FALSE(BasisFromRvalue<TrimmedCurve>);
    STATIC_REQUIRE_FALSE(BasisFromRvalue<TrimmedSurface>);
    STATIC_REQUIRE_FALSE(ProfileFromRvalue<ExtrusionSurface>);
    STATIC_REQUIRE_FALSE(ProfileFromRvalue<RevolutionSurface>);
    STATIC_REQUIRE_FALSE(AxisFromRvalue<RevolutionSurface>);

    // Owning results remain available from temporaries.
    STATIC_REQUIRE(ReversedFromRvalue<Curve>);
    STATIC_REQUIRE(ReversedFromRvalue<Surface>);
}

TEST_CASE("a named wrapper's borrow lives as long as the wrapper", "[surface][curve][borrowing]")
{
    const Line line{Point3{1.0, 2.0, 3.0}, Vec3{0.0, 0.0, 1.0}};
    const Curve curve{line};
    const Line* borrowed = curve.Target<Line>();
    REQUIRE(borrowed != nullptr);
    REQUIRE(borrowed->Origin == line.Origin);
    REQUIRE(curve.Target<Stereon::Circle>() == nullptr);

    // An owning copy taken from a temporary is inspected once it is named.
    const Curve reversed = Reversed(Curve{line});
    const Line* reversedLine = reversed.Target<Line>();
    REQUIRE(reversedLine != nullptr);
    // Read through the borrow only where it is known to be one. The model
    // behind `reversed` comes from a virtual call, so the optimiser cannot
    // rule out the null branch of Target, and REQUIRE's failure path is
    // opaque to it: without the guard, -O3 reports a null dereference.
    if (reversedLine != nullptr) {
        REQUIRE(reversedLine->Direction == -line.Direction);
    }

    const Surface surface{Cylinder{Point3{}, 2.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, true}};
    const Surface reversedSurface = Reversed(surface);
    REQUIRE(reversedSurface.Target<Cylinder>() != nullptr);
    REQUIRE(surface.Target<Plane>() == nullptr);
}
