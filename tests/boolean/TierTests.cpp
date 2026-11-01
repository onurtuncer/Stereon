// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// ADR-0003's capability tiers, checked where they bite: a user surface that
// declares only tier 1 (evaluation, derivatives, domain, bounds, projection,
// reversal, transformation) cannot take part in a Boolean, which needs the
// bounds of tier 2 and the certified intersections of tier 3. The Boolean
// must refuse with NotSupported, never guess. The second half of the ADR's
// check, the same surface with tiers 2 and 3 taking part, waits for the
// BoundedGeometry and IntersectableGeometry concepts (ADR-0010).

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

import Stereon.Boolean;

#include "Fixtures.hpp"

using BooleanTest::ctx;
using Stereon::BooleanOperation;
using Stereon::Box2;
using Stereon::Box3;
using Stereon::ErrorCode;
using Stereon::Plane;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Shape;
using Stereon::ShapeBuilder;
using Stereon::Transform;
using Stereon::Vec3;

namespace {

/// A plane behind the tier-1 operations only: a user's type, so the kernel
/// knows nothing about it beyond what these functions say.
struct UserPlane {
    Plane Base;
};
Point3 Eval(const UserPlane& s, const Point2& uv)
{
    return Eval(s.Base, uv);
}
Vec3 Deriv(const UserPlane& s, const Point2& uv, const int du, const int dv)
{
    return Deriv(s.Base, uv, du, dv);
}
Box2 Domain(const UserPlane& s)
{
    return Domain(s.Base);
}
Box3 BoundingBox(const UserPlane& s)
{
    return BoundingBox(s.Base);
}
Point2 Project(const UserPlane& s, const Point3& p)
{
    return Project(s.Base, p);
}
UserPlane Reversed(const UserPlane& s)
{
    return {Reversed(s.Base)};
}
UserPlane Transformed(const UserPlane& s, const Transform& x)
{
    return {Transformed(s.Base, x)};
}

static_assert(Stereon::SurfaceGeometry<UserPlane>);

/// The box [0, 10]^3 with its top face's plane replaced by the same plane
/// as a UserPlane: the same geometry, only its tier differs.
Shape BoxWithUserTop()
{
    const Shape box = BooleanTest::Box();
    const Stereon::FaceId top = BooleanTest::PlaneFaceAlong(box, Vec3{0.0, 0.0, 1.0});
    const Stereon::Surface surface = Stereon::SurfaceOf(box, top);
    const Plane plane = *surface.Target<Plane>();
    auto edit = ShapeBuilder::Open(box);
    REQUIRE(edit.has_value());
    const auto user = edit->AddSurface(Stereon::Surface(UserPlane{plane}));
    REQUIRE(user.has_value());
    REQUIRE(edit->SetSurface(edit->Handle(top), *user).has_value());
    auto frozen = Stereon::Freeze(std::move(*edit));
    REQUIRE(frozen.has_value());
    return std::move(*frozen);
}

} // namespace

TEST_CASE("a tier-1 user surface in a Boolean is refused, not guessed", "[boolean][adr-0003]")
{
    const Shape a = BoxWithUserTop();
    // A box across the user face, so that its intersections are needed.
    const auto b =
        MakeBoxSolid(Box3::FromCorners(Point3{5.0, 5.0, 5.0}, Point3{15.0, 15.0, 15.0}), ctx);
    REQUIRE(b.has_value());
    for (const BooleanOperation operation :
         {BooleanOperation::Union, BooleanOperation::Cut, BooleanOperation::Common}) {
        const auto result = Boolean(a, b->Shape, operation, ctx);
        REQUIRE_FALSE(result.has_value());
        INFO(std::string(result.error().GetMessage()));
        REQUIRE(result.error().GetCode() == ErrorCode::NotSupported);
    }
    // The same box with its own plane goes through: the refusal is the tier's.
    REQUIRE(Boolean(BooleanTest::Box(), b->Shape, BooleanOperation::Union, ctx).has_value());
}
