// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// One curve of a closest approach (ClosestApproach.cppm), shared by the
// primal and the products: a line O + u D, a circle or an ellipse
// C + a cos u X + b sin u Y, with its derivatives and bounds on them.
//
// Included in the purview of the closest approach's implementation units,
// after their imports, with <cmath>, <cstdint>, <limits>, <string> and
// <string_view> in their global module fragment. Everything is in an
// anonymous namespace, so each unit gets its own copy.

#pragma once

namespace Stereon::Detail {
namespace {

struct ApproachSide {
    enum class Kind : std::uint8_t { Line, Circle, Ellipse };

    Kind Shape = Kind::Line;
    Point3 Origin; // a line's origin, a conic's centre
    Vec3 D;        // a line's direction
    Vec3 X;        // a conic's frame
    Vec3 Y;
    double Ra = 0.0; // a conic's radii along X and Y
    double Rb = 0.0;
    /// Bounds on |P'|, |P''| and |P'''| over the whole curve.
    double L1 = 0.0;
    double L2 = 0.0;
    double L3 = 0.0;
    /// The largest |P - Origin| of a conic; infinite for a line.
    double Reach = std::numeric_limits<double>::infinity();

    [[nodiscard]] bool IsLine() const { return Shape == Kind::Line; }

    /// P(u) - Origin, formed from the frame.
    [[nodiscard]] Vec3 Offset(const double u) const
    {
        if (IsLine()) {
            return u * D;
        }
        return (Ra * std::cos(u)) * X + (Rb * std::sin(u)) * Y;
    }
    [[nodiscard]] Point3 At(const double u) const { return Origin + Offset(u); }
    [[nodiscard]] Vec3 D1(const double u) const
    {
        if (IsLine()) {
            return D;
        }
        return (-Ra * std::sin(u)) * X + (Rb * std::cos(u)) * Y;
    }
    [[nodiscard]] Vec3 D2(const double u) const
    {
        if (IsLine()) {
            return Vec3{};
        }
        return (-Ra * std::cos(u)) * X + (-Rb * std::sin(u)) * Y;
    }
};

[[nodiscard]] bool IsFiniteVec(const Vec3& v)
{
    return std::isfinite(v.X) && std::isfinite(v.Y) && std::isfinite(v.Z);
}

[[nodiscard]] ApproachSide ConicSideOf(const ApproachSide::Kind kind, const Point3& center,
                                       const double a, const double b, const Vec3& x, const Vec3& y)
{
    ApproachSide side;
    side.Shape = kind;
    side.Origin = center;
    side.X = x;
    side.Y = y;
    side.Ra = a;
    side.Rb = b;
    const double reach = std::max(a, b);
    side.L1 = side.L2 = side.L3 = side.Reach = reach;
    return side;
}

/// The side of a line, a circle or an ellipse; NotSupported for any other
/// curve, InvalidInput for an empty one or a degenerate line or conic.
[[nodiscard]] Result<ApproachSide> SideOf(const Curve& curve, const std::string_view where)
{
    if (curve.IsEmpty()) {
        return Fail(ErrorCode::InvalidInput, std::string(where) + ": an empty curve");
    }
    if (const Line* line = curve.Target<Line>()) {
        const double length = Length(line->Direction);
        if (!(length > 0.0) || !std::isfinite(length) || !IsFiniteVec(line->Origin - Point3{})) {
            return Fail(ErrorCode::InvalidInput,
                        std::string(where) + ": a line needs a finite origin and direction");
        }
        ApproachSide side;
        side.Origin = line->Origin;
        side.D = line->Direction;
        side.L1 = length;
        return side;
    }
    if (const Circle* circle = curve.Target<Circle>()) {
        return ConicSideOf(ApproachSide::Kind::Circle, circle->Center, circle->Radius,
                           circle->Radius, circle->XAxis, circle->YAxis);
    }
    if (const Ellipse* ellipse = curve.Target<Ellipse>()) {
        return ConicSideOf(ApproachSide::Kind::Ellipse, ellipse->Center, ellipse->MajorRadius,
                           ellipse->MinorRadius, ellipse->XAxis, ellipse->YAxis);
    }
    return Fail(ErrorCode::NotSupported,
                std::string(where) +
                    ": lines, circles and ellipses only; a B-spline waits for a subdivision "
                    "proof of its minimum (ADR-0052, addendum of 2026-11-12)");
}

} // namespace
} // namespace Stereon::Detail
