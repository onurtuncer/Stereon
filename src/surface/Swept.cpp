// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

// The contract anchor (docs/toolchain/gcc16-feature-probe.md).
#include <Stereon/Core/Contracts.hpp>
#include <Stereon/Core/Try.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <optional>
#include <type_traits>
#include <utility>
#include <vector>

module Stereon.Surface;

import :BSplineSurface;
import :Swept;

namespace Stereon::Detail {

struct ExtrusionProfile {
    /// The profile's NURBS form over its domain.
    BSplineCurve Form;
    /// The form with its poles moved onto the plane through the origin
    /// normal to the direction, P - (P . n) n, with the same knots and
    /// weights: the profile seen along the direction, on the form's
    /// parameter.
    BSplineCurve Flattened;
    CurveProjector Projector; // of Flattened
    /// True when the form is parametrised as the profile (a line, a
    /// parabola or a B-spline); otherwise a point of the form is mapped
    /// back to the profile by projecting onto it.
    bool SameParameter = false;
};

struct RevolvedPiece {
    /// The NURBS form of the full turn (MakeRevolutionSurface): its u is the
    /// rational arc parameter, its v the form's parameter.
    BSplineSurface Surface;
    SurfaceProjector Projector;
    /// As ExtrusionProfile::SameParameter, for v.
    bool SameParameter = false;
};

} // namespace Stereon::Detail

namespace Stereon {

namespace {

constexpr double Infinity = std::numeric_limits<double>::infinity();

bool Finite(const Vec3& v)
{
    return std::isfinite(v.X) && std::isfinite(v.Y) && std::isfinite(v.Z);
}

bool FiniteRange(const Interval range)
{
    return std::isfinite(range.GetLower()) && std::isfinite(range.GetUpper()) &&
           range.GetLower() < range.GetUpper() &&
           std::isfinite(range.GetUpper() - range.GetLower());
}

Result<void> CheckProfile(const BSplineCurve& profile)
{
    STEREON_TRY_VOID(
        MakeBSplineCurve(profile.Degree, profile.Knots, profile.Poles, profile.Weights));
    return {};
}

Result<BSplineSurface> CheckedSurface(BSplineSurface surface)
{
    return MakeBSplineSurface(surface.DegreeU, surface.DegreeV, std::move(surface.KnotsU),
                              std::move(surface.KnotsV), std::move(surface.Poles),
                              std::move(surface.Weights));
}

// --- The profile's NURBS form ---------------------------------------------

/// True when the profile's NURBS form keeps its parametrisation: a line, a
/// parabola or a B-spline (src/curve/Conversion.cppm); a conic's rational
/// arcs agree with it only at their knots.
bool KeepsParameter(const Curve& profile)
{
    return VisitBuiltin(profile, []<typename G>(const G&) {
        return std::is_same_v<G, Line> || std::is_same_v<G, Parabola> ||
               std::is_same_v<G, BSplineCurve>;
    });
}

/// The profile's exact NURBS form over its domain, when the kernel knows
/// one and the domain is finite.
std::optional<BSplineCurve> FormOf(const Curve& profile)
{
    if (profile.IsEmpty()) {
        return std::nullopt;
    }
    const Interval domain = Domain(profile);
    if (!FiniteRange(domain)) {
        return std::nullopt;
    }
    Result<BSplineCurve> form = ToBSpline(profile, domain);
    if (!form) {
        return std::nullopt;
    }
    return *std::move(form);
}

/// The part of v along the unit vector n, and the rest.
Vec3 Along(const Vec3& v, const Vec3& n)
{
    return Dot(v, n) * n;
}

Vec3 Across(const Vec3& v, const Vec3& n)
{
    return v - Along(v, n);
}

/// u brought into the profile's domain: by whole periods first for a
/// periodic profile, then clamped.
double IntoDomain(const double u, const Curve& profile)
{
    return NearestInRange(u, Domain(profile), Period(profile));
}

// --- Extrusion --------------------------------------------------------------

/// The profile parameter of a closest point of the extrusion to p.
double ExtrusionParameter(const ExtrusionSurface& s, const Point3& p)
{
    const Curve& profile = s.GetProfile();
    const Vec3 n = Normalized(s.GetDirection());
    if (const Detail::ExtrusionProfile* data = s.GetProjection()) {
        // Along the direction the surface is the flattened form; the closest
        // point of the surface lies over the closest point of that curve to
        // the point flattened the same way.
        const double t = data->Projector.Project(ToPoint(Across(ToVec(p), n)));
        if (data->SameParameter) {
            return IntoDomain(t, profile);
        }
        // The form's point lies on the profile, a conic that does not cross
        // itself, so projecting it back finds its own parameter.
        return IntoDomain(Project(profile, Eval(data->Form, t)), profile);
    }
    if (const Line* line = profile.Target<Line>()) {
        // A line extrudes to a plane, or to itself along the direction.
        const Vec3 e = Across(line->Direction, n);
        const double ee = Dot(e, e);
        if (ee > 0.0) {
            return Dot(Across(p - line->Origin, n), e) / ee;
        }
    }
    // Approximate: the profile's own closest point.
    return IntoDomain(Project(profile, p), profile);
}

// --- Revolution ---------------------------------------------------------------

/// (cos u, sin u) turned by `quarterTurns` quarter turns: the derivatives
/// of the turn with respect to the angle.
struct Turn {
    double C;
    double S;
};

Turn Turned(const double u, const int quarterTurns)
{
    const double c = Cos(u);
    const double s = Sin(u);
    switch (quarterTurns & 3) {
    case 0:
        return {c, s};
    case 1:
        return {-s, c};
    case 2:
        return {-c, -s};
    default:
        return {s, -c};
    }
}

/// S^(k,l) relative to the origin: the radial part of C^(l)(v) turned by u
/// and k quarter turns, plus its axial part when k is 0.
Vec3 Revolved(const RevolutionSurface& s, const Point2& uv, const int du, const int dv)
{
    const Vec3& n = s.GetAxis();
    const Vec3 vec =
        dv == 0 ? Eval(s.GetProfile(), uv.Y) - s.GetOrigin() : Deriv(s.GetProfile(), uv.Y, dv);
    const Vec3 axial = Along(vec, n);
    const Vec3 radial = vec - axial;
    const Turn t = Turned(uv.X, du);
    const Vec3 turned = t.C * radial + t.S * Cross(n, radial);
    return du == 0 ? turned + axial : turned;
}

/// The angle about the axis from the radial direction of the profile point
/// C(v) to that of `point`, in [0, 2 pi); 0 where either lies on the axis.
double AngleOf(const RevolutionSurface& s, const Point3& point, const double v)
{
    const Vec3& n = s.GetAxis();
    const Vec3 rq = Across(Eval(s.GetProfile(), v) - s.GetOrigin(), n);
    const Vec3 rp = Across(point - s.GetOrigin(), n);
    if (!(Dot(rq, rq) > 0.0) || !(Dot(rp, rp) > 0.0)) {
        return 0.0;
    }
    const double angle = Atan2(Dot(rp, Cross(n, rq)), Dot(rp, rq));
    return angle < 0.0 ? angle + CirclePeriod : angle;
}

/// p turned about the axis by -angle: into the half-plane of C(v) when
/// angle is AngleOf(s, p, v).
Point3 TurnedBack(const RevolutionSurface& s, const Point3& p, const double angle)
{
    const Vec3& n = s.GetAxis();
    const Vec3 offset = p - s.GetOrigin();
    const Vec3 radial = Across(offset, n);
    const Turn t = Turned(-angle, 0);
    return s.GetOrigin() + Along(offset, n) + t.C * radial + t.S * Cross(n, radial);
}

Point2 ProjectRevolution(const RevolutionSurface& s, const Point3& p)
{
    const Curve& profile = s.GetProfile();
    if (const Detail::RevolvedPiece* piece = s.GetPiece()) {
        const Point2 st = piece->Projector.Project(p);
        double v = st.Y;
        if (!piece->SameParameter) {
            // At angle 0 the turn's first arc pole has weight 1, so the piece
            // at (0, t) is the form's point at t, which lies on the profile.
            v = Project(profile, Eval(piece->Surface, Point2{0.0, st.Y}));
        }
        v = IntoDomain(v, profile);
        return {AngleOf(s, Eval(piece->Surface, st), v), v};
    }
    // Approximate: the profile's closest point to p turned into its
    // half-plane, a few times over, then the angle of p itself, where the
    // circle traced by that profile point comes closest.
    double v = IntoDomain(Project(profile, p), profile);
    for (int round = 0; round < 3; ++round) {
        v = IntoDomain(Project(profile, TurnedBack(s, p, AngleOf(s, p, v))), profile);
    }
    return {AngleOf(s, p, v), v};
}

/// The box of the cylinder about the axis that holds `box` turned: axial
/// extent of the box's corners, radius their largest distance from the
/// axis. Infinite where the box is.
Box3 RevolvedBox(const RevolutionSurface& s, const Box3& box)
{
    std::array<double, 3> low{};
    std::array<double, 3> high{};
    for (std::size_t d = 0; d < 3; ++d) {
        low[d] = box[d].GetLower();
        high[d] = box[d].GetUpper();
        if (!std::isfinite(low[d]) || !std::isfinite(high[d])) {
            return {Interval{-Infinity, Infinity}, Interval{-Infinity, Infinity},
                    Interval{-Infinity, Infinity}};
        }
    }
    const Vec3& n = s.GetAxis();
    double zLow = Infinity;
    double zHigh = -Infinity;
    double radius = 0.0;
    for (int corner = 0; corner < 8; ++corner) {
        const Point3 c{(corner & 1) != 0 ? high[0] : low[0], (corner & 2) != 0 ? high[1] : low[1],
                       (corner & 4) != 0 ? high[2] : low[2]};
        const Vec3 offset = c - s.GetOrigin();
        const double z = Dot(offset, n);
        zLow = std::min(zLow, z);
        zHigh = std::max(zHigh, z);
        radius = std::max(radius, Length(Across(offset, n)));
    }
    Box3 result;
    for (std::size_t d = 0; d < 3; ++d) {
        const double width = radius * Sqrt(std::max(0.0, 1.0 - n[d] * n[d]));
        const double a = s.GetOrigin()[d] + zLow * n[d];
        const double b = s.GetOrigin()[d] + zHigh * n[d];
        result[d] = Interval{std::min(a, b) - width, std::max(a, b) + width};
    }
    return result;
}

} // namespace

namespace Detail {

std::shared_ptr<const ExtrusionProfile> MakeExtrusionProfile(const Curve& profile,
                                                             const Vec3& direction)
{
    std::optional<BSplineCurve> form = FormOf(profile);
    if (!form || !Finite(direction) || !(Dot(direction, direction) > 0.0)) {
        return nullptr;
    }
    const Vec3 n = Normalized(direction);
    BSplineCurve flattened = *form;
    for (Point3& pole : flattened.Poles) {
        pole = ToPoint(Across(ToVec(pole), n));
    }
    CurveProjector projector(flattened);
    return std::make_shared<const ExtrusionProfile>(ExtrusionProfile{
        *std::move(form), std::move(flattened), std::move(projector), KeepsParameter(profile)});
}

std::shared_ptr<const RevolvedPiece> MakeRevolvedPiece(const Curve& profile, const Point3& origin,
                                                       const Vec3& axis)
{
    const std::optional<BSplineCurve> form = FormOf(profile);
    if (!form) {
        return nullptr;
    }
    Result<BSplineSurface> surface =
        MakeRevolutionSurface(*form, origin, axis, Interval{0.0, CirclePeriod});
    if (!surface) {
        return nullptr;
    }
    SurfaceProjector projector(*surface);
    return std::make_shared<const RevolvedPiece>(
        RevolvedPiece{*std::move(surface), std::move(projector), KeepsParameter(profile)});
}

} // namespace Detail

// --- Extrusion ----------------------------------------------------------------

Result<ExtrusionSurface> MakeExtrusion(const Curve& profile, const Vec3& direction)
{
    if (profile.IsEmpty()) {
        return Fail(ErrorCode::InvalidInput, "MakeExtrusion: the profile curve is empty");
    }
    STEREON_TRY_VOID(Detail::CheckFinite(direction, "the direction", "MakeExtrusion"));
    if (!(Dot(direction, direction) > 0.0)) {
        return Fail(ErrorCode::InvalidInput, "MakeExtrusion: the direction must not be zero");
    }
    return ExtrusionSurface(profile, direction);
}

Point3 Eval(const ExtrusionSurface& surface, const Point2& uv)
{
    Detail::RequireFinite(uv, "Eval(ExtrusionSurface)", "a non-finite parameter");
    Detail::RequireInDomain(Domain(surface.GetProfile()), uv.X, Period(surface.GetProfile()),
                            "Eval(ExtrusionSurface)");
    return Eval(surface.GetProfile(), uv.X) + uv.Y * surface.GetDirection();
}

Vec3 Deriv(const ExtrusionSurface& surface, const Point2& uv, const int du, const int dv)
{
    Detail::RequireFinite(uv, "Deriv(ExtrusionSurface)", "a non-finite parameter");
    Detail::RequireInDomain(Domain(surface.GetProfile()), uv.X, Period(surface.GetProfile()),
                            "Deriv(ExtrusionSurface)");
    Detail::Require(du >= 0 && dv >= 0 && du + dv >= 1, "Deriv(ExtrusionSurface)",
                    "a derivative order that does not exist");
    if (dv == 0) {
        return Deriv(surface.GetProfile(), uv.X, du);
    }
    return du == 0 && dv == 1 ? surface.GetDirection() : Vec3{};
}

Box2 Domain(const ExtrusionSurface& surface)
{
    return {Domain(surface.GetProfile()), Interval{-Infinity, Infinity}};
}

Vec2 Periods(const ExtrusionSurface& surface)
{
    return {Period(surface.GetProfile()), 0.0};
}

Box3 BoundingBox(const ExtrusionSurface& surface)
{
    const Box3 box = BoundingBox(surface.GetProfile());
    const Vec3& d = surface.GetDirection();
    Box3 swept;
    for (std::size_t i = 0; i < 3; ++i) {
        swept[i] = d[i] != 0.0 ? Interval{-Infinity, Infinity} : box[i];
    }
    return swept;
}

Point2 Project(const ExtrusionSurface& surface, const Point3& p)
{
    Detail::RequireFinite(p, "Project(ExtrusionSurface)", "a non-finite point");
    const double u = ExtrusionParameter(surface, p);
    const Vec3& d = surface.GetDirection();
    return {u, Dot(p - Eval(surface.GetProfile(), u), d) / Dot(d, d)};
}

ExtrusionSurface Reversed(const ExtrusionSurface& surface)
{
    return {Reversed(surface.GetProfile()), surface.GetDirection()};
}

ExtrusionSurface Transformed(const ExtrusionSurface& surface, const Transform& x)
{
    return {Transformed(surface.GetProfile(), x), Transformed(surface.GetDirection(), x)};
}

Result<BSplineSurface> ToBSpline(const ExtrusionSurface& surface, const Box2 range)
{
    STEREON_TRY(const BSplineCurve form, ToBSpline(surface.GetProfile(), range.X));
    return MakeExtrusionSurface(form, surface.GetDirection(), range.Y);
}

// --- Revolution -----------------------------------------------------------------

Result<RevolutionSurface> MakeRevolution(const Curve& profile, const Point3& origin,
                                         const Vec3& axis)
{
    if (profile.IsEmpty()) {
        return Fail(ErrorCode::InvalidInput, "MakeRevolution: the profile curve is empty");
    }
    STEREON_TRY_VOID(Detail::CheckFinite(origin, "the origin", "MakeRevolution"));
    STEREON_TRY_VOID(Detail::CheckFinite(axis, "the axis", "MakeRevolution"));
    if (!(Dot(axis, axis) > 0.0)) {
        return Fail(ErrorCode::InvalidInput, "MakeRevolution: the axis must not be zero");
    }
    return RevolutionSurface(profile, origin, axis);
}

Point3 Eval(const RevolutionSurface& surface, const Point2& uv)
{
    Detail::RequireFinite(uv, "Eval(RevolutionSurface)", "a non-finite parameter");
    Detail::RequireInDomain(Domain(surface.GetProfile()), uv.Y, Period(surface.GetProfile()),
                            "Eval(RevolutionSurface)");
    return surface.GetOrigin() + Revolved(surface, uv, 0, 0);
}

Vec3 Deriv(const RevolutionSurface& surface, const Point2& uv, const int du, const int dv)
{
    Detail::RequireFinite(uv, "Deriv(RevolutionSurface)", "a non-finite parameter");
    Detail::RequireInDomain(Domain(surface.GetProfile()), uv.Y, Period(surface.GetProfile()),
                            "Deriv(RevolutionSurface)");
    Detail::Require(du >= 0 && dv >= 0 && du + dv >= 1, "Deriv(RevolutionSurface)",
                    "a derivative order that does not exist");
    return Revolved(surface, uv, du, dv);
}

Box2 Domain(const RevolutionSurface& surface)
{
    return {Interval{0.0, CirclePeriod}, Domain(surface.GetProfile())};
}

Vec2 Periods(const RevolutionSurface& surface)
{
    return {CirclePeriod, Period(surface.GetProfile())};
}

Box3 BoundingBox(const RevolutionSurface& surface)
{
    const Detail::RevolvedPiece* piece = surface.GetPiece();
    return piece != nullptr ? BoundingBox(piece->Surface)
                            : RevolvedBox(surface, BoundingBox(surface.GetProfile()));
}

Point2 Project(const RevolutionSurface& surface, const Point3& p)
{
    Detail::RequireFinite(p, "Project(RevolutionSurface)", "a non-finite point");
    return ProjectRevolution(surface, p);
}

RevolutionSurface Reversed(const RevolutionSurface& surface)
{
    return {surface.GetProfile(), surface.GetOrigin(), -surface.GetAxis()};
}

RevolutionSurface Transformed(const RevolutionSurface& surface, const Transform& x)
{
    Vec3 axis = Transformed(surface.GetAxis(), x);
    if (Determinant(x.Linear) < 0.0) {
        axis = -axis; // a reflection reverses the sense of a turn
    }
    return {Transformed(surface.GetProfile(), x), Transformed(surface.GetOrigin(), x), axis};
}

Result<BSplineSurface> ToBSpline(const RevolutionSurface& surface, const Box2 range)
{
    STEREON_TRY(const BSplineCurve form, ToBSpline(surface.GetProfile(), range.Y));
    return MakeRevolutionSurface(form, surface.GetOrigin(), surface.GetAxis(), range.X);
}

// --- The exact NURBS forms ------------------------------------------------------

Result<BSplineSurface> MakeExtrusionSurface(const BSplineCurve& profile, const Vec3& direction,
                                            const Interval range)
{
    STEREON_TRY_VOID(CheckProfile(profile));
    if (!Finite(direction) || direction == Vec3{} || !FiniteRange(range)) {
        return Fail(ErrorCode::InvalidInput,
                    "MakeExtrusionSurface: expected a finite nonzero direction and finite range");
    }
    const double low = range.GetLower();
    const double high = range.GetUpper();
    BSplineSurface surface{profile.Degree, 1, profile.Knots, {low, low, high, high}, {}, {}};
    for (std::size_t i = 0; i < profile.Poles.size(); ++i) {
        surface.Poles.push_back(profile.Poles[i] + low * direction);
        surface.Poles.push_back(profile.Poles[i] + high * direction);
        if (IsRational(profile)) {
            surface.Weights.push_back(profile.Weights[i]);
            surface.Weights.push_back(profile.Weights[i]);
        }
    }
    return CheckedSurface(std::move(surface));
}

Result<BSplineSurface> MakeRevolutionSurface(const BSplineCurve& profile, const Point3& origin,
                                             const Vec3& axis, const Interval angles)
{
    STEREON_TRY_VOID(CheckProfile(profile));
    const double length = std::hypot(axis.X, axis.Y, axis.Z);
    if (!Finite(origin - Point3{}) || !Finite(axis) || !std::isfinite(length) || length == 0.0 ||
        !FiniteRange(angles) || angles.GetUpper() - angles.GetLower() > CirclePeriod) {
        return Fail(ErrorCode::InvalidInput,
                    "MakeRevolutionSurface: expected a finite axis and at most one turn");
    }
    const Vec3 normal = axis / length;
    const Detail::UnitArc arc = Detail::MakeUnitArc(angles.GetLower(), angles.GetUpper());
    BSplineSurface surface{2, profile.Degree, arc.Knots, profile.Knots, {}, {}};
    for (std::size_t i = 0; i < arc.Poles.size(); ++i) {
        for (std::size_t j = 0; j < profile.Poles.size(); ++j) {
            const Vec3 offset = profile.Poles[j] - origin;
            const Vec3 axial = Dot(offset, normal) * normal;
            const Vec3 radial = offset - axial;
            surface.Poles.push_back(origin + axial + arc.Poles[i].X * radial +
                                    arc.Poles[i].Y * Cross(normal, radial));
            const double weight = IsRational(profile) ? profile.Weights[j] : 1.0;
            surface.Weights.push_back(arc.Weights[i] * weight);
        }
    }
    return CheckedSurface(std::move(surface));
}

} // namespace Stereon
