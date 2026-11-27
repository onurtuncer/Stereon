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
#include <charconv>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <inplace_vector>
#include <limits>
#include <memory>
#include <numbers>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

module Stereon.Intersect;

import :Dispatch;
import Stereon.Robust;

#include "Analytic.hpp"
#include "BoxClip.hpp"
#include "Exact.hpp"
#include "Implicit.hpp"
#include "OvalPcurve.hpp"

/// The pairs of planes, cylinders, cones, spheres and tori that
/// PlanePlane.cpp, PlaneCylinder.cpp and the cylinder handlers do not take
/// (ADR-0010's analytic table): the closed forms first, then the
/// parametric-implicit trace (ImplicitTrace.cpp) for every configuration
/// that has none.
///
/// - Coaxial surfaces of revolution: cylinders, cones and tori sharing an
///   axis, a sphere centred on it, a plane across it. Each is the turn of a
///   meridian profile about the axis, lines and circles in the half-plane
///   (rho, z), so the surfaces meet in the circles turned by the points
///   where the profiles meet, touch along those where the profiles touch,
///   and coincide where the profiles do (Hoffmann 1989, section 7.2). Every
///   circle is exact, and its pcurve on a surface whose axis is the common
///   one is a line of constant v, exact too. Two spheres always share an
///   axis, the line of their centres, and a plane and a sphere the normal
///   through the centre. A spheroid, an ellipsoid with two equal semi-axes,
///   turns its meridian ellipse about the third; where it meets another
///   profile is isolated by certified subdivision of the ellipse's angle
///   (MeetEllipse), and its circles have line pcurves when its equal
///   semi-axes are A and B, its chart turning about its own axis.
/// - A plane through a cone's apex: two generators, one tangent generator,
///   or the apex alone.
/// - A plane containing a torus's axis: two meridian circles.
///
/// Whether axes coincide is decided exactly on the input doubles (RS-0001
/// G5). The meetings of the profiles are decided on their nominal values:
/// a cone's profile comes from the sine and cosine of its angle, so no
/// decision on it can be exact; one within the pair tolerance of a tangency
/// is declared a tangency and logged, as a coincidence within it is
/// declared (ADR-0029 section 6).

namespace Stereon::Detail {

namespace {

constexpr double Pi = std::numbers::pi;

// --- Spheroids ----------------------------------------------------------------

/// An ellipsoid with two equal semi-axes, a spheroid: the turn of its
/// meridian ellipse about the frame's axis `Along` (0 its X axis, 1 its Y
/// axis, 2 its own axis), with semi-axis `Across` across that axis and
/// `Height` along it. Equal exactly, on the input doubles: semi-axes equal
/// only within rounding make no surface of revolution. Nothing for an
/// ellipsoid with three unequal semi-axes; one with three equal ones turns
/// about every line through its centre, as a sphere does (IsSphere).
struct Spheroid {
    int Along;
    double Across;
    double Height;
};

/// The spheroid an ellipsoid is, from which semi-axes are equal.
[[nodiscard]] std::optional<Spheroid> SpheroidOf(const Ellipsoid& e)
{
    if (e.A == e.B) {
        return Spheroid{2, e.A, e.C};
    }
    if (e.B == e.C) {
        return Spheroid{0, e.B, e.A};
    }
    if (e.A == e.C) {
        return Spheroid{1, e.A, e.B};
    }
    return std::nullopt;
}

// --- The common axis --------------------------------------------------------

/// The direction of a surface's axis, exactly: its frame's cross product,
/// turned by the handedness flag; a plane's normal; a spheroid's axis of
/// revolution, its frame's X or Y axis when the equal semi-axes are the
/// other two.
template <typename G>
[[nodiscard]] ExactVector<4> ExactAxisOf(const G& g)
{
    if constexpr (std::same_as<G, Ellipsoid>) {
        const std::optional<Spheroid> s = SpheroidOf(g);
        if (s && s->Along != 2) {
            const ExactVector<1> along = Exactly(s->Along == 0 ? g.XAxis : g.YAxis);
            ExactVector<4> widened;
            widened.Hull = along.Hull;
            for (std::size_t k = 0; k < 3; ++k) {
                widened.Value[k] = [d = along.Value[k]] { return Expansion<4>(d()); };
            }
            return widened;
        }
    }
    const ExactVector<4> cross = ExactCross(Exactly(g.XAxis), Exactly(g.YAxis));
    if constexpr (std::same_as<G, Plane>) {
        return cross;
    } else {
        return g.RightHanded ? cross : ExactNegated(cross);
    }
}

[[nodiscard]] ExactVector<4> ExactAxisOf(const Elementary& e)
{
    return std::visit([](const auto& g) { return ExactAxisOf(g); }, e);
}

/// The axis's direction in doubles, as ExactAxisOf takes it.
[[nodiscard]] Vec3 AxisDirectionOf(const Elementary& e)
{
    if (const Ellipsoid* g = std::get_if<Ellipsoid>(&e)) {
        const std::optional<Spheroid> s = SpheroidOf(*g);
        if (s && s->Along != 2) {
            return Normalized(s->Along == 0 ? g->XAxis : g->YAxis);
        }
    }
    return Normalized(QuadricOf(e).Axis);
}

/// The point the frame is placed at: origin or centre.
[[nodiscard]] Point3 OriginOf(const Elementary& e)
{
    return QuadricOf(e).Origin;
}

[[nodiscard]] bool IsPlane(const Elementary& e)
{
    return std::holds_alternative<Plane>(e);
}

/// A sphere, or an ellipsoid with three equal semi-axes: round about every
/// line through its centre.
[[nodiscard]] bool IsSphere(const Elementary& e)
{
    if (const Ellipsoid* g = std::get_if<Ellipsoid>(&e)) {
        return g->A == g->B && g->B == g->C;
    }
    return std::holds_alternative<Sphere>(e);
}

/// An ellipsoid with three unequal semi-axes: no surface of revolution.
[[nodiscard]] bool IsTriaxial(const Elementary& e)
{
    const Ellipsoid* g = std::get_if<Ellipsoid>(&e);
    return g != nullptr && !SpheroidOf(*g);
}

/// A cylinder, cone, torus or spheroid: a surface with an axis of its own.
[[nodiscard]] bool HasOwnAxis(const Elementary& e)
{
    return !IsPlane(e) && !IsSphere(e) && !IsTriaxial(e);
}

[[nodiscard]] bool Parallel(const ExactVector<4>& a, const ExactVector<4>& b)
{
    return IsExactlyZero(ExactCross(a, b));
}

/// Whether point p lies on the line through o along a, exactly.
[[nodiscard]] bool OnLine(const Point3& p, const Point3& o, const ExactVector<4>& a)
{
    return IsExactlyZero(ExactCross(ExactOffset(p, o), a));
}

/// The line both surfaces turn about, if they share one: a point of it and
/// its unit direction, the first surface with an axis of its own giving
/// both, else the plane's normal through the sphere's centre, else the
/// line of the spheres' centres.
///
/// Misalignment is zero when the axes are one exactly. A pair placed by a
/// rounded rotation is coaxial in its construction but not in its doubles:
/// the cross products miss zero by units in the last place. It is taken as
/// coaxial when the misalignment moves its circles by at most the slack,
/// a sixteenth of the resolution, so that they still certify within the
/// resolution against both surfaces as given (Keep proves it, or the pair
/// fails closed); Misalignment then bounds that move.
struct Axis3 {
    Point3 Origin;
    Vec3 Direction;
    double Misalignment = 0.0;
};

/// How far a surface's points of interest reach from its frame's origin:
/// its radii and semi-axes; none for a plane, a cone's apex circle.
[[nodiscard]] double SizeOf(const Elementary& e)
{
    const Quadric q = QuadricOf(e);
    switch (q.Kind) {
    case QuadricKind::Plane:
        return 0.0;
    case QuadricKind::Torus:
        return q.Major + q.Radius;
    case QuadricKind::Ellipsoid:
        return std::max({q.Radius, q.SemiB, q.SemiC});
    default:
        return q.Radius;
    }
}

/// How far the other surface's axis is from r's, over the reach of the
/// pair: the sine between the axes times the reach, and the distance of
/// the other's origin from r's axis. A sphere has no axis of its own, a
/// plane no origin on one.
[[nodiscard]] double MisalignmentOf(const Elementary& r, const Elementary& other)
{
    const Vec3 a = AxisDirectionOf(r);
    const Vec3 offset = OriginOf(other) - OriginOf(r);
    const double reach = Length(offset) + SizeOf(r) + SizeOf(other);
    double m = 0.0;
    if (!IsSphere(other)) {
        // Axes apart by more than rounding are no common axis, however
        // small the reach: a cone's apex alone reaches nowhere.
        const double sine = Length(Cross(a, AxisDirectionOf(other)));
        if (!(sine <= 0x1p-40)) {
            return std::numeric_limits<double>::infinity();
        }
        m += sine * reach;
    }
    if (!IsPlane(other)) {
        m += Length(Cross(offset, a));
    }
    return m;
}

/// r's own axis, when the other surface shares it exactly or within the
/// slack (Axis3).
[[nodiscard]] std::optional<Axis3> OwnAxisShared(const Elementary& r, const Elementary& other,
                                                 const double slack)
{
    const ExactVector<4> axis = ExactAxisOf(r);
    const Point3 origin = OriginOf(r);
    const bool shared = IsPlane(other)    ? Parallel(axis, ExactAxisOf(other))
                        : IsSphere(other) ? OnLine(OriginOf(other), origin, axis)
                                          : Parallel(axis, ExactAxisOf(other)) &&
                                                OnLine(OriginOf(other), origin, axis);
    if (shared) {
        return Axis3{origin, AxisDirectionOf(r)};
    }
    const double misalignment = MisalignmentOf(r, other);
    if (!(misalignment <= slack)) {
        return std::nullopt;
    }
    return Axis3{origin, AxisDirectionOf(r), misalignment};
}

/// The axis both surfaces turn about, if they share one (Axis3).
[[nodiscard]] std::optional<Axis3> CommonAxis(const Elementary& a, const Elementary& b,
                                              const double slack)
{
    if (IsTriaxial(a) || IsTriaxial(b)) {
        return std::nullopt;
    }
    if (HasOwnAxis(a)) {
        return OwnAxisShared(a, b, slack);
    }
    if (HasOwnAxis(b)) {
        return OwnAxisShared(b, a, slack);
    }
    if (IsPlane(a) != IsPlane(b)) {
        const Elementary& sphere = IsPlane(a) ? b : a;
        const Elementary& plane = IsPlane(a) ? a : b;
        return Axis3{OriginOf(sphere), Normalized(QuadricOf(plane).Axis)};
    }
    // Two spheres; concentric ones about any line through their centre.
    const Vec3 between = OriginOf(b) - OriginOf(a);
    const bool concentric = LengthSquared(between) == 0.0;
    return Axis3{OriginOf(a), concentric ? Normalized(QuadricOf(a).Axis) : Normalized(between)};
}

/// Unit directions parallel within rounding, 2^-40 in sine, whose angle
/// moves a point `reach` away by at most the slack.
[[nodiscard]] bool WithinSlack(const Vec3& a, const Vec3& b, const double reach, const double slack)
{
    const double sine = Length(Cross(a, b));
    return sine <= 0x1p-40 && sine * reach <= slack;
}

/// Whether a sphere's poles lie on the common axis, so that its circles
/// about it are its parallels of latitude; whether a spheroid turns about
/// its own axis, its equal semi-axes A and B, so that its circles are too.
/// True for every other surface.
///
/// Poles within the slack of the axis, as for the axes themselves (Axis3):
/// the sine between the pole's direction and the axis's, times the radius.
[[nodiscard]] bool Aligned(const Elementary& e, const Elementary& other, const double slack)
{
    if (const Ellipsoid* g = std::get_if<Ellipsoid>(&e); g != nullptr && !IsSphere(e)) {
        return SpheroidOf(*g)->Along == 2;
    }
    if (!IsSphere(e)) {
        return true;
    }
    const bool concentric = IsSphere(other) && OriginOf(other) == OriginOf(e);
    const Vec3 pole = AxisDirectionOf(e);
    if (!IsSphere(other) || concentric) {
        return Parallel(ExactAxisOf(e), ExactAxisOf(other)) ||
               WithinSlack(pole, AxisDirectionOf(other), SizeOf(e), slack);
    }
    // Two spheres about the line of their centres.
    const ExactVector<2> between = ExactOffset(OriginOf(other), OriginOf(e));
    return IsExactlyZero(ExactCross(ExactAxisOf(e), between)) ||
           WithinSlack(pole, Normalized(OriginOf(other) - OriginOf(e)), SizeOf(e), slack);
}

// --- Meridian profiles --------------------------------------------------------

/// A line of the meridian half-plane, A rho + B z = C with (A, B) a unit
/// normal, and a circle of it.
struct ProfileLine {
    double A;
    double B;
    double C;
};

struct ProfileCircle {
    double Rho;
    double Z;
    double R;
};

/// A spheroid's meridian ellipse, centred on the axis at height Z, with
/// semi-axis Across across the axis and Height along it.
struct ProfileEllipse {
    double Z;
    double Across;
    double Height;
};

using ProfilePart = std::variant<ProfileLine, ProfileCircle, ProfileEllipse>;
using Profile = std::inplace_vector<ProfilePart, 2>;

/// The height of p along the axis.
[[nodiscard]] double HeightOf(const Axis3& axis, const Point3& p)
{
    return Dot(p - axis.Origin, axis.Direction);
}

/// The cone's two generators in the meridian plane, rho = +-k (z - z_apex)
/// with k = |tan a|, normalised: (|cos a|, -+|sin a|) . (rho, z) =
/// -+|sin a| z_apex.
[[nodiscard]] Profile ConeProfile(const Cone& c, const Axis3& axis)
{
    const double za = HeightOf(axis, Apex(c));
    const double ca = std::abs(Cos(c.SemiAngle));
    const double sa = std::abs(Sin(c.SemiAngle));
    return {ProfileLine{ca, -sa, -sa * za}, ProfileLine{ca, sa, sa * za}};
}

[[nodiscard]] Profile ProfileOf(const Elementary& e, const Axis3& axis)
{
    if (const Plane* p = std::get_if<Plane>(&e)) {
        return {ProfileLine{0.0, 1.0, HeightOf(axis, p->Origin)}};
    }
    if (const Cylinder* c = std::get_if<Cylinder>(&e)) {
        return {ProfileLine{1.0, 0.0, c->Radius}};
    }
    if (const Cone* k = std::get_if<Cone>(&e)) {
        return ConeProfile(*k, axis);
    }
    if (const Sphere* s = std::get_if<Sphere>(&e)) {
        return {ProfileCircle{0.0, HeightOf(axis, s->Center), s->Radius}};
    }
    if (const Ellipsoid* g = std::get_if<Ellipsoid>(&e)) {
        if (IsSphere(e)) {
            return {ProfileCircle{0.0, HeightOf(axis, g->Center), g->A}};
        }
        const Spheroid s = *SpheroidOf(*g);
        return {ProfileEllipse{HeightOf(axis, g->Center), s.Across, s.Height}};
    }
    const Torus& t = std::get<Torus>(e);
    return {ProfileCircle{t.MajorRadius, HeightOf(axis, t.Center), t.MinorRadius}};
}

/// Where two profiles meet: a point (rho, z), and whether they touch there
/// rather than cross.
struct Meet {
    double Rho;
    double Z;
    bool Tangent;
};

/// The profiles' meetings, and whether they coincide: exactly, Gap zero,
/// or within the tolerance, Gap then how far apart the parts declared
/// coincident are (ADR-0029 section 6).
struct Meets {
    std::vector<Meet> Points;
    bool Coincident = false;
    double Gap = 0.0;
    /// A meeting with a spheroid's ellipse left undecided (MeetEllipse).
    bool Unresolved = false;
};

/// Records two parts `gap` apart as coincident.
void DeclareCoincident(const double gap, Meets& out)
{
    out.Coincident = true;
    out.Gap = std::max(out.Gap, gap);
}

void MeetLines(const ProfileLine& p, const ProfileLine& q, const double tol, Meets& out)
{
    const double det = p.A * q.B - p.B * q.A;
    if (std::abs(det) <= 0x1p-40) {
        // Parallel: the same line within the tolerance, or none.
        const double sense = p.A * q.A + p.B * q.B > 0.0 ? 1.0 : -1.0;
        if (const double gap = std::abs(p.C - sense * q.C); gap <= tol) {
            DeclareCoincident(gap, out);
        }
        return;
    }
    out.Points.push_back({(p.C * q.B - p.B * q.C) / det, (p.A * q.C - p.C * q.A) / det, false});
}

void MeetLineCircle(const ProfileLine& l, const ProfileCircle& c, const double tol, Meets& out)
{
    const double d = l.A * c.Rho + l.B * c.Z - l.C; // The centre's signed distance.
    const double gap = std::abs(d) - c.R;
    if (gap > tol) {
        return;
    }
    const double footRho = c.Rho - d * l.A;
    const double footZ = c.Z - d * l.B;
    if (gap >= -tol) {
        out.Points.push_back({footRho, footZ, true});
        return;
    }
    const double half = std::sqrt(c.R * c.R - d * d);
    out.Points.push_back({footRho - half * l.B, footZ + half * l.A, false});
    out.Points.push_back({footRho + half * l.B, footZ - half * l.A, false});
}

/// How far along the centre line from p, d away from q, two circles touch
/// within the tolerance: the midpoint of their nearest points there, which
/// lies within half the gap of both. Externally, p's point is R_p along and
/// q's d - R_q; internally the larger one's far side meets the smaller's.
/// The radical axis's (d^2 + R_p^2 - R_q^2) / 2d is the same point only at
/// an exact tangency: near an internal one it lies gap (R_p + R_q) / 2d
/// from it, for radii 2 and 1.5 0.005 apart 0.0177, outside both spheres
/// by more than the tolerance (the OCCT corpus run of #326).
[[nodiscard]] double TouchAlong(const ProfileCircle& p, const ProfileCircle& q, const double d,
                                const double outer, const double inner)
{
    if (std::abs(outer) <= std::abs(inner)) {
        return 0.5 * (p.R + d - q.R);
    }
    return 0.5 * (d + (p.R >= q.R ? p.R + q.R : -(p.R + q.R)));
}

void MeetCircles(const ProfileCircle& p, const ProfileCircle& q, const double tol, Meets& out)
{
    const double dr = q.Rho - p.Rho;
    const double dz = q.Z - p.Z;
    const double d = std::hypot(dr, dz);
    if (d <= tol) {
        // Concentric within the tolerance: a point of one circle lies at
        // most d + |R_p - R_q| from the other.
        if (std::abs(p.R - q.R) <= tol) {
            DeclareCoincident(d + std::abs(p.R - q.R), out);
        }
        return;
    }
    const double outer = d - (p.R + q.R);
    const double inner = d - std::abs(p.R - q.R);
    if (outer > tol || inner < -tol) {
        return;
    }
    if (outer >= -tol || inner <= tol) {
        const double along = TouchAlong(p, q, d, outer, inner);
        out.Points.push_back({p.Rho + along * dr / d, p.Z + along * dz / d, true});
        return;
    }
    // Along the centre line from p, a; across it, h.
    const double a = (d * d + p.R * p.R - q.R * q.R) / (2.0 * d);
    const double baseRho = p.Rho + a * dr / d;
    const double baseZ = p.Z + a * dz / d;
    const double h = std::sqrt(std::max(0.0, p.R * p.R - a * a));
    out.Points.push_back({baseRho - h * dz / d, baseZ + h * dr / d, false});
    out.Points.push_back({baseRho + h * dz / d, baseZ - h * dr / d, false});
}

// --- A meridian ellipse against a profile part ----------------------------------

/// A profile part as a level function of the half-plane, w the height
/// above the ellipse's centre: G = Alpha (rho - Rho0)^2 + Beta (w - W0)^2 +
/// A rho + B w - C, zero on the part. A point where G = g lies within
/// |g| Reach of the part: a line has Alpha = Beta = 0 and (A, B) a unit
/// normal, so |g| is its distance; a circle Alpha = Beta = 1 and C its
/// radius R squared, and |p - q| - R = g / (|p - q| + R); an ellipse
/// centred on the axis Alpha = 1 / a^2, Beta = 1 / c^2 and C = 1, and the
/// point scaled from the centre onto it, by sqrt(1 + g), moves |g| / (1 +
/// sqrt(1 + g)) times its distance from the centre, which is at most
/// max(a, c) on the ellipse.
struct Level {
    double Alpha = 0.0;
    double Beta = 0.0;
    double Rho0 = 0.0;
    double W0 = 0.0;
    double A = 0.0;
    double B = 0.0;
    double C = 0.0;
    double Reach = 1.0;
};

/// The part's level function, heights measured from z0.
[[nodiscard]] Level LevelOf(const ProfilePart& part, const double z0)
{
    if (const ProfileLine* l = std::get_if<ProfileLine>(&part)) {
        return {0.0, 0.0, 0.0, 0.0, l->A, l->B, l->C - l->B * z0, 1.0};
    }
    if (const ProfileCircle* c = std::get_if<ProfileCircle>(&part)) {
        return {1.0, 1.0, c->Rho, c->Z - z0, 0.0, 0.0, c->R * c->R, 1.0 / c->R};
    }
    const ProfileEllipse& e = std::get<ProfileEllipse>(part);
    return {1.0 / (e.Across * e.Across), 1.0 / (e.Height * e.Height), 0.0, e.Z - z0, 0.0, 0.0, 1.0,
            std::max(e.Across, e.Height)};
}

/// g(t) = G(a cos t, c sin t) along the ellipse and its first two
/// derivatives, by the chain rule: the point's second derivative is minus
/// the point, measured from the centre.
template <typename T>
struct AlongJet {
    T Value;
    T Slope;
    T Curvature;
};

/// g, g' and g'' at angle t, in doubles or over an interval of angles.
template <typename T>
[[nodiscard]] AlongJet<T> AlongEllipse(const ProfileEllipse& e, const Level& g, const T& t)
{
    const T cosT = Cos(t);
    const T sinT = Sin(t);
    const T rho = T{e.Across} * cosT;
    const T w = T{e.Height} * sinT;
    const T dRho = -(T{e.Across} * sinT);
    const T dW = T{e.Height} * cosT;
    const T offRho = rho - T{g.Rho0};
    const T offW = w - T{g.W0};
    const T gRho = T{2.0 * g.Alpha} * offRho + T{g.A};
    const T gW = T{2.0 * g.Beta} * offW + T{g.B};
    return {T{g.Alpha} * offRho * offRho + T{g.Beta} * offW * offW + T{g.A} * rho + T{g.B} * w -
                T{g.C},
            gRho * dRho + gW * dW,
            T{2.0 * g.Alpha} * dRho * dRho + T{2.0 * g.Beta} * dW * dW - gRho * rho - gW * w};
}

/// The jet enclosed over a range of angles, or at one.
[[nodiscard]] AlongJet<Interval> JetAlong(const ProfileEllipse& e, const Level& g,
                                          const Interval& t)
{
    const AlongJet<I> j = AlongEllipse(e, g, I{t});
    return {j.Value.Value, j.Slope.Value, j.Curvature.Value};
}

/// Bisection to the last bit on the rounded sign of `f` between angles
/// whose signs differ, `low` the sign at the lower one.
template <typename F>
[[nodiscard]] double SignChange(const F& f, double lo, double hi, const int low)
{
    for (int step = 0; step < 80; ++step) {
        const double m = 0.5 * (lo + hi);
        if (!(lo < m && m < hi)) {
            break;
        }
        ((f(m) < 0.0) == (low < 0) ? lo : hi) = m;
    }
    return 0.5 * (lo + hi);
}

/// The ellipse's point at angle t, as a meeting of the profiles.
void AddEllipsePoint(const ProfileEllipse& e, const double t, const bool tangent, Meets& out)
{
    out.Points.push_back({e.Across * std::cos(t), e.Z + e.Height * std::sin(t), tangent});
}

/// An angle where the value's and the slope's signs are both certain.
struct Mark {
    double T;
    int Value;
    int Slope;
};

/// The mark at t; nothing where either sign is uncertain.
[[nodiscard]] std::optional<Mark> MarkAt(const ProfileEllipse& e, const Level& g, const double t)
{
    const AlongJet<Interval> j = JetAlong(e, g, Interval{t});
    const Mark m{t, CertainSign(j.Value), CertainSign(j.Slope)};
    if (m.Value == 0 || m.Slope == 0) {
        return std::nullopt;
    }
    return m;
}

/// One root of g between two marks whose values' signs differ, g monotone
/// between them.
void AddRootBetween(const ProfileEllipse& e, const Level& g, const Mark& lo, const Mark& hi,
                    Meets& out)
{
    const auto value = [&](const double t) { return AlongEllipse(e, g, t).Value; };
    AddEllipsePoint(e, SignChange(value, lo.T, hi.T, lo.Value), false, out);
}

/// A span over which g is strictly convex or concave, its marks' slopes of
/// opposite signs: one extremum t*, placed by bisection on the slope. Where
/// g(t*) lies within the tolerance of zero, measured as a distance
/// (Level), the profiles touch there; else they cross on each side of t*
/// whose end has g's other sign, or miss.
[[nodiscard]] bool AddExtremum(const ProfileEllipse& e, const Level& g, const Mark& lo,
                               const Mark& hi, const double tol, Meets& out)
{
    const auto slope = [&](const double t) { return AlongEllipse(e, g, t).Slope; };
    const double at = SignChange(slope, lo.T, hi.T, lo.Slope);
    const Interval value = JetAlong(e, g, Interval{at}).Value;
    const double reach = std::max(std::abs(value.GetLower()), std::abs(value.GetUpper())) * g.Reach;
    if (reach <= tol) {
        AddEllipsePoint(e, at, true, out);
        return true;
    }
    const int sign = CertainSign(value);
    if (sign == 0) {
        return false;
    }
    const auto valueAt = [&](const double t) { return AlongEllipse(e, g, t).Value; };
    if (lo.Value != sign) {
        AddEllipsePoint(e, SignChange(valueAt, lo.T, at, lo.Value), false, out);
    }
    if (hi.Value != sign) {
        AddEllipsePoint(e, SignChange(valueAt, at, hi.T, sign), false, out);
    }
    return true;
}

/// A cut of the span [lo, hi] where the marks are certain: a little off
/// its middle, so as to miss the symmetric points of symmetric profiles,
/// or further off where that one falls within rounding of a root of g or
/// of its slope.
[[nodiscard]] std::optional<Mark> CutOf(const ProfileEllipse& e, const Level& g, const Mark& lo,
                                        const Mark& hi)
{
    for (const double ratio : {0.5 + 0x1p-6, 0.5 - 0x1p-5, 0.5 + 0x1p-4, 0.5 - 0x1p-3}) {
        const double t = lo.T + ratio * (hi.T - lo.T);
        if (!(lo.T < t && t < hi.T)) {
            return std::nullopt;
        }
        if (const std::optional<Mark> m = MarkAt(e, g, t)) {
            return m;
        }
    }
    return std::nullopt;
}

/// The spans examined at most: a work budget, not a geometric tolerance.
constexpr std::size_t EllipseBudget = std::size_t{1} << 14;

/// One span decided: g keeps a sign over it, clear of the band of a touch
/// within the tolerance, no meeting; g is monotone,
/// where g' keeps its sign or g'' does with g' alike at both ends, and the
/// profiles cross exactly where g's ends differ; or g'' keeps its sign and
/// g' changes it (AddExtremum). False when the span needs splitting.
[[nodiscard]] bool ResolveSpan(const ProfileEllipse& e, const Level& g, const Mark& lo,
                               const Mark& hi, const double tol, Meets& out)
{
    const AlongJet<Interval> over = JetAlong(e, g, Interval{lo.T, hi.T});
    // Proved apart by more than the band within which a near miss is a
    // touch declared within the tolerance (AddExtremum).
    const double band = tol / g.Reach;
    if (!Intersects(over.Value, Interval{-band, band})) {
        return true;
    }
    const bool bent = !Contains(over.Curvature, 0.0);
    if (!Contains(over.Slope, 0.0) || (bent && lo.Slope == hi.Slope)) {
        if (lo.Value != hi.Value) {
            AddRootBetween(e, g, lo, hi, out);
        }
        return true;
    }
    return bent && AddExtremum(e, g, lo, hi, tol, out);
}

/// The whole turn of the ellipse's angle, from a start where the signs of
/// g and g' are certain, to the same point a period on.
[[nodiscard]] std::optional<std::pair<Mark, Mark>> WholeTurn(const ProfileEllipse& e,
                                                             const Level& g)
{
    for (int k = 0; k < 8; ++k) {
        const std::optional<Mark> start = MarkAt(e, g, -Pi + 0.375 * k);
        const std::optional<Mark> end =
            start ? MarkAt(e, g, start->T + CirclePeriod) : std::optional<Mark>{};
        if (end && end->Value == start->Value) {
            return std::pair{*start, *end};
        }
    }
    return std::nullopt;
}

/// Where the spheroid's meridian ellipse meets a profile part: certified
/// subdivision of the ellipse's angle over a whole turn, each span with the
/// signs of g and g' certain at its ends (ResolveSpan). False when a span
/// cannot be resolved within the budget: a contact of higher order, or the
/// profiles within rounding of each other over a stretch.
[[nodiscard]] bool MeetEllipse(const ProfileEllipse& e, const ProfilePart& other, const double tol,
                               Meets& out)
{
    const Level g = LevelOf(other, e.Z);
    const std::optional<std::pair<Mark, Mark>> turn = WholeTurn(e, g);
    if (!turn) {
        return false;
    }
    std::vector<std::pair<Mark, Mark>> work{*turn};
    for (std::size_t evaluated = 0; !work.empty(); ++evaluated) {
        if (evaluated >= EllipseBudget) {
            return false;
        }
        const auto [lo, hi] = work.back();
        work.pop_back();
        if (ResolveSpan(e, g, lo, hi, tol, out)) {
            continue;
        }
        const std::optional<Mark> cut = CutOf(e, g, lo, hi);
        if (!cut) {
            return false;
        }
        work.emplace_back(*cut, hi);
        work.emplace_back(lo, *cut);
    }
    return true;
}

/// Two meridian ellipses, or an ellipse and a circle about the axis, are
/// declared coincident when every point of one lies within the tolerance
/// of the other: their centres' and semi-axes' differences add up to at
/// most the tolerance.
[[nodiscard]] std::optional<double> EllipseGap(const ProfileEllipse& e, const ProfilePart& other)
{
    if (const ProfileEllipse* f = std::get_if<ProfileEllipse>(&other)) {
        return std::abs(e.Z - f->Z) + std::abs(e.Across - f->Across) +
               std::abs(e.Height - f->Height);
    }
    if (const ProfileCircle* c = std::get_if<ProfileCircle>(&other);
        c != nullptr && c->Rho == 0.0) {
        return std::abs(e.Z - c->Z) + std::abs(e.Across - c->R) + std::abs(e.Height - c->R);
    }
    return std::nullopt;
}

/// An ellipse and a part: coincident within the tolerance, or met.
void MeetEllipseParts(const ProfileEllipse& e, const ProfilePart& other, const double tol,
                      Meets& out)
{
    if (const std::optional<double> gap = EllipseGap(e, other); gap && *gap <= tol) {
        DeclareCoincident(*gap, out);
        return;
    }
    if (!MeetEllipse(e, other, tol, out)) {
        out.Unresolved = true;
    }
}

void MeetParts(const ProfilePart& p, const ProfilePart& q, const double tol, Meets& out)
{
    if (const ProfileEllipse* pe = std::get_if<ProfileEllipse>(&p)) {
        MeetEllipseParts(*pe, q, tol, out);
        return;
    }
    if (const ProfileEllipse* qe = std::get_if<ProfileEllipse>(&q)) {
        MeetEllipseParts(*qe, p, tol, out);
        return;
    }
    const ProfileLine* pl = std::get_if<ProfileLine>(&p);
    const ProfileLine* ql = std::get_if<ProfileLine>(&q);
    if (pl != nullptr && ql != nullptr) {
        MeetLines(*pl, *ql, tol, out);
    } else if (pl != nullptr) {
        MeetLineCircle(*pl, std::get<ProfileCircle>(q), tol, out);
    } else if (ql != nullptr) {
        MeetLineCircle(*ql, std::get<ProfileCircle>(p), tol, out);
    } else {
        MeetCircles(std::get<ProfileCircle>(p), std::get<ProfileCircle>(q), tol, out);
    }
}

/// The meetings in the half-plane rho >= 0, each once, in order of z then
/// rho (ADR-0029 G1); those within the tolerance of the axis are put on it.
[[nodiscard]] Meets MeetProfiles(const Profile& a, const Profile& b, const double tol)
{
    Meets all;
    for (const ProfilePart& p : a) {
        for (const ProfilePart& q : b) {
            MeetParts(p, q, tol, all);
        }
    }
    Meets kept{{}, all.Coincident, all.Gap, all.Unresolved};
    for (Meet m : all.Points) {
        if (m.Rho < -tol) {
            continue;
        }
        m.Rho = m.Rho <= tol ? 0.0 : m.Rho;
        const auto same = [&](const Meet& k) {
            return std::hypot(k.Rho - m.Rho, k.Z - m.Z) <= tol;
        };
        if (std::ranges::none_of(kept.Points, same)) {
            kept.Points.push_back(m);
        }
    }
    std::ranges::sort(kept.Points, {}, [](const Meet& m) { return std::pair{m.Z, m.Rho}; });
    return kept;
}

// --- Angles and boxes ---------------------------------------------------------

/// Whether the angle lies in the range modulo the period.
[[nodiscard]] bool AngleIn(const double angle, const Interval& range)
{
    if (!range.IsBounded() || range.GetUpper() - range.GetLower() >= CirclePeriod) {
        return true;
    }
    const double k = std::floor((angle - range.GetLower()) / CirclePeriod);
    return angle - k * CirclePeriod <= range.GetUpper();
}

/// The parameters t in `piece` at which u0 + s t lies in the range modulo
/// the period: at most two spans, since the piece is at most a turn wide.
[[nodiscard]] std::vector<Interval> AngleClip(const Interval& piece, const double u0,
                                              const double s, const Interval& range)
{
    if (!range.IsBounded() || range.GetUpper() - range.GetLower() >= CirclePeriod) {
        return {piece};
    }
    const double at0 = u0 + s * piece.GetLower();
    const double at1 = u0 + s * piece.GetUpper();
    const double low = std::min(at0, at1);
    const double first = std::floor((low - range.GetUpper()) / CirclePeriod);
    std::vector<Interval> spans;
    for (double k = first; k <= first + 2.0; k += 1.0) {
        const Interval shifted{range.GetLower() + k * CirclePeriod,
                               range.GetUpper() + k * CirclePeriod};
        const double t0 = (shifted.GetLower() - u0) / s;
        const double t1 = (shifted.GetUpper() - u0) / s;
        const Interval span = Intersection(piece, Interval{std::min(t0, t1), std::max(t0, t1)});
        if (!span.IsEmpty() && span.GetUpper() > span.GetLower()) {
            spans.push_back(span);
        }
    }
    return spans;
}

/// The parameters of p on the surface, with whether they lie in the box.
[[nodiscard]] bool ParametersIn(const Elementary& e, const Point2& uv, const Box2& box)
{
    if (IsPlane(e)) {
        return Contains(box, uv);
    }
    const bool vAngle = std::holds_alternative<Torus>(e);
    return AngleIn(uv.X, box.X) && (vAngle ? AngleIn(uv.Y, box.Y) : Contains(box.Y, uv.Y));
}

// --- Circles about the axis -----------------------------------------------------

/// The frame of the circles: XAxis from the first surface's, made square
/// to the axis.
[[nodiscard]] std::pair<Vec3, Vec3> CircleFrame(const Axis3& axis, const Elementary& first)
{
    const Vec3 x = QuadricOf(first).XAxis;
    Vec3 across = x - Dot(x, axis.Direction) * axis.Direction;
    if (LengthSquared(across) <= 0x1p-60) {
        across =
            QuadricOf(first).YAxis - Dot(QuadricOf(first).YAxis, axis.Direction) * axis.Direction;
    }
    const Vec3 cx = Normalized(across);
    return {cx, Cross(axis.Direction, cx)};
}

/// The sense a surface's angle turns about the axis.
[[nodiscard]] double SenseOf(const Elementary& e, const Axis3& axis)
{
    const Quadric q = QuadricOf(e);
    return Dot(Cross(q.XAxis, q.YAxis), axis.Direction) > 0.0 ? 1.0 : -1.0;
}

/// A circle's pcurve on one surface over [0, 2 pi]: on a plane the exact
/// plane circle; on a surface of revolution about the axis, the line
/// (u0 + s t, v0), s the sense the surface's angle turns about the axis.
struct CirclePcurve {
    bool OnPlane;
    Circle2d Planar;
    double U0;
    double V0;
    double Sense;
};

[[nodiscard]] CirclePcurve PcurveOf(const Elementary& e, const Circle& c, const Axis3& axis)
{
    if (const Plane* p = std::get_if<Plane>(&e)) {
        const auto in = [&](const Vec3& v) { return Vec2{Dot(v, p->XAxis), Dot(v, p->YAxis)}; };
        return {true, Circle2d{Project(*p, c.Center), c.Radius, in(c.XAxis), in(c.YAxis)}, 0.0, 0.0,
                1.0};
    }
    const double sense = SenseOf(e, axis);
    Point2 uv = ProjectOnto(e, Eval(c, 0.0));
    // A start a hair below the seam belongs at the seam's other end.
    if (uv.X >= CirclePeriod - 0x1p-40) {
        uv.X = sense > 0.0 ? 0.0 : CirclePeriod;
    } else if (sense < 0.0 && uv.X <= 0x1p-40) {
        uv.X = CirclePeriod;
    }
    return {false, {}, uv.X, uv.Y, sense};
}

/// The parameter where a pcurve of a surface of revolution crosses its
/// chart's seam, inside (0, 2 pi); nothing on a plane or without one.
[[nodiscard]] std::optional<double> SeamOf(const CirclePcurve& p)
{
    if (p.OnPlane) {
        return std::nullopt;
    }
    const double t = p.Sense > 0.0 ? CirclePeriod - p.U0 : p.U0;
    if (t <= 0x1p-40 || t >= CirclePeriod - 0x1p-40) {
        return std::nullopt;
    }
    return t;
}

/// The pcurve over a piece of [0, 2 pi] that crosses no seam: its angle
/// shifted by the whole periods that bring the piece into the chart.
[[nodiscard]] Curve2d PieceCurve(const CirclePcurve& p, const Interval& piece)
{
    if (p.OnPlane) {
        return Curve2d(p.Planar);
    }
    const double mid = p.U0 + p.Sense * 0.5 * (piece.GetLower() + piece.GetUpper());
    const double shift = CirclePeriod * std::floor(mid / CirclePeriod);
    return Curve2d(Line2d{Point2{p.U0 - shift, p.V0}, Vec2{p.Sense, 0.0}});
}

/// The spans of a seam-free piece inside a surface's box: an angle clip for
/// a surface of revolution, after its height; a plane keeps the circle
/// whole when its box holds it, drops it when its box misses it or lies
/// inside it, and cannot cut it otherwise (as PlaneCylinder).
[[nodiscard]] Result<std::vector<Interval>> ClipPiece(const CirclePcurve& p, const Interval& piece,
                                                      const Box2& box)
{
    if (!p.OnPlane) {
        if (!Contains(box.Y, p.V0)) {
            return std::vector<Interval>{};
        }
        const double mid = p.U0 + p.Sense * 0.5 * (piece.GetLower() + piece.GetUpper());
        const double shift = CirclePeriod * std::floor(mid / CirclePeriod);
        return AngleClip(piece, p.U0 - shift, p.Sense, box.X);
    }
    const Circle2d& c = p.Planar;
    const Box2 disk{Around(c.Center.X, c.Radius), Around(c.Center.Y, c.Radius)};
    if (!Intersects(disk, box)) {
        return std::vector<Interval>{};
    }
    if (Contains(box, disk)) {
        return std::vector<Interval>{piece};
    }
    const auto inside = [&](const double x, const double y) {
        return std::hypot(x - c.Center.X, y - c.Center.Y) < 0.5 * c.Radius;
    };
    if (IsBounded(box) && inside(box.X.GetLower(), box.Y.GetLower()) &&
        inside(box.X.GetUpper(), box.Y.GetUpper()) && inside(box.X.GetLower(), box.Y.GetUpper()) &&
        inside(box.X.GetUpper(), box.Y.GetLower())) {
        return std::vector<Interval>{};
    }
    return Fail(ErrorCode::NotImplemented,
                "Intersect: a partial circle needs certified clipping against a plane's box");
}

/// The circle's seam-free pieces of [0, 2 pi], cut where either pcurve
/// crosses its chart's seam.
[[nodiscard]] std::vector<Interval> SeamPieces(const CirclePcurve& a, const CirclePcurve& b)
{
    std::vector<double> cuts{0.0, CirclePeriod};
    for (const CirclePcurve* p : {&a, &b}) {
        if (const std::optional<double> t = SeamOf(*p)) {
            cuts.push_back(*t);
        }
    }
    std::ranges::sort(cuts);
    std::vector<Interval> pieces;
    for (std::size_t i = 1; i < cuts.size(); ++i) {
        if (cuts[i] - cuts[i - 1] > 0x1p-40) {
            pieces.emplace_back(cuts[i - 1], cuts[i]);
        }
    }
    return pieces;
}

/// The spans of [0, 2 pi] of the circle inside both boxes, each within one
/// chart of both pcurves.
[[nodiscard]] Result<std::vector<Interval>> CircleSpans(const CirclePcurve& a, const Box2& boxA,
                                                        const CirclePcurve& b, const Box2& boxB)
{
    std::vector<Interval> spans;
    for (const Interval& piece : SeamPieces(a, b)) {
        STEREON_TRY(const std::vector<Interval> onA, ClipPiece(a, piece, boxA));
        for (const Interval& spanA : onA) {
            STEREON_TRY(const std::vector<Interval> onB, ClipPiece(b, spanA, boxB));
            spans.insert(spans.end(), onB.begin(), onB.end());
        }
    }
    return spans;
}

/// The surfaces as the handlers take them, with their requested boxes.
struct Pair {
    const Surface& ErasedA;
    const Surface& ErasedB;
    Elementary A;
    Elementary B;
    Box2 BoxA;
    Box2 BoxB;
    double Tolerance;
    /// How far from coaxial a pair may be and be taken as coaxial (Axis3).
    double Slack;
};

/// The circle the profiles meet at, as branches inside both boxes.
/// The circle's branch over a span, its contact the meeting's: a tangency
/// is declared within the pair tolerance (MeetProfiles), so certified
/// within it.
[[nodiscard]] Result<void> KeepCircle(SurfaceSurfaceIntersection& result, const Pair& p,
                                      const Circle& circle, const Meet& m, Curve2d curveA,
                                      Curve2d curveB, const Interval& span, const Context& ctx)
{
    IntersectionBranch branch =
        m.Tangent ? Tangency(Curve(circle), std::move(curveA), std::move(curveB), span)
                  : Crossing(Curve(circle), std::move(curveA), std::move(curveB), span);
    return Keep(result, std::move(branch), p.ErasedA, p.ErasedB, ctx,
                m.Tangent ? p.Tolerance : 0.0);
}

[[nodiscard]] bool WholeChart(const Elementary& e, const Box2& box);

// --- Circles on a chart that does not turn about the axis -------------------------

/// The chart of a sphere or an ellipsoid, as OvalPcurve.hpp takes it.
[[nodiscard]] OvalChart OvalChartOf(const Elementary& e)
{
    if (const Sphere* s = std::get_if<Sphere>(&e)) {
        return ChartOf(*s);
    }
    return ChartOf(std::get<Ellipsoid>(e));
}

/// One side of a circle whose pcurve may be fitted: a sphere whose poles
/// are off the common axis, or a spheroid turning about its frame's X or Y
/// axis, has the certified interpolant of OvalPcurve.hpp; every other
/// surface the exact pcurve of PcurveOf.
struct CircleSide {
    const Elementary* Geometry;
    const Surface* Erased;
    const Box2* Box;
    bool Fitted;
    CirclePcurve Line;
};

/// The side of the circle on one surface: fitted or exact.
[[nodiscard]] CircleSide SideOf(const Elementary& e, const Surface& erased, const Box2& box,
                                const Elementary& other, const Circle& circle, const Axis3& axis,
                                const double slack)
{
    const bool fitted = !Aligned(e, other, slack);
    return {&e, &erased, &box, fitted, fitted ? CirclePcurve{} : PcurveOf(e, circle, axis)};
}

/// Where the side's pcurve leaves its chart: a seam crossing of either kind.
void AddSeamCuts(const CircleSide& side, const Circle& circle, std::vector<double>& cuts)
{
    if (!side.Fitted) {
        if (const std::optional<double> t = SeamOf(side.Line)) {
            cuts.push_back(*t);
        }
        return;
    }
    for (const double t : OvalSeamCrossings(OvalChartOf(*side.Geometry), OvalOf(circle))) {
        cuts.push_back(t);
    }
}

/// The side's pcurve over a piece crossing no seam of either chart.
[[nodiscard]] Result<Curve2d> SideCurve(const CircleSide& side, const Circle& circle,
                                        const Interval& piece, const Context& ctx)
{
    if (!side.Fitted) {
        return PieceCurve(side.Line, piece);
    }
    return FittedOvalPcurve(OvalChartOf(*side.Geometry), OvalOf(circle), Curve(circle),
                            *side.Erased, piece, ctx, "Intersect");
}

/// The stretches of a span inside the side's box.
[[nodiscard]] Result<std::vector<Interval>> SideSpans(const CircleSide& side, const Curve2d& curve,
                                                      const Interval& span, const Context& ctx)
{
    if (!side.Fitted) {
        return ClipPiece(side.Line, span, *side.Box);
    }
    if (WholeChart(*side.Geometry, *side.Box)) {
        return std::vector<Interval>{span};
    }
    return ClipPcurveToBox(curve, span, *side.Box, ctx, "Intersect");
}

/// The circle cut at every seam of both charts, into pieces of [0, 2 pi].
[[nodiscard]] std::vector<Interval> FittedPieces(const CircleSide& a, const CircleSide& b,
                                                 const Circle& circle)
{
    std::vector<double> cuts{0.0, CirclePeriod};
    AddSeamCuts(a, circle, cuts);
    AddSeamCuts(b, circle, cuts);
    std::ranges::sort(cuts);
    std::vector<Interval> pieces;
    for (std::size_t i = 1; i < cuts.size(); ++i) {
        if (cuts[i] - cuts[i - 1] > 0x1p-40) {
            pieces.emplace_back(cuts[i - 1], cuts[i]);
        }
    }
    return pieces;
}

/// A circle about the axis with a fitted pcurve on a side whose chart does
/// not turn about it: cut at both charts' seams, each piece's pcurves
/// computed, then clipped to both boxes.
[[nodiscard]] Result<void> AddFittedCircle(SurfaceSurfaceIntersection& result, const Pair& p,
                                           const Circle& circle, const Axis3& axis, const Meet& m,
                                           const Context& ctx)
{
    const CircleSide a = SideOf(p.A, p.ErasedA, p.BoxA, p.B, circle, axis, p.Slack);
    const CircleSide b = SideOf(p.B, p.ErasedB, p.BoxB, p.A, circle, axis, p.Slack);
    for (const Interval& piece : FittedPieces(a, b, circle)) {
        STEREON_TRY(const Curve2d curveA, SideCurve(a, circle, piece, ctx));
        STEREON_TRY(const Curve2d curveB, SideCurve(b, circle, piece, ctx));
        STEREON_TRY(const std::vector<Interval> onA, SideSpans(a, curveA, piece, ctx));
        for (const Interval& spanA : onA) {
            STEREON_TRY(const std::vector<Interval> onB, SideSpans(b, curveB, spanA, ctx));
            for (const Interval& span : onB) {
                STEREON_TRY_VOID(KeepCircle(result, p, circle, m, curveA, curveB, span, ctx));
            }
        }
    }
    return {};
}

[[nodiscard]] Result<void> AddCircle(SurfaceSurfaceIntersection& result, const Pair& p,
                                     const Axis3& axis, const Meet& m, const Context& ctx)
{
    const auto [cx, cy] = CircleFrame(axis, HasOwnAxis(p.A) || !HasOwnAxis(p.B) ? p.A : p.B);
    const Circle circle{axis.Origin + m.Z * axis.Direction, m.Rho, cx, cy};
    if (!Aligned(p.A, p.B, p.Slack) || !Aligned(p.B, p.A, p.Slack)) {
        return AddFittedCircle(result, p, circle, axis, m, ctx);
    }
    const CirclePcurve onA = PcurveOf(p.A, circle, axis);
    const CirclePcurve onB = PcurveOf(p.B, circle, axis);
    STEREON_TRY(const std::vector<Interval> spans, CircleSpans(onA, p.BoxA, onB, p.BoxB));
    for (const Interval& span : spans) {
        STEREON_TRY_VOID(KeepCircle(result, p, circle, m, PieceCurve(onA, span),
                                    PieceCurve(onB, span), span, ctx));
    }
    return {};
}

/// A meeting on the axis: a point, where both surfaces reach the axis.
void AddAxisPoint(SurfaceSurfaceIntersection& result, const Pair& p, const Axis3& axis,
                  const Meet& m)
{
    const Point3 at = axis.Origin + m.Z * axis.Direction;
    const Point2 onA = ProjectOnto(p.A, at);
    const Point2 onB = ProjectOnto(p.B, at);
    if (ParametersIn(p.A, onA, p.BoxA) && ParametersIn(p.B, onB, p.BoxB)) {
        result.Points.push_back({at, onA, onB});
    }
}

/// A pair taken as coaxial within rounding says so (Axis3).
void LogMisalignment(const Axis3& axis, const Context& ctx)
{
    if (axis.Misalignment > 0.0) {
        ctx.Log(LogLevel::Info, "Intersect: surfaces " + Number(axis.Misalignment) +
                                    " mm from coaxial, within rounding: taken as coaxial, "
                                    "their circles certified against both (ADR-0029)");
    }
}

void LogTangency(const Meets& meets, const Context& ctx)
{
    if (std::ranges::any_of(meets.Points, &Meet::Tangent)) {
        ctx.Log(LogLevel::Info, "Intersect: coaxial surfaces whose meridian profiles touch within "
                                "the pair tolerance: declared tangent along the circle (ADR-0029)");
    }
}

// --- Coincident surfaces of revolution ---------------------------------------

/// Whether a box is a whole chart of the surface: every angle and every
/// height or latitude.
[[nodiscard]] bool WholeChart(const Elementary& e, const Box2& box)
{
    const auto turn = [](const Interval& i) {
        return !i.IsBounded() || i.GetUpper() - i.GetLower() >= CirclePeriod;
    };
    if (std::holds_alternative<Torus>(e)) {
        return turn(box.X) && turn(box.Y);
    }
    if (IsSphere(e) || std::holds_alternative<Ellipsoid>(e)) {
        return turn(box.X) && Contains(box.Y, Interval{-SpherePoleLatitude, SpherePoleLatitude});
    }
    return turn(box.X) && box.Y == Interval::Everything();
}

/// A point of the chart away from the poles and the apex.
[[nodiscard]] Point2 RegularPoint(const Elementary& e)
{
    if (const Cone* c = std::get_if<Cone>(&e)) {
        return {0.3, 1.0 - c->Radius / Sin(c->SemiAngle)};
    }
    return {0.3, 0.3};
}

/// The normal Su x Sv at a point.
[[nodiscard]] Vec3 NormalAt(const Elementary& e, const Point2& uv)
{
    const Jet<double> jet = JetAt(e, uv);
    return Cross(jet.Su, jet.Sv);
}

/// Coincident surfaces of revolution, over whole charts: the overlap whose
/// map turns A's angle into B's and carries A's v onto B's. Exactly coaxial
/// and alike, so the map is affine; the sphere's only where its poles lie
/// on the other's axis. A gap above zero is a coincidence declared within
/// the pair tolerance: Tolerant evidence, and the gap as its distance.
[[nodiscard]] Result<SurfaceSurfaceIntersection> Coincident(const Pair& p, const Axis3& axis,
                                                            const double gap, const Context& ctx)
{
    if (!WholeChart(p.A, p.BoxA) || !WholeChart(p.B, p.BoxB) || !Aligned(p.A, p.B, p.Slack) ||
        !Aligned(p.B, p.A, p.Slack)) {
        return Fail(ErrorCode::NotImplemented,
                    "Intersect: coincident cones, spheres or tori over part of their charts, or "
                    "spheres whose poles differ");
    }
    const Point2 at = RegularPoint(p.A);
    const Point2 start = ProjectOnto(p.B, PointAt(p.A, Point2{0.0, at.Y}));
    const Point2 along = ProjectOnto(p.B, PointAt(p.A, Point2{0.0, at.Y + 0.25}));
    OverlapRegion overlap;
    overlap.OnA = p.BoxA;
    overlap.OnB = p.BoxB;
    overlap.AToB.U = Vec2{SenseOf(p.A, axis) * SenseOf(p.B, axis), 0.0};
    overlap.AToB.V = Vec2{0.0, along.Y >= start.Y ? 1.0 : -1.0};
    overlap.AToB.Offset = Vec2{start.X, start.Y - overlap.AToB.V.Y * at.Y};
    overlap.SameOrientation =
        Dot(NormalAt(p.A, at), NormalAt(p.B, ProjectOnto(p.B, PointAt(p.A, at)))) > 0.0;
    overlap.FullDomain = true;
    if (gap > 0.0) {
        // Declared within the pair tolerance, not exact: said so, with the
        // distance, as for parallel planes and cylinders.
        overlap.Evidence = CoincidenceEvidence::Tolerant;
        overlap.Distance = gap;
        ctx.Log(LogLevel::Info, "Intersect: coaxial surfaces of revolution " + Number(gap) +
                                    " mm apart, within the pair tolerance " + Number(p.Tolerance) +
                                    " mm: declared coincident, nothing moved (ADR-0029)");
    } else {
        ctx.Log(LogLevel::Info, "Intersect: coincident surfaces of revolution (ADR-0029)");
    }
    SurfaceSurfaceIntersection result;
    result.Overlaps.push_back(overlap);
    return result;
}

[[nodiscard]] bool UnresolvedConeAngles(const Pair& p)
{
    const auto* a = std::get_if<Cone>(&p.A);
    const auto* b = std::get_if<Cone>(&p.B);
    if (!a || !b) {
        return false;
    }
    const double difference = std::abs(a->SemiAngle) - std::abs(b->SemiAngle);
    return difference != 0.0 && std::abs(std::sin(difference)) <= 0x1p-40;
}

/// The coaxial closed form, when the surfaces share an axis and their
/// circles can be given exactly; nothing otherwise, for the trace.
[[nodiscard]] Result<std::optional<SurfaceSurfaceIntersection>> Coaxial(const Pair& p,
                                                                        const Context& ctx)
{
    const std::optional<Axis3> axis = CommonAxis(p.A, p.B, p.Slack);
    if (!axis) {
        return std::optional<SurfaceSurfaceIntersection>{};
    }
    if (UnresolvedConeAngles(p)) {
        return Fail(ErrorCode::IntersectionFailed,
                    "Nearly parallel unequal cone profiles require refinement");
    }
    const Meets meets = MeetProfiles(ProfileOf(p.A, *axis), ProfileOf(p.B, *axis), p.Tolerance);
    if (meets.Unresolved) {
        return Fail(ErrorCode::IntersectionFailed,
                    "Intersect: a coaxial spheroid's meridian ellipse meets the other profile in "
                    "a contact of higher order, or runs within rounding of it, which the "
                    "subdivision of its angle could not resolve (ADR-0029)");
    }
    if (meets.Coincident) {
        STEREON_TRY(SurfaceSurfaceIntersection overlap,
                    Coincident(p, *axis, std::max(meets.Gap, axis->Misalignment), ctx));
        return std::optional<SurfaceSurfaceIntersection>{std::move(overlap)};
    }
    LogTangency(meets, ctx);
    LogMisalignment(*axis, ctx);
    SurfaceSurfaceIntersection result;
    for (const Meet& m : meets.Points) {
        if (m.Rho > 0.0) {
            STEREON_TRY_VOID(AddCircle(result, p, *axis, m, ctx));
        } else {
            AddAxisPoint(result, p, *axis, m);
        }
    }
    return std::optional<SurfaceSurfaceIntersection>{std::move(result)};
}

// --- A sphere touching a surface off its axis ---------------------------------------

/// The nearest points of a surface of revolution to a point off its axis,
/// with their distance: by symmetry they lie in the point's meridian plane,
/// where the surface is its profile, so the nearest points of the profile's
/// lines or circle are the nearest of the surface. Two when the point is
/// as near two of them, within the tolerance.
struct Nearest {
    double Distance = std::numeric_limits<double>::infinity();
    std::inplace_vector<Point3, 2> Points;
};

/// The profile's nearest point to (rho, z) on one part, and its distance:
/// the foot on a line, the point of a circle along the ray from its centre.
[[nodiscard]] std::optional<std::pair<double, std::array<double, 2>>>
NearestOnPart(const ProfilePart& part, const double rho, const double z)
{
    if (const ProfileLine* l = std::get_if<ProfileLine>(&part)) {
        const double d = l->A * rho + l->B * z - l->C;
        return std::pair{std::abs(d), std::array{rho - d * l->A, z - d * l->B}};
    }
    const ProfileCircle& c = std::get<ProfileCircle>(part);
    const double reach = std::hypot(rho - c.Rho, z - c.Z);
    if (!(reach > 0.0)) {
        return std::nullopt; // The point at the circle's centre: every point is nearest.
    }
    const double s = c.R / reach;
    return std::pair{std::abs(reach - c.R),
                     std::array{c.Rho + s * (rho - c.Rho), c.Z + s * (z - c.Z)}};
}

/// The nearest points of a cylinder, cone or torus to the point, which lies
/// off its axis; nothing where the nearest points are not isolated.
[[nodiscard]] std::optional<Nearest> NearestPoints(const Elementary& e, const Point3& point,
                                                   const double tol)
{
    const Axis3 axis{OriginOf(e), Normalized(QuadricOf(e).Axis)};
    const double z = HeightOf(axis, point);
    const Vec3 radial = (point - axis.Origin) - z * axis.Direction;
    const double rho = Length(radial);
    if (!(rho > tol)) {
        // On the axis within the tolerance: the distance varies by at most
        // 2 rho around the circle, so a whole circle is nearest, not points.
        return std::nullopt;
    }
    const Vec3 out = radial / rho;
    Nearest nearest;
    for (const ProfilePart& part : ProfileOf(e, axis)) {
        const auto found = NearestOnPart(part, rho, z);
        if (!found) {
            return std::nullopt;
        }
        const auto& [distance, at] = *found;
        const Point3 p = axis.Origin + at[1] * axis.Direction + at[0] * out;
        if (distance < nearest.Distance - tol) {
            nearest.Distance = distance;
            nearest.Points = {p};
        } else if (distance <= nearest.Distance + tol) {
            nearest.Points.push_back(p);
        }
    }
    return nearest;
}

/// The points of contact inside both requested boxes.
void AddContacts(SurfaceSurfaceIntersection& result, const Pair& p,
                 const std::inplace_vector<Point3, 2>& points)
{
    for (const Point3& at : points) {
        const Point2 onA = ProjectOnto(p.A, at);
        const Point2 onB = ProjectOnto(p.B, at);
        if (ParametersIn(p.A, onA, p.BoxA) && ParametersIn(p.B, onB, p.BoxB)) {
            result.Points.push_back({at, onA, onB});
        }
    }
}

/// A sphere against a cylinder, cone or torus whose axis misses its centre:
/// a proved separation when the centre is further from the surface than
/// the radius and the tolerance, and the points of contact when it is the
/// radius away within the tolerance, since the open ball then misses the
/// surface. Nothing otherwise: the sphere crosses the surface, for the
/// trace.
[[nodiscard]] Result<std::optional<SurfaceSurfaceIntersection>> SphereTouch(const Pair& p,
                                                                            const Context& ctx)
{
    const bool aSphere = std::holds_alternative<Sphere>(p.A);
    const Elementary& other = aSphere ? p.B : p.A;
    if (aSphere == std::holds_alternative<Sphere>(p.B) || !HasOwnAxis(other) ||
        std::holds_alternative<Ellipsoid>(other)) {
        return std::optional<SurfaceSurfaceIntersection>{};
    }
    const Sphere& sphere = std::get<Sphere>(aSphere ? p.A : p.B);
    const std::optional<Nearest> nearest = NearestPoints(other, sphere.Center, p.Tolerance);
    if (!nearest || nearest->Distance < sphere.Radius - p.Tolerance) {
        return std::optional<SurfaceSurfaceIntersection>{};
    }
    SurfaceSurfaceIntersection result;
    if (nearest->Distance <= sphere.Radius + p.Tolerance) {
        ctx.Log(LogLevel::Info, "Intersect: a sphere touches a surface off its axis within the "
                                "pair tolerance: declared points of contact (ADR-0029)");
        AddContacts(result, p, nearest->Points);
    }
    return std::optional<SurfaceSurfaceIntersection>{std::move(result)};
}

// --- A plane through a cone's apex ---------------------------------------------

/// The generator of the cone at angle u as a line from the apex, with its
/// pcurve on each surface: on the cone the line (u, v_apex + t), since the
/// generators are unit-speed in v.
[[nodiscard]] Result<void> AddGenerator(SurfaceSurfaceIntersection& result, const Pair& p,
                                        const Plane& plane, const Cone& cone, const double u,
                                        const bool tangent, const Context& ctx)
{
    const Point3 apex = Apex(cone);
    const double vApex = -cone.Radius / Sin(cone.SemiAngle);
    const Vec3 along = Deriv(cone, Point2{u, vApex}, 0, 1);
    const Line2d onPlane{Project(plane, apex),
                         Vec2{Dot(along, plane.XAxis), Dot(along, plane.YAxis)}};
    const Line2d onCone{Point2{u, vApex}, Vec2{0.0, 1.0}};
    const bool planeFirst = IsPlane(p.A);
    const Box2& planeBox = planeFirst ? p.BoxA : p.BoxB;
    const Box2& coneBox = planeFirst ? p.BoxB : p.BoxA;
    if (!AngleIn(u, coneBox.X)) {
        return {};
    }
    const Interval coneRange = AxisRange(vApex, 1.0, coneBox.Y);
    const Interval range = Intersection(ClipLine(onPlane, planeBox), coneRange);
    if (range.IsEmpty()) {
        return {};
    }
    const Curve line{Line{apex, along}};
    Curve2d first = planeFirst ? Curve2d(onPlane) : Curve2d(onCone);
    Curve2d second = planeFirst ? Curve2d(onCone) : Curve2d(onPlane);
    IntersectionBranch branch = tangent
                                    ? Tangency(line, std::move(first), std::move(second), range)
                                    : Crossing(line, std::move(first), std::move(second), range);
    return Keep(result, std::move(branch), p.ErasedA, p.ErasedB, ctx);
}

/// The plane through the apex meets the cone in the generators whose
/// direction lies in it: n . (sin a e(u) + cos a A) = 0, that is
/// n_x cos u + n_y sin u = -cot a (n . A), with n_x and n_y the normal's
/// components along the cone's X and Y axes. Two solutions, a double one
/// within rounding of the tangent plane, or none but the apex.
[[nodiscard]] Result<SurfaceSurfaceIntersection> ApexSection(const Pair& p, const Plane& plane,
                                                             const Cone& cone, const Context& ctx)
{
    const Vec3 n = Normalized(Cross(plane.XAxis, plane.YAxis));
    const double nx = Dot(n, cone.XAxis);
    const double ny = Dot(n, cone.YAxis);
    const double across = std::hypot(nx, ny);
    const double target = -Cos(cone.SemiAngle) / Sin(cone.SemiAngle) * Dot(n, Axis(cone));
    SurfaceSurfaceIntersection result;
    const double phi = std::atan2(ny, nx);
    const double gap = across - std::abs(target);
    if (gap < -0x1p-40) {
        const Point3 apex = Apex(cone);
        const Point2 onA = ProjectOnto(p.A, apex);
        const Point2 onB = ProjectOnto(p.B, apex);
        if (ParametersIn(p.A, onA, p.BoxA) && ParametersIn(p.B, onB, p.BoxB)) {
            result.Points.push_back({apex, onA, onB});
        }
        return result;
    }
    const bool tangent = gap <= 0x1p-40;
    const double delta = tangent ? 0.0 : std::acos(std::clamp(target / across, -1.0, 1.0));
    std::array<double, 2> angles{WrappedAngle(Sin(phi - delta), Cos(phi - delta)),
                                 WrappedAngle(Sin(phi + delta), Cos(phi + delta))};
    std::ranges::sort(angles);
    for (std::size_t k = 0; k < (tangent ? 1U : 2U); ++k) {
        STEREON_TRY_VOID(AddGenerator(result, p, plane, cone, angles[k], tangent, ctx));
    }
    return result;
}

/// The plane and the cone when the plane passes through the apex within the
/// pair tolerance; nothing otherwise.
[[nodiscard]] Result<std::optional<SurfaceSurfaceIntersection>> PlaneThroughApex(const Pair& p,
                                                                                 const Context& ctx)
{
    const Plane* plane = std::get_if<Plane>(IsPlane(p.A) ? &p.A : &p.B);
    const Cone* cone = std::get_if<Cone>(IsPlane(p.A) ? &p.B : &p.A);
    if (plane == nullptr || cone == nullptr) {
        return std::optional<SurfaceSurfaceIntersection>{};
    }
    const Vec3 n = Normalized(Cross(plane->XAxis, plane->YAxis));
    if (std::abs(Dot(Apex(*cone) - plane->Origin, n)) > p.Tolerance) {
        return std::optional<SurfaceSurfaceIntersection>{};
    }
    STEREON_TRY(SurfaceSurfaceIntersection result, ApexSection(p, *plane, *cone, ctx));
    return std::optional<SurfaceSurfaceIntersection>{std::move(result)};
}

// --- A plane containing a torus's axis ------------------------------------------

/// The meridian circle of the torus at angle u, centre C + R e(u) and
/// radius r in the plane of e(u) and the axis, parametrised by the tube's
/// angle: its pcurve on the torus is the line (u, t), on the plane a circle.
[[nodiscard]] Result<void> AddMeridian(SurfaceSurfaceIntersection& result, const Pair& p,
                                       const Plane& plane, const Torus& torus, const double u,
                                       const Context& ctx)
{
    const bool planeFirst = IsPlane(p.A);
    const Box2& planeBox = planeFirst ? p.BoxA : p.BoxB;
    const Box2& torusBox = planeFirst ? p.BoxB : p.BoxA;
    if (!AngleIn(u, torusBox.X)) {
        return {};
    }
    const Vec3 e = Cos(u) * torus.XAxis + Sin(u) * torus.YAxis;
    const Vec3 a = Axis(torus);
    const Circle circle{torus.Center + torus.MajorRadius * e, torus.MinorRadius, e, a};
    const CirclePcurve onPlane =
        PcurveOf(Elementary{plane}, circle, Axis3{circle.Center, Cross(e, a)});
    STEREON_TRY(const std::vector<Interval> planeSpans,
                ClipPiece(onPlane, Interval{0.0, CirclePeriod}, planeBox));
    for (const Interval& span : planeSpans) {
        for (const Interval& piece : AngleClip(span, 0.0, 1.0, torusBox.Y)) {
            Curve2d tube(Line2d{Point2{u, 0.0}, Vec2{0.0, 1.0}});
            Curve2d flat = PieceCurve(onPlane, piece);
            IntersectionBranch branch =
                planeFirst ? Crossing(Curve(circle), std::move(flat), std::move(tube), piece)
                           : Crossing(Curve(circle), std::move(tube), std::move(flat), piece);
            STEREON_TRY_VOID(Keep(result, std::move(branch), p.ErasedA, p.ErasedB, ctx));
        }
    }
    return {};
}

/// The plane and the torus when the plane contains the torus's axis,
/// decided exactly; nothing otherwise.
[[nodiscard]] Result<std::optional<SurfaceSurfaceIntersection>> PlaneThroughAxis(const Pair& p,
                                                                                 const Context& ctx)
{
    const Plane* plane = std::get_if<Plane>(IsPlane(p.A) ? &p.A : &p.B);
    const Torus* torus = std::get_if<Torus>(IsPlane(p.A) ? &p.B : &p.A);
    if (plane == nullptr || torus == nullptr) {
        return std::optional<SurfaceSurfaceIntersection>{};
    }
    const ExactVector<4> normal = ExactAxisOf(*plane);
    const bool contains =
        SignOf(ExactDot(normal, ExactAxisOf(*torus))) == Sign::Zero &&
        SignOf(ExactDot(ExactOffset(torus->Center, plane->Origin), normal)) == Sign::Zero;
    if (!contains) {
        return std::optional<SurfaceSurfaceIntersection>{};
    }
    const Vec3 along = Cross(Cross(plane->XAxis, plane->YAxis), Axis(*torus));
    const double u = WrappedAngle(Dot(along, torus->YAxis), Dot(along, torus->XAxis));
    std::array<double, 2> angles{u, WrappedAngle(Sin(u + Pi), Cos(u + Pi))};
    std::ranges::sort(angles);
    SurfaceSurfaceIntersection result;
    for (const double angle : angles) {
        STEREON_TRY_VOID(AddMeridian(result, p, *plane, *torus, angle, ctx));
    }
    return std::optional<SurfaceSurfaceIntersection>{std::move(result)};
}

} // namespace

Result<SurfaceSurfaceIntersection> ElementaryPair(const Surface& a, const Box2& boxA,
                                                  const Surface& b, const Box2& boxB,
                                                  const double tolerance, const Context& ctx)
{
    const std::optional<Elementary> ea = ElementaryOf(a);
    const std::optional<Elementary> eb = ElementaryOf(b);
    if (!ea || !eb) {
        return Fail(ErrorCode::InvalidInput, "ElementaryPair: dispatched for other surface types");
    }
    const Pair p{a,
                 b,
                 *ea,
                 *eb,
                 boxA,
                 boxB,
                 std::max(tolerance, ctx.GetResolution()),
                 ctx.GetResolution() / 16.0};
    // An ellipsoid has the coaxial closed form when it is a spheroid on the
    // other's axis; every other configuration with one goes to the trace
    // (ADR-0048 §3.4).
    if (std::holds_alternative<Ellipsoid>(*ea) || std::holds_alternative<Ellipsoid>(*eb)) {
        STEREON_TRY(std::optional<SurfaceSurfaceIntersection> result, Coaxial(p, ctx));
        if (result) {
            return std::move(*result);
        }
        return ImplicitTrace(a, boxA, b, boxB, ctx);
    }
    using Closed =
        Result<std::optional<SurfaceSurfaceIntersection>> (*)(const Pair&, const Context&);
    for (const Closed closed : {&Coaxial, &SphereTouch, &PlaneThroughApex, &PlaneThroughAxis}) {
        STEREON_TRY(std::optional<SurfaceSurfaceIntersection> result, closed(p, ctx));
        if (result) {
            return std::move(*result);
        }
    }
    return ImplicitTrace(a, boxA, b, boxB, ctx);
}

} // namespace Stereon::Detail
