// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The elementary surfaces two ways (Hoffmann's parametric-implicit pairing,
// ADR-0010): as a parametrisation S(u, v), evaluated in doubles or over
// intervals with the geometry templates, and as an implicit equation
// F(x) = 0 in the surface's own frame, with its gradient. A pair of them
// meets where f(u, v) = F(S(u, v)) vanishes on the parametric side's
// parameter plane.
//
// With p, q and h the coordinates of x - O along the frame's XAxis, YAxis
// and axis, and rho² = p² + q²:
//
//   plane      F = h                          (the axis is the normal)
//   cylinder   F = rho² - r²
//   cone       F = rho² cos² a - (R cos a + h sin a)²   (both nappes)
//   sphere     F = rho² + h² - R²
//   torus      F = (rho² + h² + R² - r²)² - 4 R² rho²
//   ellipsoid  F = (B C p)² + (A C q)² + (A B h)² - (A B C)²
//
// The cone's form holds since rho = R + v sin a and h = v cos a give
// R cos a + h sin a = rho cos a. The local coordinates are a linear map L
// of x, so the chain rule gives f's partials as grad F . L(S_u) and
// grad F . L(S_v) whatever the frame's rounding, and every enclosure of f
// over a box of parameters is rigorous for the nominal F. The frame's own
// rounding, within FrameAllowance (ADR-0006), moves F's zero set from the
// evaluated surface by far less than the resolution; the certificate of
// every branch (ADR-0009) is taken against the evaluated surfaces.
//
// Included in the purview of the handlers' implementation units, after
// `import :Dispatch;`, with <algorithm>, <array>, <cmath>, <cstdint>,
// <optional>, <utility> and <variant> in their global module fragment.
// Everything is in an anonymous namespace, so each unit gets its own copy.

#pragma once

namespace Stereon::Detail {
namespace {

using I = IntervalScalar;
using IPoint2 = BasicPoint2<I>;
using IPoint3 = BasicPoint3<I>;
using IVec3 = BasicVec3<I>;

/// The elementary surfaces of the analytic table, as values.
using Elementary = std::variant<Plane, Cylinder, Cone, Sphere, Torus, Ellipsoid>;

/// The surface as one of the elementary types; nothing for any other.
[[nodiscard]] inline std::optional<Elementary> ElementaryOf(const Surface& s)
{
    if (const Plane* p = s.Target<Plane>()) {
        return Elementary{*p};
    }
    if (const Cylinder* c = s.Target<Cylinder>()) {
        return Elementary{*c};
    }
    if (const Cone* k = s.Target<Cone>()) {
        return Elementary{*k};
    }
    if (const Sphere* g = s.Target<Sphere>()) {
        return Elementary{*g};
    }
    if (const Torus* t = s.Target<Torus>()) {
        return Elementary{*t};
    }
    if (const Ellipsoid* e = s.Target<Ellipsoid>()) {
        return Elementary{*e};
    }
    return std::nullopt;
}

// --- Over intervals -------------------------------------------------------

[[nodiscard]] inline IVec3 WideVec(const Vec3& v)
{
    return {I{v.X}, I{v.Y}, I{v.Z}};
}

[[nodiscard]] inline IPoint3 WidePoint(const Point3& p)
{
    return {I{p.X}, I{p.Y}, I{p.Z}};
}

[[nodiscard]] inline BasicPlane<I> Widened(const Plane& p)
{
    return {WidePoint(p.Origin), WideVec(p.XAxis), WideVec(p.YAxis)};
}
[[nodiscard]] inline BasicCylinder<I> Widened(const Cylinder& c)
{
    return {WidePoint(c.Origin), I{c.Radius}, WideVec(c.XAxis), WideVec(c.YAxis), c.RightHanded};
}
[[nodiscard]] inline BasicCone<I> Widened(const Cone& c)
{
    return {WidePoint(c.Origin), I{c.Radius},      I{c.SemiAngle},
            WideVec(c.XAxis),    WideVec(c.YAxis), c.RightHanded};
}
[[nodiscard]] inline BasicSphere<I> Widened(const Sphere& s)
{
    return {WidePoint(s.Center), I{s.Radius}, WideVec(s.XAxis), WideVec(s.YAxis), s.RightHanded};
}
[[nodiscard]] inline BasicTorus<I> Widened(const Torus& t)
{
    return {WidePoint(t.Center), I{t.MajorRadius}, I{t.MinorRadius},
            WideVec(t.XAxis),    WideVec(t.YAxis), t.RightHanded};
}
[[nodiscard]] inline BasicEllipsoid<I> Widened(const Ellipsoid& e)
{
    return {WidePoint(e.Center), I{e.A},           I{e.B},       I{e.C},
            WideVec(e.XAxis),    WideVec(e.YAxis), e.RightHanded};
}

/// A point of the parametric side and its first partials.
template <Scalar T>
struct Jet {
    BasicPoint3<T> Point;
    BasicVec3<T> Su;
    BasicVec3<T> Sv;
};

/// S, S_u and S_v at (u, v), in doubles.
[[nodiscard]] inline Jet<double> JetAt(const Elementary& a, const Point2& uv)
{
    return std::visit(
        [&](const auto& g) {
            return Jet<double>{Eval(g, uv), Deriv(g, uv, 1, 0), Deriv(g, uv, 0, 1)};
        },
        a);
}

/// S, S_u and S_v enclosed over a box of parameters.
[[nodiscard]] inline Jet<I> JetOver(const Elementary& a, const IPoint2& uv)
{
    return std::visit(
        [&](const auto& g) {
            const auto w = Widened(g);
            return Jet<I>{Eval(w, uv), Deriv(w, uv, 1, 0), Deriv(w, uv, 0, 1)};
        },
        a);
}

/// The point at (u, v) in doubles.
[[nodiscard]] inline Point3 PointAt(const Elementary& a, const Point2& uv)
{
    return std::visit([&](const auto& g) { return Eval(g, uv); }, a);
}

/// The parameters of the point of the surface closest to p, in the
/// surface's own chart.
[[nodiscard]] inline Point2 ProjectOnto(const Elementary& a, const Point3& p)
{
    return std::visit([&](const auto& g) { return Project(g, p); }, a);
}

// --- The implicit side ----------------------------------------------------

enum class QuadricKind : std::uint8_t { Plane, Cylinder, Cone, Sphere, Torus, Ellipsoid };

/// The implicit side: its kind, its frame with the axis as the evaluator
/// takes it (the plane's normal), and its sizes. Radius is the cylinder's,
/// the sphere's, the cone's reference radius, the torus's minor one and the
/// ellipsoid's semi-axis A; SemiB and SemiC are the ellipsoid's other two.
struct Quadric {
    QuadricKind Kind = QuadricKind::Plane;
    Point3 Origin;
    Vec3 XAxis;
    Vec3 YAxis;
    Vec3 Axis;
    double Radius = 0.0;
    double Major = 0.0;
    double SemiAngle = 0.0;
    double SemiB = 0.0;
    double SemiC = 0.0;
};

[[nodiscard]] inline Quadric QuadricOf(const Plane& p)
{
    return {QuadricKind::Plane, p.Origin, p.XAxis, p.YAxis, Cross(p.XAxis, p.YAxis)};
}
[[nodiscard]] inline Quadric QuadricOf(const Cylinder& c)
{
    return {QuadricKind::Cylinder, c.Origin, c.XAxis, c.YAxis, Axis(c), c.Radius};
}
[[nodiscard]] inline Quadric QuadricOf(const Cone& c)
{
    return {QuadricKind::Cone, c.Origin, c.XAxis, c.YAxis, Axis(c), c.Radius, 0.0, c.SemiAngle};
}
[[nodiscard]] inline Quadric QuadricOf(const Sphere& s)
{
    return {QuadricKind::Sphere, s.Center, s.XAxis, s.YAxis, Axis(s), s.Radius};
}
[[nodiscard]] inline Quadric QuadricOf(const Torus& t)
{
    return {QuadricKind::Torus, t.Center, t.XAxis, t.YAxis, Axis(t), t.MinorRadius, t.MajorRadius};
}
[[nodiscard]] inline Quadric QuadricOf(const Ellipsoid& e)
{
    return {QuadricKind::Ellipsoid, e.Center, e.XAxis, e.YAxis, Axis(e), e.A, 0.0, 0.0, e.B, e.C};
}

[[nodiscard]] inline Quadric QuadricOf(const Elementary& e)
{
    return std::visit([](const auto& g) { return QuadricOf(g); }, e);
}

/// True for the bounded surfaces, whose whole parameter domain is a chart.
[[nodiscard]] inline bool IsCompact(const QuadricKind kind)
{
    return kind == QuadricKind::Sphere || kind == QuadricKind::Torus ||
           kind == QuadricKind::Ellipsoid;
}

/// The local coordinates (p, q, h) of a direction: the linear part of L.
template <Scalar T>
[[nodiscard]] BasicVec3<T> LocalDirection(const Quadric& s, const BasicVec3<T>& d)
{
    using V = BasicVec3<T>;
    const V x{T{s.XAxis.X}, T{s.XAxis.Y}, T{s.XAxis.Z}};
    const V y{T{s.YAxis.X}, T{s.YAxis.Y}, T{s.YAxis.Z}};
    const V a{T{s.Axis.X}, T{s.Axis.Y}, T{s.Axis.Z}};
    return {Dot(d, x), Dot(d, y), Dot(d, a)};
}

/// The local coordinates (p, q, h) of a point.
template <Scalar T>
[[nodiscard]] BasicVec3<T> LocalPoint(const Quadric& s, const BasicPoint3<T>& x)
{
    const BasicPoint3<T> o{T{s.Origin.X}, T{s.Origin.Y}, T{s.Origin.Z}};
    return LocalDirection(s, BasicVec3<T>{x - o});
}

/// F at local coordinates.
template <Scalar T>
[[nodiscard]] T ImplicitValue(const Quadric& s, const BasicVec3<T>& l)
{
    const T rho2 = l.X * l.X + l.Y * l.Y;
    const T r{s.Radius};
    switch (s.Kind) {
    case QuadricKind::Plane:
        return l.Z;
    case QuadricKind::Cylinder:
        return rho2 - r * r;
    case QuadricKind::Cone: {
        const T c = Cos(T{s.SemiAngle});
        const T k = r * c + l.Z * Sin(T{s.SemiAngle});
        return rho2 * c * c - k * k;
    }
    case QuadricKind::Sphere:
        return rho2 + l.Z * l.Z - r * r;
    case QuadricKind::Torus: {
        const T big{s.Major};
        const T k = rho2 + l.Z * l.Z + big * big - r * r;
        return k * k - T{4.0} * big * big * rho2;
    }
    case QuadricKind::Ellipsoid: {
        const T bc = T{s.SemiB} * T{s.SemiC};
        const T ac = r * T{s.SemiC};
        const T ab = r * T{s.SemiB};
        const T abc = r * bc;
        return bc * bc * l.X * l.X + ac * ac * l.Y * l.Y + ab * ab * l.Z * l.Z - abc * abc;
    }
    }
    return T{};
}

/// grad F at local coordinates, in local coordinates.
template <Scalar T>
[[nodiscard]] BasicVec3<T> ImplicitGradient(const Quadric& s, const BasicVec3<T>& l)
{
    const T r{s.Radius};
    switch (s.Kind) {
    case QuadricKind::Plane:
        return {T{}, T{}, T{1.0}};
    case QuadricKind::Cylinder:
        return {T{2.0} * l.X, T{2.0} * l.Y, T{}};
    case QuadricKind::Cone: {
        const T c = Cos(T{s.SemiAngle});
        const T n = Sin(T{s.SemiAngle});
        const T k = r * c + l.Z * n;
        return {T{2.0} * c * c * l.X, T{2.0} * c * c * l.Y, T{-2.0} * n * k};
    }
    case QuadricKind::Sphere:
        return {T{2.0} * l.X, T{2.0} * l.Y, T{2.0} * l.Z};
    case QuadricKind::Torus: {
        const T big{s.Major};
        const T k = l.X * l.X + l.Y * l.Y + l.Z * l.Z + big * big - r * r;
        const T radial = T{4.0} * (k - T{2.0} * big * big);
        return {radial * l.X, radial * l.Y, T{4.0} * k * l.Z};
    }
    case QuadricKind::Ellipsoid: {
        const T bc = T{s.SemiB} * T{s.SemiC};
        const T ac = r * T{s.SemiC};
        const T ab = r * T{s.SemiB};
        return {T{2.0} * bc * bc * l.X, T{2.0} * ac * ac * l.Y, T{2.0} * ab * ab * l.Z};
    }
    }
    return {};
}

/// f = F(S(u, v)) and its partials f_u and f_v, from the parametric side's
/// jet: the chain rule through the linear map L.
template <Scalar T>
struct Composite {
    T F;
    T Fu;
    T Fv;
};

template <Scalar T>
[[nodiscard]] Composite<T> CompositeOf(const Quadric& s, const Jet<T>& jet)
{
    const BasicVec3<T> local = LocalPoint(s, jet.Point);
    const BasicVec3<T> gradient = ImplicitGradient(s, local);
    return {ImplicitValue(s, local), Dot(gradient, LocalDirection(s, jet.Su)),
            Dot(gradient, LocalDirection(s, jet.Sv))};
}

/// f over a box, as the meet of its natural enclosure and its mean-value
/// form about the box's midpoint, which narrows with the box: the natural
/// form alone overestimates by the dependency of the terms.
struct BoxValues {
    Interval F;
    Interval Fu;
    Interval Fv;
};

[[nodiscard]] inline BoxValues ValuesOver(const Elementary& a, const Quadric& b, const Interval& u,
                                          const Interval& v)
{
    const Composite<I> over = CompositeOf(b, JetOver(a, IPoint2{I{u}, I{v}}));
    const double mu = 0.5 * (u.GetLower() + u.GetUpper());
    const double mv = 0.5 * (v.GetLower() + v.GetUpper());
    const Composite<I> mid = CompositeOf(b, JetOver(a, IPoint2{I{mu}, I{mv}}));
    const Interval du =
        Hull(Interval{u.GetLower()} - Interval{mu}, Interval{u.GetUpper()} - Interval{mu});
    const Interval dv =
        Hull(Interval{v.GetLower()} - Interval{mv}, Interval{v.GetUpper()} - Interval{mv});
    const Interval meanValue = mid.F.Value + over.Fu.Value * du + over.Fv.Value * dv;
    const Interval met = Intersection(over.F.Value, meanValue);
    return {met.IsEmpty() ? over.F.Value : met, over.Fu.Value, over.Fv.Value};
}

/// f at one parameter point, enclosed: its sign is certain when the
/// enclosure misses zero.
[[nodiscard]] inline Interval ValueAt(const Elementary& a, const Quadric& b, const Point2& uv)
{
    return CompositeOf(b, JetOver(a, IPoint2{I{uv.X}, I{uv.Y}})).F.Value;
}

/// f and its partials at a parameter point, in doubles: for placing roots,
/// never for a decision.
[[nodiscard]] inline Composite<double> NominalAt(const Elementary& a, const Quadric& b,
                                                 const Point2& uv)
{
    return CompositeOf(b, JetAt(a, uv));
}

/// -1, 0 or +1: the certain sign of an enclosure, 0 when it holds zero.
[[nodiscard]] inline int CertainSign(const Interval& x)
{
    if (x.GetLower() > 0.0) {
        return 1;
    }
    return x.GetUpper() < 0.0 ? -1 : 0;
}

} // namespace
} // namespace Stereon::Detail
