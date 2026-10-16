// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <Stereon/Core/Contracts.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <expected>
#include <span>
#include <utility>

export module Stereon.Surface:Properties;

import Stereon.Curve;
import :Concepts;

/// Local differential properties of surfaces (OCCT's GeomLProp_SLProps):
/// the unit normal, the two fundamental forms, and the principal, Gaussian
/// and mean curvatures, for any surface type with Deriv up to second order.
///
/// Curvatures are signed with respect to the unit normal n = Su x Sv /
/// |Su x Sv|: positive where the surface bends towards n. A sphere, whose
/// normal points outward in a right-handed frame, therefore has principal
/// curvatures -1/R; reversing the surface flips the normal and the signs of
/// the principal and mean curvatures, never the Gaussian curvature.
///
/// Where Su x Sv = 0, at a pole of a sphere or the apex of a cone, the
/// parametrisation is singular, the normal is undefined by this formula
/// and every property returns ErrorCode::DegenerateResult. As for curves
/// the test is exact, without a hidden tolerance (ADR-0028).

export namespace Stereon {

/// The coefficients of the first fundamental form, E du^2 + 2 F du dv +
/// G dv^2, which measures lengths on the surface, and of the second,
/// L du^2 + 2 M du dv + N dv^2, which measures how it leaves its tangent
/// plane, with the unit normal they refer to.
struct FundamentalFormsResult {
    double E = 0.0;
    double F = 0.0;
    double G = 0.0;
    double L = 0.0;
    double M = 0.0;
    double N = 0.0;
    Vec3 Normal;
};

/// The principal curvatures, Min <= Max, and their directions: unit
/// tangent vectors, orthogonal, with Cross(MinDirection, MaxDirection) =
/// Normal. At an umbilic, where Min = Max and every direction is
/// principal, the directions are an arbitrary such pair.
struct PrincipalCurvaturesResult {
    double Min = 0.0;
    double Max = 0.0;
    Vec3 MinDirection;
    Vec3 MaxDirection;
};

/// The unit normal Su x Sv / |Su x Sv|. DegenerateResult where the
/// parametrisation is singular.
template <typename S>
[[nodiscard]] Result<Vec3> UnitNormal(const S& surface, const Point2& uv)
{
    const Vec3 normal = Cross(Deriv(surface, uv, 1, 0), Deriv(surface, uv, 0, 1));
    const double length = Length(normal);
    if (!(length > 0.0)) {
        return Fail(ErrorCode::DegenerateResult,
                    "UnitNormal: the parametrisation is singular here (Su x Sv = 0)");
    }
    return normal / length;
}

} // namespace Stereon

export namespace Stereon::Detail {

/// The first and second partial derivatives at a point. A template only
/// because GCC 16.2 fails to lazy-load BasicVec3 in importers of a module
/// interface that defines a non-template struct of Vec3 members here
/// (docs/toolchain/gcc16-feature-probe.md).
template <typename V>
struct SecondPartials {
    V Su;
    V Sv;
    V Suu;
    V Suv;
    V Svv;
};

/// The partials from one combined evaluation where the type has one
/// (TotalDerivatives of a B-spline surface: one pair of basis tables, and
/// only the partials of total order at most 2), else from
/// five calls of Deriv.
template <typename S, typename V = decltype(Deriv(std::declval<const S&>(), Point2{}, 1, 0))>
[[nodiscard]] SecondPartials<V> PartialsAt(const S& surface, const Point2& uv)
{
    if constexpr (requires(std::span<V> out) { TotalDerivatives(surface, uv, 2, out); }) {
        // out[k * 3 + l] = S^(k,l).
        std::array<V, 9> out;
        TotalDerivatives(surface, uv, 2, std::span<V>{out});
        return {out[3], out[1], out[6], out[4], out[2]};
    } else {
        return {Deriv(surface, uv, 1, 0), Deriv(surface, uv, 0, 1), Deriv(surface, uv, 2, 0),
                Deriv(surface, uv, 1, 1), Deriv(surface, uv, 0, 2)};
    }
}

/// Both fundamental forms from the partials. DegenerateResult where
/// Su x Sv = 0.
template <typename Partials>
[[nodiscard]] Result<FundamentalFormsResult> FormsFrom(const Partials& d)
{
    const Vec3 cross = Cross(d.Su, d.Sv);
    const double length = Length(cross);
    if (!(length > 0.0)) {
        return Fail(ErrorCode::DegenerateResult,
                    "the parametrisation is singular here (Su x Sv = 0)");
    }
    const Vec3 normal = cross / length;
    return FundamentalFormsResult{.E = Dot(d.Su, d.Su),
                                  .F = Dot(d.Su, d.Sv),
                                  .G = Dot(d.Sv, d.Sv),
                                  .L = Dot(d.Suu, normal),
                                  .M = Dot(d.Suv, normal),
                                  .N = Dot(d.Svv, normal),
                                  .Normal = normal};
}

/// A principal direction for the curvature kappa: a null vector (du, dv) of
/// II - kappa I, from whichever row of that symmetric 2 x 2 matrix is
/// larger, as a unit tangent vector. Su when both rows vanish, at an
/// umbilic, where every direction is principal. A template only because
/// GCC 16 cannot read back an inline function of a module interface that
/// instantiates another module's templates
/// (docs/toolchain/gcc16-feature-probe.md).
template <typename Forms>
[[nodiscard]] Vec3 PrincipalDirection(const Forms& f, const double kappa, const Vec3& su,
                                      const Vec3& sv)
{
    const double a = f.L - kappa * f.E;
    const double b = f.M - kappa * f.F;
    const double c = f.N - kappa * f.G;
    const bool firstRow = a * a + b * b >= b * b + c * c;
    const double du = firstRow ? -b : -c;
    const double dv = firstRow ? a : b;
    const Vec3 direction = du * su + dv * sv;
    const double length = Length(direction);
    return length > 0.0 ? direction / length : su / Length(su);
}

} // namespace Stereon::Detail

export namespace Stereon {

/// Both fundamental forms at uv, from one evaluation of the partials.
/// DegenerateResult where the parametrisation is singular.
template <typename S>
[[nodiscard]] Result<FundamentalFormsResult> FundamentalForms(const S& surface, const Point2& uv)
{
    return Detail::FormsFrom(Detail::PartialsAt(surface, uv));
}

/// The Gaussian curvature (LN - M^2) / (EG - F^2), the product of the
/// principal curvatures. DegenerateResult where the parametrisation is
/// singular.
template <typename S>
[[nodiscard]] Result<double> GaussianCurvature(const S& surface, const Point2& uv)
{
    const Result<FundamentalFormsResult> forms = FundamentalForms(surface, uv);
    if (!forms) {
        return std::unexpected(forms.error());
    }
    const FundamentalFormsResult& f = *forms;
    return (f.L * f.N - f.M * f.M) / (f.E * f.G - f.F * f.F);
}

/// The mean curvature (EN - 2FM + GL) / (2 (EG - F^2)), the average of the
/// principal curvatures. DegenerateResult where the parametrisation is
/// singular.
template <typename S>
[[nodiscard]] Result<double> MeanCurvature(const S& surface, const Point2& uv)
{
    const Result<FundamentalFormsResult> forms = FundamentalForms(surface, uv);
    if (!forms) {
        return std::unexpected(forms.error());
    }
    const FundamentalFormsResult& f = *forms;
    return (f.E * f.N - 2.0 * f.F * f.M + f.G * f.L) / (2.0 * (f.E * f.G - f.F * f.F));
}

/// The principal curvatures and directions: the eigenvalues and
/// eigenvectors of the shape operator I^-1 II, from H +- sqrt(H^2 - K)
/// computed without cancellation, all from one evaluation of the partials.
/// DegenerateResult where the parametrisation is singular.
template <typename S>
[[nodiscard]] Result<PrincipalCurvaturesResult> PrincipalCurvatures(const S& surface,
                                                                    const Point2& uv)
{
    const auto d = Detail::PartialsAt(surface, uv);
    const Result<FundamentalFormsResult> forms = Detail::FormsFrom(d);
    if (!forms) {
        return std::unexpected(forms.error());
    }
    const FundamentalFormsResult& f = *forms;
    const double metric = f.E * f.G - f.F * f.F;
    const double gauss = (f.L * f.N - f.M * f.M) / metric;
    const double mean = (f.E * f.N - 2.0 * f.F * f.M + f.G * f.L) / (2.0 * metric);
    // H^2 - K >= 0 in exact arithmetic; rounding may push it below at an
    // umbilic. The larger root in magnitude first, the other from K = k1 k2.
    const double spread = std::sqrt(std::max(mean * mean - gauss, 0.0));
    const double far = mean >= 0.0 ? mean + spread : mean - spread;
    const double near = far != 0.0 ? gauss / far : 0.0;
    const double kMin = std::min(far, near);
    const double kMax = std::max(far, near);
    const Vec3 maxDirection = Detail::PrincipalDirection(f, kMax, d.Su, d.Sv);
    return PrincipalCurvaturesResult{.Min = kMin,
                                     .Max = kMax,
                                     .MinDirection = Cross(maxDirection, f.Normal),
                                     .MaxDirection = maxDirection};
}

} // namespace Stereon
