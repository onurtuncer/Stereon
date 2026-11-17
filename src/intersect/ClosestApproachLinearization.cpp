// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <Stereon/Core/Contracts.hpp>
#include <Stereon/Core/Try.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

module Stereon.Intersect;

import :Dispatch;
import :Linearization;
import :ClosestApproach;

#include "ApproachSide.hpp"

/// The products of a closest approach (ClosestApproach.cppm). Every
/// binding k moves each curve's point at fixed parameter by P_k and its
/// tangent by dP_k/du; the gradient G = (d . A', -d . B') then moves by
/// G_k = ((A_k - B_k) . A' + d . dA_k/ds, -(A_k - B_k) . B' - d . dB_k/dt),
/// and the free parameters by H z' = -G_k. The twelve values are linear in
/// each binding's rate, so the reverse product pairs the adjoint with each
/// binding's twelve: exactly the forward product's transpose.

namespace Stereon {

namespace {

using Detail::ApproachSide;

constexpr double Epsilon = std::numeric_limits<double>::epsilon();

/// A binding's motion of one curve at its parameter u: the point's
/// velocity and its derivative along the curve.
struct Motion {
    Vec3 Point;
    Vec3 Tangent;
};

[[nodiscard]] Vec3 Unit(const int k)
{
    return {k == 0 ? 1.0 : 0.0, k == 1 ? 1.0 : 0.0, k == 2 ? 1.0 : 0.0};
}

/// A line's origin moves its point alike; its direction by u times a unit,
/// its tangent by the unit.
[[nodiscard]] Motion LineMotion(const LineComponent c, const double u)
{
    const int k = static_cast<int>(c);
    if (k < 3) {
        return {Unit(k), Vec3{}};
    }
    return {u * Unit(k - 3), Unit(k - 3)};
}

/// A conic's centre moves its point alike; a circle's radius along
/// cos u X + sin u Y, an ellipse's radii along cos u X and sin u Y.
[[nodiscard]] Motion ConicMotion(const ApproachSide& side, const ConicComponent c, const double u)
{
    const double cu = std::cos(u);
    const double su = std::sin(u);
    switch (c) {
    case ConicComponent::Radius:
        return {cu * side.X + su * side.Y, -su * side.X + cu * side.Y};
    case ConicComponent::MajorRadius:
        return {cu * side.X, -su * side.X};
    case ConicComponent::MinorRadius:
        return {su * side.Y, cu * side.Y};
    default:
        return {Unit(static_cast<int>(c)), Vec3{}};
    }
}

/// A binding's parameter: of the context's set, and a Length.
[[nodiscard]] Result<void> CheckParameter(const ParameterSet& set, const ParameterId id,
                                          const std::string_view where)
{
    if (!set.Contains(id) || id.GetQuantity() != ParameterQuantity::Length) {
        return Fail(ErrorCode::InvalidInput,
                    std::string(where) +
                        ": a binding's parameter must be a Length of the context's set");
    }
    return {};
}

/// Whether a conic has the radius a component moves: a circle its Radius,
/// an ellipse its MajorRadius and MinorRadius, both their centre.
[[nodiscard]] bool HasRadius(const ApproachSide& side, const ConicComponent c)
{
    const bool circle = side.Shape == ApproachSide::Kind::Circle;
    if (c == ConicComponent::Radius) {
        return circle;
    }
    if (c == ConicComponent::MajorRadius || c == ConicComponent::MinorRadius) {
        return !circle;
    }
    return true;
}

/// One curve's bindings checked against what it is.
[[nodiscard]] Result<void> CheckBindings(const ApproachSide& side, const CurveBindings& bindings,
                                         const ParameterSet& set, const std::string_view where)
{
    if ((side.IsLine() && !bindings.Conic.empty()) || (!side.IsLine() && !bindings.Line.empty())) {
        return Fail(ErrorCode::InvalidInput,
                    std::string(where) +
                        ": a line binds by LineBinding, a circle or an ellipse by ConicBinding");
    }
    for (const LineBinding& b : bindings.Line) {
        STEREON_TRY_VOID(CheckParameter(set, b.Parameter, where));
    }
    for (const ConicBinding& b : bindings.Conic) {
        STEREON_TRY_VOID(CheckParameter(set, b.Parameter, where));
        if (!HasRadius(side, b.Component)) {
            return Fail(ErrorCode::InvalidInput,
                        std::string(where) +
                            ": a radius the conic does not have (a circle's Radius, an "
                            "ellipse's MajorRadius and MinorRadius)");
        }
    }
    return {};
}

/// The result lies on the current curves to rounding: a stale one gives a
/// plausible but wrong derivative.
[[nodiscard]] bool Consistent(const ApproachSide& a, const ApproachSide& b,
                              const ClosestApproachResult& r)
{
    const auto size = [](const Point3& p) {
        return std::max({std::abs(p.X), std::abs(p.Y), std::abs(p.Z)});
    };
    const double scale = std::max({1.0, size(r.OnA), size(r.OnB), size(a.Origin), size(b.Origin)});
    const double budget = 512.0 * Epsilon * scale;
    return std::isfinite(r.ParameterA) && std::isfinite(r.ParameterB) &&
           Distance(a.At(r.ParameterA), r.OnA) <= budget &&
           Distance(b.At(r.ParameterB), r.OnB) <= budget;
}

/// The system at the result: the gap, the tangents and curvatures, and the
/// Hessian of |d|^2 / 2.
struct System {
    Vec3 Gap;
    Vec3 A1;
    Vec3 B1;
    std::array<double, 3> Hess{}; // ss, st, tt
    std::array<double, 2> Grad{};
    bool FreeS = true;
    bool FreeT = true;
};

[[nodiscard]] System SystemOf(const ApproachSide& a, const ApproachSide& b,
                              const ClosestApproachResult& r)
{
    System s;
    s.Gap = (a.Origin - b.Origin) + (a.Offset(r.ParameterA) - b.Offset(r.ParameterB));
    s.A1 = a.D1(r.ParameterA);
    s.B1 = b.D1(r.ParameterB);
    s.Grad = {Dot(s.Gap, s.A1), -Dot(s.Gap, s.B1)};
    s.Hess = {Dot(s.A1, s.A1) + Dot(s.Gap, a.D2(r.ParameterA)), -Dot(s.A1, s.B1),
              Dot(s.B1, s.B1) - Dot(s.Gap, b.D2(r.ParameterB))};
    s.FreeS = r.EndA == ApproachEnd::Interior;
    s.FreeT = r.EndB == ApproachEnd::Interior;
    return s;
}

/// The smaller eigenvalue of the free block of the Hessian normalised by the
/// tangents' lengths: 1 for perpendicular lines, 0 for parallel ones.
[[nodiscard]] double NormalisedHessian(const System& s)
{
    const double na = Dot(s.A1, s.A1);
    const double nb = Dot(s.B1, s.B1);
    if (s.FreeS && s.FreeT) {
        const double hss = s.Hess[0] / na;
        const double htt = s.Hess[2] / nb;
        const double hst = s.Hess[1] / std::sqrt(na * nb);
        return 0.5 * (hss + htt) - std::hypot(0.5 * (hss - htt), hst);
    }
    if (s.FreeS) {
        return s.Hess[0] / na;
    }
    return s.FreeT ? s.Hess[2] / nb : 1.0;
}

/// An end's normalised one-sided gradient, |d . P'| / (|d| |P'|): how firmly
/// the minimum rests on the end; 1 for a parameter that is free.
[[nodiscard]] double EndFirmness(const System& s)
{
    const double gap = Length(s.Gap);
    double firm = 1.0;
    if (!s.FreeS) {
        firm = std::min(firm, std::abs(s.Grad[0]) / (gap * Length(s.A1)));
    }
    if (!s.FreeT) {
        firm = std::min(firm, std::abs(s.Grad[1]) / (gap * Length(s.B1)));
    }
    return firm;
}

[[nodiscard]] std::pair<DerivativeStatus, std::string_view> Classify(const ClosestApproachResult& r,
                                                                     const System& s)
{
    using D = DerivativeStatus;
    if (r.Uniqueness == ApproachUniqueness::NotUnique) {
        return {D::Nondifferentiable, "the closest points are a set, not a pair"};
    }
    if (r.Meets) {
        return {D::Nondifferentiable,
                "the curves meet: the distance has no derivative there (ADR-0052)"};
    }
    if (r.Uniqueness == ApproachUniqueness::Undecided) {
        return {D::Indeterminate, "the primal did not prove the minimiser unique"};
    }
    const double firm = EndFirmness(s);
    if (!(firm > 64.0 * Epsilon)) {
        return {D::Nondifferentiable,
                "the minimum rests on an end with no gradient: it may leave it"};
    }
    const double lambda = NormalisedHessian(s);
    if (!(lambda > 0.0) || !std::isfinite(lambda)) {
        return {D::Indeterminate, "the Hessian of the distance is singular at the minimum"};
    }
    if (lambda < NearEventSine || firm < NearEventSine) {
        return {D::NearEvent, "nearly parallel at the minimum, or barely held by an end"};
    }
    return {D::Valid, {}};
}

/// The twelve values of one binding moving at unit rate.
[[nodiscard]] std::array<double, ClosestApproachValues> ValuesOf(const System& s, const Motion& onA,
                                                                 const Motion& onB)
{
    const Vec3 dp = onA.Point - onB.Point;
    const std::array<double, 2> g{Dot(dp, s.A1) + Dot(s.Gap, onA.Tangent),
                                  -Dot(dp, s.B1) - Dot(s.Gap, onB.Tangent)};
    double ds = 0.0;
    double dt = 0.0;
    if (s.FreeS && s.FreeT) {
        const double det = s.Hess[0] * s.Hess[2] - s.Hess[1] * s.Hess[1];
        ds = (-s.Hess[2] * g[0] + s.Hess[1] * g[1]) / det;
        dt = (s.Hess[1] * g[0] - s.Hess[0] * g[1]) / det;
    } else if (s.FreeS) {
        ds = -g[0] / s.Hess[0];
    } else if (s.FreeT) {
        dt = -g[1] / s.Hess[2];
    }
    const Vec3 xa = onA.Point + ds * s.A1;
    const Vec3 xb = onB.Point + dt * s.B1;
    const Vec3 gap = xa - xb;
    const double distance = Dot(s.Gap, gap) / Length(s.Gap);
    return {xa.X, xa.Y, xa.Z, xb.X, xb.Y, xb.Z, ds, dt, gap.X, gap.Y, gap.Z, distance};
}

/// Each binding's rate-one values, scaled by its reference, with the index of
/// its parameter.
struct Contribution {
    std::size_t Parameter;
    double Reference;
    std::array<double, ClosestApproachValues> Values;
};

[[nodiscard]] std::vector<Contribution>
Contributions(const ApproachSide& a, const CurveBindings& bindingsA, const ApproachSide& b,
              const CurveBindings& bindingsB, const ClosestApproachResult& r, const System& s,
              const ParameterSet& set)
{
    std::vector<Contribution> out;
    const auto add = [&](const ParameterId id, const Motion& onA, const Motion& onB) {
        out.push_back({id.GetIndex(), set.GetReference(id), ValuesOf(s, onA, onB)});
    };
    const Motion none{};
    for (const LineBinding& binding : bindingsA.Line) {
        add(binding.Parameter, LineMotion(binding.Component, r.ParameterA), none);
    }
    for (const ConicBinding& binding : bindingsA.Conic) {
        add(binding.Parameter, ConicMotion(a, binding.Component, r.ParameterA), none);
    }
    for (const LineBinding& binding : bindingsB.Line) {
        add(binding.Parameter, none, LineMotion(binding.Component, r.ParameterB));
    }
    for (const ConicBinding& binding : bindingsB.Conic) {
        add(binding.Parameter, none, ConicMotion(b, binding.Component, r.ParameterB));
    }
    return out;
}

[[nodiscard]] std::uint64_t Signature(const ClosestApproachResult& r)
{
    return 0x5143A000ull ^ (static_cast<std::uint64_t>(r.EndA) << 8U) ^
           (static_cast<std::uint64_t>(r.EndB) << 16U) ^
           (static_cast<std::uint64_t>(r.Uniqueness) << 24U) ^
           (static_cast<std::uint64_t>(r.Meets) << 32U);
}

/// The common part of both products: the sides, the checks, the status.
struct Prepared {
    ApproachSide A;
    ApproachSide B;
    System Sys;
    DerivativeStatus Status = DerivativeStatus::Valid;
    std::string Reason;
};

[[nodiscard]] Result<std::optional<Prepared>>
Prepare(const Curve& a, const CurveBindings& bindingsA, const Curve& b,
        const CurveBindings& bindingsB, const ClosestApproachResult& r, const ParameterSet& set,
        const std::string_view where)
{
    const auto sideA = Detail::SideOf(a, where);
    const auto sideB = Detail::SideOf(b, where);
    if (!sideA || !sideB) {
        const KernelError& error = !sideA ? sideA.error() : sideB.error();
        if (error.GetCode() == ErrorCode::NotSupported) {
            return std::optional<Prepared>{};
        }
        return Fail(error.GetCode(), std::string(error.GetMessage()));
    }
    STEREON_TRY_VOID(CheckBindings(*sideA, bindingsA, set, where));
    STEREON_TRY_VOID(CheckBindings(*sideB, bindingsB, set, where));
    if (!Consistent(*sideA, *sideB, r)) {
        return Fail(ErrorCode::InvalidInput,
                    std::string(where) + ": a stale result: its points are not on the curves");
    }
    Prepared p{*sideA, *sideB, SystemOf(*sideA, *sideB, r), DerivativeStatus::Valid, {}};
    const auto [status, why] = Classify(r, p.Sys);
    p.Status = status;
    p.Reason = std::string(where) + ": " + std::string(why);
    return std::optional<Prepared>{std::move(p)};
}

[[nodiscard]] Result<SensitivityResult> Published(const Prepared& p, std::vector<double> values,
                                                  const ClosestApproachResult& r,
                                                  const std::string_view where)
{
    for (const double v : values) {
        if (!std::isfinite(v)) {
            return SensitivityResult::Without(DerivativeStatus::Indeterminate,
                                              std::string(where) + ": the derivative overflowed");
        }
    }
    STEREON_TRY(SensitivityResult result,
                SensitivityResult::Computed(p.Status, std::move(values), Signature(r)));
    const double residual =
        std::hypot(p.Sys.FreeS ? p.Sys.Grad[0] : 0.0, p.Sys.FreeT ? p.Sys.Grad[1] : 0.0);
    result.SetDiagnostics(residual, NormalisedHessian(p.Sys));
    return result;
}

[[nodiscard]] Result<SensitivityResult> Unsupported(const std::string_view where)
{
    return SensitivityResult::Without(DerivativeStatus::NotSupported,
                                      std::string(where) +
                                          ": lines, circles and ellipses only (ADR-0052, "
                                          "addendum of 2026-11-12)");
}

} // namespace

Result<SensitivityResult> ClosestApproachJvp(const Curve& a, const CurveBindings& bindingsA,
                                             const Curve& b, const CurveBindings& bindingsB,
                                             const ClosestApproachResult& approach,
                                             const SensitivityContext& sensitivity)
{
    constexpr std::string_view where = "ClosestApproachJvp";
    const ParameterSet& set = sensitivity.GetParameters();
    STEREON_TRY(const std::optional<Prepared> prepared,
                Prepare(a, bindingsA, b, bindingsB, approach, set, where));
    if (!prepared) {
        return Unsupported(where);
    }
    if (sensitivity.GetSeed().size() != set.GetSize()) {
        return Fail(ErrorCode::InvalidInput, std::string(where) + ": the context has no seed");
    }
    if (!HasDerivative(prepared->Status)) {
        return SensitivityResult::Without(prepared->Status, prepared->Reason);
    }
    std::vector<double> values(ClosestApproachValues, 0.0);
    const std::span<const double> seed = sensitivity.GetSeed();
    for (const Contribution& c : Contributions(prepared->A, bindingsA, prepared->B, bindingsB,
                                               approach, prepared->Sys, set)) {
        for (std::size_t i = 0; i < ClosestApproachValues; ++i) {
            values[i] += seed[c.Parameter] * c.Reference * c.Values[i];
        }
    }
    return Published(*prepared, std::move(values), approach, where);
}

Result<SensitivityResult> ClosestApproachVjp(const Curve& a, const CurveBindings& bindingsA,
                                             const Curve& b, const CurveBindings& bindingsB,
                                             const ClosestApproachResult& approach,
                                             const SensitivityContext& sensitivity)
{
    constexpr std::string_view where = "ClosestApproachVjp";
    const ParameterSet& set = sensitivity.GetParameters();
    STEREON_TRY(const std::optional<Prepared> prepared,
                Prepare(a, bindingsA, b, bindingsB, approach, set, where));
    if (!prepared) {
        return Unsupported(where);
    }
    const std::span<const double> adjoint = sensitivity.GetAdjoint();
    if (adjoint.size() != ClosestApproachValues) {
        return Fail(ErrorCode::InvalidInput,
                    std::string(where) + ": the adjoint must have twelve components");
    }
    if (!HasDerivative(prepared->Status)) {
        return SensitivityResult::Without(prepared->Status, prepared->Reason);
    }
    std::vector<double> values(set.GetSize(), 0.0);
    for (const Contribution& c : Contributions(prepared->A, bindingsA, prepared->B, bindingsB,
                                               approach, prepared->Sys, set)) {
        double pairing = 0.0;
        for (std::size_t i = 0; i < ClosestApproachValues; ++i) {
            pairing += adjoint[i] * c.Values[i];
        }
        values[c.Parameter] += c.Reference * pairing;
    }
    return Published(*prepared, std::move(values), approach, where);
}

} // namespace Stereon
