// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

// The contract anchor (docs/toolchain/gcc16-feature-probe.md); see
// src/curve/BSplineCurve.cpp.
#include <Stereon/Core/Contracts.hpp>
#include <Stereon/Core/Try.hpp>

#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

module Stereon.Curve;

import :BSplineBasis;
import :BSplineCurve;
import :PoleSensitivity;

namespace Stereon {

namespace {

/// The basis of one point: the first pole it reaches, N_j and R_j for the
/// p + 1 poles from there, and C(t) and W(t).
struct PointBasis {
    std::size_t First = 0;
    std::size_t Count = 0;
    std::array<double, Detail::BSplineBasisCount> N{};
    std::array<double, Detail::BSplineBasisCount> R{};
    Vec3 Point;
    double W = 1.0;
};

[[nodiscard]] PointBasis BasisAt(const BSplineCurve& curve, const double t)
{
    const std::size_t span = FindSpan(curve.Knots, curve.Degree, curve.Poles.size(), t);
    const auto p = static_cast<std::size_t>(curve.Degree);
    PointBasis basis;
    basis.First = span - p;
    basis.Count = p + 1;
    BasisFunctions(curve.Knots, curve.Degree, span, t, std::span{basis.N}.first(p + 1));
    double w = 0.0;
    for (std::size_t j = 0; j <= p; ++j) {
        w += basis.N[j] * Detail::WeightOf(curve, basis.First + j);
    }
    basis.W = w;
    Vec3 point{};
    for (std::size_t j = 0; j <= p; ++j) {
        basis.R[j] = basis.N[j] * Detail::WeightOf(curve, basis.First + j) / w;
        point += basis.R[j] * ToVec(curve.Poles[basis.First + j]);
    }
    basis.Point = point;
    return basis;
}

/// dC/d(component) of one binding at one point, per unit of the component.
[[nodiscard]] Vec3 Partial(const BSplineCurve& curve, const PointBasis& basis,
                           const PoleBinding& binding)
{
    if (binding.Pole < basis.First || binding.Pole >= basis.First + basis.Count) {
        return {};
    }
    const std::size_t j = binding.Pole - basis.First;
    switch (binding.Component) {
    case PoleComponent::X:
        return {basis.R[j], 0.0, 0.0};
    case PoleComponent::Y:
        return {0.0, basis.R[j], 0.0};
    case PoleComponent::Z:
        return {0.0, 0.0, basis.R[j]};
    case PoleComponent::Weight:
        return (basis.N[j] / basis.W) * (ToVec(curve.Poles[binding.Pole]) - basis.Point);
    }
    return {};
}

[[nodiscard]] Result<void> CheckInput(const BSplineCurve& curve,
                                      const std::span<const PoleBinding> bindings,
                                      const std::span<const double> parameters,
                                      const ParameterSet& set, const std::string_view function)
{
    for (const PoleBinding& binding : bindings) {
        if (!set.Contains(binding.Parameter)) {
            return Fail(ErrorCode::InvalidInput,
                        std::string(function) + ": a binding's parameter is not of the set");
        }
        if (binding.Pole >= curve.Poles.size()) {
            return Fail(ErrorCode::InvalidInput,
                        std::string(function) + ": a binding names a pole the curve does not have");
        }
        const ParameterQuantity wanted = binding.Component == PoleComponent::Weight
                                             ? ParameterQuantity::Dimensionless
                                             : ParameterQuantity::Length;
        if (binding.Parameter.GetQuantity() != wanted) {
            return Fail(ErrorCode::InvalidInput,
                        std::string(function) +
                            ": a coordinate is moved by a length and a weight by a "
                            "dimensionless parameter");
        }
    }
    const Interval domain = Domain(curve);
    for (const double t : parameters) {
        if (!std::isfinite(t) || !Contains(domain, t)) {
            return Fail(ErrorCode::InvalidInput,
                        std::string(function) + ": a parameter outside the curve's domain");
        }
    }
    return {};
}

/// The regime of the derivative: the degree, the pole count, whether the
/// curve is rational and its knots, which every evaluation holds fixed.
[[nodiscard]] std::uint64_t Signature(const BSplineCurve& curve)
{
    std::uint64_t hash = 14695981039346656037ull; // FNV-1a
    const auto mix = [&hash](const std::uint64_t value) {
        for (int byte = 0; byte < 8; ++byte) {
            hash = (hash ^ ((value >> (8 * byte)) & 0xFFu)) * 1099511628211ull;
        }
    };
    mix(static_cast<std::uint64_t>(curve.Degree));
    mix(curve.Poles.size());
    mix(IsRational(curve) ? 1u : 0u);
    for (const double u : curve.Knots) {
        mix(std::bit_cast<std::uint64_t>(u));
    }
    return hash;
}

/// The product's values as a Valid result, unless the arithmetic overflowed:
/// finite inputs can still give an infinite or NaN component, which is no
/// derivative. That is reported as Indeterminate, naming the component,
/// rather than published as Valid (SensitivityResult::Computed refuses it).
[[nodiscard]] Result<SensitivityResult>
Published(std::vector<double> values, const std::uint64_t signature, const std::string_view where)
{
    for (std::size_t k = 0; k < values.size(); ++k) {
        if (!std::isfinite(values[k])) {
            return SensitivityResult::Without(DerivativeStatus::Indeterminate,
                                              std::string(where) +
                                                  ": the derivative overflowed: component " +
                                                  std::to_string(k) + " is not finite");
        }
    }
    return SensitivityResult::Computed(DerivativeStatus::Valid, std::move(values), signature);
}

} // namespace

Result<SensitivityResult> EvalJvp(const BSplineCurve& curve,
                                  const std::span<const PoleBinding> bindings,
                                  const std::span<const double> parameters,
                                  const SensitivityContext& sensitivity)
{
    const ParameterSet& set = sensitivity.GetParameters();
    STEREON_TRY_VOID(CheckInput(curve, bindings, parameters, set, "EvalJvp"));
    const std::span<const double> seed = sensitivity.GetSeed();
    if (seed.size() != set.GetSize()) {
        return Fail(ErrorCode::InvalidInput, "EvalJvp: the context has no seed");
    }
    std::vector<double> values(3 * parameters.size(), 0.0);
    for (std::size_t k = 0; k < parameters.size(); ++k) {
        const PointBasis basis = BasisAt(curve, parameters[k]);
        Vec3 velocity{};
        for (const PoleBinding& binding : bindings) {
            const double rate =
                seed[binding.Parameter.GetIndex()] * set.GetReference(binding.Parameter);
            velocity += rate * Partial(curve, basis, binding);
        }
        values[3 * k] = velocity.X;
        values[3 * k + 1] = velocity.Y;
        values[3 * k + 2] = velocity.Z;
    }
    return Published(std::move(values), Signature(curve), "EvalJvp");
}

Result<SensitivityResult> EvalVjp(const BSplineCurve& curve,
                                  const std::span<const PoleBinding> bindings,
                                  const std::span<const double> parameters,
                                  const SensitivityContext& sensitivity)
{
    const ParameterSet& set = sensitivity.GetParameters();
    STEREON_TRY_VOID(CheckInput(curve, bindings, parameters, set, "EvalVjp"));
    const std::span<const double> adjoint = sensitivity.GetAdjoint();
    if (adjoint.size() != 3 * parameters.size()) {
        return Fail(ErrorCode::InvalidInput,
                    "EvalVjp: the adjoint must have three components per point");
    }
    std::vector<double> values(set.GetSize(), 0.0);
    for (std::size_t k = 0; k < parameters.size(); ++k) {
        const PointBasis basis = BasisAt(curve, parameters[k]);
        const Vec3 bar{adjoint[3 * k], adjoint[3 * k + 1], adjoint[3 * k + 2]};
        for (const PoleBinding& binding : bindings) {
            values[binding.Parameter.GetIndex()] +=
                set.GetReference(binding.Parameter) * Dot(bar, Partial(curve, basis, binding));
        }
    }
    return Published(std::move(values), Signature(curve), "EvalVjp");
}

} // namespace Stereon
