// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Dual.hpp"
#include "RapidCheck.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

import Stereon.Curve;

using Stereon::BSplineCurve;
using Stereon::DerivativeProduct;
using Stereon::DerivativeStatus;
using Stereon::ParameterQuantity;
using Stereon::ParameterSet;
using Stereon::Point3;
using Stereon::PoleBinding;
using Stereon::PoleComponent;
using Stereon::SensitivityContext;
using StereonTest::Dual;

// EvalJvp and EvalVjp (ADR-0021, stage 2) against an independent
// derivative: forward-mode duals pushed through the scalar-generic Eval.

namespace {

/// A random valid curve: degree 1 to 5, clamped or not, optionally rational.
BSplineCurve RandomCurve(const bool rational)
{
    const int p = *rc::gen::inRange(1, 6);
    const auto q = static_cast<std::size_t>(p);
    const std::size_t poles = q + 1 + static_cast<std::size_t>(*rc::gen::inRange(0, 6));
    const std::size_t count = poles + q + 1;
    std::vector<double> knots;
    double u = 0.0;
    while (knots.size() < count) {
        u += *rc::gen::inRange(1, 5) * 0.25;
        const auto repeat = static_cast<std::size_t>(*rc::gen::inRange(1, p + 1));
        knots.insert(knots.end(), std::min(repeat, count - knots.size()), u);
    }
    if (*rc::gen::arbitrary<bool>()) {
        std::fill(knots.begin(), knots.begin() + p + 1, knots[q]);
        std::fill(knots.end() - p - 1, knots.end(), knots[poles]);
    }
    std::vector<Point3> points;
    std::vector<double> weights;
    for (std::size_t i = 0; i < poles; ++i) {
        points.push_back({*rc::gen::inRange(-50, 50) / 10.0, *rc::gen::inRange(-50, 50) / 10.0,
                          *rc::gen::inRange(-50, 50) / 10.0});
        if (rational) {
            weights.push_back(*rc::gen::inRange(2, 40) / 10.0);
        }
    }
    const auto made = Stereon::MakeBSplineCurve(p, knots, points, weights);
    RC_PRE(made.has_value());
    return *made;
}

std::vector<double> RandomParameters(const BSplineCurve& c)
{
    const Stereon::Interval d = Domain(c);
    std::vector<double> t;
    const int count = *rc::gen::inRange(1, 6);
    for (int k = 0; k < count; ++k) {
        const double s = *rc::gen::inRange(0, 1001) / 1000.0;
        t.push_back(d.GetLower() + s * (d.GetUpper() - d.GetLower()));
    }
    return t;
}

/// A model of a few parameters with random references, and bindings of
/// each to random components of random poles.
struct Model {
    ParameterSet Parameters;
    std::vector<PoleBinding> Bindings;
};

Model RandomModel(const BSplineCurve& c)
{
    Model model;
    const int lengths = *rc::gen::inRange(1, 4);
    const int weights = *rc::gen::inRange(0, 3);
    for (int k = 0; k < lengths + weights; ++k) {
        const bool isWeight = k >= lengths;
        const auto id = model.Parameters.Add("p" + std::to_string(k),
                                             isWeight ? ParameterQuantity::Dimensionless
                                                      : ParameterQuantity::Length,
                                             *rc::gen::inRange(1, 40) / 10.0);
        RC_ASSERT(id.has_value());
        const int bound = *rc::gen::inRange(1, 4);
        for (int b = 0; b < bound; ++b) {
            const auto pole = *rc::gen::inRange<std::size_t>(0, c.Poles.size());
            const auto component = isWeight ? PoleComponent::Weight
                                            : static_cast<PoleComponent>(*rc::gen::inRange(0, 3));
            model.Bindings.push_back({*id, pole, component});
        }
    }
    return model;
}

std::vector<double> RandomVector(const std::size_t size)
{
    std::vector<double> v;
    for (std::size_t i = 0; i < size; ++i) {
        v.push_back(*rc::gen::inRange(-100, 101) / 50.0);
    }
    return v;
}

/// The curve with every pole and weight a dual whose derivative is its
/// rate under the seed: sum over its bindings of seed_k * ref_k.
Stereon::BasicBSplineCurve<Dual, 3> Lifted(const BSplineCurve& c, const Model& model,
                                           const std::vector<double>& seed)
{
    Stereon::BasicBSplineCurve<Dual, 3> lifted{c.Degree, c.Knots, {}, {}};
    for (const Point3& pole : c.Poles) {
        lifted.Poles.push_back({Dual{pole.X}, Dual{pole.Y}, Dual{pole.Z}});
    }
    const bool weighted =
        !c.Weights.empty() || std::ranges::any_of(model.Bindings, [](const PoleBinding& binding) {
            return binding.Component == PoleComponent::Weight;
        });
    if (weighted) {
        for (std::size_t i = 0; i < c.Poles.size(); ++i) {
            lifted.Weights.push_back(Dual{c.Weights.empty() ? 1.0 : c.Weights[i]});
        }
    }
    for (const PoleBinding& binding : model.Bindings) {
        const double rate =
            seed[binding.Parameter.GetIndex()] * model.Parameters.GetReference(binding.Parameter);
        auto& pole = lifted.Poles[binding.Pole];
        switch (binding.Component) {
        case PoleComponent::X:
            pole.X.Deriv += rate;
            break;
        case PoleComponent::Y:
            pole.Y.Deriv += rate;
            break;
        case PoleComponent::Z:
            pole.Z.Deriv += rate;
            break;
        case PoleComponent::Weight:
            lifted.Weights[binding.Pole].Deriv += rate;
            break;
        }
    }
    return lifted;
}

} // namespace

TEST_CASE("EvalJvp is the derivative forward-mode duals give Eval",
          "[curve][sensitivity][property]")
{
    rc::prop("for random curves, bindings, references and seeds", [] {
        const BSplineCurve c = RandomCurve(*rc::gen::arbitrary<bool>());
        const std::vector<double> t = RandomParameters(c);
        Model model = RandomModel(c);
        const std::vector<double> seed = RandomVector(model.Parameters.GetSize());
        const auto lifted = Lifted(c, model, seed);
        const std::vector<PoleBinding> bindings = model.Bindings;
        SensitivityContext sensitivity(std::move(model.Parameters));
        RC_ASSERT(sensitivity.SetSeed(seed).has_value());

        const auto jvp = Stereon::EvalJvp(c, bindings, t, sensitivity);
        RC_ASSERT(jvp.has_value());
        RC_ASSERT(jvp->GetStatus() == DerivativeStatus::Valid);
        const auto values = jvp->GetValues();
        RC_ASSERT(values.has_value());
        RC_ASSERT(values->size() == 3 * t.size());
        for (std::size_t k = 0; k < t.size(); ++k) {
            const auto point = Eval(lifted, t[k]);
            const double expected[] = {point.X.Deriv, point.Y.Deriv, point.Z.Deriv};
            for (std::size_t d = 0; d < 3; ++d) {
                RC_ASSERT(std::abs((*values)[3 * k + d] - expected[d]) <=
                          1e-11 * (1.0 + std::abs(expected[d])));
            }
        }
    });
}

TEST_CASE("EvalVjp is the transpose of EvalJvp", "[curve][sensitivity][property]")
{
    rc::prop("<x', J p'> = <J^T x', p'> for random adjoints and seeds", [] {
        const BSplineCurve c = RandomCurve(*rc::gen::arbitrary<bool>());
        const std::vector<double> t = RandomParameters(c);
        Model model = RandomModel(c);
        const std::vector<double> seed = RandomVector(model.Parameters.GetSize());
        const std::vector<double> adjoint = RandomVector(3 * t.size());
        const std::vector<PoleBinding> bindings = model.Bindings;
        SensitivityContext sensitivity(std::move(model.Parameters));
        RC_ASSERT(sensitivity.SetSeed(seed).has_value());
        RC_ASSERT(sensitivity.SetAdjoint(adjoint).has_value());

        const auto jvp = Stereon::EvalJvp(c, bindings, t, sensitivity);
        const auto vjp = Stereon::EvalVjp(c, bindings, t, sensitivity);
        RC_ASSERT(jvp.has_value() && vjp.has_value());
        RC_ASSERT(vjp->GetStatus() == DerivativeStatus::Valid);
        const auto forward = *jvp->GetValues();
        const auto reverse = *vjp->GetValues();
        RC_ASSERT(reverse.size() == seed.size());
        double left = 0.0;
        double right = 0.0;
        double scale = 0.0;
        for (std::size_t i = 0; i < forward.size(); ++i) {
            left += adjoint[i] * forward[i];
            scale += std::abs(adjoint[i] * forward[i]);
        }
        for (std::size_t k = 0; k < seed.size(); ++k) {
            right += reverse[k] * seed[k];
        }
        RC_ASSERT(std::abs(left - right) <= 1e-12 * (1.0 + scale));
        RC_ASSERT(jvp->GetSignature() == vjp->GetSignature());
    });
}

TEST_CASE("a pole's derivative on a line is its basis function", "[curve][sensitivity]")
{
    // C(t) = (1 - t) P0 + t P1: dC/dP0.x = 1 - t, and a weight on a
    // polynomial curve moves it as w0 = 1 + s would, dC/dw0 = (1 - t)(P0 - C).
    const auto c = Stereon::MakeBSplineCurve(1, {0.0, 0.0, 1.0, 1.0},
                                             {Point3{0.0, 0.0, 0.0}, Point3{2.0, 0.0, 0.0}});
    REQUIRE(c.has_value());
    ParameterSet parameters;
    const auto x0 = parameters.Add("x0", ParameterQuantity::Length);
    const auto w0 = parameters.Add("w0", ParameterQuantity::Dimensionless);
    REQUIRE((x0.has_value() && w0.has_value()));
    const std::vector<PoleBinding> bindings{{*x0, 0, PoleComponent::X},
                                            {*w0, 0, PoleComponent::Weight}};
    SensitivityContext sensitivity(std::move(parameters));
    const std::vector<double> t{0.25};

    REQUIRE(sensitivity.SetSeed(std::vector<double>{1.0, 0.0}).has_value());
    const auto dx = Stereon::EvalJvp(*c, bindings, t, sensitivity);
    REQUIRE(dx.has_value());
    REQUIRE((*dx->GetValues())[0] == 0.75);

    REQUIRE(sensitivity.SetSeed(std::vector<double>{0.0, 1.0}).has_value());
    const auto dw = Stereon::EvalJvp(*c, bindings, t, sensitivity);
    REQUIRE(dw.has_value());
    // C(0.25) = (0.5, 0, 0); (1 - t)(P0 - C) = 0.75 * (-0.5) = -0.375.
    REQUIRE((*dw->GetValues())[0] == -0.375);
}

TEST_CASE("invalid sensitivity input is refused", "[curve][sensitivity]")
{
    const auto c = Stereon::MakeBSplineCurve(1, {0.0, 0.0, 1.0, 1.0},
                                             {Point3{0.0, 0.0, 0.0}, Point3{2.0, 0.0, 0.0}});
    REQUIRE(c.has_value());
    ParameterSet parameters;
    const auto length = parameters.Add("x", ParameterQuantity::Length);
    const auto weight = parameters.Add("w", ParameterQuantity::Dimensionless);
    REQUIRE((length.has_value() && weight.has_value()));
    ParameterSet other;
    const auto foreign = other.Add("x", ParameterQuantity::Length);
    REQUIRE(foreign.has_value());
    SensitivityContext sensitivity(std::move(parameters));
    const std::vector<double> t{0.5};
    const auto refused = [&](const std::vector<PoleBinding>& bindings,
                             const std::vector<double>& at) {
        const auto result = Stereon::EvalJvp(*c, bindings, at, sensitivity);
        return !result && result.error().GetCode() == Stereon::ErrorCode::InvalidInput;
    };

    // No seed yet.
    REQUIRE(refused({{*length, 0, PoleComponent::X}}, t));
    REQUIRE(sensitivity.SetSeed(std::vector<double>{1.0, 1.0}).has_value());
    REQUIRE_FALSE(refused({{*length, 0, PoleComponent::X}}, t));
    REQUIRE(refused({{*foreign, 0, PoleComponent::X}}, t));
    REQUIRE(refused({{*length, 2, PoleComponent::X}}, t));
    REQUIRE(refused({{*length, 0, PoleComponent::Weight}}, t));
    REQUIRE(refused({{*weight, 0, PoleComponent::Y}}, t));
    REQUIRE(refused({{*length, 0, PoleComponent::X}}, {1.5}));
    REQUIRE(refused({{*length, 0, PoleComponent::X}}, {std::nan("")}));

    // An adjoint of the wrong length.
    REQUIRE(sensitivity.SetAdjoint(std::vector<double>{1.0, 0.0}).has_value());
    const auto vjp = Stereon::EvalVjp(*c, std::vector<PoleBinding>{{*length, 0, PoleComponent::X}},
                                      t, sensitivity);
    REQUIRE_FALSE(vjp.has_value());
}

TEST_CASE("Eval of a curve declares its forward and reverse products", "[curve][sensitivity]")
{
    const auto* eval = Stereon::FindCapability(Stereon::CurveCapabilities(), "Eval");
    REQUIRE(eval != nullptr);
    REQUIRE(eval->Supports(DerivativeProduct::Jvp));
    REQUIRE(eval->Supports(DerivativeProduct::Vjp));
    REQUIRE_FALSE(eval->Supports(DerivativeProduct::DenseJacobian));
    REQUIRE(eval->Correspondence.Contains(Stereon::TrackingConvention::FixedEdgeParameter));
}

TEST_CASE("a product whose arithmetic overflows is Indeterminate, not Valid",
          "[curve][sensitivity]")
{
    // Finite seeds, adjoints and references whose products overflow: the
    // audit of 2026-10-30 found these published as Valid with inf values.
    const BSplineCurve c{
        3,
        {0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0},
        {Point3{}, Point3{1.0, 1.0, 0.0}, Point3{2.0, -1.0, 0.0}, Point3{3.0, 0.0, 0.0}},
        {}};
    ParameterSet set;
    const auto id = set.Add("x", ParameterQuantity::Length, 1e300);
    REQUIRE(id.has_value());
    const std::vector<PoleBinding> bindings{{*id, 1, PoleComponent::X}};
    SensitivityContext sensitivity{std::move(set)};
    const std::vector<double> seed{1e300};
    REQUIRE(sensitivity.SetSeed(seed).has_value());
    const std::vector<double> t{0.5};
    const auto jvp = Stereon::EvalJvp(c, bindings, t, sensitivity);
    REQUIRE(jvp.has_value());
    const Stereon::SensitivityResult& forward = *jvp;
    REQUIRE(forward.GetStatus() == Stereon::DerivativeStatus::Indeterminate);
    REQUIRE_FALSE(forward.GetValues().has_value());
    REQUIRE(forward.GetReason().find("overflowed") != std::string_view::npos);

    const std::vector<double> adjoint{1e300, 0.0, 0.0};
    REQUIRE(sensitivity.SetAdjoint(adjoint).has_value());
    const auto vjp = Stereon::EvalVjp(c, bindings, t, sensitivity);
    REQUIRE(vjp.has_value());
    const Stereon::SensitivityResult& reverse = *vjp;
    REQUIRE(reverse.GetStatus() == Stereon::DerivativeStatus::Indeterminate);
    REQUIRE_FALSE(reverse.GetValues().has_value());
}
