// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Dual.hpp"
#include "RapidCheck.hpp"
#include <cppad/cppad.hpp>

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <span>
#include <string>
#include <utility>
#include <vector>

import Stereon.Curve;
import Stereon.Surface;

// CppAD, the production AD backend (ADR-0045), on the scalar-generic
// evaluators of ADR-0021: poles, weights and transform parameters are
// CppAD::AD<double>, knots and parameters stay double. Each Jacobian CppAD
// tapes is checked against the closed form or against forward-mode duals,
// a second, independent derivative.

using AD = CppAD::AD<double>;
using Stereon::BSplineCurve;
using Stereon::Point3;
using Stereon::Vec3;
using StereonTest::Dual;

static_assert(Stereon::Scalar<AD>, "CppAD's active scalar is a Stereon Scalar");

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
    const int count = *rc::gen::inRange(1, 5);
    for (int k = 0; k < count; ++k) {
        t.push_back(d.GetLower() + (*rc::gen::inRange(0, 1001) / 1000.0) * Width(d));
    }
    return t;
}

/// The independent variables of a curve: x, y, z of every pole, then every
/// weight (1 for a polynomial curve, so that weights may move).
std::vector<double> VariablesOf(const BSplineCurve& c)
{
    std::vector<double> x;
    for (const Point3& pole : c.Poles) {
        x.insert(x.end(), {pole.X, pole.Y, pole.Z});
    }
    for (std::size_t i = 0; i < c.Poles.size(); ++i) {
        x.push_back(c.Weights.empty() ? 1.0 : c.Weights[i]);
    }
    return x;
}

/// The curve's points at `t` as a function of its poles and weights, taped.
CppAD::ADFun<double> TapeCurve(const BSplineCurve& c, const std::vector<double>& t)
{
    const std::size_t n = c.Poles.size();
    std::vector<AD> x(4 * n);
    const std::vector<double> x0 = VariablesOf(c);
    std::copy(x0.begin(), x0.end(), x.begin());
    CppAD::Independent(x);
    Stereon::BasicBSplineCurve<AD, 3> active{c.Degree, c.Knots, {}, {}};
    for (std::size_t i = 0; i < n; ++i) {
        active.Poles.push_back({x[3 * i], x[3 * i + 1], x[3 * i + 2]});
        active.Weights.push_back(x[3 * n + i]);
    }
    std::vector<AD> y;
    for (const double s : t) {
        const auto point = Eval(active, s);
        y.insert(y.end(), {point.X, point.Y, point.Z});
    }
    return CppAD::ADFun<double>(x, y);
}

/// The basis functions of the curve at t: the first pole they reach and
/// N_j for the p + 1 poles from there.
std::pair<std::size_t, std::vector<double>> Basis(const BSplineCurve& c, const double t)
{
    const std::size_t span = Stereon::FindSpan(c.Knots, c.Degree, c.Poles.size(), t);
    const auto p = static_cast<std::size_t>(c.Degree);
    std::vector<double> n(p + 1);
    Stereon::BasisFunctions(c.Knots, c.Degree, span, t, n);
    return {span - p, n};
}

bool Near(const double a, const double b, const double relative)
{
    return std::abs(a - b) <= relative * (1.0 + std::abs(b));
}

} // namespace

TEST_CASE("CppAD's Jacobian of a curve's points is the basis in poles and weights",
          "[ad][cppad][property]")
{
    // dC/dP_i = R_i = N_i w_i / W per coordinate, dC/dw_i = N_i (P_i - C) / W.
    rc::prop("for random polynomial and rational curves", [] {
        const BSplineCurve c = RandomCurve(*rc::gen::arbitrary<bool>());
        const std::vector<double> t = RandomParameters(c);
        CppAD::ADFun<double> f = TapeCurve(c, t);
        const std::size_t n = c.Poles.size();
        const std::vector<double> jacobian = f.Jacobian(VariablesOf(c));
        RC_ASSERT(jacobian.size() == 3 * t.size() * 4 * n);
        for (std::size_t k = 0; k < t.size(); ++k) {
            const auto [first, basis] = Basis(c, t[k]);
            const Point3 point = Eval(c, t[k]);
            double w = 0.0;
            for (std::size_t j = 0; j < basis.size(); ++j) {
                w += basis[j] * (c.Weights.empty() ? 1.0 : c.Weights[first + j]);
            }
            for (std::size_t i = 0; i < n; ++i) {
                const bool active = i >= first && i < first + basis.size();
                const double nI = active ? basis[i - first] : 0.0;
                const double wI = c.Weights.empty() ? 1.0 : c.Weights[i];
                const Vec3 dw = (nI / w) * (c.Poles[i] - point);
                const double dwExpected[] = {dw.X, dw.Y, dw.Z};
                for (std::size_t d = 0; d < 3; ++d) {
                    const std::size_t row = 3 * k + d;
                    for (std::size_t e = 0; e < 3; ++e) {
                        const double expected = d == e ? nI * wI / w : 0.0;
                        RC_ASSERT(Near(jacobian[row * 4 * n + 3 * i + e], expected, 1e-12));
                    }
                    RC_ASSERT(Near(jacobian[row * 4 * n + 3 * n + i], dwExpected[d], 1e-12));
                }
            }
        }
    });
}

TEST_CASE("CppAD's forward and reverse modes agree with forward-mode duals",
          "[ad][cppad][property]")
{
    rc::prop("J p' from the tape equals the duals' derivative, and so does J^T x'", [] {
        const BSplineCurve c = RandomCurve(*rc::gen::arbitrary<bool>());
        const std::vector<double> t = RandomParameters(c);
        CppAD::ADFun<double> f = TapeCurve(c, t);
        const std::size_t n = c.Poles.size();
        std::vector<double> seed(4 * n);
        for (double& s : seed) {
            s = *rc::gen::inRange(-100, 101) / 50.0;
        }
        (void)f.Forward(0, VariablesOf(c));
        const std::vector<double> forward = f.Forward(1, seed);

        Stereon::BasicBSplineCurve<Dual, 3> dual{c.Degree, c.Knots, {}, {}};
        for (std::size_t i = 0; i < n; ++i) {
            const Point3& pole = c.Poles[i];
            dual.Poles.push_back({Dual{pole.X, seed[3 * i]}, Dual{pole.Y, seed[3 * i + 1]},
                                  Dual{pole.Z, seed[3 * i + 2]}});
            dual.Weights.push_back(Dual{c.Weights.empty() ? 1.0 : c.Weights[i], seed[3 * n + i]});
        }
        for (std::size_t k = 0; k < t.size(); ++k) {
            const auto point = Eval(dual, t[k]);
            RC_ASSERT(Near(forward[3 * k], point.X.Deriv, 1e-11));
            RC_ASSERT(Near(forward[3 * k + 1], point.Y.Deriv, 1e-11));
            RC_ASSERT(Near(forward[3 * k + 2], point.Z.Deriv, 1e-11));
        }

        // <x', J p'> = <J^T x', p'>.
        std::vector<double> adjoint(3 * t.size());
        for (double& a : adjoint) {
            a = *rc::gen::inRange(-100, 101) / 50.0;
        }
        const std::vector<double> reverse = f.Reverse(1, adjoint);
        double left = 0.0;
        double right = 0.0;
        double scale = 0.0;
        for (std::size_t i = 0; i < forward.size(); ++i) {
            left += adjoint[i] * forward[i];
            scale += std::abs(adjoint[i] * forward[i]);
        }
        for (std::size_t i = 0; i < seed.size(); ++i) {
            right += reverse[i] * seed[i];
        }
        RC_ASSERT(std::abs(left - right) <= 1e-12 * (1.0 + scale));
    });
}

TEST_CASE("the kernel's EvalJvp and EvalVjp equal CppAD's forward and reverse modes",
          "[ad][cppad][property]")
{
    // One design parameter per variable of the tape, reference 1: parameter
    // 3i + c moves coordinate c of pole i, parameter 3n + i its weight.
    rc::prop("for random curves, seeds and adjoints", [] {
        const BSplineCurve c = RandomCurve(*rc::gen::arbitrary<bool>());
        const std::vector<double> t = RandomParameters(c);
        CppAD::ADFun<double> f = TapeCurve(c, t);
        const std::size_t n = c.Poles.size();
        Stereon::ParameterSet parameters;
        std::vector<Stereon::PoleBinding> bindings;
        constexpr Stereon::PoleComponent Coordinates[] = {
            Stereon::PoleComponent::X, Stereon::PoleComponent::Y, Stereon::PoleComponent::Z};
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t d = 0; d < 3; ++d) {
                const auto id = parameters.Add("p" + std::to_string(3 * i + d),
                                               Stereon::ParameterQuantity::Length);
                RC_ASSERT(id.has_value());
                bindings.push_back({*id, i, Coordinates[d]});
            }
        }
        for (std::size_t i = 0; i < n; ++i) {
            const auto id =
                parameters.Add("w" + std::to_string(i), Stereon::ParameterQuantity::Dimensionless);
            RC_ASSERT(id.has_value());
            bindings.push_back({*id, i, Stereon::PoleComponent::Weight});
        }
        std::vector<double> seed(4 * n);
        for (double& s : seed) {
            s = *rc::gen::inRange(-100, 101) / 50.0;
        }
        std::vector<double> adjoint(3 * t.size());
        for (double& a : adjoint) {
            a = *rc::gen::inRange(-100, 101) / 50.0;
        }
        Stereon::SensitivityContext sensitivity(std::move(parameters));
        RC_ASSERT(sensitivity.SetSeed(seed).has_value());
        RC_ASSERT(sensitivity.SetAdjoint(adjoint).has_value());

        (void)f.Forward(0, VariablesOf(c));
        const std::vector<double> forward = f.Forward(1, seed);
        const auto jvp = Stereon::EvalJvp(c, bindings, t, sensitivity);
        RC_ASSERT(jvp.has_value());
        const auto kernelForward = *jvp->GetValues();
        RC_ASSERT(kernelForward.size() == forward.size());
        for (std::size_t i = 0; i < forward.size(); ++i) {
            RC_ASSERT(Near(kernelForward[i], forward[i], 1e-11));
        }

        const std::vector<double> reverse = f.Reverse(1, adjoint);
        const auto vjp = Stereon::EvalVjp(c, bindings, t, sensitivity);
        RC_ASSERT(vjp.has_value());
        const auto kernelReverse = *vjp->GetValues();
        RC_ASSERT(kernelReverse.size() == reverse.size());
        for (std::size_t i = 0; i < reverse.size(); ++i) {
            RC_ASSERT(Near(kernelReverse[i], reverse[i], 1e-11));
        }
    });
}

TEST_CASE("one tape evaluates the curve for any poles and weights", "[ad][cppad][property]")
{
    // Eval's branches depend on the knots and parameters only, which stay
    // double, so the tape records no comparison of active values and holds
    // for every pole and positive weight.
    rc::prop("the taped points equal Eval's for new poles and weights", [] {
        const BSplineCurve c = RandomCurve(*rc::gen::arbitrary<bool>());
        const std::vector<double> t = RandomParameters(c);
        CppAD::ADFun<double> f = TapeCurve(c, t);
        RC_ASSERT(f.compare_change_number() == 0);
        BSplineCurve moved = c;
        moved.Weights.assign(c.Poles.size(), 1.0);
        for (std::size_t i = 0; i < c.Poles.size(); ++i) {
            moved.Poles[i] = {*rc::gen::inRange(-50, 50) / 10.0, *rc::gen::inRange(-50, 50) / 10.0,
                              *rc::gen::inRange(-50, 50) / 10.0};
            moved.Weights[i] = *rc::gen::inRange(2, 40) / 10.0;
        }
        const std::vector<double> y = f.Forward(0, VariablesOf(moved));
        for (std::size_t k = 0; k < t.size(); ++k) {
            const Point3 expected = Eval(moved, t[k]);
            RC_ASSERT(Near(y[3 * k], expected.X, 1e-12));
            RC_ASSERT(Near(y[3 * k + 1], expected.Y, 1e-12));
            RC_ASSERT(Near(y[3 * k + 2], expected.Z, 1e-12));
        }
    });
}

TEST_CASE("CppAD differentiates a point through a composed rigid motion", "[ad][cppad]")
{
    // x -> R(theta)(x - c) + c + d about the unit axis k through c: the
    // derivative in theta is k x R(theta)(x - c), and in d the identity.
    const Vec3 axis{0.0, 0.6, 0.8};
    const Point3 centre{1.0, -2.0, 0.5};
    const Point3 x{3.0, 1.0, -1.0};
    const double theta = std::numbers::pi / 5.0;
    std::vector<AD> parameters{theta, 0.1, 0.2, 0.3};
    CppAD::Independent(parameters);
    const Stereon::BasicVec3<AD> k{axis.X, axis.Y, axis.Z};
    const Stereon::BasicPoint3<AD> about{centre.X, centre.Y, centre.Z};
    const auto motion = Stereon::MakeTranslation(
                            Stereon::BasicVec3<AD>{parameters[1], parameters[2], parameters[3]}) *
                        Stereon::MakeRotation(about, k, parameters[0]);
    const auto moved = Stereon::Transformed(Stereon::BasicPoint3<AD>{x.X, x.Y, x.Z}, motion);
    std::vector<AD> y{moved.X, moved.Y, moved.Z};
    CppAD::ADFun<double> f(parameters, y);
    const std::vector<double> jacobian = f.Jacobian(std::vector<double>{theta, 0.1, 0.2, 0.3});

    const auto rotation = Stereon::MakeRotation(axis, theta);
    const Vec3 dTheta = Cross(axis, Stereon::Transformed(x - centre, rotation));
    const double expected[3] = {dTheta.X, dTheta.Y, dTheta.Z};
    for (std::size_t d = 0; d < 3; ++d) {
        REQUIRE(Near(jacobian[4 * d], expected[d], 1e-14));
        for (std::size_t e = 0; e < 3; ++e) {
            REQUIRE(jacobian[4 * d + 1 + e] == (d == e ? 1.0 : 0.0));
        }
    }
}

TEST_CASE("CppAD's Jacobian of a surface point is the tensor-product basis", "[ad][cppad]")
{
    // A rational biquadratic patch: dS/dP_ij = N_i(u) M_j(v) w_ij / W.
    const std::vector<double> knots{0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0};
    std::vector<Point3> poles;
    std::vector<double> weights;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            poles.push_back({1.0 * i, 1.0 * j, 0.25 * ((i * 7 + j * 3) % 5)});
            weights.push_back(1.0 + 0.1 * ((i + 2 * j) % 4));
        }
    }
    const auto surface = Stereon::MakeBSplineSurface(2, 2, knots, knots, poles, weights);
    REQUIRE(surface.has_value());
    const double u = 0.3;
    const double v = 0.8;

    std::vector<AD> x;
    for (const Point3& pole : poles) {
        x.insert(x.end(), {pole.X, pole.Y, pole.Z});
    }
    CppAD::Independent(x);
    Stereon::BasicBSplineSurface<AD> active{2, 2, knots, knots, {}, {}};
    for (std::size_t i = 0; i < poles.size(); ++i) {
        active.Poles.push_back({x[3 * i], x[3 * i + 1], x[3 * i + 2]});
        active.Weights.push_back(AD{weights[i]});
    }
    const auto point = Eval(active, Stereon::Point2{u, v});
    std::vector<AD> y{point.X, point.Y, point.Z};
    CppAD::ADFun<double> f(x, y);
    std::vector<double> x0;
    for (const Point3& pole : poles) {
        x0.insert(x0.end(), {pole.X, pole.Y, pole.Z});
    }
    const std::vector<double> jacobian = f.Jacobian(x0);

    const auto basisOf = [&](const double s) {
        const std::size_t span = Stereon::FindSpan(knots, 2, 4, s);
        std::array<double, 3> n{};
        Stereon::BasisFunctions(knots, 2, span, s, n);
        std::array<double, 4> full{};
        for (std::size_t j = 0; j < 3; ++j) {
            full[span - 2 + j] = n[j];
        }
        return full;
    };
    const auto nu = basisOf(u);
    const auto nv = basisOf(v);
    double w = 0.0;
    for (std::size_t i = 0; i < 4; ++i) {
        for (std::size_t j = 0; j < 4; ++j) {
            w += nu[i] * nv[j] * weights[4 * i + j];
        }
    }
    const std::size_t columns = 3 * poles.size();
    for (std::size_t i = 0; i < 4; ++i) {
        for (std::size_t j = 0; j < 4; ++j) {
            const double r = nu[i] * nv[j] * weights[4 * i + j] / w;
            for (std::size_t d = 0; d < 3; ++d) {
                REQUIRE(Near(jacobian[d * columns + 3 * (4 * i + j) + d], r, 1e-13));
            }
        }
    }
}
