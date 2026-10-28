// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <cstddef>
#include <cstdint>
#include <span>

export module Stereon.Curve:PoleSensitivity;

import Stereon.Robust;
import :BSplineCurve;

/// The first derivative products of the kernel (ADR-0021, stage 2): the
/// points Eval gives a B-spline curve at fixed parameters, differentiated
/// with respect to design parameters that move its poles and weights. The
/// knots, the degree and the pole count are structural and stay fixed
/// (ADR-0021, section 4), so the map is smooth wherever the weights stay
/// positive, and every evaluation is Valid.
///
/// A parameter moves the coordinates and weights it is bound to by its
/// reference scale per unit of seed (ParameterSet), so the derivative of
/// the point at t with respect to parameter k is
///
///     dC/dp_k = ref_k * sum over k's bindings of R_i(t) e_c      (a coordinate)
///                                              or N_i(t) (P_i - C(t)) / W(t)  (a weight)
///
/// with N_i the basis functions, W = sum N_j w_j and R_i = N_i w_i / W. A
/// polynomial curve has every weight 1, and a weight bound to it is a
/// weight of 1 that may move. Samples are tracked at fixed parameter
/// (TrackingConvention::FixedEdgeParameter).

export namespace Stereon {

/// What a binding moves on a pole.
enum class PoleComponent : std::uint8_t {
    X,
    Y,
    Z,
    Weight,
};

/// A design parameter bound to one coordinate of one pole, a Length, or to
/// its weight, Dimensionless. A parameter may move several components (a
/// translation moves every pole) and a component may be moved by several
/// parameters; their contributions add.
struct PoleBinding {
    ParameterId Parameter;
    std::size_t Pole = 0;
    PoleComponent Component = PoleComponent::X;
};

/// The forward product J p' for the seed of `sensitivity`: the velocity of
/// Eval(curve, t) for each t of `parameters`, three values per point (x,
/// y, z) in order. InvalidInput for a binding whose parameter is not of the
/// context's set, whose quantity does not suit its component or whose pole
/// does not exist; a parameter that is not finite or lies outside the
/// domain; or a seed that was never set.
[[nodiscard]] Result<SensitivityResult> EvalJvp(const BSplineCurve& curve,
                                                std::span<const PoleBinding> bindings,
                                                std::span<const double> parameters,
                                                const SensitivityContext& sensitivity);

/// The reverse product J^T x' for the adjoint of `sensitivity`, a covector
/// of three components per point in the order EvalJvp gives them under the
/// Euclidean pairing: one value per parameter of the set, in its order.
/// Points are accumulated in their order, so the result does not depend on
/// anything but the input (ADR-0021, section 7). InvalidInput as EvalJvp,
/// and for an adjoint whose length is not three per point.
[[nodiscard]] Result<SensitivityResult> EvalVjp(const BSplineCurve& curve,
                                                std::span<const PoleBinding> bindings,
                                                std::span<const double> parameters,
                                                const SensitivityContext& sensitivity);

} // namespace Stereon
