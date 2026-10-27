// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <span>

export module Stereon.Curve:BatchEvaluation;

import Stereon.Robust;
import :BSplineCurve;

/// Evaluation of a B-spline curve at many parameters at once, data-parallel
/// across the parameters (ADR-0007: std::simd, through Stereon::Simd). See
/// docs/nurbs.rst, *Batch and grid evaluation*.

export namespace Stereon {

/// points[k] = Eval(curve, parameters[k]) for every k, bit for bit: a run of
/// Simd::Width<double> consecutive parameters in one knot span is evaluated
/// in the lanes of one vector, by the same operations in the same order as
/// Eval performs on each, and any other parameter by Eval itself. The result
/// therefore does not depend on the SIMD width, the order of the
/// parameters, or how they fall into spans; sorted parameters, which share
/// spans, are the fast case. Preconditions as Eval's for every parameter,
/// and as many points as parameters.
void EvalAll(const BSplineCurve& curve, std::span<const double> parameters,
             std::span<Point3> points);
void EvalAll(const BSplineCurve2d& curve, std::span<const double> parameters,
             std::span<Point2> points);

} // namespace Stereon
