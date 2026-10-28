// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <span>

export module Stereon.Surface:GridEvaluation;

import Stereon.Curve;
import :BSplineSurface;

/// Evaluation of a B-spline surface on the grid of every pair of given u and
/// v parameters (docs/nurbs.rst, *Batch and grid evaluation*).

export namespace Stereon {

/// points[i * v.size() + j] = Eval(surface, {u[i], v[j]}) for every i and j,
/// bit for bit, row by row along u. The span and basis functions of each u
/// and each v are computed once for the grid rather than once per point,
/// and within a row a run of Simd::Width<double> consecutive v parameters in
/// one knot span is evaluated in the lanes of one vector, by the operations
/// of Eval in their order; any other point by those operations one at a
/// time. The result therefore does not depend on the SIMD width or the
/// order of the parameters; sorted v parameters are the fast case.
/// Preconditions as Eval's for every pair, and u.size() * v.size() points.
void EvalGrid(const BSplineSurface& surface, std::span<const double> u, std::span<const double> v,
              std::span<Point3> points);

} // namespace Stereon
