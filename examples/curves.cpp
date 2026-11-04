// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// A cubic B-spline curve: built through the checked factory, evaluated,
// differentiated, and a point projected onto it (docs/examples.rst).

#include <cmath>
#include <cstdio>
#include <print>

import Stereon.Curve;

using namespace Stereon;

int main()
{
    // A clamped cubic on t in [0, 2]: five poles, the knots of two spans.
    const Result<BSplineCurve> curve = MakeBSplineCurve(
        3, {0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 2.0},
        {{0.0, 0.0, 0.0}, {1.0, 2.0, 0.0}, {3.0, 2.0, 1.0}, {4.0, 0.0, 1.0}, {6.0, 1.0, 0.0}});
    if (!curve) {
        std::println(stderr, "MakeBSplineCurve: {}", curve.error().GetMessage());
        return 1;
    }

    // The point and the first two derivatives at t = 0.5.
    const Point3 p = Eval(*curve, 0.5);
    const Vec3 d1 = Deriv(*curve, 0.5, 1);
    const Vec3 d2 = Deriv(*curve, 0.5, 2);
    std::println("C(0.5)   = ({:.6f}, {:.6f}, {:.6f})", p.X, p.Y, p.Z);
    std::println("C'(0.5)  = ({:.6f}, {:.6f}, {:.6f})", d1.X, d1.Y, d1.Z);
    std::println("C''(0.5) = ({:.6f}, {:.6f}, {:.6f})", d2.X, d2.Y, d2.Z);

    // The closest point of the curve to q. At an interior minimum the
    // residual q - C(t) is perpendicular to the tangent.
    const Point3 q{2.0, 3.0, 0.5};
    const double t = Project(*curve, q);
    const Vec3 residual = q - Eval(*curve, t);
    const double orthogonality = Dot(residual, Deriv(*curve, t, 1));
    std::println("closest parameter t = {:.12f}, distance {:.9f}", t, Length(residual));
    std::println("(q - C(t)) . C'(t)  = {:.3e}", orthogonality);
    return std::abs(orthogonality) < 1e-9 ? 0 : 1;
}
