// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Curve lengths computed independently of both kernels' arc-length code,
// to judge them by. Include after `import Stereon.Surface;`.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace StereonTest {

/// The exact length of a clamped polyline (a degree-1 B-spline) from its
/// start to parameter t: plain arithmetic on its points.
inline double PolylineLength(const Stereon::BSplineCurve& c, const double t)
{
    double length = 0.0;
    for (std::size_t i = 0; i + 1 < c.Poles.size(); ++i) {
        const double a = c.Knots[i + 1];
        const double b = c.Knots[i + 2];
        if (t >= b) {
            length += Distance(c.Poles[i], c.Poles[i + 1]);
        } else if (t > a) {
            length += Distance(c.Poles[i], c.Poles[i + 1]) * (t - a) / (b - a);
        }
    }
    return length;
}

/// The length of any curve between a and b by an inscribed polyline of
/// 2^k evenly spaced chords, doubled until two values agree to 1e-12 and
/// then extrapolated once, since the polyline converges as 1/n^2. It
/// shares no code with ArcLength, which makes it a fair referee.
inline double PolylineReference(const Stereon::BSplineCurve& c, const double a, const double b)
{
    double previous = 0.0;
    for (int samples = 1 << 12;; samples *= 2) {
        double length = 0.0;
        Stereon::Point3 last = Eval(c, a);
        for (int i = 1; i <= samples; ++i) {
            const Stereon::Point3 next = Eval(c, std::min(b, a + (b - a) * i / samples));
            length += Distance(last, next);
            last = next;
        }
        if (std::abs(length - previous) <= 1e-12 * length || samples >= (1 << 22)) {
            return length + (length - previous) / 3.0;
        }
        previous = length;
    }
}

} // namespace StereonTest
