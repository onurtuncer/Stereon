// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// What the analytic handler tests share: frames for the test surfaces, and
// the check that a branch's pcurves and its 3D curve agree on both surfaces
// along its range (ADR-0009's shared parameter). Included after
// `import Stereon.Intersect;`, with <algorithm>, <cmath>, <cstddef>,
// <string> and <vector> included before the import: a standard header
// first seen after the import conflicts with the module's declarations of
// it in GCC 16.2.

#pragma once

namespace IntersectTest {

inline const Stereon::Context ctx = Stereon::Context::Default();

/// A plane through `origin` with the given unit normal: its X axis is any
/// unit vector across the normal.
[[nodiscard]] inline Stereon::Plane PlaneWithNormal(const Stereon::Point3& origin,
                                                    const Stereon::Vec3& normal)
{
    const Stereon::Vec3 n = Normalized(normal);
    const Stereon::Vec3 seed =
        std::abs(n.X) < 0.9 ? Stereon::Vec3{1.0, 0.0, 0.0} : Stereon::Vec3{0.0, 1.0, 0.0};
    const Stereon::Vec3 x = Normalized(Cross(seed, n));
    return Stereon::Plane{origin, x, Cross(n, x)};
}

/// A right-handed cylinder about the given axis direction.
[[nodiscard]] inline Stereon::Cylinder CylinderAbout(const Stereon::Point3& origin,
                                                     const Stereon::Vec3& axis, const double radius)
{
    const Stereon::Vec3 a = Normalized(axis);
    const Stereon::Vec3 seed =
        std::abs(a.X) < 0.9 ? Stereon::Vec3{1.0, 0.0, 0.0} : Stereon::Vec3{0.0, 1.0, 0.0};
    const Stereon::Vec3 x = Normalized(Cross(seed, a));
    return Stereon::Cylinder{origin, radius, x, Cross(a, x), true};
}

/// Parameters across the branch's range: its ends and points between, or
/// a stretch of [-10, 10] when it is unbounded.
[[nodiscard]] inline std::vector<double> SamplesOf(const Stereon::Interval& range,
                                                   const std::size_t count)
{
    const double lower = std::isinf(range.GetLower()) ? -10.0 : range.GetLower();
    const double upper = std::isinf(range.GetUpper()) ? 10.0 : range.GetUpper();
    std::vector<double> t(count + 1);
    for (std::size_t k = 0; k <= count; ++k) {
        t[k] = lower + (upper - lower) * static_cast<double>(k) / static_cast<double>(count);
    }
    return t;
}

/// The largest distance, over `count` samples of the range, between the
/// branch's 3D point and the surfaces' points at its pcurves' parameters.
[[nodiscard]] inline double PcurveDeviation(const Stereon::IntersectionBranch& branch,
                                            const Stereon::Surface& a, const Stereon::Surface& b,
                                            const std::size_t count = 64)
{
    double worst = 0.0;
    for (const double t : SamplesOf(branch.Range, count)) {
        const Stereon::Point3 p = Eval(branch.Curve3d, t);
        worst = std::max(worst, Distance(p, Eval(a, Eval(branch.OnA, t))));
        worst = std::max(worst, Distance(p, Eval(b, Eval(branch.OnB, t))));
    }
    return worst;
}

/// A context whose log lines are appended to `lines`.
[[nodiscard]] inline Stereon::Context Logging(std::vector<std::string>& lines)
{
    return ctx.WithLogSink(
        [&lines](const Stereon::LogRecord& record) { lines.emplace_back(record.Message); });
}

} // namespace IntersectTest
