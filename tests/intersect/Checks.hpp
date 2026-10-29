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
    // The last can round a step past the range.
    t[count] = upper;
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

/// |distance from the axis - radius|: how far a point is from the cylinder.
[[nodiscard]] inline double DistanceFrom(const Stereon::Cylinder& cylinder,
                                         const Stereon::Point3& point)
{
    const Stereon::Vec3 offset = point - cylinder.Origin;
    return std::abs(Length(Cross(offset, Normalized(Axis(cylinder)))) - cylinder.Radius);
}

/// The shape of a result made of branches between their ends: how many
/// connected curves the branches form when ends that coincide are joined,
/// and how many points join four or more ends, the nodes where branches
/// meet at a tangency.
struct Connectivity {
    std::size_t Curves = 0;
    std::size_t Nodes = 0;
};

[[nodiscard]] inline Connectivity ConnectivityOf(const Stereon::SurfaceSurfaceIntersection& result,
                                                 const double coincidence = 1e-7)
{
    std::vector<Stereon::Point3> clusters;
    std::vector<std::size_t> ends; // the cluster of each branch end, two per branch
    for (const Stereon::IntersectionBranch& branch : result.Branches) {
        for (const double t : {branch.Range.GetLower(), branch.Range.GetUpper()}) {
            const Stereon::Point3 p = Eval(branch.Curve3d, t);
            std::size_t cluster = clusters.size();
            for (std::size_t i = 0; i < clusters.size(); ++i) {
                if (Distance(clusters[i], p) <= coincidence) {
                    cluster = i;
                    break;
                }
            }
            if (cluster == clusters.size()) {
                clusters.push_back(p);
            }
            ends.push_back(cluster);
        }
    }
    // Union-find over the clusters, one edge per branch.
    std::vector<std::size_t> parent(clusters.size());
    for (std::size_t i = 0; i < parent.size(); ++i) {
        parent[i] = i;
    }
    const auto find = [&](std::size_t i) {
        while (parent[i] != i) {
            i = parent[i] = parent[parent[i]];
        }
        return i;
    };
    for (std::size_t b = 0; b + 1 < ends.size(); b += 2) {
        parent[find(ends[b])] = find(ends[b + 1]);
    }
    Connectivity c;
    std::vector<std::size_t> degree(clusters.size(), 0);
    for (const std::size_t e : ends) {
        ++degree[e];
    }
    for (std::size_t i = 0; i < clusters.size(); ++i) {
        if (find(i) == i) {
            ++c.Curves;
        }
        if (degree[i] >= 4) {
            ++c.Nodes;
        }
    }
    return c;
}

/// The error's message when the result is one, for INFO before a REQUIRE.
template <typename T>
[[nodiscard]] std::string ErrorText(const Stereon::Result<T>& result)
{
    return result ? std::string() : std::string(result.error().GetMessage());
}

/// A context whose log lines are appended to `lines`.
[[nodiscard]] inline Stereon::Context Logging(std::vector<std::string>& lines)
{
    return ctx.WithLogSink(
        [&lines](const Stereon::LogRecord& record) { lines.emplace_back(record.Message); });
}

} // namespace IntersectTest
