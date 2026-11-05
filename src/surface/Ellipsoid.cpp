// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

// The contract anchor (docs/toolchain/gcc16-feature-probe.md).
#include <Stereon/Core/Contracts.hpp>

#include <array>
#include <cstddef>

module Stereon.Surface;

/// The ellipsoid's closest point (ADR-0048 §2), out of line for the reason
/// the factories are (Factories.cpp): inline non-template functions in a
/// module interface that instantiate another module's templates corrupt
/// GCC 16.2's module file.

namespace Stereon::Detail {

/// The root s >= s0 of G(s) = (n0 / (s + r0))² + (n1 / (s + r1))² +
/// (z2 / (s + 1))² - 1, bracketed and found by bisection as EllipseRoot
/// finds its own (Eberly, 2013): G decreases strictly there, so the
/// bisection needs no starting guess and ends at adjacent doubles.
[[nodiscard]] double EllipsoidRoot(const double r0, const double r1, const double z0,
                                   const double z1, const double z2, const double g)
{
    const double n0 = r0 * z0;
    const double n1 = r1 * z1;
    double s0 = z2 - 1.0;
    double s1 = g < 0.0 ? 0.0 : Hypot(Hypot(n0, n1), z2) - 1.0;
    double s = 0.0;
    for (int i = 0; i < EllipseRootSteps; ++i) {
        s = 0.5 * (s0 + s1);
        if (s == s0 || s == s1) {
            break; // the bracket holds two adjacent doubles
        }
        const double ratio0 = n0 / (s + r0);
        const double ratio1 = n1 / (s + r1);
        const double ratio2 = z2 / (s + 1.0);
        const double value = ratio0 * ratio0 + ratio1 * ratio1 + ratio2 * ratio2 - 1.0;
        if (value > 0.0) {
            s0 = s;
        } else if (value < 0.0) {
            s1 = s;
        } else {
            break;
        }
    }
    return s;
}

/// The point of the ellipsoid (x0/e0)² + (x1/e1)² + (x2/e2)² = 1 closest to
/// (y0, y1, y2), for e0 >= e1 >= e2 > 0 and a point in the first octant
/// (Eberly's case analysis). Off every coordinate plane the root of
/// EllipsoidRoot places it; on a plane through the two longer axes inside
/// the ellipse where the normals cross, it lies off that plane; otherwise
/// the problem is the ellipse's in the plane the point lies in.
[[nodiscard]] Point3 ClosestOnEllipsoidOctant(const double e0, const double e1, const double e2,
                                              const double y0, const double y1, const double y2)
{
    if (y2 > 0.0) {
        if (y1 > 0.0 && y0 > 0.0) {
            const double z0 = y0 / e0;
            const double z1 = y1 / e1;
            const double z2 = y2 / e2;
            const double g = z0 * z0 + z1 * z1 + z2 * z2 - 1.0;
            if (g == 0.0) {
                return {y0, y1, y2}; // on the ellipsoid
            }
            const double r0 = (e0 / e2) * (e0 / e2);
            const double r1 = (e1 / e2) * (e1 / e2);
            const double s = EllipsoidRoot(r0, r1, z0, z1, z2, g);
            return {r0 * y0 / (s + r0), r1 * y1 / (s + r1), y2 / (s + 1.0)};
        }
        if (y1 > 0.0) {
            const Point2 q = ClosestOnEllipseQuadrant(e1, e2, y1, y2);
            return {0.0, q.X, q.Y};
        }
        if (y0 > 0.0) {
            const Point2 q = ClosestOnEllipseQuadrant(e0, e2, y0, y2);
            return {q.X, 0.0, q.Y};
        }
        return {0.0, 0.0, e2};
    }
    const double denominator0 = e0 * e0 - e2 * e2;
    const double denominator1 = e1 * e1 - e2 * e2;
    const double numerator0 = e0 * y0;
    const double numerator1 = e1 * y1;
    if (numerator0 < denominator0 && numerator1 < denominator1) {
        const double x0 = numerator0 / denominator0;
        const double x1 = numerator1 / denominator1;
        const double rest = 1.0 - x0 * x0 - x1 * x1;
        if (rest > 0.0) {
            return {e0 * x0, e1 * x1, e2 * Sqrt(rest)};
        }
    }
    const Point2 q = ClosestOnEllipseQuadrant(e0, e1, y0, y1);
    return {q.X, q.Y, 0.0};
}

} // namespace Stereon::Detail

namespace Stereon {

/// The longitude and latitude of a point of the ellipsoid closest to p:
/// in the frame, the octant problem of ClosestOnEllipsoidOctant with the
/// semi-axes taken longest first, the point reflected back to p's octant
/// (ADR-0048 §2). p on the medial set, where two or more closest points
/// meet, gets one of them, the same one every time: the one in the
/// positive half of each coordinate p is zero in. With equal semi-axes
/// that is the sphere's Project, the centre included.
[[nodiscard]] Point2 Project(const Ellipsoid& ellipsoid, const Point3& p)
{
    Detail::RequireFinite(p, "Project(Ellipsoid)", "a non-finite point");
    const Vec3 offset = p - ellipsoid.Center;
    const std::array<double, 3> coordinate{
        Dot(offset, ellipsoid.XAxis), Dot(offset, ellipsoid.YAxis), Dot(offset, Axis(ellipsoid))};
    const std::array<double, 3> semi{ellipsoid.A, ellipsoid.B, ellipsoid.C};
    // Longest first, ties kept in frame order.
    std::array<std::size_t, 3> order{0, 1, 2};
    const auto settle = [&](const std::size_t a, const std::size_t b) {
        if (semi[order[b]] > semi[order[a]]) {
            const std::size_t kept = order[a];
            order[a] = order[b];
            order[b] = kept;
        }
    };
    settle(0, 1);
    settle(1, 2);
    settle(0, 1);
    const auto magnitude = [&](const std::size_t i) {
        return coordinate[i] < 0.0 ? -coordinate[i] : coordinate[i];
    };
    const Point3 q = Detail::ClosestOnEllipsoidOctant(semi[order[0]], semi[order[1]],
                                                      semi[order[2]], magnitude(order[0]),
                                                      magnitude(order[1]), magnitude(order[2]));
    const std::array<double, 3> sorted{q.X, q.Y, q.Z};
    std::array<double, 3> closest{};
    for (std::size_t k = 0; k < 3; ++k) {
        const std::size_t i = order[k];
        closest[i] = coordinate[i] < 0.0 ? -sorted[k] : sorted[k];
    }
    const double cosine = closest[0] / ellipsoid.A;
    const double sine = closest[1] / ellipsoid.B;
    return {WrappedAngle(sine, cosine), Atan2(closest[2] / ellipsoid.C, Hypot(cosine, sine))};
}

} // namespace Stereon
