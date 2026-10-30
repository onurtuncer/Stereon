// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The volume and area of a Shape, for the Boolean prototype's M11: what the
// OCCT differential tests compare with BRepGProp, and what the exact values
// of the fixtures check first. A test's measure, not a kernel API: it reads
// the faces' surfaces and their loops' pcurves only, and stops at the first
// face whose loops it cannot walk.
//
// Each face's loops are closed in its parameter plane (a cylinder's lateral
// face runs up and down its seam), so Green's theorem turns an integral over
// the face's region D into one around its loops:
//
//   area    = sum over D of |Su x Sv|              = loop integral of F dv,
//   volume  = 1/3 sum over faces of (S - c).N dA   = loop integral of G dv,
//
// where F and G integrate |Su x Sv| and (S - c).(Su x Sv) / 3 in u from the
// middle of the face's parameter box. A face's region is the one its loops
// enclose counterclockwise or clockwise, whichever the sign of its area
// says, and N is Su x Sv or its opposite as the face is reversed. Both
// integrals are composite Gauss-Legendre rules, run at two resolutions: the
// difference is returned as the measure's own error estimate, so a caller
// can require it to be far below what it compares.
//
// Include after `import Stereon.Boolean;`, with <cmath>, <cstddef>,
// <functional>, <numbers>, <stdexcept>, <string> and <vector> included before
// it (GCC 16.2 rejects a standard header first included after an import of a
// module that uses it).

#pragma once

namespace ShapeMeasures {

/// A shape's volume, the area of its faces and the length of its edges,
/// with how far the two resolutions of the quadrature disagreed on the
/// first two. The length is what a tolerance band along the edges scales.
struct Measures {
    double Volume = 0.0;
    double Area = 0.0;
    double EdgeLength = 0.0;
    double VolumeEstimate = 0.0;
    double AreaEstimate = 0.0;
};

/// The n-point Gauss-Legendre rule on [-1, 1], by Newton's method on the
/// Legendre polynomial from the Chebyshev guesses.
struct GaussRule {
    std::vector<double> Nodes;
    std::vector<double> Weights;

    explicit GaussRule(const int n)
    {
        for (int i = 1; i <= n; ++i) {
            double x = std::cos(std::numbers::pi * (i - 0.25) / (n + 0.5));
            double derivative = 0.0;
            for (int iteration = 0; iteration < 100; ++iteration) {
                double p0 = 1.0;
                double p1 = x;
                for (int k = 2; k <= n; ++k) {
                    const double p2 = ((2.0 * k - 1.0) * x * p1 - (k - 1.0) * p0) / k;
                    p0 = p1;
                    p1 = p2;
                }
                derivative = n * (x * p1 - p0) / (x * x - 1.0);
                const double step = p1 / derivative;
                x -= step;
                if (std::abs(step) <= 4.0 * std::numeric_limits<double>::epsilon()) {
                    break;
                }
            }
            Nodes.push_back(x);
            Weights.push_back(2.0 / ((1.0 - x * x) * derivative * derivative));
        }
    }
};

inline const GaussRule& Rule()
{
    static const GaussRule rule(16);
    return rule;
}

/// The integral of f over [a, b], cut into `pieces` equal pieces.
inline double Integrate(const std::function<double(double)>& f, const double a, const double b,
                        const int pieces)
{
    const GaussRule& rule = Rule();
    const double width = (b - a) / pieces;
    double sum = 0.0;
    for (int p = 0; p < pieces; ++p) {
        const double mid = a + (p + 0.5) * width;
        double piece = 0.0;
        for (std::size_t k = 0; k < rule.Nodes.size(); ++k) {
            piece += rule.Weights[k] * f(mid + 0.5 * width * rule.Nodes[k]);
        }
        sum += 0.5 * width * piece;
    }
    return sum;
}

/// The two loop integrals of one face at one resolution: its signed area
/// in the loops' sense, and the signed volume term, before the face's sense.
struct FaceIntegrals {
    double Area = 0.0;
    double Volume = 0.0;
};

inline FaceIntegrals Integrals(const Stereon::Shape& shape, const Stereon::FaceId face,
                               const Stereon::Point3& centre, const int pieces)
{
    using namespace Stereon;
    const Surface surface = SurfaceOf(shape, face);
    // The middle of the face's loops in the parameter plane, from which the
    // inner integrals start.
    double uLo = std::numeric_limits<double>::infinity();
    double uHi = -uLo;
    const auto loops = LoopsOf(shape, face);
    if (!loops) {
        throw std::runtime_error("ShapeMeasures: a face's loops cannot be walked");
    }
    for (const LoopId loop : *loops) {
        for (const CoedgeId coedge : *CoedgesOf(shape, loop)) {
            const std::optional<Curve2d> pcurve = PCurveOf(shape, coedge);
            if (!pcurve) {
                throw std::runtime_error(
                    "ShapeMeasures: " + Describe(shape, EdgeOf(shape, coedge)) +
                    " has a use without a pcurve");
            }
            const Interval range = RangeOf(shape, EdgeOf(shape, coedge));
            for (const double t : {range.GetLower(), range.GetUpper()}) {
                const double u = Eval(*pcurve, t).X;
                uLo = std::min(uLo, u);
                uHi = std::max(uHi, u);
            }
        }
    }
    const double u0 = 0.5 * (uLo + uHi);
    // The inner integrands at (s, v).
    const auto area = [&](const double s, const double v) {
        const Point2 uv{s, v};
        return Length(Cross(Deriv(surface, uv, 1, 0), Deriv(surface, uv, 0, 1)));
    };
    const auto volume = [&](const double s, const double v) {
        const Point2 uv{s, v};
        const Vec3 normal = Cross(Deriv(surface, uv, 1, 0), Deriv(surface, uv, 0, 1));
        return Dot(Eval(surface, uv) - centre, normal) / 3.0;
    };
    FaceIntegrals out;
    for (const LoopId loop : *loops) {
        for (const CoedgeId coedge : *CoedgesOf(shape, loop)) {
            const Curve2d pcurve = *PCurveOf(shape, coedge);
            const Interval range = RangeOf(shape, EdgeOf(shape, coedge));
            const double sense =
                StoredEdgeUseSense(shape, coedge) == EdgeUseSense::Same ? 1.0 : -1.0;
            const auto along = [&](const auto& inner) {
                return Integrate(
                    [&](const double t) {
                        const Point2 p = Eval(pcurve, t);
                        const double dv = Deriv(pcurve, t, 1).Y;
                        if (dv == 0.0) {
                            return 0.0;
                        }
                        const double v = p.Y;
                        return dv * Integrate([&](const double s) { return inner(s, v); }, u0, p.X,
                                              pieces);
                    },
                    range.GetLower(), range.GetUpper(), pieces);
            };
            out.Area += sense * along(area);
            out.Volume += sense * along(volume);
        }
    }
    return out;
}

/// The edge's length along one of its uses, the first of its radial ring:
/// |d/dt S(p(t))|, with S the use's surface and p its pcurve.
inline double LengthOf(const Stereon::Shape& shape, const Stereon::EdgeId edge, const int pieces)
{
    using namespace Stereon;
    const auto uses = RadialCoedgesOf(shape, edge);
    if (!uses || uses->begin() == uses->end()) {
        throw std::runtime_error("ShapeMeasures: " + Describe(shape, edge) + " has no use");
    }
    const CoedgeId use = *uses->begin();
    const std::optional<Curve2d> pcurve = PCurveOf(shape, use);
    if (!pcurve) {
        throw std::runtime_error("ShapeMeasures: " + Describe(shape, edge) +
                                 " has a use without a pcurve");
    }
    const Surface surface = SurfaceOf(shape, FaceOf(shape, LoopOf(shape, use)));
    const Interval range = RangeOf(shape, edge);
    return Integrate(
        [&](const double t) {
            const Point2 p = Eval(*pcurve, t);
            const Vec2 dp = Deriv(*pcurve, t, 1);
            return Length(dp.X * Deriv(surface, p, 1, 0) + dp.Y * Deriv(surface, p, 0, 1));
        },
        range.GetLower(), range.GetUpper(), pieces);
}

inline Measures MeasuresAt(const Stereon::Shape& shape, const int pieces)
{
    using namespace Stereon;
    Box3 box = Box3::Empty();
    for (const VertexId v : *AllEntities<EntityKind::Vertex>(shape)) {
        box = Hull(box, PositionOf(shape, v));
    }
    const Point3 centre =
        box.IsEmpty() ? Point3{} : Point3{Midpoint(box.X), Midpoint(box.Y), Midpoint(box.Z)};
    Measures m;
    for (const FaceId face : *AllEntities<EntityKind::Face>(shape)) {
        const FaceIntegrals f = Integrals(shape, face, centre, pieces);
        // The region is the one the loops enclose with positive area; the
        // face's normal is Su x Sv unless it is reversed.
        const double region = f.Area < 0.0 ? -1.0 : 1.0;
        const double normal = IsReversed(shape, face) ? -1.0 : 1.0;
        m.Area += region * f.Area;
        m.Volume += region * normal * f.Volume;
    }
    for (const EdgeId edge : *AllEntities<EntityKind::Edge>(shape)) {
        m.EdgeLength += LengthOf(shape, edge, pieces);
    }
    return m;
}

/// The shape's measures at 4 pieces per integral, with their difference
/// from 8 pieces as the estimate of their error.
inline Measures MeasuresOf(const Stereon::Shape& shape)
{
    Measures coarse = MeasuresAt(shape, 4);
    const Measures fine = MeasuresAt(shape, 8);
    coarse.VolumeEstimate = std::abs(fine.Volume - coarse.Volume);
    coarse.AreaEstimate = std::abs(fine.Area - coarse.Area);
    coarse.Volume = fine.Volume;
    coarse.Area = fine.Area;
    coarse.EdgeLength = fine.EdgeLength;
    return coarse;
}

} // namespace ShapeMeasures
