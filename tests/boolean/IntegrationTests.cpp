// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Intersection to Boolean, end to end: small solids whose faces meet in
// every way the plane–plane, plane–cylinder and cylinder–cylinder handlers
// and the face clip must get right (crossing, containing, sharing a face or
// an edge, tangent, oblique, skew), each through Union, Common and both
// Cuts. Nothing is checked against a count chosen by hand; every result
// must satisfy what any correct Boolean does:
//
// - it is a closed, consistently oriented solid: ValidateTopology at the
//   Solid level finds nothing but L05, which is never checked yet;
// - no two of its edges, and no two of its vertices, coincide;
// - its volume keeps the identities vol(A u B) + vol(A n B) = vol(A) +
//   vol(B) and vol(A - B) = vol(A) - vol(A n B), so no component is lost
//   or doubled;
// - a point is inside it exactly when the operation of its insides in A
//   and B says so, by winding numbers at a grid of points off every face.
//
// The pairs that fail a check today are listed with what is wrong, and a
// second test requires each of them still to fail: a fix shows up as that
// test failing, and the pair moves to the passing ones.

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

import Stereon.Boolean;

#include "BooleanCases.hpp"
#include "ShapeMeasures.hpp"

using BooleanCases::BoxFrom;
using BooleanCases::Must;
using BooleanCases::serial;
using Stereon::BooleanOperation;
using Stereon::EdgeId;
using Stereon::EntityKind;
using Stereon::Point3;
using Stereon::Shape;
using Stereon::Vec3;
using Stereon::VertexId;

namespace {

const Vec3 ex{1.0, 0.0, 0.0};
const Vec3 ey{0.0, 1.0, 0.0};
const Vec3 ez{0.0, 0.0, 1.0};

Shape Cylinder(const Point3& base, const Vec3& x, const Vec3& y, const double radius,
               const double height)
{
    return Must(MakeCylinderSolid(base, x, y, radius, height, serial), "a cylinder").Shape;
}

Shape OrientedBox(const Point3& origin, const Vec3& x, const Vec3& y, const Vec3& size)
{
    return Must(MakeBoxSolid(origin, x, y, size, serial), "a box").Shape;
}

/// A pair of solids and the box that holds both for the grid of points.
/// Known names what is wrong with the pair today; when it is empty, every
/// check must pass.
struct Pair {
    std::string Name;
    Shape A;
    Shape B;
    Point3 Lo;
    Point3 Hi;
    std::string Known;
};

std::vector<Pair> Pairs()
{
    const Shape box = BoxFrom(Point3{}, Point3{10.0, 10.0, 10.0});
    const Point3 lo{-1.0, -1.0, -1.0};
    const Point3 hi{11.0, 11.0, 11.0};
    // An exact Pythagorean frame turned about z, and one about y.
    const Vec3 turnedX{0.6, 0.8, 0.0};
    const Vec3 turnedY{-0.8, 0.6, 0.0};
    const Vec3 tiltY{-0.8, 0.0, 0.6}; // ey x tiltY = (0.6, 0, 0.8)
    std::vector<Pair> pairs;
    const auto add = [&](std::string name, Shape a, Shape b, const Point3& l, const Point3& h,
                         std::string known = {}) {
        pairs.push_back({std::move(name), std::move(a), std::move(b), l, h, std::move(known)});
    };
    // Boxes: planes against planes, and the clip of lines in faces.
    add("box corner", box, BoxFrom(Point3{5.0, 5.0, 5.0}, Point3{15.0, 15.0, 15.0}), lo,
        Point3{16.0, 16.0, 16.0});
    add("box inside", box, BoxFrom(Point3{2.0, 3.0, 4.0}, Point3{8.0, 7.0, 6.0}), lo, hi);
    add("box apart", box, BoxFrom(Point3{20.0, 0.0, 0.0}, Point3{30.0, 10.0, 10.0}), lo,
        Point3{31.0, 11.0, 11.0});
    add("box face to face", box, BoxFrom(Point3{10.0, 0.0, 0.0}, Point3{20.0, 10.0, 10.0}), lo,
        Point3{21.0, 11.0, 11.0});
    add("box on part of a face", box, BoxFrom(Point3{10.0, 2.0, 3.0}, Point3{15.0, 8.0, 7.0}), lo,
        Point3{16.0, 11.0, 11.0});
    add("box edge to edge", box, BoxFrom(Point3{10.0, 10.0, 0.0}, Point3{20.0, 20.0, 10.0}), lo,
        Point3{21.0, 21.0, 11.0});
    add("box sharing top and bottom", box, BoxFrom(Point3{5.0, 5.0, 0.0}, Point3{15.0, 15.0, 10.0}),
        lo, Point3{16.0, 16.0, 11.0},
        "refused: ImprintNetwork cannot close the network of the coplanar overlap");
    add("box turned through", box,
        OrientedBox(Point3{5.0, 2.0, -2.0}, turnedX, turnedY, Vec3{4.0, 3.0, 14.0}), lo,
        Point3{11.0, 11.0, 13.0});
    // A box and a cylinder: planes against a cylinder, circles, ellipses
    // and rulings clipped in faces.
    add("cylinder through", box, Cylinder(Point3{5.0, 5.0, -2.0}, ex, ey, 3.0, 14.0), lo,
        Point3{11.0, 11.0, 13.0});
    add("cylinder pocket", box, Cylinder(Point3{5.0, 5.0, 6.0}, ex, ey, 2.0, 8.0), lo,
        Point3{11.0, 11.0, 15.0});
    add("cylinder tangent to a side", box, Cylinder(Point3{13.0, 5.0, -2.0}, ex, ey, 3.0, 14.0), lo,
        Point3{17.0, 11.0, 13.0});
    add("cylinder about an edge", box, Cylinder(Point3{10.0, 10.0, -2.0}, ex, ey, 3.0, 14.0), lo,
        Point3{14.0, 14.0, 13.0});
    add("cylinder across", box, Cylinder(Point3{-2.0, 5.0, 5.0}, ey, ez, 2.0, 14.0),
        Point3{-3.0, -1.0, -1.0}, Point3{13.0, 11.0, 11.0});
    add("cylinder oblique", box, Cylinder(Point3{3.0, 5.0, -3.0}, ey, tiltY, 2.0, 16.0),
        Point3{-1.0, -1.0, -4.0}, Point3{15.0, 11.0, 13.0});
    // Two cylinders: coaxial, parallel, perpendicular and skew axes.
    const Shape upright = Cylinder(Point3{0.0, 0.0, -5.0}, ex, ey, 3.0, 10.0);
    const Point3 cl{-7.0, -4.0, -6.0};
    const Point3 ch{7.0, 4.0, 6.0};
    add("cylinders coaxial", upright, Cylinder(Point3{0.0, 0.0, -7.0}, ex, ey, 1.5, 14.0), cl,
        Point3{4.0, 4.0, 8.0});
    add("cylinders parallel", upright, Cylinder(Point3{2.0, 0.0, -7.0}, ex, ey, 3.0, 14.0), cl,
        Point3{6.0, 4.0, 8.0});
    add("cylinders parallel, caps shared", upright,
        Cylinder(Point3{2.0, 0.0, -5.0}, ex, ey, 3.0, 10.0), cl, Point3{6.0, 4.0, 6.0});
    add("cylinders perpendicular", upright, Cylinder(Point3{-6.0, 0.0, 0.0}, ey, ez, 1.5, 12.0), cl,
        ch);
    add("cylinders perpendicular, equal radii", upright,
        Cylinder(Point3{-6.0, 0.0, 0.0}, ey, ez, 3.0, 12.0), cl, ch,
        "refused: ImprintNetwork, at the node where the two curves of intersection cross");
    add("cylinders skew", upright, Cylinder(Point3{-6.0, 1.0, 1.0}, ey, ez, 1.5, 12.0), cl, ch);
    add("cylinders skew, oblique", upright,
        Cylinder(Point3{-5.0, -3.0, -2.0}, ez, Vec3{-0.8, 0.6, 0.0}, 1.0, 12.0), cl, ch);
    return pairs;
}

/// What the solid-level check finds besides L05.
std::string SolidFindings(const Shape& shape)
{
    const Stereon::ValidityReport report =
        ValidateTopology(shape, Stereon::ValidityLevel::Solid, serial);
    std::string text;
    for (const Stereon::ValidityFinding& f : report.GetFindings()) {
        if (f.Rule != Stereon::TopologyRule::SelfIntersectionUnchecked) {
            text += f.Message + "\n";
        }
    }
    return text;
}

/// The first use of an edge in its radial ring, if it has one.
std::optional<Stereon::CoedgeId> FirstUseOf(const Shape& shape, const EdgeId edge)
{
    const auto uses = RadialCoedgesOf(shape, edge);
    if (!uses || uses->begin() == uses->end()) {
        return std::nullopt;
    }
    return *uses->begin();
}

/// The middle of an edge, through its first use's pcurve on its face.
std::optional<Point3> MiddleOf(const Shape& shape, const EdgeId edge)
{
    const std::optional<Stereon::CoedgeId> use = FirstUseOf(shape, edge);
    if (!use) {
        return std::nullopt;
    }
    const std::optional<Stereon::Curve2d> pcurve = PCurveOf(shape, *use);
    if (!pcurve) {
        return std::nullopt;
    }
    const Stereon::Surface surface = SurfaceOf(shape, FaceOf(shape, LoopOf(shape, *use)));
    return Eval(surface, Eval(*pcurve, Midpoint(RangeOf(shape, edge))));
}

/// The shell a vertex's boundary belongs to, through its first used edge.
std::optional<Stereon::ShellId> ShellAt(const Shape& shape, const VertexId v)
{
    const auto edges = EdgesOf(shape, v);
    if (!edges) {
        return std::nullopt;
    }
    for (const EdgeId e : *edges) {
        if (const std::optional<Stereon::CoedgeId> use = FirstUseOf(shape, e)) {
            return ShellOf(shape, FaceOf(shape, LoopOf(shape, *use)));
        }
    }
    return std::nullopt;
}

/// Pairs of coincident vertices of one shell, and of edges between the same
/// vertices whose middles coincide, each within the sum of their
/// tolerances; and any edge no face uses.
std::string DuplicateVertices(const Shape& shape)
{
    std::string text;
    std::vector<VertexId> vertices;
    for (const VertexId v : *AllEntities<EntityKind::Vertex>(shape)) {
        vertices.push_back(v);
    }
    for (std::size_t i = 0; i < vertices.size(); ++i) {
        for (std::size_t j = i + 1; j < vertices.size(); ++j) {
            // Lumps touching at a point or along an edge have coincident
            // vertices, one in each shell: only a shell's own repeat one.
            if (ShellAt(shape, vertices[i]) != ShellAt(shape, vertices[j])) {
                continue;
            }
            const double d =
                Distance(PositionOf(shape, vertices[i]), PositionOf(shape, vertices[j]));
            if (d <= ToleranceOf(shape, vertices[i]) + ToleranceOf(shape, vertices[j])) {
                text += "vertices " + Describe(shape, vertices[i]) + " and " +
                        Describe(shape, vertices[j]) + " coincide\n";
            }
        }
    }
    return text;
}

std::string DuplicateEdges(const Shape& shape)
{
    std::string text;
    std::vector<EdgeId> edges;
    for (const EdgeId e : *AllEntities<EntityKind::Edge>(shape)) {
        if (!FirstUseOf(shape, e)) {
            text += "edge " + Describe(shape, e) + " has no use\n";
            continue;
        }
        edges.push_back(e);
    }
    for (std::size_t i = 0; i < edges.size(); ++i) {
        for (std::size_t j = i + 1; j < edges.size(); ++j) {
            const EdgeId e = edges[i];
            const EdgeId f = edges[j];
            const bool same =
                (StartOf(shape, e) == StartOf(shape, f) && EndOf(shape, e) == EndOf(shape, f)) ||
                (StartOf(shape, e) == EndOf(shape, f) && EndOf(shape, e) == StartOf(shape, f));
            if (!same) {
                continue;
            }
            const std::optional<Point3> m = MiddleOf(shape, e);
            const std::optional<Point3> n = MiddleOf(shape, f);
            if (m && n && Distance(*m, *n) <= ToleranceOf(shape, e) + ToleranceOf(shape, f)) {
                text +=
                    "edges " + Describe(shape, e) + " and " + Describe(shape, f) + " coincide\n";
            }
        }
    }
    return text;
}

std::string Duplicates(const Shape& shape)
{
    return DuplicateVertices(shape) + DuplicateEdges(shape);
}

/// A body's side of a point, from one evaluation of its winding number:
/// inside or outside, nothing where the number is not decided, and why
/// where it is refused. A refusal is not an undecided point: InvalidResult
/// means the faces do not add up to 0 or 1, which a closed, consistently
/// oriented body always does, so it is a problem of the body, not of the
/// point.
struct Side {
    std::optional<bool> Inside;
    std::optional<std::string> Refusal;
};

Side SideOf(const Shape& shape, const Point3& p)
{
    if (AllEntities<EntityKind::Face>(shape)->empty()) {
        return {false, std::nullopt};
    }
    const auto w = WindingNumberOf(shape, p, 8, serial);
    if (!w) {
        if (w.error().GetCode() == Stereon::ErrorCode::InvalidResult) {
            return {std::nullopt, std::string(w.error().GetMessage())};
        }
        return {};
    }
    return {w->IsDecided() ? std::optional<bool>{w->IsInside()} : std::nullopt, std::nullopt};
}

bool Expected(const BooleanOperation op, const bool a, const bool b, const bool reversed)
{
    switch (op) {
    case BooleanOperation::Union:
        return a || b;
    case BooleanOperation::Common:
        return a && b;
    case BooleanOperation::Cut:
        return reversed ? b && !a : a && !b;
    }
    return false;
}

struct Outcome {
    std::optional<Shape> Result;
    double Volume = 0.0;
    double Estimate = 0.0;
};

void CheckMembership(const Pair& pair, const Shape& shape, const BooleanOperation operation,
                     const bool reversed, const std::string& name,
                     std::vector<std::string>& problems)
{
    // Membership at a grid of points off the faces of both operands.
    constexpr int n = 4;
    const auto at = [&](const double lo, const double hi, const int s) {
        return lo + (hi - lo) * (s + 0.4142135623730951) / n;
    };
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                const Point3 p{at(pair.Lo.X, pair.Hi.X, i), at(pair.Lo.Y, pair.Hi.Y, j),
                               at(pair.Lo.Z, pair.Hi.Z, k)};
                const Side result = SideOf(shape, p);
                if (result.Refusal) {
                    problems.push_back(name +
                                       ": the winding number refuses it: " + *result.Refusal);
                    return;
                }
                const std::optional<bool> inA = SideOf(pair.A, p).Inside;
                const std::optional<bool> inB = SideOf(pair.B, p).Inside;
                const std::optional<bool>& inResult = result.Inside;
                if (inA && inB && inResult &&
                    *inResult != Expected(operation, *inA, *inB, reversed)) {
                    problems.push_back(name + " has (" + std::to_string(p.X) + ", " +
                                       std::to_string(p.Y) + ", " + std::to_string(p.Z) +
                                       ") on the wrong side");
                }
            }
        }
    }
}

/// Everything wrong with the four operations on a pair, one line each: a
/// refusal, a solid-level finding, a duplicate, a point on the wrong side,
/// a broken volume identity.
std::vector<std::string> Problems(const Pair& pair)
{
    std::vector<std::string> problems;
    const ShapeMeasures::Measures a = ShapeMeasures::MeasuresOf(pair.A);
    const ShapeMeasures::Measures b = ShapeMeasures::MeasuresOf(pair.B);
    struct Run {
        BooleanOperation Operation;
        bool Reversed;
    };
    const std::vector<Run> runs{{BooleanOperation::Union, false},
                                {BooleanOperation::Common, false},
                                {BooleanOperation::Cut, false},
                                {BooleanOperation::Cut, true}};
    std::vector<Outcome> outcomes;
    for (const Run& run : runs) {
        const std::string name =
            std::string(Stereon::ToString(run.Operation)) + (run.Reversed ? " B - A" : "");
        const auto result = run.Reversed ? Boolean(pair.B, pair.A, run.Operation, serial)
                                         : Boolean(pair.A, pair.B, run.Operation, serial);
        if (!result) {
            problems.push_back(name + " refused: " + std::string(result.error().GetMessage()));
            outcomes.push_back({});
            continue;
        }
        const Shape& shape = result->Shape;
        if (const std::string findings = SolidFindings(shape); !findings.empty()) {
            problems.push_back(name + " is not a closed solid: " + findings);
        }
        if (const std::string duplicates = Duplicates(shape); !duplicates.empty()) {
            problems.push_back(name + " has duplicates: " + duplicates);
        }
        const ShapeMeasures::Measures m = ShapeMeasures::MeasuresOf(shape);
        outcomes.push_back({shape, m.Volume, m.VolumeEstimate});

        CheckMembership(pair, shape, run.Operation, run.Reversed, name, problems);
    }
    // The volume identities, when every operation gave a result.
    if (std::ranges::all_of(outcomes, [](const Outcome& o) { return o.Result.has_value(); })) {
        const double united = outcomes[0].Volume;
        const double common = outcomes[1].Volume;
        const double cutAB = outcomes[2].Volume;
        const double cutBA = outcomes[3].Volume;
        double slack =
            1e-9 * std::max(1.0, a.Volume + b.Volume) + a.VolumeEstimate + b.VolumeEstimate;
        for (const Outcome& o : outcomes) {
            slack += o.Estimate;
        }
        const auto identity = [&](const char* what, const double lhs, const double rhs) {
            if (!(std::abs(lhs - rhs) <= slack)) {
                problems.push_back(std::string(what) + " breaks: " + std::to_string(lhs) +
                                   " against " + std::to_string(rhs));
            }
        };
        identity("vol(A u B) + vol(A n B) = vol(A) + vol(B)", united + common, a.Volume + b.Volume);
        identity("vol(A - B) = vol(A) - vol(A n B)", cutAB, a.Volume - common);
        identity("vol(B - A) = vol(B) - vol(A n B)", cutBA, b.Volume - common);
        identity("vol(A u B) = vol(A - B) + vol(B - A) + vol(A n B)", united,
                 cutAB + cutBA + common);
    }
    return problems;
}

std::string Lines(const std::vector<std::string>& lines)
{
    std::string text;
    for (const std::string& line : lines) {
        text += line + "\n";
    }
    return text;
}

} // namespace

TEST_CASE("Booleans of small solids are closed, unduplicated, and keep volume and membership",
          "[boolean][integration]")
{
    for (const Pair& pair : Pairs()) {
        if (!pair.Known.empty()) {
            continue;
        }
        const std::vector<std::string> problems = Problems(pair);
        INFO(pair.Name << ":\n" << Lines(problems));
        CHECK(problems.empty());
    }
}

TEST_CASE("the known Boolean defects and gaps are still there", "[boolean][integration][known]")
{
    // Each pair here fails some check today, for the reason it names. When
    // one passes, this fails: clear its Known so the test above holds it,
    // and update the ledger.
    for (const Pair& pair : Pairs()) {
        if (pair.Known.empty()) {
            continue;
        }
        const std::vector<std::string> problems = Problems(pair);
        INFO(pair.Name << " (" << pair.Known << "):\n" << Lines(problems));
        CHECK_FALSE(problems.empty());
    }
}
