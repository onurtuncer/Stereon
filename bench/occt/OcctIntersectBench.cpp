// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The OpenCASCADE side of B2.1 to B2.3 (docs/benchmarks.md, Phase 2), on the
// cases of bench/IntersectCases.hpp, each the twin of a row of
// bench/IntersectBench.cpp, built from the same Stereon geometry:
//
//   B2.1  Geom2dAPI_InterCurveCurve on the planar curves, over the same
//         domains, at a tolerance of 1e-9 as tests/occt/PlanarDifferentialTests.cpp
//         asks; GeomAPI_IntCS for the lines through the NURBS sphere and
//         cylinder;
//   B2.2  GeomAPI_IntSS's engine, GeomInt_IntSS with GeomAPI_IntSS's
//         arguments (GeomAPI_IntSS.lxx), which also returns the isolated
//         points GeomAPI_IntSS does not expose, at the case's tolerance; a
//         cone's or cylinder's box is a Geom_RectangularTrimmedSurface and an
//         ellipsoid the unit sphere's exact rational B-spline under its affine
//         map (GeomAPI_IntSS throws on OCCT 8's GeomEval_EllipsoidSurface);
//   B2.3  IntAna_QuadQuadGeo, the analytic service, at the case's tolerance
//         (an angular tolerance of Precision::Angular() where it takes one).
//
// An exception or IsDone() false leaves a case unresolved; IntAna_QuadQuadGeo's
// NoGeometricSolution too. Both libraries' answers are scored by the same
// rules (IntersectCases.hpp). OCCT is LGPL-2.1 and linked into this
// executable only (ADR-0008).

#include "Json.hpp"
#include <ElCLib.hxx>
#include <Geom2dAPI_InterCurveCurve.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom2d_Circle.hxx>
#include <Geom2d_Curve.hxx>
#include <Geom2d_Ellipse.hxx>
#include <Geom2d_Hyperbola.hxx>
#include <Geom2d_Line.hxx>
#include <Geom2d_Parabola.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <GeomAPI_IntCS.hxx>
#include <GeomConvert.hxx>
#include <GeomInt_IntSS.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_ConicalSurface.hxx>
#include <Geom_Curve.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_Line.hxx>
#include <Geom_Plane.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <Geom_SphericalSurface.hxx>
#include <Geom_Surface.hxx>
#include <Geom_ToroidalSurface.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <IntAna_QuadQuadGeo.hxx>
#include <IntAna_ResultType.hxx>
#include <NCollection_Array1.hxx>
#include <NCollection_Array2.hxx>
#include <Precision.hxx>
#include <benchmark/benchmark.h>
#include <gp_Ax22d.hxx>
#include <gp_Ax3.hxx>
#include <gp_Circ.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Dir.hxx>
#include <gp_Dir2d.hxx>
#include <gp_Elips.hxx>
#include <gp_Lin.hxx>
#include <gp_Pln.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Sphere.hxx>
#include <gp_Vec.hxx>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <limits>
#include <map>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

import Stereon.Intersect;

#include "../../tests/occt/OcctConversion.hpp"
#include "../IntersectCases.hpp"

namespace {

namespace B2 = StereonBench::B2;
using Stereon::Point2;
using Stereon::Point3;
using Stereon::Vec2;
using Stereon::Vec3;
using StereonTest::ToOcct;

/// Geom2dAPI_InterCurveCurve's tolerance (tests/occt/PlanarDifferentialTests.cpp).
constexpr double PlanarTolerance = 1e-9;

// --- Geometry -----------------------------------------------------------------

gp_Pnt2d ToOcct(const Point2& p)
{
    return {p.X, p.Y};
}
gp_Dir2d ToDir(const Vec2& v)
{
    return {v.X, v.Y};
}
gp_Dir ToDir(const Vec3& v)
{
    return {v.X, v.Y, v.Z};
}
gp_Ax22d FrameOf(const Point2& origin, const Vec2& x, const Vec2& y)
{
    return gp_Ax22d(ToOcct(origin), ToDir(x), ToDir(y));
}
gp_Ax3 FrameOf(const Point3& origin, const Vec3& x, const Vec3& y)
{
    return gp_Ax3(ToOcct(origin), ToDir(Cross(x, y)), ToDir(x));
}

occ::handle<Geom2d_Curve> SplineOf(const Stereon::BSplineCurve2d& c)
{
    const int count = static_cast<int>(c.Poles.size());
    NCollection_Array1<gp_Pnt2d> poles(1, count);
    NCollection_Array1<double> weights(1, count);
    for (int i = 0; i < count; ++i) {
        const auto k = static_cast<std::size_t>(i);
        poles(i + 1) = ToOcct(c.Poles[k]);
        weights(i + 1) = c.Weights.empty() ? 1.0 : c.Weights[k];
    }
    const StereonTest::OcctKnots knots = StereonTest::ToOcctKnots(c.Knots);
    return new Geom2d_BSplineCurve(poles, weights, knots.Values, knots.Multiplicities, c.Degree);
}

/// A plane curve as OCCT's, trimmed to the domain when it is bounded; OCCT's
/// line runs at unit speed, so a line's domain is scaled by its speed.
occ::handle<Geom2d_Curve> PlanarCurveOf(const Stereon::Curve2d& curve,
                                        const Stereon::Interval& domain)
{
    occ::handle<Geom2d_Curve> made;
    double scale = 1.0;
    if (const auto* l = curve.Target<Stereon::Line2d>()) {
        made = new Geom2d_Line(ToOcct(l->Origin), ToDir(l->Direction));
        scale = Length(l->Direction);
    } else if (const auto* c = curve.Target<Stereon::Circle2d>()) {
        made = new Geom2d_Circle(FrameOf(c->Center, c->XAxis, c->YAxis), c->Radius);
    } else if (const auto* e = curve.Target<Stereon::Ellipse2d>()) {
        made = new Geom2d_Ellipse(FrameOf(e->Center, e->XAxis, e->YAxis), e->MajorRadius,
                                  e->MinorRadius);
    } else if (const auto* h = curve.Target<Stereon::Hyperbola2d>()) {
        made = new Geom2d_Hyperbola(FrameOf(h->Center, h->XAxis, h->YAxis), h->MajorRadius,
                                    h->MinorRadius);
    } else if (const auto* p = curve.Target<Stereon::Parabola2d>()) {
        made = new Geom2d_Parabola(FrameOf(p->Vertex, p->XAxis, p->YAxis), p->FocalLength);
    } else {
        const auto* s = curve.Target<Stereon::BSplineCurve2d>();
        if (s == nullptr) {
            throw std::runtime_error("a plane curve OCCT's twin does not convert");
        }
        made = SplineOf(*s);
    }
    if (!domain.IsBounded()) {
        return made;
    }
    return new Geom2d_TrimmedCurve(made, domain.GetLower() * scale, domain.GetUpper() * scale);
}

/// An ellipsoid as OCCT's exact rational B-spline: the unit sphere's
/// (GeomConvert), its poles mapped by the ellipsoid's affine map.
occ::handle<Geom_Surface> EllipsoidOf(const Stereon::Ellipsoid& e)
{
    const occ::handle<Geom_BSplineSurface> unit =
        GeomConvert::SurfaceToBSplineSurface(new Geom_SphericalSurface(
            gp_Ax3(gp_Pnt(0.0, 0.0, 0.0), gp_Dir(0.0, 0.0, 1.0), gp_Dir(1.0, 0.0, 0.0)), 1.0));
    const Vec3 z = Cross(e.XAxis, e.YAxis);
    for (int i = 1; i <= unit->NbUPoles(); ++i) {
        for (int j = 1; j <= unit->NbVPoles(); ++j) {
            const gp_Pnt p = unit->Pole(i, j);
            unit->SetPole(i, j,
                          ToOcct(e.Center + (e.A * p.X()) * e.XAxis + (e.B * p.Y()) * e.YAxis +
                                 (e.C * p.Z()) * z));
        }
    }
    return unit;
}

/// A surface as OCCT's, trimmed to its box when the box is bounded.
occ::handle<Geom_Surface> SurfaceOf(const Stereon::Surface& s, const Stereon::Box2& box)
{
    occ::handle<Geom_Surface> made;
    if (const auto* p = s.Target<Stereon::Plane>()) {
        made = new Geom_Plane(FrameOf(p->Origin, p->XAxis, p->YAxis));
    } else if (const auto* sp = s.Target<Stereon::Sphere>()) {
        made = new Geom_SphericalSurface(FrameOf(sp->Center, sp->XAxis, sp->YAxis), sp->Radius);
    } else if (const auto* cy = s.Target<Stereon::Cylinder>()) {
        made = new Geom_CylindricalSurface(FrameOf(cy->Origin, cy->XAxis, cy->YAxis), cy->Radius);
    } else if (const auto* co = s.Target<Stereon::Cone>()) {
        made = new Geom_ConicalSurface(FrameOf(co->Origin, co->XAxis, co->YAxis), co->SemiAngle,
                                       co->Radius);
    } else if (const auto* t = s.Target<Stereon::Torus>()) {
        made = new Geom_ToroidalSurface(FrameOf(t->Center, t->XAxis, t->YAxis), t->MajorRadius,
                                        t->MinorRadius);
    } else if (const auto* e = s.Target<Stereon::Ellipsoid>()) {
        made = EllipsoidOf(*e);
    } else {
        const auto* spline = s.Target<Stereon::BSplineSurface>();
        if (spline == nullptr) {
            throw std::runtime_error("a surface OCCT's twin does not convert");
        }
        made = ToOcct(*spline);
    }
    if (!box.X.IsBounded() || !box.Y.IsBounded()) {
        return made;
    }
    return new Geom_RectangularTrimmedSurface(made, box.X.GetLower(), box.X.GetUpper(),
                                              box.Y.GetLower(), box.Y.GetUpper());
}

Point3 PointOf(const gp_Pnt& p)
{
    return {p.X(), p.Y(), p.Z()};
}

// --- B2.1 ---------------------------------------------------------------------

struct PlanarPair {
    occ::handle<Geom2d_Curve> A;
    occ::handle<Geom2d_Curve> B;
};

/// OCCT's points of a pair, or none if it threw.
std::optional<std::vector<Point2>> PlanarPoints(const PlanarPair& pair)
{
    try {
        const Geom2dAPI_InterCurveCurve occt(pair.A, pair.B, PlanarTolerance);
        std::vector<Point2> points;
        for (int i = 1; i <= occt.NbPoints(); ++i) {
            const gp_Pnt2d p = occt.Point(i);
            points.push_back(Point2{p.X(), p.Y()});
        }
        return points;
    } catch (const std::exception&) {
        return std::nullopt;
    }
}

void OcctCurveCurve(benchmark::State& state, const bool splines)
{
    static std::map<bool, std::pair<std::vector<PlanarPair>, B2::Score>> rows;
    if (!rows.contains(splines)) {
        const std::vector<B2::PlanarCase> cases = B2::PlanarCases(splines);
        std::pair<std::vector<PlanarPair>, B2::Score> row;
        for (const B2::PlanarCase& c : cases) {
            row.first.push_back({PlanarCurveOf(c.A, c.DomainA), PlanarCurveOf(c.B, c.DomainB)});
            row.second.Cases += 1.0;
            const std::optional<std::vector<Point2>> points = PlanarPoints(row.first.back());
            if (!points) {
                row.second.Unresolved += 1.0;
                row.second.Expected += static_cast<double>(c.Roots.size());
                continue;
            }
            B2::ScoreRoots(c.Roots, c.Within, *points, row.second);
        }
        rows.emplace(splines, std::move(row));
    }
    const auto& [pairs, score] = rows.at(splines);
    for (auto _ : state) {
        for (const PlanarPair& pair : pairs) {
            try {
                Geom2dAPI_InterCurveCurve occt(pair.A, pair.B, PlanarTolerance);
                benchmark::DoNotOptimize(occt);
            } catch (const std::exception&) {
                benchmark::ClobberMemory();
            }
        }
    }
    state.SetItemsProcessed(state.iterations() * static_cast<std::int64_t>(pairs.size()));
    B2::Report(state, score);
}

void OcctCurveSurface(benchmark::State& state, const bool sphere)
{
    static std::map<bool, std::pair<occ::handle<Geom_Surface>, B2::Score>> rows;
    static std::map<bool, std::vector<occ::handle<Geom_Curve>>> lines;
    if (!rows.contains(sphere)) {
        const Stereon::BSplineSurface spline =
            B2::CorpusSpline(sphere ? "nurbs-ssi/sphere-plane/circle-origin-x1"
                                    : "nurbs-ssi/cylinder-plane/circle-origin-x1",
                             0);
        const std::vector<B2::LineCase> cases =
            sphere ? B2::LinesThroughSphere() : B2::LinesThroughCylinder();
        std::pair<occ::handle<Geom_Surface>, B2::Score> row{ToOcct(spline), {}};
        for (const B2::LineCase& c : cases) {
            const occ::handle<Geom_Curve> line =
                new Geom_Line(gp_Lin(ToOcct(c.Line.Origin), ToDir(c.Line.Direction)));
            lines[sphere].push_back(line);
            row.second.Cases += 1.0;
            const GeomAPI_IntCS occt(line, row.first);
            if (!occt.IsDone()) {
                row.second.Unresolved += 1.0;
                row.second.Expected += static_cast<double>(c.Roots.size());
                continue;
            }
            std::vector<Point3> points;
            for (int i = 1; i <= occt.NbPoints(); ++i) {
                points.push_back(PointOf(occt.Point(i)));
            }
            B2::ScoreRoots(c.Roots, std::vector<double>(c.Roots.size(), B2::RootReach), points,
                           row.second);
        }
        rows.emplace(sphere, std::move(row));
    }
    const auto& [surface, score] = rows.at(sphere);
    const std::vector<occ::handle<Geom_Curve>>& curves = lines.at(sphere);
    for (auto _ : state) {
        for (const occ::handle<Geom_Curve>& line : curves) {
            GeomAPI_IntCS occt(line, surface);
            benchmark::DoNotOptimize(occt);
        }
    }
    state.SetItemsProcessed(state.iterations() * static_cast<std::int64_t>(curves.size()));
    B2::Report(state, score);
}

// --- B2.2 and B2.3: sampled answers ---------------------------------------------

/// The parameters at which an OCCT curve is sampled: its range, a window
/// about the anchor on an unbounded line, [-10, 10] of an unbounded conic.
std::pair<double, double> SampledRange(const occ::handle<Geom_Curve>& curve, const Point3& anchor,
                                       const double scale)
{
    double first = curve->FirstParameter();
    double last = curve->LastParameter();
    if (std::abs(first) < B2::Unbounded && std::abs(last) < B2::Unbounded) {
        return {first, last};
    }
    const Geom_Curve* basis = curve.get();
    if (const auto* trimmed = dynamic_cast<const Geom_TrimmedCurve*>(basis)) {
        basis = trimmed->BasisCurve().get();
    }
    if (const auto* line = dynamic_cast<const Geom_Line*>(basis)) {
        const gp_Lin lin = line->Lin();
        const double centre = gp_Vec(lin.Location(), ToOcct(anchor)).Dot(gp_Vec(lin.Direction()));
        const double half = B2::LineWindow * scale;
        return {std::max(first, centre - half), std::min(last, centre + half)};
    }
    return {std::max(first, -10.0), std::min(last, 10.0)};
}

template <typename Value>
std::vector<Point3> Sampled(const double first, const double last, const Value& value)
{
    std::vector<Point3> samples;
    for (int k = 0; k <= B2::CurveSamples; ++k) {
        samples.push_back(PointOf(value(first + (last - first) * k / B2::CurveSamples)));
    }
    return samples;
}

struct SurfacePair {
    occ::handle<Geom_Surface> A;
    occ::handle<Geom_Surface> B;
    double Tolerance = 0.0;
};

B2::Found FoundOf(const SurfacePair& pair, const B2::SsiCase& c)
{
    B2::Found found;
    try {
        const GeomInt_IntSS occt(pair.A, pair.B, pair.Tolerance, true);
        if (!occt.IsDone()) {
            return found;
        }
        for (int i = 1; i <= occt.NbLines(); ++i) {
            const occ::handle<Geom_Curve>& curve = occt.Line(i);
            const auto [first, last] = SampledRange(curve, c.Anchor, c.Scale);
            found.Curves.push_back(
                Sampled(first, last, [&](const double t) { return curve->Value(t); }));
        }
        for (int i = 1; i <= occt.NbPoints(); ++i) {
            found.Points.push_back(PointOf(occt.Point(i)));
        }
        found.Resolved = true;
    } catch (const std::exception&) {
        found = B2::Found{};
    }
    return found;
}

void OcctSurfaceSurface(benchmark::State& state, const std::string& name)
{
    static std::map<std::string, std::pair<std::vector<SurfacePair>, B2::Score>> rows;
    if (!rows.contains(name)) {
        const std::string kind = name.substr(0, name.find('/'));
        const std::vector<B2::SsiCase> cases =
            B2::FramedCases(kind + "-ssi", name.substr(name.find('/') + 1));
        std::pair<std::vector<SurfacePair>, B2::Score> row;
        std::vector<B2::Found> answers;
        for (const B2::SsiCase& c : cases) {
            row.first.push_back({SurfaceOf(c.A, c.BoxA), SurfaceOf(c.B, c.BoxB), c.Tolerance});
            answers.push_back(FoundOf(row.first.back(), c));
        }
        row.second = B2::ScoreSsi(cases, answers);
        rows.emplace(name, std::move(row));
    }
    const auto& [pairs, score] = rows.at(name);
    for (auto _ : state) {
        for (const SurfacePair& pair : pairs) {
            try {
                GeomInt_IntSS occt(pair.A, pair.B, pair.Tolerance, true);
                benchmark::DoNotOptimize(occt);
            } catch (const std::exception&) {
                benchmark::ClobberMemory();
            }
        }
    }
    state.SetItemsProcessed(state.iterations() * static_cast<std::int64_t>(pairs.size()));
    B2::Report(state, score);
}

// --- B2.3: IntAna_QuadQuadGeo ---------------------------------------------------

/// A pair of planes, spheres and cylinders as OCCT's elementary surfaces.
struct Quadrics {
    std::vector<gp_Pln> Planes;
    std::vector<gp_Cylinder> Cylinders;
    std::vector<gp_Sphere> Spheres;
    double Tolerance = 0.0;
};

Quadrics QuadricsOf(const B2::SsiCase& c)
{
    Quadrics q;
    q.Tolerance = c.Tolerance;
    for (const Stereon::Surface* s : {&c.A, &c.B}) {
        if (const auto* p = s->Target<Stereon::Plane>()) {
            q.Planes.emplace_back(FrameOf(p->Origin, p->XAxis, p->YAxis));
        } else if (const auto* cy = s->Target<Stereon::Cylinder>()) {
            q.Cylinders.emplace_back(FrameOf(cy->Origin, cy->XAxis, cy->YAxis), cy->Radius);
        } else if (const auto* sp = s->Target<Stereon::Sphere>()) {
            q.Spheres.emplace_back(FrameOf(sp->Center, sp->XAxis, sp->YAxis), sp->Radius);
        }
    }
    return q;
}

/// IntAna_QuadQuadGeo on the pair, through the constructor for its types.
IntAna_QuadQuadGeo Solve(const Quadrics& q)
{
    const double angular = Precision::Angular();
    if (q.Planes.size() == 2) {
        return IntAna_QuadQuadGeo(q.Planes[0], q.Planes[1], angular, q.Tolerance);
    }
    if (q.Planes.size() == 1 && q.Cylinders.size() == 1) {
        return IntAna_QuadQuadGeo(q.Planes[0], q.Cylinders[0], angular, q.Tolerance);
    }
    if (q.Planes.size() == 1) {
        return IntAna_QuadQuadGeo(q.Planes[0], q.Spheres[0]);
    }
    if (q.Cylinders.size() == 2) {
        return IntAna_QuadQuadGeo(q.Cylinders[0], q.Cylinders[1], q.Tolerance);
    }
    if (q.Cylinders.size() == 1) {
        return IntAna_QuadQuadGeo(q.Cylinders[0], q.Spheres[0], q.Tolerance);
    }
    return IntAna_QuadQuadGeo(q.Spheres[0], q.Spheres[1], q.Tolerance);
}

/// IntAna_QuadQuadGeo's answer, sampled.
B2::Found FoundOf(const IntAna_QuadQuadGeo& occt, const B2::SsiCase& c)
{
    B2::Found found;
    if (!occt.IsDone() || occt.TypeInter() == IntAna_NoGeometricSolution) {
        return found;
    }
    found.Resolved = true;
    const double period = 2.0 * std::numbers::pi;
    for (int i = 1; i <= occt.NbSolutions(); ++i) {
        switch (occt.TypeInter()) {
        case IntAna_Point:
            found.Points.push_back(PointOf(occt.Point(i)));
            break;
        case IntAna_Line: {
            const gp_Lin lin = occt.Line(i);
            const double centre =
                gp_Vec(lin.Location(), ToOcct(c.Anchor)).Dot(gp_Vec(lin.Direction()));
            const double half = B2::LineWindow * c.Scale;
            found.Curves.push_back(Sampled(centre - half, centre + half,
                                           [&](const double t) { return ElCLib::Value(t, lin); }));
            break;
        }
        case IntAna_Circle: {
            const gp_Circ circle = occt.Circle(i);
            found.Curves.push_back(
                Sampled(0.0, period, [&](const double t) { return ElCLib::Value(t, circle); }));
            break;
        }
        case IntAna_Ellipse: {
            const gp_Elips ellipse = occt.Ellipse(i);
            found.Curves.push_back(
                Sampled(0.0, period, [&](const double t) { return ElCLib::Value(t, ellipse); }));
            break;
        }
        case IntAna_Same:
            found.Overlap = true;
            break;
        default:
            break;
        }
    }
    if (occt.TypeInter() == IntAna_Same) {
        found.Overlap = true;
    }
    return found;
}

void OcctAnalyticSurfaces(benchmark::State& state, const std::string& pair)
{
    static std::map<std::string, std::pair<std::vector<Quadrics>, B2::Score>> rows;
    if (!rows.contains(pair)) {
        const std::vector<B2::SsiCase> cases = B2::PairCases(pair);
        std::pair<std::vector<Quadrics>, B2::Score> row;
        std::vector<B2::Found> answers;
        for (const B2::SsiCase& c : cases) {
            row.first.push_back(QuadricsOf(c));
            try {
                answers.push_back(FoundOf(Solve(row.first.back()), c));
            } catch (const std::exception&) {
                answers.emplace_back();
            }
        }
        row.second = B2::ScoreSsi(cases, answers);
        rows.emplace(pair, std::move(row));
    }
    const auto& [quadrics, score] = rows.at(pair);
    for (auto _ : state) {
        for (const Quadrics& q : quadrics) {
            try {
                IntAna_QuadQuadGeo occt = Solve(q);
                benchmark::DoNotOptimize(occt);
            } catch (const std::exception&) {
                benchmark::ClobberMemory();
            }
        }
    }
    state.SetItemsProcessed(state.iterations() * static_cast<std::int64_t>(quadrics.size()));
    B2::Report(state, score);
}

// --- Registration ---------------------------------------------------------------

void Register(const std::string& name, void (*row)(benchmark::State&, bool), const bool which)
{
    benchmark::RegisterBenchmark(name.c_str(), row, which);
}

void Register(const std::string& name, void (*row)(benchmark::State&, const std::string&),
              const std::string& cases)
{
    benchmark::RegisterBenchmark(name.c_str(), row, cases);
}

[[maybe_unused]] const bool registered = [] {
    Register("B2.1/CurveCurve/conics/OCCT", OcctCurveCurve, false);
    Register("B2.1/CurveCurve/splines/OCCT", OcctCurveCurve, true);
    Register("B2.1/CurveSurface/line-nurbs-sphere/OCCT", OcctCurveSurface, true);
    Register("B2.1/CurveSurface/line-nurbs-cylinder/OCCT", OcctCurveSurface, false);
    for (const std::string& family : B2::NurbsFamilies) {
        Register("B2.2/SSI/nurbs-" + family + "/OCCT", OcctSurfaceSurface, "nurbs/" + family);
    }
    for (const std::string& family : B2::CurvedFamilies) {
        Register("B2.2/SSI/curved-" + family + "/OCCT", OcctSurfaceSurface, "curved/" + family);
    }
    for (const std::string& pair : B2::AnalyticPairs) {
        Register("B2.3/AnalyticSSI/" + pair + "/OCCT", OcctAnalyticSurfaces, pair);
    }
    return true;
}();

} // namespace
