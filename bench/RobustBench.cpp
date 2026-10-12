// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Benchmark B0.1 of docs/benchmarks.md: the adaptive predicates against the
// plain floating-point determinant (OCCT has no exact predicate), on random
// points in the unit square or cube, where stage A decides everything, and
// on nearly degenerate points built in floating point, where the adaptive
// stages run. The exact stages are measured on their own as the worst case.
// ADR-0002 asks for the filtered predicates within 1.5x of the naive cost on
// non-degenerate input.

#include <benchmark/benchmark.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <random>
#include <vector>

import Stereon.Robust;

namespace {

using Stereon::Point2;
using Stereon::Point3;
using Stereon::Sign;
using Stereon::SignOf;
using Stereon::ToInt;
using Stereon::Vec2;
using Stereon::Vec3;

constexpr std::size_t CaseCount = 1 << 14;

struct Case2 {
    Point2 A;
    Point2 B;
    Point2 C;
    Point2 D;
};

struct Case3 {
    Point3 A;
    Point3 B;
    Point3 C;
    Point3 D;
    Point3 E;
};

void SetItems(benchmark::State& state)
{
    state.SetItemsProcessed(state.iterations() * static_cast<std::int64_t>(CaseCount));
}

// --- Plain floating point ---------------------------------------------------

Sign NaiveOrient2D(const Point2& a, const Point2& b, const Point2& c)
{
    return SignOf((a.X - c.X) * (b.Y - c.Y) - (a.Y - c.Y) * (b.X - c.X));
}

double NaiveDet3(const Vec3& p, const Vec3& q, const Vec3& r)
{
    return p.X * (q.Y * r.Z - q.Z * r.Y) - p.Y * (q.X * r.Z - q.Z * r.X) +
           p.Z * (q.X * r.Y - q.Y * r.X);
}

Sign NaiveOrient3D(const Point3& a, const Point3& b, const Point3& c, const Point3& d)
{
    return SignOf(NaiveDet3(a - d, b - d, c - d));
}

Sign NaiveInCircle(const Point2& a, const Point2& b, const Point2& c, const Point2& d)
{
    const Vec2 ad = a - d;
    const Vec2 bd = b - d;
    const Vec2 cd = c - d;
    return SignOf((ad.X * ad.X + ad.Y * ad.Y) * (bd.X * cd.Y - bd.Y * cd.X) +
                  (bd.X * bd.X + bd.Y * bd.Y) * (cd.X * ad.Y - cd.Y * ad.X) +
                  (cd.X * cd.X + cd.Y * cd.Y) * (ad.X * bd.Y - ad.Y * bd.X));
}

Sign NaiveInSphere(const Point3& a, const Point3& b, const Point3& c, const Point3& d,
                   const Point3& e)
{
    const Vec3 ae = a - e;
    const Vec3 be = b - e;
    const Vec3 ce = c - e;
    const Vec3 de = d - e;
    const auto lift = [](const Vec3& v) { return v.X * v.X + v.Y * v.Y + v.Z * v.Z; };
    return SignOf(lift(de) * NaiveDet3(ae, be, ce) - lift(ce) * NaiveDet3(ae, be, de) +
                  lift(be) * NaiveDet3(ae, ce, de) - lift(ae) * NaiveDet3(be, ce, de));
}

// --- Inputs -----------------------------------------------------------------

class Inputs {
public:
    explicit Inputs(const std::uint64_t seed)
        : m_Engine(seed)
    {}

    Point2 Point2InSquare() { return {m_Unit(m_Engine), m_Unit(m_Engine)}; }
    Point3 Point3InCube() { return {m_Unit(m_Engine), m_Unit(m_Engine), m_Unit(m_Engine)}; }

    /// Random points: stage A decides.
    std::vector<Case2> Random2()
    {
        std::vector<Case2> cases(CaseCount);
        for (Case2& c : cases) {
            c = {Point2InSquare(), Point2InSquare(), Point2InSquare(), Point2InSquare()};
        }
        return cases;
    }
    std::vector<Case3> Random3()
    {
        std::vector<Case3> cases(CaseCount);
        for (Case3& c : cases) {
            c = {Point3InCube(), Point3InCube(), Point3InCube(), Point3InCube(), Point3InCube()};
        }
        return cases;
    }

    /// C on the line AB and D on the circle through the origin-centred unit
    /// radius, both to within roundoff: the adaptive stages run.
    std::vector<Case2> NearlyDegenerate2()
    {
        std::vector<Case2> cases(CaseCount);
        for (Case2& c : cases) {
            c.A = Point2InSquare();
            c.B = Point2InSquare();
            const double t = m_Unit(m_Engine);
            c.C = {c.A.X + t * (c.B.X - c.A.X), c.A.Y + t * (c.B.Y - c.A.Y)};
            c.D = OnCircle();
        }
        return cases;
    }
    std::vector<Case3> NearlyDegenerate3()
    {
        std::vector<Case3> cases(CaseCount);
        for (Case3& c : cases) {
            c.A = Point3InCube();
            c.B = Point3InCube();
            c.C = Point3InCube();
            const double s = m_Unit(m_Engine);
            const double t = m_Unit(m_Engine);
            c.D = {c.A.X + s * (c.B.X - c.A.X) + t * (c.C.X - c.A.X),
                   c.A.Y + s * (c.B.Y - c.A.Y) + t * (c.C.Y - c.A.Y),
                   c.A.Z + s * (c.B.Z - c.A.Z) + t * (c.C.Z - c.A.Z)};
            c.E = OnSphere();
        }
        return cases;
    }

    /// Four points on the unit circle, for the InCircle cases.
    std::vector<Case2> Cocircular()
    {
        std::vector<Case2> cases(CaseCount);
        for (Case2& c : cases) {
            c = {OnCircle(), OnCircle(), OnCircle(), OnCircle()};
        }
        return cases;
    }
    std::vector<Case3> Cospherical()
    {
        std::vector<Case3> cases(CaseCount);
        for (Case3& c : cases) {
            c = {OnSphere(), OnSphere(), OnSphere(), OnSphere(), OnSphere()};
        }
        return cases;
    }

private:
    Point2 OnCircle()
    {
        const double angle = m_Unit(m_Engine) * std::numbers::pi;
        return {std::cos(angle), std::sin(angle)};
    }
    Point3 OnSphere()
    {
        const Vec3 v{m_Unit(m_Engine), m_Unit(m_Engine), m_Unit(m_Engine)};
        const double length = std::sqrt(v.X * v.X + v.Y * v.Y + v.Z * v.Z);
        return {v.X / length, v.Y / length, v.Z / length};
    }

    std::mt19937_64 m_Engine;
    std::uniform_real_distribution<double> m_Unit{-1.0, 1.0};
};

template <typename Predicate>
void Run2(benchmark::State& state, const std::vector<Case2>& cases, Predicate predicate)
{
    for (auto _ : state) {
        int total = 0;
        for (const Case2& c : cases) {
            total += ToInt(predicate(c));
        }
        benchmark::DoNotOptimize(total);
    }
    SetItems(state);
}

template <typename Predicate>
void Run3(benchmark::State& state, const std::vector<Case3>& cases, Predicate predicate)
{
    for (auto _ : state) {
        int total = 0;
        for (const Case3& c : cases) {
            total += ToInt(predicate(c));
        }
        benchmark::DoNotOptimize(total);
    }
    SetItems(state);
}

// --- Orient2D ---------------------------------------------------------------

void Orient2DNaive(benchmark::State& state)
{
    Run2(state, Inputs{1}.Random2(), [](const Case2& c) { return NaiveOrient2D(c.A, c.B, c.C); });
}
void Orient2DRandom(benchmark::State& state)
{
    Run2(state, Inputs{1}.Random2(),
         [](const Case2& c) { return Stereon::Orient2D(c.A, c.B, c.C); });
}
void Orient2DNearlyCollinear(benchmark::State& state)
{
    Run2(state, Inputs{2}.NearlyDegenerate2(),
         [](const Case2& c) { return Stereon::Orient2D(c.A, c.B, c.C); });
}
void Orient2DExactStage(benchmark::State& state)
{
    Run2(state, Inputs{1}.Random2(),
         [](const Case2& c) { return Stereon::Orient2DExact(c.A, c.B, c.C); });
}

// --- Orient3D ---------------------------------------------------------------

void Orient3DNaive(benchmark::State& state)
{
    Run3(state, Inputs{3}.Random3(),
         [](const Case3& c) { return NaiveOrient3D(c.A, c.B, c.C, c.D); });
}
void Orient3DRandom(benchmark::State& state)
{
    Run3(state, Inputs{3}.Random3(),
         [](const Case3& c) { return Stereon::Orient3D(c.A, c.B, c.C, c.D); });
}
void Orient3DNearlyCoplanar(benchmark::State& state)
{
    Run3(state, Inputs{4}.NearlyDegenerate3(),
         [](const Case3& c) { return Stereon::Orient3D(c.A, c.B, c.C, c.D); });
}
void Orient3DExactStage(benchmark::State& state)
{
    Run3(state, Inputs{3}.Random3(),
         [](const Case3& c) { return Stereon::Orient3DExact(c.A, c.B, c.C, c.D); });
}

// --- InCircle ---------------------------------------------------------------

void InCircleNaive(benchmark::State& state)
{
    Run2(state, Inputs{5}.Random2(),
         [](const Case2& c) { return NaiveInCircle(c.A, c.B, c.C, c.D); });
}
void InCircleRandom(benchmark::State& state)
{
    Run2(state, Inputs{5}.Random2(),
         [](const Case2& c) { return Stereon::InCircle(c.A, c.B, c.C, c.D); });
}
void InCircleNearlyCocircular(benchmark::State& state)
{
    Run2(state, Inputs{6}.Cocircular(),
         [](const Case2& c) { return Stereon::InCircle(c.A, c.B, c.C, c.D); });
}
void InCircleExactStage(benchmark::State& state)
{
    Run2(state, Inputs{5}.Random2(),
         [](const Case2& c) { return Stereon::InCircleExact(c.A, c.B, c.C, c.D); });
}

// --- InSphere ---------------------------------------------------------------

void InSphereNaive(benchmark::State& state)
{
    Run3(state, Inputs{7}.Random3(),
         [](const Case3& c) { return NaiveInSphere(c.A, c.B, c.C, c.D, c.E); });
}
void InSphereRandom(benchmark::State& state)
{
    Run3(state, Inputs{7}.Random3(),
         [](const Case3& c) { return Stereon::InSphere(c.A, c.B, c.C, c.D, c.E); });
}
void InSphereNearlyCospherical(benchmark::State& state)
{
    Run3(state, Inputs{8}.Cospherical(),
         [](const Case3& c) { return Stereon::InSphere(c.A, c.B, c.C, c.D, c.E); });
}
void InSphereExactStage(benchmark::State& state)
{
    Run3(state, Inputs{7}.Random3(),
         [](const Case3& c) { return Stereon::InSphereExact(c.A, c.B, c.C, c.D, c.E); });
}

} // namespace

BENCHMARK(Orient2DNaive);
BENCHMARK(Orient2DRandom);
BENCHMARK(Orient2DNearlyCollinear);
BENCHMARK(Orient2DExactStage);
BENCHMARK(Orient3DNaive);
BENCHMARK(Orient3DRandom);
BENCHMARK(Orient3DNearlyCoplanar);
BENCHMARK(Orient3DExactStage);
BENCHMARK(InCircleNaive);
BENCHMARK(InCircleRandom);
BENCHMARK(InCircleNearlyCocircular);
BENCHMARK(InCircleExactStage);
BENCHMARK(InSphereNaive);
BENCHMARK(InSphereRandom);
BENCHMARK(InSphereNearlyCospherical);
BENCHMARK(InSphereExactStage);
