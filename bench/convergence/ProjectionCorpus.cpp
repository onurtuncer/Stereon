// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// B1.3's convergence target (docs/benchmarks.md: point projection onto NURBS
// curves and surfaces, 99.99 % convergence) measured on a declared corpus,
// not inferred from timings. Project(BSplineCurve) and Project(BSplineSurface)
// return a parameter and no refusal, so each case is either a correct global
// closest point or an incorrect success; there is no structured-refusal
// outcome to count.
//
// The corpus, "b13-v1", is generated from a fixed seed: random curves of
// degrees 1 to 7 and surfaces of degrees 1 to 5 in each direction, clamped,
// with repeated interior knots, polynomial and rational, and points of five
// kinds per geometry: on it, just off it, well off it, far away, and near an
// end or a knot. The reference is independent of the projector: the distance
// sampled densely over the domain, and every sampled local minimum refined by
// a derivative-free search (golden sections on a curve, a shrinking compass
// on a surface) that calls only Eval. It can miss a narrow minimum, never
// invent one, so a case counts as converged when the projector's point is no
// farther than the reference's, within a rounding allowance of 1e-9 of the
// scene's size; it is an incorrect success when it is farther.
//
//     stereon_b13_convergence --cases 10000 --out convergence.json
//
// --self-check answers every case with the domain's start instead, and fails
// unless the reference flags most of them: a check that it can see a wrong answer.
//
// writes the evidence tools/check_release_benchmarks.py reads, and lists
// every incorrect success with its seed index so that it can be replayed.
// Exit status 0 when it ran, whatever it measured: the checker decides.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <numbers>
#include <random>
#include <string>
#include <utility>
#include <vector>

import Stereon.Surface;

using Stereon::BSplineCurve;
using Stereon::BSplineSurface;
using Stereon::Point2;
using Stereon::Point3;

namespace {

constexpr const char* Corpus = "b13-v1";
constexpr std::uint32_t Seed = 20261101;

/// The rounding allowance, relative to the scene's size (the poles and the point).
constexpr double Allowance = 1e-9;

/// Samples of the reference search: per curve, per surface direction.
constexpr int CurveSamples = 4000;
constexpr int SurfaceSamples = 80;

/// How many of the closest sampled local minima the reference refines.
constexpr std::size_t RefinedMinima = 8;

/// The golden ratio's conjugate, for the curve's section search.
const double InverseGolden = 1.0 / std::numbers::phi;

class Generator {
public:
    explicit Generator(const std::uint32_t seed)
        : m_Engine(seed)
    {}
    int Int(const int low, const int high) // [low, high]
    {
        return std::uniform_int_distribution<int>(low, high)(m_Engine);
    }
    double Real(const double low, const double high)
    {
        return std::uniform_real_distribution<double>(low, high)(m_Engine);
    }

private:
    std::mt19937 m_Engine;
};

/// Clamped knots for `poles` poles of degree p over [0, 1], interior knots of
/// random spacing repeated up to p times.
std::vector<double> RandomKnots(Generator& g, const int p, const std::size_t poles)
{
    const auto q = static_cast<std::size_t>(p);
    std::vector<double> interior;
    while (interior.size() < poles - q - 1) {
        const double u = g.Real(0.02, 0.98);
        const auto repeat = static_cast<std::size_t>(g.Int(1, p));
        for (std::size_t r = 0; r < repeat && interior.size() < poles - q - 1; ++r) {
            interior.push_back(u);
        }
    }
    std::ranges::sort(interior);
    std::vector<double> knots(q + 1, 0.0);
    knots.insert(knots.end(), interior.begin(), interior.end());
    knots.insert(knots.end(), q + 1, 1.0);
    return knots;
}

Point3 RandomPoint(Generator& g, const double reach)
{
    return {g.Real(-reach, reach), g.Real(-reach, reach), g.Real(-reach, reach)};
}

std::vector<double> RandomWeights(Generator& g, const std::size_t count)
{
    std::vector<double> weights;
    for (std::size_t i = 0; i < count; ++i) {
        weights.push_back(g.Real(0.3, 3.0));
    }
    return weights;
}

BSplineCurve RandomCurve(Generator& g)
{
    for (;;) {
        const int p = g.Int(1, 7);
        const auto poles = static_cast<std::size_t>(p + 1 + g.Int(0, 8));
        std::vector<Point3> points;
        for (std::size_t i = 0; i < poles; ++i) {
            points.push_back(RandomPoint(g, 10.0));
        }
        const std::vector<double> weights =
            g.Int(0, 1) == 1 ? RandomWeights(g, poles) : std::vector<double>{};
        if (auto c = Stereon::MakeBSplineCurve(p, RandomKnots(g, p, poles), points, weights)) {
            return *c;
        }
    }
}

BSplineSurface RandomSurface(Generator& g)
{
    for (;;) {
        const int p = g.Int(1, 5);
        const int q = g.Int(1, 5);
        const auto nu = static_cast<std::size_t>(p + 1 + g.Int(0, 4));
        const auto nv = static_cast<std::size_t>(q + 1 + g.Int(0, 4));
        std::vector<Point3> points;
        for (std::size_t i = 0; i < nu; ++i) {
            for (std::size_t j = 0; j < nv; ++j) {
                // A wavy sheet over a grid, so the surface is not a tangle.
                points.push_back({2.0 * static_cast<double>(i) + g.Real(-0.6, 0.6),
                                  2.0 * static_cast<double>(j) + g.Real(-0.6, 0.6),
                                  g.Real(-3.0, 3.0)});
            }
        }
        const std::vector<double> weights =
            g.Int(0, 1) == 1 ? RandomWeights(g, nu * nv) : std::vector<double>{};
        if (auto s = Stereon::MakeBSplineSurface(p, q, RandomKnots(g, p, nu), RandomKnots(g, q, nv),
                                                 points, weights)) {
            return *s;
        }
    }
}

double Size(const std::vector<Point3>& poles, const Point3& p)
{
    double size = Length(Stereon::ToVec(p));
    for (const Point3& pole : poles) {
        size = std::max(size, Length(Stereon::ToVec(pole)));
    }
    return size;
}

// --- The reference: sampled, then refined with Eval alone -----------------------

double CurveReference(const BSplineCurve& c, const Point3& p)
{
    const Stereon::Interval d = Domain(c);
    const auto at = [&](const double t) { return Distance(Eval(c, t), p); };
    std::vector<double> t;
    for (int i = 0; i <= CurveSamples; ++i) {
        t.push_back(d.GetLower() + (d.GetUpper() - d.GetLower()) * i / CurveSamples);
    }
    t.insert(t.end(), c.Knots.begin(), c.Knots.end()); // the kinks, where minima like to sit
    std::ranges::sort(t);
    t.erase(std::ranges::unique(t).begin(), t.end());
    std::vector<double> f(t.size());
    for (std::size_t i = 0; i < t.size(); ++i) {
        f[i] = at(t[i]);
    }
    std::vector<std::pair<double, std::size_t>> minima;
    for (std::size_t i = 0; i < t.size(); ++i) {
        const bool left = i == 0 || f[i] <= f[i - 1];
        const bool right = i + 1 == t.size() || f[i] <= f[i + 1];
        if (left && right) {
            minima.emplace_back(f[i], i);
        }
    }
    std::ranges::sort(minima);
    double best = std::numeric_limits<double>::infinity();
    for (const double v : f) {
        best = std::min(best, v);
    }
    for (std::size_t m = 0; m < std::min(RefinedMinima, minima.size()); ++m) {
        const std::size_t i = minima[m].second;
        double a = t[i == 0 ? 0 : i - 1];
        double b = t[i + 1 == t.size() ? i : i + 1];
        for (int step = 0; step < 200 && b - a > 1e-16 * (1.0 + std::abs(a)); ++step) {
            const double x1 = b - InverseGolden * (b - a);
            const double x2 = a + InverseGolden * (b - a);
            (at(x1) <= at(x2) ? b : a) = at(x1) <= at(x2) ? x2 : x1;
        }
        best = std::min({best, at(a), at(b)});
    }
    return best;
}

double SurfaceReference(const BSplineSurface& s, const Point3& p, Point2* where = nullptr)
{
    const Stereon::Box2 d = Domain(s);
    const auto clamp = [&](const Point2& uv) {
        return Point2{std::clamp(uv.X, d.X.GetLower(), d.X.GetUpper()),
                      std::clamp(uv.Y, d.Y.GetLower(), d.Y.GetUpper())};
    };
    const auto at = [&](const Point2& uv) { return Distance(Eval(s, clamp(uv)), p); };
    const double du = (d.X.GetUpper() - d.X.GetLower()) / SurfaceSamples;
    const double dv = (d.Y.GetUpper() - d.Y.GetLower()) / SurfaceSamples;
    const auto node = [&](const int i, const int j) {
        return Point2{d.X.GetLower() + du * i, d.Y.GetLower() + dv * j};
    };
    std::vector<double> f((SurfaceSamples + 1) * (SurfaceSamples + 1));
    const auto index = [](const int i, const int j) {
        return static_cast<std::size_t>(i * (SurfaceSamples + 1) + j);
    };
    for (int i = 0; i <= SurfaceSamples; ++i) {
        for (int j = 0; j <= SurfaceSamples; ++j) {
            f[index(i, j)] = at(node(i, j));
        }
    }
    std::vector<std::pair<double, Point2>> minima;
    for (int i = 0; i <= SurfaceSamples; ++i) {
        for (int j = 0; j <= SurfaceSamples; ++j) {
            bool lowest = true;
            for (int a = -1; a <= 1 && lowest; ++a) {
                for (int b = -1; b <= 1 && lowest; ++b) {
                    const int ii = i + a;
                    const int jj = j + b;
                    if ((a != 0 || b != 0) && ii >= 0 && jj >= 0 && ii <= SurfaceSamples &&
                        jj <= SurfaceSamples) {
                        lowest = f[index(i, j)] <= f[index(ii, jj)];
                    }
                }
            }
            if (lowest) {
                minima.emplace_back(f[index(i, j)], node(i, j));
            }
        }
    }
    std::ranges::sort(minima, {}, &std::pair<double, Point2>::first);
    double best = std::numeric_limits<double>::infinity();
    for (const double v : f) {
        best = std::min(best, v);
    }
    for (std::size_t m = 0; m < std::min(RefinedMinima, minima.size()); ++m) {
        // A compass search: try the four neighbours at the current step, move
        // to a better one or halve the step.
        Point2 x = minima[m].second;
        double fx = at(x);
        double hu = du;
        double hv = dv;
        for (int step = 0; step < 400 && (hu > 1e-15 || hv > 1e-15); ++step) {
            bool moved = false;
            for (const Point2 delta :
                 {Point2{hu, 0.0}, Point2{-hu, 0.0}, Point2{0.0, hv}, Point2{0.0, -hv}}) {
                const Point2 y = clamp(Point2{x.X + delta.X, x.Y + delta.Y});
                const double fy = at(y);
                if (fy < fx) {
                    x = y;
                    fx = fy;
                    moved = true;
                }
            }
            if (!moved) {
                hu *= 0.5;
                hv *= 0.5;
            }
        }
        if (fx < best) {
            best = fx;
            if (where != nullptr) {
                *where = x;
            }
        }
    }
    return best;
}

// --- The corpus ------------------------------------------------------------------

/// A point of kind k (0 on, 1 just off, 2 well off, 3 far, 4 near an end or a
/// knot) for a geometry whose sample points `on` gives.
template <typename OnPoint>
Point3 CasePoint(Generator& g, const int kind, const OnPoint& on)
{
    const Point3 base = on(kind == 4);
    switch (kind) {
    case 0:
        return base;
    case 1:
        return base + Stereon::Vec3{g.Real(-1e-3, 1e-3), g.Real(-1e-3, 1e-3), g.Real(-1e-3, 1e-3)};
    case 2:
        return base + Stereon::Vec3{g.Real(-2.0, 2.0), g.Real(-2.0, 2.0), g.Real(-2.0, 2.0)};
    case 3:
        return RandomPoint(g, 40.0);
    default:
        return base + Stereon::Vec3{g.Real(-0.5, 0.5), g.Real(-0.5, 0.5), g.Real(-0.5, 0.5)};
    }
}

struct Tally {
    std::size_t Cases = 0;
    std::size_t Converged = 0;
    std::size_t Incorrect = 0;
    double WorstExcess = 0.0; // the largest (projector - reference) / size
    std::vector<std::string> Failures;
};

void Record(Tally& tally, const char* what, const std::size_t index, const double projected,
            const double reference, const double size)
{
    ++tally.Cases;
    const double excess = (projected - reference) / size;
    tally.WorstExcess = std::max(tally.WorstExcess, excess);
    if (std::isfinite(projected) && excess <= Allowance) {
        ++tally.Converged;
        return;
    }
    ++tally.Incorrect;
    if (tally.Failures.size() < 50) {
        char line[200];
        std::snprintf(line, sizeof line, "%s case %zu: projector %.17g, reference %.17g", what,
                      index, projected, reference);
        tally.Failures.emplace_back(line);
    }
}

/// --self-check: answer every case with the domain's start instead of the
/// projector's parameter, which the reference must catch as incorrect.
bool WrongOnPurpose = false;

Tally Curves(const std::size_t cases)
{
    Generator g(Seed);
    Tally tally;
    for (std::size_t i = 0; i < cases;) {
        const BSplineCurve c = RandomCurve(g);
        const Stereon::CurveProjector projector(c);
        for (int kind = 0; kind < 5 && i < cases; ++kind, ++i) {
            const auto on = [&](const bool special) {
                const Stereon::Interval d = Domain(c);
                const double t = special ? c.Knots[static_cast<std::size_t>(
                                               g.Int(0, static_cast<int>(c.Knots.size()) - 1))]
                                         : g.Real(d.GetLower(), d.GetUpper());
                return Eval(c, t);
            };
            const Point3 p = CasePoint(g, kind, on);
            const double t = WrongOnPurpose ? Domain(c).GetLower() : projector.Project(p);
            const double projected = Distance(Eval(c, t), p);
            Record(tally, "curve", i, projected, CurveReference(c, p), 1.0 + Size(c.Poles, p));
        }
    }
    return tally;
}

/// --dump: print every incorrect surface case in full precision, as C++
/// initialisers, so that it can become a regression test.
bool DumpCases = false;

void DumpList(const char* name, const std::vector<double>& values)
{
    std::printf("    const std::vector<double> %s{", name);
    for (std::size_t i = 0; i < values.size(); ++i) {
        std::printf("%s%.17g", i ? ", " : "", values[i]);
    }
    std::printf("};\n");
}

void DumpSurfaceCase(const std::size_t index, const BSplineSurface& s, const Point3& p,
                     const Point2& uv, const double projected, const double reference)
{
    std::printf("  // surface case %zu: projector (%.17g, %.17g) at %.17g, reference %.17g\n",
                index, uv.X, uv.Y, projected, reference);
    std::printf("    const int degreeU = %d;\n    const int degreeV = %d;\n", s.DegreeU, s.DegreeV);
    DumpList("knotsU", s.KnotsU);
    DumpList("knotsV", s.KnotsV);
    std::printf("    const std::vector<Point3> poles{");
    for (std::size_t i = 0; i < s.Poles.size(); ++i) {
        std::printf("%s{%.17g, %.17g, %.17g}", i ? ", " : "", s.Poles[i].X, s.Poles[i].Y,
                    s.Poles[i].Z);
    }
    std::printf("};\n");
    DumpList("weights", s.Weights);
    std::printf("    const Point3 point{%.17g, %.17g, %.17g};\n", p.X, p.Y, p.Z);
}

Tally Surfaces(const std::size_t cases)
{
    Generator g(Seed + 1);
    Tally tally;
    for (std::size_t i = 0; i < cases;) {
        const BSplineSurface s = RandomSurface(g);
        const Stereon::SurfaceProjector projector(s);
        for (int kind = 0; kind < 5 && i < cases; ++kind, ++i) {
            const auto on = [&](const bool special) {
                const Stereon::Box2 d = Domain(s);
                const double u = special ? (g.Int(0, 1) == 1 ? d.X.GetUpper() : d.X.GetLower())
                                         : g.Real(d.X.GetLower(), d.X.GetUpper());
                return Eval(s, Point2{u, g.Real(d.Y.GetLower(), d.Y.GetUpper())});
            };
            const Point3 p = CasePoint(g, kind, on);
            const Point2 uv = WrongOnPurpose ? Stereon::MinCorner(Domain(s)) : projector.Project(p);
            const double projected = Distance(Eval(s, uv), p);
            Point2 found{};
            const double reference = SurfaceReference(s, p, &found);
            const std::size_t before = tally.Incorrect;
            Record(tally, "surface", i, projected, reference, 1.0 + Size(s.Poles, p));
            if (DumpCases && tally.Incorrect > before) {
                DumpSurfaceCase(i, s, p, uv, projected, reference);
                std::printf("    // reference at (%.17g, %.17g)\n", found.X, found.Y);
            }
        }
    }
    return tally;
}

/// The one-sided 95 % Clopper-Pearson lower bound of a success rate.
double LowerBound(const std::size_t successes, const std::size_t trials)
{
    if (successes == trials) {
        return std::pow(0.05, 1.0 / static_cast<double>(trials));
    }
    // Bisection on the binomial tail, enough for a report.
    const auto tail = [&](const double rate) {
        double sum = 0.0;
        double term = std::pow(rate, static_cast<double>(trials)); // k = trials
        for (std::size_t k = trials; k >= successes; --k) {
            sum += term;
            if (k == 0) {
                break;
            }
            term *=
                static_cast<double>(k) / static_cast<double>(trials - k + 1) * (1.0 - rate) / rate;
        }
        return sum;
    };
    double low = 0.0;
    double high = 1.0;
    for (int i = 0; i < 100; ++i) {
        const double mid = 0.5 * (low + high);
        (tail(mid) < 0.05 ? low : high) = mid;
    }
    return low;
}

void Write(std::ostream& out, const char* name, const Tally& t)
{
    out << "    \"" << name << "\": {\"cases\": " << t.Cases << ", \"converged\": " << t.Converged
        << ", \"incorrect_success\": " << t.Incorrect << ", \"worst_excess\": " << t.WorstExcess
        << ", \"failures\": [";
    for (std::size_t i = 0; i < t.Failures.size(); ++i) {
        out << (i ? ", " : "") << '"' << t.Failures[i] << '"';
    }
    out << "]}";
}

} // namespace

int main(int argc, char** argv)
{
    std::size_t cases = 10000;
    std::string path = "convergence.json";
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--self-check") == 0) {
            WrongOnPurpose = true;
        } else if (std::strcmp(argv[i], "--dump") == 0) {
            DumpCases = true;
        }
    }
    for (int i = 1; i + 1 < argc; ++i) {
        if (std::strcmp(argv[i], "--cases") == 0) {
            cases = static_cast<std::size_t>(std::strtoull(argv[++i], nullptr, 10));
        } else if (std::strcmp(argv[i], "--out") == 0) {
            path = argv[++i];
        }
    }
    const Tally curves = Curves(cases);
    const Tally surfaces = Surfaces(cases);
    const std::size_t total = curves.Cases + surfaces.Cases;
    const std::size_t converged = curves.Converged + surfaces.Converged;
    std::ofstream out(path);
    out.precision(17);
    out << "{\n  \"corpus\": \"" << Corpus << "\",\n  \"seed\": " << Seed
        << ",\n  \"allowance\": " << Allowance << ",\n  \"cases\": " << total
        << ",\n  \"converged\": " << converged
        << ",\n  \"incorrect_success\": " << curves.Incorrect + surfaces.Incorrect
        << ",\n  \"rate\": " << static_cast<double>(converged) / static_cast<double>(total)
        << ",\n  \"lower_bound_95\": " << LowerBound(converged, total) << ",\n  \"by_kind\": {\n";
    Write(out, "curves", curves);
    out << ",\n";
    Write(out, "surfaces", surfaces);
    out << "\n  }\n}\n";
    std::printf(
        "%s: %zu of %zu converged (curves %zu/%zu, surfaces %zu/%zu), %zu incorrect successes\n",
        Corpus, converged, total, curves.Converged, curves.Cases, surfaces.Converged,
        surfaces.Cases, curves.Incorrect + surfaces.Incorrect);
    for (const Tally* t : {&curves, &surfaces}) {
        for (const std::string& failure : t->Failures) {
            std::printf("  %s\n", failure.c_str());
        }
    }
    if (WrongOnPurpose) {
        // The reference must see through a projector that is wrong on purpose.
        const bool caught =
            curves.Incorrect > curves.Cases / 2 && surfaces.Incorrect > surfaces.Cases / 2;
        std::printf("self-check: %s\n",
                    caught ? "the reference caught the wrong answers" : "MISSED");
        return caught ? 0 : 1;
    }
    return out ? 0 : 1;
}
