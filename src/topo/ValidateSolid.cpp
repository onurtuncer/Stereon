// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include "Profiling.hpp"

// The contract anchor (docs/toolchain/gcc16-feature-probe.md).
#include <Stereon/Core/Contracts.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <format>
#include <limits>
#include <map>
#include <numbers>
#include <optional>
#include <string>
#include <utility>
#include <vector>

module Stereon.Topo;

/// The Solid level of the validity checker (ADR-0018), on structurally and
/// geometrically valid tables, for each shell of each solid:
///
/// - L01: the shell is closed, every edge its faces use used exactly twice
///   by them. A Degenerate edge, a pole with no extent, is exempt.
/// - L02: the two uses run the edge in opposite senses, so the faces on
///   either side are oriented consistently.
/// - L03: the signed volume of a closed, consistent shell is decided, and
///   the solid has exactly one outer shell, of positive volume (normals
///   outward), the others voids of negative volume (normals inward).
/// - L04: every void lies inside the outer shell: the outer shell's winding
///   number at a vertex of the void is 1.
/// - L05: face–face self-intersection is not checked: it needs ADR-0010's
///   generic surface–surface path. Each solid is reported, never passed.
///
/// Volume and winding number are surface integrals over each face, turned
/// by Green's theorem into integrals over its loops in the face's parameter
/// space: with f(u, v) the integrand times |S_u x S_v|, and F(u, v) its
/// integral in u from a fixed u0, the face's integral is the sum over its
/// coedges of the integral of F(c(t)) v'(t) dt. A loop runs counter-
/// clockwise about its face's own normal, reversed or not, so in the
/// parameters it runs counter-clockwise exactly when the normal is
/// S_u x S_v: the orientation is in the loops, and no sign for a reversed
/// face is applied on top of it. Both integrals are 16-point Gauss–Legendre
/// in t and in u, once over each range and once over its halves; their
/// difference, with the rounding of the sum, is the error estimate a sign
/// is decided against. An estimate, not a certified bound: a decision it
/// cannot make is a finding, never a pass.

namespace Stereon::Detail {

// Instantiated at namespace scope before any function uses it
// (docs/toolchain/gcc16-feature-probe.md).
static_assert(sizeof(Vec3) == 3 * sizeof(double));

namespace {

/// The relative rounding allowed on a sum of face integrals, about 1e-12 of
/// their magnitudes: several digits below the precision of the evaluated
/// geometry, as Stereon.Curve's arc length aims (ADR-0028).
constexpr double RelativeRounding = 0x1p-40;

/// A sign is decided when the value exceeds this many error estimates.
constexpr double DecisionMargin = 4.0;

/// A winding number is read as the whole number within half of it, error
/// included: 1 inside a closed shell, 0 outside.
constexpr double HalfWinding = 0.5;

/// An integral with an estimate of its error.
struct Estimate {
    double Value = 0.0;
    double Error = 0.0;
    double Magnitude = 0.0;

    void Add(const Estimate& other)
    {
        Value += other.Value;
        Error += other.Error;
        Magnitude += other.Magnitude;
    }

    /// The error with the rounding of the sum.
    [[nodiscard]] double Bound() const { return Error + RelativeRounding * Magnitude; }
};

/// The quantity integrated over a shell: at a point S of a surface with
/// partials Su and Sv, the integrand times |Su x Sv|.
using Integrand = double (*)(const Point3& s, const Vec3& su, const Vec3& sv, const Point3& p);

/// x . n / 3: its integral over a closed surface is the enclosed volume.
double VolumeIntegrand(const Point3& s, const Vec3& su, const Vec3& sv, const Point3& /*p*/)
{
    return Dot(ToVec(s), Cross(su, sv)) / 3.0;
}

/// (x - p) . n / (4 pi |x - p|^3): its integral over a closed surface is the
/// winding number at p.
double WindingIntegrand(const Point3& s, const Vec3& su, const Vec3& sv, const Point3& p)
{
    const Vec3 d = s - p;
    const double r = Length(d);
    return Dot(d, Cross(su, sv)) / (4.0 * std::numbers::pi * r * r * r);
}

/// 16-point Gauss–Legendre of g over [a, b] in `pieces` equal pieces.
template <typename G>
[[nodiscard]] double Quadrature(const G& g, const double a, const double b, const int pieces)
{
    double sum = 0.0;
    const double width = (b - a) / pieces;
    for (int k = 0; k < pieces; ++k) {
        const double lo = a + width * k;
        const double half = 0.5 * width;
        const double middle = lo + half;
        for (const auto& [node, weight] : GaussLegendre16()) {
            sum += weight * half * g(middle + half * node);
        }
    }
    return sum;
}

class SolidValidator {
public:
    explicit SolidValidator(const TopologyTables& tables, ValidityReport& report)
        : m_T(tables)
        , m_Report(report)
    {}

    void Run()
    {
        for (std::uint32_t s = 0; s < m_T.Solids.GetSize(); ++s) {
            if (m_T.IsLive<EntityKind::Solid>(s)) {
                CheckSolid(s);
            }
        }
    }

    /// The shell's integral of `integrand` with its estimate, or nothing
    /// when a coedge has no pcurve to integrate along.
    [[nodiscard]] std::optional<Estimate>
    ShellIntegral(const std::uint32_t shell, const Integrand integrand, const Point3& p) const
    {
        Estimate total;
        for (const std::uint32_t f : FacesOf(shell)) {
            const std::optional<Estimate> face = FaceIntegral(f, integrand, p);
            if (!face) {
                return std::nullopt;
            }
            total.Add(*face);
        }
        return total;
    }

private:
    void Report(const TopologyRule rule, std::vector<std::string> entities, std::string message)
    {
        m_Report.Add(
            {.Rule = rule, .Entities = std::move(entities), .Message = std::move(message)});
    }

    [[nodiscard]] std::string Name(const EntityKind kind, const std::uint32_t index) const
    {
        return DescribeEntity(m_T, kind, index);
    }

    // --- Walks, bounded by the tables they walk -------------------------

    [[nodiscard]] std::vector<std::uint32_t> ShellsOf(const std::uint32_t solid) const
    {
        std::vector<std::uint32_t> shells;
        std::uint32_t s = m_T.Solids[solid].FirstShell.GetIndex();
        while (s != ShellId::Invalid && shells.size() < m_T.Shells.GetSize()) {
            shells.push_back(s);
            s = m_T.Shells[s].NextInSolid.GetIndex();
        }
        return shells;
    }

    [[nodiscard]] std::vector<std::uint32_t> FacesOf(const std::uint32_t shell) const
    {
        std::vector<std::uint32_t> faces;
        std::uint32_t f = m_T.Shells[shell].FirstFace.GetIndex();
        while (f != FaceId::Invalid && faces.size() < m_T.Faces.GetSize()) {
            faces.push_back(f);
            f = m_T.Faces[f].NextInShell.GetIndex();
        }
        return faces;
    }

    [[nodiscard]] std::vector<std::uint32_t> LoopsOf(const std::uint32_t face) const
    {
        std::vector<std::uint32_t> loops;
        std::uint32_t l = m_T.Faces[face].FirstLoop.GetIndex();
        while (l != LoopId::Invalid && loops.size() < m_T.Loops.GetSize()) {
            loops.push_back(l);
            l = m_T.Loops[l].NextInFace.GetIndex();
        }
        return loops;
    }

    [[nodiscard]] std::vector<std::uint32_t> CoedgesOf(const std::uint32_t loop) const
    {
        std::vector<std::uint32_t> coedges;
        const std::uint32_t anchor = m_T.Loops[loop].Anchor.GetIndex();
        std::uint32_t c = anchor;
        while (c != CoedgeId::Invalid && coedges.size() < m_T.Coedges.GetSize()) {
            coedges.push_back(c);
            c = m_T.Coedges[c].Next.GetIndex();
            if (c == anchor) {
                break;
            }
        }
        return coedges;
    }

    [[nodiscard]] std::vector<std::uint32_t> CoedgesOfFace(const std::uint32_t face) const
    {
        std::vector<std::uint32_t> coedges;
        for (const std::uint32_t l : LoopsOf(face)) {
            const std::vector<std::uint32_t> of = CoedgesOf(l);
            coedges.insert(coedges.end(), of.begin(), of.end());
        }
        return coedges;
    }

    // --- L01, L02: closure and orientation --------------------------------

    /// True when the shell is closed and consistently oriented; reports
    /// what is not.
    bool CheckClosure(const std::uint32_t shell)
    {
        std::map<std::uint32_t, std::vector<std::uint32_t>> uses;
        for (const std::uint32_t f : FacesOf(shell)) {
            for (const std::uint32_t c : CoedgesOfFace(f)) {
                uses[m_T.Coedges[c].Edge.GetIndex()].push_back(c);
            }
        }
        bool sound = true;
        for (const auto& [e, of] : uses) {
            if (m_T.Edges[e].Kind == EdgeCurveKind::Degenerate) {
                continue;
            }
            if (of.size() != 2) {
                Report(TopologyRule::ShellOpen,
                       {Name(EntityKind::Shell, shell), Name(EntityKind::Edge, e)},
                       std::format("the shell uses the edge {} times, not twice: it is not "
                                   "closed",
                                   of.size()));
                sound = false;
            } else if (m_T.Coedges[of[0]].Sense == m_T.Coedges[of[1]].Sense) {
                Report(TopologyRule::ShellOrientation,
                       {Name(EntityKind::Shell, shell), Name(EntityKind::Edge, e),
                        Name(EntityKind::Coedge, of[0]), Name(EntityKind::Coedge, of[1])},
                       "both uses run the edge the same way: the faces on either side are "
                       "oriented against each other");
                sound = false;
            }
        }
        return sound;
    }

    // --- Integrals over faces ----------------------------------------------

    /// F(u, v), the integral of f in u from u0, in `pieces` pieces.
    [[nodiscard]] static double Inner(const Surface& surface, const Integrand integrand,
                                      const Point3& p, const double u0, const Point2& uv,
                                      const int pieces)
    {
        const auto f = [&](const double u) {
            const Point2 at{u, uv.Y};
            return integrand(Eval(surface, at), Deriv(surface, at, 1, 0), Deriv(surface, at, 0, 1),
                             p);
        };
        return Quadrature(f, u0, uv.X, pieces);
    }

    /// The coedge's term, F(c(t)) v'(t) along its loop's direction, in
    /// `pieces` pieces in t and in u.
    [[nodiscard]] double CoedgeTerm(const Surface& surface, const std::uint32_t c,
                                    const Integrand integrand, const Point3& p, const double u0,
                                    const int pieces) const
    {
        const CoedgeRecord& use = m_T.Coedges[c];
        const Curve2d pcurve = m_T.PCurves[use.PCurve.GetIndex()];
        const Interval range = m_T.Edges[use.Edge.GetIndex()].Range;
        const auto g = [&](const double t) {
            return Inner(surface, integrand, p, u0, Eval(pcurve, t), pieces) *
                   Deriv(pcurve, t, 1).Y;
        };
        const double along = Quadrature(g, range.GetLower(), range.GetUpper(), pieces);
        return use.Sense == EdgeUseSense::Same ? along : -along;
    }

    /// The start of every inner integral of the face, so that it stays in
    /// the surface's domain: the domain's lower u where it is bounded, the
    /// smallest u the face's pcurves reach at their ends otherwise.
    [[nodiscard]] double StartU(const Surface& surface,
                                const std::vector<std::uint32_t>& coedges) const
    {
        const double lower = Domain(surface).X.GetLower();
        if (std::isfinite(lower)) {
            return lower;
        }
        double u0 = std::numeric_limits<double>::infinity();
        for (const std::uint32_t c : coedges) {
            const CoedgeRecord& use = m_T.Coedges[c];
            const Curve2d pcurve = m_T.PCurves[use.PCurve.GetIndex()];
            const Interval range = m_T.Edges[use.Edge.GetIndex()].Range;
            u0 = std::min({u0, Eval(pcurve, range.GetLower()).X, Eval(pcurve, range.GetUpper()).X});
        }
        return u0;
    }

    [[nodiscard]] std::optional<Estimate>
    FaceIntegral(const std::uint32_t face, const Integrand integrand, const Point3& p) const
    {
        const Surface surface = m_T.Surfaces[m_T.Faces[face].Surface.GetIndex()];
        const std::vector<std::uint32_t> coedges = CoedgesOfFace(face);
        for (const std::uint32_t c : coedges) {
            if (!m_T.Coedges[c].PCurve.IsValid()) {
                return std::nullopt;
            }
        }
        const double u0 = StartU(surface, coedges);
        Estimate total;
        for (const std::uint32_t c : coedges) {
            const double coarse = CoedgeTerm(surface, c, integrand, p, u0, 1);
            const double fine = CoedgeTerm(surface, c, integrand, p, u0, 2);
            total.Value += fine;
            total.Error += std::abs(fine - coarse);
            total.Magnitude += std::abs(fine);
        }
        return total;
    }

    // --- L03, L04: volumes and nesting -------------------------------------

    /// A closed, consistent shell's volume, or nothing after reporting why
    /// its sign cannot be decided.
    [[nodiscard]] std::optional<double> DecidedVolume(const std::uint32_t shell)
    {
        const std::optional<Estimate> volume = ShellIntegral(shell, &VolumeIntegrand, Point3{});
        if (!volume) {
            Report(TopologyRule::ShellVolume, {Name(EntityKind::Shell, shell)},
                   "a coedge of the shell has no pcurve, so its volume cannot be integrated");
            return std::nullopt;
        }
        if (!std::isfinite(volume->Value) ||
            !(std::abs(volume->Value) > DecisionMargin * volume->Bound())) {
            Report(TopologyRule::ShellVolume, {Name(EntityKind::Shell, shell)},
                   std::format("the shell's volume {} ± {} is not finite or too small to decide "
                               "its orientation",
                               volume->Value, volume->Bound()));
            return std::nullopt;
        }
        return volume->Value;
    }

    /// A point of the shell: the start of its first face's first edge.
    [[nodiscard]] std::optional<Point3> PointOf(const std::uint32_t shell) const
    {
        for (const std::uint32_t f : FacesOf(shell)) {
            for (const std::uint32_t c : CoedgesOfFace(f)) {
                const VertexId v = m_T.Edges[m_T.Coedges[c].Edge.GetIndex()].Start;
                return m_T.Vertices[v.GetIndex()].Position;
            }
        }
        return std::nullopt;
    }

    void CheckNesting(const std::uint32_t solid, const std::uint32_t outer,
                      const std::uint32_t inner)
    {
        const std::optional<Point3> p = PointOf(inner);
        const std::optional<Estimate> w =
            p ? ShellIntegral(outer, &WindingIntegrand, *p) : std::nullopt;
        if (w && std::abs(w->Value - 1.0) + w->Bound() < HalfWinding) {
            return; // inside
        }
        const bool outside = w && std::abs(w->Value) + w->Bound() < HalfWinding;
        Report(TopologyRule::ShellNesting,
               {Name(EntityKind::Solid, solid), Name(EntityKind::Shell, inner),
                Name(EntityKind::Shell, outer)},
               outside ? "the void shell lies outside the outer shell"
                       : std::format("whether the void shell lies inside the outer shell could "
                                     "not be decided: winding number {} ± {}",
                                     w ? w->Value : std::nan(""), w ? w->Bound() : std::nan("")));
    }

    void CheckShells(const std::uint32_t solid, const std::vector<std::uint32_t>& shells,
                     const std::vector<double>& volumes)
    {
        std::vector<std::uint32_t> outer;
        for (std::size_t i = 0; i < shells.size(); ++i) {
            // Decided against its error estimate by DecidedVolume.
            if (volumes[i] > 0.0) { // STEREON-LINT-ALLOW(epsilon): a decided sign
                outer.push_back(shells[i]);
            }
        }
        if (outer.size() != 1) {
            Report(TopologyRule::ShellVolume, {Name(EntityKind::Solid, solid)},
                   std::format("the solid has {} shells of positive volume, not one outer "
                               "shell: {}",
                               outer.size(),
                               outer.empty() ? "its normals point inward" : "two bound it"));
            return;
        }
        for (const std::uint32_t s : shells) {
            if (s != outer.front()) {
                CheckNesting(solid, outer.front(), s);
            }
        }
    }

    void CheckSolid(const std::uint32_t solid)
    {
        const std::vector<std::uint32_t> shells = ShellsOf(solid);
        bool sound = true;
        for (const std::uint32_t s : shells) {
            sound = CheckClosure(s) && sound;
        }
        std::vector<double> volumes;
        for (const std::uint32_t s : sound ? shells : std::vector<std::uint32_t>{}) {
            const std::optional<double> v = DecidedVolume(s);
            sound = v.has_value() && sound;
            volumes.push_back(v.value_or(0.0));
        }
        if (sound && !shells.empty()) {
            CheckShells(solid, shells, volumes);
        }
        Report(TopologyRule::SelfIntersectionUnchecked, {Name(EntityKind::Solid, solid)},
               "face–face self-intersection is not checked: it needs ADR-0010's generic "
               "surface–surface path, so the solid is not reported valid at the Solid level");
    }

    const TopologyTables& m_T;
    ValidityReport& m_Report;
};

} // namespace

void ValidateSolid(const TopologyTables& tables, ValidityReport& report)
{
    STEREON_PROFILE_SCOPE("topo.ValidateSolid");
    SolidValidator(tables, report).Run();
}

std::optional<double> ShellVolume(const TopologyTables& tables, const std::uint32_t shell)
{
    ValidityReport unused;
    const std::optional<Estimate> volume =
        SolidValidator(tables, unused).ShellIntegral(shell, &VolumeIntegrand, Point3{});
    return volume ? std::optional<double>(volume->Value) : std::nullopt;
}

} // namespace Stereon::Detail
