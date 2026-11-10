// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

// The contract anchor (docs/toolchain/gcc16-feature-probe.md).
#include <Stereon/Core/Contracts.hpp>
#include <Stereon/Core/Try.hpp>

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <memory>
#include <numbers>
#include <optional>
#include <ranges>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

module Stereon.Intersect;

import :Dispatch;
import Stereon.Robust;

#include "Analytic.hpp"
#include "Bernstein.hpp"
#include "BezierArc.hpp"
#include "BezierPatch.hpp"
#include "ConicArc.hpp"
#include "CurvePieces.hpp"

/// Any built-in curve against a B-spline surface, polynomial or rational.
/// The curve is cut into enclosed Bezier pieces (CurvePieces.hpp), the
/// surface into enclosed Bezier patches (BezierPatch.hpp), and every pair of
/// a piece and a patch is settled as one of:
///
/// - apart: the boxes of the piece's poles and of the patch's net are
///   disjoint, so the convex hulls are;
/// - without a root, or with exactly one, by Krawczyk's operator on
///
///     F(s, u, v) = C(s) - S(u, v) = 0,   J = [C'(s), -S_u, -S_v],
///
///   over the pair's parameter box. A unique root is a transverse crossing:
///   Krawczyk's contraction proves J regular over the box, so the curve's
///   tangent is not in the surface's tangent plane there;
/// - or cut: the piece or the patch, whichever is larger in space.
///
/// A pair that cannot be cut further, or a search past its budget, fails
/// with IntersectionFailed: a tangent contact, a curve lying in the surface,
/// or a crossing too nearly tangential to isolate. Never a dropped or an
/// unproved root (ADR-0029). A root on a shared boundary of two pieces or
/// two patches is found by both and reported once. A root whose bracket
/// reaches a requested domain but whose middle lies outside it cannot be
/// placed on either side, and fails as an undecided root does elsewhere.

namespace Stereon::Detail {

namespace {

/// The pairs of a piece and a patch one search may visit.
constexpr std::size_t CrossingBudget = std::size_t{1} << 16;

/// A piece's point and tangent over a parameter interval of its own,
/// enclosed; nothing when its weight is not proved positive.
struct CurveJet {
    Bounds Point;
    Bounds Tangent;
};

[[nodiscard]] std::optional<CurveJet> JetOf(const Arc& piece, const Interval& s)
{
    const HPoint p = EvalH(piece.Poles, s);
    if (!(p.W.GetLower() > 0.0)) {
        return std::nullopt;
    }
    const HPoint d = EvalH(Hodograph(piece.Poles), s);
    const Bounds point{p.X / p.W, p.Y / p.W, p.Z / p.W};
    return CurveJet{point, QuotientOf(d, point, p.W)};
}

[[nodiscard]] std::optional<Bounds> BoxOfPiece(const Arc& piece)
{
    Bounds box{Interval::Empty(), Interval::Empty(), Interval::Empty()};
    for (const HPoint& h : piece.Poles) {
        if (!(h.W.GetLower() > 0.0)) {
            return std::nullopt;
        }
        box = {Hull(box[0], h.X / h.W), Hull(box[1], h.Y / h.W), Hull(box[2], h.Z / h.W)};
    }
    return box;
}

/// Whether the hulls of the piece and the patch are proved disjoint.
[[nodiscard]] bool Apart(const Arc& piece, const Patch& patch)
{
    const std::optional<Bounds> box = BoxOfPiece(piece);
    if (!box || !AllBounded(patch.Extent)) {
        return false;
    }
    for (std::size_t k = 0; k < 3; ++k) {
        if (!Intersects((*box)[k], patch.Extent[k])) {
            return true;
        }
    }
    return false;
}

/// The net's homogeneous value at fractions of its stretch, given directly:
/// mapping a fraction to the surface's parameter and back, as ValueOf
/// would, inflates its enclosure once a patch is narrow beside the
/// parameter's magnitude.
[[nodiscard]] HPoint ValueAtFractions(const Net& net, const Interval& a, const Interval& b)
{
    std::vector<HPoint> across;
    for (std::size_t i = 0; i < net.CountU; ++i) {
        across.push_back(Reduced(LineOf(net, false, i), b));
    }
    return Reduced(std::move(across), a);
}

/// The root of C(s) - S(u, v) in the pair's box. Every variable runs over
/// its own piece or patch, s the piece's parameter and (a, b) the patch's
/// fractions, u = u0 + a wu and v = v0 + b wv: Krawczyk's contraction test
/// is a row sum, so variables of very different scales, a small piece's
/// parameter against a surface's, would never pass it, however transverse
/// the crossing. The root is returned in the surface's parameters.
[[nodiscard]] Isolated<3> PairRoot(const Arc& piece, const Patch& patch)
{
    const Interval unit = Reaching(Interval{0.0, 1.0});
    const std::optional<CurveJet> jet = JetOf(piece, unit);
    if (!jet) {
        return {};
    }
    const Hulls hulls = HullsOf(patch);
    const Interval wu = WidthOf(patch.Parameters.X);
    const Interval wv = WidthOf(patch.Parameters.Y);
    const Interval u0{patch.Parameters.X.GetLower()};
    const Interval v0{patch.Parameters.Y.GetLower()};
    IntervalMatrix<3> j{};
    for (std::size_t r = 0; r < 3; ++r) {
        j[r] = {jet->Tangent[r], Slackened(-hulls.Du[r] * wu), Slackened(-hulls.Dv[r] * wv)};
    }
    Isolated<3> root = Krawczyk<3>({unit, unit, unit}, j, [&](const std::array<double, 3>& m) {
        const std::optional<CurveJet> at = JetOf(piece, Interval{m[0]});
        const HPoint h = ValueAtFractions(patch.Poles, Interval{m[1]}, Interval{m[2]});
        if (!at || !(h.W.GetLower() > 0.0)) {
            return Unbounded();
        }
        return MinusOf(at->Point, CartesianOf(h));
    });
    root.Box[1] = u0 + root.Box[1] * wu;
    root.Box[2] = v0 + root.Box[2] * wv;
    return root;
}

/// A crossing found: its brackets on the curve's piece parameter and on
/// the surface's parameters.
struct Found3 {
    Interval OnCurve;
    Interval U;
    Interval V;
};

struct CrossingSearch {
    const PieceCurve* Curve = nullptr;
    Box2 Domain;
    std::size_t Pairs = 0;
    std::vector<Found3> Found;
    std::string Handler;
    /// Whether the surface closes on itself across the ends of its u or v
    /// domain, and the whole of that domain is searched: its first and last
    /// rows of poles are the same.
    bool SeamU = false;
    bool SeamV = false;
};

/// Whether the surface's first and last rows of poles across one
/// direction, weights included, are the same: it closes on itself there,
/// as a full turn of a revolved surface does.
[[nodiscard]] bool ClosedAcross(const BSplineSurface& s, const bool alongU)
{
    const std::size_t nu = PoleCountU(s);
    const std::size_t nv = PoleCountV(s);
    const std::size_t rows = alongU ? nv : nu;
    const auto index = [&](const std::size_t first, const std::size_t k) {
        return alongU ? (first == 0 ? 0 : nu - 1) * nv + k : k * nv + (first == 0 ? 0 : nv - 1);
    };
    for (std::size_t k = 0; k < rows; ++k) {
        const std::size_t a = index(0, k);
        const std::size_t b = index(1, k);
        if (s.Poles[a] != s.Poles[b] || (!s.Weights.empty() && s.Weights[a] != s.Weights[b])) {
            return false;
        }
    }
    return true;
}

/// A surface parameter's middle placed in its domain: as it is inside; at
/// the nearer end when the bracket reaches across a seam; nothing when it
/// cannot be placed.
[[nodiscard]] std::optional<double> OnSurface(const Interval& bracket, const Interval& domain,
                                              const bool seam)
{
    const double x = Midpoint(bracket);
    if (Contains(domain, x)) {
        return x;
    }
    if (seam) {
        return x < domain.GetLower() ? domain.GetLower() : domain.GetUpper();
    }
    return std::nullopt;
}

[[nodiscard]] KernelError Unresolved(const CrossingSearch& search, const Arc& piece,
                                     const Patch& patch, const std::string_view why)
{
    const Segment& segment = search.Curve->Segments[piece.Index];
    const double t = Midpoint(OnCurve(segment, InSegment(piece, Interval{0.5})));
    return Fail(ErrorCode::IntersectionFailed,
                search.Handler + ": " + std::string(why) + " near the curve's " + Number(t) +
                    " and the surface's (" + Number(Midpoint(patch.Parameters.X)) + ", " +
                    Number(Midpoint(patch.Parameters.Y)) +
                    ") (a tangent contact, the curve along the surface or a nearly tangential "
                    "crossing; ADR-0029: an unresolved pair fails, it is never dropped)")
        .error();
}

/// Records the root of an isolated box: kept when it reaches the piece's
/// stretch in the domain and the surface's requested box; on a neighbour's
/// stretch, that neighbour proves it.
void Record(CrossingSearch& search, const Arc& piece, const Patch& patch, const IntervalsN<3>& root)
{
    const Segment& segment = search.Curve->Segments[piece.Index];
    const Interval onCurve = OnCurve(segment, InSegment(piece, root[0]));
    const Interval stretch = Intersection(Interval{segment.T0, segment.T1}, search.Curve->Domain);
    if (!Intersects(onCurve, stretch) || !Intersects(root[1], Reaching(patch.Parameters.X)) ||
        !Intersects(root[2], Reaching(patch.Parameters.Y))) {
        return;
    }
    if (!Intersects(root[1], search.Domain.X) || !Intersects(root[2], search.Domain.Y)) {
        return;
    }
    search.Found.push_back({onCurve, root[1], root[2]});
}

/// The larger of the piece and the patch, in space, cut in two; false when
/// neither can be.
[[nodiscard]] bool PushHalves(std::vector<std::pair<Arc, Patch>>& stack, const Arc& piece,
                              const Patch& patch, const int depth)
{
    const std::optional<Bounds> box = BoxOfPiece(piece);
    const bool pieceLarger =
        !box || !AllBounded(patch.Extent) || LargestWidth(*box) >= LargestWidth(patch.Extent);
    const auto cutPiece = [&] {
        if (depth >= MaxDepth) {
            return false;
        }
        const auto [left, right] = Halves(piece);
        stack.emplace_back(left, patch);
        stack.emplace_back(right, patch);
        return true;
    };
    const auto cutPatch = [&] {
        const std::optional<std::pair<Patch, Patch>> halves = Halves(patch);
        if (!halves) {
            return false;
        }
        stack.emplace_back(piece, halves->first);
        stack.emplace_back(piece, halves->second);
        return true;
    };
    return pieceLarger ? (cutPiece() || cutPatch()) : (cutPatch() || cutPiece());
}

/// The depth of a piece: how often it has been halved, from its width.
[[nodiscard]] int DepthOf(const Arc& piece)
{
    return static_cast<int>(std::round(-std::log2(piece.S1 - piece.S0)));
}

[[nodiscard]] Result<void> SearchCrossings(CrossingSearch& search,
                                           std::vector<std::pair<Arc, Patch>> stack)
{
    while (!stack.empty()) {
        const auto [piece, patch] = std::move(stack.back());
        stack.pop_back();
        if (++search.Pairs > CrossingBudget) {
            return std::unexpected(
                Unresolved(search, piece, patch, "the subdivision budget is spent"));
        }
        if (Apart(piece, patch)) {
            continue;
        }
        const Isolated<3> root = PairRoot(piece, patch);
        if (root.Outcome == Isolation::None) {
            continue;
        }
        if (root.Outcome == Isolation::Unique) {
            Record(search, piece, patch, root.Box);
            continue;
        }
        if (!PushHalves(stack, piece, patch, DepthOf(piece))) {
            return std::unexpected(
                Unresolved(search, piece, patch, "a pair neither apart nor a single crossing"));
        }
    }
    return {};
}

/// Crossings whose brackets meet on the curve and the surface are one,
/// found from two pieces or two patches sharing a boundary.
[[nodiscard]] std::vector<Found3> Distinct(std::vector<Found3> found)
{
    std::ranges::stable_sort(found, [](const Found3& p, const Found3& q) {
        return p.OnCurve.GetLower() < q.OnCurve.GetLower();
    });
    std::vector<Found3> distinct;
    for (const Found3& c : found) {
        const auto same = std::ranges::find_if(distinct, [&](const Found3& d) {
            return Intersects(Reaching(d.OnCurve), c.OnCurve) && Intersects(Reaching(d.U), c.U) &&
                   Intersects(Reaching(d.V), c.V);
        });
        if (same == distinct.end()) {
            distinct.push_back(c);
        } else {
            same->OnCurve = Hull(same->OnCurve, c.OnCurve);
            same->U = Hull(same->U, c.U);
            same->V = Hull(same->V, c.V);
        }
    }
    return distinct;
}

/// The surface's data checked as MakeBSplineSurface does.
[[nodiscard]] Result<void> CheckSurface(const BSplineSurface& s, const std::string_view where)
{
    STEREON_TRY_VOID(CheckKnotVector(s.DegreeU, s.KnotsU, PoleCountU(s), where));
    STEREON_TRY_VOID(CheckKnotVector(s.DegreeV, s.KnotsV, PoleCountV(s), where));
    if (s.Poles.size() != PoleCountU(s) * PoleCountV(s)) {
        return Fail(ErrorCode::InvalidInput, std::string(where) + ": the pole grid does not fit");
    }
    STEREON_TRY_VOID(CheckWeights(s.Weights, s.Poles.size(), where));
    for (const Point3& pole : s.Poles) {
        if (!std::isfinite(pole.X) || !std::isfinite(pole.Y) || !std::isfinite(pole.Z)) {
            return Fail(ErrorCode::InvalidInput, std::string(where) + ": poles must be finite");
        }
    }
    return {};
}

/// A ball holding the patches, widened past the resolution, as BallOf does
/// for a curve's pieces.
[[nodiscard]] Ball BallOf(const std::vector<Patch>& patches, const double resolution)
{
    Bounds box{Interval::Empty(), Interval::Empty(), Interval::Empty()};
    for (const Patch& p : patches) {
        for (std::size_t k = 0; k < 3; ++k) {
            box[k] = Hull(box[k], p.Extent[k]);
        }
    }
    const Point3 center = Middle(box);
    const double radius =
        0.5 * std::hypot(Width(box[0]), Width(box[1]), Width(box[2])) * (1.0 + 0x1p-30);
    const double scale = std::max({1.0, radius, Length(center - Point3{})});
    return {center, radius + 2.0 * resolution + 0x1p-30 * scale};
}

/// The curve's pieces over its requested domain, an unbounded one clipped
/// to the ball of the surface's patches; nothing when the clip is empty.
[[nodiscard]] Result<std::optional<PieceCurve>> CurvePiecesNear(const Curve& curve,
                                                                const Interval& domain,
                                                                const std::vector<Patch>& patches,
                                                                const double resolution)
{
    STEREON_TRY(const Interval requested, RequestedOf(curve, domain));
    Interval searched = requested;
    if (!requested.IsEmpty() && !requested.IsBounded()) {
        STEREON_TRY(const Interval reach, Reach(curve, BallOf(patches, resolution)));
        searched = Intersection(requested, reach);
    }
    if (searched.IsEmpty()) {
        return std::optional<PieceCurve>{};
    }
    STEREON_TRY(PieceCurve pieces, PiecesOf(curve, searched));
    return std::optional<PieceCurve>{std::move(pieces)};
}

/// The crossings as points, in order of the curve's parameter. A root whose
/// middle lies outside a requested domain while its bracket reaches it is
/// undecided: it fails, as RootsOf's do, except across a seam of the
/// surface, where it is placed at the end. A root at a seam, found from the
/// patches on both sides, is one point.
[[nodiscard]] Result<CurveSurfaceIntersection> PointsOf(const CrossingSearch& search,
                                                        const Curve& curve, const double resolution)
{
    CurveSurfaceIntersection result;
    std::vector<Interval> brackets;
    const PieceCurve& c = *search.Curve;
    for (const Found3& f : Distinct(search.Found)) {
        const double u = Midpoint(f.OnCurve);
        const std::optional<double> su = OnSurface(f.U, search.Domain.X, search.SeamU);
        const std::optional<double> sv = OnSurface(f.V, search.Domain.Y, search.SeamV);
        if (!Contains(c.Domain, u) || !su || !sv) {
            return Fail(ErrorCode::IntersectionFailed,
                        search.Handler +
                            ": a crossing within rounding of a requested domain's end, which "
                            "cannot be placed on either side of it");
        }
        const double t = CurveParameter(c, u);
        const Point3 p = Eval(curve, t);
        const bool repeat = std::ranges::any_of(
            std::views::iota(std::size_t{0}, brackets.size()), [&](const std::size_t i) {
                return Intersects(Reaching(brackets[i]), f.OnCurve) &&
                       Distance(result.Points[i].Position, p) <= resolution;
            });
        if (!repeat) {
            brackets.push_back(f.OnCurve);
            result.Points.push_back({t, Point2{*su, *sv}, p, ContactType::Crossing, true});
        }
    }
    std::ranges::sort(result.Points, [](const CurveSurfacePoint& p, const CurveSurfacePoint& q) {
        return p.Parameter < q.Parameter;
    });
    return result;
}

} // namespace

Result<CurveSurfaceIntersection> CurveSplineSurface(const Curve& a, const Interval& domain,
                                                    const Surface& b, const Box2& box,
                                                    const Context& ctx)
{
    const BSplineSurface* surface = b.Target<BSplineSurface>();
    if (surface == nullptr) {
        return Fail(ErrorCode::InvalidInput,
                    "CurveSplineSurface: dispatched for another surface type");
    }
    const std::string handler = "Curve x BSplineSurface";
    STEREON_TRY_VOID(CheckSurface(*surface, handler));
    const Box2 own = Domain(*surface);
    const Box2 searched{Intersection(own.X, box.X), Intersection(own.Y, box.Y)};
    if (searched.X.IsEmpty() || searched.Y.IsEmpty()) {
        return CurveSurfaceIntersection{};
    }
    STEREON_TRY(const std::vector<Patch> patches, PatchesOver(*surface, searched));
    STEREON_TRY(const std::optional<PieceCurve> pieces,
                CurvePiecesNear(a, domain, patches, ctx.GetResolution()));
    if (!pieces) {
        return CurveSurfaceIntersection{};
    }
    CrossingSearch search{&*pieces, searched, 0, {}, pieces->Name + " x BSplineSurface"};
    search.SeamU = Contains(box.X, own.X) && ClosedAcross(*surface, true);
    search.SeamV = Contains(box.Y, own.Y) && ClosedAcross(*surface, false);
    std::vector<std::pair<Arc, Patch>> stack;
    for (std::size_t i = 0; i < pieces->Segments.size(); ++i) {
        for (const Patch& patch : patches) {
            stack.emplace_back(Arc{i, 0.0, 1.0, pieces->Segments[i].Poles}, patch);
        }
    }
    STEREON_TRY_VOID(SearchCrossings(search, std::move(stack)));
    return PointsOf(search, a, ctx.GetResolution());
}

} // namespace Stereon::Detail
