// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// A transverse node on a cut of the requested domain, or on an end of a
// surface's domain that does not close, for the generic surface-surface
// intersection. Included once, in GenericSSI.cpp, after TouchPoint.hpp,
// whose helpers and those of "Branch points" it uses; it is a part of that
// unit kept in its own file.
//
// A node whose parameters on S lie on a side of the searched domain is in
// no box the proof inside a cell can be centred on. The searched cell of S
// is therefore continued past that side: its Bezier net, polynomial or
// homogeneous, is the same polynomial past the cut (de Casteljau's algorithm
// at a parameter outside [0, 1], in interval arithmetic), a surface the proof
// may use though the result may not. The node is proved on the continued
// patch exactly as inside a cell (GenericSSI.cpp, "Branch points", facts 1
// to 4): four half-branches leave s*, each ending at one of the four roots
// on the centred box's edges, Q = (s - s*)^T M (s - s*) strictly increasing
// along each away from s*. Then, for each cut line through the centred box:
//
// 5. The cut keeps clear. On the cut, by Taylor's theorem about s* with
//    grad h(s*) = 0, h(s) lies in h(s*) + (s - s*)^T H (s - s*) / 2 for H in
//    the Hessian's enclosure; this is proved to exclude zero on every piece
//    of the line farther from s* than a window delta, the window's reach
//    delta |S_s| in space within the resolution and |h| over it within the
//    tolerance.
//
// What follows. A half-branch crosses the cut only within the window. Q
// increases along it, so once past the window's largest Q it never returns:
// its far part lies on the side of its end, and its near part within the
// window, which the surfaces meet within the tolerance and the resolution
// of the node. The half-branches whose ends lie in the searched domain are
// therefore the intersection there, joined at the node: reported as one
// branch from the node, clamped onto the cut, to each such end, the node a
// candidate for an event (ADR-0054) on the domain's boundary.
//
// The search: a pair whose cell of S lies on a side of the searched domain
// that does not close (a requested cut, or an end of an open chart), a
// candidate for a node on the continued pair, the boxes narrowed about it
// as for a node inside a cell; the product of the centred boxes, cut back to
// the searched domain, recorded as a settled region where nothing was found
// yet, and carved from the search.
//
// Not settled, and failing closed: a cut along a direction the Hessian's
// form does not keep of one sign (a branch along the cut), a node on cuts of
// both surfaces at once, a node on a knot line or seam as well as a cut, and
// a window wider than the resolution allows.

namespace Stereon::Detail {
namespace {

/// The sides of the searched domain a cell lies on, in the order u low, u
/// high, v low, v high, where the domain does not close: a requested cut,
/// or an end of a chart that is not proved closed. A seam of a closed chart
/// searched whole is SeamNode.hpp's.
using CutSides = std::array<bool, 4>;

/// The cut sides of the searched domain a cell lies on.
[[nodiscard]] CutSides CutSidesOf(const Atlas& atlas, const Box2& cell)
{
    CutSides sides{};
    for (const bool alongU : {true, false}) {
        const Interval mine = alongU ? cell.X : cell.Y;
        const Interval searched = alongU ? atlas.Searched.X : atlas.Searched.Y;
        const Interval whole = alongU ? atlas.Domain.X : atlas.Domain.Y;
        const bool closes = (alongU ? atlas.ClosedU : atlas.ClosedV) && searched == whole;
        const std::size_t k = alongU ? 0 : 2;
        sides.at(k) = !closes && mine.GetLower() == searched.GetLower();
        sides.at(k + 1) = !closes && mine.GetUpper() == searched.GetUpper();
    }
    return sides;
}

/// The patch continued past one side of its stretch to `to`: de Casteljau's
/// algorithm on each line at t = (to - a) / (b - a), outside [0, 1], the part
/// over [a, to] or [to, b]. The same polynomial, or rational map with its
/// homogeneous net, past the side; a weight there not proved positive leaves
/// an extent no proof can use (ExtentOf).
[[nodiscard]] std::optional<Patch> ContinuedTo(const Patch& p, const bool alongU, const double to)
{
    const Interval stretch = StretchOf(p, alongU);
    const Interval width = Interval{stretch.GetUpper()} - Interval{stretch.GetLower()};
    const bool upper = to > stretch.GetUpper();
    if (!(width.GetLower() > 0.0) || (!upper && !(to < stretch.GetLower()))) {
        return std::nullopt;
    }
    const Interval t = (Interval{to} - Interval{stretch.GetLower()}) / width;
    std::vector<std::vector<HPoint>> lines(LineCount(p.Poles, alongU));
    for (std::size_t line = 0; line < lines.size(); ++line) {
        std::vector<HPoint> left;
        std::vector<HPoint> right;
        SplitLine(LineOf(p.Poles, alongU, line), t, left, right);
        lines[line] = upper ? std::move(left) : std::move(right);
    }
    Box2 box = p.Parameters;
    (alongU ? box.X : box.Y) =
        upper ? Interval{stretch.GetLower(), to} : Interval{to, stretch.GetUpper()};
    return MakePatch(box, WithLines(p.Poles, alongU, lines));
}

/// The cell continued past each of its cut sides by its own width there.
[[nodiscard]] std::optional<Patch> ContinuedPast(Patch p, const CutSides& sides)
{
    for (std::size_t k = 0; k < 4; ++k) {
        if (!sides.at(k)) {
            continue;
        }
        const bool alongU = k < 2;
        const Interval stretch = StretchOf(p, alongU);
        const double to =
            k % 2 == 0 ? stretch.GetLower() - Width(stretch) : stretch.GetUpper() + Width(stretch);
        std::optional<Patch> next = ContinuedTo(p, alongU, to);
        if (!next) {
            return std::nullopt;
        }
        p = std::move(*next);
    }
    return p;
}

/// How deep a cut line is halved while its pieces are tested.
constexpr int CutDepth = 12;

/// The height on a piece of a cut line, by Taylor's theorem about the
/// node's box X: h* + (s - X)^T H (s - X) / 2, s on the piece.
[[nodiscard]] Interval HeightOnCut(const NodeAt& at, const bool alongU, const double cut,
                                   const Interval& piece)
{
    const Interval du = alongU ? Interval{cut} - at.Box[0] : piece - at.Box[0];
    const Interval dv = alongU ? piece - at.Box[1] : Interval{cut} - at.Box[1];
    const Matrix2& h = at.Hessian;
    return at.Height + Interval{0.5} * (h[0][0] * Squared(du) + Interval{2.0} * h[0][1] * du * dv +
                                        h[1][1] * Squared(dv));
}

/// Fact 5 on the cut a = cut (a = u when alongU) across the centred box:
/// every piece of the line keeps h off zero, or lies within `window` of the
/// node along the line with |h| within the tolerance there.
[[nodiscard]] bool CutKeepsClear(const NodeAt& at, const bool alongU, const double cut,
                                 const Interval& line, const double window, const double tolerance)
{
    const Interval node = alongU ? at.Box[1] : at.Box[0];
    std::vector<std::pair<Interval, int>> pieces{{line, 0}};
    while (!pieces.empty()) {
        const auto [piece, depth] = pieces.back();
        pieces.pop_back();
        const Interval h = HeightOnCut(at, alongU, cut, piece);
        if (ExcludesZero(h)) {
            continue;
        }
        const bool near = piece.GetLower() >= node.GetLower() - window &&
                          piece.GetUpper() <= node.GetUpper() + window;
        if (near && h.GetLower() >= -tolerance && h.GetUpper() <= tolerance) {
            continue;
        }
        if (depth >= CutDepth) {
            return false;
        }
        const double middle = Midpoint(piece);
        pieces.push_back({Interval{piece.GetLower(), middle}, depth + 1});
        pieces.push_back({Interval{middle, piece.GetUpper()}, depth + 1});
    }
    return true;
}

/// The window along a cut: as far as delta |S_s| in space is the resolution,
/// |S_s| bounded over the centred box.
[[nodiscard]] double WindowOf(const Jet& overS, const double resolution)
{
    const Interval speed = Sqrt(DotOf(overS.Pu, overS.Pu) + DotOf(overS.Pv, overS.Pv));
    return speed.GetUpper() > 0.0 ? (Interval{resolution} / Interval{speed.GetUpper()}).GetLower()
                                  : 0.0;
}

/// Fact 5 on every cut line of the searched domain through the centred box.
[[nodiscard]] bool CutsKeepClear(const NodeAt& at, const Patch& small, const Box2& searched,
                                 const CutSides& sides, const SearchSettings& s)
{
    const std::optional<Jet> overS = JetOver(small, SecondNetsOf(small));
    if (!overS) {
        return false;
    }
    const double window = WindowOf(*overS, 0x1p7 * s.Join);
    for (std::size_t k = 0; k < 4; ++k) {
        const bool alongU = k < 2;
        const Interval side = alongU ? searched.X : searched.Y;
        const double cut = k % 2 == 0 ? side.GetLower() : side.GetUpper();
        const Interval across = StretchOf(small, alongU);
        if (!sides.at(k) || !(across.GetLower() < cut && cut < across.GetUpper())) {
            continue;
        }
        if (!CutKeepsClear(at, alongU, cut, StretchOf(small, !alongU), window, s.Tolerance)) {
            return false;
        }
    }
    return true;
}

/// The height at the node, h* = n . (S(x*) - C(c*)), in the centred form
/// about the middle m of the node's box X: n . (S(m) - C(m_C)) plus the
/// partials along n, enclosed over X, times X - m. Near a node n is nearly
/// normal to both surfaces, so the partials along it nearly vanish and the
/// box's width hardly widens h*, as it does evaluated over X directly. Met
/// with that direct enclosure; both hold h*.
[[nodiscard]] Interval CentredHeight(const PatchPair& pair, const NodeFrame& f, const NodeAt& at)
{
    const std::array<double, 4> m = MidpointsOf<4>(at.Box);
    const std::optional<Jet> s = JetAt(pair.A, SecondNetsOf(pair.A), at.Box[0], at.Box[1]);
    const std::optional<Jet> c = JetAt(pair.B, SecondNetsOf(pair.B), at.Box[2], at.Box[3]);
    const Bounds sm = XyzOver(pair.A, Interval{m[0]}, Interval{m[1]});
    const Bounds cm = XyzOver(pair.B, Interval{m[2]}, Interval{m[3]});
    if (!s || !c || !AllBounded(sm) || !AllBounded(cm)) {
        return at.Height;
    }
    const auto off = [&](const std::size_t k) { return at.Box[k] - Interval{m[k]}; };
    const Interval h = DotOf(f.N, MinusOf(sm, cm)) + DotOf(f.N, s->Pu) * off(0) +
                       DotOf(f.N, s->Pv) * off(1) - DotOf(f.N, c->Pu) * off(2) -
                       DotOf(f.N, c->Pv) * off(3);
    const Interval met = Intersection(h, at.Height);
    return met.IsEmpty() ? at.Height : met;
}

/// The node's parameters on S clamped into the searched domain, its point
/// S's there.
[[nodiscard]] Seed ClampedNode(const Seed& node, const Patch& small, const Box2& searched)
{
    const Point2 on = ClampedInto(node.OnA, searched);
    const Bounds p = XyzOver(small, Interval{on.X}, Interval{on.Y});
    return {AllBounded(p) ? Middle(p) : node.Position, on, node.OnB};
}

/// A proved node on a cut: the node, the ends of its half-branches in the
/// searched domain, and the boxes of the region it settles, in the order
/// (S, C).
struct CutNode {
    Seed Node;
    std::vector<Seed> Ends;
    Box2 OnS;
    Box2 OnC;
};

/// The node proved on the centred pair (the small patch continued past the
/// cut), its cuts kept clear, and its ends in the searched domain; nothing
/// unless at least one end is.
[[nodiscard]] std::optional<CutNode> ProvedCutNode(const PatchPair& centred, const Vec3& normal,
                                                   const Box2& searched, const CutSides& sides,
                                                   const SearchSettings& s)
{
    std::optional<NodeProof> proof = ProvedNodeAt(centred, normal, s);
    const std::optional<NodeFrame> frame = FrameAbout(normal);
    if (!proof || !frame) {
        return std::nullopt;
    }
    proof->At.Height = CentredHeight(centred, *frame, proof->At);
    if (!CutsKeepClear(proof->At, centred.A, searched, sides, s)) {
        return std::nullopt;
    }
    std::vector<Seed> ends;
    for (const Seed& end : proof->Ends.Ends) {
        if (Contains(searched, end.OnA)) {
            ends.push_back(end);
        }
    }
    const Interval u = Intersection(centred.A.Parameters.X, searched.X);
    const Interval v = Intersection(centred.A.Parameters.Y, searched.Y);
    if (ends.empty() || ends.size() == proof->Ends.Ends.size() || u.IsEmpty() || v.IsEmpty()) {
        return std::nullopt;
    }
    return CutNode{ClampedNode(proof->Ends.Node, centred.A, searched), std::move(ends), Box2{u, v},
                   centred.B.Parameters};
}

/// The node, its branches and its region into the findings, in the order
/// (A, B); the pair goes back on the stack to be carved.
void RecordCutNode(const CutNode& n, const bool smallIsA, Findings& found)
{
    const auto ordered = [&](const Seed& seed) { return smallIsA ? seed : Swapped(seed); };
    const Seed node = ordered(n.Node);
    const Box2& boxA = smallIsA ? n.OnS : n.OnC;
    const Box2& boxB = smallIsA ? n.OnC : n.OnS;
    for (const Seed& e : n.Ends) {
        const Seed end = ordered(e);
        const Vec3 d = ToVec(end.Position) - ToVec(node.Position);
        found.Arcs.push_back({node, end, d / Length(d), boxA, boxB});
    }
    found.Nodes.push_back(node.Position);
    found.Regions.push_back({{boxA}, {boxB}});
}

/// The candidate's covering cell narrowed j times, and within it the small
/// continued cell narrowed until the node is proved, or nothing.
[[nodiscard]] std::optional<CutNode> CutNodeAround(const PatchPair& continued, const Placed& at,
                                                   const Box2& searched, const CutSides& sides,
                                                   const SearchSettings& s)
{
    for (int j = 0; j < CoverNarrowings; ++j) {
        const double fc = std::ldexp(0.5, -j);
        const std::optional<Centred> cover =
            CentredOn(continued.B, Around(at.OnC, fc * at.MarginC[0], fc * at.MarginC[1]));
        if (!cover) {
            continue;
        }
        const Vec3 normal =
            NormalAt(cover->Centre, SecondNetsOf(cover->Centre), at.OnC.X, at.OnC.Y);
        for (int k = 0; k < SmallNarrowings; ++k) {
            const double fs = std::ldexp(0.5, -k);
            const std::optional<Centred> small =
                CentredOn(continued.A, Around(at.OnS, fs * at.MarginS[0], fs * at.MarginS[1]));
            std::optional<CutNode> node =
                small ? ProvedCutNode({small->Centre, cover->Centre}, normal, searched, sides, s)
                      : std::nullopt;
            if (node) {
                return node;
            }
        }
    }
    return std::nullopt;
}

/// The candidate on the pair with S's cell continued past its cut sides,
/// placed, or nothing; the node must lie in the cell, within a 2^-30th of
/// its width, not past the cut.
[[nodiscard]] std::optional<Placed> CutCandidate(const PatchPair& continued, const Box2& cell,
                                                 const Box2& cellC, const SearchSettings& s)
{
    const std::optional<std::array<double, 4>> m = NodeCandidate(continued, s);
    if (!m) {
        return std::nullopt;
    }
    const Point2 onS{(*m)[0], (*m)[1]};
    const Point2 onC{(*m)[2], (*m)[3]};
    const Box2 slack{Inflated(cell.X, 0x1p-30 * Width(cell.X)),
                     Inflated(cell.Y, 0x1p-30 * Width(cell.Y))};
    const Placed at{onS, onC, Margins(onS, continued.A.Parameters), Margins(onC, cellC)};
    if (!Contains(slack, onS) || !(std::min(at.MarginS[0], at.MarginS[1]) > 0.0) ||
        !(std::min(at.MarginC[0], at.MarginC[1]) > 0.0)) {
        return std::nullopt;
    }
    return at;
}

/// The node on a cut side of S's cell, S the surface of A when smallIsA, or
/// nothing.
[[nodiscard]] std::optional<CutNode> CutNodeOf(const PatchPair& pair, const SearchSettings& s,
                                               const bool smallIsA)
{
    const PatchPair oriented = smallIsA ? pair : PatchPair{pair.B, pair.A};
    const Atlas& atlas = smallIsA ? *s.OnA : *s.OnB;
    const CutSides sides = CutSidesOf(atlas, oriented.A.Parameters);
    if (std::ranges::none_of(sides, [](const bool b) { return b; })) {
        return std::nullopt;
    }
    std::optional<Patch> continued = ContinuedPast(oriented.A, sides);
    if (!continued) {
        return std::nullopt;
    }
    const PatchPair wide{std::move(*continued), oriented.B};
    const std::optional<Placed> at =
        CutCandidate(wide, oriented.A.Parameters, oriented.B.Parameters, s);
    return at ? CutNodeAround(wide, *at, atlas.Searched, sides, s) : std::nullopt;
}

/// A pair whose cell of S lies on a cut side of the searched domain that may
/// hold a node there, either surface taken as S; the pair goes back on the
/// stack to be carved by the node's region.
[[nodiscard]] bool SettleCutNode(const PatchPair& pair, const SearchSettings& s, Findings& found,
                                 std::vector<PatchPair>& stack)
{
    if (s.OnA == nullptr || s.OnB == nullptr) {
        return false;
    }
    for (const bool smallIsA : {true, false}) {
        const std::optional<CutNode> node = CutNodeOf(pair, s, smallIsA);
        if (node && NothingFoundIn(smallIsA ? node->OnS : node->OnC,
                                   smallIsA ? node->OnC : node->OnS, found)) {
            RecordCutNode(*node, smallIsA, found);
            stack.push_back(pair);
            return true;
        }
    }
    return false;
}

} // namespace
} // namespace Stereon::Detail
