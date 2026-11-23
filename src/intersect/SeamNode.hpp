// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Branch points on knot lines and seams, for the generic surface-surface
// intersection. Included once, in GenericSSI.cpp, after the proof of a node
// inside one Bezier cell of each surface ("Branch points" there), whose
// helpers it uses; it is a part of that unit kept in its own file.
//
// A transverse node whose parameters lie on a knot line of a surface, or on
// the seam where its chart closes, is in no single Bezier patch: the patches
// around it meet there, and across the line the surface need not be smooth
// in its parameters. A full turn of a rational quadratic circle is only
// continuous across its double knots, and with its middle weights rounded
// to doubles it has a crease of a few roundings there. No enclosure can
// prove such a crease absent, so the proof allows one, and is tolerant where
// the proof inside a cell is exact.
//
// The setting. S and C are the two surfaces near the node, each a chart of
// one, two or four clipped Bezier pieces around a corner, c_S and c_C, a
// point all of a chart's pieces share; the pieces tile a box, R for S, in
// the chart's parameters, a piece beyond a seam relabelled to lie next to
// its neighbour (an affine change of its parameter, so nothing is assumed
// of the two sides' speeds).
//
// 1. Coordinates. In a frame (e1, e2, n), C's projection q onto the plane is
//    continuous on its box (across a seam because the chart is proved
//    closed, ClosesAlong), and its Jacobian, enclosed over every piece, is
//    regular in the hull of those enclosures. By the mean value theorem
//    along a segment, piece by piece, no two parameters share a projection:
//    C is a graph over the plane.
// 2. Containment. S's projection misses the projection of the outer edges
//    of C's box, and S's corner projects into C's box (Poincare-Miranda on
//    q, preconditioned), so all of S over R lies over C: the height h(s) of
//    S above C along n is defined on R, continuous, and smooth where s and
//    the point c'(s) of C below it stay in one piece each. Its zeros are
//    exactly the points of S over R on C.
// 3. The model. For s in a piece i of S with c'(s) in a piece k of C, the
//    gradient of h is Phi(s, c') = (nu . S_u, nu . S_v), nu the frame's
//    normal less C's slope (CoNormal). Phi and E = e . (S - C) are smooth on
//    the product of the two pieces' boxes, a convex box with the corners
//    (c_S, c_C) in its closure, so by the mean value theorem there, with
//    E = 0 at (s, c'(s)),
//        grad h(s) in G + H (s - c_S),
//    H the Hessian's enclosure of the proof inside a cell and G =
//    Phi(c_S, c_C) - Dp^T Hess g E(c_S, c_C), each enclosed on every pair of
//    pieces from that pair's own jets, and the hull taken over the pairs: G
//    is the common tangent plane's and the corners' residuals, rounding and
//    crease together. Likewise h(c_S) is enclosed.
// 4. A small box D around c_S is given up: nothing is proved of the
//    intersection inside it, only that S stays within the tolerance of C
//    all over it (|h|, on the model). It is the smallest of a ladder of
//    boxes on which facts 4 to 6 hold: the corners' residuals need room
//    around them, a gap g between the surfaces about sqrt(g / curvature).
//    Outside it, for every s in R and every gradient g that h has there
//    (one per pair of pieces meeting at s),
//        (s - c_S)^T M J g > 0,
//    M the Lyapunov matrix of the proof inside a cell: the quadratic part
//    is at least lambda |y|^2 and the residual at most |w| |y| in
//    coordinates y scaled to D, and lambda > |w| is proved. So zero is not
//    in h's generalized gradient, the zero set is a Lipschitz curve through
//    each of its points outside D (Clarke's implicit function theorem), and
//    Q(s) = (s - c_S)^T M (s - c_S) is strictly monotone along it: no closed
//    loop, no arc returning to its start.
// 5. Four rays from c_S, along two directions on which the Hessian's form
//    is proved positive and negative, carry no zero of h outside D: along a
//    ray h's derivative is e . G + t e^T H e, of one sign from D on, and h
//    has that sign where the ray leaves D. They cut R less D into four
//    sectors, h changing sign across each.
// 6. Counting. h has exactly four zeros on D's boundary, counted on the
//    model (monotone pieces with proved signs at their ends), and exactly
//    four on R's, each isolated by Krawczyk's operator as an edge's root on
//    C. Each sector's stretch of either boundary holds an odd number, so
//    one each.
//
// Hence in each sector the zero set is one arc, from D's boundary to R's:
// its components are arcs along which Q is strictly monotone, so their two
// ends are distinct zeros on the sector's boundary, of which there are two.
// The intersection over R is those four arcs and whatever lies over D,
// where the surfaces are within the tolerance of each other: a gap or a
// crease within it is not told apart from a crossing (ADR-0050). The node is
// reported as the shared end of four branches, a candidate for an
// arrangement event (ADR-0054) as one inside a cell is, and no more: the
// branches are not followed through D, and each is certified against both
// surfaces like any other.
//
// The search. The first pair of patches that holds a candidate for such a
// node (Newton's method in doubles, only a guess) settles it once, for
// every pair of patches around it: the product of the two boxes is recorded
// as a region of the search, and every pair of patches met afterwards is
// cut along the region's sides and its part inside dropped (Carved). A
// region is recorded only where the search has found nothing yet
// (NothingFoundOn): no arc of a pair of patches sharing area with it, no
// earlier region. In a depth-first search that is the usual case for charts
// that are not too wide: the other pairs of Bezier patches around the node
// are still whole, or would have found it first, and where an arc was
// found beside the node already, a narrower chart is tried.

#pragma once

namespace Stereon::Detail {
namespace {

// --- Charts across knot lines and seams -------------------------------------

/// Whether a knot vector is clamped: its first and last p + 1 knots equal.
[[nodiscard]] bool Clamped(const std::vector<double>& knots, const std::size_t p)
{
    const std::size_t n = knots.size();
    return n >= 2 * (p + 1) && knots[p] == knots.front() && knots[n - 1 - p] == knots.back();
}

/// Whether two poles of a surface are the same doubles, weights included.
[[nodiscard]] bool SamePole(const BSplineSurface& s, const std::size_t a, const std::size_t b)
{
    const Point3& p = s.Poles[a];
    const Point3& q = s.Poles[b];
    return p.X == q.X && p.Y == q.Y && p.Z == q.Z &&
           (s.Weights.empty() || s.Weights[a] == s.Weights[b]);
}

/// Whether the chart closes along a direction, proved on the doubles: the
/// knots clamped, so that the surface's first and last lines of poles
/// control its two ends, and those lines the same poles and weights. The
/// two ends are then one curve, parameter for parameter.
[[nodiscard]] bool ClosesAlong(const BSplineSurface& s, const bool alongU)
{
    const auto p = static_cast<std::size_t>(alongU ? s.DegreeU : s.DegreeV);
    if (!Clamped(alongU ? s.KnotsU : s.KnotsV, p)) {
        return false;
    }
    const std::size_t nu = PoleCountU(s);
    const std::size_t nv = PoleCountV(s);
    if (alongU) {
        for (std::size_t k = 0; k < nv; ++k) {
            if (!SamePole(s, k, (nu - 1) * nv + k)) {
                return false;
            }
        }
        return true;
    }
    for (std::size_t k = 0; k < nu; ++k) {
        if (!SamePole(s, k * nv, k * nv + nv - 1)) {
            return false;
        }
    }
    return true;
}

/// The patch of the atlas whose box holds a box of parameters.
[[nodiscard]] const Patch* Holding(const Atlas& atlas, const Box2& box)
{
    for (const Patch& p : atlas.Patches) {
        if (Contains(p.Parameters.X, box.X) && Contains(p.Parameters.Y, box.Y)) {
            return &p;
        }
    }
    return nullptr;
}

/// The stretch of the Bezier patch beyond a side of a patch of the atlas,
/// along one direction: across a knot line, or, where the chart closes,
/// across its seam, at the domain's other end.
struct Beyond {
    Interval Stretch;
    bool Seam;
};

/// The end of a stretch a neighbouring patch's stretch must have to lie
/// beyond a side of `mine`: that side, or, on a side of a domain that
/// closes, the domain's other end. Nothing on a side of a domain that does
/// not close.
struct Wanted {
    double End;
    bool Seam;
};

[[nodiscard]] std::optional<Wanted> WantedBeyond(const Interval& mine, const Interval& domain,
                                                 const bool closed, const bool upper)
{
    const double side = upper ? mine.GetUpper() : mine.GetLower();
    const double end = upper ? domain.GetUpper() : domain.GetLower();
    if (side != end) {
        return Wanted{side, false};
    }
    if (!closed) {
        return std::nullopt;
    }
    return Wanted{upper ? domain.GetLower() : domain.GetUpper(), true};
}

[[nodiscard]] std::optional<Beyond> BeyondSide(const Atlas& atlas, const Patch& own,
                                               const bool alongU, const bool upper)
{
    const std::optional<Wanted> wanted =
        WantedBeyond(StretchOf(own, alongU), alongU ? atlas.Domain.X : atlas.Domain.Y,
                     alongU ? atlas.ClosedU : atlas.ClosedV, upper);
    if (!wanted) {
        return std::nullopt;
    }
    for (const Patch& p : atlas.Patches) {
        const Interval s = StretchOf(p, alongU);
        const double end = upper ? s.GetLower() : s.GetUpper();
        if (end == wanted->End && StretchOf(p, !alongU) == StretchOf(own, !alongU)) {
            return Beyond{s, wanted->Seam};
        }
    }
    return std::nullopt;
}

/// One side of a chart's corner along one direction: where it lies in the
/// surface's parameters, and in the chart's.
struct Span {
    Interval True;
    Interval Label;
};

/// A chart along one direction: one span around the corner, or two meeting
/// at it, the lower first.
struct Spans {
    std::vector<Span> Sides;
    double Corner = 0.0;
};

/// How a chart is laid along one direction: the fraction of the patch it
/// takes, and how near a side a candidate counts as on it.
struct Reach {
    double Factor;
    double Snap;
};

/// As Around's boxes: shorter on one side, so that symmetric input does not
/// end a branch on a corner.
constexpr double Shorter = 2.0 * SplitFraction;

/// An interval, or nothing when its ends are not in order.
[[nodiscard]] std::optional<Interval> Ordered(const double lower, const double upper)
{
    return lower < upper ? std::optional<Interval>{Interval{lower, upper}} : std::nullopt;
}

/// The span beyond a side at c, of width wn: where it lies in the surface's
/// parameters and in the chart's. Beyond a knot line the two are the same;
/// beyond a seam the span lies at the domain's other end, and is labelled
/// to continue the chart.
[[nodiscard]] std::optional<Span> SpanBeyond(const Beyond& beyond, const bool upper, const double c,
                                             const double wn)
{
    const Interval& b = beyond.Stretch;
    if (upper) {
        const double from = beyond.Seam ? b.GetLower() : c;
        const std::optional<Interval> truth = Ordered(from, std::min(from + wn, b.GetUpper()));
        const std::optional<Interval> label = beyond.Seam ? Ordered(c, c + wn) : truth;
        return truth && label ? std::optional<Span>{Span{*truth, *label}} : std::nullopt;
    }
    const double from = beyond.Seam ? b.GetUpper() : c;
    const std::optional<Interval> truth = Ordered(std::max(from - wn, b.GetLower()), from);
    const std::optional<Interval> label = beyond.Seam ? Ordered(c - wn, c) : truth;
    return truth && label ? std::optional<Span>{Span{*truth, *label}} : std::nullopt;
}

/// The two spans meeting on a side of the patch that found the candidate:
/// its own, a fraction of it, and one beyond, no wider than the patch
/// there.
[[nodiscard]] std::optional<Spans> SpansAcross(const Interval& stretch, const bool upper,
                                               const Beyond& beyond, const double factor)
{
    const double c = upper ? stretch.GetUpper() : stretch.GetLower();
    const double w = factor * Width(stretch);
    const std::optional<Interval> own = upper ? Ordered(c - w, c) : Ordered(c, c + w);
    const std::optional<Span> theirs =
        SpanBeyond(beyond, upper, c, std::min(Shorter * w, Width(beyond.Stretch)));
    if (!own || !theirs) {
        return std::nullopt;
    }
    const Span mine{*own, *own};
    return Spans{upper ? std::vector<Span>{mine, *theirs} : std::vector<Span>{*theirs, mine}, c};
}

/// The one span around a candidate inside the patch along a direction.
[[nodiscard]] std::optional<Spans> SpansAround(const Interval& stretch, const double at,
                                               const double factor)
{
    const double r = factor * std::min(at - stretch.GetLower(), stretch.GetUpper() - at);
    const std::optional<Interval> around = Ordered(at - Shorter * r, at + r);
    if (!(r > 0.0) || !around || !(around->GetLower() < at && at < around->GetUpper())) {
        return std::nullopt;
    }
    return Spans{{Span{*around, *around}}, at};
}

/// The chart's spans along one direction, for a candidate at `at` in the
/// patch that found it, a part of the atlas's patch `own`, which the spans
/// are sized by. A candidate on a side must have a patch beyond it: on a
/// cut of the search, or on a side of the domain that does not close, there
/// is nothing to chart across.
[[nodiscard]] std::optional<Spans> SpansAlong(const Atlas& atlas, const Patch& found,
                                              const Patch& own, const bool alongU, const double at,
                                              const Reach& reach)
{
    const Interval stretch = StretchOf(found, alongU);
    const Interval whole = StretchOf(own, alongU);
    const bool low = at - stretch.GetLower() <= reach.Snap;
    const bool high = stretch.GetUpper() - at <= reach.Snap;
    if (!low && !high) {
        return SpansAround(whole, at, reach.Factor);
    }
    const bool onOwnSide =
        high ? stretch.GetUpper() == whole.GetUpper() : stretch.GetLower() == whole.GetLower();
    if (low == high || !onOwnSide) {
        return std::nullopt;
    }
    const std::optional<Beyond> beyond = BeyondSide(atlas, own, alongU, high);
    return beyond ? SpansAcross(whole, high, *beyond, reach.Factor) : std::nullopt;
}

/// One piece of a chart: a Bezier patch clipped and labelled in the chart's
/// parameters, and the box it has in the surface's own.
struct ChartPiece {
    Patch Labelled;
    Box2 True;
};

/// The pieces around a corner, tiling Box in the chart's parameters.
struct Chart {
    std::vector<ChartPiece> Pieces;
    Point2 Corner;
    Box2 Box;
};

[[nodiscard]] std::optional<ChartPiece> PieceOver(const Atlas& atlas, const Span& u, const Span& v)
{
    const Box2 box{u.True, v.True};
    const Patch* source = Holding(atlas, box);
    if (source == nullptr) {
        return std::nullopt;
    }
    Result<std::optional<Patch>> alongU = Clipped(*source, true, u.True);
    if (!alongU || !*alongU) {
        return std::nullopt;
    }
    Result<std::optional<Patch>> both = Clipped(std::move(**alongU), false, v.True);
    if (!both || !*both) {
        return std::nullopt;
    }
    Patch labelled = std::move(**both);
    labelled.Parameters = Box2{u.Label, v.Label};
    return ChartPiece{std::move(labelled), box};
}

/// The chart around a candidate at `at` in the patch that found it, or
/// nothing where none can be laid.
[[nodiscard]] std::optional<Chart> ChartAbout(const Atlas& atlas, const Patch& found,
                                              const Point2& at, const Reach& u, const Reach& v)
{
    const Patch* own = Holding(atlas, found.Parameters);
    if (own == nullptr) {
        return std::nullopt;
    }
    const std::optional<Spans> su = SpansAlong(atlas, found, *own, true, at.X, u);
    const std::optional<Spans> sv = SpansAlong(atlas, found, *own, false, at.Y, v);
    if (!su || !sv) {
        return std::nullopt;
    }
    Chart chart{{},
                Point2{su->Corner, sv->Corner},
                Box2{Hull(su->Sides.front().Label, su->Sides.back().Label),
                     Hull(sv->Sides.front().Label, sv->Sides.back().Label)}};
    for (const Span& a : su->Sides) {
        for (const Span& b : sv->Sides) {
            std::optional<ChartPiece> piece = PieceOver(atlas, a, b);
            if (!piece) {
                return std::nullopt;
            }
            chart.Pieces.push_back(std::move(*piece));
        }
    }
    return chart;
}

/// The chart with every piece moved by -c in space, enclosed (MovedNet):
/// the same surface in coordinates about c. The quotient rule's enclosures
/// of a rational piece's partials widen with the size of its points'
/// coordinates, so both charts are taken about the node before anything is
/// proved of them; the proof is of the two moved surfaces, which meet as
/// the given ones do.
[[nodiscard]] Chart About(Chart chart, const Point3& c)
{
    for (ChartPiece& piece : chart.Pieces) {
        Patch moved = MakePatch(piece.Labelled.Parameters, MovedNet(piece.Labelled.Poles, c));
        moved.NextAlongU = piece.Labelled.NextAlongU;
        piece.Labelled = std::move(moved);
    }
    return chart;
}

/// A chart's parameter in the surface's own, on one piece.
[[nodiscard]] double TrueOf(const Interval& label, const Interval& truth, const double x)
{
    if (label == truth) {
        return x;
    }
    if (x <= label.GetLower()) {
        return truth.GetLower();
    }
    if (x >= label.GetUpper()) {
        return truth.GetUpper();
    }
    const double t = (x - label.GetLower()) / Width(label);
    return std::clamp(truth.GetLower() + t * Width(truth), truth.GetLower(), truth.GetUpper());
}

[[nodiscard]] Point2 TrueOf(const ChartPiece& piece, const Point2& p)
{
    return {TrueOf(piece.Labelled.Parameters.X, piece.True.X, p.X),
            TrueOf(piece.Labelled.Parameters.Y, piece.True.Y, p.Y)};
}

// --- The model of the height -----------------------------------------------

[[nodiscard]] Bounds HullOf(const Bounds& a, const Bounds& b)
{
    return {Hull(a[0], b[0]), Hull(a[1], b[1]), Hull(a[2], b[2])};
}

[[nodiscard]] Jet HullOf(const Jet& a, const Jet& b)
{
    return {HullOf(a.P, b.P),     HullOf(a.Pu, b.Pu),   HullOf(a.Pv, b.Pv),
            HullOf(a.Puu, b.Puu), HullOf(a.Puv, b.Puv), HullOf(a.Pvv, b.Pvv)};
}

/// A chart's pieces' nets, their jets over each piece, and the hull of
/// those.
struct ChartJets {
    std::vector<SecondNets> Nets;
    std::vector<Jet> Each;
    Jet Over;
};

[[nodiscard]] std::optional<ChartJets> JetsOf(const Chart& chart)
{
    ChartJets jets;
    for (const ChartPiece& piece : chart.Pieces) {
        jets.Nets.push_back(SecondNetsOf(piece.Labelled));
        const std::optional<Jet> over = JetOver(piece.Labelled, jets.Nets.back());
        if (!over) {
            return std::nullopt;
        }
        jets.Over = jets.Each.empty() ? *over : HullOf(jets.Over, *over);
        jets.Each.push_back(*over);
    }
    return jets;
}

/// A chart's pieces' jets at its corner, each from its own side.
[[nodiscard]] std::optional<std::vector<Jet>> CornerJets(const Chart& chart, const ChartJets& jets)
{
    std::vector<Jet> at;
    for (std::size_t i = 0; i < chart.Pieces.size(); ++i) {
        const std::optional<Jet> j = JetAt(chart.Pieces[i].Labelled, jets.Nets[i],
                                           Interval{chart.Corner.X}, Interval{chart.Corner.Y});
        if (!j) {
            return std::nullopt;
        }
        at.push_back(*j);
    }
    return at;
}

using Pair2 = std::array<Interval, 2>;

[[nodiscard]] Pair2 Applied(const Matrix2& m, const Pair2& x)
{
    return {m[0][0] * x[0] + m[0][1] * x[1], m[1][0] * x[0] + m[1][1] * x[1]};
}

/// What the proof knows of the height h of S above C around the corners
/// (fact 3): grad h(s) in Gradient + Hessian (s - c_S) over S's box, h(c_S)
/// in Height. With it the frame, Dq^-1 over each piece of C, S's point at
/// its corner, C's projected, and how far apart the corners are across the
/// plane.
struct SeamModel {
    NodeFrame Frame;
    std::vector<Matrix2> Inverses;
    Matrix2 Hessian;
    Pair2 Gradient;
    Interval Height;
    Bounds OnS;
    Pair2 CornerOfC;
    Pair2 Apart;
};

[[nodiscard]] Pair2 HullOf(const Pair2& a, const Pair2& b)
{
    return {Hull(a[0], b[0]), Hull(a[1], b[1])};
}

[[nodiscard]] Matrix2 HullOf(const Matrix2& a, const Matrix2& b)
{
    return {HullOf(a[0], b[0]), HullOf(a[1], b[1])};
}

/// A chart's jets over its pieces and at its corner.
struct JetsAround {
    const ChartJets& Over;
    const std::vector<Jet>& AtCorner;
};

/// How far apart the corners are: S's less C's, and that across the plane.
struct CornerGap {
    Bounds Apart;
    Pair2 Across;
};

/// What one piece of C adds to the model, with every piece of S: the
/// gradient's residual G and the Hessian H on each pair of pieces, each from
/// that pair's own jets, and h at S's corner if C's point below it is on
/// this piece. Enclosing a pair by its own jets matters: two pieces whose
/// parameters run at different speeds have different Hessians, and a product
/// of the hulls of their partials would be far wider than the hull of the
/// two.
struct ModelPart {
    Matrix2 Hessian;
    Pair2 Gradient;
    Interval Height;
    Matrix2 Inverse;
};

[[nodiscard]] std::optional<ModelPart> PartOfModel(const NodeFrame& f, const JetsAround& onS,
                                                   const Jet& overC, const Jet& atC,
                                                   const CornerGap& gap)
{
    const std::optional<Matrix2> inverse = InverseOf(ProjectedJacobian(f, overC.Pu, overC.Pv));
    const std::optional<Matrix2> atCorner = InverseOf(ProjectedJacobian(f, atC.Pu, atC.Pv));
    if (!inverse || !atCorner) {
        return std::nullopt;
    }
    const Pair2 bent =
        Applied(GraphHessian(CoNormal(f, *inverse, overC), overC, *inverse), gap.Across);
    // C's point below S's corner is Dq^-1 E from C's corner.
    const Pair2 below = Applied(*inverse, gap.Across);
    ModelPart part{{},
                   {},
                   DotOf(f.N, gap.Apart) -
                       (DotOf(f.N, overC.Pu) * below[0] + DotOf(f.N, overC.Pv) * below[1]),
                   *inverse};
    const Bounds nu = CoNormal(f, *atCorner, atC);
    for (std::size_t i = 0; i < onS.AtCorner.size(); ++i) {
        const Jet& s = onS.Over.Each[i];
        const Matrix2 dp = ProjectedJacobian(f, s.Pu, s.Pv);
        const Pair2 gradient{
            DotOf(nu, onS.AtCorner[i].Pu) - (dp[0][0] * bent[0] + dp[1][0] * bent[1]),
            DotOf(nu, onS.AtCorner[i].Pv) - (dp[0][1] * bent[0] + dp[1][1] * bent[1])};
        const Matrix2 hessian = HeightHessianOf(f, s, overC, *inverse);
        part.Gradient = i == 0 ? gradient : HullOf(part.Gradient, gradient);
        part.Hessian = i == 0 ? hessian : HullOf(part.Hessian, hessian);
    }
    return part;
}

/// Fact 3, from the charts' jets over their pieces and at their corners.
[[nodiscard]] std::optional<SeamModel> ModelOf(const NodeFrame& f, const JetsAround& onS,
                                               const JetsAround& onC)
{
    const Bounds& atS = onS.AtCorner.front().P;
    const Bounds& atC = onC.AtCorner.front().P;
    const Bounds apart = MinusOf(atS, atC);
    const CornerGap gap{apart, {DotOf(f.E1, apart), DotOf(f.E2, apart)}};
    std::optional<ModelPart> whole;
    std::vector<Matrix2> inverses;
    for (std::size_t k = 0; k < onC.AtCorner.size(); ++k) {
        const std::optional<ModelPart> part =
            PartOfModel(f, onS, onC.Over.Each[k], onC.AtCorner[k], gap);
        if (!part) {
            return std::nullopt;
        }
        inverses.push_back(part->Inverse);
        whole = !whole ? *part
                       : ModelPart{HullOf(whole->Hessian, part->Hessian),
                                   HullOf(whole->Gradient, part->Gradient),
                                   Hull(whole->Height, part->Height), part->Inverse};
    }
    return SeamModel{f,
                     std::move(inverses),
                     whole->Hessian,
                     whole->Gradient,
                     whole->Height,
                     atS,
                     {DotOf(f.E1, atC), DotOf(f.E2, atC)},
                     gap.Across};
}

/// h over offsets from S's corner, by the model: h(c_S) + G . x + x^T H x / 2.
[[nodiscard]] Interval ModelHeight(const SeamModel& m, const Interval& du, const Interval& dv)
{
    const Matrix2& h = m.Hessian;
    return m.Height + m.Gradient[0] * du + m.Gradient[1] * dv +
           Interval{0.5} *
               (h[0][0] * Squared(du) + (h[0][1] + h[1][0]) * du * dv + h[1][1] * Squared(dv));
}

/// h's derivative along u or v over offsets from S's corner, by the model.
[[nodiscard]] Interval ModelSlope(const SeamModel& m, const bool alongU, const Interval& du,
                                  const Interval& dv)
{
    const std::size_t r = alongU ? 0 : 1;
    return m.Gradient[r] + m.Hessian[r][0] * du + m.Hessian[r][1] * dv;
}

// --- Facts 1 and 2: C a graph under all of S --------------------------------

/// Whether a side of a piece is a side of its chart's box, not one it
/// shares with another piece: the edge running along a direction, the lower
/// or the upper across it.
[[nodiscard]] bool OuterSide(const Chart& chart, const ChartPiece& piece, const bool alongU,
                             const bool upper)
{
    const Interval across = StretchOf(piece.Labelled, !alongU);
    const Interval& whole = alongU ? chart.Box.Y : chart.Box.X;
    return upper ? across.GetUpper() == whole.GetUpper() : across.GetLower() == whole.GetLower();
}

/// S's projection keeps clear of the projection of the outer edges of C's
/// box.
[[nodiscard]] bool ProjectionClear(const Chart& s, const Chart& c, const NodeFrame& f)
{
    Interval x;
    Interval y;
    for (const ChartPiece& piece : s.Pieces) {
        x = Hull(x, Across(piece.Labelled.Poles.Poles, f.E1));
        y = Hull(y, Across(piece.Labelled.Poles.Poles, f.E2));
    }
    constexpr int depth = 6;
    for (const ChartPiece& piece : c.Pieces) {
        const Net& net = piece.Labelled.Poles;
        for (const bool alongU : {true, false}) {
            for (const bool upper : {false, true}) {
                const std::size_t line = upper ? LineCount(net, alongU) - 1 : 0;
                if (OuterSide(c, piece, alongU, upper) &&
                    !ProjectedApart(LineOf(net, alongU, line), f, x, y, depth)) {
                    return false;
                }
            }
        }
    }
    return true;
}

/// S's corner projects into C's box, by the Poincare-Miranda theorem on
/// y(c') = Y (q(c') - p), p the projection of S's corner and Y a
/// preconditioner: y = -Y E + Y M (c' - c_C) for some M in Dq's enclosure
/// (the mean value theorem along the segment from c_C, piece by piece), and
/// each component of y is proved negative on the box's lower face across it
/// and positive on its upper, so y, and with it q - p, has a zero in the box.
[[nodiscard]] bool CornerOverC(const SeamModel& m, const Jet& overC, const Chart& c)
{
    const Matrix2 dq = ProjectedJacobian(m.Frame, overC.Pu, overC.Pv);
    if (!InverseOf(dq)) {
        return false; // Fact 1: the hull of the pieces' Jacobians is regular.
    }
    const std::optional<DoubleMatrix<2>> y = Inverted<2>(
        {{{Midpoint(dq[0][0]), Midpoint(dq[0][1])}, {Midpoint(dq[1][0]), Midpoint(dq[1][1])}}});
    if (!y) {
        return false;
    }
    const Matrix2 exact{{{Interval{(*y)[0][0]}, Interval{(*y)[0][1]}},
                         {Interval{(*y)[1][0]}, Interval{(*y)[1][1]}}}};
    const Matrix2 slope = Times(exact, dq);
    const Pair2 offset = Applied(exact, m.Apart);
    const Pair2 box{c.Box.X - Interval{c.Corner.X}, c.Box.Y - Interval{c.Corner.Y}};
    for (const std::size_t r : {std::size_t{0}, std::size_t{1}}) {
        const std::size_t other = 1 - r;
        const Interval rest = slope[r][other] * box[other] - offset[r];
        const Interval below = slope[r][r] * Interval{box[r].GetLower()} + rest;
        const Interval above = slope[r][r] * Interval{box[r].GetUpper()} + rest;
        if (!(below.GetUpper() < 0.0 && above.GetLower() > 0.0)) {
            return false;
        }
    }
    return true;
}

// --- Facts 4 to 6 on the model ---------------------------------------------

/// The box D given up around S's corner, as offsets from it.
struct Inner {
    Interval Du;
    Interval Dv;
};

[[nodiscard]] double Magnitude(const Interval& x)
{
    return std::max(std::abs(x.GetLower()), std::abs(x.GetUpper()));
}

/// S within the tolerance of C all over D: the point of C below a point of
/// S is |h| from it along the frame's normal, a unit vector but for its
/// rounding, which the factor of two covers many times over.
[[nodiscard]] bool WithinTolerance(const SeamModel& m, const Inner& d, const double tolerance)
{
    return 2.0 * Magnitude(ModelHeight(m, d.Du, d.Dv)) <= tolerance;
}

/// The smaller distance from zero to an end of an interval holding it.
[[nodiscard]] double Narrower(const Interval& x)
{
    return std::min(-x.GetLower(), x.GetUpper());
}

/// Fact 4: x^T M J (G + H x) > 0 for every x outside D. In y = (x_u / a,
/// x_v / b), a and b the smaller half-widths of D, every such x has
/// |y| >= 1; the form is y^T P' y + w' . y with P' = L P L, w' = L M J G and
/// L = diag(a, b), at least lambda |y|^2 - |w'| |y| for lambda the smallest
/// eigenvalue of sym P', itself at least det / trace.
[[nodiscard]] bool LevelCurvesLeaveBeyond(const SeamModel& m, const Inner& d)
{
    const std::optional<Lyapunov> l = LyapunovOf(m.Hessian);
    if (!l) {
        return false;
    }
    const Interval a{Narrower(d.Du)};
    const Interval b{Narrower(d.Dv)};
    // J G = (-G_v, G_u).
    const Pair2 w = Applied(l->M, {-m.Gradient[1], m.Gradient[0]});
    const Interval p00 = a * a * l->P00;
    const Interval p11 = b * b * l->P11;
    const Interval s01 = a * b * l->S01;
    const Interval det = p00 * p11 - s01 * s01;
    const Interval residual =
        Sqrt(Squared(Interval{Magnitude(a * w[0])}) + Squared(Interval{Magnitude(b * w[1])}));
    return p00.GetLower() > 0.0 && det.GetLower() > 0.0 &&
           (det - (p00 + p11) * residual).GetLower() > 0.0;
}

/// A direction in S's parameters, and the Hessian's form along it.
struct Ray {
    double U;
    double V;
    Interval Form;
};

[[nodiscard]] Ray RayAlong(const SeamModel& m, const double u, const double v)
{
    const Matrix2& h = m.Hessian;
    const Interval eu{u};
    const Interval ev{v};
    return {u, v, h[0][0] * eu * eu + (h[0][1] + h[1][0]) * eu * ev + h[1][1] * ev * ev};
}

/// A lower bound of the parameter at which a ray from the corner leaves D.
[[nodiscard]] double LeavesAt(const Inner& d, const Ray& ray)
{
    const auto along = [](const Interval& side, const double e) {
        if (e == 0.0) {
            return Infinity;
        }
        return (e > 0.0 ? side.GetUpper() : -side.GetLower()) / std::abs(e);
    };
    return (1.0 - 0x1p-40) * std::min(along(d.Du, ray.U), along(d.Dv, ray.V));
}

/// Fact 5 on one ray: from where it leaves D, h's derivative along it,
/// e . G + t e^T H e, has the form's sign, and so has h there, h(c_S) +
/// t e . G + t^2 e^T H e / 2.
[[nodiscard]] bool RayKeepsSign(const SeamModel& m, const Inner& d, const Ray& ray)
{
    const double from = LeavesAt(d, ray);
    if (!ExcludesZero(ray.Form) || !std::isfinite(from)) {
        return false;
    }
    const Interval t{from};
    const Interval form{std::min(std::abs(ray.Form.GetLower()), std::abs(ray.Form.GetUpper()))};
    const Interval slope{
        Magnitude(Interval{ray.U} * m.Gradient[0] + Interval{ray.V} * m.Gradient[1])};
    const Interval height{Magnitude(m.Height)};
    return (form * t - slope).GetLower() > 0.0 &&
           (Interval{0.5} * form * t * t - slope * t - height).GetLower() > 0.0;
}

/// Fact 5: two lines through the corner, the Hessian's form positive along
/// one and negative along the other, and h keeping those signs on their
/// four rays outside D. The lines are the middle Hessian's eigenvectors,
/// rounded: any two with proved signs serve.
[[nodiscard]] bool RaysKeepSigns(const SeamModel& m, const Inner& d)
{
    const Matrix2& h = m.Hessian;
    const double off = Midpoint(h[0][1]) + Midpoint(h[1][0]);
    const double turn = 0.5 * std::atan2(off, Midpoint(h[0][0]) - Midpoint(h[1][1]));
    const double c = std::cos(turn);
    const double s = std::sin(turn);
    const std::array<Ray, 4> rays{RayAlong(m, c, s), RayAlong(m, -c, -s), RayAlong(m, -s, c),
                                  RayAlong(m, s, -c)};
    const bool opposite = (rays[0].Form.GetLower() > 0.0 && rays[2].Form.GetUpper() < 0.0) ||
                          (rays[0].Form.GetUpper() < 0.0 && rays[2].Form.GetLower() > 0.0);
    return opposite &&
           std::ranges::all_of(rays, [&](const Ray& ray) { return RayKeepsSign(m, d, ray); });
}

[[nodiscard]] int SignOfHeight(const Interval& x)
{
    if (x.GetLower() > 0.0) {
        return 1;
    }
    return x.GetUpper() < 0.0 ? -1 : 0;
}

/// One edge of D: its running direction and the other offset, fixed.
struct InnerEdge {
    bool AlongU;
    double Fixed;
};

[[nodiscard]] Interval HeightOnEdge(const SeamModel& m, const InnerEdge& e, const Interval& at)
{
    const Interval fixed{e.Fixed};
    return e.AlongU ? ModelHeight(m, at, fixed) : ModelHeight(m, fixed, at);
}

/// A walk along an edge of D, piece by piece: h's sign where it was last
/// proved, the sign of its derivative over the pieces passed since, on
/// which h is monotone (zero when there are none), and the zeros counted.
struct EdgeWalk {
    int Sign;
    int Slope = 0;
    std::size_t Zeros = 0;
};

/// h's sign is proved again, at `sign`: across a monotone stretch it has
/// one zero if the sign changed and none if not. False when the sign is not
/// proved, or changed with no stretch between, which cannot be.
[[nodiscard]] bool Arrive(EdgeWalk& walk, const int sign)
{
    if (sign == 0 || (walk.Slope == 0 && sign != walk.Sign)) {
        return false;
    }
    walk.Zeros += sign != walk.Sign ? 1U : 0U;
    walk.Sign = sign;
    walk.Slope = 0;
    return true;
}

/// One piece of the edge: h keeps off zero on it, so its sign is proved
/// there; or its derivative along the edge keeps a sign, the same as on the
/// pieces passed since h's sign was last proved. False otherwise.
[[nodiscard]] bool Step(EdgeWalk& walk, const SeamModel& m, const InnerEdge& e,
                        const Interval& piece)
{
    const int sign = SignOfHeight(HeightOnEdge(m, e, piece));
    if (sign != 0) {
        return Arrive(walk, sign);
    }
    const Interval fixed{e.Fixed};
    const int slope = SignOfHeight(e.AlongU ? ModelSlope(m, true, piece, fixed)
                                            : ModelSlope(m, false, fixed, piece));
    if (slope == 0 || (walk.Slope != 0 && slope != walk.Slope)) {
        return false;
    }
    walk.Slope = slope;
    return true;
}

/// h's zeros on an edge of D, counted on the model over 2^depth equal
/// pieces: between two places where its sign is proved, across pieces on
/// which it is monotone the same way, one zero if the signs differ and none
/// if they agree. Nothing when a piece is neither, or an end's sign is not
/// proved.
[[nodiscard]] std::optional<std::size_t> ZerosOnEdge(const SeamModel& m, const InnerEdge& e,
                                                     const Interval& range, const int depth)
{
    const int pieces = 1 << depth;
    const auto at = [&](const int k) {
        return range.GetLower() + Width(range) * static_cast<double>(k) / pieces;
    };
    EdgeWalk walk{SignOfHeight(HeightOnEdge(m, e, Interval{range.GetLower()}))};
    if (walk.Sign == 0) {
        return std::nullopt;
    }
    for (int k = 0; k < pieces; ++k) {
        const double from = k == 0 ? range.GetLower() : at(k);
        const double to = k + 1 == pieces ? range.GetUpper() : at(k + 1);
        if (!(from < to) || !Step(walk, m, e, Interval{from, to})) {
            return std::nullopt;
        }
    }
    if (!Arrive(walk, SignOfHeight(HeightOnEdge(m, e, Interval{range.GetUpper()})))) {
        return std::nullopt;
    }
    return walk.Zeros;
}

/// The finest an edge of D is cut while its zeros are counted: 2^10 pieces.
constexpr int InnerDepth = 10;

[[nodiscard]] std::optional<std::size_t> ZerosOnEdge(const SeamModel& m, const InnerEdge& e,
                                                     const Interval& range)
{
    for (int depth = 4; depth <= InnerDepth; depth += 2) {
        if (const std::optional<std::size_t> zeros = ZerosOnEdge(m, e, range, depth)) {
            return zeros;
        }
    }
    return std::nullopt;
}

/// Fact 6 on D: exactly four zeros of h on its boundary.
[[nodiscard]] bool FourZerosAround(const SeamModel& m, const Inner& d)
{
    std::size_t zeros = 0;
    for (const bool alongU : {true, false}) {
        const Interval& range = alongU ? d.Du : d.Dv;
        const Interval& across = alongU ? d.Dv : d.Du;
        for (const double fixed : {across.GetLower(), across.GetUpper()}) {
            const std::optional<std::size_t> found = ZerosOnEdge(m, {alongU, fixed}, range);
            if (!found) {
                return false;
            }
            zeros += *found;
        }
    }
    return zeros == 4;
}

/// Facts 4 to 6 on D, and S within the tolerance of C over it.
[[nodiscard]] bool InnerHolds(const SeamModel& m, const Inner& d, const double tolerance)
{
    return WithinTolerance(m, d, tolerance) && LevelCurvesLeaveBeyond(m, d) &&
           RaysKeepSigns(m, d) && FourZerosAround(m, d);
}

/// D: the smallest box around S's corner on which the facts hold, tried
/// from 2^-24 of the chart's narrower side there, four times wider each
/// time, to an eighth of it. The residuals of the corner, rounding, a
/// crease and a gap within the tolerance, need a box of some size around
/// them; the tolerance allows one only so large.
[[nodiscard]] std::optional<Inner> InnerOf(const Chart& s, const SeamModel& m,
                                           const double tolerance)
{
    const double roomU = std::min(s.Corner.X - s.Box.X.GetLower(), s.Box.X.GetUpper() - s.Corner.X);
    const double roomV = std::min(s.Corner.Y - s.Box.Y.GetLower(), s.Box.Y.GetUpper() - s.Corner.Y);
    if (!(roomU > 0.0 && roomV > 0.0)) {
        return std::nullopt;
    }
    for (int e = -24; e <= -3; e += 2) {
        const double du = std::ldexp(roomU, e);
        const double dv = std::ldexp(roomV, e);
        const Inner d{Interval{-Shorter * du, du}, Interval{-Shorter * dv, dv}};
        if (InnerHolds(m, d, tolerance)) {
            return d;
        }
    }
    return std::nullopt;
}

// --- Fact 6 on R: the four ends ---------------------------------------------

/// An end of a half-branch on R's boundary: its seed in the charts'
/// parameters, and the pieces of S and C it lies on.
struct SeamEnd {
    Seed At;
    std::size_t OnS;
    std::size_t OnC;
};

void AddEnd(std::vector<SeamEnd>& ends, const SeamEnd& end, const double join)
{
    const bool known = std::ranges::any_of(
        ends, [&](const SeamEnd& e) { return Distance(e.At.Position, end.At.Position) <= join; });
    if (!known) {
        ends.push_back(end);
    }
}

/// The two charts and what is proved of them.
struct SeamSystem {
    const Chart& S;
    const Chart& C;
    const SeamModel& Model;
    double Join;
};

/// One outer edge of a piece of S: the piece, the edge's running direction,
/// which of the two across it, and its fixed parameter.
struct OuterEdge {
    std::size_t OnS;
    bool AlongU;
    bool Upper;
    double Fixed;
};

/// The part of a piece of C that can hold a root from a stretch of C's
/// parameters; no part, or nothing settled when it cannot be cut out.
[[nodiscard]] Result<std::optional<Patch>> PartOf(const ChartPiece& piece, const Interval& u,
                                                  const Interval& v)
{
    const Box2& box = piece.Labelled.Parameters;
    const Interval ku = Intersection(u, box.X);
    const Interval kv = Intersection(v, box.Y);
    // A stretch meeting the piece only on a side: a root there is on the
    // piece beyond that side as well, whose own stretch, widened, holds
    // more than a side of it, or on an outer edge of C, where fact 2 leaves
    // none.
    if (ku.IsEmpty() || kv.IsEmpty() || !(ku.GetLower() < ku.GetUpper()) ||
        !(kv.GetLower() < kv.GetUpper())) {
        return std::optional<Patch>{};
    }
    STEREON_TRY(std::optional<Patch> part, Clipped(piece.Labelled, true, ku));
    if (!part) {
        return std::optional<Patch>{};
    }
    return Clipped(std::move(*part), false, kv);
}

/// The stretch of a piece of C's parameters that can hold a root on a line
/// of S's poles, by the mean value theorem on the piece, from the corner in
/// its closure: q(c') - q(c_C) = M (c' - c_C), M in Dq's enclosure over the
/// piece, and q(c') is the line's projection.
[[nodiscard]] Pair2 ReachOnC(const SeamSystem& sys, const std::vector<HPoint>& line,
                             const std::size_t onC)
{
    const SeamModel& m = sys.Model;
    const Pair2 from{Across(line, m.Frame.E1) - m.CornerOfC[0],
                     Across(line, m.Frame.E2) - m.CornerOfC[1]};
    const Pair2 by = Applied(m.Inverses[onC], from);
    // Widened by a 2^-16th of C's box: an edge of S that projects onto a
    // line would otherwise leave a part of C a few roundings wide, on which
    // no derivative is enclosed usefully.
    return {Inflated(Interval{sys.C.Corner.X} + by[0], 0x1p-16 * Width(sys.C.Box.X)),
            Inflated(Interval{sys.C.Corner.Y} + by[1], 0x1p-16 * Width(sys.C.Box.Y))};
}

/// A piece of an outer edge, cut out of its piece of S.
struct EdgeCarrier {
    const Patch& Carrier;
    const OuterEdge& Of;
};

/// The roots of a piece of an outer edge against every piece of C, by
/// Krawczyk's operator on the part of each that can hold one. Unsettled
/// when any is.
[[nodiscard]] Isolation EndsOnCarrier(const SeamSystem& sys, const EdgeCarrier& e,
                                      std::vector<SeamEnd>& ends)
{
    const std::vector<HPoint> line = EdgeLine(e.Carrier, e.Of.AlongU, e.Of.Upper);
    const Edge edge{e.Carrier, e.Of.AlongU, e.Of.Fixed};
    const Hulls own = LocalHullsOf(e.Carrier);
    for (std::size_t k = 0; k < sys.C.Pieces.size(); ++k) {
        const Pair2 reach = ReachOnC(sys, line, k);
        const Result<std::optional<Patch>> part = PartOf(sys.C.Pieces[k], reach[0], reach[1]);
        if (!part) {
            return Isolation::Unsettled;
        }
        if (!*part) {
            continue;
        }
        const Isolated<3> root = RootOnEdge(edge, line, own, **part, LocalHullsOf(**part));
        if (root.Outcome == Isolation::Unsettled) {
            return Isolation::Unsettled;
        }
        if (root.Outcome == Isolation::Unique) {
            AddEnd(ends, {SeedOfRoot(edge, root, **part, true), e.Of.OnS, k}, sys.Join);
        }
    }
    return Isolation::None;
}

/// The roots of a piece of an outer edge: none when the model keeps the
/// height off zero there, else EndsOnCarrier's verdict.
[[nodiscard]] Isolation EndsOnPiece(const SeamSystem& sys, const OuterEdge& e,
                                    const Interval& piece, std::vector<SeamEnd>& ends)
{
    const Point2& corner = sys.S.Corner;
    const Interval du = (e.AlongU ? piece : Interval{e.Fixed}) - Interval{corner.X};
    const Interval dv = (e.AlongU ? Interval{e.Fixed} : piece) - Interval{corner.Y};
    if (ExcludesZero(ModelHeight(sys.Model, du, dv))) {
        return Isolation::None;
    }
    const Result<std::optional<Patch>> carrier =
        Clipped(sys.S.Pieces[e.OnS].Labelled, e.AlongU, piece);
    if (!carrier || !*carrier) {
        return Isolation::Unsettled; // A failed cut cannot prove an edge root absent.
    }
    return EndsOnCarrier(sys, {**carrier, e}, ends);
}

/// The roots of one outer edge, a piece of it halved while it is not
/// settled; false when one is still not settled at the deepest halving.
[[nodiscard]] bool EndsOnEdge(const SeamSystem& sys, const OuterEdge& e, std::vector<SeamEnd>& ends)
{
    const Patch& carrier = sys.S.Pieces[e.OnS].Labelled;
    std::vector<std::pair<Interval, int>> work{{StretchOf(carrier, e.AlongU), 0}};
    while (!work.empty()) {
        const auto [piece, depth] = work.back();
        work.pop_back();
        if (EndsOnPiece(sys, e, piece, ends) != Isolation::Unsettled) {
            continue;
        }
        if (depth >= EdgeDepth) {
            return false;
        }
        const double middle = Midpoint(piece);
        work.emplace_back(Interval{piece.GetLower(), middle}, depth + 1);
        work.emplace_back(Interval{middle, piece.GetUpper()}, depth + 1);
    }
    return true;
}

/// Fact 6 on R: the roots on the outer edges of S's pieces, exactly four,
/// or nothing.
[[nodiscard]] std::optional<std::vector<SeamEnd>> FourEnds(const SeamSystem& sys)
{
    std::vector<SeamEnd> ends;
    for (std::size_t i = 0; i < sys.S.Pieces.size(); ++i) {
        const ChartPiece& piece = sys.S.Pieces[i];
        for (const bool alongU : {true, false}) {
            const Interval across = StretchOf(piece.Labelled, !alongU);
            for (const bool upper : {false, true}) {
                const double fixed = upper ? across.GetUpper() : across.GetLower();
                if (OuterSide(sys.S, piece, alongU, upper) &&
                    !EndsOnEdge(sys, {i, alongU, upper, fixed}, ends)) {
                    return std::nullopt;
                }
            }
        }
    }
    if (ends.size() != 4) {
        return std::nullopt;
    }
    return ends;
}

// --- The proof and its record ------------------------------------------------

/// A node settled on two charts: its point and the four ends, in the
/// coordinates the charts were taken in.
struct SeamNode {
    Point3 Position;
    std::vector<SeamEnd> Ends;
};

/// A frame about C at its corner with e1 along C's partial there, e2
/// across it in the tangent plane: the projections' Jacobians are then near
/// diagonal where the parameter lines are near orthogonal, and their
/// enclosed inverses far narrower than in a frame turned about the normal
/// (a box of a turned Jacobian wraps it, and the determinant's enclosure
/// loses what its entries do, twice over). Any basis serves the proof.
[[nodiscard]] std::optional<NodeFrame> FrameAlong(const Jet& atC)
{
    const Vec3 du = MiddleVec(atC.Pu);
    const Vec3 normal = Cross(du, MiddleVec(atC.Pv));
    const double length = Length(normal);
    if (!(length > 0.0) || !std::isfinite(length) || !(Length(du) > 0.0)) {
        return std::nullopt;
    }
    const Vec3 n = normal / length;
    const Vec3 e1 = Normalized(du);
    return NodeFrame{Exactly(e1), Exactly(Cross(n, e1)), Exactly(n)};
}

/// Facts 1 to 3 on two charts: the model, C a graph under all of S.
[[nodiscard]] std::optional<SeamModel> ProvedModel(const Chart& s, const ChartJets& onS,
                                                   const Chart& c, const ChartJets& onC)
{
    const std::optional<std::vector<Jet>> atS = CornerJets(s, onS);
    const std::optional<std::vector<Jet>> atC = CornerJets(c, onC);
    const std::optional<NodeFrame> frame = atC ? FrameAlong(atC->front()) : std::nullopt;
    if (!frame || !atS) {
        return std::nullopt;
    }
    std::optional<SeamModel> model = ModelOf(*frame, {onS, *atS}, {onC, *atC});
    if (!model || !ProjectionClear(s, c, *frame) || !CornerOverC(*model, onC.Over, c)) {
        return std::nullopt;
    }
    return model;
}

/// Facts 1 to 6 on a chart of S and one of C, with the tolerance D and the
/// corners' gap must be within and the join of seeds.
[[nodiscard]] std::optional<SeamNode> ProvedSeamNode(const Chart& s, const Chart& c,
                                                     const SearchSettings& settings)
{
    const std::optional<ChartJets> onS = JetsOf(s);
    const std::optional<ChartJets> onC = JetsOf(c);
    if (!onS || !onC) {
        return std::nullopt;
    }
    const std::optional<SeamModel> model = ProvedModel(s, *onS, c, *onC);
    if (!model || !InnerOf(s, *model, settings.Tolerance)) {
        return std::nullopt;
    }
    std::optional<std::vector<SeamEnd>> ends = FourEnds({s, c, *model, settings.Join});
    if (!ends) {
        return std::nullopt;
    }
    return SeamNode{Middle(model->OnS), std::move(*ends)};
}

/// The node's four arcs, each on the pieces its end lies on, with the
/// node's own parameters there; the node; and the region it settles. The
/// charts were taken about `centre`, which the points are moved back by.
void RecordSeamNode(const SeamNode& node, const Chart& s, const Chart& c, const bool smallIsA,
                    const Point3& centre, Findings& found)
{
    const auto ordered = [&](const Seed& seed) { return smallIsA ? seed : Swapped(seed); };
    const Vec3 back = ToVec(centre);
    const Point3 at = node.Position + back;
    for (const SeamEnd& e : node.Ends) {
        const ChartPiece& onS = s.Pieces[e.OnS];
        const ChartPiece& onC = c.Pieces[e.OnC];
        const Seed from = ordered({at, TrueOf(onS, s.Corner), TrueOf(onC, c.Corner)});
        const Seed to =
            ordered({e.At.Position + back, TrueOf(onS, e.At.OnA), TrueOf(onC, e.At.OnB)});
        const Vec3 d = ToVec(to.Position) - ToVec(from.Position);
        found.Arcs.push_back({from, to, d / Length(d), smallIsA ? onS.True : onC.True,
                              smallIsA ? onC.True : onS.True});
    }
    found.Nodes.push_back(at);
    NodeRegion region;
    for (const ChartPiece& piece : s.Pieces) {
        (smallIsA ? region.OnA : region.OnB).push_back(piece.True);
    }
    for (const ChartPiece& piece : c.Pieces) {
        (smallIsA ? region.OnB : region.OnA).push_back(piece.True);
    }
    found.Regions.push_back(std::move(region));
}

// --- The candidate ------------------------------------------------------------

/// The critical point's system of the proof inside a cell and its Jacobian,
/// in doubles at a point: e . (S - C), S_u . N_C, S_v . N_C in (u_S, v_S,
/// u_C, v_C).
struct GuessSystem {
    std::array<double, 4> F;
    DoubleMatrix<4> J;
};

[[nodiscard]] GuessSystem GuessSystemOf(const Jet& s, const Jet& c, const Vec3& e1, const Vec3& e2)
{
    const auto v = [](const Bounds& b) { return MiddleVec(b); };
    const Vec3 r = v(s.P) - v(c.P);
    const Vec3 n = Cross(v(c.Pu), v(c.Pv));
    const Vec3 nu = Cross(v(c.Puu), v(c.Pv)) + Cross(v(c.Pu), v(c.Puv));
    const Vec3 nv = Cross(v(c.Puv), v(c.Pv)) + Cross(v(c.Pu), v(c.Pvv));
    return {{Dot(e1, r), Dot(e2, r), Dot(v(s.Pu), n), Dot(v(s.Pv), n)},
            {{{Dot(e1, v(s.Pu)), Dot(e1, v(s.Pv)), -Dot(e1, v(c.Pu)), -Dot(e1, v(c.Pv))},
              {Dot(e2, v(s.Pu)), Dot(e2, v(s.Pv)), -Dot(e2, v(c.Pu)), -Dot(e2, v(c.Pv))},
              {Dot(v(s.Puu), n), Dot(v(s.Puv), n), Dot(v(s.Pu), nu), Dot(v(s.Pu), nv)},
              {Dot(v(s.Puv), n), Dot(v(s.Pvv), n), Dot(v(s.Pv), nu), Dot(v(s.Pv), nv)}}}};
}

/// Newton's step from x, or nothing where the system cannot be formed.
[[nodiscard]] std::optional<std::array<double, 4>>
GuessStep(const PatchPair& pair, const SecondNets& a, const SecondNets& b,
          const std::array<Vec3, 2>& e, const std::array<double, 4>& x)
{
    const std::optional<Jet> s = JetAt(pair.A, a, Interval{x[0]}, Interval{x[1]});
    const std::optional<Jet> c = JetAt(pair.B, b, Interval{x[2]}, Interval{x[3]});
    if (!s || !c) {
        return std::nullopt;
    }
    const GuessSystem system = GuessSystemOf(*s, *c, e[0], e[1]);
    const std::optional<DoubleMatrix<4>> inverse = Inverted<4>(system.J);
    if (!inverse) {
        return std::nullopt;
    }
    std::array<double, 4> step{};
    for (std::size_t r = 0; r < 4; ++r) {
        for (std::size_t k = 0; k < 4; ++k) {
            step[r] -= (*inverse)[r][k] * system.F[k];
        }
    }
    return step;
}

/// How far apart the surfaces are at a guess, in doubles; infinity where a
/// point cannot be formed.
[[nodiscard]] double GuessGap(const PatchPair& pair, const std::array<double, 4>& x)
{
    const Bounds onA = XyzOver(pair.A, Interval{x[0]}, Interval{x[1]});
    const Bounds onB = XyzOver(pair.B, Interval{x[2]}, Interval{x[3]});
    return AllBounded(onA) && AllBounded(onB) ? Distance(Middle(onA), Middle(onB)) : Infinity;
}

/// Newton's iterations for a candidate.
constexpr int GuessLimit = 12;

/// Whether every component of a step is within a fraction of its box.
[[nodiscard]] bool StepWithin(const std::array<double, 4>& step, const std::array<Interval, 4>& box,
                              const double fraction)
{
    for (std::size_t r = 0; r < 4; ++r) {
        if (!(std::abs(step[r]) <= fraction * Width(box[r]))) {
            return false;
        }
    }
    return true;
}

/// Where a pair may hold a point of common tangency, in its boxes or just
/// outside them: Newton's method in doubles from the pair's middle, each
/// iterate kept in the boxes, until a step is a 2^-30th of them. A point of
/// tangency just beyond a side leaves the iterate on that side, stepping
/// outwards by its distance each time; it is accepted there, the last step
/// within a 2^-12th of the boxes, as it lies well within a chart's reach
/// across the side (ReachOf). The surfaces must be within the tolerance of
/// each other at it. A guess: it proves nothing and is trusted with nothing.
[[nodiscard]] std::optional<std::array<double, 4>> NodeGuess(const PatchPair& pair,
                                                             const double tolerance)
{
    const SecondNets a = SecondNetsOf(pair.A);
    const SecondNets b = SecondNetsOf(pair.B);
    const std::array<Interval, 4> box{pair.A.Parameters.X, pair.A.Parameters.Y, pair.B.Parameters.X,
                                      pair.B.Parameters.Y};
    std::array<double, 4> x = MidpointsOf<4>(box);
    const std::optional<NodeFrame> frame = FrameAbout(NormalAt(pair.B, b, x[2], x[3]));
    if (!frame) {
        return std::nullopt;
    }
    const std::array<Vec3, 2> e{MiddleVec(frame->E1), MiddleVec(frame->E2)};
    bool near = false;
    for (int k = 0; k < GuessLimit; ++k) {
        const std::optional<std::array<double, 4>> step = GuessStep(pair, a, b, e, x);
        if (!step) {
            return std::nullopt;
        }
        for (std::size_t r = 0; r < 4; ++r) {
            x[r] = std::clamp(x[r] + (*step)[r], box[r].GetLower(), box[r].GetUpper());
        }
        near = StepWithin(*step, box, 0x1p-12);
        if (StepWithin(*step, box, 0x1p-30)) {
            break;
        }
    }
    return near && GuessGap(pair, x) <= tolerance ? std::optional{x} : std::nullopt;
}

/// How near a side of the patch that found it a candidate counts as on it:
/// a 2^-10th of the patch along that direction. D can take an eighth of the
/// chart, so a node so near a side still has room in it.
[[nodiscard]] Reach ReachOf(const Patch& p, const bool alongU, const double factor)
{
    return {factor, 0x1p-10 * Width(StretchOf(p, alongU))};
}

/// The pair as (small, cover), a candidate's parameters on each, and which
/// of the search's surfaces the small one is.
struct SeamCandidate {
    const PatchPair& Oriented;
    Point2 OnS;
    Point2 OnC;
    bool SmallIsA;
};

/// Whether a chart shares area with any of the boxes.
[[nodiscard]] bool SharesArea(const Chart& chart, const std::vector<Box2>& boxes)
{
    return std::ranges::any_of(chart.Pieces, [&](const ChartPiece& piece) {
        return std::ranges::any_of(
            boxes, [&](const Box2& box) { return HasAreaOverlap(piece.True, box); });
    });
}

/// Whether nothing found so far lies where the charts would settle: no arc
/// of a pair of patches, and no settled region, shares area with them on
/// both surfaces. A node may be settled only where the search has left no
/// branch, or the branch would be reported twice.
[[nodiscard]] bool NothingFoundOn(const Chart& s, const Chart& c, const bool smallIsA,
                                  const Findings& found)
{
    const Chart& onA = smallIsA ? s : c;
    const Chart& onB = smallIsA ? c : s;
    const bool arc = std::ranges::any_of(found.Arcs, [&](const Arc& a) {
        return SharesArea(onA, {a.BoxA}) && SharesArea(onB, {a.BoxB});
    });
    const bool region = std::ranges::any_of(found.Regions, [&](const NodeRegion& r) {
        return SharesArea(onA, r.OnA) && SharesArea(onB, r.OnB);
    });
    return !arc && !region;
}

/// One attempt: the charts laid with C's taking 2^-(j+1) of its Bezier
/// patch and S's 2^-(k+1), the node proved on them and recorded. A node with one
/// piece on each chart is inside a cell of both, and is SettleNode's.
[[nodiscard]] bool SettleOnCharts(const SeamCandidate& at, const int j, const int k,
                                  const SearchSettings& s, Findings& found)
{
    const Atlas& atlasS = at.SmallIsA ? *s.OnA : *s.OnB;
    const Atlas& atlasC = at.SmallIsA ? *s.OnB : *s.OnA;
    const Patch& small = at.Oriented.A;
    const Patch& cover = at.Oriented.B;
    const double fs = std::ldexp(0.5, -k);
    const double fc = std::ldexp(0.5, -j);
    const std::optional<Chart> onS =
        ChartAbout(atlasS, small, at.OnS, ReachOf(small, true, fs), ReachOf(small, false, fs));
    const std::optional<Chart> onC =
        ChartAbout(atlasC, cover, at.OnC, ReachOf(cover, true, fc), ReachOf(cover, false, fc));
    if (!onS || !onC || onS->Pieces.size() * onC->Pieces.size() == 1 ||
        !NothingFoundOn(*onS, *onC, at.SmallIsA, found)) {
        return false;
    }
    const ChartPiece& first = onS->Pieces.front();
    const Bounds corner = XyzOver(first.Labelled, Interval{onS->Corner.X}, Interval{onS->Corner.Y});
    if (!AllBounded(corner)) {
        return false;
    }
    const Point3 centre = Middle(corner);
    const Chart aboutS = About(*onS, centre);
    const Chart aboutC = About(*onC, centre);
    const std::optional<SeamNode> node = ProvedSeamNode(aboutS, aboutC, s);
    if (!node) {
        return false;
    }
    RecordSeamNode(*node, aboutS, aboutC, at.SmallIsA, centre, found);
    return true;
}

/// How often the charts are narrowed around a candidate: C's, and within
/// each of C's, S's.
constexpr int ChartNarrowings = 6;
constexpr int SmallChartNarrowings = 6;

/// Whether a side of the patch is a knot line or a closing seam with a
/// patch beyond it: where a node may lie that no single patch holds.
[[nodiscard]] bool OnInterface(const Atlas& atlas, const Patch& p)
{
    const Patch* own = Holding(atlas, p.Parameters);
    if (own == nullptr) {
        return false;
    }
    for (const bool alongU : {true, false}) {
        const Interval mine = StretchOf(p, alongU);
        const Interval whole = StretchOf(*own, alongU);
        for (const bool upper : {false, true}) {
            const bool shared =
                upper ? mine.GetUpper() == whole.GetUpper() : mine.GetLower() == whole.GetLower();
            if (shared && BeyondSide(atlas, *own, alongU, upper)) {
                return true;
            }
        }
    }
    return false;
}

/// The charts around a candidate narrowed until the node is proved on
/// them.
[[nodiscard]] bool SettleAround(const SeamCandidate& at, const SearchSettings& s, Findings& found)
{
    for (int j = 0; j < ChartNarrowings; ++j) {
        for (int k = j; k < j + SmallChartNarrowings; ++k) {
            if (SettleOnCharts(at, j, k, s, found)) {
                return true;
            }
        }
    }
    return false;
}

[[nodiscard]] bool MayHoldSeamNode(const PatchPair& pair, const SearchSettings& s)
{
    return s.OnA != nullptr && s.OnB != nullptr &&
           (OnInterface(*s.OnA, pair.A) || OnInterface(*s.OnB, pair.B));
}

/// A pair that may hold a node on a knot line or seam of either surface:
/// a candidate, then the charts around it narrowed until the node is
/// proved, either surface taken as the small one. The node's region is then
/// recorded, and the pair goes back on the stack to be cut along it.
[[nodiscard]] bool SettleSeamNode(const PatchPair& pair, const SearchSettings& s, Findings& found,
                                  std::vector<PatchPair>& stack)
{
    const std::optional<std::array<double, 4>> m =
        MayHoldSeamNode(pair, s) ? NodeGuess(pair, s.Tolerance) : std::nullopt;
    if (!m) {
        return false;
    }
    const Point2 onA{(*m)[0], (*m)[1]};
    const Point2 onB{(*m)[2], (*m)[3]};
    const PatchPair swapped{pair.B, pair.A};
    if (SettleAround({pair, onA, onB, true}, s, found) ||
        SettleAround({swapped, onB, onA, false}, s, found)) {
        stack.push_back(pair);
        return true;
    }
    return false;
}

// --- Cutting the search along settled regions ------------------------------

/// What became of a patch against a box it shares area with: cut in two
/// along a side of the box that crosses it, inside the box, or not cut
/// where a cut cannot be enclosed.
enum class Cut : std::uint8_t { Halved, Inside, Failed };

struct CutPatch {
    Cut Outcome;
    std::optional<std::pair<Patch, Patch>> Halves;
};

[[nodiscard]] CutPatch CutAlong(const Patch& p, const Box2& box)
{
    for (const bool alongU : {true, false}) {
        const Interval stretch = StretchOf(p, alongU);
        const Interval& side = alongU ? box.X : box.Y;
        for (const double at : {side.GetLower(), side.GetUpper()}) {
            if (stretch.GetLower() < at && at < stretch.GetUpper()) {
                auto halves = Split(p, alongU, at);
                return halves ? CutPatch{Cut::Halved, std::move(halves)}
                              : CutPatch{Cut::Failed, std::nullopt};
            }
        }
    }
    return {Cut::Inside, std::nullopt};
}

/// What a settled region makes of a pair.
enum class Carving : std::uint8_t { Untouched, Consumed, Failed };

/// A pair against one box of A's and one of B's of a settled region it
/// shares area with on both: cut along a side of either box, its halves
/// back on the stack, or, inside both, dropped as settled.
[[nodiscard]] Carving CarvedBy(const PatchPair& pair, const Box2& onA, const Box2& onB,
                               std::vector<PatchPair>& stack)
{
    CutPatch a = CutAlong(pair.A, onA);
    if (a.Outcome == Cut::Halved) {
        stack.push_back({std::move(a.Halves->first), pair.B});
        stack.push_back({std::move(a.Halves->second), pair.B});
        return Carving::Consumed;
    }
    CutPatch b = a.Outcome == Cut::Inside ? CutAlong(pair.B, onB) : CutPatch{Cut::Failed, {}};
    if (b.Outcome == Cut::Halved) {
        stack.push_back({pair.A, std::move(b.Halves->first)});
        stack.push_back({pair.A, std::move(b.Halves->second)});
        return Carving::Consumed;
    }
    return b.Outcome == Cut::Inside ? Carving::Consumed : Carving::Failed;
}

[[nodiscard]] Carving CarvedBy(const PatchPair& pair, const NodeRegion& region,
                               std::vector<PatchPair>& stack)
{
    for (const Box2& onA : region.OnA) {
        if (!HasAreaOverlap(pair.A.Parameters, onA)) {
            continue;
        }
        for (const Box2& onB : region.OnB) {
            if (HasAreaOverlap(pair.B.Parameters, onB)) {
                return CarvedBy(pair, onA, onB, stack);
            }
        }
    }
    return Carving::Untouched;
}

/// The pair against every settled region: what lies in one is settled, so
/// the pair is cut along the region's sides until each part is inside it,
/// and dropped, or shares no area with it, and searched.
[[nodiscard]] Carving Carved(const PatchPair& pair, const std::vector<NodeRegion>& regions,
                             std::vector<PatchPair>& stack)
{
    for (const NodeRegion& region : regions) {
        const Carving carving = CarvedBy(pair, region, stack);
        if (carving != Carving::Untouched) {
            return carving;
        }
    }
    return Carving::Untouched;
}

} // namespace
} // namespace Stereon::Detail
