// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <cstddef>
#include <vector>

export module Stereon.Boolean:Winding;

import Stereon.Build;

/// The bounded generalised winding number (the Boolean prototype's M7;
/// ADR-0013, the decisions on acceptance): the fallback that decides a
/// region's seed when every ray meets the other body ambiguously. The body's
/// faces are triangulated in their parameter spaces, every face taking the
/// same points along an edge it shares, so the triangles close up into a
/// polyhedron, and the winding number of that polyhedron at the point is
/// the sum of the triangles' solid angles over 4 pi. It equals the body's
/// own winding number wherever the point is farther from every triangle
/// than the triangles stray from the face: a surface that agrees with the
/// face inside and runs along the chords at its boundary deforms onto the
/// triangles within that distance without crossing the point, so the two
/// numbers differ by a whole number that must be zero. Where a cell of a
/// face's parameter space is not separated from the point that way it is
/// split and its loops clipped to the halves, down to the round's deviation
/// target; a round that leaves cells unseparated is followed by one with
/// the target quartered, up to a budget. The result carries its bound, and a
/// side is read from it only when the number lies further from one half
/// than the bound (ADR-0029, section 4), never by guessing.

export namespace Stereon {

/// The generalised winding number of a body's boundary at a point, as the
/// triangulated faces give it, with what certifies it.
struct WindingNumber {
    /// The polyhedron's winding number: 1 inside a closed body, 0 outside.
    double Value = 0.0;
    /// The body's own number lies within Bound of Value: the cells of faces
    /// that could not be separated from the point, each worth up to one,
    /// plus the rounding of the sum.
    double Bound = 0.0;
    /// Rounds of refinement the result took, each quartering the deviation.
    int Refinements = 0;
    /// Triangles summed in the last round.
    std::size_t Triangles = 0;

    /// True when no error within the bound could change the side.
    [[nodiscard]] bool IsDecided() const noexcept;
    /// The side of one half the value lies on; meaningful once decided.
    [[nodiscard]] bool IsInside() const noexcept;

    bool operator==(const WindingNumber&) const = default;
};

/// The winding number of `body`'s boundary at `point`, refined up to
/// `refinements` rounds (quartering the deviation target each round, from a
/// sixty-fourth of the body's extent down to the Context's resolution) or
/// until decided. An undecided result is returned with its bound, for the
/// caller to report. NotSupported for a face outside the prototype's planes
/// and cylinders, or a pcurve outside its lines, circles, ellipses and
/// non-rational B-splines; InvalidResult when a decided number is neither 0
/// nor 1, which means the faces do not close up or one is oriented against
/// its loops; InvalidInput for a negative budget or a body without vertices.
[[nodiscard]] Result<WindingNumber> WindingNumberOf(const Shape& body, const Point3& point,
                                                    int refinements, const Context& ctx);

} // namespace Stereon

export namespace Stereon::Detail {

/// Where WindingNumberOf samples an edge for a deviation target, from one
/// end of its range to the other, and how far the chords between the
/// samples may lie from the edge on any face using it: the target, or more
/// where the sample count is capped. Every knot of a B-spline pcurve inside
/// the range is a sample. Exposed for the tests.
struct EdgeSampling {
    std::vector<double> Parameters;
    double Deviation = 0.0;
};

[[nodiscard]] Result<EdgeSampling> WindingSamplingOf(const Shape& shape, EdgeId edge,
                                                     double target);

} // namespace Stereon::Detail
