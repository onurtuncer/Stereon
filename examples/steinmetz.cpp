// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The Steinmetz solids (docs/examples.rst): two cylinders of one radius
// whose axes cross at right angles, united, intersected and cut either
// way, and a third across both. The common of two is the Steinmetz solid,
// the common of three the tricylinder. Each result is counted, checked
// closed by the validity checker, and asked by its own winding number
// whether it holds points known to be in some of the cylinders, in all or
// in none.

#include <array>
#include <cstdio>
#include <print>

import Stereon.Boolean;

using namespace Stereon;

namespace {

/// A point, and whether the cylinder about z, the one about x and the one
/// about y hold it.
struct Sample {
    Point3 At;
    bool InZ;
    bool InX;
    bool InY;
};

/// Radius 3: the cylinder about z is x^2 + y^2 <= 9 for |z| <= 5, the one
/// about x is y^2 + z^2 <= 9 for |x| <= 6, the one about y is x^2 + z^2 <= 9
/// for |y| <= 6. One sample for each way of being in some and out of others.
constexpr std::array samples{
    Sample{{0.5, 0.5, 0.5}, true, true, true},   Sample{{2.2, 0.5, 2.2}, true, true, false},
    Sample{{0.5, 0.5, 4.0}, true, false, false}, Sample{{0.5, 2.5, 2.5}, true, false, true},
    Sample{{5.0, 0.5, 0.5}, false, true, false}, Sample{{2.5, 2.5, 0.5}, false, true, true},
    Sample{{0.5, 3.5, 0.5}, false, false, true}, Sample{{4.0, 0.5, 4.0}, false, false, false},
};

/// Closed and consistently oriented: at the Solid level the checker finds
/// nothing but the self-intersection check, which that level does not make.
bool Closed(const Shape& shape, const Context& ctx)
{
    const ValidityReport report = ValidateTopology(shape, ValidityLevel::Solid, ctx);
    for (const ValidityFinding& finding : report.GetFindings()) {
        if (finding.Rule != TopologyRule::SelfIntersectionUnchecked) {
            return false;
        }
    }
    return true;
}

/// Whether a sample belongs to a result.
using Membership = bool (*)(const Sample&);

/// A Boolean's result, counted: whether it is closed and holds each sample
/// exactly when `expected` says it should.
bool Checked(const char* name, const Result<OperationResult>& result, const Membership expected,
             const Context& ctx)
{
    if (!result) {
        std::println(stderr, "{}: {}", name, result.error().GetMessage());
        return false;
    }
    const Shape& shape = result->Shape;
    bool members = true;
    for (const Sample& sample : samples) {
        const Result<WindingNumber> winding = WindingNumberOf(shape, sample.At, 8, ctx);
        members =
            members && winding && winding->IsDecided() && winding->IsInside() == expected(sample);
    }
    const bool closed = Closed(shape, ctx);
    std::println(
        "{:<12} solids {}  faces {:>2}  edges {:>2}  vertices {:>2}  closed {}  members {}", name,
        shape.Count<EntityKind::Solid>(), shape.Count<EntityKind::Face>(),
        shape.Count<EntityKind::Edge>(), shape.Count<EntityKind::Vertex>(), closed, members);
    return closed && members;
}

/// What each result should hold.
bool InEither(const Sample& s)
{
    return s.InZ || s.InX;
}

bool InBoth(const Sample& s)
{
    return s.InZ && s.InX;
}

bool InZAlone(const Sample& s)
{
    return s.InZ && !s.InX;
}

bool InXAlone(const Sample& s)
{
    return s.InX && !s.InZ;
}

bool InAllThree(const Sample& s)
{
    return s.InZ && s.InX && s.InY;
}

/// One Boolean of two of the cylinders.
struct Case {
    const char* Name;
    const Shape* A;
    const Shape* B;
    BooleanOperation Operation;
    Membership Expected;
};

} // namespace

int main()
{
    const Context ctx = Context::Default();

    // Three cylinders of radius 3 whose axes cross at the origin: about z
    // from z = -5 to 5, about x and about y from -6 to 6. A cylinder's axis
    // is the cross product of its two frame vectors.
    const Vec3 ex{1.0, 0.0, 0.0};
    const Vec3 ey{0.0, 1.0, 0.0};
    const Vec3 ez{0.0, 0.0, 1.0};
    const Result<OperationResult> aboutZ =
        MakeCylinderSolid(Point3{0.0, 0.0, -5.0}, ex, ey, 3.0, 10.0, ctx);
    const Result<OperationResult> aboutX =
        MakeCylinderSolid(Point3{-6.0, 0.0, 0.0}, ey, ez, 3.0, 12.0, ctx);
    const Result<OperationResult> aboutY =
        MakeCylinderSolid(Point3{0.0, -6.0, 0.0}, ez, ex, 3.0, 12.0, ctx);
    if (!aboutZ || !aboutX || !aboutY) {
        std::println(stderr, "a cylinder could not be made");
        return 1;
    }
    const Shape& z = aboutZ->Shape;
    const Shape& x = aboutX->Shape;

    // Two cylinders: their union, their common, the Steinmetz solid, and
    // each less the other.
    const std::array cases{
        Case{"z u x", &z, &x, BooleanOperation::Union, InEither},
        Case{"z n x", &z, &x, BooleanOperation::Common, InBoth},
        Case{"z - x", &z, &x, BooleanOperation::Cut, InZAlone},
        Case{"x - z", &x, &z, BooleanOperation::Cut, InXAlone},
    };
    bool good = true;
    for (const Case& c : cases) {
        good = Checked(c.Name, Boolean(*c.A, *c.B, c.Operation, ctx), c.Expected, ctx) && good;
    }

    // Three: the tricylinder is the Steinmetz solid inside the third.
    const Result<OperationResult> common = Boolean(z, x, BooleanOperation::Common, ctx);
    if (!common) {
        return 1;
    }
    good =
        Checked("z n x n y", Boolean(common->Shape, aboutY->Shape, BooleanOperation::Common, ctx),
                InAllThree, ctx) &&
        good;
    return good ? 0 : 1;
}
