// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// A box drilled through by a cylinder (docs/examples.rst): two primitive
// solids, a Boolean cut, and the result checked by the validity checker.
// Every step returns a Result; an operation either gives a valid shape or
// says why not.

#include <cstdio>
#include <print>

import Stereon.Boolean;

using namespace Stereon;

namespace {

void Describe(const char* what, const Shape& shape)
{
    std::println("{}: {} faces, {} edges, {} vertices", what, shape.Count<EntityKind::Face>(),
                 shape.Count<EntityKind::Edge>(), shape.Count<EntityKind::Vertex>());
}

} // namespace

int main()
{
    const Context ctx = Context::Default();

    // A 40 x 30 x 20 mm box, and a cylinder of radius 6 mm standing through
    // it, one millimetre proud of either face.
    const Result<OperationResult> box =
        MakeBoxSolid(Box3::FromCorners({0.0, 0.0, 0.0}, {40.0, 30.0, 20.0}), ctx);
    const Result<OperationResult> pin = MakeCylinderSolid(
        Point3{20.0, 15.0, -1.0}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, 6.0, 22.0, ctx);
    if (!box || !pin) {
        std::println(stderr, "primitive: {}", (box ? pin : box).error().GetMessage());
        return 1;
    }

    const Result<OperationResult> drilled =
        Boolean(box->Shape, pin->Shape, BooleanOperation::Cut, ctx);
    if (!drilled) {
        std::println(stderr, "cut: {}", drilled.error().GetMessage());
        return 1;
    }
    Describe("box", box->Shape);
    Describe("box - cylinder", drilled->Shape);

    // The result passed the Geometric level when it was published; check it
    // again, with the Context's tolerance rules, as a caller may.
    const ValidityReport report = ValidateTopology(drilled->Shape, ValidityLevel::Geometric, ctx);
    std::println("valid at the Geometric level: {}", report.IsValid());
    return report.IsValid() ? 0 : 1;
}
