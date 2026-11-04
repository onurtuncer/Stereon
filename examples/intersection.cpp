// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Curves against surfaces (docs/examples.rst): a line through a sphere, and
// a B-spline curve crossing a plane. Every point comes with its parameter on
// the curve, its parameters on the surface, its position, how the two meet
// there, and whether that was decided exactly.

#include <cstdio>
// Not used directly: GCC 16.2 on MinGW defines std::make_shared's tag twice
// unless a unit that reaches it through a module includes <memory> itself
// (docs/toolchain/gcc16-feature-probe.md).
#include <memory>
#include <print>
#include <string_view>

import Stereon.Intersect;

using namespace Stereon;

namespace {

std::string_view ContactName(const ContactType contact)
{
    return contact == ContactType::Crossing ? "crossing" : "touch";
}

void Print(const char* what, const CurveSurfaceIntersection& result)
{
    std::println("{}: {} point(s), {} overlap(s)", what, result.Points.size(),
                 result.Overlaps.size());
    for (const CurveSurfacePoint& p : result.Points) {
        std::println("  t = {:.9f}  (u, v) = ({:.6f}, {:.6f})  at ({:.6f}, {:.6f}, {:.6f})  {}{}",
                     p.Parameter, p.OnSurface.X, p.OnSurface.Y, p.Position.X, p.Position.Y,
                     p.Position.Z, ContactName(p.Contact), p.Decided ? "" : ", undecided");
    }
}

} // namespace

int main()
{
    const Context ctx = Context::Default();

    // A line along x, 0.6 off the centre of the unit sphere: two crossings.
    const Line line{Point3{-2.0, 0.6, 0.0}, Vec3{1.0, 0.0, 0.0}};
    const Sphere sphere{Point3{}, 1.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    const Result<CurveSurfaceIntersection> onSphere = Intersect(Curve{line}, Surface{sphere}, ctx);
    if (!onSphere) {
        std::println(stderr, "line x sphere: {}", onSphere.error().GetMessage());
        return 1;
    }
    Print("line x sphere", *onSphere);

    // An arch, z(s) = -1 + 8 s - 8 s^2 over s in [0, 1], crossing z = 0 twice.
    const Result<BSplineCurve> arch = MakeBSplineCurve(
        2, {0.0, 0.0, 0.0, 1.0, 1.0, 1.0}, {{0.0, 0.0, -1.0}, {1.0, 0.0, 3.0}, {2.0, 0.0, -1.0}});
    const Plane ground{Point3{}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    if (!arch) {
        std::println(stderr, "MakeBSplineCurve: {}", arch.error().GetMessage());
        return 1;
    }
    const Result<CurveSurfaceIntersection> onPlane = Intersect(Curve{*arch}, Surface{ground}, ctx);
    if (!onPlane) {
        std::println(stderr, "arch x plane: {}", onPlane.error().GetMessage());
        return 1;
    }
    Print("arch x plane", *onPlane);
    return onSphere->Points.size() == 2 && onPlane->Points.size() == 2 ? 0 : 1;
}
