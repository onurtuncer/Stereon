// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The program the debugger pretty-printer tests inspect
// (tools/debugger/, tests/debugger/check_printers.py). It builds one value
// of every type the printers know, then calls StereonDebuggerStop(), where
// the scripted debugger session stops and prints them. The values are fixed
// so the expected output is too.

#include <memory>

import Stereon.Surface;

namespace {

/// A curve type that is not built in, to show how the printers describe it.
struct UserCurve {
    Stereon::Line Base;
};
Stereon::Point3 Eval(const UserCurve& c, const double t)
{
    return Eval(c.Base, t);
}
Stereon::Vec3 Deriv(const UserCurve& c, const double t, const int n)
{
    return Deriv(c.Base, t, n);
}
Stereon::Interval Domain(const UserCurve& c)
{
    return Domain(c.Base);
}
Stereon::Box3 BoundingBox(const UserCurve& c)
{
    return BoundingBox(c.Base);
}
double Project(const UserCurve& c, const Stereon::Point3& p)
{
    return Project(c.Base, p);
}
UserCurve Reversed(const UserCurve& c)
{
    return {Reversed(c.Base)};
}
UserCurve Transformed(const UserCurve& c, const Stereon::Transform& x)
{
    return {Transformed(c.Base, x)};
}

} // namespace

#if defined(_MSC_VER)
#define STEREON_DEBUGGER_NOINLINE __declspec(noinline)
#else
#define STEREON_DEBUGGER_NOINLINE [[gnu::noinline]]
#endif

namespace {

/// Written by StereonDebuggerStop() so the call, and every value whose
/// address it receives, survives optimisation.
const void* volatile s_Anchor = nullptr;

} // namespace

// Not static and not inlined, so the debugger can stop in it.
STEREON_DEBUGGER_NOINLINE void StereonDebuggerStop(const void* anchor)
{
    s_Anchor = anchor;
}

int main()
{
    using namespace Stereon;

    const Vec3 vector{1.0, 2.5, -3.0};
    const Point3 point{4.0, 5.0, 6.0};
    const Point2 uv{0.5, -1.5};
    const Interval interval{1.0, 2.0};
    const Interval emptyInterval{};
    const Box3 box = Box3::FromCorners(Point3{0.0, 0.0, 0.0}, Point3{1.0, 2.0, 3.0});
    const Box2 emptyBox{};

    const VertexId vertex = VertexId::FromIndex(12);
    const EdgeId invalidEdge{};
    const OrientedCoedgeId coedge = OrientedCoedgeId::Make(CoedgeId::FromIndex(42), true);
    GenerationTable<EntityKind::Face> faces;
    const BuilderId<EntityKind::Face> face = faces.Allocate();

    const Line line{Point3{0.0, 0.0, 1.0}, Vec3{1.0, 0.0, 0.0}};
    const Circle circle{Point3{1.0, 2.0, 3.0}, 2.5, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};
    const Cylinder cylinder{Point3{}, 1.5, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, false};

    const Curve lineCurve(line);
    const Curve circleCurve(circle);
    const Curve userCurve(UserCurve{line});
    const Curve emptyCurve;
    const Curve2d planarCurve(Line2d{Point2{0.0, 0.0}, Vec2{0.0, 1.0}});
    const Surface cylinderSurface(cylinder);

    PagedArena<int, 4> arena;
    for (int i = 0; i < 6; ++i) {
        arena.PushBack(i * 10);
    }
    const PagedArena<int, 4> snapshot = arena;
    arena.Set(5, 99);

    const JournalId journal{7};
    const Sign sign = Sign::Negative;

    StereonDebuggerStop(&sign);
    // Keep every value alive across the stop.
    StereonDebuggerStop(&vector);
    StereonDebuggerStop(&point);
    StereonDebuggerStop(&uv);
    StereonDebuggerStop(&interval);
    StereonDebuggerStop(&emptyInterval);
    StereonDebuggerStop(&box);
    StereonDebuggerStop(&emptyBox);
    StereonDebuggerStop(&vertex);
    StereonDebuggerStop(&invalidEdge);
    StereonDebuggerStop(&coedge);
    StereonDebuggerStop(&face);
    StereonDebuggerStop(&lineCurve);
    StereonDebuggerStop(&circleCurve);
    StereonDebuggerStop(&userCurve);
    StereonDebuggerStop(&emptyCurve);
    StereonDebuggerStop(&planarCurve);
    StereonDebuggerStop(&cylinderSurface);
    StereonDebuggerStop(&arena);
    StereonDebuggerStop(&snapshot);
    StereonDebuggerStop(&journal);
    return 0;
}
