// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Include after importing the geometry modules. Evaluators deliberately
// violate their numerical promises so validation must refuse their output.
namespace FaultGeometry {

struct Curve {
    Stereon::Line Base{{0, 0, 0}, {1, 0, 0}};
    double At = 0.5;
    double Value = 0.0;
};
inline Stereon::Point3 Eval(const Curve& c, const double t)
{
    auto point = Stereon::Eval(c.Base, t);
    if (t == c.At) {
        point.Y = c.Value;
    }
    return point;
}
inline Stereon::Vec3 Deriv(const Curve& c, const double t, const int n)
{
    return Stereon::Deriv(c.Base, t, n);
}
inline Stereon::Interval Domain(const Curve&)
{
    return {0.0, 1.0};
}
inline Stereon::Box3 BoundingBox(const Curve& c)
{
    return Stereon::BoundingBox(c.Base);
}
inline double Project(const Curve& c, const Stereon::Point3& p)
{
    return Stereon::Project(c.Base, p);
}
inline Curve Reversed(const Curve& c)
{
    return {Stereon::Reversed(c.Base), c.At, c.Value};
}
inline Curve Transformed(const Curve& c, const Stereon::Transform& x)
{
    return {Stereon::Transformed(c.Base, x), c.At, c.Value};
}

struct Surface {
    Stereon::Surface Base;
    double At = 0.5;
    double Value = 0.0;
};
inline Stereon::Point3 Eval(const Surface& s, const Stereon::Point2 uv)
{
    auto point = Stereon::Eval(s.Base, uv);
    if (uv.X == s.At) {
        point.Z = s.Value;
    }
    return point;
}
inline Stereon::Vec3 Deriv(const Surface& s, const Stereon::Point2 uv, const int du, const int dv)
{
    return Stereon::Deriv(s.Base, uv, du, dv);
}
inline Stereon::Box2 Domain(const Surface& s)
{
    return Stereon::Domain(s.Base);
}
inline Stereon::Box3 BoundingBox(const Surface& s)
{
    return Stereon::BoundingBox(s.Base);
}
inline Stereon::Point2 Project(const Surface& s, const Stereon::Point3& p)
{
    return Stereon::Project(s.Base, p);
}
inline Surface Reversed(const Surface& s)
{
    return {Stereon::Reversed(s.Base), s.At, s.Value};
}
inline Surface Transformed(const Surface& s, const Stereon::Transform& x)
{
    return {Stereon::Transformed(s.Base, x), s.At, s.Value};
}
inline Stereon::Vec2 Periods(const Surface& s)
{
    return Stereon::Periods(s.Base);
}

} // namespace FaultGeometry
