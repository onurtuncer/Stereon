// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// What the Boolean tests share: the Build primitives as fixtures, the face
// of a shape by its surface, and one branch of a face's surface against
// another surface. Included after `import Stereon.Boolean;` and the Catch2
// headers.

#pragma once

namespace BooleanTest {

inline const Stereon::Context ctx = Stereon::Context::Default();

/// The box [0, 10]^3.
[[nodiscard]] inline Stereon::Shape Box()
{
    const auto box = MakeBoxSolid(
        Stereon::Box3::FromCorners(Stereon::Point3{}, Stereon::Point3{10.0, 10.0, 10.0}), ctx);
    REQUIRE(box.has_value());
    return box->Shape;
}

/// A cylinder of radius 3 and height 10 about z from the origin.
[[nodiscard]] inline Stereon::Shape Drum()
{
    const auto drum = MakeCylinderSolid(Stereon::Point3{}, Stereon::Vec3{1.0, 0.0, 0.0},
                                        Stereon::Vec3{0.0, 1.0, 0.0}, 3.0, 10.0, ctx);
    REQUIRE(drum.has_value());
    return drum->Shape;
}

/// The face whose plane's normal points along n.
[[nodiscard]] inline Stereon::FaceId PlaneFaceAlong(const Stereon::Shape& shape,
                                                    const Stereon::Vec3& n)
{
    const auto faces = AllEntities<Stereon::EntityKind::Face>(shape);
    REQUIRE(faces.has_value());
    for (const Stereon::FaceId face : *faces) {
        const Stereon::Surface surface = SurfaceOf(shape, face);
        if (const Stereon::Plane* plane = surface.Target<Stereon::Plane>();
            plane != nullptr && Dot(Normalized(Cross(plane->XAxis, plane->YAxis)), n) > 0.9) {
            return face;
        }
    }
    FAIL("no face along the direction");
    return {};
}

/// The face on a cylinder.
[[nodiscard]] inline Stereon::FaceId CylinderFace(const Stereon::Shape& shape)
{
    const auto faces = AllEntities<Stereon::EntityKind::Face>(shape);
    REQUIRE(faces.has_value());
    for (const Stereon::FaceId face : *faces) {
        if (SurfaceOf(shape, face).Target<Stereon::Cylinder>() != nullptr) {
            return face;
        }
    }
    FAIL("no cylindrical face");
    return {};
}

/// The one branch of the face's surface against `other`, the face's surface
/// as side A, over the face's parameter box (planes are unbounded).
[[nodiscard]] inline Stereon::IntersectionBranch BranchOn(const Stereon::Shape& shape,
                                                          const Stereon::FaceId face,
                                                          const Stereon::Surface& other,
                                                          const Stereon::Box2& box)
{
    const auto result =
        Intersect(SurfaceOf(shape, face), box, other, Stereon::Box2::Everything(), ctx);
    REQUIRE(result.has_value());
    REQUIRE(result->Branches.size() == 1);
    return result->Branches.front();
}

/// The plane x = `x`.
[[nodiscard]] inline Stereon::Surface VerticalPlane(const double x)
{
    return Stereon::Surface{Stereon::Plane{
        Stereon::Point3{x, 0.0, 0.0}, Stereon::Vec3{0.0, 1.0, 0.0}, Stereon::Vec3{0.0, 0.0, 1.0}}};
}

/// A cylinder about z through (5, 5), a drill into the box's top face.
[[nodiscard]] inline Stereon::Surface Drill(const double radius)
{
    return Stereon::Surface{Stereon::Cylinder{Stereon::Point3{5.0, 5.0, 0.0}, radius,
                                              Stereon::Vec3{1.0, 0.0, 0.0},
                                              Stereon::Vec3{0.0, 1.0, 0.0}, true}};
}

/// The box's face parameters span [0, 10] in both directions; a margin lets
/// the lines be clipped beside the face too.
inline const Stereon::Box2 faceBox{Stereon::Interval{-1.0, 11.0}, Stereon::Interval{-1.0, 11.0}};

/// The live entities of a kind.
template <Stereon::EntityKind K>
[[nodiscard]] std::size_t CountOf(const Stereon::Shape& shape)
{
    const auto all = AllEntities<K>(shape);
    REQUIRE(all.has_value());
    return all->size();
}

} // namespace BooleanTest
