// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The dumps behind the manual's stn-view screenshots (docs/debugging.rst):
//
//     stereon_viewer_dumps <directory>
//
// writes <directory>/geometry (a B-spline curve with a kink, a B-spline
// patch, a circle and two plane curves) and <directory>/topology (the
// closed cylinder of the topology tests, as ShapeArchive stores it). The
// viewer job of .github/workflows/linux.yml renders each with
// `stn-view --screenshot` and uploads the images.

#include "../../tests/topo/BodiesStd.hpp"

#include <cmath>
#include <filesystem>
#include <iostream>
#include <numbers>
#include <string>
#include <vector>

import Stereon.Io;

#include "../../tests/topo/Bodies.hpp"

namespace {

using Stereon::Point2;
using Stereon::Point3;
using Stereon::Vec2;
using Stereon::Vec3;

bool Save(const std::filesystem::path& directory, const std::string& name,
          Stereon::GeometryArchive archive)
{
    std::filesystem::create_directories(directory);
    if (!archive.Debug) {
        archive.Debug.emplace();
    }
    archive.Debug->Name = name;
    archive.Debug->File = "docs/figures/ViewerDumps.cpp";
    const auto saved = Stereon::SaveGeometry(directory / (name + ".stn"), archive);
    if (!saved) {
        std::cerr << name << ": " << saved.error().GetMessage() << '\n';
    }
    return saved.has_value();
}

/// A cubic through a kink: the knot 0.5 has multiplicity three.
Stereon::BSplineCurve KinkedCubic()
{
    return Stereon::BSplineCurve{
        3,
        {0.0, 0.0, 0.0, 0.0, 0.25, 0.5, 0.5, 0.5, 0.75, 1.0, 1.0, 1.0, 1.0},
        {Point3{0.0, 0.0, 0.0}, Point3{1.0, 2.0, 0.5}, Point3{2.5, -1.0, 1.0},
         Point3{3.5, 1.5, 1.5}, Point3{4.5, 0.0, 0.5}, Point3{5.5, 2.0, 0.0},
         Point3{6.5, -1.5, 1.0}, Point3{7.5, 1.0, 0.0}, Point3{8.5, 0.0, 0.5}},
        {}};
}

/// A bicubic 5 x 5 patch over [0, 6] x [4, 10], gently waved.
Stereon::BSplineSurface WavedPatch()
{
    Stereon::BSplineSurface s;
    s.DegreeU = 3;
    s.DegreeV = 3;
    s.KnotsU = {0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0};
    s.KnotsV = s.KnotsU;
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            const double x = 1.5 * i;
            const double y = 4.0 + 1.5 * j;
            s.Poles.push_back(Point3{x, y, 0.8 * std::sin(0.9 * x) * std::cos(0.7 * y)});
        }
    }
    return s;
}

bool WriteGeometry(const std::filesystem::path& directory)
{
    Stereon::GeometryArchive curve;
    curve.Curves.push_back(Stereon::Curve(KinkedCubic()));
    Stereon::GeometryArchive patch;
    patch.Surfaces.push_back(Stereon::Surface(WavedPatch()));
    Stereon::GeometryArchive circle;
    circle.Curves.push_back(
        Stereon::Curve(Stereon::Circle{Point3{8.0, 7.0, 0.0}, 1.5, Vec3{1, 0, 0}, Vec3{0, 1, 0}}));
    Stereon::GeometryArchive plane;
    plane.Curves2d.push_back(
        Stereon::Curve2d(Stereon::Circle2d{Point2{0.5, 0.5}, 0.3, Vec2{1, 0}, Vec2{0, 1}}));
    plane.Curves2d.push_back(Stereon::Curve2d(Stereon::BSplineCurve2d{
        2,
        {0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0},
        {Point2{0.0, 0.0}, Point2{0.4, 1.0}, Point2{0.7, -0.2}, Point2{1.0, 0.6}},
        {}}));
    return Save(directory, "kinked-cubic", curve) && Save(directory, "waved-patch", patch) &&
           Save(directory, "circle", circle) && Save(directory, "plane-curves", plane);
}

bool WriteTopology(const std::filesystem::path& directory)
{
    const auto cylinder = Stereon::ShapeArchive(TopoTest::MakeCylinder().Freeze());
    if (!cylinder) {
        std::cerr << "cylinder: " << cylinder.error().GetMessage() << '\n';
        return false;
    }
    return Save(directory, "closed-cylinder", *cylinder);
}

} // namespace

int main(const int argc, char** argv)
{
    if (argc != 2) {
        std::cerr << "usage: stereon_viewer_dumps <directory>\n";
        return 2;
    }
    const std::filesystem::path root(argv[1]);
    return WriteGeometry(root / "geometry") && WriteTopology(root / "topology") ? 0 : 1;
}
