// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The Boolean prototype's M11: union, cut and common of the M10 fixtures
// (tests/support/BooleanCases.hpp) against OpenCASCADE's BRepAlgoAPI on the
// same operands. Volumes and areas must agree within 1e-6 relative (ADR-0013)
// plus what the results' tolerances allow: a band as wide as the largest
// tolerance, over the faces for the volume and along the edges for the
// area. Face counts must agree exactly. Stereon's
// results are measured by tests/support/ShapeMeasures.hpp, whose own error
// estimate must stay below a hundredth of that bound, and OCCT's by
// BRepGProp.
//
// The operands are built twice, from the same numbers: here with
// BRepPrimAPI and in BooleanCases with Stereon.Build. Each OCCT operand's
// volume and area are checked against the Stereon operand's first, so the
// twins cannot drift apart unnoticed. The loosened box (M9) sits 0.8 of its
// tolerance from the box's face; OCCT is given that tolerance as its fuzzy
// value, its counterpart of ADR-0001's pair tolerance.
//
// OCCT is LGPL-2.1 and linked only into this test (ADR-0008); build it with
// tools/occt/build_occt.sh and configure with -DSTEREON_WITH_OCCT=ON.

#include <BRepAlgoAPI_BooleanOperation.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepGProp.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <GProp_GProps.hxx>
#include <NCollection_IndexedMap.hxx>
#include <NCollection_List.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopExp.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS_Shape.hxx>
#include <gp_Ax2.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <limits>
#include <memory>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

import Stereon.Boolean;

#include "BooleanCases.hpp"
#include "ShapeMeasures.hpp"

using Stereon::BooleanOperation;
using Stereon::EntityKind;

namespace {

TopoDS_Shape OcctBox(const gp_Pnt& lo, const gp_Pnt& hi)
{
    return BRepPrimAPI_MakeBox(lo, hi).Shape();
}

/// A cylinder standing on (x, y, z) along +z, as BooleanCases builds them.
TopoDS_Shape OcctCylinder(const double x, const double y, const double z, const double radius,
                          const double height)
{
    return BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(x, y, z), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0)),
                                    radius, height)
        .Shape();
}

/// The OCCT twins of one BooleanCases case, by its name's first word: the
/// object, the tools and the fuzzy value.
struct Twin {
    TopoDS_Shape A;
    NCollection_List<TopoDS_Shape> B;
    double Fuzzy = 0.0;
};

Twin TwinOf(const std::string& name)
{
    const TopoDS_Shape box = OcctBox(gp_Pnt(0, 0, 0), gp_Pnt(10, 10, 10));
    Twin twin;
    twin.A = box;
    if (name.starts_with("corner ")) {
        twin.B.Append(OcctBox(gp_Pnt(5, 5, 5), gp_Pnt(15, 15, 15)));
    } else if (name.starts_with("drum ")) {
        twin.B.Append(OcctCylinder(5, 5, -2, 3, 15));
    } else if (name.starts_with("loose ")) {
        const double t = 1e-3;
        twin.B.Append(OcctBox(gp_Pnt(10 + 0.8 * t, 0, 0), gp_Pnt(20, 10, 10)));
        twin.Fuzzy = t;
    } else if (name == "drilled slab Cut") {
        const int n = 3;
        twin.A = OcctBox(gp_Pnt(0, 0, 0), gp_Pnt(10.0 * n, 10.0 * n, 10));
        for (int k = 0; k < n * n; ++k) {
            twin.B.Append(OcctCylinder(10.0 * (k % n) + 5.0, 10.0 * (k / n) + 5.0, -5, 2, 20));
        }
    } else {
        throw std::runtime_error("no OCCT twin for the case " + name);
    }
    return twin;
}

TopoDS_Shape OcctBoolean(const Twin& twin, const BooleanOperation operation)
{
    NCollection_List<TopoDS_Shape> arguments;
    arguments.Append(twin.A);
    std::unique_ptr<BRepAlgoAPI_BooleanOperation> op;
    switch (operation) {
    case BooleanOperation::Union:
        op = std::make_unique<BRepAlgoAPI_Fuse>();
        break;
    case BooleanOperation::Common:
        op = std::make_unique<BRepAlgoAPI_Common>();
        break;
    case BooleanOperation::Cut:
        op = std::make_unique<BRepAlgoAPI_Cut>();
        break;
    }
    op->SetArguments(arguments);
    op->SetTools(twin.B);
    if (twin.Fuzzy > 0.0) {
        op->SetFuzzyValue(twin.Fuzzy);
    }
    op->Build();
    if (!op->IsDone() || op->HasErrors()) {
        throw std::runtime_error("BRepAlgoAPI failed");
    }
    return op->Shape();
}

struct OcctMeasures {
    double Volume = 0.0;
    double Area = 0.0;
    std::size_t Faces = 0;
};

OcctMeasures MeasuresOf(const TopoDS_Shape& shape)
{
    OcctMeasures m;
    GProp_GProps volume;
    BRepGProp::VolumeProperties(shape, volume);
    m.Volume = volume.Mass();
    GProp_GProps area;
    BRepGProp::SurfaceProperties(shape, area);
    m.Area = area.Mass();
    NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> faces;
    TopExp::MapShapes(shape, TopAbs_FACE, faces);
    m.Faces = static_cast<std::size_t>(faces.Extent());
    return m;
}

template <EntityKind K>
double LargestTolerance(const Stereon::Shape& shape)
{
    double largest = 0.0;
    for (const auto h : *AllEntities<K>(shape)) {
        largest = std::max(largest, ToleranceOf(shape, h));
    }
    return largest;
}

double LargestTolerance(const Stereon::Shape& shape)
{
    return std::max(LargestTolerance<EntityKind::Vertex>(shape),
                    LargestTolerance<EntityKind::Edge>(shape));
}

/// ADR-0013's 1e-6 relative, plus a band of the given width over `extent`
/// (the faces' area for a volume, the edges' length for an area).
double Bound(const double value, const double extent, const double width)
{
    return 1e-6 * std::abs(value) + extent * width;
}

} // namespace

TEST_CASE("Boolean fixtures: volumes, areas and face counts agree with BRepAlgoAPI",
          "[occt][boolean]")
{
    for (const BooleanCases::Case& c : BooleanCases::Cases()) {
        INFO(c.Name);
        const Twin twin = TwinOf(c.Name);

        // The twins first: each OCCT operand measures as the Stereon one.
        {
            const ShapeMeasures::Measures a = ShapeMeasures::MeasuresOf(c.A);
            const OcctMeasures twinA = MeasuresOf(twin.A);
            CHECK(std::abs(a.Volume - twinA.Volume) <= 1e-9 * twinA.Volume);
            CHECK(std::abs(a.Area - twinA.Area) <= 1e-9 * twinA.Area);
            const ShapeMeasures::Measures b = ShapeMeasures::MeasuresOf(c.B);
            double volumeB = 0.0;
            double areaB = 0.0;
            for (const TopoDS_Shape& tool : twin.B) {
                volumeB += MeasuresOf(tool).Volume;
                areaB += MeasuresOf(tool).Area;
            }
            CHECK(std::abs(b.Volume - volumeB) <= 1e-9 * volumeB);
            CHECK(std::abs(b.Area - areaB) <= 1e-9 * areaB);
        }

        const auto result = Boolean(c.A, c.B, c.Operation, BooleanCases::serial);
        REQUIRE(result.has_value());
        const ShapeMeasures::Measures ours = ShapeMeasures::MeasuresOf(result->Shape);
        const OcctMeasures theirs = MeasuresOf(OcctBoolean(twin, c.Operation));
        const std::size_t faces = AllEntities<EntityKind::Face>(result->Shape)->size();

        const double width = std::max(LargestTolerance(result->Shape), twin.Fuzzy);
        const double volumeBound = Bound(theirs.Volume, std::max(ours.Area, theirs.Area), width);
        const double areaBound = Bound(theirs.Area, ours.EdgeLength, width);
        std::printf("%-20s volume %.12g / %.12g, area %.12g / %.12g, faces %zu / %zu\n",
                    c.Name.c_str(), ours.Volume, theirs.Volume, ours.Area, theirs.Area, faces,
                    theirs.Faces);
        INFO("Stereon volume " << ours.Volume << ", OCCT " << theirs.Volume << ", bound "
                               << volumeBound);
        INFO("Stereon area " << ours.Area << ", OCCT " << theirs.Area << ", bound " << areaBound);
        REQUIRE(ours.VolumeEstimate <= 0.01 * volumeBound);
        REQUIRE(ours.AreaEstimate <= 0.01 * areaBound);
        CHECK(std::abs(ours.Volume - theirs.Volume) <= volumeBound);
        CHECK(std::abs(ours.Area - theirs.Area) <= areaBound);
        CHECK(faces == theirs.Faces);
    }
}
