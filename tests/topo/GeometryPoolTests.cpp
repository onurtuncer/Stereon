// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The geometry pools of a topology store (ADR-0004, B0.4): built-in types
// kept by value in their own arenas, everything else as wrappers, behind one
// dense index; read back wrapped or as the stored value, and shared
// copy-on-write like the entity tables.

#include "BodiesStd.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <string>
#include <type_traits>
#include <vector>

import Stereon.Topo;

#include "Bodies.hpp"

using Stereon::BSplineCurve;
using Stereon::Circle;
using Stereon::Curve;
using Stereon::Line;
using Stereon::Point3;
using Stereon::Vec3;
using Stereon::Detail::CurvePool;
using TopoTest::Deref;

namespace {

BSplineCurve Quadratic()
{
    return BSplineCurve{
        2, {0.0, 0.0, 0.0, 1.0, 1.0, 1.0}, {Point3{0, 0, 0}, Point3{1, 2, 0}, Point3{2, 0, 0}}, {}};
}

/// What Visit sees at the index: the stored type's name.
std::string StoredAs(const CurvePool& pool, const std::size_t index)
{
    return pool.Visit(index, [](const auto& g) -> std::string {
        using G = std::remove_cvref_t<decltype(g)>;
        if constexpr (std::is_same_v<G, Line>) {
            return "Line";
        } else if constexpr (std::is_same_v<G, Circle>) {
            return "Circle";
        } else if constexpr (std::is_same_v<G, Curve>) {
            return "wrapper";
        } else {
            return "other built-in";
        }
    });
}

} // namespace

TEST_CASE("a geometry pool keeps built-in types unwrapped behind dense indices", "[topo][memory]")
{
    CurvePool pool;
    const Line line{Point3{1, 2, 3}, Vec3{0, 0, 1}};
    const Circle circle{Point3{}, 2.0, Vec3{1, 0, 0}, Vec3{0, 1, 0}};
    REQUIRE(pool.PushBack(Curve(line)) == 0);
    REQUIRE(pool.PushBack(Curve(Quadratic())) == 1);
    REQUIRE(pool.PushBack(Curve(circle)) == 2);
    REQUIRE(pool.PushBack(Curve{}) == 3);
    REQUIRE(pool.GetSize() == 4);

    REQUIRE(StoredAs(pool, 0) == "Line");
    REQUIRE(StoredAs(pool, 1) == "wrapper");
    REQUIRE(StoredAs(pool, 2) == "Circle");
    REQUIRE(StoredAs(pool, 3) == "wrapper");

    // Read back wrapped, the same geometry whichever arena holds it.
    REQUIRE(Deref(pool[0].Target<Line>()).Origin == line.Origin);
    REQUIRE(Deref(pool[1].Target<BSplineCurve>()).Poles == Quadratic().Poles);
    REQUIRE(Eval(pool[2], 0.0) == Eval(circle, 0.0));
    REQUIRE(pool[3].IsEmpty());
}

TEST_CASE("replacing a pool entry may move it to another arena", "[topo][memory]")
{
    CurvePool pool;
    static_cast<void>(pool.PushBack(Curve(Line{Point3{}, Vec3{1, 0, 0}})));
    pool.Set(0, Curve(Quadratic()));
    REQUIRE(StoredAs(pool, 0) == "wrapper");
    REQUIRE(pool[0].Target<BSplineCurve>() != nullptr);
    pool.Set(0, Curve{});
    REQUIRE(pool[0].IsEmpty());
    REQUIRE(pool.GetSize() == 1);
}

TEST_CASE("a copied pool shares its pages until it is written", "[topo][memory]")
{
    CurvePool pool;
    for (int i = 0; i < 600; ++i) {
        static_cast<void>(
            pool.PushBack(Curve(Line{Point3{0, 0, static_cast<double>(i)}, Vec3{1, 0, 0}})));
    }
    CurvePool copy = pool;
    REQUIRE(SharedPageCount(pool, copy) == pool.GetPageCount());
    copy.Set(0, Curve(Line{Point3{9, 9, 9}, Vec3{0, 1, 0}}));
    // The slot page holding index 0 and the line arena's last page are
    // copied; every other page is still shared.
    REQUIRE(SharedPageCount(pool, copy) + 2 >= pool.GetPageCount());
    REQUIRE(SharedPageCount(pool, copy) < copy.GetPageCount());
    REQUIRE(Deref(pool[0].Target<Line>()).Origin == Point3{0, 0, 0});
    REQUIRE(Deref(copy[0].Target<Line>()).Origin == Point3{9, 9, 9});
}

TEST_CASE("a shape's pools hold its lines and planes without wrappers", "[topo][memory]")
{
    const Stereon::Shape box = TopoTest::MakeBox().Freeze();
    const Stereon::Detail::TopologyTables& t = box.GetSnapshot()->GetTables();
    for (std::size_t i = 0; i < t.Curves.GetSize(); ++i) {
        REQUIRE(StoredAs(t.Curves, i) == "Line");
    }
    // A line pcurve and its slot take less than one wrapper.
    STATIC_REQUIRE(sizeof(Stereon::Line2d) + sizeof(std::uint32_t) < sizeof(Stereon::Curve2d));
    // The accessors still give every edge its curve, wrapped.
    for (const Stereon::EdgeId e :
         TopoTest::Must(Stereon::AllEntities<Stereon::EntityKind::Edge>(box))) {
        REQUIRE(Deref(Stereon::CurveOf(box, e)).Target<Line>() != nullptr);
    }
}
