// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Freeze() as the certification gate of ADR-0038: what it certifies, what it
// refuses, what it keeps, and the empty Shape.

#include "BodiesStd.hpp"

#include <catch2/catch_test_macros.hpp>

#include <set>
#include <string>
#include <utility>
#include <vector>

import Stereon.Topo;

#include "Bodies.hpp"

using namespace TopoTest;
using Stereon::EntityKind;
using Stereon::ErrorCode;
using Stereon::Shape;
using Stereon::ShapeBuilder;

TEST_CASE("Freeze certifies the snapshot it makes", "[topo][freeze]")
{
    const Shape shape = MakeBox().Freeze();
    REQUIRE_FALSE(shape.IsEmpty());
    const Stereon::ValidationCertificate& certificate = shape.GetSnapshot()->GetCertificate();
    // CI probe: the Linux GCC jobs read a certificate the test TU lays out
    // differently from the unit that wrote it. Shown only on failure.
    using Stereon::ValidationCertificate;
    const char* snapshotBytes = reinterpret_cast<const char*>(shape.GetSnapshot().get());
    INFO("sizeof(ValidationCertificate) "
         << sizeof(ValidationCertificate) << ", offsets "
         << offsetof(ValidationCertificate, ValidatorVersion) << '/'
         << offsetof(ValidationCertificate, Invariants) << '/'
         << offsetof(ValidationCertificate, Snapshot) << "; sizeof(TopologyTables) "
         << sizeof(Stereon::Detail::TopologyTables) << ", sizeof(TopologySnapshot) "
         << sizeof(Stereon::TopologySnapshot) << ", certificate at snapshot+"
         << (reinterpret_cast<const char*>(&certificate) - snapshotBytes) << "; raw: version "
         << certificate.ValidatorVersion << ", invariants size " << certificate.Invariants.size()
         << " data " << static_cast<const void*>(certificate.Invariants.data()) << ", snapshot "
         << certificate.Snapshot << ", source builder " << shape.GetSnapshot()->GetSourceBuilder());
    REQUIRE(certificate.ValidatorVersion == Stereon::Detail::StructuralValidatorVersion);
    REQUIRE(certificate.Invariants == "Structural S01-S15");
    REQUIRE(certificate.Snapshot != 0);
    REQUIRE(certificate.Snapshot == shape.GetSnapshot()->GetStore());
    // The independent validator agrees with the certificate.
    REQUIRE(Stereon::ValidateTopology(shape).IsValid());
}

TEST_CASE("every Freeze makes a snapshot of its own", "[topo][freeze]")
{
    std::set<Stereon::StoreId> ids;
    for (int i = 0; i < 20; ++i) {
        ids.insert(MakeSquareSheet().Freeze().GetSnapshot()->GetStore());
    }
    REQUIRE(ids.size() == 20);
}

TEST_CASE("Freeze consumes the builder whether it succeeds or not", "[topo][freeze]")
{
    SECTION("success")
    {
        Body body = MakeSquareSheet();
        REQUIRE(Stereon::Freeze(std::move(body.Builder)).has_value());
        REQUIRE_FALSE(body.Builder.GetToken().IsValid()); // NOLINT(bugprone-use-after-move)
        REQUIRE_FALSE(body.Builder.AddVertex({0, 0, 0}, 0.0).has_value());
        REQUIRE(body.Builder.Count<EntityKind::Vertex>() == 0);
    }
    SECTION("failure")
    {
        ShapeBuilder b;
        const BuilderFaceId f = Must(b.AddFace(b.AddSurface(XyPlane())));
        static_cast<void>(Must(b.AddLoop(f)));
        REQUIRE_FALSE(Stereon::Freeze(std::move(b)).has_value());
        REQUIRE_FALSE(b.GetToken().IsValid()); // NOLINT(bugprone-use-after-move)
        REQUIRE_FALSE(Stereon::Freeze(std::move(b)).has_value());
    }
}

TEST_CASE("a refused Freeze reports the findings by rule", "[topo][freeze]")
{
    ShapeBuilder b;
    const BuilderFaceId f = Must(b.AddFace(b.AddSurface(XyPlane())));
    static_cast<void>(Must(b.AddLoop(f)));
    const Stereon::ValidityReport report = b.Validate();
    REQUIRE_FALSE(report.IsValid());
    REQUIRE(report.GetCount() == 1);
    REQUIRE(report.Has(Stereon::TopologyRule::EmptyLoop));
    REQUIRE(report.GetFindings()[0].Entities.size() == 1);
    REQUIRE(report.GetFindings()[0].Entities[0].starts_with("Loop#0 [Face#0, 0 coedges"));

    const auto frozen = Stereon::Freeze(std::move(b));
    REQUIRE_FALSE(frozen.has_value());
    REQUIRE(frozen.error().GetCode() == ErrorCode::InvalidInput);
    const std::string message(frozen.error().GetMessage());
    REQUIRE(message.find("not structurally valid") != std::string::npos);
    REQUIRE(message.find("S14 EmptyLoop: Loop#0") != std::string::npos);
}

TEST_CASE("Freeze keeps every index", "[topo][freeze]")
{
    Body body = MakeCylinder();
    const Shape shape = std::move(body).Freeze();
    for (const BuilderVertexId v : body.Vertices) {
        REQUIRE(shape.Resolve(v).GetIndex() == v.GetIndex());
    }
    for (const BuilderCoedgeId c : body.Coedges) {
        REQUIRE(shape.Resolve(c).GetIndex() == c.GetIndex());
    }
    for (const BuilderFaceId f : body.Faces) {
        REQUIRE(shape.Resolve(f).GetIndex() == f.GetIndex());
    }
}

TEST_CASE("Resolve refuses handles that name no live entity", "[topo][freeze]")
{
    Body body = MakeSquareSheet();
    const BuilderVertexId v = body.Vertices[0];
    const Shape shape = std::move(body).Freeze();
    REQUIRE(shape.Resolve(v).IsValid());
    REQUIRE_FALSE(shape.Resolve(BuilderVertexId{}).IsValid());
    REQUIRE_FALSE(
        shape.Resolve(BuilderVertexId::Make(v.GetIndex(), v.GetGeneration() + 1, 0)).IsValid());
    REQUIRE_FALSE(shape.Resolve(BuilderVertexId::Make(400, 0, 0)).IsValid());
    REQUIRE_FALSE(Shape{}.Resolve(v).IsValid());
    if constexpr (STEREON_CHECKED_HANDLES) {
        // A live slot of this shape, but a handle another builder made.
        ShapeBuilder other;
        const BuilderVertexId foreign = Must(other.AddVertex({0, 0, 0}, 0.0));
        REQUIRE(foreign.GetIndex() == v.GetIndex());
        REQUIRE_FALSE(shape.Resolve(foreign).IsValid());
    }
}

TEST_CASE("an empty builder freezes to a shape with nothing in it", "[topo][freeze]")
{
    const Shape shape = Must(Stereon::Freeze(ShapeBuilder{}));
    REQUIRE_FALSE(shape.IsEmpty());
    REQUIRE(shape.Count<EntityKind::Vertex>() == 0);
    REQUIRE(Must(Stereon::AllEntities<EntityKind::Face>(shape)).empty());
    REQUIRE(Stereon::ValidateTopology(shape).IsValid());
}

TEST_CASE("the empty Shape answers every query with InvalidInput", "[topo][freeze]")
{
    const Shape empty;
    REQUIRE(empty.IsEmpty());
    REQUIRE(empty.Count<EntityKind::Edge>() == 0);
    REQUIRE(Stereon::ValidateTopology(empty).IsValid());
    const auto all = Stereon::AllEntities<EntityKind::Vertex>(empty);
    REQUIRE_FALSE(all.has_value());
    REQUIRE(all.error().GetCode() == ErrorCode::InvalidInput);
    REQUIRE_FALSE(Stereon::LoopsOf(empty, Stereon::FaceId::FromIndex(0)).has_value());
    REQUIRE_FALSE(Stereon::CoedgesOf(empty, Stereon::LoopId::FromIndex(0)).has_value());
    REQUIRE_FALSE(Stereon::RadialCoedgesOf(empty, Stereon::EdgeId::FromIndex(0)).has_value());
    REQUIRE_FALSE(Stereon::FacesOf(empty, Stereon::ShellId::FromIndex(0)).has_value());
    REQUIRE_FALSE(Stereon::ShellsOf(empty, Stereon::SolidId::FromIndex(0)).has_value());
    REQUIRE(Stereon::Describe(empty, Stereon::EdgeId::FromIndex(3)) ==
            "Edge#3 [not in this shape]");
}

TEST_CASE("copies of a Shape share one snapshot", "[topo][freeze]")
{
    const Shape a = MakeBox().Freeze();
    const Shape b = a; // NOLINT(performance-unnecessary-copy-initialization): the point
    REQUIRE(a.GetSnapshot() == b.GetSnapshot());
    REQUIRE(a.GetSnapshot().use_count() == 2);
}
