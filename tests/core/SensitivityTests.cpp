// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

import Stereon.Core;

using Stereon::DependencySet;
using Stereon::DerivativeProduct;
using Stereon::DerivativeStatus;
using Stereon::ErrorCode;
using Stereon::OperationCapability;
using Stereon::ParameterId;
using Stereon::ParameterQuantity;
using Stereon::ParameterSet;
using Stereon::Point2;
using Stereon::SampleOrientation;
using Stereon::SensitivityContext;
using Stereon::SensitivityResult;
using Stereon::TrackedSample;
using Stereon::TrackingConvention;

namespace {

constexpr double NaN = std::numeric_limits<double>::quiet_NaN();
constexpr double Infinity = std::numeric_limits<double>::infinity();

ParameterSet TwoParameters()
{
    ParameterSet parameters;
    REQUIRE(parameters.Add("radius", ParameterQuantity::Length, 10.0));
    REQUIRE(parameters.Add("angle", ParameterQuantity::Angle));
    return parameters;
}

constexpr OperationCapability Undeclared =
    Stereon::Declare("MakeThing", {Stereon::Dependency::DefiningValues},
                     {TrackingConvention::FixedEdgeParameter}, "smooth in its values");

} // namespace

TEST_CASE("parameters are numbered in the order they are added", "[core][sensitivity]")
{
    ParameterSet parameters;
    const auto radius = parameters.Add("radius", ParameterQuantity::Length, 10.0);
    const auto angle = parameters.Add("angle", ParameterQuantity::Angle);
    REQUIRE(radius);
    REQUIRE(angle);
    REQUIRE(radius->GetIndex() == 0);
    REQUIRE(angle->GetIndex() == 1);
    REQUIRE(parameters.GetSize() == 2);
    REQUIRE(parameters.GetName(*angle) == "angle");
    REQUIRE(parameters.GetUnit(*radius) == "mm");
    REQUIRE(parameters.GetUnit(*angle) == "rad");
    REQUIRE(parameters.GetReference(*radius) == 10.0);
    REQUIRE(parameters.GetReference(*angle) == 1.0);
    REQUIRE(parameters.Find("angle") == *angle);
    REQUIRE_FALSE(parameters.Find("height").IsValid());
    REQUIRE(Stereon::UnitOf(ParameterQuantity::Dimensionless) == "1");
}

TEST_CASE("a parameter set rejects bad names and scales", "[core][sensitivity]")
{
    ParameterSet parameters;
    REQUIRE(parameters.Add("radius", ParameterQuantity::Length));
    for (const auto& [name, reference] : std::vector<std::pair<std::string, double>>{
             {"", 1.0}, {"radius", 1.0}, {"h", 0.0}, {"h", -1.0}, {"h", NaN}, {"h", Infinity}}) {
        const auto added = parameters.Add(name, ParameterQuantity::Length, reference);
        REQUIRE_FALSE(added);
        REQUIRE(added.error().GetCode() == ErrorCode::InvalidInput);
    }
    REQUIRE(parameters.GetSize() == 1);
}

TEST_CASE("a parameter id belongs to its set", "[core][sensitivity]")
{
    const ParameterSet parameters = TwoParameters();
    REQUIRE_FALSE(ParameterId{}.IsValid());
    REQUIRE_FALSE(parameters.Contains(ParameterId{}));
    REQUIRE(parameters.Contains(parameters.Find("angle")));
    REQUIRE(parameters.GetName(ParameterId{}).empty());
    REQUIRE(std::isnan(parameters.GetReference(ParameterId{})));

    // An id cannot be made from an index: only a set makes one.
    STATIC_REQUIRE_FALSE(std::is_constructible_v<ParameterId, std::uint32_t, ParameterQuantity>);
}

TEST_CASE("a parameter id of another set is foreign, whatever its index", "[core][sensitivity]")
{
    // The same position and quantity in two sets, with different names and
    // reference scales.
    ParameterSet a;
    ParameterSet b;
    const ParameterId radius = *a.Add("radius", ParameterQuantity::Length, 1.0);
    const ParameterId thickness = *b.Add("thickness", ParameterQuantity::Length, 100.0);
    REQUIRE(radius.GetIndex() == thickness.GetIndex());
    REQUIRE(radius.GetQuantity() == thickness.GetQuantity());
    REQUIRE(radius != thickness);
    REQUIRE_FALSE(b.Contains(radius));
    REQUIRE_FALSE(a.Contains(thickness));
    // A foreign id never reads another parameter's name or scale.
    REQUIRE(b.GetName(radius).empty());
    REQUIRE(b.GetUnit(radius).empty());
    REQUIRE(std::isnan(b.GetReference(radius)));
    REQUIRE(a.GetName(radius) == "radius");
    REQUIRE(b.GetReference(thickness) == 100.0);

    // An index past a set's end is refused too.
    const ParameterId angle = *b.Add("angle", ParameterQuantity::Angle);
    REQUIRE(b.Contains(angle));
    REQUIRE_FALSE(a.Contains(angle)); // index 1 is past a's end
}

TEST_CASE("a parameter set moves with its ids and is not copied", "[core][sensitivity]")
{
    STATIC_REQUIRE_FALSE(std::is_copy_constructible_v<ParameterSet>);
    STATIC_REQUIRE_FALSE(std::is_copy_assignable_v<ParameterSet>);
    STATIC_REQUIRE(std::is_nothrow_move_constructible_v<ParameterSet>);
    STATIC_REQUIRE(std::is_nothrow_move_assignable_v<ParameterSet>);

    ParameterSet original;
    const ParameterId radius = *original.Add("radius", ParameterQuantity::Length, 10.0);
    ParameterSet moved(std::move(original));
    REQUIRE(moved.Contains(radius));
    REQUIRE(moved.GetReference(radius) == 10.0);
    // The moved-from set is empty, with an identity of its own: its new
    // first parameter is not the moved one.
    REQUIRE(original.GetSize() == 0); // NOLINT(bugprone-use-after-move)
    REQUIRE_FALSE(original.Contains(radius));
    const ParameterId again = *original.Add("radius", ParameterQuantity::Length, 10.0);
    REQUIRE(again != radius);
    REQUIRE_FALSE(moved.Contains(again));

    ParameterSet assigned;
    assigned = std::move(moved);
    REQUIRE(assigned.Contains(radius));
    REQUIRE_FALSE(moved.Contains(radius)); // NOLINT(bugprone-use-after-move)

    // A context keeps the ids of the set it was given.
    const SensitivityContext context(std::move(assigned));
    REQUIRE(context.GetParameters().Contains(radius));
}

TEST_CASE("only Valid and NearEvent carry a derivative", "[core][sensitivity]")
{
    STATIC_REQUIRE(Stereon::HasDerivative(DerivativeStatus::Valid));
    STATIC_REQUIRE(Stereon::HasDerivative(DerivativeStatus::NearEvent));
    STATIC_REQUIRE_FALSE(Stereon::HasDerivative(DerivativeStatus::Nondifferentiable));
    STATIC_REQUIRE_FALSE(Stereon::HasDerivative(DerivativeStatus::NotSupported));
    STATIC_REQUIRE_FALSE(Stereon::HasDerivative(DerivativeStatus::Indeterminate));
    STATIC_REQUIRE(Stereon::ToString(DerivativeStatus::NearEvent) == "NearEvent");
    STATIC_REQUIRE(Stereon::ToString(TrackingConvention::UniqueProjection) == "UniqueProjection");
}

TEST_CASE("a computed result keeps its values and regime", "[core][sensitivity]")
{
    auto result = SensitivityResult::Computed(DerivativeStatus::NearEvent, {1.0, -2.0}, 42);
    REQUIRE(result);
    REQUIRE(result->GetStatus() == DerivativeStatus::NearEvent);
    REQUIRE(result->IsQualified());
    const auto values = result->GetValues();
    REQUIRE(values.has_value());
    REQUIRE(values->size() == 2);
    REQUIRE((*values)[1] == -2.0);
    REQUIRE(result->GetSignature() == 42);

    // Unknown diagnostics are NaN, never a plausible zero.
    REQUIRE(std::isnan(result->GetResolution()));
    REQUIRE(std::isnan(result->GetToleranceCeiling()));
    REQUIRE(std::isnan(result->GetResidual()));
    REQUIRE(std::isnan(result->GetConditioning()));
    result->SetTolerancePolicy(1e-6, 1e-3);
    result->SetDiagnostics(1e-12, 3.0);
    result->AddEntity(7, 3);
    REQUIRE(result->GetResolution() == 1e-6);
    REQUIRE(result->GetToleranceCeiling() == 1e-3);
    REQUIRE(result->GetResidual() == 1e-12);
    REQUIRE(result->GetConditioning() == 3.0);
    REQUIRE(result->GetEntities().size() == 1);
    REQUIRE(result->GetEntities()[0] == std::pair<Stereon::StoreId, std::uint32_t>{7, 3});
}

TEST_CASE("a result without a derivative has no values, not zeros", "[core][sensitivity]")
{
    for (const DerivativeStatus status :
         {DerivativeStatus::Nondifferentiable, DerivativeStatus::NotSupported,
          DerivativeStatus::Indeterminate}) {
        const auto result = SensitivityResult::Without(status, "a kink");
        REQUIRE(result);
        REQUIRE(result->GetStatus() == status);
        REQUIRE_FALSE(result->GetValues().has_value());
        REQUIRE_FALSE(result->IsQualified());
        REQUIRE(result->GetReason() == "a kink");

        // A status without a derivative cannot be given one.
        const auto computed = SensitivityResult::Computed(status, {0.0}, 0);
        REQUIRE_FALSE(computed);
        REQUIRE(computed.error().GetCode() == ErrorCode::InvalidInput);
    }
    for (const DerivativeStatus status : {DerivativeStatus::Valid, DerivativeStatus::NearEvent}) {
        const auto without = SensitivityResult::Without(status, "no values");
        REQUIRE_FALSE(without);
        REQUIRE(without.error().GetCode() == ErrorCode::InvalidInput);
    }
}

TEST_CASE("an unimplemented product reports NotSupported by name", "[core][sensitivity]")
{
    STATIC_REQUIRE(Undeclared.Products.IsEmpty());
    STATIC_REQUIRE_FALSE(Undeclared.Supports(DerivativeProduct::Jvp));
    const SensitivityResult result =
        SensitivityResult::NotSupported(Undeclared, DerivativeProduct::Vjp);
    REQUIRE(result.GetStatus() == DerivativeStatus::NotSupported);
    REQUIRE_FALSE(result.GetValues().has_value());
    REQUIRE(result.GetReason().find("MakeThing") != std::string_view::npos);
    REQUIRE(result.GetReason().find("Vjp") != std::string_view::npos);
}

TEST_CASE("capability declarations are constant and queryable", "[core][sensitivity]")
{
    static constexpr OperationCapability Catalogue[] = {
        Undeclared,
        Stereon::Declare("Fit", {Stereon::Dependency::InputPoints, Stereon::Dependency::Tolerances},
                         {}, "fixed knots", {DerivativeProduct::Jvp}),
    };
    STATIC_REQUIRE(Stereon::FindCapability(Catalogue, "Fit") == &Catalogue[1]);
    STATIC_REQUIRE(Stereon::FindCapability(Catalogue, "Missing") == nullptr);
    STATIC_REQUIRE(Catalogue[1].Supports(DerivativeProduct::Jvp));
    STATIC_REQUIRE_FALSE(Catalogue[1].Supports(DerivativeProduct::DenseJacobian));
    STATIC_REQUIRE(Catalogue[1].Dependencies.Contains(Stereon::Dependency::Tolerances));
    STATIC_REQUIRE_FALSE(Catalogue[1].Dependencies.Contains(Stereon::Dependency::Topology));
    STATIC_REQUIRE(Catalogue[1].Correspondence.IsEmpty());
    STATIC_REQUIRE(DependencySet{} == DependencySet{});
}

TEST_CASE("a tracked sample names its snapshot, entity and convention", "[core][sensitivity]")
{
    const auto face =
        TrackedSample::OnFace(5, 2, Point2{0.25, 0.75}, TrackingConvention::FixedSurfaceCoordinates,
                              SampleOrientation::Reversed);
    REQUIRE(face);
    REQUIRE(face->GetKind() == Stereon::EntityKind::Face);
    REQUIRE(face->GetIndex() == 2);
    REQUIRE(face->GetOrientation() == SampleOrientation::Reversed);
    REQUIRE(face->GetCoordinates().Y == 0.75);
    REQUIRE(face->IsOf(5));
    // The same slot of another snapshot (a compacted or edited one) is not
    // this sample's entity.
    REQUIRE_FALSE(face->IsOf(6));

    const auto edge = TrackedSample::OnEdge(5, 1, 0.5, TrackingConvention::UniqueProjection);
    REQUIRE(edge);
    REQUIRE(edge->GetKind() == Stereon::EntityKind::Edge);
    REQUIRE(edge->GetCoordinates().X == 0.5);
    REQUIRE(edge->GetOrientation() == SampleOrientation::AsStored);
}

TEST_CASE("a tracked sample rejects what cannot name a point", "[core][sensitivity]")
{
    const Point2 uv{0.5, 0.5};
    const auto fixed = TrackingConvention::FixedSurfaceCoordinates;
    REQUIRE_FALSE(TrackedSample::OnFace(0, 0, uv, fixed));
    REQUIRE_FALSE(TrackedSample::OnFace(1, Stereon::FaceId::Invalid, uv, fixed));
    REQUIRE_FALSE(TrackedSample::OnFace(1, 0, Point2{NaN, 0.0}, fixed));
    REQUIRE_FALSE(TrackedSample::OnFace(1, 0, uv, TrackingConvention::FixedEdgeParameter));
    REQUIRE_FALSE(TrackedSample::OnEdge(0, 0, 0.5, TrackingConvention::FixedEdgeParameter));
    REQUIRE_FALSE(TrackedSample::OnEdge(1, 0, Infinity, TrackingConvention::FixedEdgeParameter));
    REQUIRE_FALSE(TrackedSample::OnEdge(1, 0, 0.5, fixed));
}

TEST_CASE("a sensitivity context owns its seeds and cannot be shared", "[core][sensitivity]")
{
    STATIC_REQUIRE_FALSE(std::is_copy_constructible_v<SensitivityContext>);
    STATIC_REQUIRE_FALSE(std::is_copy_assignable_v<SensitivityContext>);
    STATIC_REQUIRE(std::is_nothrow_move_constructible_v<SensitivityContext>);

    SensitivityContext context(TwoParameters());
    REQUIRE(context.GetParameters().GetSize() == 2);

    const std::vector<double> seed{1.0, 0.0};
    REQUIRE(context.SetSeed(seed));
    REQUIRE(context.GetSeed().size() == 2);
    const std::vector<double> shortSeed{1.0};
    REQUIRE(context.SetSeed(shortSeed).error().GetCode() == ErrorCode::InvalidInput);
    const std::vector<double> badSeed{1.0, NaN};
    REQUIRE_FALSE(context.SetSeed(badSeed));
    REQUIRE(context.GetSeed()[0] == 1.0); // the rejected seed changed nothing

    const std::vector<double> adjoint{0.5, 0.5, 0.5};
    REQUIRE(context.SetAdjoint(adjoint));
    REQUIRE(context.GetAdjoint().size() == 3);
    const std::vector<double> badAdjoint{Infinity};
    REQUIRE_FALSE(context.SetAdjoint(badAdjoint));

    REQUIRE(context.Workspace(16).size() == 16);
    REQUIRE(context.Workspace(4).size() == 4);

    SensitivityContext moved(std::move(context));
    REQUIRE(moved.GetSeed().size() == 2);
    moved.Clear();
    REQUIRE(moved.GetSeed().empty());
    REQUIRE(moved.GetAdjoint().empty());
    REQUIRE(moved.GetParameters().GetSize() == 2);
}

namespace {

// Each borrowing accessor as a template over the owner, so that a call the
// owner refuses is a false requirement rather than an error (ADR-0046). The
// audit of 2026-10-30 found GetValues() read through a temporary
// (Computed(...).value().GetValues()) and used after its storage was freed.
template <typename R>
constexpr bool ValuesFromRvalue = requires(R&& r) { std::move(r).GetValues(); };
template <typename R>
constexpr bool ValuesFromLvalue = requires(const R& r) { r.GetValues(); };
template <typename R>
constexpr bool ReasonFromRvalue = requires(R&& r) { std::move(r).GetReason(); };
template <typename R>
constexpr bool EntitiesFromRvalue = requires(R&& r) { std::move(r).GetEntities(); };
template <typename C>
constexpr bool SeedFromRvalue = requires(C&& c) { std::move(c).GetSeed(); };
template <typename C>
constexpr bool AdjointFromRvalue = requires(C&& c) { std::move(c).GetAdjoint(); };
template <typename C>
constexpr bool ParametersFromRvalue = requires(C&& c) { std::move(c).GetParameters(); };
template <typename C>
constexpr bool ParametersFromLvalue = requires(const C& c) { c.GetParameters(); };
template <typename C>
constexpr bool WorkspaceFromRvalue = requires(C&& c) { std::move(c).Workspace(1); };
template <typename C>
constexpr bool WorkspaceFromLvalue = requires(C& c) { c.Workspace(1); };
template <typename S>
constexpr bool NameFromRvalue = requires(S&& s, ParameterId id) { std::move(s).GetName(id); };
template <typename S>
constexpr bool UnitFromRvalue = requires(S&& s, ParameterId id) { std::move(s).GetUnit(id); };

} // namespace

TEST_CASE("sensitivity borrows are taken from named owners only", "[core][sensitivity][borrowing]")
{
    STATIC_REQUIRE(ValuesFromLvalue<SensitivityResult>);
    STATIC_REQUIRE_FALSE(ValuesFromRvalue<SensitivityResult>);
    STATIC_REQUIRE_FALSE(ReasonFromRvalue<SensitivityResult>);
    STATIC_REQUIRE_FALSE(EntitiesFromRvalue<SensitivityResult>);
    STATIC_REQUIRE(ParametersFromLvalue<SensitivityContext>);
    STATIC_REQUIRE_FALSE(ParametersFromRvalue<SensitivityContext>);
    STATIC_REQUIRE_FALSE(SeedFromRvalue<SensitivityContext>);
    STATIC_REQUIRE_FALSE(AdjointFromRvalue<SensitivityContext>);
    STATIC_REQUIRE(WorkspaceFromLvalue<SensitivityContext>);
    STATIC_REQUIRE_FALSE(WorkspaceFromRvalue<SensitivityContext>);
    STATIC_REQUIRE_FALSE(NameFromRvalue<ParameterSet>);
    // The unit is a literal, which outlives every set.
    STATIC_REQUIRE(UnitFromRvalue<ParameterSet>);

    // The named result's borrow is the result's own storage.
    const auto result = SensitivityResult::Computed(DerivativeStatus::Valid, {1.0, 2.0}, 7);
    REQUIRE(result.has_value());
    const SensitivityResult& named = *result;
    REQUIRE(named.GetValues()->size() == 2);
    REQUIRE((*named.GetValues())[1] == 2.0);
}

TEST_CASE("a derivative that is not finite is never published as computed", "[core][sensitivity]")
{
    constexpr double inf = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (const DerivativeStatus status : {DerivativeStatus::Valid, DerivativeStatus::NearEvent}) {
        for (const double bad : {inf, -inf, nan}) {
            const auto refused = SensitivityResult::Computed(status, {1.0, bad}, 0);
            REQUIRE_FALSE(refused.has_value());
            REQUIRE(refused.error().GetCode() == ErrorCode::InvalidInput);
        }
    }
    // What an operation reports instead: no derivative, and why.
    const auto indeterminate =
        SensitivityResult::Without(DerivativeStatus::Indeterminate, "overflowed");
    REQUIRE(indeterminate.has_value());
    const SensitivityResult& named = *indeterminate;
    REQUIRE_FALSE(named.GetValues().has_value());
}
