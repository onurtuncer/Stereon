// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

module Stereon.Core;

import :Error;
import :Handle;
import :Sensitivity;
import :Vector;

namespace Stereon {

namespace {

[[nodiscard]] bool AllFinite(const std::span<const double> values)
{
    return std::ranges::all_of(values, [](const double x) { return std::isfinite(x); });
}

[[nodiscard]] std::string_view ProductName(const DerivativeProduct product) noexcept
{
    switch (product) {
    case DerivativeProduct::Jvp:
        return "Jvp";
    case DerivativeProduct::Vjp:
        return "Vjp";
    case DerivativeProduct::DenseJacobian:
        return "DenseJacobian";
    }
    return "product";
}

} // namespace

// --- ParameterSet -----------------------------------------------------------------------

Result<ParameterId> ParameterSet::Add(const std::string_view name, const ParameterQuantity quantity,
                                      const double reference)
{
    if (name.empty()) {
        return Fail(ErrorCode::InvalidInput, "ParameterSet::Add: a parameter needs a name");
    }
    if (Find(name).IsValid()) {
        return Fail(ErrorCode::InvalidInput, "ParameterSet::Add: a parameter named \"" +
                                                 std::string(name) + "\" is already in the set");
    }
    if (!std::isfinite(reference) || !(reference > 0.0)) {
        return Fail(ErrorCode::InvalidInput,
                    "ParameterSet::Add: the reference scale must be finite and positive");
    }
    const auto index = static_cast<std::uint32_t>(m_Names.size());
    m_Names.emplace_back(name);
    m_Quantities.push_back(quantity);
    m_References.push_back(reference);
    return ParameterId(index, quantity);
}

bool ParameterSet::Contains(const ParameterId id) const noexcept
{
    return id.IsValid() && id.GetIndex() < m_Names.size() &&
           m_Quantities[id.GetIndex()] == id.GetQuantity();
}

std::string_view ParameterSet::GetName(const ParameterId id) const
{
    return m_Names.at(id.GetIndex());
}

std::string_view ParameterSet::GetUnit(const ParameterId id) const
{
    return UnitOf(m_Quantities.at(id.GetIndex()));
}

double ParameterSet::GetReference(const ParameterId id) const
{
    return m_References.at(id.GetIndex());
}

ParameterId ParameterSet::Find(const std::string_view name) const noexcept
{
    for (std::size_t i = 0; i < m_Names.size(); ++i) {
        if (m_Names[i] == name) {
            return ParameterId(static_cast<std::uint32_t>(i), m_Quantities[i]);
        }
    }
    return {};
}

// --- TrackedSample ----------------------------------------------------------------------

Result<TrackedSample> TrackedSample::OnFace(const StoreId snapshot, const std::uint32_t face,
                                            const Point2& uv, const TrackingConvention convention,
                                            const SampleOrientation orientation)
{
    if (snapshot == 0 || face == FaceId::Invalid) {
        return Fail(ErrorCode::InvalidInput,
                    "TrackedSample::OnFace: the sample needs a snapshot and a face of it");
    }
    if (!std::isfinite(uv.X) || !std::isfinite(uv.Y)) {
        return Fail(ErrorCode::InvalidInput,
                    "TrackedSample::OnFace: the coordinates must be finite");
    }
    if (convention == TrackingConvention::FixedEdgeParameter) {
        return Fail(
            ErrorCode::InvalidInput,
            "TrackedSample::OnFace: a fixed edge parameter does not track a point of a face");
    }
    TrackedSample sample;
    sample.m_Kind = EntityKind::Face;
    sample.m_Snapshot = snapshot;
    sample.m_Index = face;
    sample.m_Orientation = orientation;
    sample.m_Convention = convention;
    sample.m_Coordinates = uv;
    return sample;
}

Result<TrackedSample> TrackedSample::OnEdge(const StoreId snapshot, const std::uint32_t edge,
                                            const double t, const TrackingConvention convention,
                                            const SampleOrientation orientation)
{
    if (snapshot == 0 || edge == EdgeId::Invalid) {
        return Fail(ErrorCode::InvalidInput,
                    "TrackedSample::OnEdge: the sample needs a snapshot and an edge of it");
    }
    if (!std::isfinite(t)) {
        return Fail(ErrorCode::InvalidInput, "TrackedSample::OnEdge: the parameter must be finite");
    }
    if (convention == TrackingConvention::FixedSurfaceCoordinates) {
        return Fail(ErrorCode::InvalidInput, "TrackedSample::OnEdge: fixed surface coordinates do "
                                             "not track a point of an edge");
    }
    TrackedSample sample;
    sample.m_Kind = EntityKind::Edge;
    sample.m_Snapshot = snapshot;
    sample.m_Index = edge;
    sample.m_Orientation = orientation;
    sample.m_Convention = convention;
    sample.m_Coordinates = Point2{t, 0.0};
    return sample;
}

// --- SensitivityContext -----------------------------------------------------------------

SensitivityContext::SensitivityContext(ParameterSet parameters)
    : m_Parameters(std::move(parameters))
{}

Result<void> SensitivityContext::SetSeed(const std::span<const double> seed)
{
    if (seed.size() != m_Parameters.GetSize()) {
        return Fail(ErrorCode::InvalidInput, "SensitivityContext::SetSeed: the seed has " +
                                                 std::to_string(seed.size()) + " components for " +
                                                 std::to_string(m_Parameters.GetSize()) +
                                                 " parameters");
    }
    if (!AllFinite(seed)) {
        return Fail(ErrorCode::InvalidInput,
                    "SensitivityContext::SetSeed: every component must be finite");
    }
    m_Seed.assign(seed.begin(), seed.end());
    return {};
}

Result<void> SensitivityContext::SetAdjoint(const std::span<const double> adjoint)
{
    if (!AllFinite(adjoint)) {
        return Fail(ErrorCode::InvalidInput,
                    "SensitivityContext::SetAdjoint: every component must be finite");
    }
    m_Adjoint.assign(adjoint.begin(), adjoint.end());
    return {};
}

std::span<double> SensitivityContext::Workspace(const std::size_t size)
{
    if (m_Workspace.size() < size) {
        m_Workspace.resize(size);
    }
    return std::span<double>{m_Workspace}.first(size);
}

void SensitivityContext::Clear() noexcept
{
    m_Seed.clear();
    m_Adjoint.clear();
    m_Workspace.clear();
}

// --- SensitivityResult ------------------------------------------------------------------

Result<SensitivityResult> SensitivityResult::Computed(const DerivativeStatus status,
                                                      std::vector<double> values,
                                                      const std::uint64_t signature)
{
    if (!HasDerivative(status)) {
        return Fail(ErrorCode::InvalidInput, "SensitivityResult::Computed: a " +
                                                 std::string(ToString(status)) +
                                                 " result carries no derivative");
    }
    SensitivityResult result;
    result.m_Status = status;
    result.m_Values = std::move(values);
    result.m_Signature = signature;
    return result;
}

Result<SensitivityResult> SensitivityResult::Without(const DerivativeStatus status,
                                                     std::string reason)
{
    if (HasDerivative(status)) {
        return Fail(ErrorCode::InvalidInput, "SensitivityResult::Without: a " +
                                                 std::string(ToString(status)) +
                                                 " result needs its derivative");
    }
    SensitivityResult result;
    result.m_Status = status;
    result.m_Reason = std::move(reason);
    return result;
}

SensitivityResult SensitivityResult::NotSupported(const OperationCapability& capability,
                                                  const DerivativeProduct product)
{
    SensitivityResult result;
    result.m_Status = DerivativeStatus::NotSupported;
    result.m_Reason = std::string(capability.Operation) + " implements no " +
                      std::string(ProductName(product)) + " (ADR-0021: declared, not computed)";
    return result;
}

std::optional<std::span<const double>> SensitivityResult::GetValues() const noexcept
{
    if (!HasDerivative(m_Status)) {
        return std::nullopt;
    }
    return std::span<const double>{m_Values};
}

void SensitivityResult::AddEntity(const StoreId snapshot, const std::uint32_t index)
{
    m_Entities.emplace_back(snapshot, index);
}

void SensitivityResult::SetTolerancePolicy(const double resolution,
                                           const double toleranceCeiling) noexcept
{
    m_Resolution = resolution;
    m_ToleranceCeiling = toleranceCeiling;
}

void SensitivityResult::SetDiagnostics(const double residual, const double conditioning) noexcept
{
    m_Residual = residual;
    m_Conditioning = conditioning;
}

} // namespace Stereon
