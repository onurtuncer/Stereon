// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <cstdlib>
#include <set>
#include <span>
#include <string>
#include <string_view>
#include <vector>

import Stereon.Boolean;

using Stereon::DerivativeProduct;
using Stereon::DerivativeStatus;
using Stereon::OperationCapability;

// The catalogues (ADR-0021, section 5) as code sees them. Which public
// functions must be declared is tools/check_capabilities.py's business; this
// checks what the declarations say.

namespace {

constexpr DerivativeProduct Products[] = {DerivativeProduct::Jvp, DerivativeProduct::Vjp,
                                          DerivativeProduct::DenseJacobian};

/// The declaration of `name`; the test fails when there is none. Checking
/// the pointer in the caller lets GCC 16.2 at -O3 see a null path and warn
/// (-Wnull-dereference), as Bodies.hpp's Deref explains.
const OperationCapability& Declared(const std::span<const OperationCapability> catalogue,
                                    const std::string_view name)
{
    const OperationCapability* const capability = Stereon::FindCapability(catalogue, name);
    if (capability == nullptr) {
        FAIL("no declaration of " << name);
        std::abort();
    }
    return *capability;
}

} // namespace

TEST_CASE("each library's catalogue declares each operation once", "[boolean][sensitivity]")
{
    for (const std::span<const OperationCapability> catalogue :
         {Stereon::CurveCapabilities(), Stereon::SurfaceCapabilities(),
          Stereon::IntersectCapabilities(), Stereon::TopoCapabilities(),
          Stereon::BuildCapabilities(), Stereon::BooleanCapabilities()}) {
        REQUIRE_FALSE(catalogue.empty());
        std::set<std::string_view> names;
        for (const OperationCapability& capability : catalogue) {
            INFO(capability.Operation);
            REQUIRE(names.insert(capability.Operation).second);
            REQUIRE(Stereon::FindCapability(catalogue, capability.Operation) == &capability);
        }
    }
}

TEST_CASE("every declaration states its dependencies and assumptions", "[boolean][sensitivity]")
{
    const std::vector<OperationCapability> all = Stereon::KernelCapabilities();
    REQUIRE(all.size() ==
            Stereon::CurveCapabilities().size() + Stereon::SurfaceCapabilities().size() +
                Stereon::IntersectCapabilities().size() + Stereon::TopoCapabilities().size() +
                Stereon::BuildCapabilities().size() + Stereon::BooleanCapabilities().size());
    for (const OperationCapability& capability : all) {
        INFO(capability.Operation);
        REQUIRE_FALSE(capability.Operation.empty());
        REQUIRE_FALSE(capability.Dependencies.IsEmpty());
        REQUIRE_FALSE(capability.Assumptions.empty());
    }
}

TEST_CASE("only the curve's Eval implements products; the rest answer NotSupported",
          "[boolean][sensitivity]")
{
    // The first derivative slice (ADR-0021, stage 2): Eval of a B-spline
    // curve has JVP and VJP. An evaluation asked of any other operation is
    // NotSupported and carries no values, never a zero. A new slice changes
    // this test deliberately.
    const OperationCapability& eval = Declared(Stereon::CurveCapabilities(), "Eval");
    REQUIRE(eval.Supports(DerivativeProduct::Jvp));
    REQUIRE(eval.Supports(DerivativeProduct::Vjp));
    REQUIRE_FALSE(eval.Supports(DerivativeProduct::DenseJacobian));
    std::size_t implemented = 0;
    for (const OperationCapability& capability : Stereon::KernelCapabilities()) {
        INFO(capability.Operation);
        if (!capability.Products.IsEmpty()) {
            // KernelCapabilities copies the declarations: the curve's Eval is
            // the one whose products and correspondence are eval's.
            REQUIRE(capability.Operation == eval.Operation);
            REQUIRE(capability.Products == eval.Products);
            REQUIRE(capability.Correspondence == eval.Correspondence);
            ++implemented;
            continue;
        }
        for (const DerivativeProduct product : Products) {
            REQUIRE_FALSE(capability.Supports(product));
            const auto result = Stereon::SensitivityResult::NotSupported(capability, product);
            REQUIRE(result.GetStatus() == DerivativeStatus::NotSupported);
            REQUIRE_FALSE(result.GetValues().has_value());
            REQUIRE(result.GetReason().find(capability.Operation) != std::string_view::npos);
        }
    }
    REQUIRE(implemented == 1);
}

TEST_CASE("operations on points declare how their samples correspond", "[boolean][sensitivity]")
{
    // Projection is the implicit construction of ADR-0021, section 6: its
    // samples are tracked by the unique projection, and only it may say so.
    for (const std::string_view name : {"Project", "ProjectAll"}) {
        for (const std::span<const OperationCapability> catalogue :
             {Stereon::CurveCapabilities(), Stereon::SurfaceCapabilities()}) {
            const OperationCapability& capability = Declared(catalogue, name);
            REQUIRE(
                capability.Correspondence.Contains(Stereon::TrackingConvention::UniqueProjection));
            REQUIRE(capability.Dependencies.Contains(Stereon::Dependency::InputPoints));
        }
    }
    REQUIRE(Declared(Stereon::TopoCapabilities(), "Freeze")
                .Dependencies.Contains(Stereon::Dependency::Topology));
}

TEST_CASE("a Boolean's declaration names the trims that move with its operands",
          "[boolean][sensitivity]")
{
    // The intersection edges a Boolean makes move when its operands do, so
    // their points are tracked on moving trims, and a tolerance change is a
    // branch event (ADR-0021, section 8).
    const OperationCapability& boolean = Declared(Stereon::BooleanCapabilities(), "Boolean");
    REQUIRE(boolean.Correspondence.Contains(Stereon::TrackingConvention::MovingTrimCoordinates));
    REQUIRE(boolean.Dependencies.Contains(Stereon::Dependency::InputGeometry));
    REQUIRE(boolean.Dependencies.Contains(Stereon::Dependency::Tolerances));
}
