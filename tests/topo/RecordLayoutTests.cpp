// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The sizes behind B0.4's memory per face (ADR-0004): the records of a
// release build, where handles are four bytes, and the erased geometry
// wrappers the pools hold. A field added or reordered so that a record
// grows fails here, rather than as a benchmark regression found later.

#include <catch2/catch_test_macros.hpp>

#include <cstddef>

import Stereon.Topo;

namespace Detail = Stereon::Detail;

TEST_CASE("a release build's topology records are as small as their fields allow", "[topo][memory]")
{
    STATIC_REQUIRE(sizeof(Detail::SlotState) == 4);
#if !STEREON_CHECKED_HANDLES
    {
        STATIC_REQUIRE(sizeof(Stereon::EdgeId) == 4);
        // Slot and incidence head, position, tolerance.
        STATIC_REQUIRE(sizeof(Detail::VertexRecord) == 40);
        // Slot and six handles, two geometry refs, kind and cache status, the
        // range, the tolerance and the cache bound.
        STATIC_REQUIRE(sizeof(Detail::EdgeRecord) == 72);
        // Slot, five handles, the pcurve and the sense.
        STATIC_REQUIRE(sizeof(Detail::CoedgeRecord) == 32);
        STATIC_REQUIRE(sizeof(Detail::LoopRecord) == 16);
        STATIC_REQUIRE(sizeof(Detail::FaceRecord) == 28);
    }
#endif
}

TEST_CASE("the geometry wrappers hold their largest built-in type and no more", "[topo][memory]")
{
    STATIC_REQUIRE(sizeof(Stereon::Curve) == Stereon::CurveBufferSize + sizeof(void*));
    STATIC_REQUIRE(sizeof(Stereon::Curve2d) == Stereon::Curve2dBufferSize + sizeof(void*));
    STATIC_REQUIRE(sizeof(Stereon::Surface) == Stereon::SurfaceBufferSize + sizeof(void*));
    STATIC_REQUIRE(Detail::FitsInline<Stereon::Line2d, Stereon::Curve2dBufferSize>);
    STATIC_REQUIRE(Detail::FitsInline<Stereon::Ellipse, Stereon::CurveBufferSize>);
    STATIC_REQUIRE(Detail::FitsInline<Stereon::Torus, Stereon::SurfaceBufferSize>);
}
