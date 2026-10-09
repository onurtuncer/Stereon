// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// RapidCheck with its Catch2 adapter (rc::prop). Include this instead of the
/// RapidCheck headers: the adapter only picks Catch2 v3 when the v3 macros
/// header comes first, which include sorting would otherwise undo.

#pragma once

// clang-format off
#include <catch2/catch_test_macros.hpp>
#include <rapidcheck.h>
#include <rapidcheck/catch.h>
// clang-format on
