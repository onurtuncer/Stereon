// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// Stereon.Intersect: curve-curve, curve-surface and surface-surface
/// intersection (ADR-0010). So far the type-pair dispatch of ADR-0003 with
/// stub handlers; the analytic handlers follow in Phase 2.
export module Stereon.Intersect;

export import Stereon.Surface;

export import :Capabilities;
export import :Dispatch;
