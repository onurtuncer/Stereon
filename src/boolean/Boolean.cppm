// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// Stereon.Boolean: union, cut and common (ADR-0013), Boolean() in
/// :Assemble (M8), over the pipeline of the Boolean prototype
/// (docs/planning/boolean-prototype-plan.md): face clipping and imprinting
/// (:Clip and :Imprint, M6), the two-body driver (:Bodies, M6d) and the
/// classification with its winding-number fallback (:Classify and :Winding,
/// M7).
export module Stereon.Boolean;

export import Stereon.Build;
export import Stereon.Intersect;

export import :Assemble;
export import :Bodies;
export import :Capabilities;
export import :Classify;
export import :Clip;
export import :Imprint;
export import :Winding;
