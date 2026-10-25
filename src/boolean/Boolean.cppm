// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// Stereon.Boolean: union, cut and common (ADR-0013). So far the face
/// clipping and imprinting of the Boolean prototype (:Clip and :Imprint,
/// M6); the pipeline of docs/planning/boolean-prototype-plan.md follows
/// milestone by milestone.
export module Stereon.Boolean;

export import Stereon.Build;
export import Stereon.Intersect;

export import :Bodies;
export import :Classify;
export import :Clip;
export import :Imprint;
