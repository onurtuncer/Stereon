// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// Stereon.Curve: curves in space and in a surface's parameter plane
/// (ADR-0003): the CurveGeometry and Curve2dGeometry concepts, the built-in
/// Line and Circle, and the type-erased Curve and Curve2d.
export module Stereon.Curve;

export import Stereon.Robust;

export import :Circle;
export import :Concepts;
export import :Erased;
export import :Line;
export import :Visit;
