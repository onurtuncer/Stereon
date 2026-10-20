// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// Stereon.Topo: the B-rep store (ADR-0004): entity records, ShapeBuilder,
/// Freeze() with its structural validator, the frozen Shape and its
/// read-only ranges (ADR-0036 to ADR-0038).
export module Stereon.Topo;

export import Stereon.Intersect;
export import :Builder;
export import :Capabilities;
export import :Operation;
export import :Records;
export import :Report;
export import :Shape;
