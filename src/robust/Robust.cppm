// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// Stereon.Robust: exact predicates, expansion arithmetic and interval
/// filters (ADR-0002). Every topological decision in topo, boolean and mesh
/// goes through a predicate from here.
export module Stereon.Robust;

export import Stereon.Core;

export import :Domain;
export import :Expansion;
export import :Filter;
export import :Predicates;
export import :Sign;
