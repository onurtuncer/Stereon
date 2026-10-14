// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

export module Stereon.Core:Allowances;

/// Purpose-specific numerical allowances (ADR-0028). Each is its own type,
/// so that an approximation budget cannot be passed where a modelling
/// tolerance or a solver target is expected, and every quantity states its
/// units. The modelling tolerances of ADR-0001 stay with the Context; the
/// types here are the per-operation requests that ADR-0028 separates from
/// them.

export namespace Stereon {

/// How far a fitted or cached representation may deviate from the data or
/// the geometry it represents: ADR-0028's approximation budget, in
/// millimetres. It is a request; an operation that honours it reports the
/// deviation it achieved alongside the result, and never widens a modelling
/// tolerance because of it.
struct ApproximationBudget {
    double Deviation = 0.0;

    constexpr bool operator==(const ApproximationBudget&) const = default;
};

} // namespace Stereon
