// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Private shared clipping entry point, included in the module purview
// after :Dispatch, with <vector> and <string_view> in the global fragment.
#pragma once

namespace Stereon::Detail {

/// Complete planar intersections with the box's edge lines partition the
/// bounded range. Returns ordered closed stretches, including singleton
/// contacts, or fails if a piece cannot be classified reliably.
[[nodiscard]] Result<std::vector<Interval>> ClipPcurveToBox(const Curve2d& pcurve,
                                                            const Interval& range, const Box2& box,
                                                            const Context& ctx,
                                                            std::string_view handler);

} // namespace Stereon::Detail
