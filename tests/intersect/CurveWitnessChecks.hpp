// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0

// Include after Stereon.Intersect, with <algorithm>, <cmath> and <limits>
// included before the import.
#pragma once

namespace IntersectTest {

/// A meeting witness must lie in both domains and within the model
/// resolution of both evaluated supports. The extra allowance covers only
/// floating-point evaluation/comparison, not a larger model tolerance.
[[nodiscard]] inline bool ValidCurveWitness(const Stereon::CurveCurvePoint& p,
                                            const Stereon::Curve& a, const Stereon::Interval& da,
                                            const Stereon::Curve& b, const Stereon::Interval& db,
                                            const Stereon::Context& policy)
{
    if (!std::isfinite(p.ParameterA) || !std::isfinite(p.ParameterB) ||
        !Contains(da, p.ParameterA) || !Contains(db, p.ParameterB)) {
        return false;
    }
    const Stereon::Point3 onA = Eval(a, p.ParameterA);
    const Stereon::Point3 onB = Eval(b, p.ParameterB);
    const double scale =
        std::max({1.0, Length(onA - Stereon::Point3{}), Length(onB - Stereon::Point3{})});
    const double rounding = 64.0 * std::numeric_limits<double>::epsilon() * scale;
    const double reach = policy.GetResolution() + rounding;
    return Distance(p.Position, onA) <= reach && Distance(p.Position, onB) <= reach &&
           Distance(onA, onB) <= reach;
}

} // namespace IntersectTest
