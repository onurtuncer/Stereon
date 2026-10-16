// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Stereon's points and B-splines as OCCT's, for the tests that compare the
// two kernels. Include after `import Stereon.Surface;` and the OCCT headers
// for gp_Pnt, Geom_BSplineCurve, Geom_BSplineSurface and NCollection arrays.

#pragma once

#include <cstddef>
#include <vector>

namespace StereonTest {

inline gp_Pnt ToOcct(const Stereon::Point3& p)
{
    return {p.X, p.Y, p.Z};
}

/// Distinct values and their multiplicities, 1-based, as OCCT stores knots.
struct OcctKnots {
    NCollection_Array1<double> Values;
    NCollection_Array1<int> Multiplicities;
};

inline OcctKnots ToOcctKnots(const std::vector<double>& knots)
{
    std::vector<double> values;
    std::vector<int> counts;
    for (const double u : knots) {
        if (values.empty() || values.back() != u) {
            values.push_back(u);
            counts.push_back(0);
        }
        ++counts.back();
    }
    const int n = static_cast<int>(values.size());
    OcctKnots result{NCollection_Array1<double>(1, n), NCollection_Array1<int>(1, n)};
    for (int i = 0; i < n; ++i) {
        result.Values(i + 1) = values[static_cast<std::size_t>(i)];
        result.Multiplicities(i + 1) = counts[static_cast<std::size_t>(i)];
    }
    return result;
}

/// Stereon's B-spline curve as OCCT's: distinct knots with multiplicities,
/// and 1-based arrays.
inline occ::handle<Geom_BSplineCurve> ToOcct(const Stereon::BSplineCurve& c)
{
    const int poleCount = static_cast<int>(c.Poles.size());
    NCollection_Array1<gp_Pnt> poles(1, poleCount);
    NCollection_Array1<double> weights(1, poleCount);
    for (int i = 0; i < poleCount; ++i) {
        const auto index = static_cast<std::size_t>(i);
        poles(i + 1) = ToOcct(c.Poles[index]);
        weights(i + 1) = c.Weights.empty() ? 1.0 : c.Weights[index];
    }
    const OcctKnots knots = ToOcctKnots(c.Knots);
    return new Geom_BSplineCurve(poles, weights, knots.Values, knots.Multiplicities, c.Degree);
}

/// Stereon's B-spline surface as OCCT's; Stereon's poles are row-major with
/// u along the rows.
inline occ::handle<Geom_BSplineSurface> ToOcct(const Stereon::BSplineSurface& s)
{
    const auto countU = static_cast<int>(Stereon::PoleCountU(s));
    const auto countV = static_cast<int>(Stereon::PoleCountV(s));
    NCollection_Array2<gp_Pnt> poles(1, countU, 1, countV);
    NCollection_Array2<double> weights(1, countU, 1, countV);
    for (int i = 0; i < countU; ++i) {
        for (int j = 0; j < countV; ++j) {
            const auto k = static_cast<std::size_t>(i * countV + j);
            poles(i + 1, j + 1) = ToOcct(s.Poles[k]);
            weights(i + 1, j + 1) = s.Weights.empty() ? 1.0 : s.Weights[k];
        }
    }
    const OcctKnots u = ToOcctKnots(s.KnotsU);
    const OcctKnots v = ToOcctKnots(s.KnotsV);
    return new Geom_BSplineSurface(poles, weights, u.Values, v.Values, u.Multiplicities,
                                   v.Multiplicities, s.DegreeU, s.DegreeV);
}

} // namespace StereonTest
