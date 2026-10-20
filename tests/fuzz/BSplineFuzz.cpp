// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Fuzz target: the checked B-spline factory (ADR-0006, ADR-0016) and the
// evaluators on what it accepts. The bytes give a degree, a pole count, then
// knots, poles and weights as doubles; whatever the factory accepts must
// evaluate, differentiate, bound and project at parameters across its domain
// without a contract violation, since the factory is the boundary that turns
// outside data into geometry the evaluators may trust.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <vector>

import Stereon.Curve;

namespace {

/// Reads the input front to back; past its end, zeros.
class Reader {
public:
    Reader(const std::uint8_t* data, const std::size_t size)
        : m_Data(data)
        , m_Size(size)
    {}
    std::uint8_t Byte() { return m_Offset < m_Size ? m_Data[m_Offset++] : 0; }
    double Double()
    {
        std::uint8_t raw[8] = {};
        for (std::uint8_t& b : raw) {
            b = Byte();
        }
        double x = 0.0;
        std::memcpy(&x, raw, sizeof(x));
        return x;
    }

private:
    const std::uint8_t* m_Data;
    std::size_t m_Size;
    std::size_t m_Offset = 0;
};

/// Evaluates through the erased wrapper, whose parameter checks are compiled
/// as plain branches (GCC 16.2 drops contracts in imported templates).
void Exercise(const Stereon::BSplineCurve& bspline)
{
    const Stereon::Curve curve(bspline);
    const Stereon::Interval domain = Domain(curve);
    for (int k = 0; k <= 8; ++k) {
        // A convex combination, finite for any finite domain, where the
        // width alone could overflow.
        const double s = k / 8.0;
        const double t = std::clamp(domain.GetLower() * (1.0 - s) + domain.GetUpper() * s,
                                    domain.GetLower(), domain.GetUpper());
        static_cast<void>(Eval(curve, t));
        static_cast<void>(Deriv(curve, t, 1));
        static_cast<void>(Deriv(curve, t, 2));
    }
    static_cast<void>(BoundingBox(curve));
    static_cast<void>(Project(curve, Stereon::Point3{0.5, -1.0, 2.0}));
}

} // namespace

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, const std::size_t size)
{
    Reader in(data, size);
    const int degree = in.Byte() % 8;
    const auto poles = static_cast<std::size_t>(1 + in.Byte() % 24);
    const bool rational = (in.Byte() & 1U) != 0;
    std::vector<double> knots(poles + static_cast<std::size_t>(degree) + 1);
    for (double& k : knots) {
        k = in.Double();
    }
    // Sorting often turns noise into a valid knot vector, so the evaluators
    // are reached; unsorted input stays possible through the first byte.
    if ((data != nullptr && size > 0 && (data[0] & 0x80U) == 0)) {
        std::ranges::sort(knots, [](const double a, const double b) {
            return std::isnan(b) ? !std::isnan(a) : a < b;
        });
    }
    std::vector<Stereon::Point3> points(poles);
    for (Stereon::Point3& p : points) {
        p = Stereon::Point3{in.Double(), in.Double(), in.Double()};
    }
    std::vector<double> weights;
    if (rational) {
        weights.resize(poles);
        for (double& w : weights) {
            w = in.Double();
        }
    }
    const auto curve = Stereon::MakeBSplineCurve(degree, knots, points, weights);
    if (curve) {
        Exercise(*curve);
    }
    return 0;
}
