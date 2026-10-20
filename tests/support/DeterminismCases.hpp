// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The inputs of the determinism checks (ADR-0007): a rational cubic curve
// with several spans, a rational bicubic surface, and points from a
// counter-based generator, so that every run, process and thread count sees
// the same inputs; and the results as their bit patterns, which the checks
// compare. Shared by tests/surface/DeterminismTests.cpp and the nightly
// executable tests/determinism/DeterminismRun.cpp.
//
// Include after `import Stereon.Surface;`, with <bit>, <cstdint>, <span>
// and <vector> included before it (GCC 16.2 rejects a standard header
// first included after an import of a module that uses it).

#pragma once

namespace DeterminismTest {

/// SplitMix64: the k-th value of a counter-based stream, the same in any
/// order and on any thread (ADR-0007, seeded counter-based randomness).
inline std::uint64_t Mix(std::uint64_t k)
{
    k += 0x9E3779B97F4A7C15ULL;
    k = (k ^ (k >> 30)) * 0xBF58476D1CE4E5B9ULL;
    k = (k ^ (k >> 27)) * 0x94D049BB133111EBULL;
    return k ^ (k >> 31);
}

/// A double in [low, high) from the k-th value of the stream.
inline double Uniform(const std::uint64_t k, const double low, const double high)
{
    const double unit = static_cast<double>(Mix(k) >> 11) * 0x1p-53;
    return low + (high - low) * unit;
}

/// `count` points around the curve and the surface, in [-2, 5]^3.
inline std::vector<Stereon::Point3> Points(const std::size_t count)
{
    std::vector<Stereon::Point3> points;
    points.reserve(count);
    for (std::uint64_t i = 0; i < count; ++i) {
        points.push_back({Uniform(3 * i, -2.0, 5.0), Uniform(3 * i + 1, -2.0, 5.0),
                          Uniform(3 * i + 2, -2.0, 5.0)});
    }
    return points;
}

/// A rational cubic with four spans that winds through space.
inline Stereon::BSplineCurve Curve()
{
    return *Stereon::MakeBSplineCurve(3, {0.0, 0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0, 1.0},
                                      {{0.0, 0.0, 0.0},
                                       {1.0, 2.0, 0.5},
                                       {2.0, -1.0, 1.0},
                                       {3.0, 1.0, 2.5},
                                       {4.0, 3.0, 0.0},
                                       {2.5, 4.0, -1.0},
                                       {1.0, 3.0, 1.5}},
                                      {1.0, 0.5, 2.0, 0.75, 1.5, 1.0, 1.0});
}

/// A rational bicubic, 6 x 6 poles over [0, 3]^2, with a wave in z.
inline Stereon::BSplineSurface Surface()
{
    const std::vector<double> knots{0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0, 3.0};
    std::vector<Stereon::Point3> poles;
    std::vector<double> weights;
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            const double z = ((i + 2 * j) % 3 == 0 ? 1.0 : -0.5) + 0.25 * i;
            poles.push_back({0.6 * i, 0.6 * j, z});
            weights.push_back(1.0 + 0.125 * ((i * 7 + j * 3) % 5));
        }
    }
    return *Stereon::MakeBSplineSurface(3, 3, knots, knots, std::move(poles), std::move(weights));
}

/// The bit patterns of the curve parameters, then of the surface parameters:
/// the bytes the checks compare. Comparing doubles with == would let -0.0 and
/// 0.0 pass as equal.
inline std::vector<std::uint64_t> Bits(const std::vector<double>& curve,
                                       const std::vector<Stereon::Point2>& surface)
{
    std::vector<std::uint64_t> bits;
    bits.reserve(curve.size() + 2 * surface.size());
    for (const double t : curve) {
        bits.push_back(std::bit_cast<std::uint64_t>(t));
    }
    for (const Stereon::Point2& uv : surface) {
        bits.push_back(std::bit_cast<std::uint64_t>(uv.X));
        bits.push_back(std::bit_cast<std::uint64_t>(uv.Y));
    }
    return bits;
}

/// Projects every point onto both and returns the bits, under `context`.
inline std::vector<std::uint64_t> Run(const Stereon::CurveProjector& curve,
                                      const Stereon::SurfaceProjector& surface,
                                      std::span<const Stereon::Point3> points,
                                      const Stereon::Context& context)
{
    return Bits(ProjectAll(curve, points, context), ProjectAll(surface, points, context));
}

} // namespace DeterminismTest
