// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The cases of corpus/curve-curve/ (tools/generate_curve_corpus.py), read
// as kernel curves. Shared by the kernel's own runner and the OCCT
// comparison. Include Json.hpp, <algorithm>, <filesystem>, <string> and
// <vector> before `import Stereon.Intersect;`, and this after it.

#pragma once

namespace StereonTest {

using CurveJson = JsonValue;

inline Stereon::Point3 PointFrom(const CurveJson& v)
{
    return {v[0].AsNumber(), v[1].AsNumber(), v[2].AsNumber()};
}

inline Stereon::Vec3 VectorFrom(const CurveJson& v)
{
    return {v[0].AsNumber(), v[1].AsNumber(), v[2].AsNumber()};
}

inline Stereon::Interval DomainFrom(const CurveJson& c)
{
    return Stereon::Interval{c["domain"][0].AsNumber(), c["domain"][1].AsNumber()};
}

inline Stereon::BSplineCurve SplineFrom(const CurveJson& c)
{
    std::vector<double> knots;
    for (const CurveJson& k : c["knots"].AsArray()) {
        knots.push_back(k.AsNumber());
    }
    std::vector<Stereon::Point3> poles;
    for (const CurveJson& p : c["poles"].AsArray()) {
        poles.push_back(PointFrom(p));
    }
    return Stereon::BSplineCurve{
        static_cast<int>(c["degree"].AsNumber()), std::move(knots), std::move(poles), {}};
}

inline Stereon::Curve CurveFrom(const CurveJson& c)
{
    using namespace Stereon;
    const std::string type = c["type"].AsString();
    if (type == "line") {
        return Curve{Line{PointFrom(c["origin"]), VectorFrom(c["direction"])}};
    }
    if (type == "bspline") {
        return Curve{SplineFrom(c)};
    }
    const Vec3 x = VectorFrom(c["x_axis"]);
    const Vec3 y = VectorFrom(c["y_axis"]);
    if (type == "circle") {
        return Curve{Circle{PointFrom(c["center"]), c["radius"].AsNumber(), x, y}};
    }
    if (type == "ellipse") {
        return Curve{
            Ellipse{PointFrom(c["center"]), c["major"].AsNumber(), c["minor"].AsNumber(), x, y}};
    }
    if (type == "parabola") {
        return Curve{Parabola{PointFrom(c["vertex"]), c["focal"].AsNumber(), x, y}};
    }
    return Curve{
        Hyperbola{PointFrom(c["center"]), c["major"].AsNumber(), c["minor"].AsNumber(), x, y}};
}

/// Every case file of the suite, in a stable order.
inline std::vector<std::filesystem::path> CurveCorpusFiles()
{
    const std::filesystem::path root = std::filesystem::path(STEREON_CORPUS_DIR) / "curve-curve";
    std::vector<std::filesystem::path> files;
    for (const auto& entry : std::filesystem::recursive_directory_iterator(root)) {
        if (entry.path().extension() == ".json") {
            files.push_back(entry.path());
        }
    }
    std::ranges::sort(files);
    return files;
}

} // namespace StereonTest
