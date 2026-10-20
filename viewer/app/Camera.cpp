// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Camera.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>

namespace StereonView {

namespace {

constexpr double FieldOfView = 35.0 * std::numbers::pi / 180.0;
constexpr double OrbitPerPixel = 0.008;
constexpr double ZoomPerStep = 0.88;
constexpr double PitchLimit = 1.55;

[[nodiscard]] Point Minus(const Point& a, const Point& b)
{
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

[[nodiscard]] Point Cross(const Point& a, const Point& b)
{
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

[[nodiscard]] double Dot(const Point& a, const Point& b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

[[nodiscard]] Point Unit(const Point& a)
{
    const double length = std::sqrt(Dot(a, a));
    return length > 0.0 ? Point{a[0] / length, a[1] / length, a[2] / length} : a;
}

/// The camera's axes: right, up and backwards (from the target to the eye).
struct Axes {
    Point Right;
    Point Up;
    Point Back;
};

[[nodiscard]] Axes AxesOf(const Point& eye, const Point& target)
{
    const Point back = Unit(Minus(eye, target));
    const Point right = Unit(Cross({0.0, 0.0, 1.0}, back));
    return {right, Cross(back, right), back};
}

/// Column-major 4x4 product a * b.
[[nodiscard]] std::array<double, 16> Multiply(const std::array<double, 16>& a,
                                              const std::array<double, 16>& b)
{
    std::array<double, 16> m{};
    for (int c = 0; c < 4; ++c) {
        for (int r = 0; r < 4; ++r) {
            double sum = 0.0;
            for (int k = 0; k < 4; ++k) {
                sum +=
                    a[static_cast<std::size_t>(k * 4 + r)] * b[static_cast<std::size_t>(c * 4 + k)];
            }
            m[static_cast<std::size_t>(c * 4 + r)] = sum;
        }
    }
    return m;
}

} // namespace

void Camera::Fit(const Point& lower, const Point& upper)
{
    m_Target = {(lower[0] + upper[0]) / 2, (lower[1] + upper[1]) / 2, (lower[2] + upper[2]) / 2};
    m_Radius = std::max(0.5 * std::sqrt(Dot(Minus(upper, lower), Minus(upper, lower))), 1e-6);
    m_Distance = 1.15 * m_Radius / std::sin(FieldOfView / 2);
}

void Camera::Orbit(const double dx, const double dy)
{
    m_Yaw -= dx * OrbitPerPixel;
    m_Pitch = std::clamp(m_Pitch + dy * OrbitPerPixel, -PitchLimit, PitchLimit);
}

void Camera::Pan(const double dx, const double dy, const double viewportHeight)
{
    const Axes axes = AxesOf(Eye(), m_Target);
    const double perPixel =
        2.0 * m_Distance * std::tan(FieldOfView / 2) / std::max(viewportHeight, 1.0);
    for (std::size_t i = 0; i < 3; ++i) {
        m_Target[i] += (-dx * axes.Right[i] + dy * axes.Up[i]) * perPixel;
    }
}

void Camera::Zoom(const double steps)
{
    m_Distance = std::max(m_Distance * std::pow(ZoomPerStep, steps), 1e-6);
}

Point Camera::Eye() const
{
    return {m_Target[0] + m_Distance * std::cos(m_Pitch) * std::cos(m_Yaw),
            m_Target[1] + m_Distance * std::cos(m_Pitch) * std::sin(m_Yaw),
            m_Target[2] + m_Distance * std::sin(m_Pitch)};
}

Matrix Camera::ViewProjection(const double aspect) const
{
    const Point eye = Eye();
    const Axes a = AxesOf(eye, m_Target);
    const std::array<double, 16> view{a.Right[0],         a.Up[0],         a.Back[0],         0.0,
                                      a.Right[1],         a.Up[1],         a.Back[1],         0.0,
                                      a.Right[2],         a.Up[2],         a.Back[2],         0.0,
                                      -Dot(a.Right, eye), -Dot(a.Up, eye), -Dot(a.Back, eye), 1.0};
    const double nearPlane = std::max(m_Distance - 4.0 * m_Radius, m_Distance * 1e-3);
    const double farPlane = m_Distance + 4.0 * m_Radius;
    const double f = 1.0 / std::tan(FieldOfView / 2);
    const std::array<double, 16> projection{f / std::max(aspect, 1e-6),
                                            0.0,
                                            0.0,
                                            0.0,
                                            0.0,
                                            f,
                                            0.0,
                                            0.0,
                                            0.0,
                                            0.0,
                                            (farPlane + nearPlane) / (nearPlane - farPlane),
                                            -1.0,
                                            0.0,
                                            0.0,
                                            2.0 * farPlane * nearPlane / (nearPlane - farPlane),
                                            0.0};
    const std::array<double, 16> m = Multiply(projection, view);
    Matrix out{};
    for (std::size_t i = 0; i < 16; ++i) {
        out[i] = static_cast<float>(m[i]);
    }
    return out;
}

bool Camera::Project(const Point& p, const double width, const double height, double& x,
                     double& y) const
{
    const Matrix m = ViewProjection(width / std::max(height, 1.0));
    std::array<double, 4> clip{};
    for (std::size_t r = 0; r < 4; ++r) {
        clip[r] = static_cast<double>(m[r]) * p[0] + static_cast<double>(m[4 + r]) * p[1] +
                  static_cast<double>(m[8 + r]) * p[2] + static_cast<double>(m[12 + r]);
    }
    if (!(clip[3] > 0.0)) {
        return false;
    }
    x = (clip[0] / clip[3] * 0.5 + 0.5) * width;
    y = (0.5 - clip[1] / clip[3] * 0.5) * height;
    return true;
}

} // namespace StereonView
