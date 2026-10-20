// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// stn-view's orbit camera: a target, a distance and two angles, z up.

#pragma once

#include <array>

namespace StereonView {

using Point = std::array<double, 3>;
using Matrix = std::array<float, 16>;

class Camera {
public:
    /// Looks at the box from the default direction, the whole box in view.
    void Fit(const Point& lower, const Point& upper);
    /// Turns around the target by a mouse drag, in pixels.
    void Orbit(double dx, double dy);
    /// Moves the target by a mouse drag, in pixels of a viewport this high.
    void Pan(double dx, double dy, double viewportHeight);
    /// Moves closer (positive steps) or away.
    void Zoom(double steps);

    /// The view and projection together, column-major as OpenGL takes it.
    [[nodiscard]] Matrix ViewProjection(double aspect) const;
    /// Where the point appears in a viewport of this size, in pixels from
    /// the top left; false when it is behind the camera.
    [[nodiscard]] bool Project(const Point& p, double width, double height, double& x,
                               double& y) const;

private:
    [[nodiscard]] Point Eye() const;

    Point m_Target{0.0, 0.0, 0.0};
    double m_Distance = 10.0;
    double m_Radius = 5.0;
    double m_Yaw = -0.9;
    double m_Pitch = 0.5;
};

} // namespace StereonView
