// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// stn-view's drawing: coloured lines, points and triangles in one shader,
// OpenGL 3.3 core.

#pragma once

#include "Camera.hpp"

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace StereonView {

struct Vertex {
    float X = 0.0F;
    float Y = 0.0F;
    float Z = 0.0F;
    float R = 0.0F;
    float G = 0.0F;
    float B = 0.0F;
    float A = 1.0F;
};

/// What to draw: vertex pairs for lines, single vertices for points and
/// triples for triangles.
struct Batches {
    std::vector<Vertex> Lines;
    std::vector<Vertex> Points;
    std::vector<Vertex> Triangles;
};

class Renderer {
public:
    Renderer() = default;
    Renderer(const Renderer&) = delete;
    Renderer& operator=(const Renderer&) = delete;
    ~Renderer();

    /// Compiles the shader and creates the buffers, in the current context;
    /// false with the reason in `error` otherwise.
    [[nodiscard]] bool Init(std::string& error);
    void Upload(const Batches& batches);
    void Draw(const Matrix& viewProjection, int width, int height, float pointSize);
    /// The framebuffer's pixels, RGBA, rows from the top.
    [[nodiscard]] std::vector<std::uint8_t> ReadPixels(int width, int height) const;

private:
    unsigned int m_Program = 0;
    std::array<unsigned int, 3> m_Arrays{};
    std::array<unsigned int, 3> m_Buffers{};
    std::array<int, 3> m_Counts{};
    int m_MatrixLocation = -1;
    int m_SizeLocation = -1;
};

} // namespace StereonView
