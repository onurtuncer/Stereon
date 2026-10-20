// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Renderer.hpp"

#include "Gl.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace StereonView {

namespace {

using namespace StereonGl;

constexpr const char* VertexSource = R"(#version 330 core
layout(location = 0) in vec3 position;
layout(location = 1) in vec4 colour;
uniform mat4 viewProjection;
uniform float pointSize;
out vec4 vertexColour;
void main()
{
    gl_Position = viewProjection * vec4(position, 1.0);
    gl_PointSize = pointSize;
    vertexColour = colour;
}
)";

constexpr const char* FragmentSource = R"(#version 330 core
in vec4 vertexColour;
out vec4 fragmentColour;
void main()
{
    fragmentColour = vertexColour;
}
)";

constexpr std::array<GLenum, 3> Modes{Lines, Points, Triangles};

[[nodiscard]] GLuint Compile(const GLenum kind, const char* source, std::string& error)
{
    const GLuint shader = Gl.CreateShader(kind);
    Gl.ShaderSource(shader, 1, &source, nullptr);
    Gl.CompileShader(shader);
    GLint ok = 0;
    Gl.GetShaderiv(shader, CompileStatus, &ok);
    if (ok == 0) {
        std::string log(1024, '\0');
        Gl.GetShaderInfoLog(shader, static_cast<GLsizei>(log.size()), nullptr, log.data());
        error = "shader: " + log;
        Gl.DeleteShader(shader);
        return 0;
    }
    return shader;
}

} // namespace

Renderer::~Renderer()
{
    if (m_Program != 0) {
        Gl.DeleteBuffers(3, m_Buffers.data());
        Gl.DeleteVertexArrays(3, m_Arrays.data());
        Gl.DeleteProgram(m_Program);
    }
}

bool Renderer::Init(std::string& error)
{
    const GLuint vertex = Compile(VertexShader, VertexSource, error);
    const GLuint fragment = Compile(FragmentShader, FragmentSource, error);
    if (vertex == 0 || fragment == 0) {
        return false;
    }
    m_Program = Gl.CreateProgram();
    Gl.AttachShader(m_Program, vertex);
    Gl.AttachShader(m_Program, fragment);
    Gl.LinkProgram(m_Program);
    Gl.DeleteShader(vertex);
    Gl.DeleteShader(fragment);
    GLint ok = 0;
    Gl.GetProgramiv(m_Program, LinkStatus, &ok);
    if (ok == 0) {
        std::string log(1024, '\0');
        Gl.GetProgramInfoLog(m_Program, static_cast<GLsizei>(log.size()), nullptr, log.data());
        error = "program: " + log;
        return false;
    }
    m_MatrixLocation = Gl.GetUniformLocation(m_Program, "viewProjection");
    m_SizeLocation = Gl.GetUniformLocation(m_Program, "pointSize");
    Gl.GenVertexArrays(3, m_Arrays.data());
    Gl.GenBuffers(3, m_Buffers.data());
    for (std::size_t i = 0; i < 3; ++i) {
        Gl.BindVertexArray(m_Arrays[i]);
        Gl.BindBuffer(ArrayBuffer, m_Buffers[i]);
        constexpr auto stride = static_cast<GLsizei>(sizeof(Vertex));
        Gl.VertexAttribPointer(0, 3, FloatType, False, stride, nullptr);
        Gl.VertexAttribPointer(
            1, 4, FloatType, False, stride,
            reinterpret_cast<const void*>(3 * sizeof(float))); // NOLINT: GL offset
        Gl.EnableVertexAttribArray(0);
        Gl.EnableVertexAttribArray(1);
    }
    Gl.BindVertexArray(0);
    return true;
}

void Renderer::Upload(const Batches& batches)
{
    const std::array<const std::vector<Vertex>*, 3> lists{&batches.Lines, &batches.Points,
                                                          &batches.Triangles};
    for (std::size_t i = 0; i < 3; ++i) {
        Gl.BindBuffer(ArrayBuffer, m_Buffers[i]);
        Gl.BufferData(ArrayBuffer, static_cast<GLsizeiptr>(lists[i]->size() * sizeof(Vertex)),
                      lists[i]->data(), DynamicDraw);
        m_Counts[i] = static_cast<int>(lists[i]->size());
    }
}

void Renderer::Draw(const Matrix& viewProjection, const int width, const int height,
                    const float pointSize)
{
    Gl.Viewport(0, 0, width, height);
    Gl.ClearColor(0.11F, 0.12F, 0.14F, 1.0F);
    Gl.Clear(ColorBufferBit | DepthBufferBit);
    Gl.Enable(DepthTest);
    Gl.Enable(Blend);
    Gl.Enable(ProgramPointSize);
    Gl.BlendFunc(SrcAlpha, OneMinusSrcAlpha);
    Gl.UseProgram(m_Program);
    Gl.UniformMatrix4fv(m_MatrixLocation, 1, False, viewProjection.data());
    Gl.Uniform1f(m_SizeLocation, pointSize);
    // Triangles first, so that lines and points on them stay visible.
    for (const std::size_t i : {std::size_t{2}, std::size_t{0}, std::size_t{1}}) {
        if (m_Counts[i] > 0) {
            Gl.BindVertexArray(m_Arrays[i]);
            Gl.DrawArrays(Modes[i], 0, m_Counts[i]);
        }
    }
    Gl.BindVertexArray(0);
    Gl.UseProgram(0);
    Gl.Disable(DepthTest);
}

std::vector<std::uint8_t> Renderer::ReadPixels(const int width, const int height) const
{
    const auto w = static_cast<std::size_t>(width);
    const auto h = static_cast<std::size_t>(height);
    std::vector<std::uint8_t> bottomUp(w * h * 4);
    Gl.PixelStorei(PackAlignment, 1);
    Gl.ReadPixels(0, 0, width, height, StereonGl::Rgba, UnsignedByte, bottomUp.data());
    std::vector<std::uint8_t> topDown(bottomUp.size());
    for (std::size_t row = 0; row < h; ++row) {
        const std::size_t from = (h - 1 - row) * w * 4;
        for (std::size_t k = 0; k < w * 4; ++k) {
            topDown[row * w * 4 + k] = bottomUp[from + k];
        }
    }
    return topDown;
}

} // namespace StereonView
