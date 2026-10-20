// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The OpenGL 3.3 core functions stn-view draws with, loaded through the
// window system's procedure lookup (glfwGetProcAddress). Written here rather
// than generated, so that no GL header or loader is a dependency (ADR-0035,
// Consequences: the loader's licence): the types and constants are the
// ones the OpenGL 3.3 core specification fixes, and only what the renderer
// uses is declared. Dear ImGui's OpenGL 3 backend loads its own copies.

#pragma once

#include <cstddef>
#include <cstdint>

#if defined(_WIN32) && !defined(_WIN64)
#define STEREON_GL_API __stdcall
#else
#define STEREON_GL_API
#endif

namespace StereonGl {

using GLenum = unsigned int;
using GLbitfield = unsigned int;
using GLuint = unsigned int;
using GLint = int;
using GLsizei = int;
using GLboolean = unsigned char;
using GLfloat = float;
using GLchar = char;
using GLsizeiptr = std::ptrdiff_t;

inline constexpr GLenum Points = 0x0000;
inline constexpr GLenum Lines = 0x0001;
inline constexpr GLenum Triangles = 0x0004;
inline constexpr GLenum FloatType = 0x1406;
inline constexpr GLenum DepthTest = 0x0B71;
inline constexpr GLenum Blend = 0x0BE2;
inline constexpr GLenum SrcAlpha = 0x0302;
inline constexpr GLenum OneMinusSrcAlpha = 0x0303;
inline constexpr GLenum ArrayBuffer = 0x8892;
inline constexpr GLenum DynamicDraw = 0x88E8;
inline constexpr GLenum FragmentShader = 0x8B30;
inline constexpr GLenum VertexShader = 0x8B31;
inline constexpr GLenum CompileStatus = 0x8B81;
inline constexpr GLenum LinkStatus = 0x8B82;
inline constexpr GLenum ProgramPointSize = 0x8642;
inline constexpr GLenum Rgba = 0x1908;
inline constexpr GLenum UnsignedByte = 0x1401;
inline constexpr GLenum PackAlignment = 0x0D05;
inline constexpr GLbitfield ColorBufferBit = 0x00004000;
inline constexpr GLbitfield DepthBufferBit = 0x00000100;
inline constexpr GLboolean False = 0;

struct Functions {
    void(STEREON_GL_API* Viewport)(GLint, GLint, GLsizei, GLsizei);
    void(STEREON_GL_API* ClearColor)(GLfloat, GLfloat, GLfloat, GLfloat);
    void(STEREON_GL_API* Clear)(GLbitfield);
    void(STEREON_GL_API* Enable)(GLenum);
    void(STEREON_GL_API* Disable)(GLenum);
    void(STEREON_GL_API* BlendFunc)(GLenum, GLenum);
    GLuint(STEREON_GL_API* CreateShader)(GLenum);
    void(STEREON_GL_API* ShaderSource)(GLuint, GLsizei, const GLchar* const*, const GLint*);
    void(STEREON_GL_API* CompileShader)(GLuint);
    void(STEREON_GL_API* GetShaderiv)(GLuint, GLenum, GLint*);
    void(STEREON_GL_API* GetShaderInfoLog)(GLuint, GLsizei, GLsizei*, GLchar*);
    void(STEREON_GL_API* DeleteShader)(GLuint);
    GLuint(STEREON_GL_API* CreateProgram)();
    void(STEREON_GL_API* AttachShader)(GLuint, GLuint);
    void(STEREON_GL_API* LinkProgram)(GLuint);
    void(STEREON_GL_API* GetProgramiv)(GLuint, GLenum, GLint*);
    void(STEREON_GL_API* GetProgramInfoLog)(GLuint, GLsizei, GLsizei*, GLchar*);
    void(STEREON_GL_API* DeleteProgram)(GLuint);
    void(STEREON_GL_API* UseProgram)(GLuint);
    GLint(STEREON_GL_API* GetUniformLocation)(GLuint, const GLchar*);
    void(STEREON_GL_API* UniformMatrix4fv)(GLint, GLsizei, GLboolean, const GLfloat*);
    void(STEREON_GL_API* Uniform1f)(GLint, GLfloat);
    void(STEREON_GL_API* GenVertexArrays)(GLsizei, GLuint*);
    void(STEREON_GL_API* BindVertexArray)(GLuint);
    void(STEREON_GL_API* DeleteVertexArrays)(GLsizei, const GLuint*);
    void(STEREON_GL_API* GenBuffers)(GLsizei, GLuint*);
    void(STEREON_GL_API* BindBuffer)(GLenum, GLuint);
    void(STEREON_GL_API* BufferData)(GLenum, GLsizeiptr, const void*, GLenum);
    void(STEREON_GL_API* DeleteBuffers)(GLsizei, const GLuint*);
    void(STEREON_GL_API* VertexAttribPointer)(GLuint, GLint, GLenum, GLboolean, GLsizei,
                                              const void*);
    void(STEREON_GL_API* EnableVertexAttribArray)(GLuint);
    void(STEREON_GL_API* DrawArrays)(GLenum, GLint, GLsizei);
    void(STEREON_GL_API* PixelStorei)(GLenum, GLint);
    void(STEREON_GL_API* ReadPixels)(GLint, GLint, GLsizei, GLsizei, GLenum, GLenum, void*);
};

/// The window system's lookup: glfwGetProcAddress.
using ProcLookup = void (*(*)(const char*))();

/// The loaded functions; valid after Load returned true.
extern Functions Gl;

/// Loads every function; false, naming the first missing one in `missing`,
/// if the context lacks any.
[[nodiscard]] bool Load(ProcLookup lookup, const char** missing);

} // namespace StereonGl
