// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Gl.hpp"

namespace StereonGl {

Functions Gl{};

namespace {

/// Resolves one function into its pointer, recording the first name the
/// context lacks in `missing`.
template <typename F>
void Resolve(const ProcLookup lookup, F& function, const char* name, const char** missing)
{
    // Converting between function pointer types is how every GL loader
    // works: the lookup returns a generic pointer for the named function.
    function = reinterpret_cast<F>(lookup(name)); // NOLINT
    if (function == nullptr && *missing == nullptr) {
        *missing = name;
    }
}

} // namespace

bool Load(const ProcLookup lookup, const char** missing)
{
    *missing = nullptr;
    Functions& f = Gl;
    Resolve(lookup, f.Viewport, "glViewport", missing);
    Resolve(lookup, f.ClearColor, "glClearColor", missing);
    Resolve(lookup, f.Clear, "glClear", missing);
    Resolve(lookup, f.Enable, "glEnable", missing);
    Resolve(lookup, f.Disable, "glDisable", missing);
    Resolve(lookup, f.BlendFunc, "glBlendFunc", missing);
    Resolve(lookup, f.CreateShader, "glCreateShader", missing);
    Resolve(lookup, f.ShaderSource, "glShaderSource", missing);
    Resolve(lookup, f.CompileShader, "glCompileShader", missing);
    Resolve(lookup, f.GetShaderiv, "glGetShaderiv", missing);
    Resolve(lookup, f.GetShaderInfoLog, "glGetShaderInfoLog", missing);
    Resolve(lookup, f.DeleteShader, "glDeleteShader", missing);
    Resolve(lookup, f.CreateProgram, "glCreateProgram", missing);
    Resolve(lookup, f.AttachShader, "glAttachShader", missing);
    Resolve(lookup, f.LinkProgram, "glLinkProgram", missing);
    Resolve(lookup, f.GetProgramiv, "glGetProgramiv", missing);
    Resolve(lookup, f.GetProgramInfoLog, "glGetProgramInfoLog", missing);
    Resolve(lookup, f.DeleteProgram, "glDeleteProgram", missing);
    Resolve(lookup, f.UseProgram, "glUseProgram", missing);
    Resolve(lookup, f.GetUniformLocation, "glGetUniformLocation", missing);
    Resolve(lookup, f.UniformMatrix4fv, "glUniformMatrix4fv", missing);
    Resolve(lookup, f.Uniform1f, "glUniform1f", missing);
    Resolve(lookup, f.GenVertexArrays, "glGenVertexArrays", missing);
    Resolve(lookup, f.BindVertexArray, "glBindVertexArray", missing);
    Resolve(lookup, f.DeleteVertexArrays, "glDeleteVertexArrays", missing);
    Resolve(lookup, f.GenBuffers, "glGenBuffers", missing);
    Resolve(lookup, f.BindBuffer, "glBindBuffer", missing);
    Resolve(lookup, f.BufferData, "glBufferData", missing);
    Resolve(lookup, f.DeleteBuffers, "glDeleteBuffers", missing);
    Resolve(lookup, f.VertexAttribPointer, "glVertexAttribPointer", missing);
    Resolve(lookup, f.EnableVertexAttribArray, "glEnableVertexAttribArray", missing);
    Resolve(lookup, f.DrawArrays, "glDrawArrays", missing);
    Resolve(lookup, f.PixelStorei, "glPixelStorei", missing);
    Resolve(lookup, f.ReadPixels, "glReadPixels", missing);
    return *missing == nullptr;
}

} // namespace StereonGl
