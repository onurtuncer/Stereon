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

/// Resolves one function into its pointer; false if the context has none.
template <typename F>
bool Resolve(const ProcLookup lookup, F& function, const char* name, const char** missing)
{
    // Converting between function pointer types is how every GL loader
    // works: the lookup returns a generic pointer for the named function.
    function = reinterpret_cast<F>(lookup(name)); // NOLINT
    if (function == nullptr && *missing == nullptr) {
        *missing = name;
    }
    return function != nullptr;
}

} // namespace

bool Load(const ProcLookup lookup, const char** missing)
{
    *missing = nullptr;
    Functions& f = Gl;
    bool ok = true;
    ok = Resolve(lookup, f.Viewport, "glViewport", missing) && ok;
    ok = Resolve(lookup, f.ClearColor, "glClearColor", missing) && ok;
    ok = Resolve(lookup, f.Clear, "glClear", missing) && ok;
    ok = Resolve(lookup, f.Enable, "glEnable", missing) && ok;
    ok = Resolve(lookup, f.Disable, "glDisable", missing) && ok;
    ok = Resolve(lookup, f.BlendFunc, "glBlendFunc", missing) && ok;
    ok = Resolve(lookup, f.CreateShader, "glCreateShader", missing) && ok;
    ok = Resolve(lookup, f.ShaderSource, "glShaderSource", missing) && ok;
    ok = Resolve(lookup, f.CompileShader, "glCompileShader", missing) && ok;
    ok = Resolve(lookup, f.GetShaderiv, "glGetShaderiv", missing) && ok;
    ok = Resolve(lookup, f.GetShaderInfoLog, "glGetShaderInfoLog", missing) && ok;
    ok = Resolve(lookup, f.DeleteShader, "glDeleteShader", missing) && ok;
    ok = Resolve(lookup, f.CreateProgram, "glCreateProgram", missing) && ok;
    ok = Resolve(lookup, f.AttachShader, "glAttachShader", missing) && ok;
    ok = Resolve(lookup, f.LinkProgram, "glLinkProgram", missing) && ok;
    ok = Resolve(lookup, f.GetProgramiv, "glGetProgramiv", missing) && ok;
    ok = Resolve(lookup, f.GetProgramInfoLog, "glGetProgramInfoLog", missing) && ok;
    ok = Resolve(lookup, f.DeleteProgram, "glDeleteProgram", missing) && ok;
    ok = Resolve(lookup, f.UseProgram, "glUseProgram", missing) && ok;
    ok = Resolve(lookup, f.GetUniformLocation, "glGetUniformLocation", missing) && ok;
    ok = Resolve(lookup, f.UniformMatrix4fv, "glUniformMatrix4fv", missing) && ok;
    ok = Resolve(lookup, f.Uniform1f, "glUniform1f", missing) && ok;
    ok = Resolve(lookup, f.GenVertexArrays, "glGenVertexArrays", missing) && ok;
    ok = Resolve(lookup, f.BindVertexArray, "glBindVertexArray", missing) && ok;
    ok = Resolve(lookup, f.DeleteVertexArrays, "glDeleteVertexArrays", missing) && ok;
    ok = Resolve(lookup, f.GenBuffers, "glGenBuffers", missing) && ok;
    ok = Resolve(lookup, f.BindBuffer, "glBindBuffer", missing) && ok;
    ok = Resolve(lookup, f.BufferData, "glBufferData", missing) && ok;
    ok = Resolve(lookup, f.DeleteBuffers, "glDeleteBuffers", missing) && ok;
    ok = Resolve(lookup, f.VertexAttribPointer, "glVertexAttribPointer", missing) && ok;
    ok = Resolve(lookup, f.EnableVertexAttribArray, "glEnableVertexAttribArray", missing) && ok;
    ok = Resolve(lookup, f.DrawArrays, "glDrawArrays", missing) && ok;
    ok = Resolve(lookup, f.PixelStorei, "glPixelStorei", missing) && ok;
    ok = Resolve(lookup, f.ReadPixels, "glReadPixels", missing) && ok;
    return ok;
}

} // namespace StereonGl
