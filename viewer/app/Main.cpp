// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// stn-view, the debug viewer (ADR-0035): shows the dumps of a directory
// ($STEREON_DUMP_DIR by default) and reloads them as they change.
//
//     stn-view [directory] [--screenshot file.png] [--size WIDTHxHEIGHT]
//              [--no-labels] [--no-panel]
//
// Left drag orbits, right or middle drag pans, the wheel zooms, F fits all.
// The Objects window lists every dump with its visibility, colour, name and
// source location; selecting a row highlights the object. The Layers window
// shows or hides the overlays; the Parameter space window draws plane
// curves and pcurves. With --screenshot it renders a few frames, writes the
// window as a PNG and exits, for the manual's figures.

// Every header comes before the import (docs/toolchain/gcc16-feature-probe.md).
#include "Camera.hpp"
#include "Gl.hpp"
#include "Png.hpp"
#include "Renderer.hpp"

#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

import Stereon.Viewer.Scene;

namespace {

namespace V = Stereon::Viewer;
using StereonView::Batches;
using StereonView::Camera;
using StereonView::Renderer;
using StereonView::Vertex;

constexpr V::Rgba Highlight{1.0F, 0.95F, 0.2F, 1.0F};
constexpr double RefreshSeconds = 0.5;
constexpr int ScreenshotFrames = 4;
constexpr std::array<const char*, V::LayerCount> LayerNames{
    "Geometry",   "Control polygons and nets", "Knots", "Parameter ticks",
    "Directions", "Debug primitives"};

struct Arguments {
    std::filesystem::path Directory;
    std::optional<std::filesystem::path> Screenshot;
    int Width = 1280;
    int Height = 800;
    bool Labels = true;
    bool Panel = true;
};

struct ViewState {
    std::array<bool, V::LayerCount> Layers{true, true, true, false, true, true};
    bool Labels = true;
    bool Panel = true;
    std::optional<std::size_t> Selected;
    bool Dirty = true;
};

[[nodiscard]] std::filesystem::path DefaultDirectory()
{
    const char* dump = std::getenv("STEREON_DUMP_DIR");
    return dump != nullptr && *dump != '\0' ? std::filesystem::path(dump)
                                            : std::filesystem::path(".");
}

[[nodiscard]] bool ParseSize(const std::string_view text, Arguments& args)
{
    const std::size_t x = text.find('x');
    if (x == std::string_view::npos) {
        return false;
    }
    args.Width = std::atoi(std::string(text.substr(0, x)).c_str());
    args.Height = std::atoi(std::string(text.substr(x + 1)).c_str());
    return args.Width > 0 && args.Height > 0;
}

[[nodiscard]] std::optional<Arguments> Parse(const std::vector<std::string_view>& argv)
{
    Arguments args;
    args.Directory = DefaultDirectory();
    for (std::size_t i = 0; i < argv.size(); ++i) {
        const std::string_view arg = argv[i];
        const bool hasValue = i + 1 < argv.size();
        if (arg == "--screenshot" && hasValue) {
            args.Screenshot = std::filesystem::path(argv[++i]);
        } else if (arg == "--size" && hasValue) {
            if (!ParseSize(argv[++i], args)) {
                return std::nullopt;
            }
        } else if (arg == "--no-labels") {
            args.Labels = false;
        } else if (arg == "--no-panel") {
            args.Panel = false;
        } else if (!arg.starts_with("--")) {
            args.Directory = std::filesystem::path(arg);
        } else {
            return std::nullopt;
        }
    }
    return args;
}

[[nodiscard]] Vertex ToVertex(const Stereon::Point3& p, const V::Rgba& c)
{
    return Vertex{static_cast<float>(p.X),
                  static_cast<float>(p.Y),
                  static_cast<float>(p.Z),
                  c.R,
                  c.G,
                  c.B,
                  c.A};
}

[[nodiscard]] bool Shown(const ViewState& state, const V::Layer layer)
{
    return state.Layers[static_cast<std::size_t>(layer)];
}

void AddObject(const V::Primitives& content, const ViewState& state,
               const std::optional<V::Rgba> highlight, Batches& out)
{
    // Debug primitives first: they often trace geometry (a shape's edges in
    // their kind's colour), and at equal depth the first line drawn wins.
    for (const bool debug : {true, false}) {
        for (const V::Segment& s : content.Segments) {
            if ((s.Kind == V::Layer::Debug) == debug && Shown(state, s.Kind)) {
                const V::Rgba c = highlight.value_or(s.Colour);
                out.Lines.push_back(ToVertex(s.A, c));
                out.Lines.push_back(ToVertex(s.B, c));
            }
        }
    }
    for (const V::Marker& m : content.Markers) {
        if (Shown(state, m.Kind)) {
            out.Points.push_back(ToVertex(m.Position, highlight.value_or(m.Colour)));
        }
    }
    if (Shown(state, V::Layer::Debug)) {
        for (const V::Triangle& t : content.Triangles) {
            V::Rgba c = highlight.value_or(t.Colour);
            c.A = 0.45F;
            out.Triangles.push_back(ToVertex(t.A, c));
            out.Triangles.push_back(ToVertex(t.B, c));
            out.Triangles.push_back(ToVertex(t.C, c));
        }
    }
}

[[nodiscard]] Batches BuildBatches(const V::Scene& scene, const ViewState& state)
{
    Batches out;
    const auto objects = scene.GetObjects();
    for (std::size_t i = 0; i < objects.size(); ++i) {
        if (objects[i].Visible) {
            const bool selected = state.Selected == i;
            AddObject(objects[i].Content, state,
                      selected ? std::optional<V::Rgba>{Highlight} : std::nullopt, out);
        }
    }
    return out;
}

[[nodiscard]] ImU32 ToImGui(const V::Rgba& c)
{
    return ImGui::ColorConvertFloat4ToU32(ImVec4(c.R, c.G, c.B, c.A));
}

void DrawLabels(const V::Scene& scene, const Camera& camera, const ImVec2 size)
{
    ImDrawList* draw = ImGui::GetBackgroundDrawList();
    for (const V::SceneObject& object : scene.GetObjects()) {
        if (!object.Visible) {
            continue;
        }
        for (const V::Label& label : object.Content.Labels) {
            double x = 0.0;
            double y = 0.0;
            const StereonView::Point p{label.Position.X, label.Position.Y, label.Position.Z};
            if (camera.Project(p, static_cast<double>(size.x), static_cast<double>(size.y), x, y)) {
                draw->AddText(ImVec2(static_cast<float>(x) + 4.0F, static_cast<float>(y) - 6.0F),
                              ToImGui(label.Colour), label.Text.c_str());
            }
        }
    }
}

void ObjectRow(V::Scene& scene, ViewState& state, const std::size_t i)
{
    const V::SceneObject& object = scene.GetObjects()[i];
    ImGui::PushID(static_cast<int>(i));
    bool visible = object.Visible;
    if (ImGui::Checkbox("##visible", &visible)) {
        scene.SetVisible(i, visible);
        state.Dirty = true;
    }
    ImGui::SameLine();
    const V::Rgba c = object.Colour;
    ImGui::ColorButton("##colour", ImVec4(c.R, c.G, c.B, c.A), ImGuiColorEditFlags_NoTooltip);
    ImGui::SameLine();
    if (ImGui::Selectable(object.Name.c_str(), state.Selected == i)) {
        state.Selected = state.Selected == i ? std::nullopt : std::optional<std::size_t>{i};
        state.Dirty = true;
    }
    if (ImGui::IsItemHovered()) {
        ImGui::SetTooltip("%s\n%s:%u%s%s", object.Source.string().c_str(), object.File.c_str(),
                          object.Line, object.Error.empty() ? "" : "\n", object.Error.c_str());
    }
    if (!object.File.empty()) {
        const std::string file = std::filesystem::path(object.File).filename().string();
        ImGui::SameLine();
        if (object.Line > 0) {
            ImGui::TextDisabled("%s:%u", file.c_str(), object.Line);
        } else {
            ImGui::TextDisabled("%s", file.c_str());
        }
    }
    if (!object.Error.empty()) {
        ImGui::SameLine();
        ImGui::TextColored(ImVec4(1.0F, 0.4F, 0.4F, 1.0F), "not loaded");
    }
    ImGui::PopID();
}

void DrawObjectList(V::Scene& scene, ViewState& state)
{
    ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(330, 260), ImGuiCond_FirstUseEver);
    ImGui::Begin("Objects");
    ImGui::TextDisabled("%s", scene.GetDirectory().string().c_str());
    for (std::size_t i = 0; i < scene.GetObjects().size(); ++i) {
        ObjectRow(scene, state, i);
    }
    ImGui::End();
}

void DrawLayers(ViewState& state)
{
    ImGui::SetNextWindowPos(ImVec2(10, 280), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(330, 210), ImGuiCond_FirstUseEver);
    ImGui::Begin("Layers");
    for (std::size_t i = 0; i < V::LayerCount; ++i) {
        state.Dirty = ImGui::Checkbox(LayerNames[i], &state.Layers[i]) || state.Dirty;
    }
    ImGui::Checkbox("Labels", &state.Labels);
    ImGui::Checkbox("Parameter space", &state.Panel);
    ImGui::End();
}

/// The bounds of the visible plane segments, as (min x, min y, max x, max y).
[[nodiscard]] std::optional<std::array<double, 4>> PlaneBounds(const V::Scene& scene)
{
    std::optional<std::array<double, 4>> b;
    for (const V::SceneObject& object : scene.GetObjects()) {
        if (!object.Visible) {
            continue;
        }
        for (const V::PlaneSegment& s : object.Content.PlaneSegments) {
            for (const Stereon::Point2& p : {s.A, s.B}) {
                if (!b) {
                    b = std::array<double, 4>{p.X, p.Y, p.X, p.Y};
                }
                *b = {std::min((*b)[0], p.X), std::min((*b)[1], p.Y), std::max((*b)[2], p.X),
                      std::max((*b)[3], p.Y)};
            }
        }
    }
    return b;
}

void DrawPlanePanel(const V::Scene& scene, const ViewState& state)
{
    const std::optional<std::array<double, 4>> b = PlaneBounds(scene);
    if (!state.Panel || !b) {
        return;
    }
    ImGui::SetNextWindowPos(ImVec2(10, 500), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(330, 280), ImGuiCond_FirstUseEver);
    ImGui::Begin("Parameter space");
    const ImVec2 origin = ImGui::GetCursorScreenPos();
    const ImVec2 room = ImGui::GetContentRegionAvail();
    const double width = std::max((*b)[2] - (*b)[0], 1e-9);
    const double height = std::max((*b)[3] - (*b)[1], 1e-9);
    const double roomX = static_cast<double>(room.x);
    const double roomY = static_cast<double>(room.y);
    const double scale = 0.9 * std::min(roomX / width, roomY / height);
    const auto map = [&](const Stereon::Point2& p) {
        return ImVec2(origin.x + static_cast<float>(0.05 * roomX + (p.X - (*b)[0]) * scale),
                      origin.y + static_cast<float>(0.05 * roomY + ((*b)[3] - p.Y) * scale));
    };
    ImDrawList* draw = ImGui::GetWindowDrawList();
    draw->AddRect(map(Stereon::Point2{(*b)[0], (*b)[3]}), map(Stereon::Point2{(*b)[2], (*b)[1]}),
                  ToImGui(V::Rgba{0.35F, 0.35F, 0.4F, 1.0F}));
    for (const V::SceneObject& object : scene.GetObjects()) {
        for (const V::PlaneSegment& s : object.Content.PlaneSegments) {
            if (object.Visible && Shown(state, s.Kind)) {
                draw->AddLine(map(s.A), map(s.B), ToImGui(s.Colour), 1.5F);
            }
        }
    }
    ImGui::End();
}

void FitAll(const V::Scene& scene, Camera& camera)
{
    const Stereon::Box3 box = scene.GetBounds();
    if (box.IsEmpty()) {
        camera.Fit({-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0});
        return;
    }
    camera.Fit({box.X.GetLower(), box.Y.GetLower(), box.Z.GetLower()},
               {box.X.GetUpper(), box.Y.GetUpper(), box.Z.GetUpper()});
}

void HandleInput(const V::Scene& scene, Camera& camera, const ImVec2 size)
{
    const ImGuiIO& io = ImGui::GetIO();
    if (!io.WantCaptureMouse) {
        const auto dx = static_cast<double>(io.MouseDelta.x);
        const auto dy = static_cast<double>(io.MouseDelta.y);
        if (ImGui::IsMouseDown(ImGuiMouseButton_Left)) {
            camera.Orbit(dx, dy);
        } else if (ImGui::IsMouseDown(ImGuiMouseButton_Right) ||
                   ImGui::IsMouseDown(ImGuiMouseButton_Middle)) {
            camera.Pan(dx, dy, static_cast<double>(size.y));
        }
        if (io.MouseWheel != 0.0F) {
            camera.Zoom(static_cast<double>(io.MouseWheel));
        }
    }
    if (!io.WantCaptureKeyboard && ImGui::IsKeyPressed(ImGuiKey_F)) {
        FitAll(scene, camera);
    }
}

struct App {
    GLFWwindow* Window = nullptr;
    Renderer Gpu;
    Camera View;
    V::Scene Dumps;
    ViewState State;
    double LastRefresh = -1.0;
};

/// One frame; with `pixels`, the frame's image is read back before it is
/// shown.
void Frame(App& app, std::vector<std::uint8_t>* pixels)
{
    if (glfwGetTime() - app.LastRefresh > RefreshSeconds) {
        const bool first = app.LastRefresh < 0.0;
        app.LastRefresh = glfwGetTime();
        if (app.Dumps.Refresh() || first) {
            app.State.Dirty = true;
            if (first) {
                FitAll(app.Dumps, app.View);
            }
        }
    }
    int width = 0;
    int height = 0;
    glfwGetFramebufferSize(app.Window, &width, &height);
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
    const ImVec2 size = ImGui::GetIO().DisplaySize;
    HandleInput(app.Dumps, app.View, size);
    DrawObjectList(app.Dumps, app.State);
    DrawLayers(app.State);
    DrawPlanePanel(app.Dumps, app.State);
    if (app.State.Labels) {
        DrawLabels(app.Dumps, app.View, size);
    }
    if (app.State.Dirty) {
        app.Gpu.Upload(BuildBatches(app.Dumps, app.State));
        app.State.Dirty = false;
    }
    const double aspect = static_cast<double>(width) / std::max(height, 1);
    app.Gpu.Draw(app.View.ViewProjection(aspect), width, height, 5.0F);
    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
    if (pixels != nullptr) {
        *pixels = app.Gpu.ReadPixels(width, height);
    }
    glfwSwapBuffers(app.Window);
}

[[nodiscard]] GLFWwindow* OpenWindow(const Arguments& args)
{
    glfwSetErrorCallback([](const int code, const char* message) {
        std::cerr << "stn-view: GLFW error " << code << ": " << message << '\n';
    });
    if (glfwInit() == GLFW_FALSE) {
        return nullptr;
    }
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GLFW_TRUE); // macOS needs it for 3.3 core
    GLFWwindow* window = glfwCreateWindow(args.Width, args.Height, "stn-view", nullptr, nullptr);
    if (window != nullptr) {
        glfwMakeContextCurrent(window);
        glfwSwapInterval(1);
    }
    return window;
}

[[nodiscard]] int Run(App& app, const Arguments& args)
{
    if (!args.Screenshot) {
        while (glfwWindowShouldClose(app.Window) == GLFW_FALSE) {
            glfwPollEvents();
            Frame(app, nullptr);
        }
        return 0;
    }
    std::vector<std::uint8_t> pixels;
    for (int frame = 0; frame < ScreenshotFrames; ++frame) {
        glfwPollEvents();
        Frame(app, frame + 1 == ScreenshotFrames ? &pixels : nullptr);
    }
    int width = 0;
    int height = 0;
    glfwGetFramebufferSize(app.Window, &width, &height);
    if (!StereonView::WritePng(*args.Screenshot, width, height, pixels)) {
        std::cerr << "stn-view: cannot write " << args.Screenshot->string() << '\n';
        return 1;
    }
    return 0;
}

[[nodiscard]] int Start(const Arguments& args)
{
    GLFWwindow* window = OpenWindow(args);
    if (window == nullptr) {
        std::cerr << "stn-view: cannot open an OpenGL 3.3 core window\n";
        glfwTerminate();
        return 1;
    }
    const char* missing = nullptr;
    if (!StereonGl::Load(&glfwGetProcAddress, &missing)) {
        std::cerr << "stn-view: the OpenGL context lacks " << missing << '\n';
        glfwTerminate();
        return 1;
    }
    ImGui::CreateContext();
    ImGui::GetIO().IniFilename = nullptr;
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330 core");
    int status = 1;
    {
        App app{window, Renderer{}, Camera{}, V::Scene(args.Directory), ViewState{}, -1.0};
        app.State.Labels = args.Labels;
        app.State.Panel = args.Panel;
        std::string error;
        if (app.Gpu.Init(error)) {
            status = Run(app, args);
        } else {
            std::cerr << "stn-view: " << error << '\n';
        }
    }
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();
    return status;
}

} // namespace

int main(const int argc, char** argv)
{
    const std::vector<std::string_view> raw(argv + 1, argv + argc);
    const std::optional<Arguments> args = Parse(raw);
    if (!args) {
        std::cerr << "usage: stn-view [directory] [--screenshot file.png] "
                     "[--size WIDTHxHEIGHT] [--no-labels] [--no-panel]\n";
        return 2;
    }
    return Start(*args);
}
