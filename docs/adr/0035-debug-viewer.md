# ADR-0035: Debug viewer over dump files with GLFW, OpenGL and Dear ImGui

- **Status:** Proposed
- **Date:** 2026-10-17
- **Owner:** Onur Tuncer
- **Phase:** 0 (geometry display); topology overlays follow the topology prototype
- **Scope:** The Phase 0 debug viewer: how objects reach it, where display sampling lives, its backend and dependencies, and what v0 shows
- **Related:** ADR-0008, ADR-0009, ADR-0012, ADR-0014, ADR-0018, ADR-0020, ADR-0023, ADR-0028, ADR-0033, ADR-0034

## Context

Kernel work needs a quick way to look at an intermediate object: a NURBS curve after knot insertion, a trimmed curve and its basis, a surface with its control net, the pcurves of an edge, and, once topology exists, a body with the entities a validity report names. Numbers in a test log do not show a reversed coedge or a pcurve that leaves its domain.

ADR-0034 fixes the boundary: the kernel stays headless, presentation is an adapter over public kernel interfaces, and display tolerances never become modelling tolerances. ADR-0033 asks that reusable algorithms not live only in the inspection application's source tree. Neither ADR chooses a backend. TODO.md §8 left the choice open between "Qt6 or Vulkan".

OCCT's equivalent is DRAW. It is a Tcl interpreter with its own 2D and 3D views, and `vinit` adds the AIS/OpenGL viewer. Its most valued habit is `DBRep_Set("s", shape)` from a breakpoint, which puts an object on screen in the middle of a debugging session. Its costs are a scripting language woven into the test system, a viewer that has to run inside the process under test, and toolkits (TKDraw, TKTopTest, TKViewerTest) whose display code grew toward a presentation framework.

The decision is needed now. Curves, surfaces, NURBS and trimmed curves exist, and the topology prototype (TODO.md §7) will want the viewer from its first day.

## Options

How objects reach the viewer:

1. **In-process viewer, DRAW-style.** The code under test opens a window or talks to an interpreter that owns one. This gives the fastest round trip, but every test that wants to show something needs a GUI-capable process, and sanitizer, CI and Python runs are awkward.
2. **Dump files of serialised kernel objects.** The code under test writes files, and a separate viewer reads them. It works from any process, CI can keep the files as artifacts, and they diff. The viewer links the kernel to read and sample them.
3. **Dump files of pre-sampled display primitives.** The viewer stays trivial and kernel-free, but the dump loses everything except the drawing: no control nets, knots, pcurve domains, tolerances or entity labels unless each is sampled in advance.
4. **Export only** (OBJ, VTK, then ParaView or Blender). There is nothing to build, but there are no pcurves, labels, tolerances or validity overlays. ADR-0034 already limits this option to early inspection.

The backend:

5. **Qt 6** (widgets and `QOpenGLWidget`). A full toolkit, but LGPL-3.0 or commercial, which ADR-0008 allows only as an optional separately built adapter. It is a heavy dependency for CI and MSYS2, and moc with GCC 16 modules is untried.
6. **Raw Vulkan.** Explicit and modern, but weeks of setup code (instances, swapchains, pipelines, synchronisation) that a debug viewer gains nothing from.
7. **GLFW + OpenGL 3.3 core + Dear ImGui.** Small, permissively licensed (GLFW under Zlib, Dear ImGui under MIT), packaged by MSYS2 and the Linux distributions, and immediate-mode UI suits a tool that changes often. OpenGL is a legacy API on macOS, but 3.3 core still runs there, and the viewer needs nothing newer.

## Decision

Adopt options 2 and 7.

### 1. Objects reach the viewer as `.stn` dump files

- `Stereon.Io:Dump` writes one file per call: `<STEREON_DUMP_DIR>/<name>.stn`. A dump is an ordinary `.stn` file (ADR-0012, `Stereon.Io:GeometryArchive`). It holds the geometry chunks, and the topology chunks once a `Shape` exists, plus one **optional** `DBUG` chunk. The `DBUG` chunk records the dump's name, the source location and time of the call, and the debug primitives that are not kernel objects: points, polylines, boxes, triangle soups and text labels, each with an optional colour. The `DBUG` chunk is optional, so readers that do not know it skip it with a warning, and adding it does not change the format version.
- `STEREON_DUMP(name, obj)` in `<Stereon/Io/Dump.hpp>` records the call's source location. Exported, non-template, non-inline `Dump(name, obj)` overloads for the built-in geometry, and later for `Shape`, can be called from a GDB or LLDB breakpoint. They take the place of DRAW's `DBRep_Set`.
- When `STEREON_DUMP_DIR` is unset, a dump writes nothing and costs one environment lookup. A dump writes a temporary file and renames it, so a reader never sees a half-written file. A later dump with the same name replaces the earlier one, which is what live reload (section 4) shows. Names are limited to `[A-Za-z0-9_.-]`.
- In Phase 0, tests, tools and applications can dump. Kernel algorithms below `io` cannot. When an algorithm first needs to dump an intermediate object, that goes through a sink on the `Context`, as logging and the journal do (ADR-0023), and is decided then.
- `stn-dump` (ADR-0012) reads dumps like any other `.stn` file, so the JSON view and the viewer share one reader.

### 2. Display sampling is a headless kernel-family service

- `Stereon.Mesh:Sampling` turns curves into polylines and surfaces into iso-lines. A display tolerance (chord deviation and angle) controls the sampling. That tolerance is a purpose of its own (ADR-0028, ADR-0034 §5) and never reaches modelling code. NURBS are sampled through every distinct knot. Unbounded geometry (lines, planes) is drawn over a range that the dump states, or else over the bounding box of the bounded objects in the scene.
- Sampling depends only on public curve and surface interfaces. Display tessellation of faces arrives with the meshing work (ADR-0014) and replaces nothing here.

### 3. The viewer is an application under `viewer/`

- `viewer/scene` is a headless library. It reads dump files and builds a scene of display primitives through `Stereon.Mesh:Sampling`, and Catch2 tests exercise it without a window.
- `viewer/app` builds the executable `stn-view`: GLFW for the window and input, OpenGL 3.3 core for drawing, Dear ImGui for panels, and a generated GL loader. It is built only with `STEREON_BUILD_VIEWER=ON`, which is off by default.
- GLFW, Dear ImGui and the loader are declared in `cmake/StereonDependencies.cmake` with their licences, like every dependency. `tools/check_dependencies.py` admits them only under `viewer/`, as it admits OpenCASCADE only under `tests/` and `bench/`. Nothing under `src/` may reference them, GL headers or `viewer/` (ADR-0034 §1).

### 4. What v0 shows

- Points, polylines, boxes, triangle soups and labels from the `DBUG` chunk.
- Every built-in curve and surface. Optional overlays: parameter ticks, direction arrows, control polygons and nets, and knot positions.
- A 2D panel for `Curve2d` values and pcurves in parameter space, with the surface's parameter domain drawn (ADR-0009).
- Orbit, pan and zoom, fit all, and an object list with visibility, colour, name and source location. Selecting a row in the list highlights its object.
- Live reload: the viewer polls the dump directory with `std::filesystem` (no file-watcher dependency) and reloads changed files within a second.
- Once the topology prototype lands: edges coloured by curve kind (ADR-0009), coedge orientation arrows, vertex and edge tolerances drawn as spheres and tubes, `Describe()` labels, and a `ValidityReport` overlay that highlights the entities each rule ID names (ADR-0018).

v0 deliberately has no picking of rendered primitives, no shaded faces and no editing. It has no scripting language either: the Python package (ADR-0020) takes that role, and can call `Dump` once it is bound. Picking that resolves to validated model references is step 3 of ADR-0034 and comes when interactive modelling needs it.

## Consequences

- Any test, tool, sanitizer run or debugger session can put an object on screen, and CI can keep the dumps of a failing test as artifacts. The viewer never runs inside the code under test.
- One file format serves persistence, debugging and diffing, so the `.stn` reader is exercised every time someone debugs.
- Display sampling is reusable by other frontends and testable headlessly, so it starts the `Stereon.Mesh` library before analysis-grade meshing.
- Three new dependencies, confined to an optional target. The GL loader's licence must be confirmed against the ADR-0008 allow-list when it is added. If the generated code's licence is not on the list, the loader that Dear ImGui's OpenGL 3 backend ships is the fallback.
- Linux CI needs the X11 or Wayland development packages to build the viewer. CI builds `stn-view` and runs the `viewer/scene` tests, but does not open a window.
- Dumping from inside kernel algorithms is not yet possible. That is a deliberate gap, closed through the `Context` when first needed.
- The viewer is not the future presentation layer. When document services and validated interaction arrive (ADR-0033, ADR-0034), they may reuse `viewer/scene` or replace it. Either way the kernel is unaffected.

## Verification

- The kernel builds and its tests pass with `STEREON_BUILD_VIEWER=OFF`. A check fails if anything under `src/` references GLFW, Dear ImGui, a GL header or `viewer/`.
- `tools/check_dependencies.py` accepts the three new dependencies with their licences and rejects them outside `viewer/`.
- A dump of every built-in curve and surface type reads back bit-identical through `Stereon.Io:GeometryArchive`. A reader that does not know the `DBUG` chunk loads the geometry with one warning.
- `Stereon.Mesh:Sampling` property tests: the chord deviation stays within the display tolerance on every built-in type, NURBS polylines pass through every distinct knot, and unbounded geometry is clipped to its range.
- `viewer/scene` tests build the expected primitives from a fixture dump without a window.
- Manual acceptance before Phase 0 exit: the hand-built box and cylinder of TODO.md §7 are displayed with edge kinds, orientations, tolerances and `Describe()` labels, a seeded checker fault is highlighted by its rule ID, and a re-dump appears in a running viewer within a second.

## Relationship to existing decisions

- [ADR-0034](0034-kernel-presentation-boundary.md): this is the "thin inspection viewer" of its step 2. It chooses the graphics API and UI toolkit that ADR-0034 leaves open, and it keeps the kernel headless.
- [ADR-0033](0033-modular-sdk-product-scope.md): display sampling sits in the meshing capability, and the viewer belongs to the inspection application and uses public interfaces only.
- [ADR-0012](0012-native-file-format.md): dumps are `.stn` files, and the `DBUG` chunk is an optional chunk under its rules.
- [ADR-0008](0008-licence.md): GLFW (Zlib) and Dear ImGui (MIT) are on the allow-list, and Qt was rejected partly because it is not.
- [ADR-0014](0014-mesh-generation-boundary.md) and [ADR-0028](0028-purpose-specific-numerical-tolerances.md): display sampling and its tolerance are kept separate from analysis-grade meshing and from modelling tolerances.
- [ADR-0023](0023-logging-and-journal-sinks.md): an in-kernel dump, when needed, follows the same sink pattern as logging.

This proposal supersedes none of these decisions.
