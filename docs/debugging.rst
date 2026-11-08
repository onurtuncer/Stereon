.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Debugging
=========

Index handles, type-erased geometry and copy-on-write pages are opaque to a
debugger that shows raw members: a handle is a bare integer, a ``Curve`` is
96 bytes of buffer and a pointer, a ``PagedArena`` a vector of shared
pointers. ADR-0004 therefore makes pretty-printers part of the design, and
``tools/debugger/`` provides them for GDB and LLDB.

Loading the printers
--------------------

GDB, in a session or in ``~/.gdbinit``:

.. code-block:: text

   source /path/to/Stereon/tools/debugger/stereon_gdb.py

LLDB, in a session or in ``~/.lldbinit``:

.. code-block:: text

   command script import /path/to/Stereon/tools/debugger/stereon_lldb.py

In Visual Studio Code with the C/C++ extension and GDB, add the ``source``
line to ``setupCommands`` in ``launch.json``:

.. code-block:: json

   "setupCommands": [
       { "text": "source ${workspaceFolder}/tools/debugger/stereon_gdb.py" }
   ]

What they show
--------------

Every value prints as Stereon's own ``std::format`` output does, so a
debugger session, a log line and a test failure read alike:

=========================  =======================================================================
Type                       Example
=========================  =======================================================================
``Vec3``, ``Point2``       ``(1, 2.5, -3)``
``Interval``               ``[1, 2]``, ``[]`` when empty
``Box3``                   ``[(0, 0, 0), (1, 2, 3)]``
``VertexId``               ``Vertex#12``, ``Edge#-`` when invalid
``BuilderId``              ``Face#0@0`` (index and generation)
``Oriented``               ``Coedge#42 (reversed)``
``Curve``, ``Surface``     ``Curve(Circle) {Center = (1, 2, 3), Radius = 2.5, ...}``
user geometry              ``Curve((anonymous namespace)::UserCurve) {Base = {...}}``
large geometry             ``Curve(BSplineCurve, shared) {...}``, held behind the shared pointer
``PagedArena``             ``PagedArena of 6 entries in 2 pages, 1 shared``; GDB lists the entries
``JournalId``              ``J#7``
``Sign``, ``Expansion``    ``Negative``, ``Expansion [1, 2]``
=========================  =======================================================================

In a build with checked handles (the debug and sanitizer presets), a handle
of a live shape shows the entity it names, as ``Describe()`` does:

.. code-block:: text

   (gdb) print edge
   $1 = Edge#0 [Spatial Line, tol 1e-07, faces none]

Each frozen snapshot registers its store while it is alive, and the
printers call ``stereon_describe_handle(kind, index, store)`` in the
process. With no live process (a core file), a handle of a release build,
which carries no store, or a snapshot that is gone, the handle shows its
index, ``Edge#0``.

How they work
-------------

Two details of the GCC 16 module build shape the printers. Type names carry
the module they belong to, ``Stereon::Interval@Stereon.Core``, which the
printers' matchers allow and strip. And a ``Curve`` or ``Surface`` holds its
value in a model whose dynamic type neither debugger resolves by itself: the
vtable names template arguments as ``3ull`` where the debug information says
``3``, and LLDB's demangler does not understand GCC's module mangling at all.
The printers read the model's vtable symbol, demangling it with binutils'
``c++filt`` under LLDB, take the wrapped type's name from the model's template
arguments, look that type up, and read the value just past the vtable
pointer, or through the shared pointer for a type too large for the buffer.

Both printers implement the same formatting rules, and
``tests/debugger/check_printers.py`` holds them to one expected table: it runs
``tests/debugger/DebuggerTarget.cpp`` under each debugger in batch mode,
stops in ``StereonDebuggerStop()``, prints every variable and compares. CTest
runs it as ``tools: gdb pretty-printers`` and ``tools: lldb pretty-printers``
in Debug builds, and reports a debugger that is not installed as skipped;
the Linux CI jobs install both. A change to a type's layout that breaks a
printer fails these tests.

Contract violations and journals
--------------------------------

Checked handles (``STEREON_CHECKED_HANDLES``, on in the debug presets) report
a stale or foreign handle with the handle, its store and the slot's state
before raising a contract violation. Evaluator, arena and predicate-domain
violations print the function and the problem the same way. A breakpoint on
``handle_contract_violation`` stops at any of them. The operation journal of
ADR-0023 records what an operation was called with, so a failure can be
replayed in isolation.

Dumps
-----

A dump puts an object in a file for the debug viewer (ADR-0035). It is an
ordinary ``.stn`` archive, ``$STEREON_DUMP_DIR/<name>.stn``, with an optional
``DBUG`` chunk that records the dump's name, the file, line and time of the
call, and debug primitives: points, polylines, boxes, triangle soups and
labels, each with an optional colour. When ``STEREON_DUMP_DIR`` is unset, a
dump writes nothing and costs one environment lookup. A dump of the same
name replaces the earlier one. The file is written to a temporary name and
renamed, so a reader never sees half of it.

.. code-block:: cpp

   #include <Stereon/Io/Dump.hpp>
   import Stereon.Io;

   STEREON_DUMP("after-knot-insertion", curve);  // a Curve, Curve2d, Surface or Shape
   STEREON_DUMP("probe", primitives);            // a DebugChunk of primitives alone

``STEREON_DUMP`` records the call site and ignores errors. ``Stereon::DumpAt``
returns the path written, or the error. From a debugger, call the
out-of-line overloads, which return ``true`` when a file was written:

.. code-block:: text

   (gdb) call Stereon::Dump("edge", *curve)

A ``Shape`` is dumped through ``ShapeArchive`` until the native format has
topology chunks. The archive holds the surface of every face, the curve of
every edge trimmed to its range and the pcurve of every coedge. Its debug
primitives are the topology overlays:

- every vertex, with its tolerance drawn as a box;
- every edge as a polyline coloured by its curve kind (``EdgeCurveColour``);
- a ``Describe()`` label at every vertex and at the middle of every edge.

Display sampling, ``Stereon.Mesh:Sampling``, turns curves into polylines
and surfaces into iso-lines. Its ``DisplayTolerance`` (chord and angle)
serves display only and never reaches modelling code. B-spline polylines
pass through every distinct knot. Unbounded geometry is drawn over a range
that the caller supplies.

stn-dump
--------

``stn-dump`` prints any ``.stn`` file, dumps included, as JSON: the header,
the policy and every chunk with its records and primitives. With
``--pretty`` it prints one member per line and keeps flat arrays on one
line, so two versions of a file diff line by line:

.. code-block:: sh

   stn-dump --pretty before.stn > before.json
   stn-dump --pretty after.stn > after.json
   diff before.json after.json

It is built with the kernel unless ``STEREON_BUILD_TOOLS=OFF``.

The debug viewer
----------------

``stn-view`` shows the dumps of a directory and reloads them when they
change, so a dump made in a running test or a debugger session appears
within a second. Configure with ``-DSTEREON_BUILD_VIEWER=ON``. It needs
GLFW, OpenGL 3.3 and Dear ImGui, which only the viewer may link
(ADR-0035):

.. code-block:: sh

   STEREON_DUMP_DIR=/tmp/dumps ctest -R my-test
   stn-view /tmp/dumps        # or no argument: $STEREON_DUMP_DIR

- Left drag orbits, right or middle drag pans, the wheel zooms, and ``F``
  fits everything in view.
- The **Objects** window lists every dump with its visibility, colour, name
  and source location. Selecting a row highlights the object. A file that
  does not load is listed with its error.
- The **Layers** window shows or hides geometry, control polygons and nets,
  knots, parameter ticks, direction arrows, debug primitives and labels.
- The **Parameter space** window draws plane curves and pcurves.

Curves and surfaces are sampled for display (``Stereon.Mesh:Sampling``).
Unbounded geometry, such as a line, a plane or a cylinder, is drawn around
the bounded content of the scene. The scene behind the window,
``viewer/scene``, is headless and tested without a window, so other
frontends can reuse it.

.. figure:: figures/viewer-geometry.png
   :width: 100%
   :alt: stn-view with a B-spline curve, a B-spline patch and a circle in 3D, and two plane curves in the Parameter space window.

   Geometry dumps. The cubic has a knot of multiplicity three at its
   middle, and the polyline keeps the kink. The control polygon is grey
   and the knots are yellow. The patch is drawn as iso-lines over its
   control net. Arrows show the direction of each curve's parameter. The
   two plane curves are in the Parameter space window.

.. figure:: figures/viewer-topology.png
   :width: 100%
   :alt: stn-view with the closed cylinder of the topology tests, its edges coloured by kind and labelled with Describe().

   A ``Shape`` dump: the closed cylinder of the topology tests. Each face
   is drawn in grey as its surface's iso-lines inside its loops, so the
   disks stop at their circles. The edges are coloured by kind: the seam
   in orange (``Boundary``) and the circles in blue (``Exact``). Vertices
   and edges carry their ``Describe()`` labels.

Faces drawn inside their loops
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A face of a ``Shape`` is a piece of its surface: the region its loops bound
in the surface's parameters. A dump draws it that way. ``ShapeArchive``
samples each loop's pcurves into a closed polygon in the parameter plane,
and ``SampleTrimmedSurface`` (``Stereon.Mesh:Sampling``) cuts each iso-line
of the surface where it crosses the polygons, keeping the stretches inside
by the even-odd rule, so holes stay empty. The lines go into the dump as
debug polylines in grey (``FaceColour()``). Only where a face's loops cannot
be sampled does the dump fall back to the surface itself, trimmed to the
parameter box its boundary covers, as every face was drawn before.

The difference shows on the results of Booleans, whose faces are cut out of
their surfaces. The Steinmetz solid of two cylinders
(:doc:`examples`, *The Steinmetz solids*) has four faces, each a lens of
a cylinder's lateral face:

.. list-table::
   :widths: 50 50

   * - .. figure:: figures/viewer-steinmetz-before.png
          :width: 100%
          :alt: The Steinmetz solid drawn as the two whole cylinders its faces lie on.

          Before: each face drawn over the parameter box of its boundary,
          which for these lenses takes in most of each cylinder.
     - .. figure:: figures/viewer-steinmetz-after.png
          :width: 100%
          :alt: The Steinmetz solid drawn as its four lenses, bounded by the red intersection edges.

          After: each face drawn inside its loops, the four lenses bounded
          by the red ``Intersection`` edges.

A box drilled through by a cylinder (:doc:`examples`, *A drilled box*):
the top and the bottom are squares with a hole, which the lines now skip.

.. list-table::
   :widths: 50 50

   * - .. figure:: figures/viewer-drilled-box-before.png
          :width: 100%
          :alt: A drilled box whose top face's lines run across the hole.

          Before: the top's lines run across the hole.
     - .. figure:: figures/viewer-drilled-box-after.png
          :width: 100%
          :alt: A drilled box whose top face's lines stop at the hole's circle.

          After: they stop at the hole's circle.

Validity findings on the shape
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A ``ValidityReport`` names the entities each finding is about, but as text:
``Edge#8 [Exact Line, tol 1e-05, faces 3|0]``. ``ValidityArchive(shape,
report)`` draws them where they are (ADR-0035, section 4). It writes the
shape as ``ShapeArchive`` does, then every vertex, edge and face a finding
names in red (``FindingColour()``), a coedge as its edge, with a label at the
first of them giving the rule, the entities' short names and the message.
The findings' lines are drawn before the shape's own, so they show on top,
and their labels replace the shape's ``Describe()`` labels, which would bury
them. ``SelfIntersectionUnchecked``, which every Solid-level report holds,
is left out. ``DumpAt(name, shape, report)`` writes the same as a dump:

.. code-block:: cpp

   const Stereon::ValidityReport report =
       Stereon::ValidateTopology(shape, Stereon::ValidityLevel::Solid, ctx);
   if (!report.IsValid()) {
       Stereon::DumpAt("before-commit", shape, report);
   }

.. figure:: figures/viewer-validity.png
   :width: 100%
   :alt: An open box: its four top edges drawn in red, each labelled ShellOpen with the shell and the edge it names.

   A box without its top, checked at the Solid level: the shell uses each
   of the four top edges once, not twice (``ShellOpen``, L01). Each
   finding draws its edge in red and labels it.

Intersection results
~~~~~~~~~~~~~~~~~~~~

``IntersectionArchive`` turns an intersection result into a dump, so a
traced branch can be looked at where it lies:

- for a surface–surface result, both surfaces, each trimmed to its
  requested box when the box is bounded; every branch's 3D curve over its
  range in its contact's colour (``ContactColour()``: red for a crossing,
  orange for a touch, green for an overlap), labelled with its index,
  contact, side transition and certified bound; every isolated point,
  labelled with its bound; and every branch's pcurves on A and on B, which
  the **Parameter space** window draws;
- for a curve–surface result, the curve over its domain and the surface
  over its box, every point in its contact's colour labelled with its
  parameter and contact, and every overlap as the curve's polyline over its
  range with its pcurve.

.. code-block:: cpp

   const auto result = Stereon::Intersect(a, boxA, b, boxB, ctx);
   if (result) {
       if (const auto archive = Stereon::IntersectionArchive(a, boxA, b, boxB, *result)) {
           Stereon::DumpAt("cylinder-sphere", *archive);
       }
   }

.. figure:: figures/viewer-intersection.png
   :width: 100%
   :alt: A cylinder through a sphere off its axis, the two traced branches in red with their labels, and their pcurves on both surfaces in the Parameter space window.

   A cylinder of radius 1 through a sphere of radius 2 off its axis: two
   closed curves, each a traced and certified branch, labelled with its
   bound. The Parameter space window shows the four pcurves, two on each
   surface.

Drawing this figure exposed a gap in the display sampler: a whole sphere,
or ellipsoid, could not be drawn, since the line of constant latitude at a
pole is a circle of radius about :math:`10^{-16}` whose halves turn
sharply, and its refinement never ended. A segment within rounding of a
point now passes the angle test; the angle still bounds every segment with
a length of its own.

``stn-view --screenshot file.png`` renders the window once, writes it as a
PNG and exits; ``--no-labels`` hides the labels. The figures above come
from the dumps of ``docs/figures/ViewerDumps.cpp``, which the Linux viewer
CI job renders under Xvfb with Mesa's software OpenGL. The *before* images
were rendered once, from the same dumps written by the earlier
``ShapeArchive``, and are not regenerated.

Related design decisions
------------------------

- `ADR-0003: Geometry polymorphism via external polymorphism <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0003-geometry-polymorphism.md>`_
- `ADR-0004: Topology storage and handles <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0004-topology-storage-and-handles.md>`_
