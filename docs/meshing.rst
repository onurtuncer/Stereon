.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Meshing
=======

A mesh approximates exact geometry by simple pieces: polylines for curves,
triangles for surfaces. A kernel needs meshes for two quite different
purposes, and Stereon keeps them apart (ADR-0014, ADR-0034, ADR-0035):

- **Display sampling** draws geometry on a screen. It must look right, and
  it is free to be coarse wherever the eye cannot tell.
- **Analysis meshes** feed a CFD or FEA solver. They must be watertight
  across face boundaries, carry the tags of the model's regions, and
  approximate the geometry to a stated, certified deviation.

Neither kind is ever the model. A mesh is derived from the B-rep and never
changes its geometry or topology. Its tolerances are display or
discretisation tolerances of their own, which no modelling decision reads
(ADR-0028; :doc:`tolerances`).

The ``mesh`` library today provides display sampling of curves and surfaces.
Analysis-grade, watertight and tagged surface meshes are planned for v0.7,
and certified discretisation for machining for v0.3 (ADR-0044; see the
roadmap).

Display sampling
----------------

``SampleCurve`` turns a curve into a polyline, and ``SampleSurface`` turns a
surface into a net of iso-parameter lines. Both are steered by a
``DisplayTolerance``:

.. list-table::
   :header-rows: 1
   :widths: 20 16 64

   * - Field
     - Default
     - Meaning
   * - ``Chord``
     - 0.01 mm
     - The largest distance allowed between the curve and a segment of its
       polyline.
   * - ``Angle``
     - 0.1 rad
     - The largest turn allowed between consecutive half-segments.
   * - ``MaxPoints``
     - 65 536
     - The most points one polyline may have. Sampling that would need more
       is refused, never truncated.

Adaptive chordal subdivision
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A polyline is built by adaptive chordal subdivision, the standard method of
curve flattening :cite:p:`piegl1997nurbs`. The parameter range is first
divided at its **breakpoints**, the ends and every distinct knot of a
B-spline (or of a trimmed curve's B-spline basis) inside the range. A knot
is where a spline may have a kink, and a kink must be drawn as a vertex,
not smoothed over. Each span between breakpoints starts as four equal
segments. Then each segment :math:`[t_0, t_1]` is tested
(:numref:`fig-mesh-chord-test`):

.. math::

   \max_{k = 1, \dots, 7}
   \operatorname{dist}\Bigl(C\bigl(t_0 + \tfrac{k}{8}(t_1 - t_0)\bigr),\;
   \overline{C(t_0)\,C(t_1)}\Bigr) \le \texttt{Chord}
   \quad\text{and}\quad
   \angle\bigl(C(t_m) - C(t_0),\; C(t_1) - C(t_m)\bigr) \le \texttt{Angle},

with :math:`t_m = \tfrac12(t_0 + t_1)`. A segment that passes is kept. One
that fails is split at :math:`t_m`, and both halves are tested in turn, to
a depth of 40 halvings below the span's first four segments: :math:`2^{-40}`
of the span, far below anything a display resolves.

.. _fig-mesh-chord-test:

.. figure:: _static/tikz/mesh-chord-test.png
   :width: 92%
   :alt: A curve over one segment, its chord with a tolerance band, the
         curve's points at the interior eighths, and the turn angle at the
         midpoint.

   The test of one segment. The curve's points at the seven interior eighths
   (blue) are measured against the chord (red). Here the curve leaves the
   band by :math:`\delta`, and the chords of the two halves (green) turn by
   :math:`\pi - \theta`, so the segment is split at :math:`t_m` and each half
   tested again.

Four details make the method dependable:

- Starting each span with four segments, not one, catches the S-shaped cubic
  whose midpoint happens to lie on its chord. A single chord would accept it.
- The seven check points, not only the midpoint, catch a deviation that the
  midpoint misses.
- The angle test refines where the curve turns sharply but stays close to its
  chord, near a cusp for instance.
- Failure is explicit. If the tolerance cannot be met within ``MaxPoints``
  points or the depth limit, because the parameter's precision runs out or
  the curve evaluates to a non-finite point, ``SampleCurve`` returns
  ``NotSupported`` instead of a polyline that breaks its tolerance.

The method guarantees the tolerance at the check points, not between them:
it is a display, not a certificate. A certified discretisation, with a
proved bound on the deviation everywhere, is the subject of ADR-0044.

Surfaces as iso-lines
~~~~~~~~~~~~~~~~~~~~~

``SampleSurface`` draws ``linesPerDirection`` evenly spaced iso-lines in each
parameter direction, nine by default, ends included. It adds one more at
every distinct interior knot of a B-spline surface, or of a trimmed
surface's B-spline basis, so that patch boundaries are always drawn. Each
iso-line is a curve sampled as above. The lines along :math:`u` come first,
by ascending :math:`v`.

Unbounded geometry
~~~~~~~~~~~~~~~~~~

A line, a plane or the height of a cylinder has no natural extent. It is
sampled over the range the caller gives, clipped to the domain. Without a
range it is refused, since the scene, not the geometry, knows how far to
draw it. A periodic curve or surface is sampled over its period unless a
range says otherwise.

Using it
--------

.. code-block:: cpp

   import Stereon.Mesh;

   using namespace Stereon;

   const Curve circle(*MakeCircle(Point3{}, 5.0, Vec3{1, 0, 0}, Vec3{0, 1, 0}));
   const DisplayTolerance display{.Chord = 1e-3, .Angle = 0.05};
   const Result<Polyline> polyline = SampleCurve(circle, display);
   // polyline->Points[k] is the curve at polyline->Parameters[k].

The debug viewer draws everything through these two functions
(:doc:`debugging`).

Implementation
--------------

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Code
     - Role
   * - ``src/mesh/Sampling.cppm``
     - ``DisplayTolerance``, ``Polyline``, ``IsoLine``, ``Breakpoints``,
       ``SampleCurve``, ``SampleSurface``.
   * - ``src/mesh/Sampling.cpp``
     - The subdivision: ``InitialSegments = 4``, ``ChordChecks = 8``,
       ``MaxDepth = 40``.
   * - ``tests/mesh/``
     - Tolerance met at the check points, knots as vertices, refusal of
       unbounded and impossible requests.

Related design decisions
------------------------

- `ADR-0014: Mesh generation boundary <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0014-mesh-generation-boundary.md>`_
- `ADR-0034: Kernel presentation boundary <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0034-kernel-presentation-boundary.md>`_
- `ADR-0035: Debug viewer <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0035-debug-viewer.md>`_
- `ADR-0028: Purpose-specific numerical tolerances <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0028-purpose-specific-numerical-tolerances.md>`_
- `ADR-0044: B-rep to tri-dexel conversion <https://github.com/onurtuncer/Stereon/blob/main/docs/adr/0044-brep-to-tri-dexel-conversion.md>`_
