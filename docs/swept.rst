.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Extruded and revolved geometry
==============================

A profile curve moved along a direction is an extrusion; turned about an
axis, a revolution (STEP's ``surface_of_linear_extrusion`` and
``surface_of_revolution``, OCCT's ``Geom_SurfaceOfLinearExtrusion`` and
``Geom_SurfaceOfRevolution``). Stereon has them in two forms.

``ExtrusionSurface`` and ``RevolutionSurface`` wrap any profile in an
erased ``Curve`` and evaluate the sweep itself: an extrusion is unbounded
along its direction, a revolution runs a full turn and is periodic in the
angle. They are the surfaces a face carries before it is trimmed, and the
ones data exchange produces. *The wrapper types* below describe them.

``MakeExtrusionSurface`` and ``MakeRevolutionSurface`` construct finite
swept patches as exact tensor-product ``BSplineSurface`` values from a
``BSplineCurve`` profile, the form the wrappers convert to with
``ToBSpline`` over a range. They do not construct faces or solids. Convert
an analytic or trimmed profile with ``ToBSpline`` first; that conversion
may change its parameter. Both factories validate the profile before
accessing its poles or weights, then validate the generated surface. They
return ``InvalidInput`` for malformed data or non-finite generated values.
There is no fitting budget: these are exact rational constructions up to
floating-point rounding.

Extrusion
---------

``MakeExtrusionSurface(profile, direction, range)`` constructs
``S(u,v) = Eval(profile,u) + v * direction``. The profile parameter is
preserved in u. The v range must be finite and increasing. The finite,
nonzero direction is a displacement per unit v, not a unit vector; it is
not normalised. Rational profile weights are duplicated along the linear
sweep direction. Polynomial profiles remain polynomial.

Revolution
----------

``MakeRevolutionSurface(profile, origin, axis, angles)`` rotates the profile
about the oriented line through origin along axis, using the right-hand
rule. The axis is normalised; its length must be finite and nonzero.
The increasing angular interval must be finite and at most one turn.
Intervals may cross the conventional angular seam.

The profile parameter is preserved in v. The u parameter is the rational
quadratic circle parameter: it agrees with angle at the arc knots, but
not between them. Do not interpret u derivatives as derivatives with
respect to angle. Each angular segment spans at most a quarter turn.
For one segment starting at a with span d, its normalised parameter t
corresponds to angle ``a + 2*atan2(t*sin(d/2), 1-t+t*cos(d/2))``.

Each profile pole is decomposed into an axial and a radial component.
The circle control polygon rotates the radial component and leaves the
axial component fixed; tensor-product weights are the products of profile
and circle weights. This supports rational, non-planar profiles as well
as conventional planar meridians (The NURBS Book, section 8.5).

Limits and validation
---------------------

The wrapper types
-----------------

``MakeExtrusion(profile, direction)`` gives an ``ExtrusionSurface``,
:math:`S(u, v) = C(u) + v\,\mathbf d`: the profile's parameter in
:math:`u`, with the profile's domain and period, and the whole real line
in :math:`v`, the direction being a displacement per unit :math:`v`, not
normalised. :math:`S^{(k,0)} = C^{(k)}`, :math:`S^{(0,1)} = \mathbf d`,
every other partial is zero. Its bounding box is the profile's swept
along the direction, unbounded in every axis the direction has a
component along.

``MakeRevolution(profile, origin, axis)`` gives a ``RevolutionSurface``,
:math:`S(u, v) = O + \mathbf a(v) + \cos u\,\mathbf r(v) + \sin u\,
(\mathbf n \times \mathbf r(v))`, the profile point
:math:`C(v) = O + \mathbf a(v) + \mathbf r(v)` split into its part
:math:`\mathbf a` along the unit axis :math:`\mathbf n` and its radial
part :math:`\mathbf r`, turned by the angle :math:`u` about the axis with
the right-hand rule: :math:`u \in [0, 2\pi]`, periodic, and :math:`v`
the profile's parameter. :math:`S^{(k,l)}` is the radial part of
:math:`C^{(l)}(v)` turned by :math:`u` and :math:`k` quarter turns, plus
the axial part of :math:`C^{(l)}` when :math:`k = 0`; the turn's
derivatives cost nothing beyond the profile's. Unlike the NURBS patch, the
wrapper's :math:`u` *is* the angle, so :math:`u` derivatives are
derivatives with respect to angle. The profile need not be planar.

Reversing an extrusion reverses its profile, which flips :math:`S_u` and
with it the normal; reversing a revolution turns the other way about the
axis, ``Eval(Reversed(s), (2 pi - u, v))`` being ``Eval(s, (u, v))``.
Transforming moves the profile and the direction or the origin and axis;
a reflection also turns a revolution's axis round, so that the turn's
sense follows the map and ``Eval(Transformed(s, x), uv)`` stays
``Transformed(Eval(s, uv), x)``.

**Projection.** Where the kernel knows the profile's exact NURBS form over
a finite domain (a circle, an ellipse, a B-spline, or a trimmed curve of
any built-in curve, a line segment included, reversed or nested:
:doc:`trimmed`, *Conversion over a subrange*), projection is exact.
An extrusion keeps that form flattened along the direction, its poles
moved onto the plane through the origin normal to it, with a
``CurveProjector``: the closest point of the surface to :math:`\mathbf p`
lies over the closest point of the flattened curve to :math:`\mathbf p`
flattened the same way, and :math:`v` follows from the displacement along
the direction. A revolution keeps the NURBS form of its full turn with a
``SurfaceProjector``, and maps the result back to the angle and the
profile's parameter, projecting the form's point onto the profile where
the form does not keep its parameter (a conic's rational arcs). A line
profile, which has no finite domain, extrudes to a plane and is handled as
one. For other profiles, those of a user's type or an unbounded conic,
``Project`` is approximate: the profile's own closest point, turned into
the point's half-plane for a revolution a few times over, with the sweep
parameter chosen for it. ``GetProjection()`` and ``GetPiece()`` are null
in that case, and ``BoundingBox`` of such a revolution is the box of the
cylinder about the axis that holds the profile's box turned.

``ToBSpline(extrusion, range)`` is ``MakeExtrusionSurface`` of the
profile's form over ``range.X``, which must lie in the profile's domain or
span at most one period, over the finite ``range.Y``;
``ToBSpline(revolution, range)`` is ``MakeRevolutionSurface`` through the
angles ``range.X``, at most one turn, of the form over ``range.Y``.
``NotSupported`` for a profile without a known form. The wrappers are not
built-in types for ``VisitBuiltin``, as a trimmed curve is not; they are
wrapped in a ``Surface`` like any user type, found by
``Surface::Target<ExtrusionSurface>()``, and the erased ``ToBSpline``
reaches them.

Limits and validation
---------------------

A revolution's axis may have any finite, non-zero length: it is scaled by
the power of two that brings its largest component into :math:`[1, 2)`
before it is normalised, so an axis of :math:`10^{200}` or
:math:`10^{-200}` has a direction like any other, and :math:`2^k` times an
axis gives the same unit axis bit for bit. An extrusion keeps its direction
as given and divides by its squared length, so ``MakeExtrusion`` requires
that squared length to be a finite, non-zero double.

Collapsed patches, axis crossings and poles on the axis are allowed.
Construction does not certify regularity, injectivity, orientation of a
solid, or absence of self-intersections. Singular local-property queries
retain the existing structured-error behaviour. Projection has the same
limitations as other B-spline surfaces, including uncertified final
Newton leaves, and is approximate for the profiles named above. The
geometry archive stores a wrapper as its profile with the direction, or
the origin and axis, and rebuilds the projection data on load
(:doc:`data_exchange`, format version 3).

Tests compare extrusion with direct translation and revolution with
Rodrigues rotation at interior parameters, the wrappers with their
profiles, with finite differences, with a torus and with dense sampling of
their surfaces, cover derivatives, reversal, reflection and invalid inputs,
and check native round trips of the patches. The optional OCCT suite
compares the patches and the wrappers against its extrusion and
revolution evaluators and the rectangular trimmed surface against its
own.
