.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Extruded and revolved geometry
==============================

``MakeExtrusionSurface`` and ``MakeRevolutionSurface`` construct finite
swept patches as exact tensor-product ``BSplineSurface`` values. They do
not construct faces or solids. The result uses the existing surface
operations, type erasure and version-1 geometry archive; construction
history is not retained by the archive.

Pass a valid ``BSplineCurve`` profile. Convert an analytic or trimmed
profile with ``ToBSpline`` first; that conversion may change its parameter.
Both factories validate the profile before accessing its poles or weights,
then validate the generated surface. They return ``InvalidInput`` for
malformed data or non-finite generated values. There is no fitting budget:
these are exact rational constructions up to floating-point rounding.

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

Collapsed patches, axis crossings and poles on the axis are allowed.
Construction does not certify regularity, injectivity, orientation of a
solid, or absence of self-intersections. Singular local-property queries
retain the existing structured-error behaviour. Projection has the same
limitations as other B-spline surfaces, including uncertified final
Newton leaves. Unbounded sweeps and evaluator-only custom profiles are
not supported by these factories. Dedicated extrusion/revolution wrapper
types and rectangular trimmed-surface wrappers remain future work.

Tests compare extrusion with direct translation and revolution with
Rodrigues rotation at interior parameters, cover derivatives and invalid
inputs, and check native round trips. The optional OCCT suite compares
against its extrusion and revolution evaluators.
