// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// Stereon.Surface: surfaces (ADR-0003): the SurfaceGeometry concept, the
/// built-in Plane, Cylinder, Cone, Sphere and Torus, B-spline and NURBS
/// surfaces, and the type-erased Surface.
export module Stereon.Surface;

export import Stereon.Curve;

export import :BSplineSurface;
export import :Capabilities;
export import :Certify;
export import :Concepts;
export import :Cone;
export import :Conversion;
export import :Cylinder;
export import :Erased;
export import :Fitting;
export import :GridEvaluation;
export import :Plane;
export import :Properties;
export import :Sphere;
export import :Swept;
export import :Torus;
export import :Trimmed;
export import :Visit;
