.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Camera projection, scene generation and image encoding
======================================================

The viewer derives camera axes from the eye and target, composes view and projection transforms, and builds display primitives from sampled geometry. The PNG path writes image bytes as stored DEFLATE blocks inside a zlib stream. PNG chunks carry CRCs and the zlib stream carries an Adler checksum. For the encoding algorithms consult the normative specifications: `PNG <https://www.w3.org/TR/png-3/>`__, `zlib (RFC 1950) <https://www.rfc-editor.org/rfc/rfc1950>`__ and `DEFLATE (RFC 1951) <https://www.rfc-editor.org/rfc/rfc1951>`__.

For camera transforms and rendering, see :cite:p:`akenine2018rendering`.

Detailed treatment: :doc:`/debugging`.

Implementation notes
--------------------

.. algorithm-sources::
   docs/figures/FigureData.cpp
   docs/figures/ViewerDumps.cpp
   docs/figures/make_fitting_figures.py
   viewer/app/Camera.cpp
   viewer/app/Camera.hpp
   viewer/app/Gl.cpp
   viewer/app/Gl.hpp
   viewer/app/Main.cpp
   viewer/app/Png.cpp
   viewer/app/Png.hpp
   viewer/app/Renderer.cpp
   viewer/app/Renderer.hpp
   viewer/scene/Primitives.cpp
   viewer/scene/Scene.cpp
   viewer/scene/Scene.cppm
