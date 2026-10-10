.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

API reference
=============

The API reference is generated from the source comments by Doxygen and
rendered here through Breathe.

.. stereon-api-begin

.. only:: have_doxygen

   .. todo::

      Replace the full index below with curated pages per library
      (Core, Robust, Curve, Surface, Intersect, Topo, Build, Boolean,
      Blend, Mesh, Io) once the public API settles.

   .. doxygenindex::
      :project: Stereon

.. stereon-api-end

.. only:: not have_doxygen

   .. note::

      This build has no Doxygen XML, so the API reference is empty.
      Configure with ``-DSTEREON_BUILD_DOCS=ON`` and Doxygen installed,
      or set ``STEREON_RUN_DOXYGEN=1`` when running ``sphinx-build``.
