.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Algorithm reference
===================

Start with a method family below. Each page links to the mathematical
chapter, identifies the reference behind the method, and exposes implementation
notes alongside the revision-specific source. A reference to a numerical
method does not certify a particular implementation. Supported domains and
failure conditions are defined by the detailed chapters and regression tests.

.. toctree::
   :maxdepth: 1

   arithmetic
   core
   spline
   fitting
   geometry
   certification
   sensitivity
   planar
   space
   curve-surface
   ssi
   topology
   boolean
   winding
   mesh
   io
   viewer
   verification
   benchmarks
   integration

Maintaining this reference
--------------------------

List each first-party implementation, header, script or test with the family
it implements or exercises. For a new algorithm in an existing file,
document its inputs, construction, decision or stopping rule, limitations,
and primary literature in that family's narrative or linked chapter.
Use a specification or project design record for a project-specific format
or policy. Attribute adaptations separately from the original method.

The source inventory check detects missing and stale files. On pull requests
it also compares named C++ and Python definitions against the base revision.
A new definition needs a triple-slash comment, a Python docstring, or an
explicit description naming the routine in its family's narrative. Review
changes inside existing function bodies for new mathematical methods and
verify that the cited source supports the particular technique. C++ triple-slash
comments and Python docstrings are rendered directly from the checkout.

To audit existing kernel definitions, install the pinned analysis tools and run:

.. code-block:: console

   python -m pip install -r tools/requirements.txt
   python tools/check_algorithm_docs.py --audit-existing

The audit lists definitions where the scanner found neither an adjacent
documentation comment nor an exact routine name in the manual. Some entries
are wrappers or have their contract on a declaration. Conversely, a name
match does not establish that a whole method or its citations are correct.
Use the list to direct source review alongside the family pages.
