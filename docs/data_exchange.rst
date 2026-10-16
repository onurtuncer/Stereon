.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Data exchange
=============

STEP, STL, OBJ and glTF import and export, and the native file format with
its versioning rules. So far the native format stores geometry; shapes,
with their topology, history and attributes, follow with the topology
library (ADR-0012 places native shape persistence in v0.3).

.. todo::

   STEP, STL, OBJ and glTF import and export.

Saving and loading geometry
---------------------------

``Stereon.Io`` saves curves, plane curves and surfaces of the built-in
types, together with the numerical policy they were built under, and reads
them back bit for bit:

.. code-block:: cpp

   import Stereon.Io;
   using namespace Stereon;

   GeometryArchive archive;
   archive.Policy = PolicyOf(context);           // the model's numerical policy
   archive.Curves.push_back(Curve(circle));
   archive.Surfaces.push_back(Surface(nurbs));

   if (const Result<void> saved = SaveGeometry("part.stn", archive); !saved) {
       // InvalidInput (cannot write, bad policy) or NotSupported (a user type)
   }
   const Result<LoadedGeometry> loaded = LoadGeometry("part.stn");
   if (loaded) {
       // loaded->Archive.Policy is the stored policy, not this kernel's defaults;
       // loaded->Warnings lists the optional chunks the reader skipped.
   }

``WriteGeometry`` and ``ReadGeometry`` do the same with bytes, and
``GeometryToJson`` gives a JSON view of an archive for debugging and diffs,
the ``stn-dump`` view of ADR-0012.

Layout
~~~~~~

Every value is little-endian, doubles as their IEEE 754 bit pattern, so the
same geometry gives the same bytes on every platform (ADR-0007).

.. list-table::
   :header-rows: 1
   :widths: 22 12 66

   * - Field
     - Bytes
     - Content
   * - Magic
     - 4
     - ``STN\x1A``
   * - Format version
     - 4
     - 1 for the layout described here
   * - Writer version
     - 6 + 2
     - The kernel version that wrote the file (three u16), two zero bytes
   * - Length unit
     - 4
     - 1 = millimetre (ADR-0001)
   * - Resolution
     - 8
     - The model's linear resolution, in that unit
   * - Reserved
     - 4
     - Zero
   * - Chunks
     - …
     - Each: a tag of four ASCII letters, flags (u32, bit 0 = required), the
       payload length (u64), the CRC-32 of the payload (u32), the payload

Format version 1 has four chunks, all required:

- ``POLI``, the rest of the numerical policy: the angular resolution, the
  tolerance ceiling and the frame allowance under which the geometry was
  accepted;
- ``CRVS``, ``PCRV`` and ``SRFS``: curves in space, curves in a surface's
  parameter plane, and surfaces. A payload is a record count (u64) and the
  records. A record is the type's name (a u16 length and its bytes, such as
  ``Circle`` or ``BSplineSurface``), the length of the value (u64) and the
  value as the reflection serialiser writes it: its members in declaration
  order, a ``std::vector`` as its element count and its elements
  (:doc:`geometry`, *Serialisation*).

Records are named, not numbered, so adding a built-in type never changes
what an existing name means.

Reading
~~~~~~~

The reader treats a file as data from outside the kernel (ADR-0006,
ADR-0016). Everything it reads is checked, and every failure is a
``ParseError`` that says where:

- the magic, and a format version this reader knows; a newer version is
  refused, never guessed at;
- each chunk's length against what is left, and its CRC-32: a damaged chunk
  is refused;
- every record's value through its type's checked factory (``MakeCircle``,
  ``MakeBSplineSurface`` and so on), so a negative radius or a bad knot
  vector in a file is caught exactly as it would be from any other source;
- every element count against the bytes left, before anything is allocated,
  so a hostile count cannot exhaust memory.

Following ADR-0012, an unknown **required** chunk is refused, because a
model read without part of its geometry is not the model that was saved,
and an unknown **optional** chunk is skipped with a line in ``Warnings``. A
newer writer that adds required data therefore bumps the format version.

The numerical policy is part of the model's meaning (ADR-0028). The reader
returns the stored values, refuses a file that has none, and never
substitutes its own defaults: a model saved with a resolution of
:math:`2 \cdot 10^{-5}` mm loads with :math:`2 \cdot 10^{-5}` mm, whatever
the reading program's ``Context`` says. What to do when the two differ is
the caller's decision.

Compatibility
~~~~~~~~~~~~~

``tests/io/data/geometry-v1.stn`` holds one of every built-in type written
in format version 1. A test reads it and requires the kernel to write the
same bytes back, so every later version keeps reading version-1 files
exactly. When the format changes, the version goes up, a migration reads
the old one, and the old file stays in the repository: it is regenerated
only for a new format version, never to make a test pass.

Limits
~~~~~~

- Only the built-in types have records. Saving a user's curve or surface
  type fails with ``NotSupported``; a registration mechanism for user types
  belongs with the plug-in work of ADR-0003.
- Files are read whole into memory. The chunked layout allows streaming
  and optional zstd compression (ADR-0012); neither is implemented yet.
- Compilers without static reflection (Clang and MSVC so far) report
  ``NotSupported`` for every call.

Related design decisions
------------------------


