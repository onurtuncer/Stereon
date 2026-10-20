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
       // InvalidInput (cannot write, bad policy) or NotSupported (a user type);
       // an existing part.stn is left as it was
   }
   const Result<LoadedGeometry> loaded = LoadGeometry("part.stn");
   if (loaded) {
       // loaded->Archive.Policy is the stored policy, not this kernel's defaults;
       // loaded->Warnings lists the optional chunks the reader skipped.
   }

``SaveGeometry`` writes the archive to ``part.stn.partial`` beside the
target, flushes and closes it, checking each step, and only then renames it
over ``part.stn``. A failure at any step, including one the stream reports
only when its buffer reaches the disk, is returned as ``InvalidInput``, and
leaves an existing file untouched.

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

The format has four chunks, all required, each exactly once and in any
order. Every writer stores all four, empty or not, so the reader refuses
an archive that lacks one (a file cut after a complete chunk, whose
checksums are all intact) or repeats one (which would replace what the
first decoded), with a ``ParseError`` naming the chunk:

- ``POLI``, the rest of the numerical policy: the angular resolution, the
  tolerance ceiling and the frame allowance under which the geometry was
  accepted. The resolution must not exceed the tolerance ceiling, as in a
  ``Context``. The frame allowance must be this kernel's
  ``FrameAllowance``, the one its factories check frames against when they
  validate the geometry on write and on read; the writer refuses another
  with ``InvalidInput``, and the reader with ``NotSupported``, rather than
  claim to restore a policy it did not apply;
- ``CRVS``, ``PCRV`` and ``SRFS``: curves in space, curves in a surface's
  parameter plane, and surfaces. A payload is a record count (u64) and the
  records. A record is the type's name (a u16 length and its bytes, such as
  ``Circle`` or ``BSplineSurface``), the length of the value (u64) and the
  value as the reflection serialiser writes it: its members in declaration
  order, a ``std::vector`` as its element count and its elements
  (:doc:`geometry`, *Serialisation*).

Records are named, not numbered, so adding a built-in type never changes
what an existing name means. A record of a type that came in a later
format version than the file states (a ``TrimmedCurve`` in a version-1
file, a sweep in a version-2 one) is refused with a ``ParseError``: no
writer of that version produced it.

Format version 2 adds trimmed curves (:doc:`trimmed`) to ``CRVS`` and
``PCRV``, named ``TrimmedCurve`` and ``TrimmedCurve2d``. The value of a
trimmed-curve record is its basis, written as a record of the same chunk
(name, length and value), followed by the trim: the range's ends and the
sense, ``Lower`` and ``Upper`` (f64) and ``Reversed`` (one byte), as the
reflection serialiser writes them. A basis may itself be a trimmed curve,
nested at most ``MaxTrimmedNesting`` (8) deep. The writer refuses deeper
nesting with ``NotSupported``, and the reader refuses it with
``ParseError``, so a hostile file cannot exhaust the stack by recursion.
The piece's NURBS form and projector are not stored: they are rebuilt from
the basis on load.

Format version 3 adds trimmed surfaces and the sweep wrappers
(:doc:`trimmed`, :doc:`swept`) to ``SRFS``. The value of a
``TrimmedSurface`` record is its basis as a record of the same chunk, then
the trim: ``LowerU``, ``UpperU``, ``LowerV`` and ``UpperV`` (f64) and
``Reversed`` (one byte). The value of an ``ExtrusionSurface`` record is
its profile as a record of the ``CRVS`` kind, so a trimmed profile is
stored as one, then ``Direction`` (three f64); a ``RevolutionSurface``
record is its profile likewise, then ``Origin`` and ``Axis`` (three f64
each, the axis as the surface keeps it, of unit length: the reader keeps it
as stored, since normalising it again could move it by an ulp, and refuses
one that is not of unit length). Bases and
profiles nest under the one bound of ``MaxTrimmedNesting``, trims and
sweeps counted together, and the exact patches and projectors are rebuilt
on load as for trimmed curves. Each version is the next without its new
records, so a version-1 or version-2 file needs no migration.

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
  a trimmed curve goes through ``MakeTrimmedCurve`` after its basis, once
  its range has been checked finite and non-empty;
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
in format version 1, ``geometry-v2.stn`` the same with trimmed curves
added, among them an arc across a circle's seam, a reversed piece of a
B-spline and a trimmed curve of a trimmed curve, and ``geometry-v3.stn``
that with trimmed surfaces and sweeps added: a cylinder patch across its
seam, a reversed corner of a B-spline surface, an extrusion of a trimmed
profile, a revolution and a band of it. A test reads each and
requires the kernel to write the same bytes back, apart from the versions
in the header, so every later version keeps reading the files of every
earlier one exactly. When the format changes, the version goes up, a
migration reads the old one if it needs one, and the old files stay in the
repository unchanged: only the current version's file is generated, once,
when that version is new, never to make a test pass.

Limits
~~~~~~

- Only the built-in types and trimmed curves of them have records. Saving
  a user's curve or surface type, or a trimmed curve of one, fails with
  ``NotSupported``; a registration mechanism for user types
  belongs with the plug-in work of ADR-0003.
- Files are read whole into memory. The chunked layout allows streaming
  and optional zstd compression (ADR-0012); neither is implemented yet.
- Compilers without static reflection (Clang and MSVC so far) report
  ``NotSupported`` for every call.

Related design decisions
------------------------


