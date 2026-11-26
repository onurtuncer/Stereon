.. SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
.. SPDX-License-Identifier: MPL-2.0
..
.. This Source Code Form is subject to the terms of the Mozilla Public
.. License, v. 2.0. If a copy of the MPL was not distributed with this
.. file, You can obtain one at https://mozilla.org/MPL/2.0/.

Native records, replay and language bindings
============================================

Stereon's native archive is defined by its versioned record layout. Read the record type and length, validate payload boundaries and reconstruct geometry through the checked factories. Reject incompatible numerical policies. Debug dumps carry extra information for replay and inspection. The format specification and the code are the authorities for this project-specific algorithm. :doc:`/python` describes the Python interface to the same kernel.

Detailed treatment: :doc:`/data_exchange`.

Implementation notes
--------------------

.. algorithm-sources::
   python/StereonModule.cpp
   src/io/Debug.cpp
   src/io/Debug.cppm
   src/io/Dump.cpp
   src/io/Dump.cppm
   src/io/GeometryArchive.cpp
   src/io/GeometryArchive.cppm
   src/io/Io.cppm
   src/io/Records.cpp
   src/io/Records.cppm
   src/io/include/Stereon/Io/Dump.hpp
   tools/debugger/stereon_gdb.py
   tools/debugger/stereon_lldb.py
   tools/stn-dump/StnDump.cpp
