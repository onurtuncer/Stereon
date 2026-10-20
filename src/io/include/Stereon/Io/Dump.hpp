// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef STEREON_IO_DUMP_HPP
#define STEREON_IO_DUMP_HPP

/// Dumps an object for the debug viewer (ADR-0035) as
/// `$STEREON_DUMP_DIR/<name>.stn`, recording where the call was made; does
/// nothing when STEREON_DUMP_DIR is unset. The object is anything
/// Stereon::DumpAt takes: a Curve, a Curve2d, a Surface, a Shape, a
/// GeometryArchive or a DebugChunk of primitives. Errors are ignored: a
/// dump is a debugging aid and never changes what the code does. Use it in
/// a unit that imports Stereon.Io:
///
///     #include <Stereon/Io/Dump.hpp>
///     import Stereon.Io;
///     ...
///     STEREON_DUMP("after-knot-insertion", curve);
#define STEREON_DUMP(name, object)                                                                 \
    static_cast<void>(::Stereon::DumpAt((name), (object), ::Stereon::SourceLocation::Current()))

#endif // STEREON_IO_DUMP_HPP
