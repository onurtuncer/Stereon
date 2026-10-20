// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// Stereon.Io: STEP, STL, OBJ, glTF and the native format (ADR-0012). So
/// far the geometry subset of the native format (:GeometryArchive), its
/// debug chunk (:Debug) and dumps for the debug viewer (:Dump, ADR-0035).
export module Stereon.Io;

export import Stereon.Mesh;
export import :Debug;
export import :Dump;
export import :GeometryArchive;
export import :Records;
