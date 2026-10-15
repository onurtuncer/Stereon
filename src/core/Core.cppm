// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// Stereon.Core: math types, intervals, boxes, handles, arenas, errors, the
/// Context with its log and journal sinks, plus the Stereon::Simd polyfill
/// wrapper (ADR-0016). Stereon::Exec is a
/// header, <Stereon/Core/Exec.hpp>, because of compiler bugs (see there).
export module Stereon.Core;

export import :Allowances;
export import :Arena;
export import :Box;
export import :Context;
export import :Enum;
export import :EvaluatorChecks;
export import :Error;
export import :Frame;
export import :GeometryChecks;
export import :Handle;
export import :Interval;
export import :Journal;
export import :Log;
export import :Matrix;
export import :Scalar;
export import :Serialize;
export import :Simd;
export import :SmallBuffer;
export import :Trigonometry;
export import :Vector;
export import :Version;
