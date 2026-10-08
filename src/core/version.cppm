// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

export module stereon.core:version;

export namespace stereon {

/// Library version, following semantic versioning.
struct Version {
    int major;
    int minor;
    int patch;

    // Member, not hidden friend: a defaulted friend operator== exported from a
    // module crashes GCC 16.2 (docs/toolchain/gcc16-feature-probe.md).
    constexpr bool operator==(const Version&) const = default;
};

/// Set from the CMake project version.
inline constexpr Version version{STEREON_VERSION_MAJOR, STEREON_VERSION_MINOR,
                                 STEREON_VERSION_PATCH};

/// True when the library was built with checked handles (ADR-0004).
inline constexpr bool checked_handles =
#if defined(STEREON_CHECKED_HANDLES)
    true;
#else
    false;
#endif

} // namespace stereon
