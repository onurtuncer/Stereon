// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <cstdint>
#include <expected>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

export module Stereon.Core:Error;

/// Error handling (ADR-0006). Operations that can fail return
/// Result<T> = std::expected<T, KernelError>; the kernel does not throw,
/// although std::bad_alloc may still propagate. Invalid data from outside the
/// kernel is a KernelError; a broken programmer obligation is a contract
/// violation (ADR-0016), never a KernelError.
///
/// STEREON_TRY and STEREON_TRY_VOID, in <Stereon/Core/Try.hpp>, propagate
/// errors from one Result to the caller's.
export namespace Stereon {

/// Why an operation failed. The C ABI exposes the numbers (ADR-0020), so
/// values never change: new codes are appended, and 0 is left for success.
enum class ErrorCode : std::uint32_t {
    InvalidInput = 1,    ///< Input data failed validation.
    InvalidResult,       ///< The result failed its validity post-check (ADR-0018).
    ToleranceExceeded,   ///< A tolerance would grow above the ceiling (ADR-0001).
    IntersectionFailed,  ///< An intersection could not be computed.
    CertificationFailed, ///< A deviation bound could not be certified (ADR-0009).
    DegenerateResult,    ///< The result would be degenerate.
    NotSupported,        ///< The input is valid but outside what the kernel handles.
    NotImplemented,      ///< The case is planned but not yet implemented.
    ParseError,          ///< A file or stream is malformed.
};

/// The enumerator's name, e.g. "ToleranceExceeded"; "Unknown" for a value
/// that is not an enumerator.
[[nodiscard]] constexpr std::string_view ToString(const ErrorCode code) noexcept
{
    switch (code) {
    case ErrorCode::InvalidInput:
        return "InvalidInput";
    case ErrorCode::InvalidResult:
        return "InvalidResult";
    case ErrorCode::ToleranceExceeded:
        return "ToleranceExceeded";
    case ErrorCode::IntersectionFailed:
        return "IntersectionFailed";
    case ErrorCode::CertificationFailed:
        return "CertificationFailed";
    case ErrorCode::DegenerateResult:
        return "DegenerateResult";
    case ErrorCode::NotSupported:
        return "NotSupported";
    case ErrorCode::NotImplemented:
        return "NotImplemented";
    case ErrorCode::ParseError:
        return "ParseError";
    }
    return "Unknown";
}

/// A failed operation: what went wrong, a message for people and the
/// tolerance values at stake.
///
/// ADR-0006 also attaches the entities involved, a ValidityReport and a
/// journal ID for replay. They arrive with handles, the validity checker
/// and the operation journal.
class KernelError {
public:
    KernelError(const ErrorCode code, std::string message)
        : m_Code(code)
        , m_Message(std::move(message))
    {}

    [[nodiscard]] ErrorCode GetCode() const noexcept { return m_Code; }
    [[nodiscard]] std::string_view GetMessage() const noexcept { return m_Message; }

    /// Tolerances at stake in millimetres, in the order they were added.
    [[nodiscard]] std::span<const double> GetTolerances() const noexcept { return m_Tolerances; }

    /// Adds a tolerance value; chains on a temporary:
    /// KernelError(ErrorCode::ToleranceExceeded, "...").WithTolerance(0.02).
    [[nodiscard]] KernelError WithTolerance(const double tolerance) &&
    {
        m_Tolerances.push_back(tolerance);
        return std::move(*this);
    }

private:
    ErrorCode m_Code;
    std::string m_Message;
    std::vector<double> m_Tolerances;
};

/// The return type of every public operation that can fail (ADR-0006).
template <typename T>
using Result = std::expected<T, KernelError>;

/// An error to return from a function that returns Result:
/// return Fail(ErrorCode::InvalidInput, "radius must be positive");
[[nodiscard]] inline std::unexpected<KernelError> Fail(const ErrorCode code, std::string message)
{
    return std::unexpected<KernelError>(std::in_place, code, std::move(message));
}

} // namespace Stereon
