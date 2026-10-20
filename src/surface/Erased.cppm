// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <Stereon/Core/Contracts.hpp>

#include <concepts>
#include <cstddef>
#include <string_view>
#include <type_traits>
#include <utility>

export module Stereon.Surface:Erased;

import Stereon.Curve;
import :Concepts;

/// Surface: any type satisfying SurfaceGeometry behind one copyable value
/// type (ADR-0003), stored like Curve: in place when small, shared and
/// immutable when large.

export namespace Stereon {

class Surface;

} // namespace Stereon

export namespace Stereon::Detail {

template <typename T, std::size_t Size>
class SurfaceModel;

/// The erased interface; member names differ from the free functions, as in
/// CurveConcept.
struct SurfaceConcept : ErasedConcept<SurfaceConcept> {
    template <typename T, std::size_t Size>
    using Model = SurfaceModel<T, Size>;

    [[nodiscard]] virtual Point3 EvalAt(const Point2& uv) const = 0;
    [[nodiscard]] virtual Vec3 DerivAt(const Point2& uv, int du, int dv) const = 0;
    [[nodiscard]] virtual Box2 GetDomain() const = 0;
    /// Periods(value) for a PeriodicSurface, (0, 0) otherwise.
    [[nodiscard]] virtual Vec2 GetPeriods() const = 0;
    [[nodiscard]] virtual Box3 GetBoundingBox() const = 0;
    [[nodiscard]] virtual Point2 ProjectPoint(const Point3& p) const = 0;
    [[nodiscard]] virtual Surface MakeReversed() const = 0;
    [[nodiscard]] virtual Surface MakeTransformed(const Transform& x) const = 0;
};

} // namespace Stereon::Detail

export namespace Stereon {

class Surface {
public:
    /// The empty surface; nothing may be evaluated on it.
    Surface() = default;

    /// Wraps a copy of `geometry`.
    template <typename G>
        requires(!std::same_as<std::remove_cvref_t<G>, Surface>) &&
                SurfaceGeometry<std::remove_cvref_t<G>>
    explicit Surface(G&& geometry)
    {
        using Model = Detail::SurfaceModel<std::remove_cvref_t<G>, SurfaceBufferSize>;
        m_Storage.Emplace<Model>(std::forward<G>(geometry));
    }

    [[nodiscard]] bool IsEmpty() const noexcept { return m_Storage.IsEmpty(); }

    /// The wrapped value if it is a T, otherwise null; see Curve::Target.
    template <typename T>
    [[nodiscard]] const T* Target() const noexcept
    {
        const Detail::SurfaceConcept* model = m_Storage.Get();
        if (model == nullptr || model->GetTypeId() != Detail::TypeIdOf<T>()) {
            return nullptr;
        }
        return static_cast<const T*>(model->GetValue());
    }

    [[nodiscard]] const Detail::SurfaceConcept& GetModel() const STEREON_HOT_PRE(!IsEmpty())
    {
        return *m_Storage.Get();
    }

private:
    Detail::SmallBuffer<Detail::SurfaceConcept, SurfaceBufferSize> m_Storage;
};

} // namespace Stereon

namespace Stereon::Detail {

// The checked-build preconditions of the Surface functions below, defined in
// Surface.cpp. Out of line, because instantiating core's check templates in
// this partition's inline functions makes GCC 16.2 write a module file that
// cannot be read back ("Bad file data"; docs/toolchain/gcc16-feature-probe.md).
void CheckSurfaceParameter(const SurfaceConcept& model, const Point2& uv,
                           std::string_view function) noexcept;
void CheckDerivativeOrder(int du, int dv, std::string_view function) noexcept;
void CheckSurfacePoint(const Point3& p, std::string_view function) noexcept;

} // namespace Stereon::Detail

export namespace Stereon {

[[nodiscard]] inline Point3 Eval(const Surface& surface, const Point2& uv)
{
    if constexpr (Detail::EvaluatorChecks) {
        Detail::CheckSurfaceParameter(surface.GetModel(), uv, "Eval(Surface)");
    }
    return surface.GetModel().EvalAt(uv);
}

[[nodiscard]] inline Vec3 Deriv(const Surface& surface, const Point2& uv, const int du,
                                const int dv)
{
    if constexpr (Detail::EvaluatorChecks) {
        Detail::CheckSurfaceParameter(surface.GetModel(), uv, "Deriv(Surface)");
        Detail::CheckDerivativeOrder(du, dv, "Deriv(Surface)");
    }
    return surface.GetModel().DerivAt(uv, du, dv);
}

/// The wrapped surface's periods in u and v, 0 where it is not periodic.
[[nodiscard]] inline Vec2 Periods(const Surface& surface)
{
    return surface.GetModel().GetPeriods();
}

[[nodiscard]] inline Box2 Domain(const Surface& surface)
{
    return surface.GetModel().GetDomain();
}

[[nodiscard]] inline Box3 BoundingBox(const Surface& surface)
{
    return surface.GetModel().GetBoundingBox();
}

[[nodiscard]] inline Point2 Project(const Surface& surface, const Point3& p)
{
    if constexpr (Detail::EvaluatorChecks) {
        Detail::CheckSurfacePoint(p, "Project(Surface)");
    }
    return surface.GetModel().ProjectPoint(p);
}

[[nodiscard]] inline Surface Reversed(const Surface& surface)
{
    return surface.GetModel().MakeReversed();
}

[[nodiscard]] inline Surface Transformed(const Surface& surface, const Transform& x)
{
    return surface.GetModel().MakeTransformed(x);
}

} // namespace Stereon

export namespace Stereon::Detail {

template <typename T, std::size_t Size>
class SurfaceModel final : public ModelBase<SurfaceConcept, T, Size> {
    using Base = ModelBase<SurfaceConcept, T, Size>;

public:
    using Base::Base;

    [[nodiscard]] Point3 EvalAt(const Point2& uv) const override { return Eval(this->Get(), uv); }
    [[nodiscard]] Vec3 DerivAt(const Point2& uv, const int du, const int dv) const override
    {
        return Deriv(this->Get(), uv, du, dv);
    }
    [[nodiscard]] Box2 GetDomain() const override { return Domain(this->Get()); }
    [[nodiscard]] Vec2 GetPeriods() const override
    {
        if constexpr (PeriodicSurface<T>) {
            return Periods(this->Get());
        } else {
            return {};
        }
    }
    [[nodiscard]] Box3 GetBoundingBox() const override { return BoundingBox(this->Get()); }
    [[nodiscard]] Point2 ProjectPoint(const Point3& p) const override
    {
        return Project(this->Get(), p);
    }
    [[nodiscard]] Surface MakeReversed() const override { return Surface(Reversed(this->Get())); }
    [[nodiscard]] Surface MakeTransformed(const Transform& x) const override
    {
        return Surface(Transformed(this->Get(), x));
    }
};

} // namespace Stereon::Detail

export namespace Stereon {

static_assert(SurfaceGeometry<Surface>);
static_assert(std::semiregular<Surface>);
// The buffer and the model pointer, nothing more (B0.4).
static_assert(sizeof(Surface) == SurfaceBufferSize + sizeof(void*));

} // namespace Stereon
