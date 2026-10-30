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

export module Stereon.Curve:Erased;

import Stereon.Robust;
import :Concepts;

/// Curve and Curve2d: any type satisfying CurveGeometry or Curve2dGeometry,
/// behind one copyable value type (ADR-0003). A small geometry type lives in
/// the wrapper's own buffer; a large one is shared, immutable, between
/// copies. Calls go through one virtual function each. Hot loops over
/// built-in types will dispatch once with VisitBuiltin instead (a later
/// step of TODO section 6).
///
/// Curve and Curve2d are double only (ADR-0021): scalar-generic code works on
/// the built-in types directly.

export namespace Stereon {

template <std::size_t D>
class BasicErasedCurve;

/// A curve in space, of any CurveGeometry type.
using Curve = BasicErasedCurve<3>;
/// A curve in a surface's parameter plane, of any Curve2dGeometry type.
using Curve2d = BasicErasedCurve<2>;

} // namespace Stereon

export namespace Stereon::Detail {

template <std::size_t D, typename T, std::size_t Size>
class CurveModel;

/// The erased interface. Member names differ from the free functions, so that
/// a model's call to Eval(value, t) finds the geometry type's function by
/// argument-dependent lookup rather than the member.
template <std::size_t D>
struct CurveConcept : ErasedConcept<CurveConcept<D>> {
    template <typename T, std::size_t Size>
    using Model = CurveModel<D, T, Size>;

    [[nodiscard]] virtual BasicPointN<double, D> EvalAt(double t) const = 0;
    [[nodiscard]] virtual BasicVecN<double, D> DerivAt(double t, int n) const = 0;
    [[nodiscard]] virtual Interval GetDomain() const = 0;
    /// Period(value) for a PeriodicCurve, 0 otherwise.
    [[nodiscard]] virtual double GetPeriod() const = 0;
    [[nodiscard]] virtual BoxN<D> GetBoundingBox() const = 0;
    [[nodiscard]] virtual double ProjectPoint(const BasicPointN<double, D>& p) const = 0;
    [[nodiscard]] virtual BasicErasedCurve<D> MakeReversed() const = 0;
    /// Called only for D == 3; see Transformed(const Curve&, ...).
    [[nodiscard]] virtual BasicErasedCurve<D> MakeTransformed(const Transform& x) const = 0;
};

/// The geometry concept a type must satisfy to be wrapped in dimension D.
template <typename G, std::size_t D>
concept WrappableCurve = (D == 2 && Curve2dGeometry<G>) || (D == 3 && CurveGeometry<G>);

} // namespace Stereon::Detail

export namespace Stereon {

template <std::size_t D>
class BasicErasedCurve {
public:
    /// The empty curve; nothing may be evaluated on it.
    BasicErasedCurve() = default;

    /// Wraps a copy of `geometry`.
    template <typename G>
        requires(!std::same_as<std::remove_cvref_t<G>, BasicErasedCurve>) &&
                Detail::WrappableCurve<std::remove_cvref_t<G>, D>
    explicit BasicErasedCurve(G&& geometry)
    {
        using Model = Detail::CurveModel<D, std::remove_cvref_t<G>, CurveBufferSizeN<D>>;
        m_Storage.template Emplace<Model>(std::forward<G>(geometry));
    }

    [[nodiscard]] bool IsEmpty() const noexcept { return m_Storage.IsEmpty(); }

    /// The wrapped value if it is a T, otherwise null. A borrow of this
    /// wrapper (ADR-0046): the pointer is into the wrapper's own buffer, or
    /// into the value the copies of a large type share, and stays valid until
    /// the wrapper is assigned, moved from or destroyed. Lvalues only: a
    /// temporary curve would leave it dangling, so name the curve first.
    /// Copies of a curve whose type is too large for the buffer share one
    /// value, so their Target<T>() pointers are equal.
    template <typename T>
    [[nodiscard]] const T* Target() const& noexcept
    {
        const Detail::CurveConcept<D>* model = m_Storage.Get();
        if (model == nullptr || model->GetTypeId() != Detail::TypeIdOf<T>()) {
            return nullptr;
        }
        return static_cast<const T*>(model->GetValue());
    }
    template <typename T>
    const T* Target() const&& = delete;

    /// The model, for the free functions below: a borrow of this wrapper,
    /// as Target<T>() is, and lvalue-only for the same reason.
    [[nodiscard]] const Detail::CurveConcept<D>& GetModel() const& STEREON_HOT_PRE(!IsEmpty())
    {
        return *m_Storage.Get();
    }
    const Detail::CurveConcept<D>& GetModel() const&& = delete;

private:
    Detail::SmallBuffer<Detail::CurveConcept<D>, CurveBufferSizeN<D>> m_Storage;
};

} // namespace Stereon

namespace Stereon::Detail {

/// The checks every erased curve evaluation makes in checked builds: a
/// finite parameter inside the domain of a non-periodic curve. They run
/// here, at the boundary to user geometry, because a user type need not
/// check anything itself.
template <std::size_t D>
void CheckCurveParameter(const CurveConcept<D>& model, const double t,
                         const std::string_view function) noexcept
{
    if constexpr (EvaluatorChecks) {
        RequireFinite(t, function, "a non-finite parameter");
        RequireInDomain(model.GetDomain(), t, model.GetPeriod(), function);
    }
}

} // namespace Stereon::Detail

export namespace Stereon {

template <std::size_t D>
[[nodiscard]] BasicPointN<double, D> Eval(const BasicErasedCurve<D>& curve, const double t)
{
    Detail::CheckCurveParameter(curve.GetModel(), t, "Eval(Curve)");
    return curve.GetModel().EvalAt(t);
}

template <std::size_t D>
[[nodiscard]] BasicVecN<double, D> Deriv(const BasicErasedCurve<D>& curve, const double t,
                                         const int n)
{
    Detail::CheckCurveParameter(curve.GetModel(), t, "Deriv(Curve)");
    Detail::Require(n >= 1, "Deriv(Curve)", "a derivative order below 1");
    return curve.GetModel().DerivAt(t, n);
}

/// The wrapped curve's period, or 0 when it is not periodic.
template <std::size_t D>
[[nodiscard]] double Period(const BasicErasedCurve<D>& curve)
{
    return curve.GetModel().GetPeriod();
}

template <std::size_t D>
[[nodiscard]] Interval Domain(const BasicErasedCurve<D>& curve)
{
    return curve.GetModel().GetDomain();
}

template <std::size_t D>
[[nodiscard]] BoxN<D> BoundingBox(const BasicErasedCurve<D>& curve)
{
    return curve.GetModel().GetBoundingBox();
}

template <std::size_t D>
[[nodiscard]] double Project(const BasicErasedCurve<D>& curve, const BasicPointN<double, D>& p)
{
    Detail::RequireFinite(p, "Project(Curve)", "a non-finite point");
    return curve.GetModel().ProjectPoint(p);
}

template <std::size_t D>
[[nodiscard]] BasicErasedCurve<D> Reversed(const BasicErasedCurve<D>& curve)
{
    return curve.GetModel().MakeReversed();
}

template <std::size_t D>
    requires(D == 3)
[[nodiscard]] BasicErasedCurve<D> Transformed(const BasicErasedCurve<D>& curve, const Transform& x)
{
    return curve.GetModel().MakeTransformed(x);
}

} // namespace Stereon

export namespace Stereon::Detail {

template <std::size_t D, typename T, std::size_t Size>
class CurveModel final : public ModelBase<CurveConcept<D>, T, Size> {
    using Base = ModelBase<CurveConcept<D>, T, Size>;

public:
    using Base::Base;

    [[nodiscard]] BasicPointN<double, D> EvalAt(const double t) const override
    {
        return Eval(this->Get(), t);
    }
    [[nodiscard]] BasicVecN<double, D> DerivAt(const double t, const int n) const override
    {
        return Deriv(this->Get(), t, n);
    }
    [[nodiscard]] Interval GetDomain() const override { return Domain(this->Get()); }
    [[nodiscard]] double GetPeriod() const override
    {
        if constexpr (PeriodicCurve<T>) {
            return Period(this->Get());
        } else {
            return 0.0;
        }
    }
    [[nodiscard]] BoxN<D> GetBoundingBox() const override { return BoundingBox(this->Get()); }
    [[nodiscard]] double ProjectPoint(const BasicPointN<double, D>& p) const override
    {
        return Project(this->Get(), p);
    }
    [[nodiscard]] BasicErasedCurve<D> MakeReversed() const override
    {
        return BasicErasedCurve<D>(Reversed(this->Get()));
    }
    [[nodiscard]] BasicErasedCurve<D> MakeTransformed(const Transform& x) const override
    {
        if constexpr (D == 3) {
            return BasicErasedCurve<D>(Transformed(this->Get(), x));
        } else {
            // Curve2d has no Transformed, so nothing calls this.
            static_cast<void>(x);
            std::unreachable();
        }
    }
};

} // namespace Stereon::Detail

export namespace Stereon {

static_assert(CurveGeometry<Curve>);
static_assert(Curve2dGeometry<Curve2d>);
static_assert(std::semiregular<Curve>);
// The buffer plus the model pointer, padded to the buffer's alignment.
// The buffer and the model pointer, nothing more (B0.4).
static_assert(sizeof(Curve) == CurveBufferSize + sizeof(void*));
static_assert(sizeof(Curve2d) == Curve2dBufferSize + sizeof(void*));

} // namespace Stereon
