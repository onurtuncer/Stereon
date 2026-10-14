// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <array>
#include <cstddef>
#include <functional>
#include <memory>
#include <new>
#include <optional>
#include <type_traits>
#include <utility>

export module Stereon.Core:SmallBuffer;

/// Storage for the type-erased geometry wrappers of ADR-0003: a model object
/// placed in a fixed buffer inside the wrapper, so that wrapping a small
/// geometry type allocates nothing. A geometry type too large for the buffer
/// is held through a shared pointer to an immutable value instead, so that
/// copies of the wrapper share it and stay cheap.
///
/// The interface a wrapper erases derives from ErasedConcept, and its models
/// hold the value in an InlineValue or a SharedValue. SmallBuffer copies,
/// moves and destroys models through the concept's virtual functions.

export namespace Stereon {

/// Bytes of in-place storage in Curve, Curve2d and Surface (ADR-0003). A
/// tuning constant fixed by the copy benchmark in bench/GeometryBench.cpp;
/// not part of the API or ABI. 128 holds Cone and Torus (96 bytes plus the
/// model's table pointer), whose copies cost about 15 ns inline against 43 ns
/// shared, for about 2.5 ns more on a Line or Circle copy than at 96.
inline constexpr std::size_t GeometryBufferSize = 128;

/// A list of types, for the closed sets of built-in geometry that
/// VisitBuiltin and the intersection tables are generated over.
template <typename... Ts>
struct TypeList {
    static constexpr std::size_t Size = sizeof...(Ts);
};

} // namespace Stereon

export namespace Stereon::Detail {

/// The base of every erased interface: the operations SmallBuffer needs.
template <typename Concept>
struct ErasedConcept {
    ErasedConcept() = default;
    ErasedConcept(const ErasedConcept&) = default;
    ErasedConcept& operator=(const ErasedConcept&) = delete;
    virtual ~ErasedConcept() = default;

    /// A copy of the model constructed at `buffer`.
    virtual Concept* CloneInto(void* buffer) const = 0;
    /// The model moved to `buffer`; this one is left to be destroyed.
    virtual Concept* MoveInto(void* buffer) noexcept = 0;
    /// Identifies the erased value's type; see TypeIdOf.
    [[nodiscard]] virtual const void* GetTypeId() const noexcept = 0;
    /// The erased value, for Target<T>().
    [[nodiscard]] virtual const void* GetValue() const noexcept = 0;
};

/// One distinct address per type, so that the erased type can be compared
/// without RTTI. Inline variables have one definition in the program.
template <typename T>
inline constexpr char TypeTag = 0;

template <typename T>
[[nodiscard]] constexpr const void* TypeIdOf() noexcept
{
    return &TypeTag<T>;
}

/// A geometry value stored in the buffer.
template <typename T>
class InlineValue {
public:
    explicit InlineValue(T value)
        : m_Value(std::move(value))
    {}
    [[nodiscard]] const T& Get() const noexcept { return m_Value; }

private:
    T m_Value;
};

/// A geometry value too large for the buffer, held immutable and shared
/// between copies.
template <typename T>
class SharedValue {
public:
    explicit SharedValue(T value)
        : m_Value(std::make_shared<const T>(std::move(value)))
    {}
    [[nodiscard]] const T& Get() const noexcept { return *m_Value; }

private:
    std::shared_ptr<const T> m_Value;
};

/// True when a model holding T in place, plus its virtual table pointer,
/// fits in `Size` bytes and T is trivially copyable. A type whose copy
/// allocates, such as a B-spline with its pole and knot vectors, is shared
/// however small its own footprint, so that copying a wrapper never copies
/// its geometry (ADR-0003: "larger types store an immutable shared payload
/// so copies stay cheap").
template <typename T, std::size_t Size>
inline constexpr bool FitsInline =
    sizeof(T) + sizeof(void*) <= Size && alignof(T) <= alignof(std::max_align_t) &&
    std::is_trivially_copyable_v<T>;

/// InlineValue<T> when it fits, SharedValue<T> otherwise.
template <typename T, std::size_t Size>
using ValueHolder = std::conditional_t<FitsInline<T, Size>, InlineValue<T>, SharedValue<T>>;

/// The common model plumbing: the holder, the clone and move operations and
/// the type identity. Models of a concrete interface derive from it.
template <typename Concept, typename T, std::size_t Size>
class ModelBase : public Concept {
public:
    explicit ModelBase(T value)
        : m_Holder(std::move(value))
    {}

    Concept* CloneInto(void* buffer) const override
    {
        return ::new (buffer) ModelOf(static_cast<const ModelOf&>(*this));
    }
    Concept* MoveInto(void* buffer) noexcept override
    {
        return ::new (buffer) ModelOf(std::move(static_cast<ModelOf&>(*this)));
    }
    [[nodiscard]] const void* GetTypeId() const noexcept override { return TypeIdOf<T>(); }
    [[nodiscard]] const void* GetValue() const noexcept override { return &Get(); }

protected:
    /// The most derived model, which the concrete interface defines.
    using ModelOf = typename Concept::template Model<T, Size>;

    [[nodiscard]] const T& Get() const noexcept { return m_Holder.Get(); }

private:
    ValueHolder<T, Size> m_Holder;
};

/// A model of `Concept` placed in `Size` bytes, or nothing.
template <typename Concept, std::size_t Size>
class SmallBuffer {
public:
    /// Empty. User-provided, so that a const empty wrapper can be declared:
    /// the storage is raw bytes and stays uninitialised until Emplace().
    SmallBuffer() noexcept {} // NOLINT(modernize-use-equals-default)

    SmallBuffer(const SmallBuffer& other)
    {
        if (other.m_Model != nullptr) {
            m_Model = other.m_Model->CloneInto(&m_Storage);
        }
    }
    SmallBuffer(SmallBuffer&& other) noexcept
    {
        if (other.m_Model != nullptr) {
            m_Model = other.m_Model->MoveInto(&m_Storage);
            other.Reset();
        }
    }
    SmallBuffer& operator=(const SmallBuffer& other)
    {
        if (this != &other) {
            SmallBuffer copy(other);
            *this = std::move(copy);
        }
        return *this;
    }
    SmallBuffer& operator=(SmallBuffer&& other) noexcept
    {
        if (this != &other) {
            Reset();
            if (other.m_Model != nullptr) {
                m_Model = other.m_Model->MoveInto(&m_Storage);
                other.Reset();
            }
        }
        return *this;
    }
    ~SmallBuffer() { Reset(); }

    /// Places a model of type M, constructed from `args`.
    template <typename M, typename... Args>
    void Emplace(Args&&... args)
    {
        static_assert(sizeof(M) <= Size, "the model does not fit the buffer");
        static_assert(alignof(M) <= alignof(std::max_align_t));
        Reset();
        m_Model = ::new (&m_Storage) M(std::forward<Args>(args)...);
    }

    [[nodiscard]] bool IsEmpty() const noexcept { return m_Model == nullptr; }
    [[nodiscard]] const Concept* Get() const noexcept { return m_Model; }

private:
    void Reset() noexcept
    {
        if (m_Model != nullptr) {
            m_Model->~Concept();
            m_Model = nullptr;
        }
    }

    alignas(std::max_align_t) std::byte m_Storage[Size];
    Concept* m_Model = nullptr;
};

/// 0, 1, ..., N - 1, for expansion statements over a type pack.
template <std::size_t N>
inline constexpr std::array<std::size_t, N> Indices = [] {
    std::array<std::size_t, N> indices{};
    for (std::size_t i = 0; i < N; ++i) {
        indices[i] = i;
    }
    return indices;
}();

/// The position of the wrapped value's type in Ts..., or Size when it is not
/// one of them. One virtual call; the comparisons are unrolled.
///
/// Both functions below are generated with an expansion statement and pack
/// indexing (ADR-0016) where the compiler has them, and with a fold
/// expression otherwise: MSVC and Clang lack expansion statements so far,
/// and they build in non-blocking CI. The two forms make the same
/// comparisons in the same order.
template <typename Wrapper, typename... Ts>
[[nodiscard]] std::size_t IndexIn(TypeList<Ts...>, const Wrapper& wrapper) noexcept
{
    const void* id = wrapper.GetModel().GetTypeId();
#if defined(__cpp_expansion_statements)
    template for (constexpr std::size_t i : Indices<sizeof...(Ts)>)
    {
        if (id == TypeIdOf<Ts...[i]>()) {
            return i;
        }
    }
    return sizeof...(Ts);
#else
    std::size_t index = 0;
    const bool found = ((id == TypeIdOf<Ts>() || (++index, false)) || ...);
    return found ? index : sizeof...(Ts);
#endif
}

/// Calls f with the wrapped value as its own type when that is one of Ts...,
/// and with the wrapper otherwise (ADR-0003). The type test runs once, so a
/// loop inside f runs on the concrete type with no virtual calls. Every call
/// of f must return the same type, by value.
template <typename Wrapper, typename F, typename... Ts>
decltype(auto) VisitErased(TypeList<Ts...>, const Wrapper& wrapper, F&& f)
{
    const auto& model = wrapper.GetModel();
    const void* id = model.GetTypeId();
#if defined(__cpp_expansion_statements)
    template for (constexpr std::size_t i : Indices<sizeof...(Ts)>)
    {
        using T = Ts...[i];
        if (id == TypeIdOf<T>()) {
            return std::invoke(f, *static_cast<const T*>(model.GetValue()));
        }
    }
    return std::invoke(f, wrapper);
#else
    using Result = std::invoke_result_t<F&, const Wrapper&>;
    if constexpr (std::is_void_v<Result>) {
        const bool done = ((id == TypeIdOf<Ts>() &&
                            (std::invoke(f, *static_cast<const Ts*>(model.GetValue())), true)) ||
                           ...);
        if (!done) {
            std::invoke(f, wrapper);
        }
    } else {
        std::optional<Result> result;
        const bool done =
            ((id == TypeIdOf<Ts>() &&
              (result.emplace(std::invoke(f, *static_cast<const Ts*>(model.GetValue()))), true)) ||
             ...);
        if (!done) {
            return Result(std::invoke(f, wrapper));
        }
        return Result(std::move(*result));
    }
#endif
}

} // namespace Stereon::Detail
