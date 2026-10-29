// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Detail::SmallBuffer, the storage behind the erased geometry wrappers
// (ADR-0003), against ADR-0046's rules for manual object lifetime (section
// 6, AC-0046-06): the empty state, inline and shared payloads, the alignment
// limit, copy, move and self-assignment, a model whose construction throws,
// and exactly one destruction per model. The models count their live
// instances, so a leak or a double destruction shows as a count that does
// not return to zero.

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

import Stereon.Core;

namespace {

using Stereon::Detail::ErasedConcept;
using Stereon::Detail::ModelBase;
using Stereon::Detail::SmallBuffer;

/// How many models are alive, over every type below.
int liveModels = 0;

/// A small erased interface: a value, and the plumbing ErasedConcept asks for.
struct Numbered : ErasedConcept<Numbered> {
    template <typename T, std::size_t Size>
    class Model;

    [[nodiscard]] virtual int Value() const noexcept = 0;
};

template <typename T, std::size_t Size>
class Numbered::Model final : public ModelBase<Numbered, T, Size> {
public:
    explicit Model(T value)
        : ModelBase<Numbered, T, Size>(std::move(value))
    {
        ++liveModels;
    }
    Model(const Model& other)
        : ModelBase<Numbered, T, Size>(other)
    {
        ++liveModels;
    }
    Model(Model&& other) noexcept
        : ModelBase<Numbered, T, Size>(std::move(other))
    {
        ++liveModels;
    }
    Model& operator=(const Model&) = delete;
    Model& operator=(Model&&) = delete;
    ~Model() override { --liveModels; }

    [[nodiscard]] int Value() const noexcept override { return this->Get().Number; }
};

/// A trivially copyable value that fits inline.
struct Small {
    int Number = 0;
};

/// A value too large for the buffer, so that its model shares it.
struct Large {
    int Number = 0;
    std::array<double, 64> Padding{};
};

/// A value whose copies can be told to fail: building it never throws, so
/// the throw comes from inside Emplace(), where the model copies or moves
/// the value into place, not from the caller's argument.
struct Fragile {
    Fragile(const int number, const bool fail) noexcept
        : Number(number)
        , Fail(fail)
    {}
    Fragile(const Fragile& other)
        : Number(other.Number)
        , Fail(other.Fail)
    {
        if (Fail) {
            throw std::runtime_error("Fragile refused to be copied");
        }
    }
    Fragile(Fragile&& other) // NOLINT(performance-noexcept-move-constructor): it throws on purpose
        : Fragile(std::as_const(other))
    {}
    Fragile& operator=(const Fragile&) = default;
    Fragile& operator=(Fragile&&) = default;
    ~Fragile() = default;

    int Number;
    bool Fail;
};

constexpr std::size_t BufferSize = 64;
using Buffer = SmallBuffer<Numbered, BufferSize>;
using SmallModel = Numbered::Model<Small, BufferSize>;
using LargeModel = Numbered::Model<Large, BufferSize>;
using FragileModel = Numbered::Model<Fragile, BufferSize>;

int ValueOf(const Buffer& buffer)
{
    return buffer.Get()->Value();
}

/// Whether Get() is callable on an lvalue or an rvalue buffer, as templates
/// so that the rejected call is a false requirement, not an error.
template <typename B>
constexpr bool GetFromLvalue = requires(const B& b) { b.Get(); };
template <typename B>
constexpr bool GetFromRvalue = requires(B&& b) { std::move(b).Get(); };

} // namespace

TEST_CASE("an empty SmallBuffer holds nothing and may be copied and moved", "[core][erased]")
{
    liveModels = 0;
    {
        const Buffer empty;
        REQUIRE(empty.IsEmpty());
        REQUIRE(empty.Get() == nullptr);
        Buffer copy(empty);
        REQUIRE(copy.IsEmpty());
        Buffer moved(std::move(copy));
        REQUIRE(moved.IsEmpty());
        moved = empty;
        REQUIRE(moved.IsEmpty());
    }
    REQUIRE(liveModels == 0);
}

TEST_CASE("inline and shared payloads are placed, copied, moved and destroyed once",
          "[core][erased]")
{
    STATIC_REQUIRE(Stereon::Detail::FitsInline<Small, BufferSize>);
    STATIC_REQUIRE_FALSE(Stereon::Detail::FitsInline<Large, BufferSize>);
    STATIC_REQUIRE(sizeof(SmallModel) <= BufferSize);
    STATIC_REQUIRE(sizeof(LargeModel) <= BufferSize); // the model holds a shared_ptr, not Large
    STATIC_REQUIRE(alignof(SmallModel) <= Stereon::GeometryAlignment);
    STATIC_REQUIRE(alignof(LargeModel) <= Stereon::GeometryAlignment);

    liveModels = 0;
    {
        Buffer small;
        small.Emplace<SmallModel>(Small{7});
        Buffer large;
        large.Emplace<LargeModel>(Large{11, {}});
        REQUIRE(liveModels == 2);
        REQUIRE(ValueOf(small) == 7);
        REQUIRE(ValueOf(large) == 11);

        // Copies hold models of their own; the large one shares the value.
        const Buffer smallCopy = small;
        const Buffer largeCopy = large;
        REQUIRE(liveModels == 4);
        REQUIRE(ValueOf(smallCopy) == 7);
        REQUIRE(ValueOf(largeCopy) == 11);
        REQUIRE(smallCopy.Get() != small.Get()); // each in its own buffer
        REQUIRE(largeCopy.Get() != large.Get());

        // A move leaves the source empty and the model count unchanged.
        Buffer movedSmall = std::move(small);
        Buffer movedLarge = std::move(large);
        REQUIRE(small.IsEmpty()); // NOLINT(bugprone-use-after-move): the point of the test
        REQUIRE(large.IsEmpty()); // NOLINT(bugprone-use-after-move)
        REQUIRE(liveModels == 4);
        REQUIRE(ValueOf(movedSmall) == 7);
        REQUIRE(ValueOf(movedLarge) == 11);

        // Assignment replaces the model, destroying the old one once.
        movedSmall = largeCopy;
        REQUIRE(liveModels == 4);
        REQUIRE(ValueOf(movedSmall) == 11);
        movedLarge = std::move(movedSmall);
        REQUIRE(movedSmall.IsEmpty()); // NOLINT(bugprone-use-after-move)
        REQUIRE(liveModels == 3);
        REQUIRE(ValueOf(movedLarge) == 11);

        // Re-placing destroys the previous model first.
        movedLarge.Emplace<SmallModel>(Small{3});
        REQUIRE(liveModels == 3);
        REQUIRE(ValueOf(movedLarge) == 3);
    }
    REQUIRE(liveModels == 0);
}

TEST_CASE("self-assignment keeps the model", "[core][erased]")
{
    liveModels = 0;
    {
        Buffer buffer;
        buffer.Emplace<SmallModel>(Small{5});
        Buffer& alias = buffer;
        buffer = alias;
        REQUIRE(ValueOf(buffer) == 5);
        buffer = std::move(alias);
        REQUIRE(ValueOf(buffer) == 5);
        REQUIRE(liveModels == 1);
    }
    REQUIRE(liveModels == 0);
}

TEST_CASE("a model whose construction throws leaves the buffer empty and leaks nothing",
          "[core][erased]")
{
    liveModels = 0;
    {
        Buffer buffer;
        buffer.Emplace<SmallModel>(Small{9});
        REQUIRE(liveModels == 1);
        REQUIRE_THROWS_AS(buffer.Emplace<FragileModel>(Fragile{1, true}), std::runtime_error);
        // The previous model was destroyed before the failed construction;
        // nothing was counted for the one that did not come to be.
        REQUIRE(buffer.IsEmpty());
        REQUIRE(buffer.Get() == nullptr);
        REQUIRE(liveModels == 0);
        buffer.Emplace<FragileModel>(Fragile{2, false});
        REQUIRE(liveModels == 1);
        REQUIRE(ValueOf(buffer) == 2);
    }
    REQUIRE(liveModels == 0);
}

TEST_CASE("the model is borrowed from an lvalue buffer only", "[core][erased]")
{
    // ADR-0046, section 3 (AC-0046-02): Get() on a temporary buffer would
    // return a pointer into storage that is about to go.
    STATIC_REQUIRE(GetFromLvalue<Buffer>);
    STATIC_REQUIRE_FALSE(GetFromRvalue<Buffer>);
}

TEST_CASE("a vector of buffers moves its models with it", "[core][erased]")
{
    // Reallocation moves every model into new storage and destroys the old
    // ones; the pointers into the old buffers are gone with them.
    liveModels = 0;
    {
        std::vector<Buffer> buffers;
        for (int i = 0; i < 16; ++i) {
            Buffer b;
            if (i % 2 == 0) {
                b.Emplace<SmallModel>(Small{i});
            } else {
                b.Emplace<LargeModel>(Large{i, {}});
            }
            buffers.push_back(std::move(b));
        }
        REQUIRE(liveModels == 16);
        for (int i = 0; i < 16; ++i) {
            REQUIRE(ValueOf(buffers[static_cast<std::size_t>(i)]) == i);
        }
    }
    REQUIRE(liveModels == 0);
}
