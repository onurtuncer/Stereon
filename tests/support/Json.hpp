// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// A small JSON reader for the corpus cases (corpus/README.md): the values a
// case file holds, parsed by recursive descent. Malformed input throws
// std::runtime_error naming the offset, which fails the test that reads it.
// Tests only; the kernel reads no JSON.

#pragma once

#include <charconv>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

namespace StereonTest {

class JsonValue {
public:
    using Array = std::vector<JsonValue>;
    using Object = std::vector<std::pair<std::string, JsonValue>>;

    std::variant<std::nullptr_t, bool, double, std::string, Array, Object> Value;

    [[nodiscard]] bool IsNull() const { return std::holds_alternative<std::nullptr_t>(Value); }
    [[nodiscard]] bool IsObject() const { return std::holds_alternative<Object>(Value); }
    [[nodiscard]] bool AsBool() const { return std::get<bool>(Value); }
    [[nodiscard]] double AsNumber() const { return std::get<double>(Value); }
    [[nodiscard]] const std::string& AsString() const { return std::get<std::string>(Value); }
    [[nodiscard]] const Array& AsArray() const { return std::get<Array>(Value); }
    [[nodiscard]] const Object& AsObject() const { return std::get<Object>(Value); }

    /// The member `key`, or null when there is none.
    [[nodiscard]] const JsonValue* Find(const std::string_view key) const
    {
        for (const auto& [name, value] : AsObject()) {
            if (name == key) {
                return &value;
            }
        }
        return nullptr;
    }

    [[nodiscard]] bool Has(const std::string_view key) const { return Find(key) != nullptr; }

    [[nodiscard]] const JsonValue& operator[](const std::string_view key) const
    {
        if (const JsonValue* value = Find(key)) {
            return *value;
        }
        throw std::runtime_error("JSON object has no member \"" + std::string(key) + "\"");
    }

    [[nodiscard]] const JsonValue& operator[](const std::size_t index) const
    {
        return AsArray().at(index);
    }
};

namespace Detail {

class JsonParser {
public:
    explicit JsonParser(const std::string_view text)
        : m_Text(text)
    {}

    JsonValue ParseDocument()
    {
        JsonValue value = ParseValue();
        SkipSpace();
        if (m_Pos != m_Text.size()) {
            Fail("trailing characters");
        }
        return value;
    }

private:
    [[noreturn]] void Fail(const std::string& what) const
    {
        throw std::runtime_error("JSON: " + what + " at offset " + std::to_string(m_Pos));
    }

    void SkipSpace()
    {
        while (m_Pos < m_Text.size() && (m_Text[m_Pos] == ' ' || m_Text[m_Pos] == '\n' ||
                                         m_Text[m_Pos] == '\r' || m_Text[m_Pos] == '\t')) {
            ++m_Pos;
        }
    }

    [[nodiscard]] char Peek() const { return m_Pos < m_Text.size() ? m_Text[m_Pos] : '\0'; }

    void Expect(const char c)
    {
        if (Peek() != c) {
            Fail(std::string("expected '") + c + "'");
        }
        ++m_Pos;
    }

    bool Consume(const std::string_view word)
    {
        if (m_Text.substr(m_Pos, word.size()) == word) {
            m_Pos += word.size();
            return true;
        }
        return false;
    }

    JsonValue ParseValue()
    {
        SkipSpace();
        switch (Peek()) {
        case '{':
            return ParseObject();
        case '[':
            return ParseArray();
        case '"':
            return JsonValue{ParseString()};
        case 't':
            return Consume("true") ? JsonValue{true} : (Fail("bad literal"), JsonValue{});
        case 'f':
            return Consume("false") ? JsonValue{false} : (Fail("bad literal"), JsonValue{});
        case 'n':
            return Consume("null") ? JsonValue{nullptr} : (Fail("bad literal"), JsonValue{});
        default:
            return ParseNumber();
        }
    }

    JsonValue ParseNumber()
    {
        double value = 0.0;
        const char* begin = m_Text.data() + m_Pos;
        const char* end = m_Text.data() + m_Text.size();
        const std::from_chars_result r = std::from_chars(begin, end, value);
        if (r.ec != std::errc{} || r.ptr == begin) {
            Fail("expected a value");
        }
        m_Pos = static_cast<std::size_t>(r.ptr - m_Text.data());
        return JsonValue{value};
    }

    std::string ParseString()
    {
        Expect('"');
        std::string s;
        while (m_Pos < m_Text.size() && m_Text[m_Pos] != '"') {
            const char c = m_Text[m_Pos++];
            if (c != '\\') {
                s += c;
                continue;
            }
            s += Unescaped();
        }
        Expect('"');
        return s;
    }

    std::string Unescaped()
    {
        if (m_Pos >= m_Text.size()) {
            Fail("unterminated escape");
        }
        const char e = m_Text[m_Pos++];
        switch (e) {
        case 'b':
            return "\b";
        case 'f':
            return "\f";
        case 'n':
            return "\n";
        case 'r':
            return "\r";
        case 't':
            return "\t";
        case 'u': {
            unsigned code = 0;
            const char* begin = m_Text.data() + m_Pos;
            if (m_Pos + 4 > m_Text.size() ||
                std::from_chars(begin, begin + 4, code, 16).ptr != begin + 4) {
                Fail("bad \\u escape");
            }
            m_Pos += 4;
            return code < 0x80 ? std::string(1, static_cast<char>(code)) : std::string("?");
        }
        default:
            return std::string(1, e); // \" \\ \/
        }
    }

    JsonValue ParseArray()
    {
        Expect('[');
        JsonValue::Array items;
        SkipSpace();
        if (Peek() == ']') {
            ++m_Pos;
            return JsonValue{std::move(items)};
        }
        for (;;) {
            items.push_back(ParseValue());
            SkipSpace();
            if (Peek() == ',') {
                ++m_Pos;
                continue;
            }
            Expect(']');
            return JsonValue{std::move(items)};
        }
    }

    JsonValue ParseObject()
    {
        Expect('{');
        JsonValue::Object members;
        SkipSpace();
        if (Peek() == '}') {
            ++m_Pos;
            return JsonValue{std::move(members)};
        }
        for (;;) {
            SkipSpace();
            std::string key = ParseString();
            SkipSpace();
            Expect(':');
            members.emplace_back(std::move(key), ParseValue());
            SkipSpace();
            if (Peek() == ',') {
                ++m_Pos;
                continue;
            }
            Expect('}');
            return JsonValue{std::move(members)};
        }
    }

    std::string_view m_Text;
    std::size_t m_Pos = 0;
};

} // namespace Detail

[[nodiscard]] inline JsonValue ParseJson(const std::string_view text)
{
    return Detail::JsonParser(text).ParseDocument();
}

[[nodiscard]] inline JsonValue ReadJsonFile(const std::filesystem::path& path)
{
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::runtime_error("cannot open " + path.string());
    }
    std::ostringstream text;
    text << in.rdbuf();
    return ParseJson(text.str());
}

} // namespace StereonTest
