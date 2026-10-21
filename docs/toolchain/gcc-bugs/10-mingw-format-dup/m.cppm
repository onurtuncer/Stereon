module;
#include <format>
#include <string_view>
export module m;
export template <class T> unsigned fmt(const T* s) { return std::format("{}", std::string_view(s)).size(); }
