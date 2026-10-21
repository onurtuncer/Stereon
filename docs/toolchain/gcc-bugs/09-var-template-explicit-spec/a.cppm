export module a;
export struct L {};
export template <typename T> inline constexpr int N = 0;
template <> inline constexpr int N<L> = 1;
