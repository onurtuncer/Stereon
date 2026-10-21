module;
#include <source_location>
export module m;
export template <class T> int f() { return std::source_location::current().line(); }
