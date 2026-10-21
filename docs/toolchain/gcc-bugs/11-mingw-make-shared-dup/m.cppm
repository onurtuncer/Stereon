module;
#include <memory>
export module m;
export template <class T> std::shared_ptr<T> make(T v) { return std::make_shared<T>(v); }
