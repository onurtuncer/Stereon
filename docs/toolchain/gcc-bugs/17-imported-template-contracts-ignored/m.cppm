export module m;
export template <class T> T g(T x) pre(x > 0) { return x; }
