export module a;
export template <typename T> struct V { T x; };
export template <typename T> struct P {
  T x;
  friend V<T> operator-(const P& p, const P& q) { return {p.x - q.x}; }
};
export inline int g(const P<int>& p) { return p.x; }
