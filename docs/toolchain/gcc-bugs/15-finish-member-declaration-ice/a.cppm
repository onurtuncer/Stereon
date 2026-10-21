export module a;
export template <typename T> struct V {
  T x;
  friend V operator*(T s, V v) { return {s * v.x}; }
};
export using V1 = V<int>;
