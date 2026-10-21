export module m;
export struct S {
  int x;
  friend bool operator==(const S&, const S&) = default;
};
