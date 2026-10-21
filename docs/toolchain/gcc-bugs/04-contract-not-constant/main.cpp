struct S {
  int n = 1;
  constexpr int size() const { return n; }
  constexpr int get() const pre(size() > 0) { return n; }
};

constexpr int f(const S &s) { return s.size() + s.get(); }

static_assert(f(S{}) == 2);

int main() { return f(S{}) - 2; }
