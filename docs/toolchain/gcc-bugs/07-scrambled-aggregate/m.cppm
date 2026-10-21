export module m;
export struct Empty {};
export struct Big { long a, b, c; };
export struct S {
  static long f(Empty, Big b) { return b.a * 100 + b.b * 10 + b.c; }
};
