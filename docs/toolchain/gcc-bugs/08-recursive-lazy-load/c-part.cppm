export module c:part;
import b;
export template <typename S>
V1 g(S) {
  A<V1> a{};
  return a.e[2];
}
