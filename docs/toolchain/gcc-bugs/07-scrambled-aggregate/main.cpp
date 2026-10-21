#include <cstdio>
import m;
int main() {
  long r = S::f(Empty{}, Big{1, 2, 3});
  std::printf("%ld\n", r);
  return r == 123 ? 0 : 1;
}
