import b;
void g(V1& r) { r = 2 * r; }
int main() { V1 v{1}; g(v); return v.x == 2 ? 0 : 1; }
