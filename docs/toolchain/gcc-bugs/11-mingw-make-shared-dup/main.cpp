#include <memory>
int fa();
int main() { return *std::make_shared<int>(1) + fa() - 2; }
