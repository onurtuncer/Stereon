#include <format>
#include <string_view>
unsigned fa();
int main() { return std::format("{}", std::string_view("y")).size() + fa() - 2; }
