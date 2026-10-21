module;
#include <source_location>
export module m;
export struct Record {
  std::source_location loc;
};
export inline Record make(std::source_location loc = std::source_location::current()) {
  return {loc};
}
