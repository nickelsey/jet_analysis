// string_util.hh

#include <string>
#include <sstream>

inline void MakeStringInternal(std::stringstream&) {}

template <typename T>
inline void MakeStringInternal(std::stringstream& ss, const T& t) {
  ss << t;
}

template <typename T, typename... Args>
inline void
MakeStringInternal(std::stringstream& ss, const T& t, const Args&... args) {
  MakeStringInternal(ss, t);
  MakeStringInternal(ss, args...);
}

// given a set of inputs, use a stringstream to consruct a string
template <typename... Args>
std::string MakeString(const Args&... args) {
  std::stringstream ss;
  MakeStringInternal(ss, args...);
  return std::string(ss.str());
}

template <>
inline std::string MakeString(const std::string& str) {
  return str;
}

inline std::string MakeString(const char* c_str) {
  return std::string(c_str);
}


// create a string by contatenating the entries of an standard library
// style container, with a specific delimiter
template <class Container>
inline std::string Join(const std::string& delim, const Container& container) {
  std::stringstream ss;
  auto count = container.size() - 1;
  for (auto i = container.begin(); i != container.end(); ++i, --count) {
    ss << *i << (count ? delim : "");
  }
  return ss.str();
}
