// arg_helper.hh

// helper functions for parsing command line
// options in c++ code

#ifndef ARG_HELPER_HH
#define ARG_HELPER_HH

#include <string>
#include <sstream>
#include <set>


bool Consume(std::string& arg, const std::string& seq) {
  if (seq.size() > arg.size()) return false;
  if (arg.compare( 0, seq.length(), seq ) == 0) {
    arg = std::string(arg.begin() + seq.length(), arg.end());
    return true;
  }
  return false;
}

bool Consume(std::string& arg, char c) {
  if (arg.size() == 0)
    return false;
  if (arg[0] == c) {
    arg = std::string(arg.begin()+1, arg.end());
    return true;
  }
  return false;
}

template<typename T>
bool CanCast(std::string s) {
  std::istringstream iss(s);
  T dummy;
  iss >> std::skipws >> dummy;
  return iss && iss.eof();
}

template<typename T>
T CastTo(std::string s) {
  std::istringstream iss(s);
  T dummy;
  iss >> std::skipws >> dummy;
  return dummy;
}

bool ParseBoolFlag(std::string arg, const std::string &flag, bool* opt) {
  if (Consume(arg, flag)) {
    if (arg.empty() || arg == "=true") {
      *opt = true;
      return true;
    }
    else if (arg == "=false") {
      *opt = false;
      return true;
    }
  }
  return false;
}

bool ParseStrFlag(std::string arg, const std::string &flag, std::string* opt) {
  if (Consume(arg, flag) && Consume(arg, "=")) {
    *opt = arg;
    return true;
  }
  return false;
}

bool ParseFloatFlag(std::string arg, const std::string &flag, double* opt) {
  if (Consume(arg, flag) && Consume(arg, "=")) {
    if (CanCast<double>(arg)){
      *opt = atof(arg.c_str());
      return true;
    }
  }
  return false;
}

bool ParseIntFlag(std::string arg, const std::string &flag, int* opt) {
  if (Consume(arg, flag) && Consume(arg, "=")) {
    if (CanCast<int>(arg)){
      *opt = atof(arg.c_str());
      return true;
    }
  }
  return false;
}

bool HasEnding(const std::string &full_string, const std::string &ending) {
  if (full_string.length() >= ending.length()) {
    return (0 == full_string.compare (full_string.length() - ending.length(), ending.length(), ending) );
  } else {
    return false;
  }
}

namespace patch {
  template<typename T, typename P>
  T remove_if(T beg, T end, P pred)
  {
    T dest = beg;
    for (T itr = beg;itr != end; ++itr)
      if (!pred(*itr))
        *(dest++) = *itr;
    return dest;
  }
}

template<typename T>
std::set<T> ParseArgString(std::string str) {
  std::set<T> ret;
  
  // remove all spaces
  str.erase(patch::remove_if(str.begin(), str.end(), ::isspace), str.end());
  
  std::string token;
  while ( str.find(",") != std::string::npos ) {
    size_t pos = str.find(",");
    token = str.substr(0, pos);
    if (CanCast<T>(token)) {
      ret.insert(CastTo<T>(token));
    }
    str.erase(0, pos + 1);
  }
  if (CanCast<T>(str))
    ret.insert(CastTo<T>(str));
  
  return ret;
}

template<typename T>
std::vector<T> ParseArgStringToVec(std::string str) {
  std::vector<T> ret;
  
  // remove all spaces
  str.erase(patch::remove_if(str.begin(), str.end(), ::isspace), str.end());
  
  std::string token;
  while ( str.find(",") != std::string::npos ) {
    size_t pos = str.find(",");
    token = str.substr(0, pos);
    if (CanCast<T>(token)) {
      ret.push_back(CastTo<T>(token));
    }
    str.erase(0, pos + 1);
  }
  if (CanCast<T>(str))
    ret.push_back(CastTo<T>(str));
  
  return ret;
}

#endif // ARG_HELPER_HH
