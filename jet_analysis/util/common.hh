#ifndef JET_ANALYSIS_UTILS_COMMON_HH
#define JET_ANALYSIS_UTILS_COMMON_HH

#include <string>
#include <vector>
#include <memory>
#include <map>
#include <unordered_map>
#include <set>

#include "gflags/gflags.h"
#include "glog/stl_logging.h"
#include "glog/logging.h"

// implementation of the make_unique in case we are compiling with c++11
namespace std {
  template<class T> struct _Unique_if {
    typedef unique_ptr<T> _Single_object;
  };
  
  template<class T> struct _Unique_if<T[]> {
    typedef unique_ptr<T[]> _Unknown_bound;
  };
  
  template<class T, size_t N> struct _Unique_if<T[N]> {
    typedef void _Known_bound;
  };
  
  template<class T, class... Args>
  typename _Unique_if<T>::_Single_object
  make_unique(Args&&... args) {
    return unique_ptr<T>(new T(std::forward<Args>(args)...));
  }
  
  template<class T>
  typename _Unique_if<T>::_Unknown_bound
  make_unique(size_t n) {
    typedef typename remove_extent<T>::type U;
    return unique_ptr<T>(new U[n]());
  }
  
  template<class T, class... Args>
  typename _Unique_if<T>::_Known_bound
  make_unique(Args&&...) = delete;
}

// commonly used classes from std, use statements
// in scaffold namespace so that including scaffold headers
// won't modify the global namespace
using std::pair;
using std::string;
using std::vector;
using std::unique_ptr;
using std::shared_ptr;
using std::make_unique;
using std::map;
using std::unordered_map;
using std::set;

// if using enum classes and pairs as keys, need to specify
// the hash function
struct EnumClassHash {
  template <typename T>
  std::size_t operator()(T t) const {
    return static_cast<std::size_t>(t);
  }
};

struct PairHash {
  template <typename T, typename U>
  std::size_t operator()(const std::pair<T, U> &x) const
  {
    return std::hash<T>()(x.first) ^ std::hash<U>()(x.second);
  }
};

struct PairEnumHash {
  template <typename T, typename U>
  std::size_t operator()(const std::pair<T, U> &x) const
  {
    return EnumClassHash()(x.first) ^ EnumClassHash()(x.second);
  }
};

#endif // JET_ANALYSIS_UTILS_COMMON_HH
