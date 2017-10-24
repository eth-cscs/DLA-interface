#ifndef DLA_INTERFACE_TEST_INCLUDE_UTIL_COMPLEX_H
#define DLA_INTERFACE_TEST_INCLUDE_UTIL_COMPLEX_H

#include <complex>
#include <type_traits>

namespace testing {

  // Returns r for floating point type,
  // Returns r * exp(2*Pi*i*arg/12)
  template <class ElType>
  std::enable_if_t<std::is_floating_point<ElType>::value, ElType> value(ElType r, int arg) {
    return r;
  }

  template <class ElType>
  std::enable_if_t<std::is_floating_point<typename ElType::value_type>::value, ElType> value(
      typename ElType::value_type r, int arg) {
    // clang-format off
    constexpr typename ElType::value_type vals[12] = {
      0, 0.5, 0.86602540378443864676372317075293618, 1, 0.86602540378443864676372317075293618, 0.5,
      0, -0.5, -0.86602540378443864676372317075293618, -1, -0.86602540378443864676372317075293618, -0.5
    };
    // clang-format on
    int index_sin = (arg % 12 + 12) % 12;
    int index_cos = (index_sin + 3) % 12;

    return r * ElType(vals[index_cos], vals[index_sin]);
  }
}

#endif  // DLA_INTERFACE_TEST_INCLUDE_UTIL_COMPLEX_H
