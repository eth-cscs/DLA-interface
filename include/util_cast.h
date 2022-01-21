//
// Distributed Linear Algebra Interface (DLAI)
//
// Copyright (c) 2018-2021, ETH Zurich
// All rights reserved.
//
// Please, refer to the LICENSE file in the root directory.
// SPDX-License-Identifier: BSD-3-Clause
//

#ifndef DLA_INTERFACE_UTIL_CAST_H
#define DLA_INTERFACE_UTIL_CAST_H

#include <complex>

namespace dla_interface {
  namespace util {

    template <class Type>
    inline Type& castToC(Type& x) {
      return x;
    }

    template <class Type>
    inline const Type& castToC(const Type& x) {
      return x;
    }

    inline _Complex float castToC(std::complex<float>& x) {
      return *reinterpret_cast<float(&)[2]>(x);
    }

    inline const _Complex float castToC(const std::complex<float>& x) {
      return *reinterpret_cast<const float(&)[2]>(x);
    }

    inline _Complex double castToC(std::complex<double>& x) {
      return *reinterpret_cast<double(&)[2]>(x);
    }

    inline const _Complex double castToC(const std::complex<double>& x) {
      return *reinterpret_cast<const double(&)[2]>(x);
    }
  }
}

#endif  // DLA_INTERFACE_UTIL_CAST_H
