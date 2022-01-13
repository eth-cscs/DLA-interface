//
// Distributed Linear Algebra Interface (DLAI)
//
// Copyright (c) 2018-2021, ETH Zurich
// All rights reserved.
//
// Please, refer to the LICENSE file in the root directory.
// SPDX-License-Identifier: BSD-3-Clause
//

#ifndef DLA_INTERFACE_UTIL_DEBUG_OUTPUT_H
#define DLA_INTERFACE_UTIL_DEBUG_OUTPUT_H

#ifdef DLAI_PRINT_DEBUG_CALL_PARAM
#include <iostream>
#include "util_output.h"
#include "util_types.h"
#endif

namespace dla_interface {
  namespace util {

#ifdef DLAI_PRINT_DEBUG_CALL_PARAM
    template <class T>
    void debug_print(std::stringstream& out, T arg1) {
      out << arg1;
    }

    template <class T, class... Args>
    void debug_print(std::stringstream& out, T arg1, Args... args) {
      out << arg1 << ", ";
      debug_print(out, args...);
    }

    template <class... Args>
    void debug_print_call_param(const char* func, Args... args) {
      std::stringstream out;
      out << func << "(";
      debug_print(out, args...);
      out << ")" << std::endl;
      std::cout << out.str();
    }

#define dlai__util__debug_print_call_param(...) \
  dla_interface::util::debug_print_call_param(__func__, __VA_ARGS__)

#else
#define dlai__util__debug_print_call_param(...)
#endif
  }
}

#endif  // DLA_INTERFACE_UTIL_DEBUG_OUTPUT_H
