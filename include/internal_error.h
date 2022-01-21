//
// Distributed Linear Algebra Interface (DLAI)
//
// Copyright (c) 2018-2021, ETH Zurich
// All rights reserved.
//
// Please, refer to the LICENSE file in the root directory.
// SPDX-License-Identifier: BSD-3-Clause
//

#ifndef DLA_INTERFACE_INTERNAL_ERROR_H
#define DLA_INTERFACE_INTERNAL_ERROR_H

#include <stdexcept>

namespace dla_interface {
  namespace error {
    class InternalError : public std::runtime_error {
      public:
      explicit InternalError(const std::string& what_arg) : std::runtime_error(what_arg) {}
      explicit InternalError(const char* what_arg) : std::runtime_error(what_arg) {}
    };
  }
}

#endif  // DLA_INTERFACE_INTERNAL_ERROR_H
