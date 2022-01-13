//
// Distributed Linear Algebra Interface (DLAI)
//
// Copyright (c) 2018-2021, ETH Zurich
// All rights reserved.
//
// Please, refer to the LICENSE file in the root directory.
// SPDX-License-Identifier: BSD-3-Clause
//

#ifndef DLA_INTERFACE_UTIL_OUTPUT_H
#define DLA_INTERFACE_UTIL_OUTPUT_H

#include <ostream>
#include <iomanip>
#include "matrix_index.h"

namespace dla_interface {
  namespace util {

    template <class t1, class t2>
    std::ostream& operator<<(std::ostream& s, const std::pair<t1, t2>& value) {
      return s << "(" << value.first << ", " << value.second << ")";
    }

    inline std::ostream& operator<<(std::ostream& s, const Local2DIndex& index) {
      return s << "(" << index.row << ", " << index.col << ")";
    }
    inline std::ostream& operator<<(std::ostream& s, const Global2DIndex& index) {
      return s << "(" << index.row << ", " << index.col << ")";
    }

    template <class Out, class ElType>
    void dumpDistributedMatrixElement(Out& out, const Global2DIndex& gindex, Local2DIndex& lindex,
                                      std::size_t storage_index, std::size_t storage_base_index,
                                      ElType val) {
      out << std::setw(16) << gindex << std::setw(16) << lindex << std::setw(12) << storage_index
          << " (" << std::setw(12) << storage_base_index << ") " << val << std::endl;
    }
  }
}

using dla_interface::util::operator<<;

#endif  // DLA_INTERFACE_UTIL_OUTPUT_H
