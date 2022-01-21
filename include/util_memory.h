//
// Distributed Linear Algebra Interface (DLAI)
//
// Copyright (c) 2018-2021, ETH Zurich
// All rights reserved.
//
// Please, refer to the LICENSE file in the root directory.
// SPDX-License-Identifier: BSD-3-Clause
//

#ifndef DLA_INTERFACE_UTIL_MEMORY_H
#define DLA_INTERFACE_UTIL_MEMORY_H

#include <cstring>
#include <type_traits>
#include <utility>
#include "types.h"
#include "util_math.h"

namespace dla_interface {
  namespace util {
    namespace memory {

      /// Copies src_ptr[i] to dest_ptr[i] for 0 <= i < n.<br>
      ///
      /// **Precondition:**
      /// - src_ptr[i] and dest_ptr[i] are valid for 0 <= i < n.
      /// - ElType is trivially copyable.
      template <class ElType>
      void copy(std::size_t n, const ElType* src_ptr, ElType* dest_ptr) {
        static_assert(std::is_trivially_copyable<ElType>(), "Type is not trivially copyable!");
        std::memcpy(dest_ptr, src_ptr, sizeof(ElType) * n);
      }

      /// Copies src_ptr[k1] to dest_ptr[k2]
      /// - for k = i + j * ld_src and k = i + j * ld_dest,
      /// - for 0 <= i < size.first and 0 <= j < size.second.
      ///
      /// **Precondition:**
      /// - src_ptr[k1] is valid for 0 <= k1 < ld_src * size.second,
      /// - dest_ptr[k2] is valid for 0 <= k1 < ld_src * size.second,
      /// - ElType is trivially copyable.
      template <class ElType>
      void copy2D(const std::pair<SizeType, SizeType>& size, const ElType* src_ptr, int ld_src,
                  ElType* dest_ptr, int ld_dest) {
        static_assert(std::is_trivially_copyable<ElType>(), "Type is not trivially copyable!");
        std::size_t len = multSize(size.first, size.second);
        if (len > 0) {
          if (size.first == ld_src && size.first == ld_dest) {
            copy(len, src_ptr, dest_ptr);
          }
          else {
            for (int j = 0; j < size.second; ++j) {
              copy(size.first, src_ptr + multSize(ld_src, j), dest_ptr + multSize(ld_dest, j));
            }
          }
        }
      }
    }
  }
}

#endif  // DLA_INTERFACE_UTIL_MEMORY_H
