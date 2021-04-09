//
// Distributed Linear Algebra Interface (DLAI)
//
// Copyright (c) 2018-2019, ETH Zurich
// All rights reserved.
//
// Please, refer to the LICENSE file in the root directory.
// SPDX-License-Identifier: BSD-3-Clause
//

#ifndef DLA_INTERFACE_MATRIX_INDEX_H
#define DLA_INTERFACE_MATRIX_INDEX_H

#include "types.h"

namespace dla_interface {

  struct Global2DIndex {
    Global2DIndex() : row(0), col(0) {}
    Global2DIndex(IndexType row_index, IndexType col_index) : row(row_index), col(col_index) {}

    bool operator==(const Global2DIndex& rhs) const {
      return row == rhs.row && col == rhs.col;
    }
    bool operator!=(const Global2DIndex& rhs) const {
      return row != rhs.row || col != rhs.col;
    }

    IndexType row;
    IndexType col;
  };

  inline Global2DIndex operator+(const Global2DIndex& a, const Global2DIndex& b) {
    return Global2DIndex(a.row + b.row, a.col + b.col);
  }
  inline Global2DIndex add(const Global2DIndex& a, IndexType row, IndexType col) {
    return Global2DIndex(a.row + row, a.col + col);
  }
  inline Global2DIndex operator-(const Global2DIndex& a, const Global2DIndex& b) {
    return Global2DIndex(a.row - b.row, a.col - b.col);
  }

  struct Local2DIndex {
    Local2DIndex() : row(0), col(0) {}
    Local2DIndex(IndexType row_index, IndexType col_index) : row(row_index), col(col_index) {}

    bool operator==(const Local2DIndex& rhs) const {
      return row == rhs.row && col == rhs.col;
    }
    bool operator!=(const Local2DIndex& rhs) const {
      return row != rhs.row || col != rhs.col;
    }

    IndexType row;
    IndexType col;
  };

  inline Local2DIndex operator+(const Local2DIndex& a, const Local2DIndex& b) {
    return Local2DIndex(a.row + b.row, a.col + b.col);
  }
  inline Local2DIndex operator-(const Local2DIndex& a, const Local2DIndex& b) {
    return Local2DIndex(a.row - b.row, a.col - b.col);
  }
}

#endif  // DLA_INTERFACE_MATRIX_INDEX_H
