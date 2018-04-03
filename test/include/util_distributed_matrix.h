#ifndef DLA_INTERFACE_TEST_INCLUDE_UTIL_DISTRIBUTED_MATRIX_H
#define DLA_INTERFACE_TEST_INCLUDE_UTIL_DISTRIBUTED_MATRIX_H

#include <cmath>
#include <cassert>
#include <iostream>
#include "types.h"

namespace testing {
  using namespace dla_interface;
  template <class Type, class F>
  void fillDistributedMatrix(DistributedMatrix<Type>& mat, F& el_val) {
    for (int j = 0; j < mat.size().second; ++j) {
      for (int i = 0; i < mat.size().first; ++i) {
        if (mat.commGrid().id2D() == mat.getRankId2D(i, j))
          mat(i, j) = el_val(i, j);
      }
    }
  }

  template <class Type, class F, class Out>
  bool checkDistributedMatrix(const DistributedMatrix<Type>& mat, F& el_val, Out& out) {
    for (int j = 0; j < mat.size().second; ++j) {
      for (int i = 0; i < mat.size().first; ++i) {
        if (mat.commGrid().id2D() == mat.getRankId2D(i, j)) {
          if (mat(i, j) != el_val(i, j)) {
            out << "Element (" << i << ", " << j << ") is wrong." << std::endl;
            out << "Expected " << el_val(i, j) << "," << std::endl;
            out << "got " << mat(i, j) << "." << std::endl;
            return false;
          }
        }
      }
    }
    return true;
  }

  template <class Type, class F>
  bool checkDistributedMatrix(const DistributedMatrix<Type>& mat, F& el_val) {
    return checkDistributedMatrix(mat, el_val, std::cout);
  }

  template <class Type, class F, class Out>
  bool checkNearDistributedMatrix(const DistributedMatrix<Type>& mat, F& el_val,
                                  BaseType<Type> diff, BaseType<Type> threshold_rel_abs, Out& out) {
    assert(diff > 0);
    assert(threshold_rel_abs > 0);

    for (int j = 0; j < mat.size().second; ++j) {
      for (int i = 0; i < mat.size().first; ++i) {
        if (mat.commGrid().id2D() == mat.getRankId2D(i, j)) {
          BaseType<Type> el_diff =
              std::abs(mat(i, j) - el_val(i, j)) /
              std::max(std::max(threshold_rel_abs, std::abs(mat(i, j))), std::abs(el_val(i, j)));
          if (el_diff > diff) {
            out << "Element (" << i << ", " << j << ") is wrong." << std::endl;
            out << "Expected " << el_val(i, j) << "," << std::endl;
            out << "got " << mat(i, j) << "." << std::endl;
            out << "Difference abs: " << std::abs(mat(i, j) - el_val(i, j)) << "." << std::endl;
            out << "Difference rel/abs: " << el_diff << "." << std::endl;
            return false;
          }
        }
      }
    }
    return true;
  }

  template <class Type, class F>
  bool checkNearDistributedMatrix(const DistributedMatrix<Type>& mat, F& el_val,
                                  BaseType<Type> diff, BaseType<Type> threshold_rel_abs = 1e-3) {
    return checkNearDistributedMatrix(mat, el_val, diff, threshold_rel_abs, std::cout);
  }

  template <class Type, class Out>
  bool checkDistributedMatrixSame(const DistributedMatrix<Type>& mat1,
                                  const DistributedMatrix<Type>& mat2, Out& out) {
    if (mat1.size() != mat2.size()) {
      out << "Sizes are different";
      out << " (" << mat1.size().first << ", " << mat1.size().second << ")";
      out << " != (" << mat2.size().first << ", " << mat2.size().second << ")" << std::endl;
      return false;
    }

    for (int j = 0; j < mat1.size().second; ++j) {
      for (int i = 0; i < mat1.size().first; ++i) {
        if (mat1.getRankId2D(i, j) != mat2.getRankId2D(i, j)) {
          out << "Different element location (" << i << ", " << j << ")" << std::endl;
          return false;
        }
        if (mat1.commGrid().id2D() == mat2.getRankId2D(i, j)) {
          if (mat1(i, j) != mat2(i, j)) {
            out << "Element (" << i << ", " << j << ") is wrong." << std::endl;
            out << "Expected " << mat1(i, j) << "," << std::endl;
            out << "got " << mat2(i, j) << "." << std::endl;
            return false;
          }
        }
      }
    }
    return true;
  }

  template <class Type>
  bool checkDistributedMatrixSame(const DistributedMatrix<Type>& mat1,
                                  const DistributedMatrix<Type>& mat2) {
    return checkDistributedMatrixSame(mat1, mat2, std::cout);
  }

  template <class Type, class Out>
  bool checkIfDistributedSubMatrix(const DistributedMatrix<Type>& mat,
                                   std::pair<SizeType, SizeType> subsize,
                                   std::pair<IndexType, IndexType> subindex,
                                   const DistributedMatrix<Type>& submat, Out& out) {
    if (subsize != submat.size()) {
      out << "Sizes are different";
      out << " (" << subsize.first << ", " << subsize.second << ")";
      out << " != (" << submat.size().first << ", " << submat.size().second << ")" << std::endl;
      return false;
    }

    for (int j = 0; j < subsize.second; ++j) {
      for (int i = 0; i < subsize.first; ++i) {
        if (mat.getRankId2D(subindex.first + i, subindex.second + j) != submat.getRankId2D(i, j)) {
          out << "Different element location (" << i << ", " << j << ")" << std::endl;
          return false;
        }
        if (mat.commGrid().id2D() == submat.getRankId2D(i, j)) {
          if (mat.ptr(subindex.first + i, subindex.second + j) != submat.ptr(i, j)) {
            out << "Pointers of element (" << i << ", " << j << ") are not the same."
                      << std::endl;
            return false;
          }
        }
      }
    }
    return true;
  }

  template <class Type>
  bool checkIfDistributedSubMatrix(const DistributedMatrix<Type>& mat,
                                   std::pair<SizeType, SizeType> subsize,
                                   std::pair<IndexType, IndexType> subindex,
                                   const DistributedMatrix<Type>& submat) {
    return checkIfDistributedSubMatrix(mat, subsize, subindex, submat, std::cout);
  }
}

#endif  // DLA_INTERFACE_TEST_INCLUDE_UTIL_DISTRIBUTED_MATRIX_H
