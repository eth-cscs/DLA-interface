#ifndef DLA_INTERFACE_TEST_INCLUDE_UTIL_DISTRIBUTED_MATRIX_H
#define DLA_INTERFACE_TEST_INCLUDE_UTIL_DISTRIBUTED_MATRIX_H

#include <cmath>
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

  template <class Type, class F>
  bool checkDistributedMatrix(const DistributedMatrix<Type>& mat, F& el_val) {
    for (int j = 0; j < mat.size().second; ++j) {
      for (int i = 0; i < mat.size().first; ++i) {
        if (mat.commGrid().id2D() == mat.getRankId2D(i, j)) {
          if (mat(i, j) != el_val(i, j)) {
            std::cout << "Element (" << i << ", " << j << ") is wrong." << std::endl;
            std::cout << "Expected " << el_val(i, j) << "," << std::endl;
            std::cout << "Got " << mat(i, j) << "." << std::endl;
            return false;
          }
        }
      }
    }
    return true;
  }

  template <class Type, class F>
  bool checkNearDistributedMatrix(const DistributedMatrix<Type>& mat, F& el_val, BaseType<Type> diff) {
    for (int j = 0; j < mat.size().second; ++j) {
      for (int i = 0; i < mat.size().first; ++i) {
        if (mat.commGrid().id2D() == mat.getRankId2D(i, j)) {
          if (std::abs(mat(i, j) - el_val(i, j)) > diff) {
            std::cout << "Element (" << i << ", " << j << ") is wrong." << std::endl;
            std::cout << "Expected " << el_val(i, j) << "," << std::endl;
            std::cout << "Got " << mat(i, j) << "." << std::endl;
            std::cout << "Difference " << std::abs(mat(i, j) - el_val(i, j)) << "." << std::endl;
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
    if (mat1.size() != mat2.size()) {
      std::cout << "Sizes are different";
      std::cout << " (" << mat1.size().first << ", " << mat1.size().second << ")";
      std::cout << " != (" << mat2.size().first << ", " << mat2.size().second << ")" << std::endl;
      return false;
    }
    if (mat1.localSize() != mat2.localSize()) {
      std::cout << "Local sizes are different";
      std::cout << " (" << mat1.localSize().first << ", " << mat1.localSize().second << ")";
      std::cout << " != (" << mat2.localSize().first << ", " << mat2.localSize().second << ")"
                << std::endl;
      return false;
    }
    if (&mat1.commGrid() != &mat2.commGrid()) {
      std::cout << "Communicator grids are different" << std::endl;
      return false;
    }

    for (int j = 0; j < mat1.size().second; ++j) {
      for (int i = 0; i < mat1.size().first; ++i) {
        if (mat1.getRankId2D(i, j) != mat2.getRankId2D(i, j)) {
          std::cout << "Different element location (" << i << ", " << j << ")" << std::endl;
          return false;
        }
        if (mat1.commGrid().id2D() == mat2.getRankId2D(i, j)) {
          if (mat1(i, j) != mat2(i, j)) {
            std::cout << "Element (" << i << ", " << j << ") is wrong." << std::endl;
            std::cout << "Expected " << mat1(i, j) << "," << std::endl;
            std::cout << "Got " << mat2(i, j) << "." << std::endl;
            return false;
          }
        }
      }
    }
    return true;
  }

  template <class Type>
  bool checkDistributedSubMatrixSame(const DistributedMatrix<Type>& mat,
                                     std::pair<SizeType, SizeType> subsize,
                                     std::pair<IndexType, IndexType> subindex,
                                     const DistributedMatrix<Type>& submat) {
    if (subsize != submat.size()) {
      std::cout << "Sizes are different";
      std::cout << " (" << subsize.first << ", " << subsize.second << ")";
      std::cout << " != (" << submat.size().first << ", " << submat.size().second << ")" << std::endl;
      return false;
    }

    if (&mat.commGrid() != &submat.commGrid()) {
      std::cout << "Communicator grids are different" << std::endl;
      return false;
    }

    for (int j = 0; j < subsize.second; ++j) {
      for (int i = 0; i < subsize.first; ++i) {
        if (mat.getRankId2D(subindex.first + i, subindex.second + j) != submat.getRankId2D(i, j)) {
          std::cout << "Different element location (" << i << ", " << j << ")" << std::endl;
          return false;
        }
        if (mat.commGrid().id2D() == submat.getRankId2D(i, j)) {
          if (mat.ptr(subindex.first + i, subindex.second + j) != submat.ptr(i, j)) {
            std::cout << "Pointers of element (" << i << ", " << j << ") are not the same."
                      << std::endl;
            return false;
          }
        }
      }
    }
    return true;
  }
}

#endif  // DLA_INTERFACE_TEST_INCLUDE_UTIL_DISTRIBUTED_MATRIX_H
