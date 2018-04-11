#ifndef DLA_INTERFACE_UTIL_CHECK_H
#define DLA_INTERFACE_UTIL_CHECK_H

#include "error_message.h"

namespace dla_interface {
  namespace util {
    template <class Matrix>
    void checkIsSquare(const Matrix& matrix, const char* file, int line, const char* func,
                       const char* matrix_name) {
      if (matrix.size().first != matrix.size().second)
        throw std::invalid_argument(errorMessageInternal(
            file, line, func, matrix_name, " is not square. (Size is ", matrix.size(), ".)"));
    }

    template <class DistMatrix>
    void checkBlocksAreSquare(const DistMatrix& matrix, const char* file, int line,
                              const char* func, const char* matrix_name) {
      if (matrix.blockSize().first != matrix.blockSize().second)
        throw std::invalid_argument(errorMessageInternal(file, line, func, matrix_name,
                                                         " blocks are not square. (Block size is ",
                                                         matrix.blockSize(), ".)"));
    }

    template <class DistMatrix>
    void checkBaseIndexAtBlock(const DistMatrix& matrix, const char* file, int line,
                               const char* func, const char* matrix_name) {
      if (matrix.baseIndex().row % matrix.blockSize().first != 0 ||
          matrix.baseIndex().col % matrix.blockSize().second != 0)
        throw std::invalid_argument(
            errorMessageInternal(file, line, func, matrix_name,
                                 " base index is not aligned with blocks. (Block size is ",
                                 matrix.blockSize(), ", base index is ", matrix.baseIndex(), ")"));
    }

    // Returns 0 if all the matrices have the same communicator or the position (0-based) of the
    // first matrix with a different communicator.
    template <class DistMatrix1, class DistMatrix2>
    int compareComm2D(int n, const DistMatrix1& mat_0, const DistMatrix2& mat_n) {
      if (&mat_0.commGrid() != &mat_n.commGrid())
        return n;
      return 0;
    }
    template <class DistMatrix1, class DistMatrix2, class... Args>
    int compareComm2D(int n, const DistMatrix1& mat_0, const DistMatrix2& mat_n, Args&... args) {
      if (&mat_0.commGrid() != &mat_n.commGrid())
        return n;
      return compareComm2D(++n, mat_0, args...);
    }

    template <class... Args>
    void checkSameComm2D(const char* file, int line, const char* func, Args&... args) {
      int n = compareComm2D(1, args...);

      if (n)
        throw std::invalid_argument(errorMessageInternal(file, line, func, "Matrices 0 and ", n,
                                                         " do not have the same communicator."));
    }
  }
}

// Throws a std::invalid_argument if matrix is not square.
#define dlai__util__checkIsSquare(matrix) \
  dla_interface::util::checkIsSquare(matrix, __FILE__, __LINE__, __func__, #matrix)

#define dlai__util__checkBlocksAreSquare(matrix) \
  dla_interface::util::checkBlocksAreSquare(matrix, __FILE__, __LINE__, __func__, #matrix)

#define dlai__util__checkBaseIndexAtBlock(matrix) \
  dla_interface::util::checkBaseIndexAtBlock(matrix, __FILE__, __LINE__, __func__, #matrix)

#define dlai__util__checkSameComm2D(matrix, ...) \
  dla_interface::util::checkSameComm2D(__FILE__, __LINE__, __func__, matrix, __VA_ARGS__)

#endif  // DLA_INTERFACE_UTIL_CHECK_H
