#ifndef DLA_INTERFACE_TEST_INCLUDE_UTIL_LOCAL_MATRIX_H
#define DLA_INTERFACE_TEST_INCLUDE_UTIL_LOCAL_MATRIX_H

namespace testing {
  using namespace dla_interface;
  template <class Type, class F>
  void fillLocalMatrix(LocalMatrix<Type>& mat, F& el_val) {
    for (int j = 0; j < mat.size().second; ++j) {
      for (int i = 0; i < mat.size().first; ++i) {
        mat(i, j) = el_val(i, j);
      }
    }
  }

  template <class Type, class F>
  bool checkLocalMatrix(const LocalMatrix<Type>& mat, F& el_val) {
    for (int j = 0; j < mat.size().second; ++j) {
      for (int i = 0; i < mat.size().first; ++i) {
        if (mat(i, j) != el_val(i, j)) {
          std::cout << "Element (" << i << ", " << j << ") is wrong." << std::endl;
          std::cout << "Expected " << el_val(i, j) << "," << std::endl;
          std::cout << "Got " << mat(i, j) << "." << std::endl;
          return false;
        }
      }
    }
    return true;
  }
}

#endif  // DLA_INTERFACE_TEST_INCLUDE_UTIL_LOCAL_MATRIX_H
