#include "local_matrix.h"

#include <memory>
#include <stdexcept>
#include "gtest/gtest.h"
#include "util_local_matrix.h"

using namespace dla_interface;
using namespace testing;

TEST(LocalMatrixTest, DefaultConstructors) {
  using Type = double;

  LocalMatrix<Type> mat;
  EXPECT_EQ(std::make_pair(0, 0), mat.size());
  EXPECT_EQ(1, mat.leadingDimension());
  EXPECT_EQ(nullptr, mat.ptr());
}

TEST(LocalMatrixTest, Constructor1) {
  using Type = double;

  int m = 8;
  int n = 5;

  LocalMatrix<Type> mat1(m, n);
  EXPECT_EQ(std::make_pair(m, n), mat1.size());
  EXPECT_LE(m, mat1.leadingDimension());
  EXPECT_NE(nullptr, mat1.ptr());

  EXPECT_THROW(LocalMatrix<Type>(-1, n), std::invalid_argument);
  EXPECT_THROW(LocalMatrix<Type>(m, -1), std::invalid_argument);

  // Empty matrix constructors
  LocalMatrix<Type> mat2(0, n);
  EXPECT_EQ(std::make_pair(0, 0), mat2.size());
  EXPECT_EQ(1, mat2.leadingDimension());
  EXPECT_EQ(nullptr, mat2.ptr());

  LocalMatrix<Type> mat3(m, 0);
  EXPECT_EQ(std::make_pair(0, 0), mat3.size());
  EXPECT_EQ(1, mat3.leadingDimension());
  EXPECT_EQ(nullptr, mat3.ptr());
}

TEST(LocalMatrixTest, Constructor2) {
  using Type = double;

  int m = 8;
  int n = 5;
  int ld = 10;

  LocalMatrix<Type> mat1(m, n, m);
  EXPECT_EQ(std::make_pair(m, n), mat1.size());
  EXPECT_EQ(m, mat1.leadingDimension());
  EXPECT_NE(nullptr, mat1.ptr());

  LocalMatrix<Type> mat2(m, n, ld);
  EXPECT_EQ(std::make_pair(m, n), mat2.size());
  EXPECT_EQ(ld, mat2.leadingDimension());
  EXPECT_NE(nullptr, mat2.ptr());

  EXPECT_THROW(LocalMatrix<Type>(m, n, m - 1), std::invalid_argument);
  EXPECT_THROW(LocalMatrix<Type>(-1, n, ld), std::invalid_argument);
  EXPECT_THROW(LocalMatrix<Type>(m, -1, ld), std::invalid_argument);

  // Empty matrix constructors
  LocalMatrix<Type> mat3(0, n, ld);
  EXPECT_EQ(std::make_pair(0, 0), mat3.size());
  EXPECT_EQ(ld, mat3.leadingDimension());
  EXPECT_EQ(nullptr, mat3.ptr());

  LocalMatrix<Type> mat4(m, 0, ld);
  EXPECT_EQ(std::make_pair(0, 0), mat4.size());
  EXPECT_EQ(ld, mat4.leadingDimension());
  EXPECT_EQ(nullptr, mat4.ptr());

  LocalMatrix<Type> mat5(m, 0, 1);
  EXPECT_EQ(std::make_pair(0, 0), mat5.size());
  EXPECT_EQ(1, mat5.leadingDimension());
  EXPECT_EQ(nullptr, mat5.ptr());

  EXPECT_THROW(LocalMatrix<Type>(0, n, 0), std::invalid_argument);
}

TEST(LocalMatrixTest, Constructor3) {
  using Type = double;

  int m = 8;
  int n = 5;
  int ld = 10;
  std::vector<Type> buf(ld * n);

  LocalMatrix<Type> mat1(m, n, &buf[0], ld);
  EXPECT_EQ(std::make_pair(m, n), mat1.size());
  EXPECT_EQ(ld, mat1.leadingDimension());
  EXPECT_EQ(&buf[0], mat1.ptr());

  LocalMatrix<Type> mat2(m, n, &buf[0], m);
  EXPECT_EQ(std::make_pair(m, n), mat2.size());
  EXPECT_EQ(m, mat2.leadingDimension());
  EXPECT_EQ(&buf[0], mat2.ptr());

  EXPECT_THROW(LocalMatrix<Type>(m, n, nullptr, m), std::invalid_argument);
  EXPECT_THROW(LocalMatrix<Type>(m, n, &buf[0], m - 1), std::invalid_argument);
  EXPECT_THROW(LocalMatrix<Type>(-1, n, &buf[0], ld), std::invalid_argument);
  EXPECT_THROW(LocalMatrix<Type>(m, -1, &buf[0], ld), std::invalid_argument);

  // Empty matrix constructors
  for (Type* ptr : {(Type*)nullptr, &buf[0]}) {
    LocalMatrix<Type> mat3(0, n, ptr, ld);
    EXPECT_EQ(std::make_pair(0, 0), mat3.size());
    EXPECT_EQ(ld, mat3.leadingDimension());
    EXPECT_EQ(nullptr, mat3.ptr());

    LocalMatrix<Type> mat4(m, 0, ptr, ld);
    EXPECT_EQ(std::make_pair(0, 0), mat4.size());
    EXPECT_LE(ld, mat4.leadingDimension());
    EXPECT_EQ(nullptr, mat4.ptr());

    LocalMatrix<Type> mat5(m, 0, ptr, 1);
    EXPECT_EQ(std::make_pair(0, 0), mat5.size());
    EXPECT_LE(1, mat5.leadingDimension());
    EXPECT_EQ(nullptr, mat5.ptr());

    EXPECT_THROW(LocalMatrix<Type>(0, n, ptr, 0), std::invalid_argument);
  }
}

TEST(LocalMatrixTest, ElementAndPointer) {
  using Type = double;

  std::vector<LocalMatrix<Type>> vmat;
  vmat.emplace_back(8, 5);
  vmat.emplace_back(12, 7, 14);
  vmat.emplace_back(0, 7);
  vmat.emplace_back(7, 0, 3);
  vmat.emplace_back();

  for (auto& mat : vmat) {
    auto el_val = [](int i, int j) { return j + .001 * i; };
    const auto& const_mat = mat;
    int m = mat.size().first;
    int n = mat.size().second;

    for (int j = 0; j < mat.size().second; ++j) {
      for (int i = 0; i < mat.size().first; ++i) {
        mat(i, j) = el_val(i, j);
      }
    }
    auto ptr = mat.ptr();
    auto const_ptr = const_mat.ptr();
    ASSERT_EQ(ptr, const_ptr);

    for (int j = 0; j < n; ++j) {
      for (int i = 0; i < m; ++i) {
        EXPECT_EQ(el_val(i, j), mat(i, j));
        EXPECT_EQ(ptr + i + j * mat.leadingDimension(), mat.ptr(i, j));
        EXPECT_EQ(mat.ptr(i, j), &mat(i, j));
        EXPECT_EQ(el_val(i, j), const_mat(i, j));
        EXPECT_EQ(const_ptr + i + j * const_mat.leadingDimension(), const_mat.ptr(i, j));
        EXPECT_EQ(const_mat.ptr(i, j), &const_mat(i, j));
      }
    }
    EXPECT_THROW(mat(-1, 0), std::invalid_argument);
    EXPECT_THROW(mat(0, -1), std::invalid_argument);
    EXPECT_THROW(mat(m + 1, 0), std::invalid_argument);
    EXPECT_THROW(mat(0, n + 1), std::invalid_argument);
    EXPECT_THROW(mat.ptr(-1, 0), std::invalid_argument);
    EXPECT_THROW(mat.ptr(0, -1), std::invalid_argument);
    EXPECT_THROW(mat.ptr(m + 1, 0), std::invalid_argument);
    EXPECT_THROW(mat.ptr(0, n + 1), std::invalid_argument);
    EXPECT_THROW(const_mat(-1, 0), std::invalid_argument);
    EXPECT_THROW(const_mat(0, -1), std::invalid_argument);
    EXPECT_THROW(const_mat(m + 1, 0), std::invalid_argument);
    EXPECT_THROW(const_mat(0, n + 1), std::invalid_argument);
    EXPECT_THROW(const_mat.ptr(-1, 0), std::invalid_argument);
    EXPECT_THROW(const_mat.ptr(0, -1), std::invalid_argument);
    EXPECT_THROW(const_mat.ptr(m + 1, 0), std::invalid_argument);
    EXPECT_THROW(const_mat.ptr(0, n + 1), std::invalid_argument);
  }
}

TEST(LocalMatrixTest, CopyConstructor) {
  using Type = double;
  auto el_val = [](int i, int j) { return j + .001 * i; };
  auto el_val2 = [](int i, int j) { return -j + .05 * i; };

  LocalMatrix<Type> mat(8, 5, 13);
  fillLocalMatrix(mat, el_val);
  const auto& const_mat(mat);
  auto size = mat.size();
  auto ld = mat.leadingDimension();

  LocalMatrix<Type> mat2(mat);
  LocalMatrix<Type> mat3(const_mat);

  EXPECT_EQ(size, mat.size());
  EXPECT_EQ(ld, mat.leadingDimension());
  EXPECT_TRUE(checkLocalMatrix(mat, el_val));

  EXPECT_EQ(size, mat2.size());
  EXPECT_EQ(ld, mat2.leadingDimension());
  EXPECT_TRUE(checkLocalMatrix(mat2, el_val));

  EXPECT_EQ(size, mat3.size());
  EXPECT_EQ(ld, mat3.leadingDimension());
  EXPECT_TRUE(checkLocalMatrix(mat3, el_val));

  EXPECT_NE(mat.ptr(), mat2.ptr());
  EXPECT_NE(mat.ptr(), mat3.ptr());
  EXPECT_NE(mat2.ptr(), mat3.ptr());

  // Change the elements of mat and check that the elements of mat2 and mat3 are not changed.
  fillLocalMatrix(mat, el_val2);

  EXPECT_TRUE(checkLocalMatrix(mat2, el_val));
  EXPECT_TRUE(checkLocalMatrix(mat3, el_val));
}

TEST(LocalMatrixTest, AssignementOperator) {
  using Type = double;
  auto el_val = [](int i, int j) { return j + .001 * i; };
  auto el_val2 = [](int i, int j) { return -j + .05 * i; };

  int m = 5;
  int n = 7;
  int ld = 7;

  LocalMatrix<Type> mat1(m, n, ld);
  fillLocalMatrix(mat1, el_val);
  auto size1 = mat1.size();
  auto ld1 = mat1.leadingDimension();

  LocalMatrix<Type> mat2(m, n, ld);
  auto mat2_oldptr = mat2.ptr();
  mat2 = mat1;
  EXPECT_EQ(size1, mat1.size());
  EXPECT_EQ(ld1, mat1.leadingDimension());
  EXPECT_TRUE(checkLocalMatrix(mat1, el_val));
  EXPECT_EQ(ld1, mat2.leadingDimension());
  EXPECT_EQ(size1, mat2.size());
  EXPECT_TRUE(checkLocalMatrix(mat2, el_val));
  EXPECT_NE(mat2_oldptr, mat2.ptr());

  // Change the elements of mat1 and check that the elements of mat2 are not changed.
  fillLocalMatrix(mat1, el_val2);

  EXPECT_EQ(ld1, mat2.leadingDimension());
  EXPECT_EQ(size1, mat2.size());
  EXPECT_TRUE(checkLocalMatrix(mat2, el_val));

  // Check that a reference to the object is returned.
  LocalMatrix<Type> mat3;
  auto pmat3 = &mat3;
  EXPECT_EQ(pmat3, &(mat3 = mat1));
  EXPECT_EQ(ld1, mat3.leadingDimension());
  EXPECT_EQ(size1, mat3.size());
  EXPECT_TRUE(checkLocalMatrix(mat3, el_val2));
}

TEST(LocalMatrixTest, MoveConstructor) {
  using Type = double;
  auto el_val = [](int i, int j) { return j + .001 * i; };
  auto el_val2 = [](int i, int j) { return -j + .05 * i; };

  LocalMatrix<Type> mat1(8, 5, 13);
  fillLocalMatrix(mat1, el_val);
  auto size1 = mat1.size();
  auto ld1 = mat1.leadingDimension();
  auto mat1_oldptr = mat1.ptr();

  LocalMatrix<Type> mat2(std::move(mat1));
  EXPECT_EQ(size1, mat2.size());
  EXPECT_EQ(ld1, mat2.leadingDimension());
  EXPECT_TRUE(checkLocalMatrix(mat2, el_val));
  EXPECT_EQ(mat1_oldptr, mat2.ptr());

  // Re-assign mat1 and check that mat2 do not change.
  mat1 = LocalMatrix<Type>(3, 4, 7);
  fillLocalMatrix(mat1, el_val2);

  EXPECT_EQ(size1, mat2.size());
  EXPECT_EQ(ld1, mat2.leadingDimension());
  EXPECT_TRUE(checkLocalMatrix(mat2, el_val));
}

TEST(LocalMatrixTest, MoveAssignementOperator) {
  using Type = double;
  auto el_val = [](int i, int j) { return j + .001 * i; };
  auto el_val2 = [](int i, int j) { return -j + .05 * i; };

  LocalMatrix<Type> mat1(8, 5, 13);
  fillLocalMatrix(mat1, el_val);
  auto size1 = mat1.size();
  auto ld1 = mat1.leadingDimension();
  auto mat1_oldptr = mat1.ptr();

  LocalMatrix<Type> mat2;
  mat2 = std::move(mat1);
  EXPECT_EQ(ld1, mat2.leadingDimension());
  EXPECT_EQ(size1, mat2.size());
  EXPECT_TRUE(checkLocalMatrix(mat2, el_val));
  EXPECT_EQ(mat1_oldptr, mat2.ptr());

  // Re-assign mat1 and check that mat2 do not change.
  mat1 = LocalMatrix<Type>(3, 4, 7);
  fillLocalMatrix(mat1, el_val2);

  EXPECT_EQ(ld1, mat2.leadingDimension());
  EXPECT_EQ(size1, mat2.size());
  EXPECT_TRUE(checkLocalMatrix(mat2, el_val));

  // Check that a reference to the object is returned.
  LocalMatrix<Type> mat3;
  auto pmat3 = &mat3;
  EXPECT_EQ(pmat3, &(mat3 = std::move(mat1)));
}

TEST(LocalMatrixTest, Copy) {
  using Type = double;
  auto el_val = [](int i, int j) { return j + .001 * i; };
  auto el_val2 = [](int i, int j) { return -j + .05 * i; };

  int m = 8;
  int n = 5;
  LocalMatrix<Type> mat1(m, n, 13);
  std::vector<LocalMatrix<Type>> vmat;
  vmat.emplace_back(m, n, 13);
  vmat.emplace_back(m, n, 17);
  vmat.emplace_back(m, n);

  for (auto& mat : vmat) {
    fillLocalMatrix(mat1, el_val);
    fillLocalMatrix(mat, el_val2);

    mat.copy(mat1);
    fillLocalMatrix(mat1, el_val2);
    EXPECT_TRUE(checkLocalMatrix(mat, el_val));

    // Check that a reference to the object is returned.
    auto pmat = &mat;
    EXPECT_EQ(pmat, &(mat.copy(mat1)));
  }

  LocalMatrix<Type> mat2(m + 1, n);
  EXPECT_THROW(mat2.copy(mat1), std::invalid_argument);
  EXPECT_THROW(mat1.copy(mat2), std::invalid_argument);
  LocalMatrix<Type> mat3(m, n + 1);
  EXPECT_THROW(mat3.copy(mat1), std::invalid_argument);
  EXPECT_THROW(mat1.copy(mat3), std::invalid_argument);
}

TEST(LocalMatrixTest, SubMatrix) {
  using Type = double;
  auto el_val = [](int i, int j) { return j + .001 * i; };
  auto el_val2 = [](int i, int j) { return -j + .05 * i; };

  std::vector<LocalMatrix<Type>> vmat;
  vmat.emplace_back(7, 3);
  vmat.emplace_back(8, 5, 13);

  for (auto& mat : vmat) {
    fillLocalMatrix(mat, el_val);
    int m = mat.size().first;
    int n = mat.size().second;
    const auto& const_mat = mat;
    for (int si = 0; si < m; ++si) {
      for (int sj = 0; sj < n; ++sj) {
        for (int sm = 0; sm < m - si; ++sm) {
          for (int sn = 0; sn < n - sj; ++sn) {
            auto psubmat = mat.subMatrix(sm, sn, si, sj);
            auto pconst_submat = const_mat.subMatrix(sm, sn, si, sj);
            if (sm > 0 && sn > 0) {
              EXPECT_EQ(std::make_pair(sm, sn), psubmat->size());
              EXPECT_EQ(mat.leadingDimension(), psubmat->leadingDimension());
              EXPECT_EQ(std::make_pair(sm, sn), pconst_submat->size());
              EXPECT_EQ(mat.leadingDimension(), pconst_submat->leadingDimension());

              EXPECT_EQ(mat.ptr(si, sj), psubmat->ptr());
              EXPECT_EQ(mat.ptr(si, sj), psubmat->ptr(0, 0));
              EXPECT_EQ(mat.ptr(si, sj), pconst_submat->ptr());
              EXPECT_EQ(mat.ptr(si, sj), pconst_submat->ptr(0, 0));
              Type val = -.123;
              (*psubmat)(0, 0) = val;
              EXPECT_EQ(val, (*pconst_submat)(0, 0));
              EXPECT_EQ(val, mat(si, sj));
            }
            else {
              EXPECT_EQ(std::make_pair(0, 0), psubmat->size());
              EXPECT_EQ(mat.leadingDimension(), psubmat->leadingDimension());
              EXPECT_EQ(nullptr, psubmat->ptr());
              EXPECT_EQ(std::make_pair(0, 0), pconst_submat->size());
              EXPECT_EQ(mat.leadingDimension(), pconst_submat->leadingDimension());
              EXPECT_EQ(nullptr, pconst_submat->ptr());
            }
            fillLocalMatrix(mat, el_val2);
            auto sub_el_val2 = [&el_val2, si, sj](int i, int j) { return el_val2(i + si, j + sj); };
            EXPECT_TRUE(checkLocalMatrix(*psubmat, sub_el_val2));
            EXPECT_TRUE(checkLocalMatrix(*pconst_submat, sub_el_val2));
          }
        }
      }
    }
    EXPECT_THROW(mat.subMatrix(m + 1, n, 0, 0), std::invalid_argument);
    EXPECT_THROW(mat.subMatrix(-1, n, 1, 0), std::invalid_argument);
    EXPECT_THROW(mat.subMatrix(2, n, m - 1, 0), std::invalid_argument);
    EXPECT_THROW(mat.subMatrix(1, n, -1, 0), std::invalid_argument);
    EXPECT_THROW(mat.subMatrix(m, n + 1, 0, 0), std::invalid_argument);
    EXPECT_THROW(mat.subMatrix(m, -1, 0, 1), std::invalid_argument);
    EXPECT_THROW(mat.subMatrix(m, 2, 0, n - 1), std::invalid_argument);
    EXPECT_THROW(mat.subMatrix(m, 1, 0, -1), std::invalid_argument);
  }

  // Check that the submatrix doesn't change if the original matrix is re-assigned.
  int m = 5;
  int n = 3;

  LocalMatrix<Type> mat1(m, n);
  fillLocalMatrix(mat1, el_val);
  auto psubmat1 = mat1.subMatrix(m - 2, n - 2, 0, 0);
  EXPECT_TRUE(checkLocalMatrix(*psubmat1, el_val));

  mat1 = LocalMatrix<Type>(m, n);
  fillLocalMatrix(mat1, el_val2);

  EXPECT_TRUE(checkLocalMatrix(*psubmat1, el_val));

  LocalMatrix<Type> mat2(m, n);
  fillLocalMatrix(mat2, el_val);
  auto psubmat2 = mat2.subMatrix(m - 2, n - 2, 0, 0);
  EXPECT_TRUE(checkLocalMatrix(*psubmat2, el_val));

  LocalMatrix<Type> tmp(m, n);
  mat2 = std::move(tmp);
  fillLocalMatrix(mat2, el_val);

  EXPECT_TRUE(checkLocalMatrix(*psubmat2, el_val));

  // Check that the submatrix is still valid after original matrix goes out of scope.
  std::shared_ptr<LocalMatrix<Type>> psubmat3;
  {
    LocalMatrix<Type> mat3(m, n);
    fillLocalMatrix(mat3, el_val);
    psubmat3 = mat3.subMatrix(m - 2, n - 2, 0, 0);
    EXPECT_TRUE(checkLocalMatrix(*psubmat3, el_val));
  }
  EXPECT_TRUE(checkLocalMatrix(*psubmat3, el_val));
}
