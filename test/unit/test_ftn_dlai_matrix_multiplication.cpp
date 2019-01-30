#include "dla_interface.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include "gtest/gtest.h"
#include "communicator_grid.h"
#include "communicator_manager.h"
#include "types.h"
#include "util_complex.h"
#include "util_distributed_matrix.h"
#include "test_ftn_dlai_main.h"

using namespace dla_interface;
using namespace testing;

#define TEST_FTN_MATRIX_MULTIPLICATION(function_name)                                           \
  void function_name(const char* transa, const char* transb, const int* m, const int* n,        \
                     const int* k, void* alpha, void* a, const int* ia, const int* ja,          \
                     const int* desca, void* b, const int* ib, const int* jb, const int* descb, \
                     void* beta, void* c, const int* ic, const int* jc, const int* descc,       \
                     const char* solver, int* info);
extern "C" {
TEST_FTN_MATRIX_MULTIPLICATION(test_ftn_s_matrix_multiplication);
TEST_FTN_MATRIX_MULTIPLICATION(test_ftn_d_matrix_multiplication);
TEST_FTN_MATRIX_MULTIPLICATION(test_ftn_c_matrix_multiplication);
TEST_FTN_MATRIX_MULTIPLICATION(test_ftn_z_matrix_multiplication);
}

bool matrixMultiplicationTestThrows(SolverType solver) {
#ifdef DLA_HAVE_SCALAPACK
  if (solver == ScaLAPACK)
    return false;
#endif
#ifdef DLA_HAVE_DPLASMA
  if (solver == DPlasma)
    return false;
#endif
  return true;
}

template <typename T>
class FtnDLATypedTest : public ::testing::Test {
  public:
  BaseType<T> epsilon() {
    return std::numeric_limits<BaseType<T>>::epsilon();
  }

  T value(BaseType<T> r, int arg) {
    return testing::value<T>(r, arg);
  }

  void testMatrixMultiplication(OpTrans trans_a, OpTrans trans_b, DistributedMatrix<T>& a,
                                DistributedMatrix<T>& b, DistributedMatrix<T>& c, SolverType solver) {
    T alpha(2);
    T beta(-1);

    bool swap_index_a = trans_a == NoTrans ? false : true;
    int exp_complex_a = trans_a == ConjTrans ? -1 : 1;

    auto el_val_a = [this, swap_index_a, exp_complex_a](int i, int j) {
      if (swap_index_a)
        std::swap(i, j);
      return this->value(static_cast<BaseType<T>>(i + 1) / (j + 1), exp_complex_a * (-i + j));
    };

    bool swap_index_b = trans_b == NoTrans ? false : true;
    int exp_complex_b = trans_b == ConjTrans ? -1 : 1;

    auto el_val_b = [this, swap_index_b, exp_complex_b](int i, int j) {
      if (swap_index_b)
        std::swap(i, j);
      return this->value(static_cast<BaseType<T>>(i + 1) / (j + 1), exp_complex_b * (-i + j));
    };

    auto el_val_c = [this](int i, int j) {
      return this->value(static_cast<BaseType<T>>(i + 1) / (j + 1), -i + j);
    };

    int k = trans_a == NoTrans ? a.size().second : a.size().first;
    auto el_val_c_expected = [this, k, alpha, beta](int i, int j) {
      return this->value(static_cast<BaseType<T>>(i + 1) / (j + 1), -i + j) *
             (alpha * static_cast<BaseType<T>>(k) + beta);
    };

    fillDistributedMatrix(a, el_val_a);
    fillDistributedMatrix(b, el_val_b);
    fillDistributedMatrix(c, el_val_c);

    auto a_scalapack = a.getScalapackDescription();
    auto b_scalapack = b.getScalapackDescription();
    auto c_scalapack = c.getScalapackDescription();
    int m = c.size().first;
    int n = c.size().second;
    const char c_trans_a = static_cast<char>(trans_a);
    const char c_trans_b = static_cast<char>(trans_b);
    auto solver_s = util::getSolverString(solver);
    int info = 0;

    if (std::is_same<T, float>::value)
      test_ftn_s_matrix_multiplication(
          &c_trans_a, &c_trans_b, &m, &n, &k, &alpha, std::get<0>(a_scalapack),
          &std::get<1>(a_scalapack), &std::get<2>(a_scalapack), &std::get<3>(a_scalapack)[0],
          std::get<0>(b_scalapack), &std::get<1>(b_scalapack), &std::get<2>(b_scalapack),
          &std::get<3>(b_scalapack)[0], &beta, std::get<0>(c_scalapack), &std::get<1>(c_scalapack),
          &std::get<2>(c_scalapack), &std::get<3>(c_scalapack)[0], solver_s.c_str(), &info);
    if (std::is_same<T, double>::value)
      test_ftn_d_matrix_multiplication(
          &c_trans_a, &c_trans_b, &m, &n, &k, &alpha, std::get<0>(a_scalapack),
          &std::get<1>(a_scalapack), &std::get<2>(a_scalapack), &std::get<3>(a_scalapack)[0],
          std::get<0>(b_scalapack), &std::get<1>(b_scalapack), &std::get<2>(b_scalapack),
          &std::get<3>(b_scalapack)[0], &beta, std::get<0>(c_scalapack), &std::get<1>(c_scalapack),
          &std::get<2>(c_scalapack), &std::get<3>(c_scalapack)[0], solver_s.c_str(), &info);
    if (std::is_same<T, std::complex<float>>::value)
      test_ftn_c_matrix_multiplication(
          &c_trans_a, &c_trans_b, &m, &n, &k, &alpha, std::get<0>(a_scalapack),
          &std::get<1>(a_scalapack), &std::get<2>(a_scalapack), &std::get<3>(a_scalapack)[0],
          std::get<0>(b_scalapack), &std::get<1>(b_scalapack), &std::get<2>(b_scalapack),
          &std::get<3>(b_scalapack)[0], &beta, std::get<0>(c_scalapack), &std::get<1>(c_scalapack),
          &std::get<2>(c_scalapack), &std::get<3>(c_scalapack)[0], solver_s.c_str(), &info);
    if (std::is_same<T, std::complex<double>>::value)
      test_ftn_z_matrix_multiplication(
          &c_trans_a, &c_trans_b, &m, &n, &k, &alpha, std::get<0>(a_scalapack),
          &std::get<1>(a_scalapack), &std::get<2>(a_scalapack), &std::get<3>(a_scalapack)[0],
          std::get<0>(b_scalapack), &std::get<1>(b_scalapack), &std::get<2>(b_scalapack),
          &std::get<3>(b_scalapack)[0], &beta, std::get<0>(c_scalapack), &std::get<1>(c_scalapack),
          &std::get<2>(c_scalapack), &std::get<3>(c_scalapack)[0], solver_s.c_str(), &info);

    if (matrixMultiplicationTestThrows(solver)) {
      EXPECT_NE(0, info);
      return;
    }

    EXPECT_TRUE(
        checkNearDistributedMatrix(c, el_val_c_expected, k * this->epsilon(), 1e-3, *outstream));
  }
};

typedef ::testing::Types<float, double, std::complex<float>, std::complex<double>> MyTypes;
TYPED_TEST_CASE(FtnDLATypedTest, MyTypes);

TYPED_TEST(FtnDLATypedTest, MatrixMultiplication) {
  using ElType = TypeParam;
  int m = 49;
  int n = 59;
  int k = 37;
  int nb = 3;

  for (auto comm_ptr : comms) {
    for (auto dist : {scalapack_dist}) {
      for (auto solver : SOLVER_SET) {
        for (auto trans_a : OPTRANS_SET) {
          for (auto trans_b : OPTRANS_SET) {
            int a_m = trans_a == NoTrans ? m : k;
            int a_n = trans_a == NoTrans ? k : m;
            int b_m = trans_b == NoTrans ? k : n;
            int b_n = trans_b == NoTrans ? n : k;

            DistributedMatrix<ElType> a(a_m, a_n, nb, nb, *comm_ptr, dist);
            DistributedMatrix<ElType> b(b_m, b_n, nb, nb, *comm_ptr, dist);
            DistributedMatrix<ElType> c(m, n, nb, nb, *comm_ptr, dist);

            std::stringstream s;
            s << "Case: " << trans_a << trans_b << " Solver: " << solver << " Distribution: " << dist;
            SCOPED_TRACE(s.str());

            this->testMatrixMultiplication(trans_a, trans_b, a, b, c, solver);
          }
        }
      }
    }
  }
}

TYPED_TEST(FtnDLATypedTest, MatrixMultiplicationDiffBlockSize) {
  using ElType = TypeParam;
  int m = 49;
  int n = 59;
  int k = 37;
  int opa_mb = 5;
  int opa_nb = 3;
  int opb_mb = 5;
  int opb_nb = 4;
  int c_mb = 4;
  int c_nb = 3;

  for (auto comm_ptr : comms) {
    for (auto dist : {scalapack_dist}) {
      for (auto solver : SOLVER_SET) {
        for (auto trans_a : OPTRANS_SET) {
          for (auto trans_b : OPTRANS_SET) {
            int a_m = trans_a == NoTrans ? m : k;
            int a_n = trans_a == NoTrans ? k : m;
            int a_mb = trans_a == NoTrans ? opa_mb : opa_nb;
            int a_nb = trans_a == NoTrans ? opa_nb : opa_mb;
            int b_m = trans_b == NoTrans ? k : n;
            int b_n = trans_b == NoTrans ? n : k;
            int b_mb = trans_b == NoTrans ? opb_mb : opb_nb;
            int b_nb = trans_b == NoTrans ? opb_nb : opb_mb;

            DistributedMatrix<ElType> a(a_m, a_n, a_mb, a_nb, *comm_ptr, dist);
            DistributedMatrix<ElType> b(b_m, b_n, b_mb, b_nb, *comm_ptr, dist);
            DistributedMatrix<ElType> c(m, n, c_mb, c_nb, *comm_ptr, dist);

            std::stringstream s;
            s << "Case: " << trans_a << trans_b << " Solver: " << solver << " Distribution: " << dist;
            SCOPED_TRACE(s.str());

            this->testMatrixMultiplication(trans_a, trans_b, a, b, c, solver);
          }
        }
      }
    }
  }
}

TYPED_TEST(FtnDLATypedTest, MatrixMultiplicationSubMatrix) {
  using ElType = TypeParam;
  int m = 49;
  int n = 59;
  int k = 37;
  int nb = 3;

  for (auto comm_ptr : comms) {
    for (auto dist : {scalapack_dist}) {
      for (auto solver : SOLVER_SET) {
        for (auto trans_a : OPTRANS_SET) {
          for (auto trans_b : OPTRANS_SET) {
            int a_m = trans_a == NoTrans ? m : k;
            int a_n = trans_a == NoTrans ? k : m;
            int b_m = trans_b == NoTrans ? k : n;
            int b_n = trans_b == NoTrans ? n : k;

            auto a_ptr = DistributedMatrix<ElType>(a_m + 2 * nb, a_n + nb, nb, nb, *comm_ptr, dist)
                             .subMatrix(a_m, a_n, 2 * nb, nb);
            auto b_ptr =
                DistributedMatrix<ElType>(b_m + 3 * nb, b_n + 2 * nb, nb, nb, *comm_ptr, dist)
                    .subMatrix(b_m, b_n, 3 * nb, 2 * nb);
            auto c_ptr =
                DistributedMatrix<ElType>(m + nb, n, nb, nb, *comm_ptr, dist).subMatrix(m, n, nb, 0);

            std::stringstream s;
            s << "Case: " << trans_a << trans_b << " Solver: " << solver << " Distribution: " << dist;
            SCOPED_TRACE(s.str());

            this->testMatrixMultiplication(trans_a, trans_b, *a_ptr, *b_ptr, *c_ptr, solver);
          }
        }
      }
    }
  }
}
