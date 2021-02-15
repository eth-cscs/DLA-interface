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

#define TEST_FTN_TRIANGULAR_INVERSE(function_name)                                             \
  void function_name(const char* uplo, const char* diag, const int* n, void* a, const int* ia, \
                     const int* ja, const int* desca, const char* solver, int* info)

extern "C" {
TEST_FTN_TRIANGULAR_INVERSE(test_ftn_s_triangular_inverse);
TEST_FTN_TRIANGULAR_INVERSE(test_ftn_d_triangular_inverse);
TEST_FTN_TRIANGULAR_INVERSE(test_ftn_c_triangular_inverse);
TEST_FTN_TRIANGULAR_INVERSE(test_ftn_z_triangular_inverse);
}

bool triangularInverseTestThrows(SolverType solver) {
#ifdef DLAI_WITH_SCALAPACK
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

  void testTriangularInverse(UpLo uplo, Diag diag, DistributedMatrix<T>& a, SolverType solver) {
    auto el_val = [this, diag](int i, int j) {
      if (diag == Unit)
        return this->value(i == j ? 10 : std::exp2(-std::abs(i - j)), -i + j);
      return this->value(-std::exp2(-std::abs(i - j)), -i + j);
    };

    auto el_val_expected = [this, diag, &el_val](int i, int j) {
      if (diag == Unit) {
        if (i == j)
          return el_val(i, j);
        if (i - j == 1)
          return this->value(-0.5, -1);
        if (i - j == -1)
          return this->value(-0.5, 1);
        return this->value(0, 0);
      }
      if (i == j)
        return this->value(-1, 0);
      if (i - j == 1)
        return this->value(0.5, -1);
      if (i - j == -1)
        return this->value(0.5, 1);
      return this->value(0, 0);
    };

    auto el_val_expected_lower = [&el_val, &el_val_expected](int i, int j) {
      return i < j ? el_val(i, j) : el_val_expected(i, j);
    };
    auto el_val_expected_upper = [&el_val, &el_val_expected](int i, int j) {
      return i > j ? el_val(i, j) : el_val_expected(i, j);
    };

    fillDistributedMatrix(a, el_val);

    int n = a.size().first;
    auto a_scalapack = a.getScalapackDescription();
    const char c_uplo = static_cast<char>(uplo);
    const char c_diag = static_cast<char>(diag);
    auto solver_s = util::getSolverString(solver);
    int info = 0;

    if (std::is_same<T, float>::value)
      test_ftn_s_triangular_inverse(&c_uplo, &c_diag, &n, std::get<0>(a_scalapack),
                                    &std::get<1>(a_scalapack), &std::get<2>(a_scalapack),
                                    &std::get<3>(a_scalapack)[0], solver_s.c_str(), &info);
    if (std::is_same<T, double>::value)
      test_ftn_d_triangular_inverse(&c_uplo, &c_diag, &n, std::get<0>(a_scalapack),
                                    &std::get<1>(a_scalapack), &std::get<2>(a_scalapack),
                                    &std::get<3>(a_scalapack)[0], solver_s.c_str(), &info);
    if (std::is_same<T, std::complex<float>>::value)
      test_ftn_c_triangular_inverse(&c_uplo, &c_diag, &n, std::get<0>(a_scalapack),
                                    &std::get<1>(a_scalapack), &std::get<2>(a_scalapack),
                                    &std::get<3>(a_scalapack)[0], solver_s.c_str(), &info);
    if (std::is_same<T, std::complex<double>>::value)
      test_ftn_z_triangular_inverse(&c_uplo, &c_diag, &n, std::get<0>(a_scalapack),
                                    &std::get<1>(a_scalapack), &std::get<2>(a_scalapack),
                                    &std::get<3>(a_scalapack)[0], solver_s.c_str(), &info);

    if (triangularInverseTestThrows(solver)) {
      EXPECT_NE(0, info);
      return;
    }

    EXPECT_EQ(0, info);
    if (uplo == Lower)
      EXPECT_TRUE(checkNearDistributedMatrix(a, el_val_expected_lower, 100 * this->epsilon(), 1e-1,
                                             *outstream));
    if (uplo == Upper)
      EXPECT_TRUE(checkNearDistributedMatrix(a, el_val_expected_upper, 100 * this->epsilon(), 1e-1,
                                             *outstream));
  }
};

typedef ::testing::Types<float, double, std::complex<float>, std::complex<double>> MyTypes;
TYPED_TEST_CASE(FtnDLATypedTest, MyTypes);

TYPED_TEST(FtnDLATypedTest, TriangularInverse) {
  using ElType = TypeParam;
  int n = 59;
  int nb = 3;

  for (auto comm_ptr : comms) {
    // Fortran interface only for ScaLAPACK matrices.
    for (auto dist : {scalapack_dist}) {
      for (auto solver : SOLVER_SET) {
        for (auto uplo : UPLO_SET) {
          for (auto diag : DIAG_SET) {
            DistributedMatrix<ElType> a(n, n, nb, nb, *comm_ptr, dist);

            this->testTriangularInverse(uplo, diag, a, solver);
          }
        }
      }
    }
  }
}

TYPED_TEST(FtnDLATypedTest, TriangularInverseSubmatrixSame) {
  using ElType = TypeParam;
  int n = 59;
  int nb = 3;

  for (auto comm_ptr : comms) {
    for (auto dist : {scalapack_dist}) {
      for (auto solver : SOLVER_SET) {
        for (auto uplo : UPLO_SET) {
          for (auto diag : DIAG_SET) {
            auto a_ptr = DistributedMatrix<ElType>(n + nb, n + nb, nb, nb, *comm_ptr, dist)
                             .subMatrix(n, n, nb, nb);

            this->testTriangularInverse(uplo, diag, *a_ptr, solver);
          }
        }
      }
    }
  }
}

TYPED_TEST(FtnDLATypedTest, TriangularInverseSubmatrix) {
  using ElType = TypeParam;
  int n = 59;
  int nb = 3;

  for (auto comm_ptr : comms) {
    for (auto dist : {scalapack_dist}) {
      for (auto solver : SOLVER_SET) {
        for (auto uplo : UPLO_SET) {
          for (auto diag : DIAG_SET) {
            auto a_ptr = DistributedMatrix<ElType>(n + 2 * nb, n + nb, nb, nb, *comm_ptr, dist)
                             .subMatrix(n, n, 2 * nb, nb);

            this->testTriangularInverse(uplo, diag, *a_ptr, solver);
          }
        }
      }
    }
  }
}
