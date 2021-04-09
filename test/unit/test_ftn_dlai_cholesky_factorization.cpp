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

#define TEST_FTN_CHOLESKY_FACTORIZATION(function_name)                                      \
  void function_name(const char* uplo, const int* n, void* a, const int* ia, const int* ja, \
                     const int* desca, const char* solver, int* info)

extern "C" {
TEST_FTN_CHOLESKY_FACTORIZATION(test_ftn_s_cholesky_factorization);
TEST_FTN_CHOLESKY_FACTORIZATION(test_ftn_d_cholesky_factorization);
TEST_FTN_CHOLESKY_FACTORIZATION(test_ftn_c_cholesky_factorization);
TEST_FTN_CHOLESKY_FACTORIZATION(test_ftn_z_cholesky_factorization);
}

bool choleskyFactorizationTestThrows(SolverType solver) {
#ifdef DLAI_WITH_SCALAPACK
  if (solver == ScaLAPACK)
    return false;
#endif
#ifdef DLAI_WITH_DPLASMA
  if (solver == DPlasma)
    return false;
#endif
#ifdef DLAI_WITH_HPX_LINALG
  if (solver == HPX_LINALG)
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

  void testCholeskyFactorization(UpLo uplo, DistributedMatrix<T>& a, SolverType solver) {
    auto el_val = [this](int i, int j) {
      return this->value(
          std::exp2(-(i + j) + 2 * (std::min(i, j) + 1)) / 3 - std::exp2(-(i + j)) / 3, -i + j);
    };

    auto el_val_expected_lower = [this, &el_val](int i, int j) {
      return i < j ? el_val(i, j) : this->value(std::exp2(-std::abs(i - j)), -i + j);
    };
    auto el_val_expected_upper = [this, &el_val](int i, int j) {
      return i > j ? el_val(i, j) : this->value(std::exp2(-std::abs(i - j)), -i + j);
    };

    fillDistributedMatrix(a, el_val);

    auto a_scalapack = a.getScalapackDescription();
    int n = a.size().first;
    const char c_uplo = static_cast<char>(uplo);
    auto solver_s = util::getSolverString(solver);
    int info = 0;

    if (std::is_same<T, float>::value)
      test_ftn_s_cholesky_factorization(&c_uplo, &n, std::get<0>(a_scalapack),
                                        &std::get<1>(a_scalapack), &std::get<2>(a_scalapack),
                                        &std::get<3>(a_scalapack)[0], solver_s.c_str(), &info);
    if (std::is_same<T, double>::value)
      test_ftn_d_cholesky_factorization(&c_uplo, &n, std::get<0>(a_scalapack),
                                        &std::get<1>(a_scalapack), &std::get<2>(a_scalapack),
                                        &std::get<3>(a_scalapack)[0], solver_s.c_str(), &info);
    if (std::is_same<T, std::complex<float>>::value)
      test_ftn_c_cholesky_factorization(&c_uplo, &n, std::get<0>(a_scalapack),
                                        &std::get<1>(a_scalapack), &std::get<2>(a_scalapack),
                                        &std::get<3>(a_scalapack)[0], solver_s.c_str(), &info);
    if (std::is_same<T, std::complex<double>>::value)
      test_ftn_z_cholesky_factorization(&c_uplo, &n, std::get<0>(a_scalapack),
                                        &std::get<1>(a_scalapack), &std::get<2>(a_scalapack),
                                        &std::get<3>(a_scalapack)[0], solver_s.c_str(), &info);

    if (choleskyFactorizationTestThrows(solver)) {
      EXPECT_NE(0, info);
      return;
    }

    EXPECT_EQ(0, info);
    if (uplo == Lower)
      EXPECT_TRUE(checkNearDistributedMatrix(a, el_val_expected_lower, 10 * this->epsilon(), 1e-3,
                                             *outstream));
    if (uplo == Upper)
      EXPECT_TRUE(checkNearDistributedMatrix(a, el_val_expected_upper, 10 * this->epsilon(), 1e-3,
                                             *outstream));
  }
};

typedef ::testing::Types<float, double, std::complex<float>, std::complex<double>> MyTypes;
TYPED_TEST_CASE(FtnDLATypedTest, MyTypes);

TYPED_TEST(FtnDLATypedTest, CholeskyFactorization) {
  using ElType = TypeParam;
  int n = 59;
  int nb = 3;

  for (auto comm_ptr : comms) {
    // Fortran interface only for ScaLAPACK matrices.
    for (auto dist : {scalapack_dist}) {
      for (auto solver : SOLVER_SET) {
        for (auto uplo : UPLO_SET) {
          DistributedMatrix<ElType> a(n, n, nb, nb, *comm_ptr, dist);

          this->testCholeskyFactorization(uplo, a, solver);
        }
      }
    }
  }
}

TYPED_TEST(FtnDLATypedTest, CholeskyFactorizationSubmatrixSame) {
  using ElType = TypeParam;
  int n = 59;
  int nb = 3;

  for (auto comm_ptr : comms) {
    for (auto dist : {scalapack_dist}) {
      for (auto solver : SOLVER_SET) {
        for (auto uplo : UPLO_SET) {
          auto a_ptr =
              DistributedMatrix<ElType>(n + nb, n + nb, nb, nb, *comm_ptr, dist).subMatrix(n, n, nb, nb);

          this->testCholeskyFactorization(uplo, *a_ptr, solver);
        }
      }
    }
  }
}

TYPED_TEST(FtnDLATypedTest, CholeskyFactorizationSubmatrix) {
  using ElType = TypeParam;
  int n = 59;
  int nb = 3;

  for (auto comm_ptr : comms) {
    for (auto dist : {scalapack_dist}) {
      for (auto solver : SOLVER_SET) {
        for (auto uplo : UPLO_SET) {
          auto a_ptr = DistributedMatrix<ElType>(n + 2 * nb, n + nb, nb, nb, *comm_ptr, dist)
                           .subMatrix(n, n, 2 * nb, nb);

          this->testCholeskyFactorization(uplo, *a_ptr, solver);
        }
      }
    }
  }
}
