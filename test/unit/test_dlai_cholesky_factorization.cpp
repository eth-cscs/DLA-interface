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
#include "test_dlai_main.h"

using namespace dla_interface;
using namespace testing;

bool choleskyFactorizationTestThrows(SolverType solver) {
#ifdef DLA_HAVE_SCALAPACK
  if (solver == ScaLAPACK)
    return false;
#endif
#ifdef DLA_HAVE_DPLASMA
  if (solver == DPlasma)
    return false;
#endif
#ifdef DLA_HAVE_DLAF
  if (solver == DLAF)
    return false;
#endif
  return true;
}

template <typename T>
class DLATypedTest : public ::testing::Test {
  public:
  BaseType<T> epsilon() {
    return std::numeric_limits<BaseType<T>>::epsilon();
  }

  T value(BaseType<T> r, int arg) {
    return testing::value<T>(r, arg);
  }

  void testCholeskyFactorization(UpLo uplo, DistributedMatrix<T>& a, SolverType solver) {
    if (choleskyFactorizationTestThrows(solver)) {
      EXPECT_THROW(choleskyFactorization(uplo, a, solver), std::invalid_argument);
      return;
    }

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

    choleskyFactorization(uplo, a, solver);

    if (uplo == Lower)
      EXPECT_TRUE(checkNearDistributedMatrix(a, el_val_expected_lower, 10 * this->epsilon(), 1e-3,
                                             *outstream));
    if (uplo == Upper)
      EXPECT_TRUE(checkNearDistributedMatrix(a, el_val_expected_upper, 10 * this->epsilon(), 1e-3,
                                             *outstream));
  }
};

typedef ::testing::Types<float, double, std::complex<float>, std::complex<double>> MyTypes;
TYPED_TEST_CASE(DLATypedTest, MyTypes);

TYPED_TEST(DLATypedTest, CholeskyFactorization) {
  using ElType = TypeParam;
  int n = 59;
  int nb = 3;

  for (auto comm_ptr : comms) {
    for (auto dist : DISTRIBUTION_SET) {
      for (auto solver : SOLVER_SET) {
        for (auto uplo : UPLO_SET) {
          DistributedMatrix<ElType> a(n, n, nb, nb, *comm_ptr, dist);

          this->testCholeskyFactorization(uplo, a, solver);
        }
      }
    }
  }
}

TYPED_TEST(DLATypedTest, CholeskyFactorizationSubmatrixSame) {
  using ElType = TypeParam;
  int n = 59;
  int nb = 3;

  for (auto comm_ptr : comms) {
    for (auto dist : DISTRIBUTION_SET) {
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

TYPED_TEST(DLATypedTest, CholeskyFactorizationSubmatrix) {
  using ElType = TypeParam;
  int n = 59;
  int nb = 3;

  for (auto comm_ptr : comms) {
    for (auto dist : DISTRIBUTION_SET) {
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
