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

bool choleskyInverseTestThrows(SolverType solver) {
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
class DLATypedTest : public ::testing::Test {
  public:
  BaseType<T> epsilon() {
    return std::numeric_limits<BaseType<T>>::epsilon();
  }

  T value(BaseType<T> r, int arg) {
    return testing::value<T>(r, arg);
  }

  void testCholeskyInverse(UpLo uplo, DistributedMatrix<T>& a, SolverType solver) {
    if (choleskyInverseTestThrows(solver)) {
      EXPECT_THROW(choleskyInverse(uplo, a, solver), std::invalid_argument);
      return;
    }

    auto el_val = [this](int i, int j) { return this->value(std::exp2(-std::abs(i - j)), -i + j); };

    int n = a.size().first;
    auto el_val_expected = [this, n](int i, int j) {
      if (i == j) {
        if (i == n - 1)
          return this->value(1, 0);
        return this->value(1.25, 0);
      }
      if (i - j == 1)
        return this->value(-0.5, -1);
      if (i - j == -1)
        return this->value(-0.5, 1);
      return this->value(0, 0);
    };

    auto el_val_expected_lower = [&el_val, &el_val_expected](int i, int j) {
      return i < j ? el_val(i, j) : el_val_expected(i, j);
    };
    auto el_val_expected_upper = [&el_val, &el_val_expected](int i, int j) {
      return i > j ? el_val(i, j) : el_val_expected(i, j);
    };

    fillDistributedMatrix(a, el_val);

    choleskyInverse(uplo, a, solver);

    if (uplo == Lower)
      EXPECT_TRUE(checkNearDistributedMatrix(a, el_val_expected_lower, 100 * this->epsilon(), 1e-1,
                                             *outstream));
    if (uplo == Upper)
      EXPECT_TRUE(checkNearDistributedMatrix(a, el_val_expected_upper, 100 * this->epsilon(), 1e-1,
                                             *outstream));
  }
};

typedef ::testing::Types<float, double, std::complex<float>, std::complex<double>> MyTypes;
TYPED_TEST_CASE(DLATypedTest, MyTypes);

TYPED_TEST(DLATypedTest, CholeskyInverse) {
  using ElType = TypeParam;
  int n = 59;
  int nb = 3;

  for (auto comm_ptr : comms) {
    for (auto dist : DISTRIBUTION_SET) {
      for (auto solver : SOLVER_SET) {
        for (auto uplo : UPLO_SET) {
          DistributedMatrix<ElType> a(n, n, nb, nb, *comm_ptr, dist);

          this->testCholeskyInverse(uplo, a, solver);
        }
      }
    }
  }
}

TYPED_TEST(DLATypedTest, CholeskyInverseSubmatrixSame) {
  using ElType = TypeParam;
  int n = 59;
  int nb = 3;

  for (auto comm_ptr : comms) {
    for (auto dist : DISTRIBUTION_SET) {
      for (auto solver : SOLVER_SET) {
        for (auto uplo : UPLO_SET) {
          auto a_ptr =
              DistributedMatrix<ElType>(n + nb, n + nb, nb, nb, *comm_ptr, dist).subMatrix(n, n, nb, nb);

          this->testCholeskyInverse(uplo, *a_ptr, solver);
        }
      }
    }
  }
}

TYPED_TEST(DLATypedTest, CholeskyInverseSubmatrix) {
  using ElType = TypeParam;
  int n = 59;
  int nb = 3;

  for (auto comm_ptr : comms) {
    for (auto dist : DISTRIBUTION_SET) {
      for (auto solver : SOLVER_SET) {
        for (auto uplo : UPLO_SET) {
          auto a_ptr = DistributedMatrix<ElType>(n + 2 * nb, n + nb, nb, nb, *comm_ptr, dist)
                           .subMatrix(n, n, 2 * nb, nb);

          this->testCholeskyInverse(uplo, *a_ptr, solver);
        }
      }
    }
  }
}
