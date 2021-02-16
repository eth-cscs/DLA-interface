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

bool triangularInverseTestThrows(SolverType solver) {
#ifdef DLAI_WITH_SCALAPACK
  if (solver == ScaLAPACK)
    return false;
#endif
#ifdef DLAI_WITH_DPLASMA
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

  void testTriangularInverse(UpLo uplo, Diag diag, DistributedMatrix<T>& a, SolverType solver) {
    if (triangularInverseTestThrows(solver)) {
      EXPECT_THROW(triangularInverse(uplo, diag, a, solver), std::invalid_argument);
      return;
    }

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

    triangularInverse(uplo, diag, a, solver);

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

TYPED_TEST(DLATypedTest, TriangularInverse) {
  using ElType = TypeParam;
  int n = 59;
  int nb = 3;

  for (auto comm_ptr : comms) {
    for (auto dist : DISTRIBUTION_SET) {
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

TYPED_TEST(DLATypedTest, TriangularInverseSubmatrixSame) {
  using ElType = TypeParam;
  int n = 59;
  int nb = 3;

  for (auto comm_ptr : comms) {
    for (auto dist : DISTRIBUTION_SET) {
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

TYPED_TEST(DLATypedTest, TriangularInverseSubmatrix) {
  using ElType = TypeParam;
  int n = 59;
  int nb = 3;

  for (auto comm_ptr : comms) {
    for (auto dist : DISTRIBUTION_SET) {
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
