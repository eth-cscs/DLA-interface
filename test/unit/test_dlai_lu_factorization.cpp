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

bool LUFactorizationTestThrows(SolverType solver) {
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

  void testLUFactorization(DistributedMatrix<T>& a, std::vector<int>& ipiv, SolverType solver) {
    if (LUFactorizationTestThrows(solver)) {
      EXPECT_THROW(LUFactorization(a, ipiv, solver), std::invalid_argument);
      return;
    }

    auto el_val = [this](int i, int j) {
      return this->value(
          std::exp2(-(i + 2 * j) + 3 * (std::min(i, j) + 1)) / 7 - std::exp2(-(i + 2 * j)) / 7,
          -i + j);
    };

    auto el_val_expected = [this](int i, int j) {
      return i < j ? this->value(std::exp2(-2 * std::abs(i - j)), -i + j)
                   : this->value(std::exp2(-std::abs(i - j)), -i + j);
    };

    fillDistributedMatrix(a, el_val);

    LUFactorization(a, ipiv, solver);

    EXPECT_TRUE(
        checkNearDistributedMatrix(a, el_val_expected, 10 * this->epsilon(), 1e-3, *outstream));

    bool ipiv_check = true;
    for (int i = 0; i < a.localSize().first; ++i) {
      int global_row_index = a.getGlobal2DIndex(Local2DIndex(i, 0)).row;
      if (global_row_index < a.size().second) {
        int expected = global_row_index + a.baseIndex().row + 1;
        if (expected != ipiv[a.localBaseIndex().row + i]) {
          *outstream << "ipiv element (" << i << ") is wrong." << std::endl;
          *outstream << "Expected " << expected << "," << std::endl;
          *outstream << "got " << ipiv[a.localBaseIndex().row + i] << "." << std::endl;
          ipiv_check = false;
          break;
        }
      }
    }
    EXPECT_TRUE(ipiv_check);
  }
};

typedef ::testing::Types<float, double, std::complex<float>, std::complex<double>> MyTypes;
TYPED_TEST_CASE(DLATypedTest, MyTypes);

TYPED_TEST(DLATypedTest, LUFactorization) {
  using ElType = TypeParam;
  int n = 59;
  int nb = 3;

  for (auto comm_ptr : comms) {
    for (auto dist : DISTRIBUTION_SET) {
      for (auto solver : SOLVER_SET) {
        for (int m : {47, 59, 79}) {
          DistributedMatrix<ElType> a(m, n, nb, nb, *comm_ptr, dist);
          std::vector<int> ipiv;

          this->testLUFactorization(a, ipiv, solver);
        }
      }
    }
  }
}

TYPED_TEST(DLATypedTest, LUFactorizationSubmatrixSame) {
  using ElType = TypeParam;
  int n = 59;
  int nb = 3;

  for (auto comm_ptr : comms) {
    for (auto dist : DISTRIBUTION_SET) {
      for (auto solver : SOLVER_SET) {
        for (int m : {47, 59, 79}) {
          auto a_ptr =
              DistributedMatrix<ElType>(m + nb, n + nb, nb, nb, *comm_ptr, dist).subMatrix(m, n, nb, nb);
          std::vector<int> ipiv;

          this->testLUFactorization(*a_ptr, ipiv, solver);
        }
      }
    }
  }
}

TYPED_TEST(DLATypedTest, LUFactorizationSubmatrix) {
  using ElType = TypeParam;
  int n = 59;
  int nb = 3;

  for (auto comm_ptr : comms) {
    for (auto dist : DISTRIBUTION_SET) {
      for (auto solver : SOLVER_SET) {
        for (int m : {47, 59, 79}) {
          auto a_ptr = DistributedMatrix<ElType>(m + 2 * nb, n + nb, nb, nb, *comm_ptr, dist)
                           .subMatrix(m, n, 2 * nb, nb);
          std::vector<int> ipiv;

          this->testLUFactorization(*a_ptr, ipiv, solver);
        }
      }
    }
  }
}
