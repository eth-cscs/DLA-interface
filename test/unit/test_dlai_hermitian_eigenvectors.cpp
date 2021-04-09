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

bool hermitianEigenvectorsTestThrows(SolverType solver) {
#ifdef DLAI_WITH_SCALAPACK
  if (solver == ScaLAPACK)
    return false;
#endif
#ifdef DLAI_WITH_ELPA
  if (solver == ELPA)
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

  void testHermitianEigenvectors(UpLo uplo, DistributedMatrix<T>& a, DistributedMatrix<T>& evec,
                                 DistributedMatrix<T>& tmp, SolverType solver) {
    if (hermitianEigenvectorsTestThrows(solver)) {
      std::vector<BaseType<T>> eval;
      EXPECT_THROW(hermitianEigenvectors(uplo, a, eval, evec, solver), std::invalid_argument);
      return;
    }

    int n = a.size().first;
    auto eval_expected = [n](int i) {
      // E.g. n == 4: {-1, -.3333, .3333, 1}
      //      n == 5: {-1, -.5, 0, .5, 1}
      return 2. * i / (n - 1) - 1.;
    };

    auto el_val = [this, n = a.size().first, &eval_expected](int i, int j) {
      // Antidiagonal symmetric matrix.
      if (i + j == n - 1)
        return this->value(eval_expected(std::max(i, j)), 0);
      return this->value(0, 0);
    };
    auto el_val_uplo = [this, uplo, &el_val](int i, int j) {
      // Set the non used elements to -1.
      if (i < j && uplo == Lower)
        return this->value(-1, 0);
      if (i > j && uplo == Upper)
        return this->value(-1, 0);
      return el_val(i, j);
    };

    fillDistributedMatrix(a, el_val_uplo);

    std::vector<BaseType<T>> eval;
    hermitianEigenvectors(uplo, a, eval, evec, solver);

    EXPECT_EQ(a.size().first, static_cast<int>(eval.size()));
    bool eval_check = true;
    for (size_t i = 0; i < eval.size(); ++i) {
      if (std::abs(eval_expected(i) - eval[i]) > 1000 * this->epsilon()) {
        *outstream << "eval element (" << i << ") is wrong." << std::endl;
        *outstream << "Expected " << eval_expected(i) << "," << std::endl;
        *outstream << "got " << eval[i] << "." << std::endl;
        *outstream << "Difference: " << std::abs(eval[i] - eval_expected(i)) << "." << std::endl;
        eval_check = false;
        break;
      }
    }
    EXPECT_TRUE(eval_check);
    // Skip eigenvectors check if the eigenvalues are wrong.
    if (!eval_check)
      return;

    // reset a to the original matrix.
    fillDistributedMatrix(a, el_val);

    matrixMultiplication(NoTrans, NoTrans, this->value(1, 0), a, evec, this->value(0, 0), tmp,
                         ScaLAPACK);
    auto ev_val = [this, &eval](int i, int j) {
      // diagonal matrix with eigenvalues.
      if (i == j)
        return this->value(eval[i], 0);
      return this->value(0, 0);
    };
    // set a to the eigenvalue matrix.
    fillDistributedMatrix(a, ev_val);

    matrixMultiplication(NoTrans, NoTrans, this->value(-1, 0), evec, a, this->value(1, 0), tmp,
                         ScaLAPACK);

    auto expected = [this](int, int) { return this->value(0, 0); };
    EXPECT_TRUE(checkNearDistributedMatrix(tmp, expected, n * this->epsilon(), 1., *outstream));
  }
};

typedef ::testing::Types<float, double, std::complex<float>, std::complex<double>> MyTypes;
TYPED_TEST_CASE(DLATypedTest, MyTypes);

TYPED_TEST(DLATypedTest, HermitianEigenvectors) {
  using ElType = TypeParam;
  int nb = 3;

  for (auto n : {54, 59, 62})
    for (auto comm_ptr : comms) {
      for (auto dist : DISTRIBUTION_SET) {
        DistributedMatrix<ElType> a(n, n, nb, nb, *comm_ptr, dist);
        DistributedMatrix<ElType> evec(n, n, nb, nb, *comm_ptr, dist);
        DistributedMatrix<ElType> tmp(n, n, nb, nb, *comm_ptr, dist);
        for (auto solver : SOLVER_SET) {
          for (auto uplo : UPLO_SET) {
            this->testHermitianEigenvectors(uplo, a, evec, tmp, solver);
          }
        }
      }
    }
}
