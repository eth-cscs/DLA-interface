#include "dla_interface.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include "gtest/gtest.h"
#include "mpi_listener.h"
#include "communicator_grid.h"
#include "communicator_manager.h"
#include "types.h"
#include "util_complex.h"
#include "util_distributed_matrix.h"

using namespace dla_interface;
using namespace testing;

std::vector<comm::Communicator2DGrid*> comms;
std::ofstream* outstream;

bool choleskyFactorizationTestThrows(SolverType solver) {
#ifdef DLA_HAVE_SCALAPACK
  if (solver == ScaLAPACK)
    return false;
#endif
#ifdef DLA_HAVE_DPLASMA
  if (solver == DPlasma)
    return false;
#endif
#ifdef DLA_HAVE_HPX_LINALG
  if (solver == HPX_LINALG)
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

int main(int argc, char** argv) {
  comm::CommunicatorManager::initialize(2, &argc, &argv, true);

  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 6) {
    std::cout << "This test need 6 MPI ranks (" << size << " provided)!" << std::endl;
    return 1;
  }

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Create communicators used in the tests.
  for (auto order : {RowMajor, ColMajor}) {
    comms.push_back(&comm::CommunicatorManager::createCommunicator2DGrid(MPI_COMM_WORLD, 2, 3, order));
    comms.push_back(&comm::CommunicatorManager::createCommunicator2DGrid(MPI_COMM_WORLD, 3, 2, order));
    comms.push_back(&comm::CommunicatorManager::createCommunicator2DGrid(MPI_COMM_WORLD, 1, 6, order));
  }

  comm::CommunicatorManager::getFallbackInfo().setFullReport(true);

  ::testing::InitGoogleTest(&argc, argv);

  outstream = &::testing::setMPIListener("results_test_dlai_cholesky_factorization");

  auto ret = RUN_ALL_TESTS();

  comm::CommunicatorManager::finalize();

  return ret;
}
