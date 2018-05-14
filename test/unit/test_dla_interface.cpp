#include "dla_interface.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include "gtest/gtest.h"
#include "mpi_listener.h"
#include "communicator_grid.h"
#include "communicator_manager.h"
#include "util_complex.h"
#include "util_distributed_matrix.h"

using namespace dla_interface;
using namespace testing;

constexpr auto dists = {scalapack_dist, tile_dist};
constexpr auto solvers = {ScaLAPACK, ELPA, DPlasma, Chameleon, HPX_LINALG};
std::vector<comm::Communicator2DGrid*> comms;

template <typename T>
class DLATypedTest : public ::testing::Test {
  public:
  BaseType<T> epsilon() {
    return std::numeric_limits<BaseType<T>>::epsilon();
  }

  T value(BaseType<T> r, int arg) {
    return testing::value<T>(r, arg);
  }
};

typedef ::testing::Types<float, double, std::complex<float>, std::complex<double>> MyTypes;
TYPED_TEST_CASE(DLATypedTest, MyTypes);

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

TYPED_TEST(DLATypedTest, LUFactorization) {
  using ElType = TypeParam;
  int n = 59;
  int nb = 3;

  auto el_val = [this](int i, int j) {
    return this->value(
        std::exp2(-(i + 2 * j) + 3 * (std::min(i, j) + 1)) / 7 - std::exp2(-(i + 2 * j)) / 7, -i + j);
  };
  auto el_val_expected = [this, &el_val](int i, int j) {
    return i < j ? this->value(std::exp2(-2 * std::abs(i - j)), -i + j)
                 : this->value(std::exp2(-std::abs(i - j)), -i + j);
  };
  for (int m : {47, 59, 79}) {
    for (auto comm_ptr : comms) {
      for (auto dist : dists) {
        for (auto solver : solvers) {
          auto A1 = std::make_shared<DistributedMatrix<ElType>>(m, n, nb, nb, *comm_ptr, dist);
          auto A2 = DistributedMatrix<ElType>(m + nb, n + 2 * nb, nb, nb, *comm_ptr, dist)
                        .subMatrix(m, n, nb, 2 * nb);
          for (auto A_ptr : {A1, A2}) {
            auto& A = *A_ptr;
            fillDistributedMatrix(A, el_val);
            std::vector<int> ipiv;

            if (LUFactorizationTestThrows(solver)) {
              EXPECT_THROW(LUFactorization(A, ipiv, solver), std::invalid_argument);
            }
            else {
              LUFactorization(A, ipiv, solver);

              EXPECT_TRUE(checkNearDistributedMatrix(A, el_val_expected, 10 * this->epsilon()));
              for (int i = 0; i < A.localSize().first; ++i) {
                int global_row_index = A.getGlobal2DIndex(Local2DIndex(i, 0)).row;
                if (global_row_index < A.size().second) {
                  int expected = global_row_index + A.baseIndex().row + 1;
                  EXPECT_EQ(expected, ipiv[A.localBaseIndex().row + i]);
                }
              }
            }
          }
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

  ::testing::setMPIListener("results_test_dla_interface");

  auto ret = RUN_ALL_TESTS();

  comm::CommunicatorManager::finalize();

  return ret;
}
