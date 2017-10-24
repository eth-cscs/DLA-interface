#include "dla_interface.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include "gtest/gtest.h"
#include "util_complex.h"
#include "util_distributed_matrix.h"
#include "communicator_grid.h"
#include "communicator_manager.h"

using namespace dla_interface;
using namespace testing;

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

TYPED_TEST(DLATypedTest, CholeskyFactorization) {
  using ElType = TypeParam;
  int n = 59;
  int nb = 3;

  auto el_val = [this](int i, int j) {
    return this->value(std::exp2(-(i + j) + 2 * (std::min(i, j) + 1)) / 3 - std::exp2(-(i + j)) / 3,
                       -i + j);
  };
  auto el_val_expected_lower = [this, &el_val](int i, int j) {
    return i < j ? el_val(i, j) : this->value(std::exp2(-std::abs(i - j)), -i + j);
  };
  auto el_val_expected_upper = [this, &el_val](int i, int j) {
    return i > j ? el_val(i, j) : this->value(std::exp2(-std::abs(i - j)), -i + j);
  };

  for (auto comm_ptr : comms) {
    for (auto dist : {scalapack_dist, tile_dist}) {
      for (auto solver : {ScaLAPACK}) {
        for (auto uplo : {Lower, Upper}) {
          DistributedMatrix<ElType> A(n, n, nb, nb, *comm_ptr, dist);
          fillDistributedMatrix(A, el_val);

          choleskyFactorization(uplo, A, solver);
          if (uplo == Lower)
            EXPECT_TRUE(checkNearDistributedMatrix(A, el_val_expected_lower, 10 * this->epsilon()));
          if (uplo == Upper)
            EXPECT_TRUE(checkNearDistributedMatrix(A, el_val_expected_upper, 10 * this->epsilon()));
        }
      }
    }
  }
}

int main(int argc, char** argv) {
  comm::CommunicatorManager::initialize(true);

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
  }

  ::testing::InitGoogleTest(&argc, argv);

  auto ret = RUN_ALL_TESTS();

  comm::CommunicatorManager::finalize();

  return ret;
}
