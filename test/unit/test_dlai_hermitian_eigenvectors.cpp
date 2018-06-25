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

bool hermitianEigenvectorsTestThrows(SolverType solver) {
#ifdef DLA_HAVE_SCALAPACK
  if (solver == ScaLAPACK)
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
                                 SolverType solver) {
    if (hermitianEigenvectorsTestThrows(solver)) {
      std::vector<BaseType<T>> eval;
      EXPECT_THROW(hermitianEigenvectors(uplo, a, eval, evec, solver), std::invalid_argument);
      return;
    }

    auto eval_expected = [n = a.size().first](int i) {
      // E.g. n == 4: {1, .3333, -.3333, 1}
      //      n == 5: {1, .5, 0, -.5, 1}
      return 2. * i / (n - 1) - 1.;
    };

    auto el_val = [ this, n = a.size().first, &eval_expected ](int i, int j) {
      // Antidiagonal symmetric matrix.
      if (i + j == n - 1)
        return this->value(eval_expected(std::max(i, j)), 0);
      return this->value(0, 0);
    };

    fillDistributedMatrix(a, el_val);

    std::vector<BaseType<T>> eval;
    hermitianEigenvectors(uplo, a, eval, evec, solver);

    EXPECT_EQ(a.size().first, eval.size());
    bool eval_check = true;
    for (int i = 0; i < eval.size(); ++i) {
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
  }

  // TODO Eigenvectors check. (Note: may be opposite vector.)
};

typedef ::testing::Types<float, double, std::complex<float>, std::complex<double>> MyTypes;
TYPED_TEST_CASE(DLATypedTest, MyTypes);

TYPED_TEST(DLATypedTest, HermitianEigenvectors) {
  using ElType = TypeParam;
  int nb = 3;

  for (auto n : {54, 59, 62})
    for (auto comm_ptr : comms) {
      for (auto dist : DISTRIBUTION_SET) {
        for (auto solver : SOLVER_SET) {
          for (auto uplo : UPLO_SET) {
            DistributedMatrix<ElType> a(n, n, nb, nb, *comm_ptr, dist);
            DistributedMatrix<ElType> evec(n, n, nb, nb, *comm_ptr, dist);

            this->testHermitianEigenvectors(uplo, a, evec, solver);
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

  outstream = &::testing::setMPIListener("results_test_dlai_hermitian_eigenvectors");

  auto ret = RUN_ALL_TESTS();

  comm::CommunicatorManager::finalize();

  return ret;
}
