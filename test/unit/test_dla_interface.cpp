#include "dla_interface.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include "gtest/gtest.h"
#include "mpi_listener.h"
#include "util_complex.h"
#include "util_distributed_matrix.h"
#include "communicator_grid.h"
#include "communicator_manager.h"

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
    for (auto dist : dists) {
      for (auto solver : solvers) {
        for (auto uplo : {Lower, Upper}) {
          auto A1 = std::make_shared<DistributedMatrix<ElType>>(n, n, nb, nb, *comm_ptr, dist);
          auto A2 =
              DistributedMatrix<ElType>(n + nb, n + nb, nb, nb, *comm_ptr, dist).subMatrix(n, n, nb, nb);
          std::vector<std::shared_ptr<DistributedMatrix<ElType>>> mats{A1};
          if (solver != HPX_LINALG) mats.push_back(A2);
          if (solver == HPX_LINALG && uplo != Lower) continue;
          for (auto A_ptr : mats) {
            auto& A = *A_ptr;
            fillDistributedMatrix(A, el_val);

            if (choleskyFactorizationTestThrows(solver)) {
              EXPECT_THROW(choleskyFactorization(uplo, A, solver), std::invalid_argument);
            }
            else {
              choleskyFactorization(uplo, A, solver);

              if (uplo == Lower)
                EXPECT_TRUE(
                    checkNearDistributedMatrix(A, el_val_expected_lower, 10 * this->epsilon()));
              if (uplo == Upper)
                EXPECT_TRUE(
                    checkNearDistributedMatrix(A, el_val_expected_upper, 10 * this->epsilon()));
            }
          }
        }
      }
    }
  }
}

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

bool matrixMultiplicationTestThrows(SolverType solver) {
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

TYPED_TEST(DLATypedTest, Gemm) {
  using ElType = TypeParam;
  int m = 49;
  int n = 59;
  int k = 37;
  int nb = 3;

  ElType alpha(2);
  ElType beta(-1);

  auto el_val_c = [this](int i, int j) {
    return this->value(static_cast<BaseType<ElType>>(i + 1) / (j + 1), -i + j);
  };
  auto el_val_c_expected = [this, k, alpha, beta](int i, int j) {
    return this->value(static_cast<BaseType<ElType>>(i + 1) / (j + 1), -i + j) *
           (alpha * static_cast<BaseType<ElType>>(k) + beta);
  };

  for (auto trans_a : {NoTrans, Trans, ConjTrans}) {
    int a_m = trans_a == NoTrans ? m : k;
    int a_n = trans_a == NoTrans ? k : m;
    bool swap_index_a = trans_a == NoTrans ? false : true;
    int exp_complex_a = trans_a == ConjTrans ? -1 : 1;

    auto el_val_a = [this, swap_index_a, exp_complex_a](int i, int j) {
      if (swap_index_a)
        std::swap(i, j);
      return this->value(static_cast<BaseType<ElType>>(i + 1) / (j + 1), exp_complex_a * (-i + j));
    };
    for (auto trans_b : {NoTrans, Trans, ConjTrans}) {
      int b_m = trans_b == NoTrans ? k : n;
      int b_n = trans_b == NoTrans ? n : k;
      bool swap_index_b = trans_b == NoTrans ? false : true;
      int exp_complex_b = trans_b == ConjTrans ? -1 : 1;

      auto el_val_b = [this, swap_index_b, exp_complex_b](int i, int j) {
        if (swap_index_b)
          std::swap(i, j);
        return this->value(static_cast<BaseType<ElType>>(i + 1) / (j + 1), exp_complex_b * (-i + j));
      };

      for (auto comm_ptr : comms) {
        for (auto dist : dists) {
          for (auto solver : solvers) {
            auto a_ptr = DistributedMatrix<ElType>(a_m + nb, a_n + nb, nb, nb, *comm_ptr, dist)
                             .subMatrix(a_m, a_n, nb, nb);
            auto b_ptr = DistributedMatrix<ElType>(b_m + 2 * nb, b_n, nb, nb, *comm_ptr, dist)
                             .subMatrix(b_m, b_n, 2 * nb, 0);
            auto c_ptr =
                DistributedMatrix<ElType>(m, n + nb, nb, nb, *comm_ptr, dist).subMatrix(m, n, 0, nb);

            auto& a = *a_ptr;
            auto& b = *b_ptr;
            auto& c = *c_ptr;

            fillDistributedMatrix(a, el_val_a);
            fillDistributedMatrix(b, el_val_b);
            fillDistributedMatrix(c, el_val_c);

            if (matrixMultiplicationTestThrows(solver)) {
              EXPECT_THROW(matrixMultiplication(trans_a, trans_b, alpha, a, b, beta, c, solver),
                           std::invalid_argument);
            }
            else {
              matrixMultiplication(trans_a, trans_b, alpha, a, b, beta, c, solver);

              EXPECT_TRUE(checkNearDistributedMatrix(c, el_val_c_expected, k * this->epsilon()));
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
