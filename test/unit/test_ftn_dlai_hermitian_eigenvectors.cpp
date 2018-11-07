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
#include "test_ftn_dlai_main.h"

using namespace dla_interface;
using namespace testing;

#define TEST_FTN_HERMITIAN_EIGENVECTORS(function_name)                                      \
  void function_name(const char* uplo, const int* n, void* a, const int* ia, const int* ja, \
                     const int* desca, void* evals, void* v, const int* iv, const int* jv,  \
                     const int* descv, const char* solver, int* info)

extern "C" {
TEST_FTN_HERMITIAN_EIGENVECTORS(test_ftn_s_hermitian_eigenvectors);
TEST_FTN_HERMITIAN_EIGENVECTORS(test_ftn_d_hermitian_eigenvectors);
TEST_FTN_HERMITIAN_EIGENVECTORS(test_ftn_c_hermitian_eigenvectors);
TEST_FTN_HERMITIAN_EIGENVECTORS(test_ftn_z_hermitian_eigenvectors);
}

bool hermitianEigenvectorsTestThrows(SolverType solver) {
#ifdef DLA_HAVE_SCALAPACK
  if (solver == ScaLAPACK)
    return false;
#endif
#ifdef DLA_HAVE_ELPA
  if (solver == ELPA)
    return false;
#endif
  return true;
}

template <typename T>
class FtnDLATypedTest : public ::testing::Test {
  public:
  BaseType<T> epsilon() {
    return std::numeric_limits<BaseType<T>>::epsilon();
  }

  T value(BaseType<T> r, int arg) {
    return testing::value<T>(r, arg);
  }

  void testHermitianEigenvectors(UpLo uplo, DistributedMatrix<T>& a, DistributedMatrix<T>& evec,
                                 DistributedMatrix<T>& tmp, SolverType solver) {
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

    auto a_scalapack = a.getScalapackDescription();
    auto v_scalapack = evec.getScalapackDescription();
    std::vector<BaseType<T>> eval(n);
    const char c_uplo = static_cast<char>(uplo);
    auto solver_s = util::getSolverString(solver);
    int info = 0;

    if (std::is_same<T, float>::value)
      test_ftn_s_hermitian_eigenvectors(
          &c_uplo, &n, std::get<0>(a_scalapack), &std::get<1>(a_scalapack),
          &std::get<2>(a_scalapack), &std::get<3>(a_scalapack)[0], &eval[0],
          std::get<0>(v_scalapack), &std::get<1>(v_scalapack), &std::get<2>(v_scalapack),
          &std::get<3>(v_scalapack)[0], solver_s.c_str(), &info);
    if (std::is_same<T, double>::value)
      test_ftn_d_hermitian_eigenvectors(
          &c_uplo, &n, std::get<0>(a_scalapack), &std::get<1>(a_scalapack),
          &std::get<2>(a_scalapack), &std::get<3>(a_scalapack)[0], &eval[0],
          std::get<0>(v_scalapack), &std::get<1>(v_scalapack), &std::get<2>(v_scalapack),
          &std::get<3>(v_scalapack)[0], solver_s.c_str(), &info);
    if (std::is_same<T, std::complex<float>>::value)
      test_ftn_c_hermitian_eigenvectors(
          &c_uplo, &n, std::get<0>(a_scalapack), &std::get<1>(a_scalapack),
          &std::get<2>(a_scalapack), &std::get<3>(a_scalapack)[0], &eval[0],
          std::get<0>(v_scalapack), &std::get<1>(v_scalapack), &std::get<2>(v_scalapack),
          &std::get<3>(v_scalapack)[0], solver_s.c_str(), &info);
    if (std::is_same<T, std::complex<double>>::value)
      test_ftn_z_hermitian_eigenvectors(
          &c_uplo, &n, std::get<0>(a_scalapack), &std::get<1>(a_scalapack),
          &std::get<2>(a_scalapack), &std::get<3>(a_scalapack)[0], &eval[0],
          std::get<0>(v_scalapack), &std::get<1>(v_scalapack), &std::get<2>(v_scalapack),
          &std::get<3>(v_scalapack)[0], solver_s.c_str(), &info);

    if (hermitianEigenvectorsTestThrows(solver)) {
      EXPECT_NE(0, info);
      return;
    }

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
TYPED_TEST_CASE(FtnDLATypedTest, MyTypes);

TYPED_TEST(FtnDLATypedTest, HermitianEigenvectors) {
  using ElType = TypeParam;
  int nb = 3;

  for (auto n : {54, 59, 62})
    for (auto comm_ptr : comms) {
      for (auto dist : {scalapack_dist}) {
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
