#include <mpi.h>
#include "c_dla_interface.h"
#include "dla_interface.h"
#include "util_types.h"
#include "distributed_matrix.h"

// Note the routines for the complex case are declared using the pointer to the real part of the
// first element.

using namespace dla_interface;

int dlai_print_timer_value = 0;

extern "C" void dlai_initialize(const int* nr_cores, const int* initialize_mpi) {
  dlai_initialize_arg(nr_cores, nullptr, nullptr, initialize_mpi);
}
extern "C" void dlai_initialize_arg(const int* nr_cores, int* argc, char*** argv,
                                    const int* initialize_mpi) {
  comm::CommunicatorManager::initialize(*nr_cores, argc, argv, *initialize_mpi);
}

extern "C" void dlai_finalize() {
  comm::CommunicatorManager::finalize();
}

extern "C" int dlai_create_2d_grid(const MPI_Fint* base_comm, const int* row_size,
                                   const int* col_size, const char* ordering) {
  return comm::CommunicatorManager::createCommunicator2DGrid(
             MPI_Comm_f2c(*base_comm), *row_size, *col_size, util::getOrdering(*ordering))
      .blacsContext();
}

extern "C" int dlai_create_2d_grid_blacs(int* blacs_handle, const int* row_size,
                                         const int* col_size, const char* ordering) {
  return comm::CommunicatorManager::createCommunicator2DGridBlacs(
             *blacs_handle, *row_size, *col_size, util::getOrdering(*ordering))
      .blacsContext();
}

extern "C" void dlai_free_2d_grid_blacs(int* blacs_context) {
  return comm::CommunicatorManager::free2DGridFromBlacsContext(*blacs_context);
}

extern "C" int dlai_get_print_timer_option() {
  return dlai_print_timer_value;
}

extern "C" void dlai_set_print_timer_option(const int* print_timer) {
  dlai_print_timer_value = *print_timer;
}

#define DLA_DEFINE_CHOLESKY_FACTORIZATION(function_name, CType, CppType)                \
  extern "C" int function_name(const char* uplo, const int* n, CType* a, const int* ia, \
                               const int* ja, const int* desca, const char* solver) {   \
    CppType* a_ = reinterpret_cast<CppType*>(a);                                        \
    UpLo uplo_ = util::getUpLo(*uplo);                                                  \
    DistributedMatrix<CppType> mat_a(scalapack_dist, *n, *n, a_, *ia, *ja, desca);      \
    SolverType solver_ = util::getSolverType(solver);                                   \
    try {                                                                               \
      choleskyFactorization(uplo_, mat_a, solver_, dlai_print_timer_value);             \
    }                                                                                   \
    catch (std::invalid_argument & exc) {                                               \
      return -1;                                                                        \
    }                                                                                   \
    return 0;                                                                           \
  }

DLA_DEFINE_CHOLESKY_FACTORIZATION(dlai_s_cholesky_factorization, float, float)
DLA_DEFINE_CHOLESKY_FACTORIZATION(dlai_d_cholesky_factorization, double, double)
DLA_DEFINE_CHOLESKY_FACTORIZATION(dlai_c_cholesky_factorization, float, std::complex<float>)
DLA_DEFINE_CHOLESKY_FACTORIZATION(dlai_z_cholesky_factorization, double, std::complex<double>)

#define DLA_DEFINE_LU_FACTORIZATION(function_name, CType, CppType)                                 \
  extern "C" int function_name(const int* m, const int* n, CType* a, const int* ia, const int* ja, \
                               const int* desca, int* ipiv, const char* solver) {                  \
    CppType* a_ = reinterpret_cast<CppType*>(a);                                                   \
    DistributedMatrix<CppType> mat_a(scalapack_dist, *m, *n, a_, *ia, *ja, desca);                 \
    SolverType solver_ = util::getSolverType(solver);                                              \
    try {                                                                                          \
      LUFactorization(mat_a, ipiv, solver_, dlai_print_timer_value);                               \
    }                                                                                              \
    catch (std::invalid_argument & exc) {                                                          \
      return -1;                                                                                   \
    }                                                                                              \
    return 0;                                                                                      \
  }

DLA_DEFINE_LU_FACTORIZATION(dlai_s_lu_factorization, float, float)
DLA_DEFINE_LU_FACTORIZATION(dlai_d_lu_factorization, double, double)
DLA_DEFINE_LU_FACTORIZATION(dlai_c_lu_factorization, float, std::complex<float>)
DLA_DEFINE_LU_FACTORIZATION(dlai_z_lu_factorization, double, std::complex<double>)

#define DLA_DEFINE_MATRIX_MULTIPLICATION(function_name, CType, CppType)                            \
  extern "C" int function_name(                                                                    \
      const char* trans_a, const char* trans_b, const int* m, const int* n, const int* k,          \
      const CType* alpha, const CType* a, const int* ia, const int* ja, const int* desca,          \
      const CType* b, const int* ib, const int* jb, const int* descb, const CType* beta, CType* c, \
      const int* ic, const int* jc, const int* descc, const char* solver) {                        \
    const CppType* alpha_ = reinterpret_cast<const CppType*>(alpha);                               \
    const CppType* beta_ = reinterpret_cast<const CppType*>(beta);                                 \
    const CppType* a_ = reinterpret_cast<const CppType*>(a);                                       \
    const CppType* b_ = reinterpret_cast<const CppType*>(b);                                       \
    CppType* c_ = reinterpret_cast<CppType*>(c);                                                   \
    OpTrans trans_a_ = util::getOpTrans(*trans_a);                                                 \
    OpTrans trans_b_ = util::getOpTrans(*trans_b);                                                 \
    int m_a = trans_a_ == NoTrans ? *m : *k;                                                       \
    int n_a = trans_a_ == NoTrans ? *k : *m;                                                       \
    int m_b = trans_b_ == NoTrans ? *k : *n;                                                       \
    int n_b = trans_b_ == NoTrans ? *n : *k;                                                       \
    auto mat_a_ptr =                                                                               \
        DistributedMatrix<CppType>::convertConst(scalapack_dist, m_a, n_a, a_, *ia, *ja, desca);   \
    auto mat_b_ptr =                                                                               \
        DistributedMatrix<CppType>::convertConst(scalapack_dist, m_b, n_b, b_, *ib, *jb, descb);   \
    DistributedMatrix<CppType> mat_c(scalapack_dist, *m, *n, c_, *ic, *jc, descc);                 \
    SolverType solver_ = util::getSolverType(solver);                                              \
    try {                                                                                          \
      matrixMultiplication(trans_a_, trans_b_, *alpha_, *mat_a_ptr, *mat_b_ptr, *beta_, mat_c,     \
                           solver_, dlai_print_timer_value);                                       \
    }                                                                                              \
    catch (std::invalid_argument & exc) {                                                          \
      return -1;                                                                                   \
    }                                                                                              \
    return 0;                                                                                      \
  }

DLA_DEFINE_MATRIX_MULTIPLICATION(dlai_s_matrix_multiplication, float, float)
DLA_DEFINE_MATRIX_MULTIPLICATION(dlai_d_matrix_multiplication, double, double)
DLA_DEFINE_MATRIX_MULTIPLICATION(dlai_c_matrix_multiplication, float, std::complex<float>)
DLA_DEFINE_MATRIX_MULTIPLICATION(dlai_z_matrix_multiplication, double, std::complex<double>)
