#include <mpi.h>
#include "c_dla_interface.h"
#include "dla_interface.h"
#include "util_types.h"
#include "distributed_matrix.h"

// Note the routines for the complex case are declared using the pointer to the real part of the
// first element.

using namespace dla_interface;

int dlai_print_timer_value = 0;

extern "C" void dlai_initialize_(const int* nr_cores, const int* initialize_mpi) {
  dlai_initialize_arg_(nr_cores, nullptr, nullptr, initialize_mpi);
}
extern "C" void dlai_initialize_arg_(const int* nr_cores, int* argc, char*** argv,
                                     const int* initialize_mpi) {
  comm::CommunicatorManager::initialize(*nr_cores, argc, argv, *initialize_mpi);
}

extern "C" void dlai_finalize_() {
  comm::CommunicatorManager::finalize();
}

extern "C" int dlai_create_2d_grid_(const MPI_Fint* base_comm, const int* row_size,
                                    const int* col_size, const char* ordering) {
  return comm::CommunicatorManager::createCommunicator2DGrid(
             MPI_Comm_f2c(*base_comm), *row_size, *col_size, util::getOrdering(*ordering))
      .blacsContext();
}

extern "C" int dlai_get_print_timer_option_() {
  return dlai_print_timer_value;
}

extern "C" void dlai_set_print_timer_option_(int print_timer) {
  dlai_print_timer_value = print_timer;
}

#define DLA_DEFINE_CHOLESKY_FACTORIZATION(function_name, CType, CppType)                \
  extern "C" int function_name(const char* uplo, const int* n, CType* a, const int* ia, \
                               const int* ja, int* desca, const char* solver) {         \
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

DLA_DEFINE_CHOLESKY_FACTORIZATION(dlai_s_cholesky_factorization_, float, float)
DLA_DEFINE_CHOLESKY_FACTORIZATION(dlai_d_cholesky_factorization_, double, double)
DLA_DEFINE_CHOLESKY_FACTORIZATION(dlai_c_cholesky_factorization_, float, std::complex<float>)
DLA_DEFINE_CHOLESKY_FACTORIZATION(dlai_z_cholesky_factorization_, double, std::complex<double>)

#define DLA_DEFINE_MATRIX_MULTIPLY(function_name, CType, CppType)                                  \
  extern "C" int function_name(                                                                    \
      const char* trans_a, const char* trans_b, const int* m, const int* n, const int* k,          \
      const CType* alpha, /* const */ CType* a, const int* ia, const int* ja, int* desca,          \
      /* const */ CType* b, const int* ib, const int* jb, int* descb, const CType* beta, CType* c, \
      const int* ic, const int* jc, int* descc, const char* solver) {                              \
    const CppType* alpha_ = reinterpret_cast<const CppType*>(alpha);                               \
    const CppType* beta_ = reinterpret_cast<const CppType*>(beta);                                 \
    CppType* a_ = reinterpret_cast<CppType*>(a);                                                   \
    CppType* b_ = reinterpret_cast<CppType*>(b);                                                   \
    CppType* c_ = reinterpret_cast<CppType*>(c);                                                   \
    OpTrans trans_a_ = util::getOpTrans(*trans_a);                                                 \
    OpTrans trans_b_ = util::getOpTrans(*trans_b);                                                 \
    int m_a = trans_a_ == NoTrans ? *m : *k;                                                       \
    int n_a = trans_a_ == NoTrans ? *k : *m;                                                       \
    int m_b = trans_b_ == NoTrans ? *k : *n;                                                       \
    int n_b = trans_b_ == NoTrans ? *n : *k;                                                       \
    DistributedMatrix<CppType> mat_a(scalapack_dist, m_a, n_a, a_, *ia, *ja, desca);               \
    DistributedMatrix<CppType> mat_b(scalapack_dist, m_b, n_b, b_, *ib, *jb, descb);               \
    DistributedMatrix<CppType> mat_c(scalapack_dist, *m, *n, c_, *ic, *jc, descc);                 \
    SolverType solver_ = util::getSolverType(solver);                                              \
    try {                                                                                          \
      matrixMultiply(trans_a_, trans_b_, *alpha_, mat_a, mat_b, *beta_, mat_c, solver_,            \
                     dlai_print_timer_value);                                                      \
    }                                                                                              \
    catch (std::invalid_argument & exc) {                                                          \
      return -1;                                                                                   \
    }                                                                                              \
    return 0;                                                                                      \
  }

DLA_DEFINE_MATRIX_MULTIPLY(dlai_s_matrix_multiply_, float, float)
DLA_DEFINE_MATRIX_MULTIPLY(dlai_d_matrix_multiply_, double, double)
DLA_DEFINE_MATRIX_MULTIPLY(dlai_c_matrix_multiply_, float, std::complex<float>)
DLA_DEFINE_MATRIX_MULTIPLY(dlai_z_matrix_multiply_, double, std::complex<double>)
