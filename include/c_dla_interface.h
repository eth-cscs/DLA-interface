#ifndef DLA_INTERFACE_C_DLA_INTERFACE_H
#define DLA_INTERFACE_C_DLA_INTERFACE_H

#include <mpi.h>

// Note the routines for the complex case are declared using the pointer to the real part of the
// first element.

#ifdef __cplusplus
extern "C" {
#endif
void dlai_initialize_(const int* nr_cores, const int* initialize_mpi);
void dlai_initialize_arg_(const int* nr_cores, int* argc, char*** argv, const int* initialize_mpi);

void dlai_finalize_();

int dlai_create_2d_grid_(const MPI_Fint* base_comm, const int* row_size, const int* col_size,
                         const char* ordering);

int dlai_create_2d_grid_blacs_(int* blacs_handle, const int* row_size, const int* col_size,
                               const char* ordering);

void dlai_free_2d_grid_blacs_(int* blacs_context);

int dlai_get_print_timer_option_();
void dlai_set_print_timer_option_(int print_timer);

#define DLA_DECLARE_CHOLESKY_FACTORIZATION(function_name, Type)                            \
  int function_name(const char* uplo, const int* n, Type* a, const int* ia, const int* ja, \
                    const int* desca, const char* solver)

DLA_DECLARE_CHOLESKY_FACTORIZATION(dlai_s_cholesky_factorization_, float);
DLA_DECLARE_CHOLESKY_FACTORIZATION(dlai_d_cholesky_factorization_, double);
DLA_DECLARE_CHOLESKY_FACTORIZATION(dlai_c_cholesky_factorization_, float);
DLA_DECLARE_CHOLESKY_FACTORIZATION(dlai_z_cholesky_factorization_, double);

#define DLA_DECLARE_MATRIX_MULTIPLICATION(function_name, Type)                                    \
  int function_name(const char* transa, const char* transb, const int* m, const int* n,           \
                    const int* k, const Type* alpha, const Type* a, const int* ia, const int* ja, \
                    const int* desca, const Type* b, const int* ib, const int* jb,                \
                    const int* descb, const Type* beta, Type* c, const int* ic, const int* jc,    \
                    const int* descc, const char* solver)

DLA_DECLARE_MATRIX_MULTIPLICATION(dlai_s_matrix_multiplication_, float);
DLA_DECLARE_MATRIX_MULTIPLICATION(dlai_d_matrix_multiplication_, double);
DLA_DECLARE_MATRIX_MULTIPLICATION(dlai_c_matrix_multiplication_, float);
DLA_DECLARE_MATRIX_MULTIPLICATION(dlai_z_matrix_multiplication_, double);

#ifdef __cplusplus
}
#endif

#endif  // DLA_INTERFACE_C_DLA_INTERFACE_H
