//
// Distributed Linear Algebra Interface (DLAI)
//
// Copyright (c) 2018-2019, ETH Zurich
// All rights reserved.
//
// Please, refer to the LICENSE file in the root directory.
// SPDX-License-Identifier: BSD-3-Clause
//

#ifndef DLA_INTERFACE_C_DLA_INTERFACE_H
#define DLA_INTERFACE_C_DLA_INTERFACE_H

#include <mpi.h>

	/** DLAINT C Interface.
	 *
	 * The DLAINT C Interface is a layer in DLAINT to enable communication with C compilers.
	 *
	 * Note the routines for the complex case are declared using the pointer
	 * to the real part of the first element.
	 */

#ifdef __cplusplus
extern "C" {
#endif
void dlai_initialize(const int* nr_cores, const int* initialize_mpi);
void dlai_initialize_arg(const int* nr_cores, int* argc, char*** argv, const int* initialize_mpi);

void dlai_finalize();

int dlai_create_2d_grid(const MPI_Fint* base_comm, const int* row_size, const int* col_size,
                        const char* ordering);

int dlai_create_2d_grid_blacs(int* blacs_handle, const int* row_size, const int* col_size,
                              const char* ordering);

void dlai_free_2d_grid_blacs(int* blacs_context);

int dlai_get_print_timer_option();
void dlai_set_print_timer_option(const int* print_timer);

#define DLA_DECLARE_CHOLESKY_FACTORIZATION(function_name, Type)                            \
  int function_name(const char* uplo, const int* n, Type* a, const int* ia, const int* ja, \
                    const int* desca, const char* solver)

DLA_DECLARE_CHOLESKY_FACTORIZATION(dlai_s_cholesky_factorization, float);
DLA_DECLARE_CHOLESKY_FACTORIZATION(dlai_d_cholesky_factorization, double);
DLA_DECLARE_CHOLESKY_FACTORIZATION(dlai_c_cholesky_factorization, float);
DLA_DECLARE_CHOLESKY_FACTORIZATION(dlai_z_cholesky_factorization, double);

#define DLA_DECLARE_CHOLESKY_INVERSE(function_name, Type)                                  \
  int function_name(const char* uplo, const int* n, Type* a, const int* ia, const int* ja, \
                    const int* desca, const char* solver)

DLA_DECLARE_CHOLESKY_INVERSE(dlai_s_cholesky_inverse, float);
DLA_DECLARE_CHOLESKY_INVERSE(dlai_d_cholesky_inverse, double);
DLA_DECLARE_CHOLESKY_INVERSE(dlai_c_cholesky_inverse, float);
DLA_DECLARE_CHOLESKY_INVERSE(dlai_z_cholesky_inverse, double);

#define DLA_DECLARE_LU_FACTORIZATION(function_name, Type)                              \
  int function_name(const int* m, const int* n, Type* a, const int* ia, const int* ja, \
                    const int* desca, int* ipiv, const char* solver)

DLA_DECLARE_LU_FACTORIZATION(dlai_s_lu_factorization, float);
DLA_DECLARE_LU_FACTORIZATION(dlai_d_lu_factorization, double);
DLA_DECLARE_LU_FACTORIZATION(dlai_c_lu_factorization, float);
DLA_DECLARE_LU_FACTORIZATION(dlai_z_lu_factorization, double);

#define DLA_DECLARE_MATRIX_MULTIPLICATION(function_name, Type)                                    \
  int function_name(const char* transa, const char* transb, const int* m, const int* n,           \
                    const int* k, const Type* alpha, const Type* a, const int* ia, const int* ja, \
                    const int* desca, const Type* b, const int* ib, const int* jb,                \
                    const int* descb, const Type* beta, Type* c, const int* ic, const int* jc,    \
                    const int* descc, const char* solver)

DLA_DECLARE_MATRIX_MULTIPLICATION(dlai_s_matrix_multiplication, float);
DLA_DECLARE_MATRIX_MULTIPLICATION(dlai_d_matrix_multiplication, double);
DLA_DECLARE_MATRIX_MULTIPLICATION(dlai_c_matrix_multiplication, float);
DLA_DECLARE_MATRIX_MULTIPLICATION(dlai_z_matrix_multiplication, double);

#define DLA_DECLARE_HERMITIAN_EIGENVECTORS(function_name, Type, RealType)                     \
  int function_name(const char* uplo, const int* n, Type* a, const int* ia, const int* ja,    \
                    const int* desca, RealType* evals, Type* v, const int* iv, const int* jv, \
                    const int* descv, const char* solver);

DLA_DECLARE_HERMITIAN_EIGENVECTORS(dlai_s_hermitian_eigenvectors, float, float)
DLA_DECLARE_HERMITIAN_EIGENVECTORS(dlai_d_hermitian_eigenvectors, double, double)
DLA_DECLARE_HERMITIAN_EIGENVECTORS(dlai_c_hermitian_eigenvectors, float, float)
DLA_DECLARE_HERMITIAN_EIGENVECTORS(dlai_z_hermitian_eigenvectors, double, double)

#define DLA_DECLARE_TRIANGULAR_INVERSE(function_name, Type)                                   \
  int function_name(const char* uplo, const char* diag, const int* n, Type* a, const int* ia, \
                    const int* ja, const int* desca, const char* solver)

DLA_DECLARE_TRIANGULAR_INVERSE(dlai_s_triangular_inverse, float);
DLA_DECLARE_TRIANGULAR_INVERSE(dlai_d_triangular_inverse, double);
DLA_DECLARE_TRIANGULAR_INVERSE(dlai_c_triangular_inverse, float);
DLA_DECLARE_TRIANGULAR_INVERSE(dlai_z_triangular_inverse, double);
#ifdef __cplusplus
}
#endif

#endif  // DLA_INTERFACE_C_DLA_INTERFACE_H
