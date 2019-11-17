//
// Distributed Linear Algebra Interface (DLAI)
//
// Copyright (c) 2018-2019, ETH Zurich
// All rights reserved.
//
// Please, refer to the LICENSE file in the root directory.
// SPDX-License-Identifier: BSD-3-Clause
//

#ifndef DLA_INTERFACE_DLA_INTERFACE_H
#define DLA_INTERFACE_DLA_INTERFACE_H

#include <stdexcept>
#include <vector>
#include "distributed_matrix.h"
#include "dla_scalapack.h"
#include "dla_dplasma.h"
#include "dla_elpa.h"
#include "dla_hpx_linalg.h"
#include "error_message.h"
#include "util_cast.h"
#include "util_check.h"
#include "util_debug_output.h"
#include "util_fallback.h"
#include "util_types.h"
#include "timer.h"

namespace dla_interface {

  // print_timers <= 0: do not print any timer info,
  // print_timers >= 1: print routine total time info,
  // print_timers >= 2: print convertion time and solver time as well.
  // Note: enabled timers add MPI barriers.

  template <class ElType>
  void matrixMultiplication(OpTrans trans_a, OpTrans trans_b, ElType alpha,
                            const DistributedMatrix<ElType>& mat_a,
                            const DistributedMatrix<ElType>& mat_b, ElType beta,
                            DistributedMatrix<ElType>& mat_c, SolverType solver,
                            int print_timers = 0);

  template <class ElType>
  void choleskyFactorization(UpLo uplo, DistributedMatrix<ElType>& mat, SolverType solver,
                             int print_timers = 0);

  template <class ElType>
  void choleskyInverse(UpLo uplo, DistributedMatrix<ElType>& mat, SolverType solver,
                       int print_timers = 0);

  template <class ElType>
  void LUFactorization(DistributedMatrix<ElType>& mat, std::vector<int>& ipiv, SolverType solver,
                       int print_timers = 0);

  template <class ElType>
  void LUFactorization(DistributedMatrix<ElType>& mat, int* ipiv, SolverType solver,
                       int print_timers = 0);

  template <class ElType>
  void hermitianEigenvectors(UpLo uplo, DistributedMatrix<ElType>& mat,
                             std::vector<BaseType<ElType>>& evalues,
                             DistributedMatrix<ElType>& evectors, SolverType solver,
                             int print_timers = 0);

  template <class ElType>
  void hermitianEigenvectors(UpLo uplo, DistributedMatrix<ElType>& mat, BaseType<ElType>* evalues,
                             DistributedMatrix<ElType>& evectors, SolverType solver,
                             int print_timers = 0);

  template <class ElType>
  void triangularInverse(UpLo uplo, Diag diag, DistributedMatrix<ElType>& mat, SolverType solver,
                         int print_timers = 0);

#include "dla_interface/cholesky_factorization.ipp"
#include "dla_interface/cholesky_inverse.ipp"
#include "dla_interface/lu_factorization.ipp"
#include "dla_interface/matrix_multiplication.ipp"
#include "dla_interface/hermitian_eigenvectors.ipp"
#include "dla_interface/triangular_inverse.ipp"
}

#endif  // DLA_INTERFACE_DLA_INTERFACE_H
