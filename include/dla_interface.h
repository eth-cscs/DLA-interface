#ifndef DLA_INTERFACE_DLA_INTERFACE_H
#define DLA_INTERFACE_DLA_INTERFACE_H

#include <stdexcept>
#include "distributed_matrix.h"
#include "dla_scalapack.h"
#include "error_message.h"

namespace dla_interface {

  template <class ElType>
  void choleskyFactorization(UpLo uplo, DistributedMatrix<ElType>& mat, SolverType solver) {
    // check size (square), blocksize?
    switch (solver) {
#ifdef DLA_HAVE_SCALAPACK
      case ScaLAPACK: {
        DistributedMatrix<ElType> mat_scalapack(scalapack_dist, mat);
        auto info = mat_scalapack.getScalapackDescription();
        scalapack_wrappers::ppotrf(uplo, mat_scalapack.size().first, std::get<0>(info),
                                   std::get<1>(info), std::get<2>(info), &std::get<3>(info)[0]);
        break;
      }
#endif
      default:
        throw std::invalid_argument(
            errorMessage("Cholesky factorization is not available for solver ", solver));
    }
  }
}

#endif  // DLA_INTERFACE_DLA_INTERFACE_H
