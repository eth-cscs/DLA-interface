#ifndef DLA_INTERFACE_DLA_INTERFACE_H
#define DLA_INTERFACE_DLA_INTERFACE_H

#include <stdexcept>
#include "distributed_matrix.h"
#include "dla_scalapack.h"
#include "dla_dplasma.h"
#include "error_message.h"

namespace dla_interface {

  template <class ElType>
  void choleskyFactorization(UpLo uplo, DistributedMatrix<ElType>& mat, SolverType solver) {
    // check size (square), blocksize?
    switch (solver) {
#ifdef DLA_HAVE_SCALAPACK
      case ScaLAPACK: {
        int info = 0;
        DistributedMatrix<ElType> mat_scalapack(scalapack_dist, mat);
        auto matrix_info = mat_scalapack.getScalapackDescription();

        scalapack_wrappers::ppotrf(uplo, mat_scalapack.size().first, std::get<0>(matrix_info),
                                   std::get<1>(matrix_info), std::get<2>(matrix_info),
                                   &std::get<3>(matrix_info)[0], info);

        if (info != 0) {
          if (info < 0)
            throw std::invalid_argument(errorMessage("Argument ", -info, " is wrong."));
          else
            throw std::invalid_argument(
                errorMessage("Matrix is not positive definite (", info, ")"));
        }
        break;
      }
#endif

#ifdef DLA_HAVE_DPLASMA
      case DPlasma: {
        int info = 0;
        DistributedMatrix<ElType> mat_tile(tile_dist, mat);
        auto matrix_info = mat_tile.getDPlasmaDescription();
        parsec_tiled_matrix_dc_t* dp_mat =
            reinterpret_cast<parsec_tiled_matrix_dc_t*>(&std::get<0>(matrix_info));

        dplasma_wrappers::dplasma_run<dplasma_wrappers::ppotrf<ElType>>(
            std::get<1>(matrix_info), dplasma_wrappers::plasmaUpLo(uplo), dp_mat, &info);

        if (info != 0) {
          throw std::invalid_argument(errorMessage("Matrix is not positive definite (", info, ")"));
        }
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
