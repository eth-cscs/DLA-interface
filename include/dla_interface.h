#ifndef DLA_INTERFACE_DLA_INTERFACE_H
#define DLA_INTERFACE_DLA_INTERFACE_H

#include <stdexcept>
#include "distributed_matrix.h"
#include "dla_scalapack.h"
#include "dla_dplasma.h"
#include "error_message.h"
#include "util_cast.h"

namespace dla_interface {

  // TODO: Major change to matrix interface for change of distribution for const matrices.
  template <class ElType>
  void matrixMultiply(OpTrans trans_a, OpTrans trans_b, ElType alpha,
                      /* const */ DistributedMatrix<ElType>& mat_a,
                      /* const */ DistributedMatrix<ElType>& mat_b, ElType beta,
                      DistributedMatrix<ElType>& mat_c, SolverType solver) {
    // TODO: check blocksize, communicator?
    int m = mat_c.size().first;
    int n = mat_c.size().second;
    int m2 = trans_a == NoTrans ? mat_a.size().first : mat_a.size().second;
    int k = trans_a == NoTrans ? mat_a.size().second : mat_a.size().first;
    int k2 = trans_b == NoTrans ? mat_b.size().first : mat_b.size().second;
    int n2 = trans_b == NoTrans ? mat_b.size().second : mat_b.size().first;
    if (m != m2 || n != n2 || k != k2)
      throw std::invalid_argument(errorMessage("Incompatible sizes in matrix multiplication: ", m,
                                               " != ", m2, " || ", n, " != ", n2, " || ", k, " != ",
                                               k2));

    switch (solver) {
#ifdef DLA_HAVE_SCALAPACK
      case ScaLAPACK: {
        DistributedMatrix<ElType> mat_a_scalapack(scalapack_dist, mat_a);
        DistributedMatrix<ElType> mat_b_scalapack(scalapack_dist, mat_b);
        DistributedMatrix<ElType> mat_c_scalapack(scalapack_dist, mat_c);
        auto matrix_a_info = mat_a_scalapack.getScalapackDescription();
        auto matrix_b_info = mat_b_scalapack.getScalapackDescription();
        auto matrix_c_info = mat_c_scalapack.getScalapackDescription();

        scalapack_wrappers::pgemm(
            trans_a, trans_b, m, n, k, alpha, std::get<0>(matrix_a_info),
            std::get<1>(matrix_a_info), std::get<2>(matrix_a_info), &std::get<3>(matrix_a_info)[0],
            std::get<0>(matrix_b_info), std::get<1>(matrix_b_info), std::get<2>(matrix_b_info),
            &std::get<3>(matrix_b_info)[0], beta, std::get<0>(matrix_c_info),
            std::get<1>(matrix_c_info), std::get<2>(matrix_c_info), &std::get<3>(matrix_c_info)[0]);

        break;
      }
#endif

#ifdef DLA_HAVE_DPLASMA
      case DPlasma: {
        DistributedMatrix<ElType> mat_a_tile(tile_dist, mat_a);
        DistributedMatrix<ElType> mat_b_tile(tile_dist, mat_b);
        DistributedMatrix<ElType> mat_c_tile(tile_dist, mat_c);
        auto matrix_a_info = mat_a_tile.getDPlasmaDescription();
        auto matrix_b_info = mat_b_tile.getDPlasmaDescription();
        auto matrix_c_info = mat_c_tile.getDPlasmaDescription();
        const parsec_tiled_matrix_dc_t* dp_mat_a =
            reinterpret_cast<const parsec_tiled_matrix_dc_t*>(&std::get<0>(matrix_a_info));
        const parsec_tiled_matrix_dc_t* dp_mat_b =
            reinterpret_cast<const parsec_tiled_matrix_dc_t*>(&std::get<0>(matrix_b_info));
        parsec_tiled_matrix_dc_t* dp_mat_c =
            reinterpret_cast<parsec_tiled_matrix_dc_t*>(&std::get<0>(matrix_c_info));

        dplasma_wrappers::dplasma_run<dplasma_wrappers::pgemm<ElType>>(
            std::get<1>(matrix_c_info), dplasma_wrappers::plasmaTrans(trans_a),
            dplasma_wrappers::plasmaTrans(trans_b), util::castToC(alpha), dp_mat_a, dp_mat_b,
            util::castToC(beta), dp_mat_c);

        break;
      }
#endif

      default:
        throw std::invalid_argument(
            errorMessage("Matrix multiplication is not available for solver ", solver));
    }
  }
  template <class ElType>
  void choleskyFactorization(UpLo uplo, DistributedMatrix<ElType>& mat, SolverType solver) {
    // TODO: check size (square), blocksize?
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
