#ifndef DLA_INTERFACE_DLA_DPLASMA_H
#define DLA_INTERFACE_DLA_DPLASMA_H

#ifdef DLA_HAVE_DPLASMA

#include <mpi.h>
#include "communicator_manager.h"
#include "ordered_dplasma.h"
#include "util_thread.h"

namespace dla_interface {
  namespace dplasma_wrappers {

    inline PLASMA_enum plasmaUpLo(UpLo uplo) {
      switch (uplo) {
        case Lower:
          return PlasmaLower;
        case Upper:
          return PlasmaUpper;
        default:
          throw(std::invalid_argument(errorMessage("Invalid UpLo element.", uplo)));
      }
    }

    inline PLASMA_enum plasmaTrans(OpTrans trans) {
      switch (trans) {
        case NoTrans:
          return PlasmaNoTrans;
        case Trans:
          return PlasmaTrans;
        case ConjTrans:
          return PlasmaConjTrans;
        default:
          throw(std::invalid_argument(errorMessage("Invalid OpTrans element.", trans)));
      }
    }

    inline PLASMA_enum plasmaDiag(Diag diag) {
      switch (diag) {
        case Unit:
          return PlasmaUnit;
        case NonUnit:
          return PlasmaNonUnit;
        default:
          throw(std::invalid_argument(errorMessage("Invalid Diag element.", diag)));
      }
    }

    template <class Routine, class... Args>
    inline int dplasma_run(MPI_Comm comm, Args... args) {
      auto parsec = comm::CommunicatorManager::getParsecContext();
      // Sets the MPI communicator.
      // TODO: uncomment this after parsec issue #135 is fixed
      //       and remove fallback to ScaLAPACK.
      // parsec_remote_dep_set_ctx(parsec, &comm);
      // Abort if the communicator is not equivalent to MPI_COMM_WORLD
      int size1, size2, rank1, rank2;
      MPI_Comm_size(MPI_COMM_WORLD, &size1);
      MPI_Comm_rank(MPI_COMM_WORLD, &rank1);
      MPI_Comm_size(comm, &size2);
      MPI_Comm_rank(comm, &rank2);
      if (size1 != size2 || rank1 != rank2) {
        std::cerr << "Unsupported communicator for DPlasma" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 10000);
      }

      // Set number of BLAS threads and CPU binding. (Reset on object destruction)
      util::SetNumThreadsAndCpuBind config(comm::CommunicatorManager::getDPlasmaConfigInfo());

      // Creates and schedules the operation
      return Routine::run(parsec, args...);
    }

    template <class ElType>
    struct pgemm {};

    template <>
    struct pgemm<float> {
      template <class... Args>
      static auto run(Args... args) {
        return dplasma_sgemm(args...);
      }
    };

    template <>
    struct pgemm<double> {
      template <class... Args>
      static auto run(Args... args) {
        return dplasma_dgemm(args...);
      }
    };

    template <>
    struct pgemm<std::complex<float>> {
      template <class... Args>
      static auto run(Args... args) {
        return dplasma_cgemm(args...);
      }
    };

    template <>
    struct pgemm<std::complex<double>> {
      template <class... Args>
      static auto run(Args... args) {
        return dplasma_zgemm(args...);
      }
    };

    template <class ElType>
    struct pgetrf {};

    template <>
    struct pgetrf<float> {
      template <class... Args>
      static auto run(Args... args) {
        return dplasma_sgetrf(args...);
      }
    };

    template <>
    struct pgetrf<double> {
      template <class... Args>
      static auto run(Args... args) {
        return dplasma_dgetrf(args...);
      }
    };

    template <>
    struct pgetrf<std::complex<float>> {
      template <class... Args>
      static auto run(Args... args) {
        return dplasma_cgetrf(args...);
      }
    };

    template <>
    struct pgetrf<std::complex<double>> {
      template <class... Args>
      static auto run(Args... args) {
        return dplasma_zgetrf(args...);
      }
    };

    template <class ElType>
    struct ppotrf {};

    template <>
    struct ppotrf<float> {
      template <class... Args>
      static auto run(Args... args) {
        return dplasma_spotrf(args...);
      }
    };

    template <>
    struct ppotrf<double> {
      template <class... Args>
      static auto run(Args... args) {
        return dplasma_dpotrf(args...);
      }
    };

    template <>
    struct ppotrf<std::complex<float>> {
      template <class... Args>
      static auto run(Args... args) {
        return dplasma_cpotrf(args...);
      }
    };

    template <>
    struct ppotrf<std::complex<double>> {
      template <class... Args>
      static auto run(Args... args) {
        return dplasma_zpotrf(args...);
      }
    };

    template <class ElType>
    struct ppotri {};

    template <>
    struct ppotri<float> {
      template <class... Args>
      static auto run(Args... args) {
        return dplasma_spotri(args...);
      }
    };

    template <>
    struct ppotri<double> {
      template <class... Args>
      static auto run(Args... args) {
        return dplasma_dpotri(args...);
      }
    };

    template <>
    struct ppotri<std::complex<float>> {
      template <class... Args>
      static auto run(Args... args) {
        return dplasma_cpotri(args...);
      }
    };

    template <>
    struct ppotri<std::complex<double>> {
      template <class... Args>
      static auto run(Args... args) {
        return dplasma_zpotri(args...);
      }
    };

    template <class ElType>
    struct ptrtri {};

    template <>
    struct ptrtri<float> {
      template <class... Args>
      static auto run(Args... args) {
        return dplasma_strtri(args...);
      }
    };

    template <>
    struct ptrtri<double> {
      template <class... Args>
      static auto run(Args... args) {
        return dplasma_dtrtri(args...);
      }
    };

    template <>
    struct ptrtri<std::complex<float>> {
      template <class... Args>
      static auto run(Args... args) {
        return dplasma_ctrtri(args...);
      }
    };

    template <>
    struct ptrtri<std::complex<double>> {
      template <class... Args>
      static auto run(Args... args) {
        return dplasma_ztrtri(args...);
      }
    };
  }
}

#endif

#endif  // DLA_INTERFACE_DLA_DPLASMA_H
