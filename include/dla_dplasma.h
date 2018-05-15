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
      }
      throw(std::invalid_argument(errorMessage("Invalid UpLo element.", uplo)));
    }

    inline PLASMA_enum plasmaTrans(OpTrans trans) {
      switch (trans) {
        case NoTrans:
          return PlasmaNoTrans;
        case Trans:
          return PlasmaTrans;
        case ConjTrans:
          return PlasmaConjTrans;
      }
      throw(std::invalid_argument(errorMessage("Invalid OpTrans element.", trans)));
    }

    template <class Routine, class... Args>
    inline void dplasma_run(MPI_Comm comm, Args... args) {
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
      parsec_taskpool_t* taskpool = Routine::create(args...);
      parsec_enqueue(parsec, taskpool);
      // Parsec starts the computation
      parsec_context_start(parsec);
      // Wait for the completion
      parsec_context_wait(parsec);
    }

    template <class ElType>
    struct pgemm {};

    template <>
    struct pgemm<float> {
      template <class... Args>
      static auto create(Args... args) {
        return dplasma_sgemm_New(args...);
      }
    };

    template <>
    struct pgemm<double> {
      template <class... Args>
      static auto create(Args... args) {
        return dplasma_dgemm_New(args...);
      }
    };

    template <>
    struct pgemm<std::complex<float>> {
      template <class... Args>
      static auto create(Args... args) {
        return dplasma_cgemm_New(args...);
      }
    };

    template <>
    struct pgemm<std::complex<double>> {
      template <class... Args>
      static auto create(Args... args) {
        return dplasma_zgemm_New(args...);
      }
    };

    template <class ElType>
    struct pgetrf {};

    template <>
    struct pgetrf<float> {
      template <class... Args>
      static auto create(Args... args) {
        return dplasma_sgetrf_New(args...);
      }
    };

    template <>
    struct pgetrf<double> {
      template <class... Args>
      static auto create(Args... args) {
        return dplasma_dgetrf_New(args...);
      }
    };

    template <>
    struct pgetrf<std::complex<float>> {
      template <class... Args>
      static auto create(Args... args) {
        return dplasma_cgetrf_New(args...);
      }
    };

    template <>
    struct pgetrf<std::complex<double>> {
      template <class... Args>
      static auto create(Args... args) {
        return dplasma_zgetrf_New(args...);
      }
    };

    template <class ElType>
    struct ppotrf {};

    template <>
    struct ppotrf<float> {
      template <class... Args>
      static auto create(Args... args) {
        return dplasma_spotrf_New(args...);
      }
    };

    template <>
    struct ppotrf<double> {
      template <class... Args>
      static auto create(Args... args) {
        return dplasma_dpotrf_New(args...);
      }
    };

    template <>
    struct ppotrf<std::complex<float>> {
      template <class... Args>
      static auto create(Args... args) {
        return dplasma_cpotrf_New(args...);
      }
    };

    template <>
    struct ppotrf<std::complex<double>> {
      template <class... Args>
      static auto create(Args... args) {
        return dplasma_zpotrf_New(args...);
      }
    };
  }
}

#endif

#endif  // DLA_INTERFACE_DLA_DPLASMA_H
