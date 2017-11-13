#ifndef DLA_INTERFACE_DLA_DPLASMA_H
#define DLA_INTERFACE_DLA_DPLASMA_H

#ifdef DLA_HAVE_DPLASMA

#include <mpi.h>
#include "communicator_manager.h"
#include "ordered_dplasma.h"

namespace dla_interface {
  namespace dplasma_wrappers {

    inline PLASMA_enum plasmaUpLo(UpLo uplo) {
      switch (uplo) {
        case Lower:
          return PlasmaLower;
        case Upper:
          return PlasmaUpper;
      }
    }

    template <class Routine, class... Args>
    inline void dplasma_run(MPI_Comm comm, Args... args) {
      auto parsec = comm::CommunicatorManager::getParsecContext();
      // Sets the MPI communicator.
      // TODO: uncomment this after parsec issue #135 is fixed
      //       and remove continue in test_dla_interface.cpp:66.
      // parsec_remote_dep_set_ctx(parsec, &comm);
      // Creates and schedules the operation
      parsec_taskpool_t* taskpool = Routine::create(args...);
      parsec_enqueue(parsec, taskpool);
      // Parsec starts the computation
      parsec_context_start(parsec);
      // Wait for the completion
      parsec_context_wait(parsec);
    }

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
