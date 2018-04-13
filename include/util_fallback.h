#ifndef DLA_INTERFACE_UTIL_FALLBACK_H
#define DLA_INTERFACE_UTIL_FALLBACK_H

#include <mpi.h>
#include "communicator_grid.h"
#include "communicator_manager.h"
#include "types.h"

namespace dla_interface {
  namespace util {

    inline SolverType fallbackCommunicator(const char* func, const comm::Communicator2DGrid& comm,
                                           SolverType solver) {
#ifdef DLA_HAVE_DPLASMA
      if (solver == DPlasma) {
        int res;
        MPI_Comm_compare(comm.rowOrderedMPICommunicator(), MPI_COMM_WORLD, &res);

        if (res == MPI_UNEQUAL || res == MPI_SIMILAR) {
          // only rank (0, 0) of comm report the fallback;
          if (comm.id2D() == std::make_pair(0, 0))
            comm::CommunicatorManager::getFallbackInfo().report(func, solver, ScaLAPACK,
                                                                "Non compatible Communicator.");
          return ScaLAPACK;
        }
      }
#endif
      return solver;
    }
  }
}

#define dlai__util__fallbackCommunicator(comm, solver) \
  dla_interface::util::fallbackCommunicator(__func__, comm, solver)

#endif  // DLA_INTERFACE_UTIL_FALLBACK_H
