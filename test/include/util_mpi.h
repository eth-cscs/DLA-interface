#ifndef DLA_INTERFACE_TEST_INCLUDE_UTIL_MPI_H
#define DLA_INTERFACE_TEST_INCLUDE_UTIL_MPI_H

#include <mpi.h>
#include <utility>

namespace testing {

  inline std::pair<int, int> commInfo(MPI_Comm comm = MPI_COMM_WORLD) {
    int id = -1;
    int size = -1;
    MPI_Comm_rank(comm, &id);
    MPI_Comm_size(comm, &size);
    return std::make_pair(id, size);
  }
}

#endif  // DLA_INTERFACE_TEST_INCLUDE_UTIL_MPI_H
