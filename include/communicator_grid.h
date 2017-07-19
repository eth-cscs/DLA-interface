#ifndef DLA_INTERFACE_COMMUNICATOR_GRID_H
#define DLA_INTERFACE_COMMUNICATOR_GRID_H

#include <mpi.h>
#include <utility>
#include "types.h"

namespace dla_interface {
  namespace comm {

    class Communicator2DGrid {
      public:
      // Initializes the 2D communication grid:
      // - it initializes the MPI communicators,
      // - if ScaLAPACK is available it initializes the BLACS context
      // Throws a std::invalid_argument exception
      // if the number of ranks in base_comm is different from nr_rows * nr_cols.
      Communicator2DGrid(MPI_Comm base_comm, int nr_rows, int nr_cols, Ordering comm_ordering);

      Communicator2DGrid(const Communicator2DGrid& rhs) = delete;
      Communicator2DGrid operator=(const Communicator2DGrid& rhs) = delete;

      ~Communicator2DGrid();

      // Returns the row communicator (the communicator size is nr_cols).
      MPI_Comm rowMPICommunicator() const;
      // Returns the column communicator (the communicator size is nr_rows).
      MPI_Comm colMPICommunicator() const;

      // Returns the MPI communicator with the id of the ranks of the 2D grid in row major order.
      MPI_Comm rowOrderedMPICommunicator() const;
// TODO: Needed ?
// MPI_Comm MPI_ColOrderedCommunicator() const;

#ifdef DLA_HAVE_SCALAPACK
      BlacsContextType blacsContext() const;
#endif

      // Returns the 2D id of this rank.
      std::pair<int, int> id2D() const;
      // Returns the 2D size of this grid.
      std::pair<int, int> size2D() const;

      // Returns the rank order of this grid.
      Ordering rankOrder() const;
    };
  }
}

#endif  // DLA_INTERFACE_COMMUNICATOR_GRID_H
