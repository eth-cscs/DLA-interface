//
// Distributed Linear Algebra Interface (DLAI)
//
// Copyright (c) 2018-2019, ETH Zurich
// All rights reserved.
//
// Please, refer to the LICENSE file in the root directory.
// SPDX-License-Identifier: BSD-3-Clause
//

#ifndef DLA_INTERFACE_COMMUNICATOR_GRID_H
#define DLA_INTERFACE_COMMUNICATOR_GRID_H

#include <mpi.h>
#include <utility>
#include "types.h"

namespace dla_interface {
  namespace comm {

  /** Communicator2DGrid class.
   *
   * A 2D Grid communicator class for MPI layer.
   *
   */

    class Communicator2DGrid {
      public:

    	/** Initializes the 2D communication grid.
    	 *
    	 * Initializes the 2D communication grid:
         * - it initializes the MPI communicators,
         * - if ScaLAPACK is available it initializes the BLACS context
         *
         * @throws std::invalid_argument if the number of ranks in base_comm is different from nr_rows * nr_cols.
         *
    	 * @param base_comm MPI base communicator.
    	 * @param nr_rows Number of rows.
    	 * @param nr_cols Number of columns.
    	 * @param comm_ordering Communication order.
    	 */
      Communicator2DGrid(MPI_Comm base_comm, int nr_rows, int nr_cols, Ordering comm_ordering);

      /// Copy constructor is prohibited.
      Communicator2DGrid(const Communicator2DGrid& rhs) = delete;
      /// Copy function is prohibited;
      Communicator2DGrid operator=(const Communicator2DGrid& rhs) = delete;

      /// Communicator2DGrid destructor.
      ~Communicator2DGrid();

      /// Returns the row communicator (the communicator size is nr_cols).
      MPI_Comm rowMPICommunicator() const {
        return row_comm_;
      }

      /// Returns the column communicator (the communicator size is nr_rows).
      MPI_Comm colMPICommunicator() const {
        return col_comm_;
      }

      /// Returns the MPI communicator with the id of the ranks of the 2D grid in row major order.
      MPI_Comm rowOrderedMPICommunicator() const {
        return row_ordered_comm_;
      }

/// @TODO: Needed ? MPI_Comm MPI_ColOrderedCommunicator() const;

#ifdef DLA_HAVE_SCALAPACK
      BlacsContextType blacsContext() const {
        return blacs_ictxt_;
      }
#endif

      /// Returns the 2D id of this rank.
      std::pair<int, int> id2D() const {
        return id2D_;
      }

      /// Returns the 2D size of this grid.
      std::pair<int, int> size2D() const {
        return size2D_;
      }

      /// Returns the the total size of this grid.
      int size() const {
        return size2D_.first * size2D_.second;
      }

      /// Returns the rank order of this grid.
      Ordering rankOrder() const {
        return order_;
      }

      private:
      std::pair<int, int> id2D_;
      std::pair<int, int> size2D_;
      Ordering order_;
      MPI_Comm row_comm_;
      MPI_Comm col_comm_;
      MPI_Comm row_ordered_comm_;

#ifdef DLA_HAVE_SCALAPACK
      BlacsContextType blacs_ictxt_;
#endif
    };
  }
}

#endif  // DLA_INTERFACE_COMMUNICATOR_GRID_H
