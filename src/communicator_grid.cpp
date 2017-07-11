#include <array>
#include <mpi.h>
#include <stdexcept>
#include <string>
#include <utility>
#include "blacs.h"
#include "types.h"
#include "communicator_grid.h"

namespace dla_interface {
  namespace comm {

    Communicator2DGrid::Communicator2DGrid(MPI_Comm base_comm, int nr_rows, int nr_cols,
                                           Ordering comm_ordering)
     : id2D_(-1, -1), size2D_(nr_rows, nr_cols), order_(comm_ordering), row_comm_(MPI_COMM_NULL),
       col_comm_(MPI_COMM_NULL), row_ordered_comm_(MPI_COMM_NULL) {
      int rank = -1;
      MPI_Comm_rank(base_comm, &rank);
      int size = -1;
      MPI_Comm_size(base_comm, &size);

      if (size != nr_rows * nr_cols) {
        std::string msg("The size of the communicator is wrong:");
        msg += "comm_size (" + std::to_string(size) + ") != ";
        msg += "nr_row (" + std::to_string(nr_rows) + ") * nr_col (" + std::to_string(nr_cols) + ")";

        throw std::invalid_argument(msg);
      }
      if (order_ == ColMajor) {
        id2D_ = std::make_pair(rank % nr_rows, rank / nr_rows);
      }
      else {
        id2D_ = std::make_pair(rank / nr_cols, rank % nr_cols);
      }

      MPI_Comm tmp_comm;
      if (order_ == ColMajor) {
        // Create a temporary communicator with RowMajor ordered ranks.
        int row_ordered_rank = id2D_.first * nr_cols + id2D_.second;
        MPI_Comm_split(base_comm, 1, row_ordered_rank, &tmp_comm);
      }
      else {
        // We already have a RowMajor ordered communicator.
        tmp_comm = base_comm;
      }

      std::array<int, 2> dims({nr_rows, nr_cols});
      std::array<int, 2> periods({1, 1});

      MPI_Cart_create(tmp_comm, 2, &dims[0], &periods[0], 0, &row_ordered_comm_);

      if (order_ == ColMajor) {
        // Free the temporary communicator created before.
        MPI_Comm_free(&tmp_comm);
      }

      std::array<int, 2> free_coords;
      free_coords = {0, 1};
      MPI_Cart_sub(row_ordered_comm_, &free_coords[0], &row_comm_);
      free_coords = {1, 0};
      MPI_Cart_sub(row_ordered_comm_, &free_coords[0], &col_comm_);

#ifdef DLA_HAVE_SCALAPACK
      blacs_ictxt_ = blacs::Csys2blacs_handle(base_comm);
      char ordering = 'R';
      if (comm_ordering == ColMajor)
        ordering = 'C';
      blacs::Cblacs_gridinit(&blacs_ictxt_, &ordering, nr_rows, nr_cols);
#endif
    }

    Communicator2DGrid::~Communicator2DGrid() {
#ifdef DLA_HAVE_SCALAPACK
      blacs::Cblacs_gridexit(blacs_ictxt_);
#endif
      MPI_Comm_free(&col_comm_);
      MPI_Comm_free(&row_comm_);
      MPI_Comm_free(&row_ordered_comm_);
    }
  }
}
