#ifndef DLA_INTERFACE_COMMUNICATOR_MANAGER_H
#define DLA_INTERFACE_COMMUNICATOR_MANAGER_H

#include <map>
#include <memory>
#include <mpi.h>
#include "communicator_grid.h"
#include "types.h"

namespace dla_interface {
  namespace comm {

    class CommunicatorManager {
      public:
      // Initializes the CommunicatorManager.
      // If DLA_HAVE_DPLASMA is true parsec is initialized with nr_cores (default -1).
      // If initialize_mpi is true:
      // - MPI_Init_thread is called using MPI_THREAD_SERIALIZED,
      // - MPI_Finalize is called when finalize is called.
      // Precondition: initialize must have not been called before. (It can be called only once.)
      // If specified argc and argv are used in MPI and DPLASMA (if DLA_HAVE_DPLASMA is true)
      // initialization (arguments after -- are passed to Parsec).
      static void initialize(bool initialize_mpi = true);
      static void initialize(int nr_cores, bool initialize_mpi = true);
      static void initialize(int nr_cores, int* argc, char*** argv, bool initialize_mpi = true);
      // Destructs all the 2Dgrids created with createCommunicator2DGrid() and
      // if initialize was called with initialize_mpi == true this function finalize MPI as well.
      // Precondition: initialize must have been called before and
      //               finalize must have not been called before. (It can be called only once.)
      static void finalize();
#ifdef DLA_HAVE_DPLASMA
      // Returns the parsec context.
      static ParsecContext getParsecContext();
#endif

      // Creates a 2D grid based on the given MPI communicator, with nr_rows rows
      // and nr_cols columns, where the ranks are ordered in comm_ordering order.
      // The created Communicator2DGrid object is stored and a reference to it is returned.
      // Throws a std::invalid_argument exception
      // if the number of ranks in base_comm is different from nr_rows * nr_cols.
      // Precondition: initialize must have been called before and
      //               finalize must have not been called before.
      static Communicator2DGrid& createCommunicator2DGrid(MPI_Comm base_comm, int nr_rows,
                                                          int nr_cols, Ordering comm_ordering);

#ifdef DLA_HAVE_SCALAPACK
      // Returns a reference of the Communicator2DGrid created with createCommunicator2DGrid
      // whose BLACS context index is id.
      // Throws a std::invalid_argument exception if no Communicator2DGrid is found.
      static Communicator2DGrid& getCommunicator2DGridFromBlacsContext(BlacsContextType id);
#endif

      // Returns a reference of the Communicator2DGrid created with createCommunicator2DGrid
      // whose row, col or row_ordered (TODO: or col_ordered) communicator is comm.
      // Throws a std::invalid_argument exception if no Communicator2DGrid is found.
      static Communicator2DGrid& getCommunicator2DGridFromMPIComm(MPI_Comm comm);

      protected:
      Communicator2DGrid& communicator2DGrid(MPI_Comm base_comm, int row_size, int col_size,
                                             Ordering comm_ordering);
      Communicator2DGrid& communicator2DGridFromMPIComm(MPI_Comm comm) const;
#ifdef DLA_HAVE_SCALAPACK
      Communicator2DGrid& communicator2DGridFromBlacsContext(BlacsContextType id) const;
#endif

      private:
      CommunicatorManager(int nr_cores = -1, int* argc = nullptr, char*** argv = nullptr,
                          bool initialize_mpi = true);

      CommunicatorManager(const CommunicatorManager&) = delete;
      CommunicatorManager operator=(const CommunicatorManager&) = delete;

      public:
      ~CommunicatorManager();

      private:
      enum Status { non_initialized = 0, initialized = 1, finalized = 2 };
      static Status status_;
      static std::unique_ptr<CommunicatorManager> comm_manager_;

      bool init_mpi_;
#ifdef DLA_HAVE_DPLASMA
      ParsecContext parsec_handle_;
#endif

      std::map<MPI_Comm, std::shared_ptr<Communicator2DGrid>> comm_grid_map;
#ifdef DLA_HAVE_SCALAPACK
      std::map<BlacsContextType, std::shared_ptr<Communicator2DGrid>> ictxt_grid_map;
#endif
    };
  }
}

#endif  // DLA_INTERFACE_COMMUNICATOR_MANAGER_H
