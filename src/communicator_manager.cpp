#include <map>
#include <memory>
#include <mpi.h>
#include <stdexcept>
#include <string>
#include "types.h"
#include "communicator_grid.h"
#include "communicator_manager.h"
#include "internal_error.h"

namespace dla_interface {
  namespace comm {
    // static members:
    CommunicatorManager::Status CommunicatorManager::status_(CommunicatorManager::non_initialized);
    std::unique_ptr<CommunicatorManager> CommunicatorManager::comm_manager_(nullptr);

    void CommunicatorManager::initialize(bool initialize_mpi) {
      initialize(-1, initialize_mpi);
    }
    void CommunicatorManager::initialize(int nr_cores, bool initialize_mpi) {
      initialize(nr_cores, nullptr, nullptr, initialize_mpi);
    }
    void CommunicatorManager::initialize(int nr_cores, int* argc, char*** argv, bool initialize_mpi) {
      comm_manager_.reset(new CommunicatorManager(nr_cores, argc, argv, initialize_mpi));
      status_ = initialized;
    }

    void CommunicatorManager::finalize() {
      status_ = finalized;
      comm_manager_ = nullptr;
    }

#ifdef DLA_HAVE_DPLASMA
    ParsecContext CommunicatorManager::getParsecContext() {
      return comm_manager_->parsec_handle_;
    }
#endif

    Communicator2DGrid& CommunicatorManager::createCommunicator2DGrid(  //
        MPI_Comm base_comm, int nr_rows, int nr_cols, Ordering comm_ordering) {
      return comm_manager_->communicator2DGrid(base_comm, nr_rows, nr_cols, comm_ordering);
    }

#ifdef DLA_HAVE_SCALAPACK
    Communicator2DGrid& CommunicatorManager::getCommunicator2DGridFromBlacsContext(BlacsContextType id) {
      return comm_manager_->communicator2DGridFromBlacsContext(id);
    }
#endif

    Communicator2DGrid& CommunicatorManager::getCommunicator2DGridFromMPIComm(MPI_Comm comm) {
      return comm_manager_->communicator2DGridFromMPIComm(comm);
    }

    // non-static members:
    Communicator2DGrid& CommunicatorManager::communicator2DGrid(MPI_Comm base_comm, int row_size,
                                                                int col_size, Ordering comm_ordering) {
      std::shared_ptr<Communicator2DGrid> comm_grid_(
          new Communicator2DGrid(base_comm, row_size, col_size, comm_ordering));

      if (!comm_grid_map.insert(std::make_pair(comm_grid_->rowMPICommunicator(), comm_grid_)).second)
        throw error::InternalError("Cannot insert row commmunicator in the map");
      if (!comm_grid_map.insert(std::make_pair(comm_grid_->colMPICommunicator(), comm_grid_)).second)
        throw error::InternalError("Cannot insert column commmunicator in the map");
      if (!comm_grid_map.insert(std::make_pair(comm_grid_->rowOrderedMPICommunicator(), comm_grid_)).second)
        throw error::InternalError("Cannot insert row ordered commmunicator in the map");
#ifdef DLA_HAVE_SCALAPACK
      if (!ictxt_grid_map.insert(std::make_pair(comm_grid_->blacsContext(), comm_grid_)).second)
        throw error::InternalError("Cannot insert BLACS context id in the map");
#endif
      return *comm_grid_;
    }
    Communicator2DGrid& CommunicatorManager::communicator2DGridFromMPIComm(MPI_Comm comm) const {
      try {
        return *comm_grid_map.at(comm);
      }
      catch (std::out_of_range) {
        throw std::invalid_argument("No communicator2DGrid found with the given MPI_Comm");
      }
    }
#ifdef DLA_HAVE_SCALAPACK
    Communicator2DGrid& CommunicatorManager::communicator2DGridFromBlacsContext(BlacsContextType id) const {
      try {
        return *ictxt_grid_map.at(id);
      }
      catch (std::out_of_range) {
        throw std::invalid_argument("No communicator2DGrid found with the given BLACS context id");
      }
    }
#endif

    CommunicatorManager::CommunicatorManager(int nr_cores, int* argc, char*** argv,
                                             bool initialize_mpi)
     : init_mpi_(initialize_mpi) {
      if (init_mpi_) {
        int provided;
        MPI_Init_thread(argc, argv, MPI_THREAD_SERIALIZED, &provided);
        if (provided != MPI_THREAD_SERIALIZED)
          throw std::runtime_error(
              "MPI cannot be initializd with MPI_THREAD_SERIALIZED. (provided = " +
              std::to_string(provided) + ")");
      }
#ifdef DLA_HAVE_DPLASMA
      parsec_handle_ = parsec_init(nr_cores, argc, argv);
#endif
    }

    CommunicatorManager::~CommunicatorManager() {
      // Make sure that the grids are destroyed before calling MPI_Finalize
      comm_grid_map.clear();
#ifdef DLA_HAVE_SCALAPACK
      ictxt_grid_map.clear();
#endif
#ifdef DLA_HAVE_DPLASMA
      parsec_fini(&parsec_handle_);
#endif

      if (init_mpi_) {
        MPI_Finalize();
      }
    }
  }
}
