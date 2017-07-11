#include <map>
#include <memory>
#include <mpi.h>
#include <stdexcept>
#include "types.h"
#include "communicator_grid.h"
#include "communicator_manager.h"

namespace dla_interface {
  namespace comm {
    // static members:
    CommunicatorManager::Status CommunicatorManager::status_(CommunicatorManager::non_initialized);
    std::unique_ptr<CommunicatorManager> CommunicatorManager::comm_manager_(nullptr);

    void CommunicatorManager::initialize(bool initialize_mpi) {
      comm_manager_.reset(new CommunicatorManager(initialize_mpi));
      status_ = initialized;
    }

    void CommunicatorManager::finalize() {
      status_ = finalized;
      comm_manager_ = nullptr;
    }

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
        ;  // TODO: abort
      if (!comm_grid_map.insert(std::make_pair(comm_grid_->colMPICommunicator(), comm_grid_)).second)
        ;  // TODO: abort
      if (!comm_grid_map.insert(std::make_pair(comm_grid_->rowOrderedMPICommunicator(), comm_grid_)).second)
        ;  // TODO: abort
#ifdef DLA_HAVE_SCALAPACK
      if (!ictxt_grid_map.insert(std::make_pair(comm_grid_->blacsContext(), comm_grid_)).second)
        ;  // TODO: abort
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

    CommunicatorManager::CommunicatorManager(bool initialize_mpi) : init_mpi_(initialize_mpi) {
      if (init_mpi_) {
        int provided;
        MPI_Init_thread(nullptr, nullptr, MPI_THREAD_SERIALIZED, &provided);
        if (provided != MPI_THREAD_SERIALIZED)
          ;  // TODO: abort
      }
    }

    CommunicatorManager::~CommunicatorManager() {
      // Make sure that the grids are destroyed before calling MPI_Finalize
      comm_grid_map.clear();
#ifdef DLA_HAVE_SCALAPACK
      ictxt_grid_map.clear();
#endif

      if (init_mpi_) {
        MPI_Finalize();
      }
    }
  }
}
