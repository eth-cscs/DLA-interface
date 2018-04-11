#include <algorithm>
#include <cassert>
#include <cstring>
#include <map>
#include <memory>
#include <mpi.h>
#include <stdexcept>
#include <string>
#include "blacs.h"
#include "communicator_grid.h"
#include "communicator_manager.h"
#include "internal_error.h"
#include "types.h"

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
    Communicator2DGrid& CommunicatorManager::createCommunicator2DGridBlacs(  //
        int blacs_handle, int nr_rows, int nr_cols, Ordering comm_ordering) {
      return comm_manager_->communicator2DGrid(blacs::Cblacs2sys_handle(blacs_handle), nr_rows,
                                               nr_cols, comm_ordering);
    }
#endif

    Communicator2DGrid& CommunicatorManager::getCommunicator2DGridFromMPIComm(MPI_Comm comm) {
      return comm_manager_->communicator2DGridFromMPIComm(comm);
    }

#ifdef DLA_HAVE_SCALAPACK
    Communicator2DGrid& CommunicatorManager::getCommunicator2DGridFromBlacsContext(BlacsContextType id) {
      return comm_manager_->communicator2DGridFromBlacsContext(id);
    }
#endif

    void CommunicatorManager::free2DGrid(Communicator2DGrid& grid) {
      return comm_manager_->destroy2DGrid(grid);
    }
    void CommunicatorManager::free2DGridFromMPIComm(MPI_Comm comm) {
      free2DGrid(getCommunicator2DGridFromMPIComm(comm));
    }
#ifdef DLA_HAVE_SCALAPACK
    void CommunicatorManager::free2DGridFromBlacsContext(BlacsContextType id) {
      free2DGrid(getCommunicator2DGridFromBlacsContext(id));
    }
#endif

    // non-static members:
    Communicator2DGrid& CommunicatorManager::communicator2DGrid(MPI_Comm base_comm, int row_size,
                                                                int col_size, Ordering comm_ordering) {
      std::shared_ptr<Communicator2DGrid> comm_grid(
          new Communicator2DGrid(base_comm, row_size, col_size, comm_ordering));

      if (!comm_grid_map_.insert(std::make_pair(comm_grid->rowMPICommunicator(), comm_grid)).second)
        throw error::InternalError("Cannot insert row commmunicator in the map");
      if (!comm_grid_map_.insert(std::make_pair(comm_grid->colMPICommunicator(), comm_grid)).second)
        throw error::InternalError("Cannot insert column commmunicator in the map");
      if (!comm_grid_map_.insert(std::make_pair(comm_grid->rowOrderedMPICommunicator(), comm_grid)).second)
        throw error::InternalError("Cannot insert row ordered commmunicator in the map");
#ifdef DLA_HAVE_SCALAPACK
      if (!ictxt_grid_map_.insert(std::make_pair(comm_grid->blacsContext(), comm_grid)).second)
        throw error::InternalError("Cannot insert BLACS context id in the map");
#endif
      return *comm_grid;
    }
    Communicator2DGrid& CommunicatorManager::communicator2DGridFromMPIComm(MPI_Comm comm) const {
      try {
        return *comm_grid_map_.at(comm);
      }
      catch (std::out_of_range) {
        throw std::invalid_argument("No communicator2DGrid found with the given MPI_Comm");
      }
    }
#ifdef DLA_HAVE_SCALAPACK
    Communicator2DGrid& CommunicatorManager::communicator2DGridFromBlacsContext(BlacsContextType id) const {
      try {
        return *ictxt_grid_map_.at(id);
      }
      catch (std::out_of_range) {
        throw std::invalid_argument("No communicator2DGrid found with the given BLACS context id");
      }
    }
#endif

    void CommunicatorManager::destroy2DGrid(Communicator2DGrid& grid) {
#ifdef DLA_HAVE_SCALAPACK
      // Internal Check of the number of shared pointer
      assert(comm_grid_map_.at(grid.rowOrderedMPICommunicator()).use_count() == 4);
      ictxt_grid_map_.erase(grid.blacsContext());
#endif
      // Internal Check of the number of shared pointer
      assert(comm_grid_map_.at(grid.rowOrderedMPICommunicator()).use_count() == 3);
      comm_grid_map_.erase(grid.rowMPICommunicator());
      comm_grid_map_.erase(grid.colMPICommunicator());
      // Internal Check of the number of shared pointer (only one remaining)
      assert(comm_grid_map_.at(grid.rowOrderedMPICommunicator()).use_count() == 1);
      comm_grid_map_.erase(grid.rowOrderedMPICommunicator());
    }

    CommunicatorManager::CommunicatorManager(int nr_cores, int* argc, char*** argv,
                                             bool initialize_mpi)
     : init_mpi_(initialize_mpi), topo_() {
      if (init_mpi_) {
        int provided;
        MPI_Init_thread(argc, argv, MPI_THREAD_SERIALIZED, &provided);
        if (provided != MPI_THREAD_SERIALIZED)
          throw std::runtime_error(
              "MPI cannot be initializd with MPI_THREAD_SERIALIZED. (provided = " +
              std::to_string(provided) + ")");
      }

      thread::CpuSet application_cpuset = topo_.getCpuBind();
#ifdef DLA_HAVE_SCALAPACK
      scalapack_cpuset_ = application_cpuset;
      scalapack_nr_threads_ = thread::getOmpBlasThreads();
#endif

#ifdef DLA_HAVE_DPLASMA
      if (argc == nullptr) {
        parsec_handle_ = parsec_init(nr_cores, nullptr, nullptr);
      }
      else {
        char** p_argv = std::find_if(*argv, (*argv) + (*argc),
                                     [](const char* val) { return std::strcmp(val, "--") == 0; });
        char* tmp = nullptr;
        int p_argc = *argc - (p_argv - *argv);
        if (p_argc > 0) {
          // "--" found. Replace it with the application name and pass to parsec_init.
          tmp = *p_argv;
          *p_argv = (*argv)[0];
        }
        else {
          // No arguments for parsec found, just pass to parsec_init the application name.
          p_argv = *argv;
          p_argc = 1;
        }

        parsec_handle_ = parsec_init(nr_cores, &p_argc, &p_argv);

        // Restore the "--" argument if it was replaced.
        if (tmp)
          *p_argv = tmp;
      }

      dplasma_cpuset_ = topo_.getCpuBind();
      dplasma_nr_threads_ = thread::NumThreads(1);
#endif
      topo_.setCpuBind(application_cpuset);
    }

    CommunicatorManager::~CommunicatorManager() {
      // Make sure that the grids are destroyed before calling MPI_Finalize
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      std::stringstream s;
      fall_back_info_.finalReport(rank, s);
      std::cerr << s.str() << std::endl;

      comm_grid_map_.clear();
#ifdef DLA_HAVE_SCALAPACK
      ictxt_grid_map_.clear();
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
