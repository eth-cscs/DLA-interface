#ifndef DLA_INTERFACE_COMMUNICATOR_MANAGER_H
#define DLA_INTERFACE_COMMUNICATOR_MANAGER_H

#include <map>
#include <memory>
#include <mpi.h>
#include <tuple>
#include "communicator_grid.h"
#include "fallback_info.h"
#include "thread_binding.h"
#include "thread_blas.h"
#include "types.h"

namespace dla_interface {
  namespace comm {

  /// A CommunicatorManager class is a manager for MPI communication.

    class CommunicatorManager {
      public:
      /// Initializes the CommunicatorManager.
      /// If DLA_HAVE_DPLASMA is true parsec is initialized with nr_cores (default -1).
      /// If initialize_mpi is true:
      ///- MPI_Init_thread is called using MPI_THREAD_SERIALIZED,
      ///- MPI_Finalize is called when finalize is called.
      /// **Precondition:** initialize must have not been called before. (It can be called only once.)
      /// If specified argc and argv are used in MPI and DPLASMA (if DLA_HAVE_DPLASMA is true)
      /// initialization (arguments after -- are passed to Parsec).
      /// 
      /// @param initialize_mpi Is MPI initialized.
       
      static void initialize(bool initialize_mpi = true);

      ///  Initializes the CommunicatorManager.
      /// If DLA_HAVE_DPLASMA is true parsec is initialized with nr_cores (default -1).
      /// If initialize_mpi is true:
      /// - MPI_Init_thread is called using MPI_THREAD_SERIALIZED,
      /// - MPI_Finalize is called when finalize is called.
      /// **Precondition:** initialize must have not been called before. (It can be called only once.)
      /// If specified argc and argv are used in MPI and DPLASMA (if DLA_HAVE_DPLASMA is true)
      /// initialization (arguments after -- are passed to Parsec).
      /// 
      /// @param nr_cores Set number of cores to be used.
      /// @param initialize_mpi Is MPI initialized.
       
      static void initialize(int nr_cores, bool initialize_mpi = true);

      /// Initializes the CommunicatorManager.
      /// If DLA_HAVE_DPLASMA is true parsec is initialized with nr_cores (default -1).
      /// If initialize_mpi is true:
      /// - MPI_Init_thread is called using MPI_THREAD_SERIALIZED,
      /// - MPI_Finalize is called when finalize is called.
      /// Precondition: initialize must have not been called before. (It can be called only once.)
      /// If specified argc and argv are used in MPI and DPLASMA (if DLA_HAVE_DPLASMA is true)
      /// initialization (arguments after -- are passed to Parsec).
      /// 
      /// @param nr_cores Set number of cores to be used.
      /// @param argc Number of arguments.
      /// @param argv Arguments values.
      /// @param initialize_mpi Is MPI initialized.
    
      static void initialize(int nr_cores, int* argc, char*** argv, bool initialize_mpi = true);

      /// Finalize communication. Destructs all the 2Dgrids created with createCommunicator2DGrid() and
      /// if initialize was called with initialize_mpi == true this function finalize MPI as well.
      /// **Precondition:** initialize must have been called before and finalize must have not been
      /// called before. (It can be called only once.)
       
      static void finalize();
#ifdef DLAI_WITH_DPLASMA
      /// Returns the parsec context.
      static ParsecContext getParsecContext();
#endif

      /// Creates a 2D grid based on the given MPI communicator.
      /// Creates a 2D grid based on the given MPI communicator, with nr_rows rows
      /// and nr_cols columns, where the ranks are ordered in comm_ordering order.
      /// The created Communicator2DGrid object is stored and a reference to it is returned.
      /// **Precondition:** initialize must have been called before and finalize must have not been called before.
      /// 
      /// @throws std::invalid_argument if the number of ranks in base_comm is different from nr_rows * nr_cols.
      /// 
      /// @param base_comm MPI base communicator.
      /// @param nr_rows Number of rows.
      /// @param nr_cols Number of columns.
      /// @param comm_ordering Ordering type (Row or Column based).
      /// @return Communicator2DGrid&
       
      static Communicator2DGrid& createCommunicator2DGrid(MPI_Comm base_comm, int nr_rows,
                                                          int nr_cols, Ordering comm_ordering);

#ifdef DLAI_WITH_SCALAPACK
      /// Creates a 2D grid based on the MPI communicator base_comm = Cblacs2sys_handle(blacs_handle).
      /// Calls createCommunicator2DGrid with base_comm = Cblacs2sys_handle(blacs_handle),
      /// with nr_rows rows and nr_cols columns, where the ranks are ordered in comm_ordering order.
      /// The created Communicator2DGrid object is stored and a reference to it is returned.
      /// **Precondition:** initialize must have been called before and finalize must have not been called before.
      /// 
      /// @throws std::invalid_argument if the number of ranks in base_comm is different from nr_rows * nr_cols.
      /// 
      /// @param blacs_handle Pointer to BLACS.
      /// @param nr_rows Number of rows.
      /// @param nr_cols Number of columns.
      /// @param comm_ordering Ordering type (Row or Column based).
      /// @return Communicator2DGrid&
       
      static Communicator2DGrid& createCommunicator2DGridBlacs(int blacs_handle, int nr_rows,
                                                               int nr_cols, Ordering comm_ordering);
#endif

      /// Returns a reference of the Communicator2DGrid.
      /// Returns a reference of the Communicator2DGrid created with createCommunicator2DGrid
      /// whose row, col or row_ordered (TODO: or col_ordered) communicator is comm.
      /// 
      /// @throws std::invalid_argument if no Communicator2DGrid is found.
      /// 
      /// @param comm MPI base communicator.
      /// @return Communicator2DGrid&
       
      static Communicator2DGrid& getCommunicator2DGridFromMPIComm(MPI_Comm comm);

#ifdef DLAI_WITH_SCALAPACK
      /// 
      /// Returns a reference of the Communicator2DGrid created with createCommunicator2DGrid
      /// whose BLACS context index is id.
      /// 
      /// @throws std::invalid_argument if no Communicator2DGrid is found.
      /// 
      /// @param id BlacsContextType
    
      static Communicator2DGrid& getCommunicator2DGridFromBlacsContext(BlacsContextType id);
#endif

      /// Free the Communicator2DGrid object and the related MPI communicators and BLACS grids.
      static void free2DGrid(Communicator2DGrid& grid);

      /// Free the Communicator2DGrid object and the related MPI communicators and BLACS grids.
      static void free2DGridFromMPIComm(MPI_Comm comm);
#ifdef DLAI_WITH_SCALAPACK
      /// Free the Communicator2DGrid object and the related MPI communicators and BLACS grids.
      static void free2DGridFromBlacsContext(BlacsContextType id);
#endif

      /// Returns the number of threads and cpuset for Scalapack.
#ifdef DLAI_WITH_SCALAPACK
      static std::tuple<const thread::NumThreads&, const thread::CpuSet&> getScalapackConfigInfo() {
        return std::make_tuple(std::cref(commManager()->scalapack_nr_threads_),
                               std::cref(commManager()->scalapack_cpuset_));
      }
#endif
#ifdef DLAI_WITH_DPLASMA
      /// Returns the number of threads and cpuset for DPlasma.
      static std::tuple<const thread::NumThreads&, const thread::CpuSet&> getDPlasmaConfigInfo() {
        return std::make_tuple(std::cref(commManager()->dplasma_nr_threads_),
                               std::cref(commManager()->dplasma_cpuset_));
      }
#endif
#ifdef DLA_HAVE_DLAF
      /// Returns the number of threads and cpuset for DLA-Future.
      static std::tuple<const thread::NumThreads&, const thread::CpuSet&> getDLAFConfigInfo() {
        return std::make_tuple(std::cref(commManager()->dlaf_nr_threads_),
                               std::cref(commManager()->dlaf_cpuset_));
      }
#endif

      /// Returns the thread CpuSet binding.
      static thread::CpuSet getCpuBind() {
        return commManager()->getCpuBindInternal();
      }

      /// Sets the thread CpuSet binding.
      static void setCpuBind(const thread::CpuSet& cpuset) {
        commManager()->setCpuBindInternal(cpuset);
      }

      /// Returns the Fall back information.
      static FallbackInfo& getFallbackInfo() {
        return commManager()->fall_back_info_;
      }

      protected:
      Communicator2DGrid& communicator2DGrid(MPI_Comm base_comm, int row_size, int col_size,
                                             Ordering comm_ordering);
      Communicator2DGrid& communicator2DGridFromMPIComm(MPI_Comm comm) const;
#ifdef DLAI_WITH_SCALAPACK
      Communicator2DGrid& communicator2DGridFromBlacsContext(BlacsContextType id) const;
#endif
      void destroy2DGrid(Communicator2DGrid& grid);

      thread::CpuSet getCpuBindInternal() {
        return topo_.getCpuBind();
      }

      void setCpuBindInternal(const thread::CpuSet& cpuset) {
        topo_.setCpuBind(cpuset);
      }

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
      static std::unique_ptr<CommunicatorManager>& commManager() {
        static std::unique_ptr<CommunicatorManager> comm_manager_;
        return comm_manager_;
      }

      bool init_mpi_;
#ifdef DLAI_WITH_DPLASMA
      ParsecContext parsec_handle_;
#endif

      std::map<MPI_Comm, std::shared_ptr<Communicator2DGrid>> comm_grid_map_;
#ifdef DLAI_WITH_SCALAPACK
      std::map<BlacsContextType, std::shared_ptr<Communicator2DGrid>> ictxt_grid_map_;
#endif

      thread::SystemTopology topo_;
#ifdef DLAI_WITH_SCALAPACK
      thread::NumThreads scalapack_nr_threads_;
      thread::CpuSet scalapack_cpuset_;
#endif
      
#ifdef DLAI_WITH_DPLASMA
      thread::NumThreads dplasma_nr_threads_;
      thread::CpuSet dplasma_cpuset_;
#endif
      
#ifdef DLAI_WITH_DLAF
      thread::NumThreads dlaf_nr_threads_;
      thread::CpuSet dlaf_cpuset_;
#endif
      
      FallbackInfo fall_back_info_;
    };
  }
}

#endif  // DLA_INTERFACE_COMMUNICATOR_MANAGER_H
