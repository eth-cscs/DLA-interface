//
// Distributed Linear Algebra Interface (DLAI)
//
// Copyright (c) 2018-2019, ETH Zurich
// All rights reserved.
//
// Please, refer to the LICENSE file in the root directory.
// SPDX-License-Identifier: BSD-3-Clause
//

#ifndef DLA_INTERFACE_DLA_DLAF_H_
#define DLA_INTERFACE_DLA_DLAF_H_

#ifdef DLA_HAVE_DLAF

#include <string>
#include <vector>

#include <mpi.h>
#include <hpx/hpx.hpp>
#include <hpx/hpx_start.hpp>
#include <hpx/hpx_suspend.hpp>

#include <dlaf/common/index2d.h>
#include <dlaf/communication/communicator.h>
#include <dlaf/communication/communicator_grid.h>
#include <dlaf/matrix/matrix.h>
#include <dlaf/factorization/cholesky.h>
#include "dlaf/types.h"
#include "dlaf/util_matrix.h"

#include "communicator_grid.h"
#include "communicator_manager.h"
#include "distributed_matrix.h"

namespace dla_interface {

  namespace hpx_wrappers {
    inline void start(int argc, char** argv, std::vector<std::string> cfg) {
      using namespace hpx::program_options;

      hpx::init_params p;
      p.cfg = std::move(cfg);
      p.rp_callback = [](auto& rp, auto) {
        int ntasks;
        DLAF_MPI_CALL(MPI_Comm_size(MPI_COMM_WORLD, &ntasks));
        // if the user has asked for special thread pools for communication
        // then set them up
        if (ntasks > 1) {
          // Create a thread pool with a single core that we will use for all
          // communication related tasks
          rp.create_thread_pool("mpi", hpx::resource::scheduling_policy::local_priority_fifo);
          rp.add_resource(rp.numa_domains()[0].cores()[0].pus()[0], "mpi");
        }
      };

      hpx::start(nullptr, argc, argv, p);
      hpx::suspend();
    }

    inline void stop() {
      hpx::resume();
      hpx::apply([]() { hpx::finalize(); });
      hpx::stop();
    }
  }

  namespace dlaf_wrappers {

    typedef dlaf::comm::Communicator Communicator;
    typedef dlaf::comm::CommunicatorGrid CommunicatorGrid;

    // Wrapper class to DLA-Future CommunicatorGrid class
    inline CommunicatorGrid comm_grid(const comm::Communicator2DGrid& comm_grid)
    {
      dlaf::common::Ordering dlaf_ordering = dlaf::common::Ordering::RowMajor;
      Communicator dlaf_communicator;

      if(comm_grid.rankOrder() == Ordering::RowMajor) {
        dlaf_communicator = comm_grid.rowMPICommunicator();
        dlaf_ordering = dlaf::common::Ordering::RowMajor;
      }
      else if(comm_grid.rankOrder() == Ordering::ColMajor) {
        dlaf_communicator = comm_grid.colMPICommunicator();
        dlaf_ordering = dlaf::common::Ordering::ColumnMajor;
      } else {
        // Throw something! Ask what!
        std::cout << "Ordering not recognized!  " << std::endl;
        exit(1);
      }

      dlaf::comm::IndexT_MPI grid_rows = comm_grid.size2D().first;
      dlaf::comm::IndexT_MPI grid_cols = comm_grid.size2D().second;

      return CommunicatorGrid(dlaf_communicator, grid_rows, grid_cols, dlaf_ordering);
    }

    inline CommunicatorGrid comm_grid(MPI_Comm base_comm, int grid_rows, int grid_cols, Ordering comm_ordering)
    {
      dlaf::common::Ordering dlaf_ordering = dlaf::common::Ordering::RowMajor;
      Communicator dlaf_communicator(base_comm);

      if(comm_ordering == Ordering::RowMajor)
        dlaf_ordering = dlaf::common::Ordering::RowMajor;
      else if(comm_ordering == Ordering::ColMajor)
        dlaf_ordering = dlaf::common::Ordering::ColumnMajor;
      else {
        // Throw something! Ask what!
        std::cout << "Ordering not recognized!  " << std::endl;
        exit(1);
      }


      return CommunicatorGrid(dlaf_communicator, grid_rows, grid_cols, dlaf_ordering);
    }

    // Distributed matrix wrapper
    // Device:
    //   - dlaf::Device::CPU
    //   - dlaf::Device::GPU
    template <class Type, dlaf::Device Device>
    inline dlaf::Matrix<Type, Device> matrix(DistributedMatrix<Type>& mat, const CommunicatorGrid& comm_grid){

      // Allocate memory for the matrix
      dlaf::GlobalElementSize matrix_size(mat.size().first, mat.size().second);
      dlaf::TileElementSize block_size(mat.blockSize().first, mat.blockSize().second);
      dlaf::SizeType ld = mat.size().first;
      Type* mem_ptr = mat.ptr();

      //dlaf::Matrix<Type, Device> dlaf_matrix(matrix_size, block_size, comm_grid);

      dlaf::Matrix<Type, Device> dlaf_matrix = dlaf::matrix::createMatrixFromColMajor<Device, Type>(
          matrix_size,
          block_size,
          ld,
          comm_grid,
          mem_ptr);

      return dlaf_matrix;
    }

    // TODO: if needed implement also dlaf::Matrix -> dlai::DistributedMatrix
    /*
    template <class Type, dlaf::Device Device>
    inline DistributedMatrix<Type> matrix(const dlaf::Matrix<Type, Device>& mat, const comm::Communicator2DGrid& comm_grid){

      // Allocate memory for the matrix
      int n = mat.size().rows();
      int m = mat.size().cols();
      int nb = mat.blockSize().rows();
      int mb = mat.blockSize().cols();

      DistributedMatrix<Type> dlai_matrix(n, m, nb, mb, comm_grid);

      // Create a shallow copy of dlaf::Matrix with elements from dla_interface::DistributedMatrix
      // Here needed help

      return dlai_matrix;
    }
    */

  } // dlaf_wrappers
} // dla_interface

#endif // DLA_HAVE_DLAF
#endif /* DLA_INTERFACE__DLA_DLAF_H_ */
