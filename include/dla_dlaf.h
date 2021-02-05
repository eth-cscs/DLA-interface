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

#include <mpi.h>
#include <dlaf/communication/communicator_grid.h>
#include <dlaf/matrix/matrix.h>

namespace dla_interface {
	namespace dlaf_wrappers {

	/*
	typedef dlaf::comm::Communicator Communicator;
	typedef dlaf::comm::CommunicatorGrid CommunicatorGrid;

	CommunicatorGrid comm_grid(int grid_rows, int grid_cols, dlaf::common::Ordering ordering) {

		Communicator world(MPI_COMM_WORLD);

		return CommunicatorGrid(world, grid_rows, grid_cols, ordering);
	}

	class CommmunicatorGrid{

	};

	class DistributedMatrix {

	};
	*/

	} // dlaf_wrappers
} // dla_interface

#endif // DLA_HAVE_DLAF
#endif /* DLA_INTERFACE__DLA_DLAF_H_ */
