#include <mpi.h>
#include <vector>
#include "dla_interface.h"
#include "distributed_matrix.h"

using namespace dla_interface;

int main(int argc, char** argv) {

	const int n = 4096;
	const int nb = 256;

	const int p = 1;
	const int q = 1;
	const int nr_threads = 1;

	comm::CommunicatorManager::initialize(nr_threads, &argc, &argv, true);

	auto& comm_grid = comm::CommunicatorManager::createCommunicator2DGrid(MPI_COMM_WORLD, p, q, RowMajor);

	DistributedMatrix<double> mat(n, n, nb, nb, comm_grid, scalapack_dist);
	std::unique_ptr<DistributedMatrix<double>> mat_copy;

#ifdef DLA_WITH_SCALAPACK
	std::cout << "solver: SCALAPACK" << std::endl;
#endif

#ifdef DLA_WITH_DLAF
	std::cout << "solver: DLAF" << std::endl;
#endif

	return 0;
}
