#include <fstream>
#include <iostream>
#include "gtest/gtest.h"
#include "mpi_listener.h"
#include "communicator_grid.h"
#include "communicator_manager.h"
#include "test_dlai_main.h"

using namespace dla_interface;
using namespace testing;

std::vector<comm::Communicator2DGrid*> comms;
std::ofstream* outstream;

int main(int argc, char** argv) {
  comm::CommunicatorManager::initialize(2, &argc, &argv, true);

  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 6) {
    std::cout << "This test need 6 MPI ranks (" << size << " provided)!" << std::endl;
    return 1;
  }

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Create communicators used in the tests.
  for (auto order : {RowMajor, ColMajor}) {
    comms.push_back(&comm::CommunicatorManager::createCommunicator2DGrid(MPI_COMM_WORLD, 2, 3, order));
    comms.push_back(&comm::CommunicatorManager::createCommunicator2DGrid(MPI_COMM_WORLD, 3, 2, order));
    comms.push_back(&comm::CommunicatorManager::createCommunicator2DGrid(MPI_COMM_WORLD, 1, 6, order));
  }

  comm::CommunicatorManager::getFallbackInfo().setFullReport(true);

  ::testing::InitGoogleTest(&argc, argv);

  outstream = &::testing::setMPIListener("results_test_dlai_matrix_multiplication");

  auto ret = RUN_ALL_TESTS();

  comm::CommunicatorManager::finalize();

  return ret;
}
