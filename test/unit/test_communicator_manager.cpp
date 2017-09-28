#include "communicator_manager.h"

#include <mpi.h>
#include "gtest/gtest.h"
#include "communicator_grid.h"
#include "util_mpi.h"

// This tests have to be executed using 6 MPI ranks.

std::vector<MPI_Comm> mpi_comms;
using namespace dla_interface;
using namespace testing;

TEST(CommunicatorManager, Test) {
  std::vector<comm::Communicator2DGrid*> comms;
  for (auto order : {RowMajor, ColMajor}) {
    auto id_size_pair = commInfo();
    const int size = id_size_pair.second;
    int size1 = (size % 2 == 0 ? 2 : 1);
    int size2 = size / size1;

    EXPECT_THROW(
        comm::CommunicatorManager::createCommunicator2DGrid(MPI_COMM_WORLD, size1 + 1, size2, order),
        std::invalid_argument);
    {
      auto& comm(
          comm::CommunicatorManager::createCommunicator2DGrid(MPI_COMM_WORLD, size1, size2, order));
      EXPECT_EQ(std::make_pair(size1, size2), comm.size2D());
      EXPECT_EQ(order, comm.rankOrder());
      comms.push_back(&comm);
    }
    if (size1 != size2) {
      auto& comm(
          comm::CommunicatorManager::createCommunicator2DGrid(MPI_COMM_WORLD, size2, size1, order));
      EXPECT_EQ(std::make_pair(size2, size1), comm.size2D());
      EXPECT_EQ(order, comm.rankOrder());
      comms.push_back(&comm);
    }
  }

#ifdef DLA_HAVE_SCALAPACK
  int ictxt_max = 0;
#endif
  for (auto comm : comms) {
    EXPECT_EQ(comm, &(comm::CommunicatorManager::getCommunicator2DGridFromMPIComm(
                        comm->colMPICommunicator())));
    EXPECT_EQ(comm, &(comm::CommunicatorManager::getCommunicator2DGridFromMPIComm(
                        comm->rowMPICommunicator())));
    EXPECT_EQ(comm, &(comm::CommunicatorManager::getCommunicator2DGridFromMPIComm(
                        comm->rowOrderedMPICommunicator())));
#ifdef DLA_HAVE_SCALAPACK
    EXPECT_EQ(comm, &(comm::CommunicatorManager::getCommunicator2DGridFromBlacsContext(
                        comm->blacsContext())));
    ictxt_max = std::max(ictxt_max, comm->blacsContext());
#endif
  }
  EXPECT_THROW(comm::CommunicatorManager::getCommunicator2DGridFromMPIComm(MPI_COMM_NULL),
               std::invalid_argument);
#ifdef DLA_HAVE_SCALAPACK
  // No BLACS context with id ictxt_max + 1 exists.
  EXPECT_THROW(comm::CommunicatorManager::getCommunicator2DGridFromBlacsContext(ictxt_max + 1),
               std::invalid_argument);
#endif
}

int main(int argc, char** argv) {
#ifdef COMM_INITS_MPI
  comm::CommunicatorManager::initialize(true);
#else
  int provided = 0;
  MPI_Init_thread(nullptr, nullptr, MPI_THREAD_SERIALIZED, &provided);
  if (MPI_THREAD_SERIALIZED != provided) {
    std::cout << "MPI_Init_thread error!" << std::endl;
    MPI_Finalize();
    return 1;
  }

  comm::CommunicatorManager::initialize(false);
#endif

  auto id_size_pair = commInfo();
  int size = id_size_pair.second;
  if (size != 6) {
    std::cout << "This test need 6 MPI ranks (" << size << " provided)!" << std::endl;
#ifndef COMM_INITS_MPI
    MPI_Finalize();
#endif
    return 1;
  }

  ::testing::InitGoogleTest(&argc, argv);

  auto ret = RUN_ALL_TESTS();

  comm::CommunicatorManager::finalize();
#ifndef COMM_INITS_MPI
  MPI_Finalize();
#endif
  return ret;
}
