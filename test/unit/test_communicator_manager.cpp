#include "communicator_manager.h"

#include <mpi.h>
#include <stdexcept>
#include <vector>
#include "gtest/gtest.h"
#include "mpi_listener.h"
#include "communicator_grid.h"
#include "util_mpi.h"

#ifdef DLA_HAVE_SCALAPACK
#include "blacs.h"
#endif

// This tests have to be executed using 6 MPI ranks.

std::vector<MPI_Comm> mpi_comms;
using namespace dla_interface;
using namespace testing;

TEST(CommunicatorManager, Test) {
  std::vector<comm::Communicator2DGrid*> comms;
  for (const auto& mpi_comm : mpi_comms) {
    auto id_size_pair = commInfo(mpi_comm);
    const int id = id_size_pair.first;
    const int size = id_size_pair.second;
    for (auto order : {RowMajor, ColMajor}) {
      int size1 = (size % 2 == 0 ? 2 : 1);
      int size2 = size / size1;
      auto ref_id2D(order == ColMajor
                        ? [](int id, int ld) { return std::make_pair(id % ld, id / ld); }
                        : [](int id, int ld) { return std::make_pair(id / ld, id % ld); });

      EXPECT_THROW(
          comm::CommunicatorManager::createCommunicator2DGrid(mpi_comm, size1 + 1, size2, order),
          std::invalid_argument);
      {
        auto& comm(comm::CommunicatorManager::createCommunicator2DGrid(mpi_comm, size1, size2, order));
        EXPECT_EQ(std::make_pair(size1, size2), comm.size2D());
        EXPECT_EQ(size1 * size2, comm.size());
        EXPECT_EQ(ref_id2D(id, (order == ColMajor ? size1 : size2)), comm.id2D());
        EXPECT_EQ(order, comm.rankOrder());
        comms.push_back(&comm);
      }
      if (size1 != size2) {
        auto& comm(comm::CommunicatorManager::createCommunicator2DGrid(mpi_comm, size2, size1, order));
        EXPECT_EQ(std::make_pair(size2, size1), comm.size2D());
        EXPECT_EQ(size1 * size2, comm.size());
        EXPECT_EQ(ref_id2D(id, (order == ColMajor ? size2 : size1)), comm.id2D());
        EXPECT_EQ(order, comm.rankOrder());
        comms.push_back(&comm);
      }
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

  for (size_t i = 0; i < comms.size(); ++i) {
    auto comm = comms[i];
    MPI_Comm row_comm = comm->rowMPICommunicator();
    MPI_Comm col_comm = comm->colMPICommunicator();
    MPI_Comm row_ordered_comm = comm->rowOrderedMPICommunicator();
#ifdef DLA_HAVE_SCALAPACK
    int ictxt = comm->blacsContext();
#endif
    switch (i % 5) {
      case 1:
        comm::CommunicatorManager::free2DGridFromMPIComm(row_comm);
        break;
      case 2:
        comm::CommunicatorManager::free2DGridFromMPIComm(col_comm);
        break;
      case 3:
        comm::CommunicatorManager::free2DGridFromMPIComm(row_ordered_comm);
        break;
#ifdef DLA_HAVE_SCALAPACK
      case 4:
        comm::CommunicatorManager::free2DGridFromBlacsContext(ictxt);
        break;
#endif
      default:
        comm::CommunicatorManager::free2DGrid(*comm);
    }

    EXPECT_THROW(comm::CommunicatorManager::getCommunicator2DGridFromMPIComm(row_comm),
                 std::invalid_argument);
#ifdef DLA_HAVE_SCALAPACK
    EXPECT_THROW(comm::CommunicatorManager::getCommunicator2DGridFromBlacsContext(ictxt),
                 std::invalid_argument);
#endif
  }
}

#ifdef DLA_HAVE_SCALAPACK
TEST(CommunicatorManager, TestBlacsConstr) {
  std::vector<comm::Communicator2DGrid*> comms;
  for (const auto& mpi_comm : mpi_comms) {
    auto id_size_pair = commInfo(mpi_comm);
    const int id = id_size_pair.first;
    const int size = id_size_pair.second;
    for (auto order : {RowMajor, ColMajor}) {
      int size1 = (size % 2 == 0 ? 2 : 1);
      int size2 = size / size1;
      auto ref_id2D(order == ColMajor
                        ? [](int id, int ld) { return std::make_pair(id % ld, id / ld); }
                        : [](int id, int ld) { return std::make_pair(id / ld, id % ld); });

      int blacs_handle = blacs::Csys2blacs_handle(mpi_comm);
      EXPECT_THROW(comm::CommunicatorManager::createCommunicator2DGridBlacs(blacs_handle, size1 + 1,
                                                                            size2, order),
                   std::invalid_argument);
      {
        auto& comm(comm::CommunicatorManager::createCommunicator2DGridBlacs(blacs_handle, size1,
                                                                            size2, order));
        EXPECT_EQ(std::make_pair(size1, size2), comm.size2D());
        EXPECT_EQ(size1 * size2, comm.size());
        EXPECT_EQ(ref_id2D(id, (order == ColMajor ? size1 : size2)), comm.id2D());
        EXPECT_EQ(order, comm.rankOrder());
        comms.push_back(&comm);
      }
      if (size1 != size2) {
        auto& comm(comm::CommunicatorManager::createCommunicator2DGridBlacs(blacs_handle, size2,
                                                                            size1, order));
        EXPECT_EQ(std::make_pair(size2, size1), comm.size2D());
        EXPECT_EQ(size1 * size2, comm.size());
        EXPECT_EQ(ref_id2D(id, (order == ColMajor ? size2 : size1)), comm.id2D());
        EXPECT_EQ(order, comm.rankOrder());
        comms.push_back(&comm);
      }
    }
  }

  int ictxt_max = 0;
  for (auto comm : comms) {
    EXPECT_EQ(comm, &(comm::CommunicatorManager::getCommunicator2DGridFromMPIComm(
                        comm->colMPICommunicator())));
    EXPECT_EQ(comm, &(comm::CommunicatorManager::getCommunicator2DGridFromMPIComm(
                        comm->rowMPICommunicator())));
    EXPECT_EQ(comm, &(comm::CommunicatorManager::getCommunicator2DGridFromMPIComm(
                        comm->rowOrderedMPICommunicator())));
    EXPECT_EQ(comm, &(comm::CommunicatorManager::getCommunicator2DGridFromBlacsContext(
                        comm->blacsContext())));
    ictxt_max = std::max(ictxt_max, comm->blacsContext());
  }
  EXPECT_THROW(comm::CommunicatorManager::getCommunicator2DGridFromMPIComm(MPI_COMM_NULL),
               std::invalid_argument);
  // No BLACS context with id ictxt_max + 1 exists.
  EXPECT_THROW(comm::CommunicatorManager::getCommunicator2DGridFromBlacsContext(ictxt_max + 1),
               std::invalid_argument);
}
#endif

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
  int rank = id_size_pair.first;
  int size = id_size_pair.second;
  if (size != 6) {
    std::cout << "This test need 6 MPI ranks (" << size << " provided)!" << std::endl;
#ifndef COMM_INITS_MPI
    MPI_Finalize();
#endif
    return 1;
  }

  mpi_comms.push_back(MPI_COMM_WORLD);

  MPI_Comm tmp;
  MPI_Comm_split(MPI_COMM_WORLD, (rank == 1 || rank == 5 ? 1 : 2), size - rank, &tmp);
  mpi_comms.push_back(tmp);
  MPI_Comm_split(MPI_COMM_WORLD, rank % 2, rank, &tmp);
  mpi_comms.push_back(tmp);

  ::testing::InitGoogleTest(&argc, argv);

#ifdef COMM_INITS_MPI
  ::testing::setMPIListener("results_test_communicator_manager_init");
#else
  ::testing::setMPIListener("results_test_communicator_manager");
#endif

  auto ret = RUN_ALL_TESTS();

  for (auto& comm : mpi_comms)
    if (comm != MPI_COMM_WORLD)
      MPI_Comm_free(&comm);

  comm::CommunicatorManager::finalize();
#ifndef COMM_INITS_MPI
  MPI_Finalize();
#endif
  return ret;
}
