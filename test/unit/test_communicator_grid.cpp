#include "communicator_grid.h"

#include <mpi.h>
#include <stdexcept>
#include "blacs.h"
#include "gtest/gtest.h"
#include "util_mpi.h"

// This tests have to be execuuted using 6 MPI ranks.

std::vector<MPI_Comm> mpi_comms;
using namespace dla_interface;
using namespace testing;

TEST(Communicato2DGridTest, ConstructorColMajor) {
  for (const auto& mpi_comm : mpi_comms) {
    auto id_size_pair = commInfo(mpi_comm);
    const int id = id_size_pair.first;
    const int size = id_size_pair.second;
    int size1 = (size % 2 == 0 ? 2 : 1);
    int size2 = size / size1;
    EXPECT_THROW(comm::Communicator2DGrid(mpi_comm, size1 - 1, size2, ColMajor),
                 std::invalid_argument);
    EXPECT_THROW(comm::Communicator2DGrid(mpi_comm, size1 + 1, size2, ColMajor),
                 std::invalid_argument);
    EXPECT_THROW(comm::Communicator2DGrid(mpi_comm, size1, size2 - 1, ColMajor),
                 std::invalid_argument);
    EXPECT_THROW(comm::Communicator2DGrid(mpi_comm, size1, size2 + 1, ColMajor),
                 std::invalid_argument);

    auto ref_id2D = [](int id, int ld) { return std::make_pair(id % ld, id / ld); };
    {
      comm::Communicator2DGrid comm(mpi_comm, size1, size2, ColMajor);
      EXPECT_EQ(std::make_pair(size1, size2), comm.size2D());
      EXPECT_EQ(ref_id2D(id, size1), comm.id2D());
      EXPECT_EQ(ColMajor, comm.rankOrder());
    }
    if (size1 != size2) {
      comm::Communicator2DGrid comm(mpi_comm, size2, size1, ColMajor);
      EXPECT_EQ(std::make_pair(size2, size1), comm.size2D());
      EXPECT_EQ(ref_id2D(id, size2), comm.id2D());
      EXPECT_EQ(ColMajor, comm.rankOrder());
    }
  }
}

TEST(Communicato2DGridTest, ConstructorRowMajor) {
  for (const auto& mpi_comm : mpi_comms) {
    auto id_size_pair = commInfo(mpi_comm);
    const int id = id_size_pair.first;
    const int size = id_size_pair.second;
    int size1 = (size % 2 == 0 ? 2 : 1);
    int size2 = size / size1;
    EXPECT_THROW(comm::Communicator2DGrid(mpi_comm, size1 - 1, size2, RowMajor),
                 std::invalid_argument);
    EXPECT_THROW(comm::Communicator2DGrid(mpi_comm, size1 + 1, size2, RowMajor),
                 std::invalid_argument);
    EXPECT_THROW(comm::Communicator2DGrid(mpi_comm, size1, size2 - 1, RowMajor),
                 std::invalid_argument);
    EXPECT_THROW(comm::Communicator2DGrid(mpi_comm, size1, size2 + 1, RowMajor),
                 std::invalid_argument);

    auto ref_id2D = [](int id, int ld) { return std::make_pair(id / ld, id % ld); };
    {
      comm::Communicator2DGrid comm(mpi_comm, size1, size2, RowMajor);
      EXPECT_EQ(std::make_pair(size1, size2), comm.size2D());
      EXPECT_EQ(ref_id2D(id, size2), comm.id2D());
      EXPECT_EQ(RowMajor, comm.rankOrder());
    }
    if (size1 != size2) {
      comm::Communicator2DGrid comm(mpi_comm, size2, size1, RowMajor);
      EXPECT_EQ(std::make_pair(size2, size1), comm.size2D());
      EXPECT_EQ(ref_id2D(id, size1), comm.id2D());
      EXPECT_EQ(RowMajor, comm.rankOrder());
    }
  }
}

TEST(Communicato2DGridTest, MPIColCommunicator) {
  for (const auto& mpi_comm : mpi_comms) {
    auto id_size_pair = commInfo(mpi_comm);
    const int size = id_size_pair.second;
    int size1 = (size % 2 == 0 ? 2 : 1);
    int size2 = size / size1;
    for (auto order : {RowMajor, ColMajor}) {
      {
        comm::Communicator2DGrid comm(mpi_comm, size1, size2, order);
        auto id_size = commInfo(comm.colMPICommunicator());
        EXPECT_EQ(comm.id2D().first, id_size.first);
        EXPECT_EQ(comm.size2D().first, id_size.second);
      }

      if (size1 != size2) {
        comm::Communicator2DGrid comm(mpi_comm, size2, size1, order);
        auto id_size = commInfo(comm.colMPICommunicator());
        EXPECT_EQ(comm.id2D().first, id_size.first);
        EXPECT_EQ(comm.size2D().first, id_size.second);
      }
    }
  }
}

TEST(Communicato2DGridTest, MPIRowCommunicator) {
  for (const auto& mpi_comm : mpi_comms) {
    auto id_size_pair = commInfo(mpi_comm);
    const int size = id_size_pair.second;
    int size1 = (size % 2 == 0 ? 2 : 1);
    int size2 = size / size1;
    for (auto order : {RowMajor, ColMajor}) {
      {
        comm::Communicator2DGrid comm(mpi_comm, size1, size2, order);
        auto id_size = commInfo(comm.rowMPICommunicator());
        EXPECT_EQ(comm.id2D().second, id_size.first);
        EXPECT_EQ(comm.size2D().second, id_size.second);
      }

      if (size1 != size2) {
        comm::Communicator2DGrid comm(mpi_comm, size2, size1, order);
        auto id_size = commInfo(comm.rowMPICommunicator());
        EXPECT_EQ(comm.id2D().second, id_size.first);
        EXPECT_EQ(comm.size2D().second, id_size.second);
      }
    }
  }
}

TEST(Communicato2DGridTest, MPIRowOrderedCommunicator) {
  for (const auto& mpi_comm : mpi_comms) {
    auto id_size_pair = commInfo(mpi_comm);
    const int size = id_size_pair.second;
    int size1 = (size % 2 == 0 ? 2 : 1);
    int size2 = size / size1;
    for (auto order : {RowMajor, ColMajor}) {
      {
        comm::Communicator2DGrid comm(mpi_comm, size1, size2, order);
        auto id_size = commInfo(comm.rowOrderedMPICommunicator());
        EXPECT_EQ(comm.id2D().first * comm.size2D().second + comm.id2D().second, id_size.first);
        EXPECT_EQ(comm.size2D().first * comm.size2D().second, id_size.second);
      }

      if (size1 != size2) {
        comm::Communicator2DGrid comm(mpi_comm, size2, size1, order);
        auto id_size = commInfo(comm.rowOrderedMPICommunicator());
        EXPECT_EQ(comm.id2D().first * comm.size2D().second + comm.id2D().second, id_size.first);
        EXPECT_EQ(comm.size2D().first * comm.size2D().second, id_size.second);
      }
    }
  }
}

#ifdef DLA_HAVE_SCALAPACK
TEST(Communicato2DGridTest, BlacsContext) {
  for (const auto& mpi_comm : mpi_comms) {
    auto id_size_pair = commInfo(mpi_comm);
    const int size = id_size_pair.second;
    int size1 = (size % 2 == 0 ? 2 : 1);
    int size2 = size / size1;
    for (auto order : {RowMajor, ColMajor}) {
      {
        comm::Communicator2DGrid comm(mpi_comm, size1, size2, order);
        auto ictxt = comm.blacsContext();
        int np_col = -1;
        int np_row = -1;
        int my_col = -1;
        int my_row = -1;
        blacs::Cblacs_gridinfo(ictxt, &np_row, &np_col, &my_row, &my_col);
        EXPECT_EQ(comm.id2D(), std::make_pair(my_row, my_col));
        EXPECT_EQ(comm.size2D(), std::make_pair(np_row, np_col));
      }

      if (size1 != size2) {
        comm::Communicator2DGrid comm(mpi_comm, size2, size1, order);
        auto ictxt = comm.blacsContext();
        int np_col = -1;
        int np_row = -1;
        int my_col = -1;
        int my_row = -1;
        blacs::Cblacs_gridinfo(ictxt, &np_row, &np_col, &my_row, &my_col);
        EXPECT_EQ(comm.id2D(), std::make_pair(my_row, my_col));
        EXPECT_EQ(comm.size2D(), std::make_pair(np_row, np_col));
      }
    }
  }
}
#endif

int main(int argc, char** argv) {
  int provided = 0;
  MPI_Init_thread(nullptr, nullptr, MPI_THREAD_SERIALIZED, &provided);
  if (MPI_THREAD_SERIALIZED != provided) {
    std::cout << "MPI_Init_thread error!" << std::endl;
    MPI_Finalize();
    return 1;
  }

  auto id_size_pair = commInfo();
  int rank = id_size_pair.first;
  int size = id_size_pair.second;
  if (size != 6) {
    std::cout << "This test need 6 MPI ranks (" << size << " provided)!" << std::endl;
    MPI_Finalize();
    return 1;
  }

  mpi_comms.push_back(MPI_COMM_WORLD);

  MPI_Comm tmp;
  MPI_Comm_split(MPI_COMM_WORLD, (rank == 1 || rank == 5 ? 1 : 2), size - rank, &tmp);
  mpi_comms.push_back(tmp);
  MPI_Comm_split(MPI_COMM_WORLD, rank % 2, rank, &tmp);
  mpi_comms.push_back(tmp);

  ::testing::InitGoogleTest(&argc, argv);

  auto ret = RUN_ALL_TESTS();

  for (auto& comm : mpi_comms)
    if (comm != MPI_COMM_WORLD)
      MPI_Comm_free(&comm);

  MPI_Finalize();
  return ret;
}
