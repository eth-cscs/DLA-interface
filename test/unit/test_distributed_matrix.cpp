#include "distributed_matrix.h"

#include <memory>
#include <stdexcept>
#include "communicator_grid.h"
#include "communicator_manager.h"
#include "local_matrix.h"
#include "gtest/gtest.h"
#include "mpi_listener.h"
#include "null_stream.h"
#include "util_local_matrix.h"
#include "util_distributed_matrix.h"
#include "ref_scalapack_tools.h"
#include "tile_matrix_tools.h"

using namespace dla_interface;
using namespace testing;

std::vector<comm::Communicator2DGrid*> comms;

TEST(DistributedMatrixTest, DefaultConstructor) {
  using Type = double;

  DistributedMatrix<Type> mat;
  EXPECT_EQ(std::make_pair(0, 0), mat.size());
  EXPECT_EQ(std::make_pair(0, 0), mat.localSize());
  EXPECT_EQ(std::make_pair(1, 1), mat.blockSize());
  EXPECT_EQ(Global2DIndex(0, 0), mat.baseIndex());
  EXPECT_EQ(Local2DIndex(0, 0), mat.localBaseIndex());
  EXPECT_LE(1, mat.leadingDimension());
  EXPECT_EQ(1, mat.leadingNumberOfBlocks());
  EXPECT_EQ(nullptr, mat.ptr());
}

TEST(DistributedMatrixTest, Constructor1) {
  using Type = double;

  int m = 19;
  int n = 17;
  int mb = 2;
  int nb = 3;

  for (auto comm_ptr : comms) {
    auto id_2D = comm_ptr->id2D();
    int local_m = reference_scalapack_tools::numroc(m, mb, comm_ptr->id2D().first, 0,
                                                    comm_ptr->size2D().first);
    int local_n = reference_scalapack_tools::numroc(n, nb, comm_ptr->id2D().second, 0,
                                                    comm_ptr->size2D().second);
    for (auto dist : {scalapack_dist, tile_dist}) {
      DistributedMatrix<Type> mat(m, n, mb, nb, *comm_ptr, dist);
      EXPECT_EQ(std::make_pair(m, n), mat.size());
      EXPECT_EQ(std::make_pair(local_m, local_n), mat.localSize());
      EXPECT_EQ(std::make_pair(mb, nb), mat.blockSize());
      EXPECT_EQ(Global2DIndex(0, 0), mat.baseIndex());
      EXPECT_EQ(Local2DIndex(0, 0), mat.localBaseIndex());
      EXPECT_EQ(id_2D, mat.commGrid().id2D());
      EXPECT_EQ(dist, mat.distribution());
      int ld_min;
      int lnrbl;
      switch (dist) {
        case scalapack_dist:
          ld_min = local_m;
          lnrbl = 1;
          break;
        case tile_dist:
          ld_min = mb;
          lnrbl = (local_m + mb - 1) / mb;
          break;
        default:
          ASSERT_TRUE(false);
      }
      EXPECT_LE(ld_min, mat.leadingDimension());
      EXPECT_EQ(lnrbl, mat.leadingNumberOfBlocks());
      if (local_m * local_n == 0)
        EXPECT_EQ(nullptr, mat.ptr());
      else
        EXPECT_NE(nullptr, mat.ptr());

      EXPECT_THROW(DistributedMatrix<Type>(-1, n, mb, nb, *comm_ptr, dist), std::invalid_argument);
      EXPECT_THROW(DistributedMatrix<Type>(m, -1, mb, nb, *comm_ptr, dist), std::invalid_argument);
      EXPECT_THROW(DistributedMatrix<Type>(m, n, 0, nb, *comm_ptr, dist), std::invalid_argument);
      EXPECT_THROW(DistributedMatrix<Type>(m, n, mb, 0, *comm_ptr, dist), std::invalid_argument);
    }
  }

  // Empty matrix constructors
  for (auto comm_ptr : comms) {
    auto id_2D = comm_ptr->id2D();
    for (auto dist : {scalapack_dist, tile_dist}) {
      DistributedMatrix<Type> mat1(0, n, mb, nb, *comm_ptr, dist);
      EXPECT_EQ(std::make_pair(0, 0), mat1.size());
      EXPECT_EQ(std::make_pair(0, 0), mat1.localSize());
      EXPECT_EQ(std::make_pair(mb, nb), mat1.blockSize());
      EXPECT_EQ(Global2DIndex(0, 0), mat1.baseIndex());
      EXPECT_EQ(Local2DIndex(0, 0), mat1.localBaseIndex());
      EXPECT_LE(1, mat1.leadingDimension());
      EXPECT_EQ(1, mat1.leadingNumberOfBlocks());
      EXPECT_EQ(nullptr, mat1.ptr());
      EXPECT_EQ(id_2D, mat1.commGrid().id2D());
      EXPECT_EQ(dist, mat1.distribution());

      DistributedMatrix<Type> mat2(m, 0, mb, nb, *comm_ptr, dist);
      EXPECT_EQ(std::make_pair(0, 0), mat2.size());
      EXPECT_EQ(std::make_pair(0, 0), mat2.localSize());
      EXPECT_EQ(std::make_pair(mb, nb), mat2.blockSize());
      EXPECT_EQ(Global2DIndex(0, 0), mat2.baseIndex());
      EXPECT_EQ(Local2DIndex(0, 0), mat2.localBaseIndex());
      EXPECT_LE(1, mat2.leadingDimension());
      EXPECT_EQ(1, mat2.leadingNumberOfBlocks());
      EXPECT_EQ(nullptr, mat2.ptr());
      EXPECT_EQ(id_2D, mat2.commGrid().id2D());
      EXPECT_EQ(dist, mat2.distribution());
    }
  }
}

TEST(DistributedMatrixTest, Constructor2) {
  using Type = double;

  int m = 19;
  int n = 17;
  int mb = 2;
  int nb = 3;

  for (auto comm_ptr : comms) {
    auto id_2D = comm_ptr->id2D();
    int local_m = reference_scalapack_tools::numroc(m, mb, comm_ptr->id2D().first, 0,
                                                    comm_ptr->size2D().first);
    int local_n = reference_scalapack_tools::numroc(n, nb, comm_ptr->id2D().second, 0,
                                                    comm_ptr->size2D().second);
    for (auto dist : {scalapack_dist, tile_dist}) {
      int ld_min;
      int lnrbl;
      switch (dist) {
        case scalapack_dist:
          ld_min = std::max(1, local_m);
          lnrbl = 1;
          break;
        case tile_dist:
          ld_min = mb;
          lnrbl = (local_m + mb - 1) / mb;
          break;
        default:
          ASSERT_TRUE(false);
      }
      for (int ld : {ld_min, ld_min + 3}) {
        DistributedMatrix<Type> mat(m, n, mb, nb, ld, *comm_ptr, dist);
        EXPECT_EQ(std::make_pair(m, n), mat.size());
        EXPECT_EQ(std::make_pair(local_m, local_n), mat.localSize());
        EXPECT_EQ(std::make_pair(mb, nb), mat.blockSize());
        EXPECT_EQ(Global2DIndex(0, 0), mat.baseIndex());
        EXPECT_EQ(Local2DIndex(0, 0), mat.localBaseIndex());
        EXPECT_EQ(ld, mat.leadingDimension());
        EXPECT_EQ(lnrbl, mat.leadingNumberOfBlocks());
        if (local_m * local_n == 0)
          EXPECT_EQ(nullptr, mat.ptr());
        else
          EXPECT_NE(nullptr, mat.ptr());
        EXPECT_EQ(id_2D, mat.commGrid().id2D());
        EXPECT_EQ(dist, mat.distribution());
      }

      int ld = ld_min + 3;
      EXPECT_THROW(DistributedMatrix<Type>(m, n, mb, nb, ld_min - 1, *comm_ptr, dist),
                   std::invalid_argument);
      EXPECT_THROW(DistributedMatrix<Type>(-1, n, mb, nb, ld, *comm_ptr, dist),
                   std::invalid_argument);
      EXPECT_THROW(DistributedMatrix<Type>(m, -1, mb, nb, ld, *comm_ptr, dist),
                   std::invalid_argument);
      EXPECT_THROW(DistributedMatrix<Type>(m, n, 0, nb, ld, *comm_ptr, dist), std::invalid_argument);
      EXPECT_THROW(DistributedMatrix<Type>(m, n, mb, 0, ld, *comm_ptr, dist), std::invalid_argument);
    }
  }

  // Empty matrix constructors
  for (auto comm_ptr : comms) {
    auto id_2D = comm_ptr->id2D();
    int ld = 10;
    for (auto dist : {scalapack_dist, tile_dist}) {
      DistributedMatrix<Type> mat1(0, n, mb, nb, ld, *comm_ptr, dist);
      EXPECT_EQ(std::make_pair(0, 0), mat1.size());
      EXPECT_EQ(std::make_pair(0, 0), mat1.localSize());
      EXPECT_EQ(std::make_pair(mb, nb), mat1.blockSize());
      EXPECT_EQ(Global2DIndex(0, 0), mat1.baseIndex());
      EXPECT_EQ(Local2DIndex(0, 0), mat1.localBaseIndex());
      EXPECT_EQ(ld, mat1.leadingDimension());
      EXPECT_EQ(1, mat1.leadingNumberOfBlocks());
      EXPECT_EQ(nullptr, mat1.ptr());
      EXPECT_EQ(id_2D, mat1.commGrid().id2D());
      EXPECT_EQ(dist, mat1.distribution());

      DistributedMatrix<Type> mat2(m, 0, mb, nb, ld, *comm_ptr, dist);
      EXPECT_EQ(std::make_pair(0, 0), mat2.size());
      EXPECT_EQ(std::make_pair(0, 0), mat2.localSize());
      EXPECT_EQ(std::make_pair(mb, nb), mat2.blockSize());
      EXPECT_EQ(Global2DIndex(0, 0), mat2.baseIndex());
      EXPECT_EQ(Local2DIndex(0, 0), mat2.localBaseIndex());
      EXPECT_EQ(ld, mat2.leadingDimension());
      EXPECT_EQ(1, mat2.leadingNumberOfBlocks());
      EXPECT_EQ(nullptr, mat2.ptr());
      EXPECT_EQ(id_2D, mat2.commGrid().id2D());
      EXPECT_EQ(dist, mat2.distribution());

      EXPECT_THROW(DistributedMatrix<Type>(0, n, mb, nb, 0, *comm_ptr, dist), std::invalid_argument);
    }

    DistributedMatrix<Type> mat3(m, 0, mb, nb, 1, *comm_ptr, scalapack_dist);
    EXPECT_EQ(std::make_pair(0, 0), mat3.size());
    EXPECT_EQ(std::make_pair(0, 0), mat3.localSize());
    EXPECT_EQ(std::make_pair(mb, nb), mat3.blockSize());
    EXPECT_EQ(Global2DIndex(0, 0), mat3.baseIndex());
    EXPECT_EQ(Local2DIndex(0, 0), mat3.localBaseIndex());
    EXPECT_EQ(1, mat3.leadingDimension());
    EXPECT_EQ(1, mat3.leadingNumberOfBlocks());
    EXPECT_EQ(nullptr, mat3.ptr());
    EXPECT_EQ(id_2D, mat3.commGrid().id2D());
    EXPECT_EQ(scalapack_dist, mat3.distribution());
  }
}

TEST(DistributedMatrixTest, Constructor3) {
  using Type = double;

  int m = 19;
  int n = 17;
  int mb = 2;
  int nb = 3;

  auto el_val = [](int i, int j) { return j + .001 * i; };
  auto el_val2 = [](int i, int j) { return -j + .05 * i; };
  for (auto comm_ptr : comms) {
    auto id_2D = comm_ptr->id2D();
    int local_m = reference_scalapack_tools::numroc(m, mb, comm_ptr->id2D().first, 0,
                                                    comm_ptr->size2D().first);
    int local_n = reference_scalapack_tools::numroc(n, nb, comm_ptr->id2D().second, 0,
                                                    comm_ptr->size2D().second);
    std::size_t max_len = 0;
    for (auto dist1 : {scalapack_dist, tile_dist}) {
      for (auto dist2 : {scalapack_dist, tile_dist}) {
        int ld_min_dist1;
        int lnrbl_dist1;
        switch (dist1) {
          case scalapack_dist:
            ld_min_dist1 = std::max(1, local_m);
            lnrbl_dist1 = 1;
            break;
          case tile_dist:
            ld_min_dist1 = mb;
            lnrbl_dist1 = (local_m + mb - 1) / mb;
            break;
          default:
            ASSERT_TRUE(false);
        }
        int ld_min_dist2;
        int lnrbl_dist2;
        switch (dist2) {
          case scalapack_dist:
            ld_min_dist2 = std::max(1, local_m);
            lnrbl_dist2 = 1;
            break;
          case tile_dist:
            ld_min_dist2 = mb;
            lnrbl_dist2 = (local_m + mb - 1) / mb;
            break;
          default:
            ASSERT_TRUE(false);
        }
        for (int ld : {ld_min_dist1, ld_min_dist1 + 3}) {
          DistributedMatrix<Type> mat_ref(m, n, mb, nb, ld, *comm_ptr, dist1);
          fillDistributedMatrix(mat_ref, el_val);
          ASSERT_EQ(Global2DIndex(0, 0), mat_ref.baseIndex());
          std::size_t len = mat_ref.localStorageSize();
          max_len = std::max(max_len, len);
          {
            DistributedMatrix<Type> mat(m, n, mb, nb, *comm_ptr, dist2, mat_ref.ptr(), len,
                                        mat_ref.leadingDimension(), mat_ref.leadingNumberOfBlocks(),
                                        mat_ref.distribution());
            EXPECT_EQ(std::make_pair(m, n), mat.size());
            EXPECT_EQ(std::make_pair(local_m, local_n), mat.localSize());
            EXPECT_EQ(std::make_pair(mb, nb), mat.blockSize());
            EXPECT_EQ(Global2DIndex(0, 0), mat.baseIndex());
            EXPECT_EQ(Local2DIndex(0, 0), mat.localBaseIndex());
            if (dist1 == dist2) {
              EXPECT_EQ(mat_ref.leadingDimension(), mat.leadingDimension());
              EXPECT_EQ(mat_ref.leadingNumberOfBlocks(), mat.leadingNumberOfBlocks());
            }
            else {
              EXPECT_LE(ld_min_dist2, mat.leadingDimension());
              EXPECT_EQ(lnrbl_dist2, mat.leadingNumberOfBlocks());
            }
            EXPECT_TRUE(checkDistributedMatrix(mat, el_val));
            if (local_m * local_n == 0)
              EXPECT_EQ(nullptr, mat.ptr());
            else if (dist1 == dist2)
              EXPECT_EQ(mat_ref.ptr(), mat.ptr());
            else
              EXPECT_NE(mat_ref.ptr(), mat.ptr());
            EXPECT_EQ(id_2D, mat.commGrid().id2D());
            EXPECT_EQ(dist2, mat.distribution());

            fillDistributedMatrix(mat, el_val2);
          }
          EXPECT_EQ(dist1, mat_ref.distribution());
          EXPECT_TRUE(checkDistributedMatrix(mat_ref, el_val2));
        }

        int ld = ld_min_dist1 + 3;
        std::vector<Type> buf(max_len);
        if (local_m * local_n == 0)
          EXPECT_NO_THROW(DistributedMatrix<Type>(m, n, mb, nb, *comm_ptr, dist2, nullptr, 0, ld,
                                                  lnrbl_dist1, dist1));
        else
          EXPECT_THROW(DistributedMatrix<Type>(m, n, mb, nb, *comm_ptr, dist2, nullptr, max_len, ld,
                                               lnrbl_dist1, dist1),
                       std::invalid_argument);
        EXPECT_THROW(DistributedMatrix<Type>(m, n, mb, nb, *comm_ptr, dist2, &buf[0], max_len,
                                             ld_min_dist1 - 1, lnrbl_dist1, dist1),
                     std::invalid_argument);
        if (dist1 == scalapack_dist)  // ld_nr_bl is ignored by the constructor
          EXPECT_NO_THROW(DistributedMatrix<Type>(m, n, mb, nb, *comm_ptr, dist2, &buf[0], max_len,
                                                  ld_min_dist1, lnrbl_dist1 - 1, dist1));
        else
          EXPECT_THROW(DistributedMatrix<Type>(m, n, mb, nb, *comm_ptr, dist2, &buf[0], max_len,
                                               ld_min_dist1, lnrbl_dist1 - 1, dist1),
                       std::invalid_argument);
        EXPECT_THROW(DistributedMatrix<Type>(-1, n, mb, nb, *comm_ptr, dist2, &buf[0], max_len, ld,
                                             lnrbl_dist1, dist1),
                     std::invalid_argument);
        EXPECT_THROW(DistributedMatrix<Type>(m, -1, mb, nb, *comm_ptr, dist2, &buf[0], max_len, ld,
                                             lnrbl_dist1, dist1),
                     std::invalid_argument);
        EXPECT_THROW(DistributedMatrix<Type>(m, n, 0, nb, *comm_ptr, dist2, &buf[0], max_len, ld,
                                             lnrbl_dist1, dist1),
                     std::invalid_argument);
        EXPECT_THROW(DistributedMatrix<Type>(m, n, mb, 0, *comm_ptr, dist2, &buf[0], max_len, ld,
                                             lnrbl_dist1, dist1),
                     std::invalid_argument);
      }
    }
  }

  // Empty matrix constructors
  for (auto comm_ptr : comms) {
    auto id_2D = comm_ptr->id2D();
    int ld = 10;
    std::vector<Type> buf(13);
    for (Type* ptr : {(Type*)nullptr, &buf[0]}) {
      for (auto dist1 : {scalapack_dist, tile_dist}) {
        for (auto dist2 : {scalapack_dist, tile_dist}) {
          int ldmin_dist1 = dist1 == tile_dist ? mb : 1;
          int ldmin_dist2 = dist2 == tile_dist ? mb : 1;
          DistributedMatrix<Type> mat1(0, n, mb, nb, *comm_ptr, dist2, ptr, 0, ld, 1, dist1);
          EXPECT_EQ(std::make_pair(0, 0), mat1.size());
          EXPECT_EQ(std::make_pair(0, 0), mat1.localSize());
          EXPECT_EQ(std::make_pair(mb, nb), mat1.blockSize());
          EXPECT_EQ(Global2DIndex(0, 0), mat1.baseIndex());
          EXPECT_EQ(Local2DIndex(0, 0), mat1.localBaseIndex());
          if (dist1 == dist2)
            EXPECT_EQ(ld, mat1.leadingDimension());
          else
            EXPECT_LE(ldmin_dist2, mat1.leadingDimension());
          EXPECT_EQ(1, mat1.leadingNumberOfBlocks());
          EXPECT_EQ(nullptr, mat1.ptr());
          EXPECT_EQ(id_2D, mat1.commGrid().id2D());
          EXPECT_EQ(dist2, mat1.distribution());

          DistributedMatrix<Type> mat2(m, 0, mb, nb, *comm_ptr, dist2, ptr, 0, ld, 1, dist1);
          EXPECT_EQ(std::make_pair(0, 0), mat2.size());
          EXPECT_EQ(std::make_pair(0, 0), mat2.localSize());
          EXPECT_EQ(std::make_pair(mb, nb), mat2.blockSize());
          EXPECT_EQ(Global2DIndex(0, 0), mat2.baseIndex());
          EXPECT_EQ(Local2DIndex(0, 0), mat2.localBaseIndex());
          if (dist1 == dist2)
            EXPECT_EQ(ld, mat2.leadingDimension());
          else
            EXPECT_LE(ldmin_dist2, mat2.leadingDimension());
          EXPECT_EQ(1, mat2.leadingNumberOfBlocks());
          EXPECT_EQ(nullptr, mat2.ptr());
          EXPECT_EQ(id_2D, mat2.commGrid().id2D());
          EXPECT_EQ(dist2, mat2.distribution());

          DistributedMatrix<Type> mat3(m, 0, mb, nb, *comm_ptr, dist2, ptr, 0, ldmin_dist1, 1, dist1);
          EXPECT_EQ(std::make_pair(0, 0), mat3.size());
          EXPECT_EQ(std::make_pair(0, 0), mat3.localSize());
          EXPECT_EQ(std::make_pair(mb, nb), mat3.blockSize());
          EXPECT_EQ(Global2DIndex(0, 0), mat3.baseIndex());
          EXPECT_EQ(Local2DIndex(0, 0), mat3.localBaseIndex());
          if (dist1 == dist2)
            EXPECT_EQ(ldmin_dist1, mat3.leadingDimension());
          else
            EXPECT_LE(ldmin_dist2, mat3.leadingDimension());
          EXPECT_EQ(1, mat3.leadingNumberOfBlocks());
          EXPECT_EQ(nullptr, mat3.ptr());
          EXPECT_EQ(id_2D, mat3.commGrid().id2D());
          EXPECT_EQ(dist2, mat3.distribution());

          EXPECT_THROW(DistributedMatrix<Type>(0, n, mb, nb, *comm_ptr, dist2, ptr, 0,
                                               ldmin_dist1 - 1, 1, dist1),
                       std::invalid_argument);
        }
      }
    }
  }
}

TEST(DistributedMatrixTest, Constructor4) {
  using Type = double;

  int m = 19;
  int n = 17;
  int mb = 2;
  int nb = 3;

  auto el_val = [](int i, int j) { return j + .001 * i; };
  auto el_val2 = [](int i, int j) { return -j + .05 * i; };
  for (auto comm_ptr : comms) {
    auto id_2D = comm_ptr->id2D();
    int local_m = reference_scalapack_tools::numroc(m, mb, comm_ptr->id2D().first, 0,
                                                    comm_ptr->size2D().first);
    int local_n = reference_scalapack_tools::numroc(n, nb, comm_ptr->id2D().second, 0,
                                                    comm_ptr->size2D().second);
    for (auto dist1 : {scalapack_dist, tile_dist}) {
      for (auto dist2 : {scalapack_dist, tile_dist}) {
        int ld_min_dist1;
        switch (dist1) {
          case scalapack_dist:
            ld_min_dist1 = std::max(1, local_m);
            break;
          case tile_dist:
            ld_min_dist1 = mb;
            break;
          default:
            ASSERT_TRUE(false);
        }
        int ld_min_dist2;
        int lnrbl_dist2;
        switch (dist2) {
          case scalapack_dist:
            ld_min_dist2 = std::max(1, local_m);
            lnrbl_dist2 = 1;
            break;
          case tile_dist:
            ld_min_dist2 = mb;
            lnrbl_dist2 = (local_m + mb - 1) / mb;
            break;
          default:
            ASSERT_TRUE(false);
        }
        for (int ld : {ld_min_dist1, ld_min_dist1 + 3}) {
          auto mat1_ptr =
              std::make_shared<DistributedMatrix<Type>>(m, n, mb, nb, ld, *comm_ptr, dist1);
          auto mat2_ptr = DistributedMatrix<Type>(m, n, mb, nb, ld, *comm_ptr, dist1)
                              .subMatrix(m - mb - 1, n - 3 * nb + 1, mb + 1, 3 * nb - 1);
          for (auto mat_ref_ptr : {mat1_ptr, mat2_ptr}) {
            auto& mat_ref = *mat_ref_ptr;
            fillDistributedMatrix(mat_ref, el_val);
            {
              DistributedMatrix<Type> mat(dist2, mat_ref);
              EXPECT_EQ(mat_ref.size(), mat.size());
              EXPECT_EQ(mat_ref.localSize(), mat.localSize());
              EXPECT_EQ(mat_ref.blockSize(), mat.blockSize());
              EXPECT_EQ(mat_ref.baseIndex(), mat.baseIndex());
              EXPECT_EQ(mat_ref.localBaseIndex(), mat.localBaseIndex());
              EXPECT_LE(ld_min_dist2, mat.leadingDimension());
              if (dist1 == dist2)
                EXPECT_EQ(ld, mat.leadingDimension());
              EXPECT_EQ(lnrbl_dist2, mat.leadingNumberOfBlocks());
              EXPECT_TRUE(checkDistributedMatrix(mat, el_val));
              if (local_m * local_n == 0)
                EXPECT_EQ(nullptr, mat.ptr());
              else {
                if (dist1 == dist2)
                  EXPECT_EQ(mat_ref.ptr(), mat.ptr());
                else
                  EXPECT_NE(mat_ref.ptr(), mat.ptr());
              }
              EXPECT_EQ(id_2D, mat.commGrid().id2D());
              EXPECT_EQ(dist2, mat.distribution());

              fillDistributedMatrix(mat, el_val2);
            }
            EXPECT_EQ(dist1, mat_ref.distribution());
            EXPECT_TRUE(checkDistributedMatrix(mat_ref, el_val2));
          }
        }
      }
    }
  }
}

TEST(DistributedMatrixTest, CopyConvert1) {
  using Type = double;

  int m = 19;
  int n = 17;
  int mb = 2;
  int nb = 3;

  auto el_val = [](int i, int j) { return j + .001 * i; };
  for (auto comm_ptr : comms) {
    auto id_2D = comm_ptr->id2D();
    int local_m = reference_scalapack_tools::numroc(m, mb, comm_ptr->id2D().first, 0,
                                                    comm_ptr->size2D().first);
    int local_n = reference_scalapack_tools::numroc(n, nb, comm_ptr->id2D().second, 0,
                                                    comm_ptr->size2D().second);
    std::size_t max_len = 0;
    for (auto dist1 : {scalapack_dist, tile_dist}) {
      for (auto dist2 : {scalapack_dist, tile_dist}) {
        int ld_min_dist1;
        int lnrbl_dist1;
        switch (dist1) {
          case scalapack_dist:
            ld_min_dist1 = std::max(1, local_m);
            lnrbl_dist1 = 1;
            break;
          case tile_dist:
            ld_min_dist1 = mb;
            lnrbl_dist1 = (local_m + mb - 1) / mb;
            break;
          default:
            ASSERT_TRUE(false);
        }
        int ld_min_dist2;
        int lnrbl_dist2;
        switch (dist2) {
          case scalapack_dist:
            ld_min_dist2 = std::max(1, local_m);
            lnrbl_dist2 = 1;
            break;
          case tile_dist:
            ld_min_dist2 = mb;
            lnrbl_dist2 = (local_m + mb - 1) / mb;
            break;
          default:
            ASSERT_TRUE(false);
        }
        for (int ld : {ld_min_dist1, ld_min_dist1 + 3}) {
          DistributedMatrix<Type> mat_ref(m, n, mb, nb, ld, *comm_ptr, dist1);
          fillDistributedMatrix(mat_ref, el_val);
          ASSERT_EQ(Global2DIndex(0, 0), mat_ref.baseIndex());
          std::size_t len = mat_ref.localStorageSize();
          max_len = std::max(max_len, len);
          {
            const auto mat_ref_const_ptr = mat_ref.ptr();
            std::shared_ptr<const DistributedMatrix<Type>> mat = DistributedMatrix<Type>::convertConst(
                m, n, mb, nb, *comm_ptr, dist2, mat_ref_const_ptr, len, mat_ref.leadingDimension(),
                mat_ref.leadingNumberOfBlocks(), mat_ref.distribution());
            EXPECT_EQ(std::make_pair(m, n), mat->size());
            EXPECT_EQ(std::make_pair(local_m, local_n), mat->localSize());
            EXPECT_EQ(std::make_pair(mb, nb), mat->blockSize());
            EXPECT_EQ(Global2DIndex(0, 0), mat->baseIndex());
            EXPECT_EQ(Local2DIndex(0, 0), mat->localBaseIndex());
            if (dist1 == dist2) {
              EXPECT_EQ(mat_ref.leadingDimension(), mat->leadingDimension());
              EXPECT_EQ(mat_ref.leadingNumberOfBlocks(), mat->leadingNumberOfBlocks());
            }
            else {
              EXPECT_LE(ld_min_dist2, mat->leadingDimension());
              EXPECT_EQ(lnrbl_dist2, mat->leadingNumberOfBlocks());
            }
            EXPECT_TRUE(checkDistributedMatrix(*mat, el_val));
            if (local_m * local_n == 0)
              EXPECT_EQ(nullptr, mat->ptr());
            else if (dist1 == dist2)
              EXPECT_EQ(mat_ref_const_ptr, mat->ptr());
            else
              EXPECT_NE(mat_ref_const_ptr, mat->ptr());
            EXPECT_EQ(id_2D, mat->commGrid().id2D());
            EXPECT_EQ(dist2, mat->distribution());
          }
        }

        int ld = ld_min_dist1 + 3;
        std::vector<Type> buf(max_len);
        const Type* buf_ptr = &buf[0];
        if (local_m * local_n == 0)
          EXPECT_NO_THROW(DistributedMatrix<Type>::convertConst(
              m, n, mb, nb, *comm_ptr, dist2, nullptr, 0, ld, lnrbl_dist1, dist1));
        else
          EXPECT_THROW(DistributedMatrix<Type>::convertConst(m, n, mb, nb, *comm_ptr, dist2, nullptr,
                                                             max_len, ld, lnrbl_dist1, dist1),
                       std::invalid_argument);
        EXPECT_THROW(DistributedMatrix<Type>::convertConst(m, n, mb, nb, *comm_ptr, dist2, buf_ptr,
                                                           max_len, ld_min_dist1 - 1, lnrbl_dist1,
                                                           dist1),
                     std::invalid_argument);
        if (dist1 == scalapack_dist)  // ld_nr_bl is ignored by the constructor
          EXPECT_NO_THROW(DistributedMatrix<Type>::convertConst(m, n, mb, nb, *comm_ptr, dist2,
                                                                buf_ptr, max_len, ld_min_dist1,
                                                                lnrbl_dist1 - 1, dist1));
        else
          EXPECT_THROW(DistributedMatrix<Type>::convertConst(m, n, mb, nb, *comm_ptr, dist2,
                                                             buf_ptr, max_len, ld_min_dist1,
                                                             lnrbl_dist1 - 1, dist1),
                       std::invalid_argument);
        EXPECT_THROW(DistributedMatrix<Type>::convertConst(-1, n, mb, nb, *comm_ptr, dist2, buf_ptr,
                                                           max_len, ld, lnrbl_dist1, dist1),
                     std::invalid_argument);
        EXPECT_THROW(DistributedMatrix<Type>::convertConst(m, -1, mb, nb, *comm_ptr, dist2, buf_ptr,
                                                           max_len, ld, lnrbl_dist1, dist1),
                     std::invalid_argument);
        EXPECT_THROW(DistributedMatrix<Type>::convertConst(m, n, 0, nb, *comm_ptr, dist2, buf_ptr,
                                                           max_len, ld, lnrbl_dist1, dist1),
                     std::invalid_argument);
        EXPECT_THROW(DistributedMatrix<Type>::convertConst(m, n, mb, 0, *comm_ptr, dist2, buf_ptr,
                                                           max_len, ld, lnrbl_dist1, dist1),
                     std::invalid_argument);
      }
    }
  }

  // Empty matrix conversion
  for (auto comm_ptr : comms) {
    auto id_2D = comm_ptr->id2D();
    int ld = 10;
    std::vector<Type> buf(13);
    for (const Type* ptr : {(const Type*)nullptr, (const Type*)&buf[0]}) {
      for (auto dist1 : {scalapack_dist, tile_dist}) {
        for (auto dist2 : {scalapack_dist, tile_dist}) {
          int ldmin_dist1 = dist1 == tile_dist ? mb : 1;
          int ldmin_dist2 = dist2 == tile_dist ? mb : 1;
          auto mat1_ptr = DistributedMatrix<Type>::convertConst(0, n, mb, nb, *comm_ptr, dist2, ptr,
                                                                0, ld, 1, dist1);
          EXPECT_EQ(std::make_pair(0, 0), mat1_ptr->size());
          EXPECT_EQ(std::make_pair(0, 0), mat1_ptr->localSize());
          EXPECT_EQ(std::make_pair(mb, nb), mat1_ptr->blockSize());
          EXPECT_EQ(Global2DIndex(0, 0), mat1_ptr->baseIndex());
          EXPECT_EQ(Local2DIndex(0, 0), mat1_ptr->localBaseIndex());
          if (dist1 == dist2)
            EXPECT_EQ(ld, mat1_ptr->leadingDimension());
          else
            EXPECT_LE(ldmin_dist2, mat1_ptr->leadingDimension());
          EXPECT_EQ(1, mat1_ptr->leadingNumberOfBlocks());
          EXPECT_EQ(nullptr, mat1_ptr->ptr());
          EXPECT_EQ(id_2D, mat1_ptr->commGrid().id2D());
          EXPECT_EQ(dist2, mat1_ptr->distribution());

          auto mat2_ptr = DistributedMatrix<Type>::convertConst(m, 0, mb, nb, *comm_ptr, dist2, ptr,
                                                                0, ld, 1, dist1);
          EXPECT_EQ(std::make_pair(0, 0), mat2_ptr->size());
          EXPECT_EQ(std::make_pair(0, 0), mat2_ptr->localSize());
          EXPECT_EQ(std::make_pair(mb, nb), mat2_ptr->blockSize());
          EXPECT_EQ(Global2DIndex(0, 0), mat2_ptr->baseIndex());
          EXPECT_EQ(Local2DIndex(0, 0), mat2_ptr->localBaseIndex());
          if (dist1 == dist2)
            EXPECT_EQ(ld, mat2_ptr->leadingDimension());
          else
            EXPECT_LE(ldmin_dist2, mat2_ptr->leadingDimension());
          EXPECT_EQ(1, mat2_ptr->leadingNumberOfBlocks());
          EXPECT_EQ(nullptr, mat2_ptr->ptr());
          EXPECT_EQ(id_2D, mat2_ptr->commGrid().id2D());
          EXPECT_EQ(dist2, mat2_ptr->distribution());

          auto mat3_ptr = DistributedMatrix<Type>::convertConst(m, 0, mb, nb, *comm_ptr, dist2, ptr,
                                                                0, ldmin_dist1, 1, dist1);
          EXPECT_EQ(std::make_pair(0, 0), mat3_ptr->size());
          EXPECT_EQ(std::make_pair(0, 0), mat3_ptr->localSize());
          EXPECT_EQ(std::make_pair(mb, nb), mat3_ptr->blockSize());
          EXPECT_EQ(Global2DIndex(0, 0), mat3_ptr->baseIndex());
          EXPECT_EQ(Local2DIndex(0, 0), mat3_ptr->localBaseIndex());
          if (dist1 == dist2)
            EXPECT_EQ(ldmin_dist1, mat3_ptr->leadingDimension());
          else
            EXPECT_LE(ldmin_dist2, mat3_ptr->leadingDimension());
          EXPECT_EQ(1, mat3_ptr->leadingNumberOfBlocks());
          EXPECT_EQ(nullptr, mat3_ptr->ptr());
          EXPECT_EQ(id_2D, mat3_ptr->commGrid().id2D());
          EXPECT_EQ(dist2, mat3_ptr->distribution());

          EXPECT_THROW(DistributedMatrix<Type>::convertConst(0, n, mb, nb, *comm_ptr, dist2, ptr, 0,
                                                             ldmin_dist1 - 1, 1, dist1),
                       std::invalid_argument);
        }
      }
    }
  }
}

TEST(DistributedMatrixTest, ConvertConst2) {
  using Type = double;

  int m = 19;
  int n = 17;
  int mb = 2;
  int nb = 3;

  auto el_val = [](int i, int j) { return j + .001 * i; };
  for (auto comm_ptr : comms) {
    auto id_2D = comm_ptr->id2D();
    int local_m = reference_scalapack_tools::numroc(m, mb, comm_ptr->id2D().first, 0,
                                                    comm_ptr->size2D().first);
    int local_n = reference_scalapack_tools::numroc(n, nb, comm_ptr->id2D().second, 0,
                                                    comm_ptr->size2D().second);
    for (auto dist1 : {scalapack_dist, tile_dist}) {
      for (auto dist2 : {scalapack_dist, tile_dist}) {
        int ld_min_dist1;
        switch (dist1) {
          case scalapack_dist:
            ld_min_dist1 = std::max(1, local_m);
            break;
          case tile_dist:
            ld_min_dist1 = mb;
            break;
          default:
            ASSERT_TRUE(false);
        }
        int ld_min_dist2;
        int lnrbl_dist2;
        switch (dist2) {
          case scalapack_dist:
            ld_min_dist2 = std::max(1, local_m);
            lnrbl_dist2 = 1;
            break;
          case tile_dist:
            ld_min_dist2 = mb;
            lnrbl_dist2 = (local_m + mb - 1) / mb;
            break;
          default:
            ASSERT_TRUE(false);
        }
        for (int ld : {ld_min_dist1, ld_min_dist1 + 3}) {
          auto mat1_ptr =
              std::make_shared<DistributedMatrix<Type>>(m, n, mb, nb, ld, *comm_ptr, dist1);
          auto mat2_ptr = DistributedMatrix<Type>(m, n, mb, nb, ld, *comm_ptr, dist1)
                              .subMatrix(m - mb - 1, n - 3 * nb + 1, mb + 1, 3 * nb - 1);
          for (auto mat_ref_ptr : {mat1_ptr, mat2_ptr}) {
            fillDistributedMatrix(*mat_ref_ptr, el_val);
            const auto& mat_ref = *mat_ref_ptr;
            {
              auto mat_ptr = mat_ref.convertConst(dist2);
              const auto& mat = *mat_ptr;
              EXPECT_EQ(mat_ref.size(), mat.size());
              EXPECT_EQ(mat_ref.localSize(), mat.localSize());
              EXPECT_EQ(mat_ref.blockSize(), mat.blockSize());
              EXPECT_EQ(mat_ref.baseIndex(), mat.baseIndex());
              EXPECT_EQ(mat_ref.localBaseIndex(), mat.localBaseIndex());
              EXPECT_LE(ld_min_dist2, mat.leadingDimension());
              if (dist1 == dist2)
                EXPECT_EQ(ld, mat.leadingDimension());
              EXPECT_EQ(lnrbl_dist2, mat.leadingNumberOfBlocks());
              EXPECT_TRUE(checkDistributedMatrix(mat, el_val));
              if (local_m * local_n == 0)
                EXPECT_EQ(nullptr, mat.ptr());
              else {
                if (dist1 == dist2)
                  EXPECT_EQ(mat_ref.ptr(), mat.ptr());
                else
                  EXPECT_NE(mat_ref.ptr(), mat.ptr());
              }
              EXPECT_EQ(id_2D, mat.commGrid().id2D());
              EXPECT_EQ(dist2, mat.distribution());
            }
          }
        }
      }
    }
  }
}

TEST(DistributedMatrixTest, IndexConversion) {
  using Type = double;

  int m = 19;
  int n = 17;
  int mb = 2;
  int nb = 3;

  for (auto comm_ptr : comms) {
    auto id_2D = comm_ptr->id2D();
    auto comm_size = comm_ptr->size2D();
    int local_m = reference_scalapack_tools::numroc(m, mb, comm_ptr->id2D().first, 0,
                                                    comm_ptr->size2D().first);
    for (auto dist : {scalapack_dist, tile_dist}) {
      int ld_min;
      switch (dist) {
        case scalapack_dist:
          ld_min = std::max(1, local_m);
          break;
        case tile_dist:
          ld_min = mb;
          break;
        default:
          ASSERT_TRUE(false);
      }

      DistributedMatrix<Type> mat(m, n, mb, nb, ld_min + 3, *comm_ptr, dist);
      const auto& const_mat(mat);
      auto size = mat.size();
      auto block_size = mat.blockSize();
      int ld = mat.leadingDimension();

      auto g_base_index = mat.baseIndex();
      auto l_base_index = mat.localBaseIndex();

      for (int j = 0; j < size.second; ++j) {
        for (int i = 0; i < size.first; ++i) {
          int scalapack_gi = i + g_base_index.row + 1;
          int scalapack_gj = j + g_base_index.col + 1;
          std::pair<int, int> rank;

          int scalapack_li;
          int scalapack_lj;
          reference_scalapack_tools::infog1l(scalapack_gi, block_size.first, comm_size.first,
                                             id_2D.first, 0, scalapack_li, rank.first);
          reference_scalapack_tools::infog1l(scalapack_gj, block_size.second, comm_size.second,
                                             id_2D.second, 0, scalapack_lj, rank.second);

          Local2DIndex l_index_ref(scalapack_li - l_base_index.row - 1,
                                   scalapack_lj - l_base_index.col - 1);
          Global2DIndex g_index_ref(i, j);

          EXPECT_EQ(l_index_ref, mat.getLocal2DIndex(i, j));
          EXPECT_EQ(l_index_ref, mat.getLocal2DIndex(g_index_ref));
          EXPECT_EQ(l_index_ref, const_mat.getLocal2DIndex(i, j));
          EXPECT_EQ(l_index_ref, const_mat.getLocal2DIndex(g_index_ref));

          EXPECT_EQ(rank, mat.getRankId2D(i, j));
          EXPECT_EQ(rank, mat.getRankId2D(g_index_ref));
          EXPECT_EQ(rank, const_mat.getRankId2D(i, j));
          EXPECT_EQ(rank, const_mat.getRankId2D(g_index_ref));

          if (id_2D == rank) {
            EXPECT_EQ(g_index_ref, mat.getGlobal2DIndex(l_index_ref));
            EXPECT_EQ(g_index_ref, const_mat.getGlobal2DIndex(l_index_ref));

            std::size_t storage_index_ref = 0;
            switch (dist) {
              case scalapack_dist:
                storage_index_ref = scalapack_li - 1 + ld * (scalapack_lj - 1);
                break;
              case tile_dist:
                storage_index_ref = tile_matrix_tools::arrayIndex(
                    l_index_ref, l_base_index, block_size, mat.leadingNumberOfBlocks(), ld);
                break;
              default:
                ASSERT_TRUE(false);
            }
            EXPECT_EQ(storage_index_ref, mat.getLocalStorageIndex(i, j));
            EXPECT_EQ(storage_index_ref, mat.getLocalStorageIndex(g_index_ref));
            ASSERT_EQ(storage_index_ref, mat.getLocalStorageIndex(l_index_ref));
            EXPECT_EQ(storage_index_ref, const_mat.getLocalStorageIndex(i, j));
            EXPECT_EQ(storage_index_ref, const_mat.getLocalStorageIndex(g_index_ref));
            EXPECT_EQ(storage_index_ref, const_mat.getLocalStorageIndex(l_index_ref));
          }
          else {
            EXPECT_THROW(mat.getLocalStorageIndex(i, j), std::invalid_argument);
            EXPECT_THROW(mat.getLocalStorageIndex(g_index_ref), std::invalid_argument);
            EXPECT_THROW(const_mat.getLocalStorageIndex(i, j), std::invalid_argument);
            EXPECT_THROW(const_mat.getLocalStorageIndex(g_index_ref), std::invalid_argument);
          }
        }
      }
    }
  }
}

TEST(DistributedMatrixTest, ElementAndPointer) {
  using Type = double;

  std::vector<DistributedMatrix<Type>> vmat;
  for (auto comm_ptr : comms) {
    for (auto dist : {scalapack_dist, tile_dist}) {
      vmat.emplace_back(17, 13, 2, 2, *comm_ptr, dist);
      vmat.emplace_back(13, 7, 3, 3, 16, *comm_ptr, dist);
      vmat.emplace_back(5, 2, 4, 4, *comm_ptr, dist);
      vmat.emplace_back(0, 7, 2, 2, *comm_ptr, dist);
      vmat.emplace_back(7, 0, 2, 2, 2, *comm_ptr, dist);
      vmat.emplace_back();
    }
  }

  for (auto& mat : vmat) {
    auto el_val = [](int i, int j) { return j + .001 * i; };
    const auto& const_mat = mat;
    int m = mat.size().first;
    int n = mat.size().second;

    for (int j = 0; j < mat.size().second; ++j) {
      for (int i = 0; i < mat.size().first; ++i) {
        if (mat.commGrid().id2D() == mat.getRankId2D(i, j)) {
          mat(i, j) = el_val(i, j);
        }
      }
    }
    auto ptr = mat.ptr();
    auto const_ptr = const_mat.ptr();
    ASSERT_EQ(ptr, const_ptr);

    for (int j = 0; j < n; ++j) {
      for (int i = 0; i < m; ++i) {
        Global2DIndex g_index(i, j);
        if (mat.commGrid().id2D() == mat.getRankId2D(i, j)) {
          Local2DIndex l_index = mat.getLocal2DIndex(i, j);
          EXPECT_EQ(el_val(i, j), mat(i, j));
          EXPECT_EQ(ptr + mat.getLocalStorageIndex(i, j), mat.ptr(i, j));
          EXPECT_EQ(mat.ptr(i, j), mat.ptr(g_index));
          EXPECT_EQ(mat.ptr(i, j), mat.ptr(l_index));
          EXPECT_EQ(mat.ptr(i, j), &mat(i, j));
          EXPECT_EQ(mat.ptr(i, j), &mat(g_index));
          EXPECT_EQ(mat.ptr(i, j), &mat(l_index));
          EXPECT_EQ(el_val(i, j), const_mat(i, j));
          EXPECT_EQ(const_ptr + const_mat.getLocalStorageIndex(i, j), const_mat.ptr(i, j));
          EXPECT_EQ(const_mat.ptr(i, j), const_mat.ptr(g_index));
          EXPECT_EQ(const_mat.ptr(i, j), const_mat.ptr(l_index));
          EXPECT_EQ(const_mat.ptr(i, j), &const_mat(i, j));
          EXPECT_EQ(const_mat.ptr(i, j), &const_mat(g_index));
          EXPECT_EQ(const_mat.ptr(i, j), &const_mat(l_index));
        }
        else {
          EXPECT_THROW(mat(i, j), std::invalid_argument);
          EXPECT_THROW(mat(g_index), std::invalid_argument);
          EXPECT_THROW(const_mat(i, j), std::invalid_argument);
          EXPECT_THROW(const_mat(g_index), std::invalid_argument);
        }
      }
    }
    int ml = mat.localSize().first;
    int nl = mat.localSize().second;
    EXPECT_THROW(mat(-1, 0), std::invalid_argument);
    EXPECT_THROW(mat(0, -1), std::invalid_argument);
    EXPECT_THROW(mat(m + 1, 0), std::invalid_argument);
    EXPECT_THROW(mat(0, n + 1), std::invalid_argument);
    EXPECT_THROW(mat.ptr(-1, 0), std::invalid_argument);
    EXPECT_THROW(mat.ptr(0, -1), std::invalid_argument);
    EXPECT_THROW(mat.ptr(m + 1, 0), std::invalid_argument);
    EXPECT_THROW(mat.ptr(0, n + 1), std::invalid_argument);

    EXPECT_THROW(mat(Global2DIndex(-1, 0)), std::invalid_argument);
    EXPECT_THROW(mat(Global2DIndex(0, -1)), std::invalid_argument);
    EXPECT_THROW(mat(Global2DIndex(m + 1, 0)), std::invalid_argument);
    EXPECT_THROW(mat(Global2DIndex(0, n + 1)), std::invalid_argument);
    EXPECT_THROW(mat.ptr(Global2DIndex(-1, 0)), std::invalid_argument);
    EXPECT_THROW(mat.ptr(Global2DIndex(0, -1)), std::invalid_argument);
    EXPECT_THROW(mat.ptr(Global2DIndex(m + 1, 0)), std::invalid_argument);
    EXPECT_THROW(mat.ptr(Global2DIndex(0, n + 1)), std::invalid_argument);

    EXPECT_THROW(mat(Local2DIndex(-1, 0)), std::invalid_argument);
    EXPECT_THROW(mat(Local2DIndex(0, -1)), std::invalid_argument);
    EXPECT_THROW(mat(Local2DIndex(ml + 1, 0)), std::invalid_argument);
    EXPECT_THROW(mat(Local2DIndex(0, nl + 1)), std::invalid_argument);
    EXPECT_THROW(mat.ptr(Local2DIndex(-1, 0)), std::invalid_argument);
    EXPECT_THROW(mat.ptr(Local2DIndex(0, -1)), std::invalid_argument);
    EXPECT_THROW(mat.ptr(Local2DIndex(ml + 1, 0)), std::invalid_argument);
    EXPECT_THROW(mat.ptr(Local2DIndex(0, nl + 1)), std::invalid_argument);

    EXPECT_THROW(const_mat(-1, 0), std::invalid_argument);
    EXPECT_THROW(const_mat(0, -1), std::invalid_argument);
    EXPECT_THROW(const_mat(m + 1, 0), std::invalid_argument);
    EXPECT_THROW(const_mat(0, n + 1), std::invalid_argument);
    EXPECT_THROW(const_mat.ptr(-1, 0), std::invalid_argument);
    EXPECT_THROW(const_mat.ptr(0, -1), std::invalid_argument);
    EXPECT_THROW(const_mat.ptr(m + 1, 0), std::invalid_argument);
    EXPECT_THROW(const_mat.ptr(0, n + 1), std::invalid_argument);

    EXPECT_THROW(const_mat(Global2DIndex(-1, 0)), std::invalid_argument);
    EXPECT_THROW(const_mat(Global2DIndex(0, -1)), std::invalid_argument);
    EXPECT_THROW(const_mat(Global2DIndex(m + 1, 0)), std::invalid_argument);
    EXPECT_THROW(const_mat(Global2DIndex(0, n + 1)), std::invalid_argument);
    EXPECT_THROW(const_mat.ptr(Global2DIndex(-1, 0)), std::invalid_argument);
    EXPECT_THROW(const_mat.ptr(Global2DIndex(0, -1)), std::invalid_argument);
    EXPECT_THROW(const_mat.ptr(Global2DIndex(m + 1, 0)), std::invalid_argument);
    EXPECT_THROW(const_mat.ptr(Global2DIndex(0, n + 1)), std::invalid_argument);

    EXPECT_THROW(const_mat(Local2DIndex(-1, 0)), std::invalid_argument);
    EXPECT_THROW(const_mat(Local2DIndex(0, -1)), std::invalid_argument);
    EXPECT_THROW(const_mat(Local2DIndex(ml + 1, 0)), std::invalid_argument);
    EXPECT_THROW(const_mat(Local2DIndex(0, nl + 1)), std::invalid_argument);
    EXPECT_THROW(const_mat.ptr(Local2DIndex(-1, 0)), std::invalid_argument);
    EXPECT_THROW(const_mat.ptr(Local2DIndex(0, -1)), std::invalid_argument);
    EXPECT_THROW(const_mat.ptr(Local2DIndex(ml + 1, 0)), std::invalid_argument);
    EXPECT_THROW(const_mat.ptr(Local2DIndex(0, nl + 1)), std::invalid_argument);
  }
}

TEST(DistributedMatrixTest, CopyConstructor) {
  using Type = double;

  std::vector<DistributedMatrix<Type>> vmat;
  for (auto comm_ptr : comms) {
    for (auto dist : {scalapack_dist, tile_dist}) {
      vmat.emplace_back(17, 13, 2, 2, *comm_ptr, dist);
      vmat.emplace_back(13, 7, 3, 3, 16, *comm_ptr, dist);
      vmat.emplace_back(5, 2, 4, 4, *comm_ptr, dist);
      vmat.emplace_back(0, 7, 2, 2, *comm_ptr, dist);
      vmat.emplace_back(7, 0, 2, 2, 2, *comm_ptr, dist);
    }
  }

  auto el_val = [](int i, int j) { return j + .001 * i; };
  auto el_val2 = [](int i, int j) { return -j + .05 * i; };
  for (auto& mat : vmat) {
    const auto& const_mat(mat);
    fillDistributedMatrix(mat, el_val);
    auto size = mat.size();
    auto local_size = mat.localSize();
    auto distribution = mat.distribution();
    auto& comm_grid = mat.commGrid();

    DistributedMatrix<Type> mat2(mat);
    DistributedMatrix<Type> mat3(const_mat);

    EXPECT_EQ(size, mat.size());
    EXPECT_EQ(local_size, mat.localSize());
    EXPECT_EQ(distribution, mat.distribution());
    EXPECT_EQ(&comm_grid, &mat.commGrid());
    EXPECT_TRUE(checkDistributedMatrix(mat, el_val));

    EXPECT_EQ(size, mat2.size());
    EXPECT_EQ(local_size, mat2.localSize());
    EXPECT_EQ(distribution, mat2.distribution());
    EXPECT_EQ(&comm_grid, &mat2.commGrid());
    EXPECT_TRUE(checkDistributedMatrixSame(mat, mat2));

    EXPECT_EQ(size, mat3.size());
    EXPECT_EQ(local_size, mat3.localSize());
    EXPECT_EQ(distribution, mat3.distribution());
    EXPECT_EQ(&comm_grid, &mat3.commGrid());
    EXPECT_TRUE(checkDistributedMatrixSame(mat, mat3));

    // Change the elements of mat and check that the elements of mat2 and mat3 are not changed.
    fillDistributedMatrix(mat, el_val2);

    EXPECT_TRUE(checkDistributedMatrix(mat2, el_val));
    EXPECT_TRUE(checkDistributedMatrix(mat3, el_val));
  }
}

TEST(DistributedMatrixTest, AssignementOperator) {
  using Type = double;

  std::vector<DistributedMatrix<Type>> vmat;
  for (auto comm_ptr : comms) {
    for (auto dist : {scalapack_dist, tile_dist}) {
      vmat.emplace_back(17, 13, 2, 2, *comm_ptr, dist);
      vmat.emplace_back(13, 7, 3, 3, 16, *comm_ptr, dist);
      vmat.emplace_back(5, 2, 4, 4, *comm_ptr, dist);
      vmat.emplace_back(0, 7, 2, 2, *comm_ptr, dist);
      vmat.emplace_back(7, 0, 2, 2, 2, *comm_ptr, dist);
    }
  }

  auto el_val = [](int i, int j) { return j + .001 * i; };
  auto el_val2 = [](int i, int j) { return -j + .05 * i; };
  for (auto& mat : vmat) {
    const auto& const_mat(mat);
    fillDistributedMatrix(mat, el_val);
    auto size = mat.size();
    auto local_size = mat.localSize();
    auto distribution = mat.distribution();
    auto& comm_grid = mat.commGrid();

    DistributedMatrix<Type> mat2(7, 5, 1, 1, comm_grid, distribution);
    DistributedMatrix<Type> mat3(7, 5, 1, 1, comm_grid, distribution);

    auto oldptr_mat2 = mat2.ptr();
    auto oldptr_mat3 = mat3.ptr();

    mat2 = mat;
    mat3 = const_mat;

    EXPECT_EQ(size, mat.size());
    EXPECT_EQ(local_size, mat.localSize());
    EXPECT_EQ(distribution, mat.distribution());
    EXPECT_EQ(&comm_grid, &mat.commGrid());
    EXPECT_TRUE(checkDistributedMatrix(mat, el_val));

    EXPECT_EQ(size, mat2.size());
    EXPECT_EQ(local_size, mat2.localSize());
    EXPECT_EQ(distribution, mat2.distribution());
    EXPECT_EQ(&comm_grid, &mat2.commGrid());
    EXPECT_NE(oldptr_mat2, mat2.ptr());
    EXPECT_TRUE(checkDistributedMatrixSame(mat, mat2));

    EXPECT_EQ(size, mat3.size());
    EXPECT_EQ(local_size, mat3.localSize());
    EXPECT_EQ(distribution, mat3.distribution());
    EXPECT_EQ(&comm_grid, &mat3.commGrid());
    EXPECT_NE(oldptr_mat3, mat3.ptr());
    EXPECT_TRUE(checkDistributedMatrixSame(mat, mat3));

    // Change the elements of mat and check that the elements of mat2 and mat3 are not changed.
    fillDistributedMatrix(mat, el_val2);

    EXPECT_TRUE(checkDistributedMatrix(mat2, el_val));
    EXPECT_TRUE(checkDistributedMatrix(mat3, el_val));

    // Check that a reference to the object is returned.
    DistributedMatrix<Type> mat4;
    auto pmat4 = &mat4;
    EXPECT_EQ(pmat4, &(mat4 = mat));
  }
}

TEST(DistributedMatrixTest, MoveConstructor) {
  using Type = double;

  std::vector<DistributedMatrix<Type>> vmat;
  for (auto comm_ptr : comms) {
    for (auto dist : {scalapack_dist, tile_dist}) {
      vmat.emplace_back(17, 13, 2, 2, *comm_ptr, dist);
      vmat.emplace_back(13, 7, 3, 3, 16, *comm_ptr, dist);
      vmat.emplace_back(5, 2, 4, 4, *comm_ptr, dist);
      vmat.emplace_back(0, 7, 2, 2, *comm_ptr, dist);
      vmat.emplace_back(7, 0, 2, 2, 2, *comm_ptr, dist);
    }
  }

  auto el_val = [](int i, int j) { return j + .001 * i; };
  auto el_val2 = [](int i, int j) { return -j + .05 * i; };
  for (auto& mat : vmat) {
    fillDistributedMatrix(mat, el_val);
    auto size = mat.size();
    auto block_size = mat.blockSize();
    auto local_size = mat.localSize();
    auto ptr = mat.ptr();
    auto distribution = mat.distribution();
    auto& comm_grid = mat.commGrid();

    DistributedMatrix<Type> mat2(std::move(mat));

    EXPECT_EQ(size, mat2.size());
    EXPECT_EQ(local_size, mat2.localSize());
    EXPECT_EQ(distribution, mat2.distribution());
    EXPECT_EQ(&comm_grid, &mat2.commGrid());
    EXPECT_EQ(ptr, mat2.ptr());
    EXPECT_TRUE(checkDistributedMatrix(mat2, el_val));

    // Re-assign mat and check that mat2 do not change.
    mat = DistributedMatrix<Type>(size.first, size.second, block_size.first, block_size.second,
                                  comm_grid, distribution);
    fillDistributedMatrix(mat, el_val2);

    EXPECT_TRUE(checkDistributedMatrix(mat2, el_val));
  }
}

TEST(DistributedMatrixTest, MoveAssignementOperator) {
  using Type = double;

  std::vector<DistributedMatrix<Type>> vmat;
  for (auto comm_ptr : comms) {
    for (auto dist : {scalapack_dist, tile_dist}) {
      vmat.emplace_back(17, 13, 2, 2, *comm_ptr, dist);
      vmat.emplace_back(13, 7, 3, 3, 16, *comm_ptr, dist);
      vmat.emplace_back(5, 2, 4, 4, *comm_ptr, dist);
      vmat.emplace_back(0, 7, 2, 2, *comm_ptr, dist);
      vmat.emplace_back(7, 0, 2, 2, 2, *comm_ptr, dist);
    }
  }

  auto el_val = [](int i, int j) { return j + .001 * i; };
  auto el_val2 = [](int i, int j) { return -j + .05 * i; };
  for (auto& mat : vmat) {
    fillDistributedMatrix(mat, el_val);
    auto size = mat.size();
    auto block_size = mat.blockSize();
    auto local_size = mat.localSize();
    auto ptr = mat.ptr();
    auto distribution = mat.distribution();
    auto& comm_grid = mat.commGrid();

    DistributedMatrix<Type> mat2(7, 5, 1, 1, comm_grid, distribution);

    auto oldptr_mat2 = mat2.ptr();

    mat2 = std::move(mat);

    EXPECT_EQ(size, mat2.size());
    EXPECT_EQ(local_size, mat2.localSize());
    EXPECT_EQ(distribution, mat2.distribution());
    EXPECT_EQ(&comm_grid, &mat2.commGrid());
    EXPECT_EQ(ptr, mat2.ptr());
    EXPECT_NE(oldptr_mat2, mat2.ptr());
    EXPECT_TRUE(checkDistributedMatrix(mat2, el_val));

    // Re-assign mat and check that mat2 do not change.
    mat = DistributedMatrix<Type>(size.first, size.second, block_size.first, block_size.second,
                                  comm_grid, distribution);
    fillDistributedMatrix(mat, el_val2);

    EXPECT_TRUE(checkDistributedMatrix(mat2, el_val));

    // Check that a reference to the object is returned.
    DistributedMatrix<Type> mat3;
    auto pmat3 = &mat3;
    EXPECT_EQ(pmat3, &(mat3 = std::move(mat)));
  }
}

TEST(DistributedMatrixTest, SubMatrix) {
  using Type = double;

  std::vector<DistributedMatrix<Type>> vmat;
  for (auto comm_ptr : comms) {
    for (auto dist : {scalapack_dist, tile_dist}) {
      vmat.emplace_back(13, 8, 2, 3, *comm_ptr, dist);
      vmat.emplace_back(11, 7, 3, 2, 16, *comm_ptr, dist);
      vmat.emplace_back(5, 2, 4, 4, *comm_ptr, dist);
      vmat.emplace_back(7, 0, 2, 2, 2, *comm_ptr, dist);
    }
  }

  auto el_val = [](int i, int j) { return j + .001 * i; };
  auto el_val2 = [](int i, int j) { return -j + .05 * i; };
  for (auto& mat : vmat) {
    int m = mat.size().first;
    int n = mat.size().second;
    int mb = mat.blockSize().first;
    int nb = mat.blockSize().second;
    auto& comm = mat.commGrid();
    const auto& const_mat = mat;
    for (int si = 0; si < m; ++si) {
      for (int sj = 0; sj < n; ++sj) {
        for (int sm = 0; sm <= m - si; ++sm) {
          for (int sn = 0; sn <= n - sj; ++sn) {
            fillDistributedMatrix(mat, el_val);
            int local_si =
                reference_scalapack_tools::numroc(si, mb, comm.id2D().first, 0, comm.size2D().first);
            int local_sj = reference_scalapack_tools::numroc(sj, nb, comm.id2D().second, 0,
                                                             comm.size2D().second);
            int local_sm = reference_scalapack_tools::numroc(si + sm, mb, comm.id2D().first, 0,
                                                             comm.size2D().first) -
                           local_si;
            int local_sn = reference_scalapack_tools::numroc(sj + sn, nb, comm.id2D().second, 0,
                                                             comm.size2D().second) -
                           local_sj;

            auto psubmat = mat.subMatrix(sm, sn, si, sj);
            auto pconst_submat = const_mat.subMatrix(sm, sn, si, sj);

            auto size = std::make_pair(sm, sn);
            Global2DIndex base_index(si, sj);
            auto local_size = std::make_pair(local_sm, local_sn);
            Local2DIndex local_base_index(local_si, local_sj);
            if (sm == 0 || sn == 0) {
              size = std::make_pair(0, 0);
              base_index = Global2DIndex(0, 0);
            }
            if (local_sm == 0 || local_sn == 0) {
              EXPECT_EQ(nullptr, psubmat->ptr());
              EXPECT_EQ(nullptr, pconst_submat->ptr());
              local_size = std::make_pair(0, 0);
              local_base_index = Local2DIndex(0, 0);
            }
            EXPECT_EQ(size, psubmat->size());
            EXPECT_EQ(size, pconst_submat->size());
            EXPECT_EQ(base_index, psubmat->baseIndex());
            EXPECT_EQ(base_index, pconst_submat->baseIndex());
            EXPECT_EQ(local_size, psubmat->localSize());
            EXPECT_EQ(local_size, pconst_submat->localSize());
            EXPECT_EQ(local_base_index, psubmat->localBaseIndex());
            EXPECT_EQ(local_base_index, pconst_submat->localBaseIndex());
            EXPECT_TRUE(checkIfDistributedSubMatrix(mat, size, std::make_pair(si, sj), *psubmat));
            EXPECT_EQ(mat.distribution(), psubmat->distribution());
            EXPECT_TRUE(
                checkIfDistributedSubMatrix(mat, size, std::make_pair(si, sj), *pconst_submat));
            EXPECT_EQ(mat.distribution(), pconst_submat->distribution());
            auto sub_el_val = [&el_val, si, sj](int i, int j) { return el_val(i + si, j + sj); };
            EXPECT_TRUE(checkDistributedMatrix(*psubmat, sub_el_val));
            EXPECT_TRUE(checkDistributedMatrix(*pconst_submat, sub_el_val));

            fillDistributedMatrix(mat, el_val2);
            auto sub_el_val2 = [&el_val2, si, sj](int i, int j) { return el_val2(i + si, j + sj); };
            EXPECT_TRUE(checkDistributedMatrix(*psubmat, sub_el_val2));
            EXPECT_TRUE(checkDistributedMatrix(*pconst_submat, sub_el_val2));
          }
        }
      }
    }

    EXPECT_THROW(mat.subMatrix(m + 1, n, 0, 0), std::invalid_argument);
    EXPECT_THROW(mat.subMatrix(-1, n, 1, 0), std::invalid_argument);
    EXPECT_THROW(mat.subMatrix(2, n, m - 1, 0), std::invalid_argument);
    EXPECT_THROW(mat.subMatrix(1, n, -1, 0), std::invalid_argument);
    EXPECT_THROW(mat.subMatrix(m, n + 1, 0, 0), std::invalid_argument);
    EXPECT_THROW(mat.subMatrix(m, -1, 0, 1), std::invalid_argument);
    EXPECT_THROW(mat.subMatrix(m, 2, 0, n - 1), std::invalid_argument);
    EXPECT_THROW(mat.subMatrix(m, 1, 0, -1), std::invalid_argument);
  }
  vmat.clear();

  auto& comm_grid = *(comms[0]);
  auto comm_size = comm_grid.size2D();
  int mb = 2;
  int nb = 3;
  NullStream nullout;

  int sm = mb * comm_size.first;
  int sn = nb * comm_size.second;
  // 4 block per rank.
  DistributedMatrix<Type> mat(2 * sm, 2 * sn, mb, nb, comm_grid, scalapack_dist);
  // 1 block per rank.
  auto psubmat = mat.subMatrix(sm, sn, 0, 0);
  EXPECT_TRUE(
      checkIfDistributedSubMatrix(mat, std::make_pair(sm, sn), std::make_pair(0, 0), *psubmat));
  // Wrong size.
  EXPECT_FALSE(checkIfDistributedSubMatrix(mat, std::make_pair(sm - 1, sn), std::make_pair(0, 0),
                                           *psubmat, nullout));
  EXPECT_FALSE(checkIfDistributedSubMatrix(mat, std::make_pair(sm, sn - 1), std::make_pair(0, 0),
                                           *psubmat, nullout));
  // Wrong position.
  EXPECT_FALSE(checkIfDistributedSubMatrix(mat, std::make_pair(sm, sn), std::make_pair(mb, nb),
                                           *psubmat, nullout));
  EXPECT_FALSE(checkIfDistributedSubMatrix(mat, std::make_pair(sm, sn), std::make_pair(sm, sn),
                                           *psubmat, nullout));

  // Check that the submatrix doesn't change if the original matrix is re-assigned.
  int m = 17;
  int n = 13;
  DistributedMatrix<Type> mat1(m, n, mb, nb, comm_grid, scalapack_dist);
  fillDistributedMatrix(mat1, el_val);
  auto psubmat1 = mat1.subMatrix(m - 2, n - 2, 0, 0);
  EXPECT_TRUE(checkDistributedMatrix(*psubmat1, el_val));

  mat1 = DistributedMatrix<Type>(m, n, mb, nb, comm_grid, scalapack_dist);
  fillDistributedMatrix(mat1, el_val2);

  EXPECT_TRUE(checkDistributedMatrix(*psubmat1, el_val));

  DistributedMatrix<Type> mat2(m, n, mb, nb, comm_grid, scalapack_dist);
  fillDistributedMatrix(mat2, el_val);
  auto psubmat2 = mat2.subMatrix(m - 2, n - 2, 0, 0);
  EXPECT_TRUE(checkDistributedMatrix(*psubmat2, el_val));

  DistributedMatrix<Type> tmp(m, n, mb, nb, comm_grid, scalapack_dist);
  mat2 = std::move(tmp);
  fillDistributedMatrix(mat2, el_val);

  EXPECT_TRUE(checkDistributedMatrix(*psubmat2, el_val));

  // Check that the submatrix is still valid after original matrix goes out of scope.
  std::shared_ptr<DistributedMatrix<Type>> psubmat3;
  {
    DistributedMatrix<Type> mat3(m, n, mb, nb, comm_grid, scalapack_dist);
    fillDistributedMatrix(mat3, el_val);
    psubmat3 = mat3.subMatrix(m - 2, n - 2, 0, 0);
    EXPECT_TRUE(checkDistributedMatrix(*psubmat3, el_val));
  }
  EXPECT_TRUE(checkDistributedMatrix(*psubmat3, el_val));
}

TEST(DistributedMatrixTest, LocalSubMatrix) {
  using Type = double;

  std::vector<DistributedMatrix<Type>> vmat;
  for (auto comm_ptr : comms) {
    for (auto dist : {scalapack_dist, tile_dist}) {
      vmat.emplace_back(13, 8, 2, 3, *comm_ptr, dist);
      vmat.emplace_back(11, 7, 3, 2, 16, *comm_ptr, dist);
      vmat.emplace_back(5, 2, 4, 4, *comm_ptr, dist);
      vmat.emplace_back(7, 0, 2, 2, 2, *comm_ptr, dist);
    }
  }

  auto el_val = [](int i, int j) { return j + .001 * i; };
  auto el_val2 = [](int i, int j) { return -j + .05 * i; };
  for (auto& mat : vmat) {
    int m = mat.size().first;
    int n = mat.size().second;
    int mb = mat.blockSize().first;
    int nb = mat.blockSize().second;
    const auto& const_mat = mat;
    for (int si = 0; si < m; ++si) {
      for (int sj = 0; sj < n; ++sj) {
        int max_m = std::min(m - si, mb - si % mb);
        int max_n = std::min(n - sj, nb - sj % nb);
        for (int sm = 0; sm < max_m; ++sm) {
          for (int sn = 0; sn < max_n; ++sn) {
            fillDistributedMatrix(mat, el_val);

            auto psubmat = mat.localSubMatrix(sm, sn, si, sj);
            auto pconst_submat = const_mat.localSubMatrix(sm, sn, si, sj);
            if (mat.commGrid().id2D() == mat.getRankId2D(si, sj) && sm > 0 && sn > 0) {
              EXPECT_EQ(std::make_pair(sm, sn), psubmat->size());
              EXPECT_EQ(mat.leadingDimension(), psubmat->leadingDimension());
              EXPECT_EQ(std::make_pair(sm, sn), pconst_submat->size());
              EXPECT_EQ(mat.leadingDimension(), pconst_submat->leadingDimension());

              EXPECT_EQ(mat.ptr(si, sj), psubmat->ptr());
              EXPECT_EQ(mat.ptr(si, sj), psubmat->ptr(0, 0));
              EXPECT_EQ(mat.ptr(si, sj), pconst_submat->ptr());
              EXPECT_EQ(mat.ptr(si, sj), pconst_submat->ptr(0, 0));
            }
            else {
              EXPECT_EQ(std::make_pair(0, 0), psubmat->size());
              EXPECT_EQ(mat.leadingDimension(), psubmat->leadingDimension());
              EXPECT_EQ(nullptr, psubmat->ptr());
              EXPECT_EQ(std::make_pair(0, 0), pconst_submat->size());
              EXPECT_EQ(mat.leadingDimension(), pconst_submat->leadingDimension());
              EXPECT_EQ(nullptr, pconst_submat->ptr());
            }
            auto sub_el_val = [&el_val, si, sj](int i, int j) { return el_val(i + si, j + sj); };
            EXPECT_TRUE(checkLocalMatrix(*psubmat, sub_el_val));
            EXPECT_TRUE(checkLocalMatrix(*pconst_submat, sub_el_val));

            fillDistributedMatrix(mat, el_val2);
            auto sub_el_val2 = [&el_val2, si, sj](int i, int j) { return el_val2(i + si, j + sj); };
            EXPECT_TRUE(checkLocalMatrix(*psubmat, sub_el_val2));
            EXPECT_TRUE(checkLocalMatrix(*pconst_submat, sub_el_val2));
          }
        }
        EXPECT_THROW(mat.localSubMatrix(max_m + 1, max_n, si, sj), std::invalid_argument);
        EXPECT_THROW(mat.localSubMatrix(-1, max_n, si + 1, sj), std::invalid_argument);
        EXPECT_THROW(mat.localSubMatrix(2, max_n, si + max_m - 1, sj), std::invalid_argument);
        EXPECT_THROW(mat.localSubMatrix(max_m, max_n + 1, si, sj), std::invalid_argument);
        EXPECT_THROW(mat.localSubMatrix(max_m, -1, si, sj + 1), std::invalid_argument);
        EXPECT_THROW(mat.localSubMatrix(max_m, 2, si, sj + max_n - 1), std::invalid_argument);
      }
    }
    EXPECT_THROW(mat.localSubMatrix(1, 1, -1, 0), std::invalid_argument);
    EXPECT_THROW(mat.localSubMatrix(1, 1, 0, -1), std::invalid_argument);
  }

  // Check that the submatrix doesn't change if the original matrix is re-assigned.
  int m = 17;
  int n = 13;
  int mb = 2;
  int nb = 3;
  auto& comm_grid = *(comms[0]);

  DistributedMatrix<Type> mat1(m, n, mb, nb, comm_grid, scalapack_dist);
  fillDistributedMatrix(mat1, el_val);
  auto psubmat1 = mat1.localSubMatrix(mb, nb, 0, 0);
  EXPECT_TRUE(checkLocalMatrix(*psubmat1, el_val));

  mat1 = DistributedMatrix<Type>(m, n, mb, nb, comm_grid, scalapack_dist);
  fillDistributedMatrix(mat1, el_val2);

  EXPECT_TRUE(checkLocalMatrix(*psubmat1, el_val));

  DistributedMatrix<Type> mat2(m, n, mb, nb, comm_grid, scalapack_dist);
  fillDistributedMatrix(mat2, el_val);
  auto psubmat2 = mat2.localSubMatrix(mb, nb, 0, 0);
  EXPECT_TRUE(checkLocalMatrix(*psubmat2, el_val));

  DistributedMatrix<Type> tmp(m, n, mb, nb, comm_grid, scalapack_dist);
  mat2 = std::move(tmp);
  fillDistributedMatrix(mat2, el_val);

  EXPECT_TRUE(checkLocalMatrix(*psubmat2, el_val));

  // Check that the submatrix is still valid after original matrix goes out of scope.
  std::shared_ptr<LocalMatrix<Type>> psubmat3;
  {
    DistributedMatrix<Type> mat3(m, n, mb, nb, comm_grid, scalapack_dist);
    fillDistributedMatrix(mat3, el_val);
    psubmat3 = mat3.localSubMatrix(mb, nb, 0, 0);
    EXPECT_TRUE(checkLocalMatrix(*psubmat3, el_val));
  }
  EXPECT_TRUE(checkLocalMatrix(*psubmat3, el_val));
}

TEST(DistributedMatrixTest, isSameMatrix) {
  using Type = double;

  auto& comm1 = *comms[0];
  auto& comm2 = *comms[1];

  int m = 17;
  int n = 13;
  int mb = 2;
  int nb = 3;
  int ld = 32;

  for (auto dist1 : {scalapack_dist, tile_dist}) {
    DistributedMatrix<Type> mat1(m, n, mb, nb, comm1, dist1);
    EXPECT_TRUE(mat1.isSameMatrix(mat1));

    // different memory
    DistributedMatrix<Type> mat2(m, n, mb, nb, comm1, dist1);
    EXPECT_FALSE(mat1.isSameMatrix(mat2));

    DistributedMatrix<Type> mat3(m, n, mb, nb, comm1, mat1.distribution(), mat1.ptr(),
                                 mat1.localStorageSize(), mat1.leadingDimension(),
                                 mat1.leadingNumberOfBlocks(), mat1.distribution());
    EXPECT_TRUE(mat1.isSameMatrix(mat3));

    // different size
    DistributedMatrix<Type> mat4(m - 1, n, mb, nb, comm1, mat1.distribution(), mat1.ptr(),
                                 mat1.localStorageSize(), mat1.leadingDimension(),
                                 mat1.leadingNumberOfBlocks(), mat1.distribution());
    EXPECT_FALSE(mat1.isSameMatrix(mat4));

    // different size
    DistributedMatrix<Type> mat5(m, n - 1, mb, nb, comm1, mat1.distribution(), mat1.ptr(),
                                 mat1.localStorageSize(), mat1.leadingDimension(),
                                 mat1.leadingNumberOfBlocks(), mat1.distribution());
    EXPECT_FALSE(mat1.isSameMatrix(mat5));
  }

  std::size_t len = (4 * (mb + 1) + ld) * (n + nb + 1);
  std::vector<Type> buf1(len);
  std::vector<Type> buf2(len);
  for (auto dist1 : {scalapack_dist, tile_dist}) {
    int ld1 = dist1 == scalapack_dist ? ld : mb + 1;
    int ld_nr_blocks1 = dist1 == scalapack_dist ? 1 : m / mb + 1;
    DistributedMatrix<Type> mat1(m, n, mb, nb, comm1, dist1, &buf1[0], len, ld1, ld_nr_blocks1,
                                 dist1);

    DistributedMatrix<Type> mat2(m, n, mb, nb, comm1, dist1, &buf1[0], len, ld1, ld_nr_blocks1,
                                 dist1);
    EXPECT_TRUE(mat1.isSameMatrix(mat2));

    // different block size
    int ld_nr_blocks = dist1 == scalapack_dist ? 1 : m / (mb + 1) + 1;
    DistributedMatrix<Type> mat3(m, n, mb + 1, nb, comm1, dist1, &buf1[0], len, ld1, ld_nr_blocks,
                                 dist1);
    EXPECT_FALSE(mat1.isSameMatrix(mat3));

    // different block size
    DistributedMatrix<Type> mat4(m, n, mb, nb + 1, comm1, dist1, &buf1[0], len, ld1, ld_nr_blocks1,
                                 dist1);
    EXPECT_FALSE(mat1.isSameMatrix(mat4));

    // different communicator
    DistributedMatrix<Type> mat5(m, n, mb, nb, comm2, mat1.distribution(), mat1.ptr(),
                                 mat1.localStorageSize(), mat1.leadingDimension(),
                                 mat1.leadingNumberOfBlocks(), mat1.distribution());
    EXPECT_FALSE(mat1.isSameMatrix(mat5));

    // different distribution
    for (auto dist2 : {scalapack_dist, tile_dist}) {
      int ld2 = dist2 == scalapack_dist ? ld : mb + 1;
      int ld_nr_blocks2 = dist2 == scalapack_dist ? 1 : m / mb + 1;
      DistributedMatrix<Type> mat(m, n, mb, nb, comm1, dist2, &buf1[0], len, ld2, ld_nr_blocks2,
                                  dist2);
      if (dist1 == dist2)
        EXPECT_TRUE(mat1.isSameMatrix(mat));
      else
        EXPECT_FALSE(mat1.isSameMatrix(mat));
    }

    // different memory
    DistributedMatrix<Type> mat6(m, n, mb, nb + 1, comm1, dist1, &buf2[0], len, ld1, ld_nr_blocks1,
                                 dist1);
    EXPECT_FALSE(mat1.isSameMatrix(mat6));

    // different leading dimension
    DistributedMatrix<Type> mat7(m, n, mb, nb, comm1, dist1, &buf2[0], len, ld1 + 1, ld_nr_blocks1,
                                 dist1);
    EXPECT_FALSE(mat1.isSameMatrix(mat7));

    // different leading number of blocks
    DistributedMatrix<Type> mat8(m, n, mb, nb, comm1, dist1, &buf2[0], len, ld1, ld_nr_blocks1 + 1,
                                 dist1);
    EXPECT_FALSE(mat1.isSameMatrix(mat8));
  }
}

TEST(DistributedMatrixTest, Copy) {
  // This test handle the case when the distribution is the same and when is different.
  using Type = double;
  int m = 17;
  int n = 13;
  int mb = 2;
  int nb = 3;

  auto el_val = [](int i, int j) { return j + .001 * i; };
  auto el_val2 = [](int i, int j) { return -j + .05 * i; };

  for (auto comm_ptr : comms) {
    for (auto dist1 : {scalapack_dist, tile_dist}) {
      DistributedMatrix<Type> mat1(m, n, mb, nb, *comm_ptr, dist1);
      auto pmat2 =
          DistributedMatrix<Type>(m + 1, n + 1, mb, nb, *comm_ptr, dist1).subMatrix(m, n, 1, 1);
      auto& mat2 = *pmat2;

      for (auto dist2 : {scalapack_dist, tile_dist}) {
        std::vector<DistributedMatrix<Type>> vmat;
        vmat.emplace_back(m, n, mb, nb, *comm_ptr, dist2);
        vmat.emplace_back(m, n, mb, nb, m + 2, *comm_ptr, dist2);

        for (auto& mat : vmat) {
          fillDistributedMatrix(mat1, el_val);
          fillDistributedMatrix(mat, el_val2);

          mat.copy(mat1);
          fillDistributedMatrix(mat1, el_val2);
          EXPECT_TRUE(checkDistributedMatrix(mat, el_val));

          EXPECT_THROW(mat.copy(mat2), std::invalid_argument);

          // Check that a reference to the object is returned.
          auto pmat = &mat;
          EXPECT_EQ(pmat, &(mat.copy(mat1)));
        }

        vmat.clear();
        vmat.emplace_back(m + 1, n + 1, mb, nb, *comm_ptr, dist2);
        vmat.emplace_back(m + 1, n + 1, mb, nb, m + 4, *comm_ptr, dist2);

        for (auto& full_mat : vmat) {
          auto pmat = full_mat.subMatrix(m, n, 1, 1);
          auto& mat = *pmat;

          EXPECT_THROW(mat.copy(mat1), std::invalid_argument);

          fillDistributedMatrix(mat2, el_val);
          fillDistributedMatrix(mat, el_val2);

          mat.copy(mat2);
          fillDistributedMatrix(mat2, el_val2);
          EXPECT_TRUE(checkDistributedMatrix(mat, el_val));
        }

        vmat.clear();
        vmat.emplace_back(m + 1, n, mb, nb, *comm_ptr, dist2);
        vmat.emplace_back(m, n + 1, mb, nb, *comm_ptr, dist2);
        vmat.emplace_back(m, n, mb + 1, nb, *comm_ptr, dist2);
        vmat.emplace_back(m, n, mb, nb + 1, *comm_ptr, dist2);

        for (auto& mat : vmat) {
          EXPECT_THROW(mat.copy(mat1), std::invalid_argument);
          EXPECT_THROW(mat1.copy(mat), std::invalid_argument);
        }
      }
    }
  }
}

#ifdef DLA_HAVE_SCALAPACK
TEST(DistributedMatrixTest, ScalapackDescriptor) {
  using Type = double;
  int m = 17;
  int n = 13;
  int mb = 3;
  int nb = 2;

  for (auto comm_ptr : comms) {
    DistributedMatrix<Type> mat(m, n, mb, nb, *comm_ptr, scalapack_dist);
    auto data_scalapack = mat.getScalapackDescription();
    EXPECT_EQ(mat.ptr(), std::get<0>(data_scalapack));
    EXPECT_EQ(1, std::get<1>(data_scalapack));
    EXPECT_EQ(1, std::get<2>(data_scalapack));
    EXPECT_EQ(1, std::get<3>(data_scalapack)[0]);
    EXPECT_EQ(comm_ptr->blacsContext(), std::get<3>(data_scalapack)[1]);
    EXPECT_EQ(m, std::get<3>(data_scalapack)[2]);
    EXPECT_EQ(n, std::get<3>(data_scalapack)[3]);
    EXPECT_EQ(mb, std::get<3>(data_scalapack)[4]);
    EXPECT_EQ(nb, std::get<3>(data_scalapack)[5]);
    EXPECT_EQ(0, std::get<3>(data_scalapack)[6]);
    EXPECT_EQ(0, std::get<3>(data_scalapack)[7]);
    EXPECT_EQ(mat.leadingDimension(), std::get<3>(data_scalapack)[8]);

    auto psubmat = mat.subMatrix(m - 1, n - 1, 1, 1);
    DistributedMatrix<Type>& submat = *psubmat;
    data_scalapack = submat.getScalapackDescription();
    EXPECT_EQ(mat.ptr(), std::get<0>(data_scalapack));
    EXPECT_EQ(2, std::get<1>(data_scalapack));
    EXPECT_EQ(2, std::get<2>(data_scalapack));
    EXPECT_EQ(1, std::get<3>(data_scalapack)[0]);
    EXPECT_EQ(comm_ptr->blacsContext(), std::get<3>(data_scalapack)[1]);
    EXPECT_EQ(m, std::get<3>(data_scalapack)[2]);
    EXPECT_EQ(n, std::get<3>(data_scalapack)[3]);
    EXPECT_EQ(mb, std::get<3>(data_scalapack)[4]);
    EXPECT_EQ(nb, std::get<3>(data_scalapack)[5]);
    EXPECT_EQ(0, std::get<3>(data_scalapack)[6]);
    EXPECT_EQ(0, std::get<3>(data_scalapack)[7]);
    EXPECT_EQ(mat.leadingDimension(), std::get<3>(data_scalapack)[8]);

    DistributedMatrix<Type> tile_mat(m, n, mb, nb, *comm_ptr, tile_dist);
    EXPECT_THROW(tile_mat.getScalapackDescription(), std::invalid_argument);
  }
}

TEST(DistributedMatrixTest, ConstructorScalapack) {
  using Type = double;
  int m = 11;
  int n = 7;
  int mb = 3;
  int nb = 2;

  auto el_val0 = [](int i, int j) { return 0; };
  auto el_val = [](int i, int j) { return j + .001 * i; };
  auto el_val2 = [](int i, int j) { return -j + .05 * i; };

  for (auto comm_ptr : comms) {
    DistributedMatrix<Type> mat_ref(m, n, mb, nb, *comm_ptr, scalapack_dist);
    int ld = mat_ref.leadingDimension();

    for (int si = 0; si < m; ++si) {
      for (int sj = 0; sj < n; ++sj) {
        for (int sm = 0; sm < m - si; ++sm) {
          for (int sn = 0; sn < n - sj; ++sn) {
            fillDistributedMatrix(mat_ref, el_val0);
            auto psubmat = mat_ref.subMatrix(sm, sn, si, sj);
            DistributedMatrix<Type>& submat_ref = *psubmat;
            auto data_scalapack = submat_ref.getScalapackDescription();
            for (auto dist : {scalapack_dist, tile_dist}) {
              fillDistributedMatrix(submat_ref, el_val);

              {
                DistributedMatrix<Type> mat(dist, sm, sn, std::get<0>(data_scalapack),
                                            std::get<1>(data_scalapack), std::get<2>(data_scalapack),
                                            &std::get<3>(data_scalapack)[0]);
                auto size = std::make_pair(sm, sn);
                Global2DIndex base_index(si, sj);
                if (sm == 0 || sn == 0) {
                  size = std::make_pair(0, 0);
                  base_index = Global2DIndex(0, 0);
                }
                EXPECT_EQ(size, mat.size());
                EXPECT_EQ(std::make_pair(mb, nb), mat.blockSize());
                EXPECT_EQ(base_index, mat.baseIndex());
                if (dist == scalapack_dist)
                  EXPECT_EQ(ld, mat.leadingDimension());
                else
                  EXPECT_EQ(mb, mat.leadingDimension());
                EXPECT_TRUE(checkDistributedMatrix(mat, el_val));
                EXPECT_EQ(dist, mat.distribution());

                fillDistributedMatrix(mat, el_val2);
              }
              EXPECT_TRUE(checkDistributedMatrix(submat_ref, el_val2));
            }
          }
        }
      }
    }

    auto data_scalapack = mat_ref.getScalapackDescription();
    auto ptr = std::get<0>(data_scalapack);
    auto desc = &std::get<3>(data_scalapack)[0];
    EXPECT_THROW(DistributedMatrix<Type> mat(scalapack_dist, -1, n, ptr, 1, 1, desc),
                 std::invalid_argument);
    EXPECT_THROW(DistributedMatrix<Type> mat(scalapack_dist, m, -1, ptr, 1, 1, desc),
                 std::invalid_argument);
    EXPECT_THROW(DistributedMatrix<Type> mat(scalapack_dist, m, n, ptr, 0, 1, desc),
                 std::invalid_argument);
    EXPECT_THROW(DistributedMatrix<Type> mat(scalapack_dist, m, n, ptr, 1, 0, desc),
                 std::invalid_argument);
    EXPECT_THROW(DistributedMatrix<Type> mat(scalapack_dist, m - 1, n, ptr, 3, 1, desc),
                 std::invalid_argument);
    EXPECT_THROW(DistributedMatrix<Type> mat(scalapack_dist, m, n - 1, ptr, 1, 3, desc),
                 std::invalid_argument);
  }
}

TEST(DistributedMatrixTest, ConvertConstScalapack) {
  using Type = double;
  int m = 11;
  int n = 7;
  int mb = 3;
  int nb = 2;

  auto el_val0 = [](int i, int j) { return 0; };
  auto el_val = [](int i, int j) { return j + .001 * i; };

  for (auto comm_ptr : comms) {
    DistributedMatrix<Type> mat_ref(m, n, mb, nb, *comm_ptr, scalapack_dist);
    int ld = mat_ref.leadingDimension();

    for (int si = 0; si < m; ++si) {
      for (int sj = 0; sj < n; ++sj) {
        for (int sm = 0; sm < m - si; ++sm) {
          for (int sn = 0; sn < n - sj; ++sn) {
            fillDistributedMatrix(mat_ref, el_val0);
            auto psubmat = mat_ref.subMatrix(sm, sn, si, sj);
            DistributedMatrix<Type>& submat_ref = *psubmat;
            auto data_scalapack = submat_ref.getScalapackDescription();
            const Type* const_ptr = std::get<0>(data_scalapack);
            for (auto dist : {scalapack_dist, tile_dist}) {
              fillDistributedMatrix(submat_ref, el_val);

              {
                std::shared_ptr<const DistributedMatrix<Type>> mat =
                    DistributedMatrix<Type>::convertConst(
                        dist, sm, sn, const_ptr, std::get<1>(data_scalapack),
                        std::get<2>(data_scalapack), &std::get<3>(data_scalapack)[0]);
                auto size = std::make_pair(sm, sn);
                Global2DIndex base_index(si, sj);
                if (sm == 0 || sn == 0) {
                  size = std::make_pair(0, 0);
                  base_index = Global2DIndex(0, 0);
                }
                EXPECT_EQ(size, mat->size());
                EXPECT_EQ(std::make_pair(mb, nb), mat->blockSize());
                EXPECT_EQ(base_index, mat->baseIndex());
                if (dist == scalapack_dist)
                  EXPECT_EQ(ld, mat->leadingDimension());
                else
                  EXPECT_EQ(mb, mat->leadingDimension());
                EXPECT_TRUE(checkDistributedMatrix(*mat, el_val));
                EXPECT_EQ(dist, mat->distribution());
              }
            }
          }
        }
      }
    }

    auto data_scalapack = mat_ref.getScalapackDescription();
    const Type* ptr = std::get<0>(data_scalapack);
    auto desc = &std::get<3>(data_scalapack)[0];
    EXPECT_THROW(DistributedMatrix<Type>::convertConst(scalapack_dist, -1, n, ptr, 1, 1, desc),
                 std::invalid_argument);
    EXPECT_THROW(DistributedMatrix<Type>::convertConst(scalapack_dist, m, -1, ptr, 1, 1, desc),
                 std::invalid_argument);
    EXPECT_THROW(DistributedMatrix<Type>::convertConst(scalapack_dist, m, n, ptr, 0, 1, desc),
                 std::invalid_argument);
    EXPECT_THROW(DistributedMatrix<Type>::convertConst(scalapack_dist, m, n, ptr, 1, 0, desc),
                 std::invalid_argument);
    EXPECT_THROW(DistributedMatrix<Type>::convertConst(scalapack_dist, m - 1, n, ptr, 3, 1, desc),
                 std::invalid_argument);
    EXPECT_THROW(DistributedMatrix<Type>::convertConst(scalapack_dist, m, n - 1, ptr, 1, 3, desc),
                 std::invalid_argument);
  }
}
#endif

#ifdef DLA_HAVE_DPLASMA
TEST(DistributedMatrixTest, DPlasmaDescriptor) {
  using Type = double;
  int m = 17;
  int n = 13;
  int mb = 3;
  int nb = 2;

  for (auto comm_ptr : comms) {
    DistributedMatrix<Type> mat(m, n, mb, nb, *comm_ptr, tile_dist);
    auto psubmat1 = mat.subMatrix(m, n, 0, 0);
    auto psubmat2 = mat.subMatrix(m - 1, n - 1, 1, 1);
    auto psubmat3 = mat.subMatrix(3 * mb, 4 * nb, 2 * mb - 1, 2 * nb - 1);
    for (auto pmat : {psubmat1, psubmat2, psubmat3}) {
      auto data_dplasma = pmat->getDPlasmaDescription();
      // The pointer of DPlasma matrix points to (0, 0) of the original matrix.
      EXPECT_EQ(mat.ptr(), std::get<0>(data_dplasma).mat);
      EXPECT_EQ(pmat->leadingNumberOfBlocks(), std::get<0>(data_dplasma).nb_elem_r);
      int nb_elem_c = dla_interface::util::ceilDiv(
          pmat->localBaseIndex().col + pmat->localSize().second, mat.blockSize().second);

      EXPECT_EQ(nb_elem_c, std::get<0>(data_dplasma).nb_elem_c);

      matrix_type type = TypeInfo<Type>::dplasma_type;
      EXPECT_EQ(type, std::get<0>(data_dplasma).super.mtype);
      EXPECT_EQ(matrix_Tile, std::get<0>(data_dplasma).super.storage);
      EXPECT_TRUE(parsec_tiled_matrix_dc_type & std::get<0>(data_dplasma).super.dtype);
      EXPECT_TRUE(two_dim_block_cyclic_type & std::get<0>(data_dplasma).super.dtype);
      EXPECT_EQ(pmat->leadingDimension(), std::get<0>(data_dplasma).super.tileld);
      EXPECT_EQ(mb, std::get<0>(data_dplasma).super.mb);
      EXPECT_EQ(nb, std::get<0>(data_dplasma).super.nb);
      EXPECT_EQ(pmat->leadingDimension() * nb, std::get<0>(data_dplasma).super.bsiz);
      EXPECT_GE(util::ceilDiv(m, mb) * mb, std::get<0>(data_dplasma).super.lm);
      EXPECT_GE(util::ceilDiv(n, nb) * nb, std::get<0>(data_dplasma).super.ln);
      EXPECT_EQ(pmat->baseIndex().row, std::get<0>(data_dplasma).super.i);
      EXPECT_EQ(pmat->baseIndex().col, std::get<0>(data_dplasma).super.j);
      EXPECT_EQ(pmat->size().first, std::get<0>(data_dplasma).super.m);
      EXPECT_EQ(pmat->size().second, std::get<0>(data_dplasma).super.n);
      EXPECT_EQ(comm_ptr->rowOrderedMPICommunicator(), std::get<1>(data_dplasma));
    }

    DistributedMatrix<Type> scalapack_mat(m, n, mb, nb, *comm_ptr, scalapack_dist);
    EXPECT_THROW(scalapack_mat.getDPlasmaDescription(), std::invalid_argument);
  }
}
#endif

TEST(DistributedMatrixTest, TestUtilities) {
  using Type = double;
  auto el_val = [](int i, int j) { return j + .001 * i + 1e-4; };

  DistributedMatrix<Type> mat1(19, 21, 3, 2, *comms[0], scalapack_dist);
  DistributedMatrix<Type> mat2(19, 21, 3, 2, *comms[0], tile_dist);
  fillDistributedMatrix(mat1, el_val);
  fillDistributedMatrix(mat2, el_val);

  auto size = mat1.size();
  EXPECT_EQ(size, mat2.size());
  for (int j = 0; j < size.second; ++j)
    for (int i = 0; i < size.first; ++i)
      if (comms[0]->id2D() == mat1.getRankId2D(i, j)) {
        EXPECT_EQ(mat1(i, j), el_val(i, j));
        EXPECT_EQ(mat2(i, j), el_val(i, j));
      }
      else {
        EXPECT_THROW(mat1(i, j), std::invalid_argument);
        EXPECT_THROW(mat2(i, j), std::invalid_argument);
      }

  NullStream nullout;
  for (int j = 0; j < size.second; ++j)
    for (int i = 0; i < size.first; ++i) {
      if (comms[0]->id2D() == mat1.getRankId2D(i, j)) {
        fillDistributedMatrix(mat1, el_val);
        fillDistributedMatrix(mat2, el_val);
        EXPECT_TRUE(checkDistributedMatrix(mat1, el_val));
        EXPECT_TRUE(checkDistributedMatrix(mat2, el_val));
        EXPECT_TRUE(checkDistributedMatrixSame(mat1, mat2));
        mat1(i, j) *= (1 - 1e-4);
        EXPECT_TRUE(checkNearDistributedMatrix(mat1, el_val, 1.01e-4, el_val(i, j) / 10));
        EXPECT_TRUE(checkNearDistributedMatrix(mat1, el_val, 1.01e-5, 10 * el_val(i, j)));
        EXPECT_FALSE(checkNearDistributedMatrix(mat1, el_val, .99e-4, el_val(i, j) / 10, nullout));
        EXPECT_FALSE(checkNearDistributedMatrix(mat1, el_val, .99e-5, 10 * el_val(i, j), nullout));
        EXPECT_FALSE(checkDistributedMatrix(mat1, el_val, nullout));
        EXPECT_FALSE(checkDistributedMatrixSame(mat1, mat2, nullout));
        mat2(i, j) *= (1 - 1e-4);
        EXPECT_TRUE(checkNearDistributedMatrix(mat2, el_val, 1.01e-4, el_val(i, j) / 10));
        EXPECT_TRUE(checkNearDistributedMatrix(mat2, el_val, 1.01e-5, 10 * el_val(i, j)));
        EXPECT_FALSE(checkNearDistributedMatrix(mat2, el_val, .99e-4, el_val(i, j) / 10, nullout));
        EXPECT_FALSE(checkNearDistributedMatrix(mat2, el_val, .99e-5, 10 * el_val(i, j), nullout));
        EXPECT_FALSE(checkDistributedMatrix(mat2, el_val, nullout));
      }
    }
  DistributedMatrix<Type> mat3(20, 21, 3, 2, *comms[0], tile_dist);
  EXPECT_FALSE(checkDistributedMatrixSame(mat1, mat3, nullout));
  DistributedMatrix<Type> mat4(19, 21, 2, 2, *comms[0], tile_dist);
  EXPECT_FALSE(checkDistributedMatrixSame(mat1, mat4, nullout));
}

int main(int argc, char** argv) {
  comm::CommunicatorManager::initialize(true);

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
  }

  ::testing::InitGoogleTest(&argc, argv);

  ::testing::setMPIListener("results_test_distributed_matrix");

  auto ret = RUN_ALL_TESTS();

  comm::CommunicatorManager::finalize();

  return ret;
}
