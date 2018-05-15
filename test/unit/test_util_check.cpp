#include "util_check.h"

#include <utility>
#include <stdexcept>
#include "gtest/gtest.h"
#include "matrix_index.h"

using namespace dla_interface;
using namespace testing;

class MockMatrix {
  public:
  MockMatrix(std::pair<SizeType, SizeType> size) : size_(size) {}
  std::pair<SizeType, SizeType> size() const {
    return size_;
  }

  private:
  std::pair<SizeType, SizeType> size_;
};

class MockComm {
  public:
  MockComm() {}
  MockComm(const MockComm&) = delete;
  MockComm(MockComm&&) = delete;
  MockComm& operator=(const MockComm&) = delete;
  MockComm& operator=(MockComm&&) = delete;
};

class MockDistMatrix {
  public:
  MockDistMatrix(std::pair<SizeType, SizeType> size, std::pair<SizeType, SizeType> block_size,
                 Global2DIndex base_index, const MockComm& comm)
   : size_(size), block_size_(block_size), base_index_(base_index), comm_(&comm) {}
  std::pair<SizeType, SizeType> size() const {
    return size_;
  }
  std::pair<SizeType, SizeType> blockSize() const {
    return block_size_;
  }
  Global2DIndex baseIndex() const {
    return base_index_;
  }
  const MockComm& commGrid() const {
    return *comm_;
  }

  private:
  std::pair<SizeType, SizeType> size_;
  std::pair<SizeType, SizeType> block_size_;
  Global2DIndex base_index_;
  const MockComm* comm_;
};

TEST(CheckUtilTest, CheckIsSquare) {
  MockComm comm;
  for (int m : {3, 7, 14, 15}) {
    for (int n : {3, 7, 14, 15}) {
      MockMatrix mat_test(std::make_pair(m, n));
      MockDistMatrix dist_mat_test(std::make_pair(m, n), std::make_pair(m / 2, n / 2),
                                   Global2DIndex(m, n), comm);
      if (m == n) {
        EXPECT_NO_THROW(dlai__util__checkIsSquare(mat_test));
        EXPECT_NO_THROW(dlai__util__checkIsSquare(dist_mat_test));
      }
      else {
        EXPECT_THROW(dlai__util__checkIsSquare(mat_test), std::invalid_argument);
        EXPECT_THROW(dlai__util__checkIsSquare(dist_mat_test), std::invalid_argument);
      }
    }
  }
}

TEST(CheckUtilTest, CheckBlocksAreSquare) {
  MockComm comm;
  for (int mb : {3, 7, 14, 15}) {
    for (int nb : {3, 7, 14, 15}) {
      MockDistMatrix dist_mat_test(std::make_pair(7, 5), std::make_pair(mb, nb),
                                   Global2DIndex(0, 2), comm);
      if (mb == nb)
        EXPECT_NO_THROW(dlai__util__checkBlocksAreSquare(dist_mat_test));
      else
        EXPECT_THROW(dlai__util__checkBlocksAreSquare(dist_mat_test), std::invalid_argument);
    }
  }
}

TEST(CheckUtilTest, CheckBaseIndexAtBlock) {
  MockComm comm;
  for (int mb : {3, 7, 14, 15}) {
    for (int nb : {3, 7, 14, 15}) {
      for (int ia : {3, 7, 14, 15}) {
        for (int ja : {3, 7, 14, 15}) {
          MockDistMatrix dist_mat_test(std::make_pair(7, 5), std::make_pair(mb, nb),
                                       Global2DIndex(ia, ja), comm);
          if (ia % mb == 0 && ja % nb == 0)
            EXPECT_NO_THROW(dlai__util__checkBaseIndexAtBlock(dist_mat_test));
          else
            EXPECT_THROW(dlai__util__checkBaseIndexAtBlock(dist_mat_test), std::invalid_argument);
        }
      }
    }
  }
}

TEST(CheckUtilTest, CheckSameComm2D) {
  MockComm comm1;
  MockComm comm2;
  MockComm comm3;

  std::vector<MockDistMatrix> mat1;
  std::vector<MockDistMatrix> mat2;
  std::vector<MockDistMatrix> mat3;
  for (int m : {13, 17}) {
    for (int n : {18, 23}) {
      for (int mb : {3, 5}) {
        for (int nb : {3, 5}) {
          for (int ia : {3, 7}) {
            for (int ja : {3, 7}) {
              mat1.emplace_back(std::make_pair(m, n), std::make_pair(mb, nb), Global2DIndex(ia, ja),
                                comm1);
              mat2.emplace_back(std::make_pair(m, n), std::make_pair(mb, nb), Global2DIndex(ia, ja),
                                comm2);
              mat3.emplace_back(std::make_pair(m, n), std::make_pair(mb, nb), Global2DIndex(ia, ja),
                                comm3);
            }
          }
        }
      }
    }
  }
  int i0 = 0;
  for (size_t i1 = 0; i1 < mat1.size(); ++i1) {
    EXPECT_NO_THROW(dlai__util__checkSameComm2D(mat1[i0], mat1[i1]));
    EXPECT_NO_THROW(dlai__util__checkSameComm2D(mat2[i0], mat2[i1]));
    EXPECT_NO_THROW(dlai__util__checkSameComm2D(mat3[i0], mat3[i1]));
    EXPECT_THROW(dlai__util__checkSameComm2D(mat1[i0], mat2[i1]), std::invalid_argument);
    EXPECT_THROW(dlai__util__checkSameComm2D(mat1[i0], mat3[i1]), std::invalid_argument);
    for (size_t i2 = 0; i2 < mat1.size(); ++i2) {
      EXPECT_NO_THROW(dlai__util__checkSameComm2D(mat1[i0], mat1[i1], mat1[i2]));
      EXPECT_NO_THROW(dlai__util__checkSameComm2D(mat2[i0], mat2[i1], mat2[i2]));
      EXPECT_NO_THROW(dlai__util__checkSameComm2D(mat3[i0], mat3[i1], mat3[i2]));
      EXPECT_THROW(dlai__util__checkSameComm2D(mat1[i0], mat1[i1], mat2[i2]), std::invalid_argument);
      EXPECT_THROW(dlai__util__checkSameComm2D(mat1[i0], mat2[i1], mat1[i2]), std::invalid_argument);
      EXPECT_THROW(dlai__util__checkSameComm2D(mat1[i0], mat2[i1], mat2[i2]), std::invalid_argument);
      EXPECT_THROW(dlai__util__checkSameComm2D(mat1[i0], mat2[i1], mat3[i2]), std::invalid_argument);
    }
  }
}
