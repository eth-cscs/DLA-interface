#include "util_memory.h"

#include <limits>
#include "gtest/gtest.h"

using namespace dla_interface;
using namespace testing;

template <typename Type>
class MemoryUtilTest : public ::testing::Test {};

typedef ::testing::Types<int, float, double, unsigned long long> TestTypes;
TYPED_TEST_CASE(MemoryUtilTest, TestTypes);

TYPED_TEST(MemoryUtilTest, Copy) {
  using Type = TypeParam;
  constexpr std::size_t sz = 10;

  std::vector<Type> src(sz);
  for (std::size_t i = 0; i < sz; ++i)
    src[i] = static_cast<Type>(i);

  std::size_t start = 2;
  for (std::size_t len = 0; len < 4; ++len) {
    std::vector<Type> dst(sz, 0);
    const Type* src_ptr = &src[start];
    util::memory::copy(len, src_ptr, &dst[start]);

    for (size_t i = 0; i < sz; ++i) {
      Type exp = i >= start && i < start + len ? src[i] : 0;
      EXPECT_EQ(exp, dst[i]);
    }
  }
}

TYPED_TEST(MemoryUtilTest, Copy2D) {
  using Type = TypeParam;
  int ld_src = 10;
  constexpr int n_src = 10;

  std::vector<Type> src(ld_src * n_src);
  for (int j = 0; j < n_src; ++j) {
    for (int i = 0; i < ld_src; ++i) {
      src[i + ld_src * j] = static_cast<Type>(100 * i + j);
    }
  }

  int i_start = 1;
  int j_start = 1;
  int ld_dest = 8;
  for (int m = 0; m < 3; ++m) {
    for (int n = 0; n < 3; ++n) {
      std::vector<Type> dst(ld_dest * n_src, 0);
      const Type* src_ptr = &src[i_start + ld_src * j_start];
      util::memory::copy2D(std::make_pair(m, n), src_ptr, ld_src, &dst[i_start + ld_dest * j_start],
                           ld_dest);

      for (int j = 0; j < n_src; ++j) {
        for (int i = 0; i < ld_dest; ++i) {
          bool i_ok = i >= i_start && i < i_start + m;
          bool j_ok = j >= j_start && j < j_start + n;

          Type exp = i_ok && j_ok ? src[i + ld_src * j] : 0;
          EXPECT_EQ(exp, dst[i + ld_dest * j]);
        }
      }
    }
  }

  // one copy case: ld_src = ld_dest = m
  ld_dest = ld_src;
  i_start = 0;
  int m = ld_src;
  int n = 5;
  std::vector<Type> dst(ld_dest * n_src, 0);
  const Type* src_ptr = &src[i_start + ld_src * j_start];
  util::memory::copy2D(std::make_pair(m, n), src_ptr, ld_src, &dst[i_start + ld_dest * j_start],
                       ld_dest);

  for (int j = 0; j < n_src; ++j) {
    for (int i = 0; i < ld_dest; ++i) {
      bool i_ok = i >= i_start && i < i_start + m;
      bool j_ok = j >= j_start && j < j_start + n;

      Type exp = i_ok && j_ok ? src[i + ld_src * j] : 0;
      EXPECT_EQ(exp, dst[i + ld_dest * j]);
    }
  }
}
