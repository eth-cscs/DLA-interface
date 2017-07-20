#include "util_math.h"

#include <limits>
#include "gtest/gtest.h"

using namespace dla_interface;
using namespace testing;

TEST(MathUtilTest, CeilDiv) {
  EXPECT_EQ(0, util::ceilDiv(0, 1));
  EXPECT_EQ(0, util::ceilDiv(0, 4));
  EXPECT_EQ(3, util::ceilDiv(3, 1));
  EXPECT_EQ(6, util::ceilDiv(42, 7));
  EXPECT_EQ(7, util::ceilDiv(43, 7));
  EXPECT_EQ(7, util::ceilDiv(49, 7));
  EXPECT_EQ(6l, util::ceilDiv(42l, 7l));
  EXPECT_EQ(7l, util::ceilDiv(43l, 7l));
  EXPECT_EQ(7l, util::ceilDiv(49l, 7l));
  EXPECT_EQ(6ll, util::ceilDiv(42ll, 7ll));
  EXPECT_EQ(7ll, util::ceilDiv(43ll, 7ll));
  EXPECT_EQ(7ll, util::ceilDiv(49ll, 7ll));
  EXPECT_EQ(6u, util::ceilDiv(42u, 7u));
  EXPECT_EQ(7u, util::ceilDiv(43u, 7u));
  EXPECT_EQ(7u, util::ceilDiv(49u, 7u));
  EXPECT_EQ(6ul, util::ceilDiv(42ul, 7ul));
  EXPECT_EQ(7ul, util::ceilDiv(43ul, 7ul));
  EXPECT_EQ(7ul, util::ceilDiv(49ul, 7ul));
  EXPECT_EQ(6ull, util::ceilDiv(42ull, 7ull));
  EXPECT_EQ(7ull, util::ceilDiv(43ull, 7ull));
  EXPECT_EQ(7ull, util::ceilDiv(49ull, 7ull));

  EXPECT_TRUE((std::is_same<int, decltype(util::ceilDiv(3, 1))>()));
  EXPECT_TRUE((std::is_same<unsigned int, decltype(util::ceilDiv(3u, 1u))>()));
  short a;
  EXPECT_TRUE((std::is_same<short, decltype(util::ceilDiv(a, a))>()));
  unsigned short ua;
  EXPECT_TRUE((std::is_same<unsigned short, decltype(util::ceilDiv(ua, ua))>()));
  EXPECT_TRUE((std::is_same<long, decltype(util::ceilDiv(3l, 1l))>()));
  EXPECT_TRUE((std::is_same<unsigned long, decltype(util::ceilDiv(3ul, 1ul))>()));
  EXPECT_TRUE((std::is_same<long long, decltype(util::ceilDiv(3ll, 1ll))>()));
  EXPECT_TRUE((std::is_same<unsigned long long, decltype(util::ceilDiv(3ull, 1ull))>()));
}

TEST(MathUtilTest, MultSize) {
  std::size_t res = 56;
  EXPECT_EQ(res, util::multSize(7, 8));
  EXPECT_EQ(res, util::multSize(7u, 8u));
  EXPECT_EQ(res, util::multSize(7l, 8l));
  EXPECT_EQ(res, util::multSize(7ul, 8ul));
  EXPECT_EQ(res, util::multSize(7ll, 8ll));
  EXPECT_EQ(res, util::multSize(7ull, 8ull));

  if (std::numeric_limits<int>::max() < std::numeric_limits<std::size_t>::max()) {
    int b = 2;
    int a = std::numeric_limits<int>::max() / b + 1;
    std::size_t res = static_cast<std::size_t>(a) * static_cast<std::size_t>(b);

    EXPECT_EQ(res, util::multSize(a, b));
    EXPECT_NE(static_cast<std::size_t>(a * b), util::multSize(a, b));
  }

  EXPECT_TRUE((std::is_same<std::size_t, decltype(util::multSize(3, 1))>()));
  EXPECT_TRUE((std::is_same<std::size_t, decltype(util::multSize(3u, 1u))>()));
  short a;
  EXPECT_TRUE((std::is_same<std::size_t, decltype(util::multSize(a, a))>()));
  unsigned short ua;
  EXPECT_TRUE((std::is_same<std::size_t, decltype(util::multSize(ua, ua))>()));
  EXPECT_TRUE((std::is_same<std::size_t, decltype(util::multSize(3l, 1l))>()));
  EXPECT_TRUE((std::is_same<std::size_t, decltype(util::multSize(3ul, 1ul))>()));
  EXPECT_TRUE((std::is_same<std::size_t, decltype(util::multSize(3ll, 1ll))>()));
  EXPECT_TRUE((std::is_same<std::size_t, decltype(util::multSize(3ull, 1ull))>()));
}
