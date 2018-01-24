#include "thread_binding.h"

#include "gtest/gtest.h"

using namespace dla_interface;
using namespace testing;

TEST(ThreadBindingTest, ThreadBindingRunTest) {
  thread::SystemTopology topo;

  thread::CpuSet cpuset = topo.getCpuBind();

  thread::CpuSet cpuset2;
  cpuset2.add(0);

  topo.setCpuBind(cpuset2);
  ASSERT_EQ(cpuset2, topo.getCpuBind());
  topo.setCpuBind(cpuset);
  ASSERT_EQ(cpuset, topo.getCpuBind());
}

#ifdef DLA_HAVE_HWLOC
std::string expectedSet(const char* str) {
  return str;
}
#else
std::string expectedSet(const char* str) {
  return "<Warning: Cpu binding not Enabled>";
}
#endif

TEST(ThreadBindingTest, CpuSetTest) {
  thread::CpuSet set;
  ASSERT_EQ(expectedSet("0x0"), set.str());

  set.add(0);
  ASSERT_EQ(expectedSet("0x00000001"), set.str());

  set.add(3);
  ASSERT_EQ(expectedSet("0x00000009"), set.str());

  set.add(16);
  ASSERT_EQ(expectedSet("0x00010009"), set.str());

  set.remove(0);
  ASSERT_EQ(expectedSet("0x00010008"), set.str());

  thread::CpuSet set2(set);
  set2.remove(16);
  ASSERT_EQ(expectedSet("0x00010008"), set.str());
  ASSERT_EQ(expectedSet("0x00000008"), set2.str());

  set = set2;
  ASSERT_EQ(expectedSet("0x00000008"), set.str());

  thread::CpuSet set3 = std::move(set2);
  ASSERT_EQ(expectedSet("0x00000008"), set3.str());

  set2 = std::move(set3);
  ASSERT_EQ(expectedSet("0x00000008"), set2.str());

  ASSERT_TRUE(set == set2);
  ASSERT_FALSE(set != set2);
}
