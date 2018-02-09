#include "util_thread.h"

#include <tuple>
#include "gtest/gtest.h"
#include "communicator_manager.h"
#include "thread_binding.h"
#include "thread_blas.h"

using namespace dla_interface;
using namespace testing;

TEST(ThreadUtilTest, SetNumThreadAndCpuBind) {
  auto nr_threads1 = thread::getOmpBlasThreads();
  thread::setOmpBlasThreads(2);
  auto nr_threads2 = thread::getOmpBlasThreads();

  auto cpuset1 = comm::CommunicatorManager::getCpuBind();
  thread::CpuSet cpuset2;
  cpuset2.add(3);
  comm::CommunicatorManager::setCpuBind(cpuset2);
  cpuset2 = comm::CommunicatorManager::getCpuBind();

  {
    util::SetNumThreadsAndCpuBind tmp(nr_threads1, cpuset1);

    EXPECT_EQ(nr_threads1, thread::getOmpBlasThreads());
    EXPECT_EQ(cpuset1, comm::CommunicatorManager::getCpuBind());
  }
  EXPECT_EQ(nr_threads2, thread::getOmpBlasThreads());
  EXPECT_EQ(cpuset2, comm::CommunicatorManager::getCpuBind());

  {
    util::SetNumThreadsAndCpuBind tmp(std::make_tuple(std::cref(nr_threads1), std::cref(cpuset1)));

    EXPECT_EQ(nr_threads1, thread::getOmpBlasThreads());
    EXPECT_EQ(cpuset1, comm::CommunicatorManager::getCpuBind());
  }
  EXPECT_EQ(nr_threads2, thread::getOmpBlasThreads());
  EXPECT_EQ(cpuset2, comm::CommunicatorManager::getCpuBind());
}

int main(int argc, char** argv) {
  comm::CommunicatorManager::initialize(true);
  ::testing::InitGoogleTest(&argc, argv);
  auto ret = RUN_ALL_TESTS();
  comm::CommunicatorManager::finalize();
  return ret;
}
