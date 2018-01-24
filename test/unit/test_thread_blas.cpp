#include "thread_blas.h"

#include <omp.h>
#include "gtest/gtest.h"
#ifdef DLA_HAVE_MKL_NUM_THREADS_UTIL
#include <mkl.h>
#endif

using namespace dla_interface;
using namespace testing;

TEST(ThreadBlasTest, GetNumThreadsTest) {
  omp_set_num_threads(4);
#if DLA_HAVE_MKL_NUM_THREADS_UTIL
  mkl_set_num_threads(4);
#endif

  thread::NumThreads n_threads = thread::getOmpBlasThreads();

  EXPECT_EQ(omp_get_max_threads(), n_threads.omp_num_threads);
#if DLA_HAVE_MKL_NUM_THREADS_UTIL
  EXPECT_EQ(mkl_get_max_threads(), n_threads.mkl_num_threads);
#endif

  omp_set_num_threads(2);
#if DLA_HAVE_MKL_NUM_THREADS_UTIL
  mkl_set_num_threads(2);
#endif

  n_threads = thread::getOmpBlasThreads();

  EXPECT_EQ(omp_get_max_threads(), n_threads.omp_num_threads);
#if DLA_HAVE_MKL_NUM_THREADS_UTIL
  EXPECT_EQ(mkl_get_max_threads(), n_threads.mkl_num_threads);
#endif
}

TEST(ThreadBlasTest, SetNumThreadsTest1) {
  omp_set_num_threads(2);
  int ref_omp_threads = omp_get_max_threads();

  omp_set_num_threads(4);

#if DLA_HAVE_MKL_NUM_THREADS_UTIL
  mkl_set_num_threads(2);
  int ref_mkl_threads = mkl_get_max_threads();

  mkl_set_num_threads(4);
#endif

  thread::setOmpBlasThreads(2);

  EXPECT_EQ(ref_omp_threads, omp_get_max_threads());
#if DLA_HAVE_MKL_NUM_THREADS_UTIL
  EXPECT_EQ(ref_mkl_threads, mkl_get_max_threads());
#endif
}

TEST(ThreadBlasTest, SetNumThreadsTest2) {
  thread::NumThreads n_threads;

  n_threads.omp_num_threads = 3;
  omp_set_num_threads(n_threads.omp_num_threads);
  int ref_omp_threads = omp_get_max_threads();

  omp_set_num_threads(4);

#if DLA_HAVE_MKL_NUM_THREADS_UTIL
  n_threads.mkl_num_threads = 3;
  mkl_set_num_threads(n_threads.mkl_num_threads);
  int ref_mkl_threads = mkl_get_max_threads();

  mkl_set_num_threads(2);
#endif

  thread::setOmpBlasThreads(n_threads);

  EXPECT_EQ(ref_omp_threads, omp_get_max_threads());
#if DLA_HAVE_MKL_NUM_THREADS_UTIL
  EXPECT_EQ(ref_mkl_threads, mkl_get_max_threads());
#endif
}
