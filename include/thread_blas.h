//
// Distributed Linear Algebra Interface (DLAI)
//
// Copyright (c) 2018-2021, ETH Zurich
// All rights reserved.
//
// Please, refer to the LICENSE file in the root directory.
// SPDX-License-Identifier: BSD-3-Clause
//

#ifndef DLA_INTERFACE_THREAD_BLAS_H
#define DLA_INTERFACE_THREAD_BLAS_H

#include <iostream>
#include <omp.h>

#ifdef DLA_HAVE_MKL_NUM_THREADS_UTIL
#include <mkl_service.h>
#endif

namespace dla_interface {
  namespace thread {

    struct NumThreads {
      NumThreads(int nr_threads = 0)
       : omp_num_threads(nr_threads)
#ifdef DLA_HAVE_MKL_NUM_THREADS_UTIL
         ,
         mkl_num_threads(nr_threads)
#endif
      {
      }

      bool operator==(const NumThreads& rhs) const {
        bool res = omp_num_threads == rhs.omp_num_threads;
#ifdef DLA_HAVE_MKL_NUM_THREADS_UTIL
        res = res && mkl_num_threads == rhs.mkl_num_threads;
#endif
        return res;
      }

      bool operator!=(const NumThreads& rhs) const {
        return !operator==(rhs);
      }

      int omp_num_threads;
#ifdef DLA_HAVE_MKL_NUM_THREADS_UTIL
      int mkl_num_threads;
#endif
    };

    /// Returns the number of threads used by BLAS (OpenMP number of threads and MKL number of
    /// threads if DLA_HAVE_MKL_NUM_THREADS_UTIL is defined).
    inline NumThreads getOmpBlasThreads() {
      NumThreads num_threads;
#ifdef DLA_HAVE_MKL_NUM_THREADS_UTIL
      num_threads.mkl_num_threads = mkl_get_max_threads();
#endif
      num_threads.omp_num_threads = omp_get_max_threads();
      return num_threads;
    }

    inline void printOmpBlasThreadDebugInfo(const char* str) {
#ifdef DLA_THREAD_DEBUG_INFO
      auto tmp = getOmpBlasThreads();
      std::cout << str << " " << tmp.omp_num_threads
#ifdef DLA_HAVE_MKL_NUM_THREADS_UTIL
                << ", " << tmp.mkl_num_threads
#endif
                << std::endl;
#endif
    }

    /// Sets the number of threads used by BLAS using
    ///
    /// - omp_set_num_threads and,
    /// - mkl_set_num_threads if DLA_HAVE_MKL_NUM_THREADS_UTIL is defined.
    inline void setOmpBlasThreads(NumThreads nr_threads) {
      printOmpBlasThreadDebugInfo("setOmpBlasThreads called:\n  Previous value:");

      if (nr_threads.omp_num_threads > 0) {
        omp_set_num_threads(nr_threads.omp_num_threads);
      }
#ifdef DLA_HAVE_MKL_NUM_THREADS_UTIL
      if (nr_threads.mkl_num_threads > 0) {
        mkl_set_num_threads(nr_threads.mkl_num_threads);
      }
#endif

      printOmpBlasThreadDebugInfo("  New      value:");
    }

    /// Sets the number of threads used by BLAS if nr_threads > 0 using
    ///
    /// - omp_set_num_threads and,
    /// - mkl_set_num_threads if DLA_HAVE_MKL_NUM_THREADS_UTIL is defined.
    inline void setOmpBlasThreads(int nr_threads) {
      setOmpBlasThreads(NumThreads(nr_threads));
    }
  }
}

#endif  // DLA_INTERFACE_THREAD_BLAS_H
