//
// Distributed Linear Algebra Interface (DLAI)
//
// Copyright (c) 2018-2019, ETH Zurich
// All rights reserved.
//
// Please, refer to the LICENSE file in the root directory.
// SPDX-License-Identifier: BSD-3-Clause
//

#ifndef DLA_INTERFACE_TIMER_H
#define DLA_INTERFACE_TIMER_H

#include <chrono>
#include <iostream>
#include <vector>
#include <mpi.h>

namespace dla_interface {
  namespace util {

    template <class clock = std::chrono::high_resolution_clock>
    class Timer {
      public:
      // Initialize an enabled timer if enabled == true or a disabled timer.
      Timer(bool enabled = true) : enabled_(enabled), do_barrier_(false), time_points_(1, now()) {}

      // Initialize an enabled timer if enabled == true or a disabled timer, which performs
      // MPI_Barrier
      // on comm before measuring time.
      Timer(MPI_Comm comm, bool enabled = true)
       : enabled_(enabled), do_barrier_(true), comm_(comm), time_points_(1, now()) {}

      // If enabled stores the timepoint and return its index,
      // Otherwise returns -1.
      int save_time() {
        if (!enabled_)
          return -1;
        time_points_.push_back(now());
        return time_points_.size() - 1;
      }

      // If enabled returns the amount of time that has passed (in seconds)
      // from time point with index start_index
      // to time_point with index end_index.
      // Otherwise returns 0;
      // Precondition: the indeces are either 0 or indeces returned by save_time.
      double elapsed(int start_index, int end_index) const {
        if (!enabled_)
          return 0;
        return std::chrono::duration_cast<std::chrono::duration<double>>(time_points_[end_index] -
                                                                         time_points_[start_index])
            .count();
      }

      // If enabled outputs to std::cout
      // - "<str><el> s, <nr_ops/el/1e9> GFlops" if nr_ops >= 0,
      // - "<str><el> s" otherwise,
      // where el is elapsed(start_index, end_index)
      void print_elapsed(int start_index, int end_index, std::string str, double nr_ops = -1) const {
        print_elapsed(start_index, end_index, str, std::cout, nr_ops);
      }

      // If enabled output (using operator<<(Out& out, const T& obj))
      // - "<str><el> s, <nr_ops/el/1e9> GFlops" if nr_ops >= 0,
      // - "<str><el> s" otherwise,
      // where el is elapsed(start_index, end_index)
      template <class Out>
      void print_elapsed(int start_index, int end_index, const std::string str, Out& out,
                         double nr_ops = -1) const {
        if (!enabled_)
          return;

        double el = elapsed(start_index, end_index);
        out << str << el << " s";
        if (nr_ops > 0)
          out << ", " << nr_ops / el / 1e9 << " GFlop/s";
        out << std::endl;
      }

      private:
      using time_point = std::chrono::time_point<clock>;

      bool enabled_;
      bool do_barrier_;
      MPI_Comm comm_;
      std::vector<time_point> time_points_;

      time_point now() const {
        if (do_barrier_) {
          MPI_Barrier(comm_);
        }
        return clock::now();
      }
    };
  }
}

#endif  // DLA_INTERFACE_TIMER_H
