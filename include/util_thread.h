#ifndef DLA_INTERFACE_UTIL_THREAD_H
#define DLA_INTERFACE_UTIL_THREAD_H

#include <tuple>
#include "communicator_manager.h"
#include "thread_binding.h"
#include "thread_blas.h"

namespace dla_interface {
  namespace util {
    class SetNumThreadsAndCpuBind {
      // On costruction sets the BLAS number of threads and cpu binding to the given value.
      // On destruction resets the BLAS number of threads and cpu binding to the original value.

      public:
      SetNumThreadsAndCpuBind(std::tuple<const thread::NumThreads&, const thread::CpuSet&> settings)
       : SetNumThreadsAndCpuBind(std::get<0>(settings), std::get<1>(settings)) {}

      SetNumThreadsAndCpuBind(const thread::NumThreads& nr_threads, const thread::CpuSet& cpuset) {
        nr_threads_ = std::make_unique<thread::NumThreads>(thread::getOmpBlasThreads());
        if (nr_threads != *nr_threads_) {
          thread::setOmpBlasThreads(nr_threads);
        }
        else {
          nr_threads_ = nullptr;
        }

        cpuset_ = std::make_unique<thread::CpuSet>(comm::CommunicatorManager::getCpuBind());
        if (cpuset != *cpuset_) {
          comm::CommunicatorManager::setCpuBind(cpuset);
        }
        else {
          cpuset_ = nullptr;
        }
      }

      ~SetNumThreadsAndCpuBind() {
        if (nr_threads_) {
          thread::setOmpBlasThreads(*nr_threads_);
        }
        if (cpuset_) {
          comm::CommunicatorManager::setCpuBind(*cpuset_);
        }
      }

      SetNumThreadsAndCpuBind(const SetNumThreadsAndCpuBind&) = delete;
      SetNumThreadsAndCpuBind(SetNumThreadsAndCpuBind&&) = delete;
      SetNumThreadsAndCpuBind operator=(const SetNumThreadsAndCpuBind&) = delete;
      SetNumThreadsAndCpuBind operator=(SetNumThreadsAndCpuBind&&) = delete;

      private:
      std::unique_ptr<thread::NumThreads> nr_threads_;
      std::unique_ptr<thread::CpuSet> cpuset_;
    };
  }
}

#endif  // DLA_INTERFACE_UTIL_THREAD_H
