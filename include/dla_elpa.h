//
// Distributed Linear Algebra Interface (DLAI)
//
// Copyright (c) 2018-2021, ETH Zurich
// All rights reserved.
//
// Please, refer to the LICENSE file in the root directory.
// SPDX-License-Identifier: BSD-3-Clause
//

#ifndef DLA_INTERFACE_DLA_ELPA_H
#define DLA_INTERFACE_DLA_ELPA_H

#include <unistd.h>

#ifdef DLAI_WITH_ELPA

#include <complex>
#include <mpi.h>
#include "communicator_manager.h"
extern "C" {
// include only some ELPA headers.
#include <limits.h>

#include <elpa/elpa_version.h>

struct elpa_struct;
typedef struct elpa_struct* elpa_t;

struct elpa_autotune_struct;
typedef struct elpa_autotune_struct* elpa_autotune_t;


#include <elpa/elpa_constants.h>
#include <elpa/elpa_generated_c_api.h>
#define complex _Complex
#include <elpa/elpa_generated.h>
#undef complex

const char* elpa_strerr(int elpa_error);
}

#if ELPA_API_VERSION != 20200417
#warning ELPA VERSION ELPA_API_VERSION MIGHT NOT BE SUPPORTED
#endif

namespace dla_interface {
  namespace elpa {
    namespace generic {
      /// C++ Generic equivalent of elpa_set.
      template <class T>
      void set(elpa_t, const char*, T, int*);
      template <>
      inline void set(elpa_t handle, const char* name, int value, int* error) {
        elpa_set_integer(handle, name, value, error);
      }
      template <>
      inline void set(elpa_t handle, const char* name, double value, int* error) {
        elpa_set_double(handle, name, value, error);
      }

      template <class T, class BT>
      void eigenvectors(elpa_t handle, T* a, BT* ev, T* q, int* error);

      template <>
      inline void eigenvectors(elpa_t handle, float* a, float* ev, float* q, int* error) {
        elpa_eigenvectors_f(handle, a, ev, q, error);
      }
      template <>
      inline void eigenvectors(elpa_t handle, double* a, double* ev, double* q, int* error) {
        elpa_eigenvectors_d(handle, a, ev, q, error);
      }
      template <>
      inline void eigenvectors(elpa_t handle, std::complex<float>* a, float* ev,
                               std::complex<float>* q, int* error) {
        using CT = float _Complex;
        elpa_eigenvectors_fc(handle, (CT*)a, ev, (CT*)q, error);
      }
      template <>
      inline void eigenvectors(elpa_t handle, std::complex<double>* a, double* ev,
                               std::complex<double>* q, int* error) {
        using CT = double _Complex;
        elpa_eigenvectors_dc(handle, (CT*)a, ev, (CT*)q, error);
      }
    }

    inline void init() {
      if (elpa_init(ELPA_API_VERSION) != ELPA_OK)
        throw std::invalid_argument("Error: ELPA API version not supported");
    }

    inline void uninit() {
      int error;
      elpa_uninit(&error);
    }

    class Handle {
      public:
      template <class DistMatrix>
      Handle(DistMatrix& a, int nev) {
        auto& comm = a.commGrid();
        int error;

        handle_ = elpa_allocate(&error);
        generic::set(handle_, "na", a.size().first, &error);
        generic::set(handle_, "nev", nev, &error);
        generic::set(handle_, "local_nrows", a.leadingDimension(), &error);
        generic::set(handle_, "local_ncols", a.localSize().second, &error);
        generic::set(handle_, "nblk", a.blockSize().first, &error);
        generic::set(handle_, "mpi_comm_parent", MPI_Comm_c2f(comm.rowOrderedMPICommunicator()),
                     &error);
        generic::set(handle_, "process_row", comm.id2D().first, &error);
        generic::set(handle_, "process_col", comm.id2D().second, &error);

        generic::set(handle_, "solver", (int)ELPA_SOLVER_2STAGE, &error);

        elpa_setup(handle_);
      }

      Handle(const Handle&) = delete;
      Handle(Handle&&) = delete;
      Handle& operator=(const Handle&) = delete;
      Handle& operator=(Handle&&) = delete;

      ~Handle() {
        int error;
        elpa_deallocate(handle_, &error);
      }

      auto& getElpaHandle() {
        return handle_;
      }

      private:
      elpa_t handle_;
    };

    template <class T>
    auto eigenvectors(Handle& handle, T* a, BaseType<T>* ev, T* q, int* error) {
      generic::eigenvectors(handle.getElpaHandle(), a, ev, q, error);
    }
  }
}

#endif

#endif  // DLA_INTERFACE_DLA_ELPA_H
