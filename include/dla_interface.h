#ifndef DLA_INTERFACE_DLA_INTERFACE_H
#define DLA_INTERFACE_DLA_INTERFACE_H

#include "distributed_matrix.h"

namespace dla_interface {

  template <class ElType>
  void choleskyFactorization(UpLo uplo, DistributedMatrix<ElType>& mat, SolverType solver);
}

#endif  // DLA_INTERFACE_DLA_INTERFACE_H
