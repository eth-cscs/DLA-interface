#ifndef DLA_INTERFACE_SCALAPACK_H
#define DLA_INTERFACE_SCALAPACK_H

#include <mpi.h>
#include <blacs.h>

#ifdef DLA_HAVE_SCALAPACK
namespace scalapack {

  extern "C" {
  // Descriptor setup
  void descinit_(int* desc, const int* m, const int* n, const int* mb, const int* nb,
                 const int* irsrc, const int* icsrc, const int* ictxt, const int* lld, int* info);
  void descset_(int* desc, const int* m, const int* n, const int* mb, const int* nb,
                const int* irsrc, const int* icsrc, const int* ictxt, const int* lld);
  }
}
#endif

#endif  // DLA_INTERFACE_SCALAPACK_H
