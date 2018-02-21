#ifndef DLA_INTERFACE_SCALAPACK_H
#define DLA_INTERFACE_SCALAPACK_H

#include <mpi.h>
#include "blacs.h"

#ifdef DLA_HAVE_SCALAPACK
namespace scalapack {

  extern "C" {
  // Descriptor setup
  void descinit_(int* desc, const int* m, const int* n, const int* mb, const int* nb,
                 const int* irsrc, const int* icsrc, const int* ictxt, const int* lld, int* info);
  void descset_(int* desc, const int* m, const int* n, const int* mb, const int* nb,
                const int* irsrc, const int* icsrc, const int* ictxt, const int* lld);

  void psgemm_(const char* trans_a, const char* transb, const int* m, const int* n, const int* k,
               const float* alpha, const float* a, const int* ia, const int* ja, const int* desca,
               const float* b, const int* ib, const int* jb, const int* descb, const float* beta,
               float* c, const int* ic, const int* jc, const int* descc);
  void pdgemm_(const char* trans_a, const char* transb, const int* m, const int* n, const int* k,
               const double* alpha, const double* a, const int* ia, const int* ja, const int* desca,
               const double* b, const int* ib, const int* jb, const int* descb, const double* beta,
               double* c, const int* ic, const int* jc, const int* descc);
  void pcgemm_(const char* trans_a, const char* transb, const int* m, const int* n, const int* k,
               const std::complex<float>* alpha, const std::complex<float>* a, const int* ia,
               const int* ja, const int* desca, const std::complex<float>* b, const int* ib,
               const int* jb, const int* descb, const std::complex<float>* beta,
               std::complex<float>* c, const int* ic, const int* jc, const int* descc);
  void pzgemm_(const char* trans_a, const char* transb, const int* m, const int* n, const int* k,
               const std::complex<double>* alpha, const std::complex<double>* a, const int* ia,
               const int* ja, const int* desca, const std::complex<double>* b, const int* ib,
               const int* jb, const int* descb, const std::complex<double>* beta,
               std::complex<double>* c, const int* ic, const int* jc, const int* descc);

  void psgetrf_(const int* m, const int* n, float* a, const int* ia, const int* ja,
                const int* desca, int* ipiv, int* info);
  void pdgetrf_(const int* m, const int* n, double* a, const int* ia, const int* ja,
                const int* desca, int* ipiv, int* info);
  void pcgetrf_(const int* m, const int* n, std::complex<float>* a, const int* ia, const int* ja,
                const int* desca, int* ipiv, int* info);
  void pzgetrf_(const int* m, const int* n, std::complex<double>* a, const int* ia, const int* ja,
                const int* desca, int* ipiv, int* info);

  void pspotrf_(const char* uplo, const int* n, float* a, const int* ia, const int* ja,
                const int* desca, int* info);
  void pdpotrf_(const char* uplo, const int* n, double* a, const int* ia, const int* ja,
                const int* desca, int* info);
  void pcpotrf_(const char* uplo, const int* n, std::complex<float>* a, const int* ia,
                const int* ja, const int* desca, int* info);
  void pzpotrf_(const char* uplo, const int* n, std::complex<double>* a, const int* ia,
                const int* ja, const int* desca, int* info);
  }
}
#endif

#endif  // DLA_INTERFACE_SCALAPACK_H
