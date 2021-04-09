#ifndef DLA_INTERFACE_SCALAPACK_H
#define DLA_INTERFACE_SCALAPACK_H

#include <mpi.h>
#include "blacs.h"

#ifdef DLAI_WITH_SCALAPACK
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

  void pspotri_(const char* uplo, const int* n, float* a, const int* ia, const int* ja,
                const int* desca, int* info);
  void pdpotri_(const char* uplo, const int* n, double* a, const int* ia, const int* ja,
                const int* desca, int* info);
  void pcpotri_(const char* uplo, const int* n, std::complex<float>* a, const int* ia,
                const int* ja, const int* desca, int* info);
  void pzpotri_(const char* uplo, const int* n, std::complex<double>* a, const int* ia,
                const int* ja, const int* desca, int* info);

  void pssyevd_(const char* jobz, const char* uplo, const int* n, float* a, const int* ia,
                const int* ja, const int* desca, float* w, float* z, const int* iz, const int* jz,
                const int* descz, float* work, const int* lwork, int* iwork, const int* liwork,
                int* info);
  void pdsyevd_(const char* jobz, const char* uplo, const int* n, double* a, const int* ia,
                const int* ja, const int* desca, double* w, double* z, const int* iz, const int* jz,
                const int* descz, double* work, const int* lwork, int* iwork, const int* liwork,
                int* info);
  void pcheevd_(const char* jobz, const char* uplo, const int* n, std::complex<float>* a,
                const int* ia, const int* ja, const int* desca, float* w, std::complex<float>* z,
                const int* iz, const int* jz, const int* descz, std::complex<float>* work,
                const int* lwork, float* rwork, const int* lrwork, int* iwork, const int* liwork,
                int* info);
  void pzheevd_(const char* jobz, const char* uplo, const int* n, std::complex<double>* a,
                const int* ia, const int* ja, const int* desca, double* w, std::complex<double>* z,
                const int* iz, const int* jz, const int* descz, std::complex<double>* work,
                const int* lwork, double* rwork, const int* lrwork, int* iwork, const int* liwork,
                int* info);

  void pstradd_(const char* uplo, const char* trans, const int* m, const int* n, const float* alpha,
                const float* a, const int* ia, const int* ja, const int* desca, const float* beta,
                float* c, const int* ic, const int* jc, const int* descc);
  void pdtradd_(const char* uplo, const char* trans, const int* m, const int* n, const double* alpha,
                const double* a, const int* ia, const int* ja, const int* desca, const double* beta,
                double* c, const int* ic, const int* jc, const int* descc);
  void pctradd_(const char* uplo, const char* trans, const int* m, const int* n,
                const std::complex<float>* alpha, const std::complex<float>* a, const int* ia,
                const int* ja, const int* desca, const std::complex<float>* beta,
                std::complex<float>* c, const int* ic, const int* jc, const int* descc);
  void pztradd_(const char* uplo, const char* trans, const int* m, const int* n,
                const std::complex<double>* alpha, const std::complex<double>* a, const int* ia,
                const int* ja, const int* desca, const std::complex<double>* beta,
                std::complex<double>* c, const int* ic, const int* jc, const int* descc);

  void pstrtri_(const char* uplo, const char* diag, const int* n, float* a, const int* ia,
                const int* ja, const int* desca, int* info);
  void pdtrtri_(const char* uplo, const char* diag, const int* n, double* a, const int* ia,
                const int* ja, const int* desca, int* info);
  void pctrtri_(const char* uplo, const char* diag, const int* n, std::complex<float>* a,
                const int* ia, const int* ja, const int* desca, int* info);
  void pztrtri_(const char* uplo, const char* diag, const int* n, std::complex<double>* a,
                const int* ia, const int* ja, const int* desca, int* info);
  }
}
#endif

#endif  // DLA_INTERFACE_SCALAPACK_H
