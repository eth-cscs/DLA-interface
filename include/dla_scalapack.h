#ifndef DLA_INTERFACE_DLA_SCALAPACK_H
#define DLA_INTERFACE_DLA_SCALAPACK_H

#include "scalapack.h"
#include "error_message.h"

#ifdef DLA_HAVE_SCALAPACK
namespace dla_interface {
  namespace scalapack_wrappers {

    inline void pgemm(OpTrans trans_a, OpTrans trans_b, int m, int n, int k, float alpha, float* a,
                      int ia, int ja, int* desca, float* b, int ib, int jb, int* descb, float beta,
                      float* c, int ic, int jc, int* descc) {
      const char char_trans_a = static_cast<char>(trans_a);
      const char char_trans_b = static_cast<char>(trans_b);
      scalapack::psgemm_(&char_trans_a, &char_trans_b, &m, &n, &k, &alpha, a, &ia, &ja, desca, b,
                         &ib, &jb, descb, &beta, c, &ic, &jc, descc);
    }
    inline void pgemm(OpTrans trans_a, OpTrans trans_b, int m, int n, int k, double alpha,
                      double* a, int ia, int ja, int* desca, double* b, int ib, int jb, int* descb,
                      double beta, double* c, int ic, int jc, int* descc) {
      const char char_trans_a = static_cast<char>(trans_a);
      const char char_trans_b = static_cast<char>(trans_b);
      scalapack::pdgemm_(&char_trans_a, &char_trans_b, &m, &n, &k, &alpha, a, &ia, &ja, desca, b,
                         &ib, &jb, descb, &beta, c, &ic, &jc, descc);
    }
    inline void pgemm(OpTrans trans_a, OpTrans trans_b, int m, int n, int k,
                      std::complex<float> alpha, std::complex<float>* a, int ia, int ja, int* desca,
                      std::complex<float>* b, int ib, int jb, int* descb, std::complex<float> beta,
                      std::complex<float>* c, int ic, int jc, int* descc) {
      const char char_trans_a = static_cast<char>(trans_a);
      const char char_trans_b = static_cast<char>(trans_b);
      scalapack::pcgemm_(&char_trans_a, &char_trans_b, &m, &n, &k, &alpha, a, &ia, &ja, desca, b,
                         &ib, &jb, descb, &beta, c, &ic, &jc, descc);
    }
    inline void pgemm(OpTrans trans_a, OpTrans trans_b, int m, int n, int k,
                      std::complex<double> alpha, std::complex<double>* a, int ia, int ja,
                      int* desca, std::complex<double>* b, int ib, int jb, int* descb,
                      std::complex<double> beta, std::complex<double>* c, int ic, int jc, int* descc) {
      const char char_trans_a = static_cast<char>(trans_a);
      const char char_trans_b = static_cast<char>(trans_b);
      scalapack::pzgemm_(&char_trans_a, &char_trans_b, &m, &n, &k, &alpha, a, &ia, &ja, desca, b,
                         &ib, &jb, descb, &beta, c, &ic, &jc, descc);
    }

    inline void ppotrf(UpLo uplo, int n, float* a, int ia, int ja, int* desca, int& info) {
      const char char_uplo = static_cast<char>(uplo);
      scalapack::pspotrf_(&char_uplo, &n, a, &ia, &ja, desca, &info);
    }
    inline void ppotrf(UpLo uplo, int n, double* a, int ia, int ja, int* desca, int& info) {
      const char char_uplo = static_cast<char>(uplo);
      scalapack::pdpotrf_(&char_uplo, &n, a, &ia, &ja, desca, &info);
    }
    inline void ppotrf(UpLo uplo, int n, std::complex<float>* a, int ia, int ja, int* desca,
                       int& info) {
      const char char_uplo = static_cast<char>(uplo);
      scalapack::pcpotrf_(&char_uplo, &n, a, &ia, &ja, desca, &info);
    }
    inline void ppotrf(UpLo uplo, int n, std::complex<double>* a, int ia, int ja, int* desca,
                       int& info) {
      const char char_uplo = static_cast<char>(uplo);
      scalapack::pzpotrf_(&char_uplo, &n, a, &ia, &ja, desca, &info);
    }
  }
}
#endif

#endif  // DLA_INTERFACE_DLA_SCALAPACK_H
