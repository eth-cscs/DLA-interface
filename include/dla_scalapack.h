#ifndef DLA_INTERFACE_DLA_SCALAPACK_H
#define DLA_INTERFACE_DLA_SCALAPACK_H

#include "scalapack.h"
#include "error_message.h"

#ifdef DLA_HAVE_SCALAPACK
namespace dla_interface {
  namespace scalapack_wrappers {

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
