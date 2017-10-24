#ifndef DLA_INTERFACE_DLA_SCALAPACK_H
#define DLA_INTERFACE_DLA_SCALAPACK_H

#include "scalapack.h"
#include "error_message.h"

#ifdef DLA_HAVE_SCALAPACK
namespace dla_interface {
  namespace scalapack_wrappers {

    inline void ppotrf(UpLo uplo, int n, float* a, int ia, int ja, int* desca) {
      const char char_uplo = static_cast<char>(uplo);
      int info = 0;
      scalapack::pspotrf_(&char_uplo, &n, a, &ia, &ja, desca, &info);
      if (info != 0) {
        if (info < 0)
          throw std::invalid_argument(errorMessage("Argument ", -info, " is wrong."));
        else
          throw std::invalid_argument(errorMessage("Matrix is not positive definite (", info, ")"));
      }
    }
    inline void ppotrf(UpLo uplo, int n, double* a, int ia, int ja, int* desca) {
      const char char_uplo = static_cast<char>(uplo);
      int info = 0;
      scalapack::pdpotrf_(&char_uplo, &n, a, &ia, &ja, desca, &info);
      return;
      if (info != 0) {
        if (info < 0)
          throw std::invalid_argument(errorMessage("Argument ", -info, " is wrong."));
        else
          throw std::invalid_argument(errorMessage("Matrix is not positive definite (", info, ")"));
      }
    }
    inline void ppotrf(UpLo uplo, int n, std::complex<float>* a, int ia, int ja, int* desca) {
      const char char_uplo = static_cast<char>(uplo);
      int info = 0;
      scalapack::pcpotrf_(&char_uplo, &n, a, &ia, &ja, desca, &info);
      if (info != 0) {
        if (info < 0)
          throw std::invalid_argument(errorMessage("Argument ", -info, " is wrong."));
        else
          throw std::invalid_argument(errorMessage("Matrix is not positive definite (", info, ")"));
      }
    }
    inline void ppotrf(UpLo uplo, int n, std::complex<double>* a, int ia, int ja, int* desca) {
      const char char_uplo = static_cast<char>(uplo);
      int info = 0;
      scalapack::pzpotrf_(&char_uplo, &n, a, &ia, &ja, desca, &info);
      if (info != 0) {
        if (info < 0)
          throw std::invalid_argument(errorMessage("Argument ", -info, " is wrong."));
        else
          throw std::invalid_argument(errorMessage("Matrix is not positive definite (", info, ")"));
      }
    }
  }
}
#endif

#endif  // DLA_INTERFACE_DLA_SCALAPACK_H
