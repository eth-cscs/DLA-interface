//
// Distributed Linear Algebra Interface (DLAI)
//
// Copyright (c) 2018-2021, ETH Zurich
// All rights reserved.
//
// Please, refer to the LICENSE file in the root directory.
// SPDX-License-Identifier: BSD-3-Clause
//

#ifndef DLA_INTERFACE_DLA_SCALAPACK_H
#define DLA_INTERFACE_DLA_SCALAPACK_H

#include <complex>
#include <vector>
#include "scalapack.h"
#include "error_message.h"
#include "util_thread.h"

#ifdef DLAI_WITH_SCALAPACK
namespace dla_interface {
  namespace scalapack_wrappers {

    inline void pgemm(OpTrans trans_a, OpTrans trans_b, int m, int n, int k, float alpha,
                      const float* a, int ia, int ja, int* desca, const float* b, int ib, int jb,
                      int* descb, float beta, float* c, int ic, int jc, int* descc) {
      const char char_trans_a = static_cast<char>(trans_a);
      const char char_trans_b = static_cast<char>(trans_b);

      util::SetNumThreadsAndCpuBind config(comm::CommunicatorManager::getScalapackConfigInfo());
      scalapack::psgemm_(&char_trans_a, &char_trans_b, &m, &n, &k, &alpha, a, &ia, &ja, desca, b,
                         &ib, &jb, descb, &beta, c, &ic, &jc, descc);
    }
    inline void pgemm(OpTrans trans_a, OpTrans trans_b, int m, int n, int k, double alpha,
                      const double* a, int ia, int ja, int* desca, const double* b, int ib, int jb,
                      int* descb, double beta, double* c, int ic, int jc, int* descc) {
      const char char_trans_a = static_cast<char>(trans_a);
      const char char_trans_b = static_cast<char>(trans_b);

      util::SetNumThreadsAndCpuBind config(comm::CommunicatorManager::getScalapackConfigInfo());
      scalapack::pdgemm_(&char_trans_a, &char_trans_b, &m, &n, &k, &alpha, a, &ia, &ja, desca, b,
                         &ib, &jb, descb, &beta, c, &ic, &jc, descc);
    }
    inline void pgemm(OpTrans trans_a, OpTrans trans_b, int m, int n, int k,
                      std::complex<float> alpha, const std::complex<float>* a, int ia, int ja,
                      int* desca, const std::complex<float>* b, int ib, int jb, int* descb,
                      std::complex<float> beta, std::complex<float>* c, int ic, int jc, int* descc) {
      const char char_trans_a = static_cast<char>(trans_a);
      const char char_trans_b = static_cast<char>(trans_b);

      util::SetNumThreadsAndCpuBind config(comm::CommunicatorManager::getScalapackConfigInfo());
      scalapack::pcgemm_(&char_trans_a, &char_trans_b, &m, &n, &k, &alpha, a, &ia, &ja, desca, b,
                         &ib, &jb, descb, &beta, c, &ic, &jc, descc);
    }
    inline void pgemm(OpTrans trans_a, OpTrans trans_b, int m, int n, int k,
                      std::complex<double> alpha, const std::complex<double>* a, int ia, int ja,
                      int* desca, const std::complex<double>* b, int ib, int jb, int* descb,
                      std::complex<double> beta, std::complex<double>* c, int ic, int jc, int* descc) {
      const char char_trans_a = static_cast<char>(trans_a);
      const char char_trans_b = static_cast<char>(trans_b);

      util::SetNumThreadsAndCpuBind config(comm::CommunicatorManager::getScalapackConfigInfo());
      scalapack::pzgemm_(&char_trans_a, &char_trans_b, &m, &n, &k, &alpha, a, &ia, &ja, desca, b,
                         &ib, &jb, descb, &beta, c, &ic, &jc, descc);
    }

    inline void pgetrf(int m, int n, float* a, int ia, int ja, int* desca, int* ipiv, int& info) {
      util::SetNumThreadsAndCpuBind config(comm::CommunicatorManager::getScalapackConfigInfo());
      scalapack::psgetrf_(&m, &n, a, &ia, &ja, desca, ipiv, &info);
    }
    inline void pgetrf(int m, int n, double* a, int ia, int ja, int* desca, int* ipiv, int& info) {
      util::SetNumThreadsAndCpuBind config(comm::CommunicatorManager::getScalapackConfigInfo());
      scalapack::pdgetrf_(&m, &n, a, &ia, &ja, desca, ipiv, &info);
    }
    inline void pgetrf(int m, int n, std::complex<float>* a, int ia, int ja, int* desca, int* ipiv,
                       int& info) {
      util::SetNumThreadsAndCpuBind config(comm::CommunicatorManager::getScalapackConfigInfo());
      scalapack::pcgetrf_(&m, &n, a, &ia, &ja, desca, ipiv, &info);
    }
    inline void pgetrf(int m, int n, std::complex<double>* a, int ia, int ja, int* desca, int* ipiv,
                       int& info) {
      util::SetNumThreadsAndCpuBind config(comm::CommunicatorManager::getScalapackConfigInfo());
      scalapack::pzgetrf_(&m, &n, a, &ia, &ja, desca, ipiv, &info);
    }

    inline void ppotrf(UpLo uplo, int n, float* a, int ia, int ja, int* desca, int& info) {
      const char char_uplo = static_cast<char>(uplo);

      util::SetNumThreadsAndCpuBind config(comm::CommunicatorManager::getScalapackConfigInfo());
      scalapack::pspotrf_(&char_uplo, &n, a, &ia, &ja, desca, &info);
    }
    inline void ppotrf(UpLo uplo, int n, double* a, int ia, int ja, int* desca, int& info) {
      const char char_uplo = static_cast<char>(uplo);

      util::SetNumThreadsAndCpuBind config(comm::CommunicatorManager::getScalapackConfigInfo());
      scalapack::pdpotrf_(&char_uplo, &n, a, &ia, &ja, desca, &info);
    }
    inline void ppotrf(UpLo uplo, int n, std::complex<float>* a, int ia, int ja, int* desca,
                       int& info) {
      const char char_uplo = static_cast<char>(uplo);

      util::SetNumThreadsAndCpuBind config(comm::CommunicatorManager::getScalapackConfigInfo());
      scalapack::pcpotrf_(&char_uplo, &n, a, &ia, &ja, desca, &info);
    }
    inline void ppotrf(UpLo uplo, int n, std::complex<double>* a, int ia, int ja, int* desca,
                       int& info) {
      const char char_uplo = static_cast<char>(uplo);

      util::SetNumThreadsAndCpuBind config(comm::CommunicatorManager::getScalapackConfigInfo());
      scalapack::pzpotrf_(&char_uplo, &n, a, &ia, &ja, desca, &info);
    }

    inline void ppotri(UpLo uplo, int n, float* a, int ia, int ja, int* desca, int& info) {
      const char char_uplo = static_cast<char>(uplo);

      util::SetNumThreadsAndCpuBind config(comm::CommunicatorManager::getScalapackConfigInfo());
      scalapack::pspotri_(&char_uplo, &n, a, &ia, &ja, desca, &info);
    }
    inline void ppotri(UpLo uplo, int n, double* a, int ia, int ja, int* desca, int& info) {
      const char char_uplo = static_cast<char>(uplo);

      util::SetNumThreadsAndCpuBind config(comm::CommunicatorManager::getScalapackConfigInfo());
      scalapack::pdpotri_(&char_uplo, &n, a, &ia, &ja, desca, &info);
    }
    inline void ppotri(UpLo uplo, int n, std::complex<float>* a, int ia, int ja, int* desca,
                       int& info) {
      const char char_uplo = static_cast<char>(uplo);

      util::SetNumThreadsAndCpuBind config(comm::CommunicatorManager::getScalapackConfigInfo());
      scalapack::pcpotri_(&char_uplo, &n, a, &ia, &ja, desca, &info);
    }
    inline void ppotri(UpLo uplo, int n, std::complex<double>* a, int ia, int ja, int* desca,
                       int& info) {
      const char char_uplo = static_cast<char>(uplo);

      util::SetNumThreadsAndCpuBind config(comm::CommunicatorManager::getScalapackConfigInfo());
      scalapack::pzpotri_(&char_uplo, &n, a, &ia, &ja, desca, &info);
    }

    inline void pheevd(UpLo uplo, int n, float* a, int ia, int ja, int* desca, float* w, float* z,
                       int iz, int jz, int* descz, int& info) {
      const char char_uplo = static_cast<char>(uplo);
      util::SetNumThreadsAndCpuBind config(comm::CommunicatorManager::getScalapackConfigInfo());
      // Workspace query
      int lwork = -1;
      int liwork = -1;
      {
        float work;
        int iwork;
        scalapack::pssyevd_("V", &char_uplo, &n, a, &ia, &ja, desca, w, z, &iz, &jz, descz, &work,
                            &lwork, &iwork, &liwork, &info);
        lwork = static_cast<int>(work);
        liwork = iwork;
      }
      // Allocate workspace
      std::vector<float> work(lwork);
      std::vector<int> iwork(liwork);

      // Compute evs
      scalapack::pssyevd_("V", &char_uplo, &n, a, &ia, &ja, desca, w, z, &iz, &jz, descz, &work[0],
                          &lwork, &iwork[0], &liwork, &info);
    }
    inline void pheevd(UpLo uplo, int n, double* a, int ia, int ja, int* desca, double* w,
                       double* z, int iz, int jz, int* descz, int& info) {
      const char char_uplo = static_cast<char>(uplo);
      util::SetNumThreadsAndCpuBind config(comm::CommunicatorManager::getScalapackConfigInfo());
      // Workspace query
      int lwork = -1;
      int liwork = -1;
      {
        double work;
        int iwork;
        scalapack::pdsyevd_("V", &char_uplo, &n, a, &ia, &ja, desca, w, z, &iz, &jz, descz, &work,
                            &lwork, &iwork, &liwork, &info);
        lwork = static_cast<int>(work);
        liwork = iwork;
      }
      // Allocate workspace
      std::vector<double> work(lwork);
      std::vector<int> iwork(liwork);

      // Compute evs
      scalapack::pdsyevd_("V", &char_uplo, &n, a, &ia, &ja, desca, w, z, &iz, &jz, descz, &work[0],
                          &lwork, &iwork[0], &liwork, &info);
    }
    inline void pheevd(UpLo uplo, int n, std::complex<float>* a, int ia, int ja, int* desca,
                       float* w, std::complex<float>* z, int iz, int jz, int* descz, int& info) {
      const char char_uplo = static_cast<char>(uplo);
      util::SetNumThreadsAndCpuBind config(comm::CommunicatorManager::getScalapackConfigInfo());
      // Workspace query
      int lwork = -1;
      int lrwork = -1;
      int liwork = -1;
      {
        std::complex<float> work;
        float rwork;
        int iwork;
        scalapack::pcheevd_("V", &char_uplo, &n, a, &ia, &ja, desca, w, z, &iz, &jz, descz, &work,
                            &lwork, &rwork, &lrwork, &iwork, &liwork, &info);
        lwork = static_cast<int>(work.real());
        lrwork = static_cast<int>(rwork);
        liwork = iwork;
      }
      // Allocate workspace
      std::vector<std::complex<float>> work(lwork);
      std::vector<float> rwork(lrwork);
      std::vector<int> iwork(liwork);

      // Compute evs
      scalapack::pcheevd_("V", &char_uplo, &n, a, &ia, &ja, desca, w, z, &iz, &jz, descz, &work[0],
                          &lwork, &rwork[0], &lrwork, &iwork[0], &liwork, &info);
    }
    inline void pheevd(UpLo uplo, int n, std::complex<double>* a, int ia, int ja, int* desca,
                       double* w, std::complex<double>* z, int iz, int jz, int* descz, int& info) {
      const char char_uplo = static_cast<char>(uplo);
      util::SetNumThreadsAndCpuBind config(comm::CommunicatorManager::getScalapackConfigInfo());
      // Workspace query
      int lwork = -1;
      int lrwork = -1;
      int liwork = -1;
      {
        std::complex<double> work;
        double rwork;
        int iwork;
        scalapack::pzheevd_("V", &char_uplo, &n, a, &ia, &ja, desca, w, z, &iz, &jz, descz, &work,
                            &lwork, &rwork, &lrwork, &iwork, &liwork, &info);
        lwork = static_cast<int>(work.real());
        lrwork = static_cast<int>(rwork);
        liwork = iwork;
      }
      // Allocate workspace
      std::vector<std::complex<double>> work(lwork);
      std::vector<double> rwork(lrwork);
      std::vector<int> iwork(liwork);

      // Compute evs
      scalapack::pzheevd_("V", &char_uplo, &n, a, &ia, &ja, desca, w, z, &iz, &jz, descz, &work[0],
                          &lwork, &rwork[0], &lrwork, &iwork[0], &liwork, &info);
    }

    inline void ptradd(UpLo uplo, OpTrans trans, int m, int n, float alpha, const float* a, int ia,
                       int ja, int* desca, float beta, float* c, int ic, int jc, int* descc) {
      const char char_uplo = static_cast<char>(uplo);
      const char char_trans = static_cast<char>(trans);

      util::SetNumThreadsAndCpuBind config(comm::CommunicatorManager::getScalapackConfigInfo());
      scalapack::pstradd_(&char_uplo, &char_trans, &m, &n, &alpha, a, &ia, &ja, desca, &beta, c,
                          &ic, &jc, descc);
    }
    inline void ptradd(UpLo uplo, OpTrans trans, int m, int n, double alpha, const double* a, int ia,
                       int ja, int* desca, double beta, double* c, int ic, int jc, int* descc) {
      const char char_uplo = static_cast<char>(uplo);
      const char char_trans = static_cast<char>(trans);

      util::SetNumThreadsAndCpuBind config(comm::CommunicatorManager::getScalapackConfigInfo());
      scalapack::pdtradd_(&char_uplo, &char_trans, &m, &n, &alpha, a, &ia, &ja, desca, &beta, c,
                          &ic, &jc, descc);
    }
    inline void ptradd(UpLo uplo, OpTrans trans, int m, int n, std::complex<float> alpha,
                       const std::complex<float>* a, int ia, int ja, int* desca,
                       std::complex<float> beta, std::complex<float>* c, int ic, int jc, int* descc) {
      const char char_uplo = static_cast<char>(uplo);
      const char char_trans = static_cast<char>(trans);

      util::SetNumThreadsAndCpuBind config(comm::CommunicatorManager::getScalapackConfigInfo());
      scalapack::pctradd_(&char_uplo, &char_trans, &m, &n, &alpha, a, &ia, &ja, desca, &beta, c,
                          &ic, &jc, descc);
    }
    inline void ptradd(UpLo uplo, OpTrans trans, int m, int n, std::complex<double> alpha,
                       const std::complex<double>* a, int ia, int ja, int* desca,
                       std::complex<double> beta, std::complex<double>* c, int ic, int jc,
                       int* descc) {
      const char char_uplo = static_cast<char>(uplo);
      const char char_trans = static_cast<char>(trans);

      util::SetNumThreadsAndCpuBind config(comm::CommunicatorManager::getScalapackConfigInfo());
      scalapack::pztradd_(&char_uplo, &char_trans, &m, &n, &alpha, a, &ia, &ja, desca, &beta, c,
                          &ic, &jc, descc);
    }

    inline void ptrtri(UpLo uplo, Diag diag, int n, float* a, int ia, int ja, int* desca, int& info) {
      const char char_uplo = static_cast<char>(uplo);
      const char char_diag = static_cast<char>(diag);

      util::SetNumThreadsAndCpuBind config(comm::CommunicatorManager::getScalapackConfigInfo());
      scalapack::pstrtri_(&char_uplo, &char_diag, &n, a, &ia, &ja, desca, &info);
    }
    inline void ptrtri(UpLo uplo, Diag diag, int n, double* a, int ia, int ja, int* desca, int& info) {
      const char char_uplo = static_cast<char>(uplo);
      const char char_diag = static_cast<char>(diag);

      util::SetNumThreadsAndCpuBind config(comm::CommunicatorManager::getScalapackConfigInfo());
      scalapack::pdtrtri_(&char_uplo, &char_diag, &n, a, &ia, &ja, desca, &info);
    }
    inline void ptrtri(UpLo uplo, Diag diag, int n, std::complex<float>* a, int ia, int ja,
                       int* desca, int& info) {
      const char char_uplo = static_cast<char>(uplo);
      const char char_diag = static_cast<char>(diag);

      util::SetNumThreadsAndCpuBind config(comm::CommunicatorManager::getScalapackConfigInfo());
      scalapack::pctrtri_(&char_uplo, &char_diag, &n, a, &ia, &ja, desca, &info);
    }
    inline void ptrtri(UpLo uplo, Diag diag, int n, std::complex<double>* a, int ia, int ja,
                       int* desca, int& info) {
      const char char_uplo = static_cast<char>(uplo);
      const char char_diag = static_cast<char>(diag);

      util::SetNumThreadsAndCpuBind config(comm::CommunicatorManager::getScalapackConfigInfo());
      scalapack::pztrtri_(&char_uplo, &char_diag, &n, a, &ia, &ja, desca, &info);
    }
  }
}
#endif

#endif  // DLA_INTERFACE_DLA_SCALAPACK_H
