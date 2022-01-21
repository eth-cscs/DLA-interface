//
// Distributed Linear Algebra Interface (DLAI)
//
// Copyright (c) 2018-2021, ETH Zurich
// All rights reserved.
//
// Please, refer to the LICENSE file in the root directory.
// SPDX-License-Identifier: BSD-3-Clause
//

#ifndef DLA_INTERFACE_TYPES_H
#define DLA_INTERFACE_TYPES_H

#include <complex>
#include <map>
#include <mpi.h>
#include "util_macro.h"
#ifdef DLAI_WITH_DPLASMA
#include "ordered_dplasma.h"
#endif

namespace dla_interface {

  /// @note
  /// When an element is added to a enum the corresponding
  /// SET and name map (if it exists) have to be updated as well.
  enum SolverType { ScaLAPACK = 1, ELPA = 2, DPlasma = 3, Chameleon = 4, HPX_LINALG = 5 };
  MAYBE_UNUSED constexpr auto SOLVER_SET = {ScaLAPACK, ELPA, DPlasma, Chameleon, HPX_LINALG};

  const std::map<SolverType, std::string> solver_names = {{ScaLAPACK, "ScaLAPACK"},
                                                          {ELPA, "ELPA"},
                                                          {DPlasma, "DPLASMA"},
                                                          {Chameleon, "Chameleon"},
                                                          {HPX_LINALG, "HPX_LINALG"}};

  enum OpTrans { NoTrans = 'N', Trans = 'T', ConjTrans = 'C' };
  MAYBE_UNUSED constexpr auto OPTRANS_SET = {NoTrans, Trans, ConjTrans};

  enum UpLo { All = 'A', Lower = 'L', Upper = 'U' };
  MAYBE_UNUSED constexpr auto UPLO_SET = {All, Lower, Upper};

  enum Diag { Unit = 'U', NonUnit = 'N' };
  MAYBE_UNUSED constexpr auto DIAG_SET = {Unit, NonUnit};

  /// Distributed matrix distribution (see documentation for details):<br>
  /// Scalapack = 2D block cyclic distribution,<br>
  /// Tile = 2D block cyclic tile distribution.
  enum DistributionType { scalapack_dist = 1, tile_dist = 2 };
  MAYBE_UNUSED constexpr auto DISTRIBUTION_SET = {scalapack_dist, tile_dist};

  const std::map<DistributionType, std::string> dist_names = {{scalapack_dist, "ScaLAPACK"},
                                                              {tile_dist, "Tile"}};

  using IndexType = int;
  using SizeType = int;

#ifdef DLAI_WITH_SCALAPACK
  using ScalapackIndex = int;  // Base 1 indeces.
  using ScalapackDescriptor = int*;
  using constScalapackDescriptor = const int*;
  using BlacsContextType = int;
#endif
#ifdef DLAI_WITH_DPLASMA
  using ParsecContext = parsec_context_t*;
  using DPlasmaDescriptor = two_dim_block_cyclic_t;
#endif

  enum Ordering { RowMajor = 'R', ColMajor = 'C' };
  MAYBE_UNUSED constexpr auto ORDER_SET = {RowMajor, ColMajor};

  /// structure to access BaseType and ComplexType defined in the following way:<br>
  ///
  /// - ElType, BaseType: ComplexType<br>
  /// - float, float:    std::complex<float>
  /// - double, double:   std::complex<double>
  /// - std::complex<float>, float: std::complex<float>
  /// - std::complex<double>, double: std::complex<double>
  template <class ElType>
  struct TypeInfo {};

  template <>
  struct TypeInfo<int> {
    using ElementType = int;
#ifdef DLAI_WITH_DPLASMA
    static constexpr matrix_type dplasma_type = matrix_Integer;
#endif
  };
  template <>
  struct TypeInfo<float> {
    using ElementType = float;
    using BaseType = float;
    using ComplexType = std::complex<float>;
    static constexpr int ops_add = 1;
    static constexpr int ops_mult = 1;
#ifdef DLAI_WITH_DPLASMA
    static constexpr matrix_type dplasma_type = matrix_RealFloat;
#endif
  };
  template <>
  struct TypeInfo<double> {
    using ElementType = double;
    using BaseType = double;
    using ComplexType = std::complex<double>;
    static constexpr int ops_add = 1;
    static constexpr int ops_mult = 1;
#ifdef DLAI_WITH_DPLASMA
    static constexpr matrix_type dplasma_type = matrix_RealDouble;
#endif
  };
  template <>
  struct TypeInfo<std::complex<float>> {
    using ElementType = std::complex<float>;
    using BaseType = float;
    using ComplexType = std::complex<float>;
    static constexpr int ops_add = 2;
    static constexpr int ops_mult = 6;
#ifdef DLAI_WITH_DPLASMA
    static constexpr matrix_type dplasma_type = matrix_ComplexFloat;
#endif
  };
  template <>
  struct TypeInfo<std::complex<double>> {
    using ElementType = std::complex<double>;
    using BaseType = double;
    using ComplexType = std::complex<double>;
    static constexpr int ops_add = 2;
    static constexpr int ops_mult = 6;
#ifdef DLAI_WITH_DPLASMA
    static constexpr matrix_type dplasma_type = matrix_ComplexDouble;
#endif
  };

  template <class ElType>
  using BaseType = typename TypeInfo<ElType>::BaseType;

  template <class ElType>
  using ComplexType = typename TypeInfo<ElType>::ComplexType;

  const int I_ZERO = 0;
  const float S_ZERO = 0;
  const double D_ZERO = 0;
  const std::complex<float> C_ZERO(0);
  const std::complex<double> Z_ZERO(0);
#define DLA_I_ZERO (&dla_interface::I_ZERO)
#define DLA_S_ZERO (&dla_interface::S_ZERO)
#define DLA_D_ZERO (&dla_interface::D_ZERO)
#define DLA_C_ZERO (&dla_interface::C_ZERO)
#define DLA_Z_ZERO (&dla_interface::Z_ZERO)
  const int I_ONE = 1;
  const float S_ONE = 1;
  const double D_ONE = 1;
  const std::complex<float> C_ONE(1);
  const std::complex<double> Z_ONE(1);
#define DLA_I_ONE (&dla_interface::I_ONE)
#define DLA_S_ONE (&dla_interface::S_ONE)
#define DLA_D_ONE (&dla_interface::D_ONE)
#define DLA_C_ONE (&dla_interface::C_ONE)
#define DLA_Z_ONE (&dla_interface::Z_ONE)
  const int I_M_ONE = -1;
  const float S_M_ONE = -1;
  const double D_M_ONE = -1;
  const std::complex<float> C_M_ONE(-1);
  const std::complex<double> Z_M_ONE(-1);
#define DLA_I_M_ONE (&dla_interface::I_M_ONE)
#define DLA_S_M_ONE (&dla_interface::S_M_ONE)
#define DLA_D_M_ONE (&dla_interface::D_M_ONE)
#define DLA_C_M_ONE (&dla_interface::C_M_ONE)
#define DLA_Z_M_ONE (&dla_interface::Z_M_ONE)
}

#endif  // DLA_INTERFACE_TYPES_H
