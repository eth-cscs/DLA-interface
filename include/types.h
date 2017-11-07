#ifndef DLA_INTERFACE_TYPES_H
#define DLA_INTERFACE_TYPES_H

#include <complex>
#include <mpi.h>
#ifdef DLA_HAVE_DPLASMA
#include "ordered_dplasma.h"
#endif

namespace dla_interface {
  enum SolverType { ScaLAPACK = 1, ELPA = 2, DPlasma = 3, Chameleon = 4 };

  enum UpLo { Lower = 'L', Upper = 'U' };
  // Distributed matrix distribution (see documentation for details):
  // Scalapack = 2D block cyclic distribution,
  // Tile = 2D block cyclic tile distribution.
  enum DistributionType { scalapack_dist = 1, tile_dist = 2 };

  using IndexType = int;
  using SizeType = int;

#ifdef DLA_HAVE_SCALAPACK
  using ScalapackIndex = int;  // Base 1 indeces.
  using ScalapackDescriptor = int*;
  using BlacsContextType = int;
#endif
#ifdef DLA_HAVE_DPLASMA
  using ParsecContext = parsec_context_t*;
  using DPlasmaDescriptor = two_dim_block_cyclic_t;
#endif

  enum Ordering { RowMajor = 'R', ColMajor = 'C' };

  // structure to access BaseType and ComplexType
  // defined in the following way:
  // ElType,               BaseType, ComplexType
  // float,                float,    std::complex<float>
  // double,               double,   std::complex<double>
  // std::complex<float>,  float,    std::complex<float>
  // std::complex<double>, double,   std::complex<double>
  template <class ElType>
  struct TypeInfo {};

  template <>
  struct TypeInfo<float> {
    using ElementType = float;
    using BaseType = float;
    using ComplexType = std::complex<float>;
  };
  template <>
  struct TypeInfo<double> {
    using ElementType = double;
    using BaseType = double;
    using ComplexType = std::complex<double>;
  };
  template <>
  struct TypeInfo<std::complex<float>> {
    using ElementType = std::complex<float>;
    using BaseType = float;
    using ComplexType = std::complex<float>;
  };
  template <>
  struct TypeInfo<std::complex<double>> {
    using ElementType = std::complex<double>;
    using BaseType = double;
    using ComplexType = std::complex<double>;
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
}

#endif  // DLA_INTERFACE_TYPES_H
