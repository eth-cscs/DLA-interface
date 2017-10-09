#ifndef DLA_INTERFACE_TYPES_H
#define DLA_INTERFACE_TYPES_H

namespace dla_interface {
  enum SolverType { scalapack = 1, elpa = 2, dplasma = 3, chameleon = 4 };

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
  using DPlasmaDescriptor = two_dim_block_cyclic_t;
#endif

  enum Ordering { RowMajor = 'R', ColMajor = 'C' };

  const int I_ZERO = 0;
  const float S_ZERO = 0;
  const double D_ZERO = 0;
#define DLA_I_ZERO (&dla_interface::I_ZERO)
#define DLA_S_ZERO (&dla_interface::S_ZERO)
#define DLA_D_ZERO (&dla_interface::D_ZERO)
  const int I_ONE = 1;
  const float S_ONE = 1;
  const double D_ONE = 1;
#define DLA_I_ONE (&dla_interface::I_ONE)
#define DLA_S_ONE (&dla_interface::S_ONE)
#define DLA_D_ONE (&dla_interface::D_ONE)
}

#endif  // DLA_INTERFACE_TYPES_H
