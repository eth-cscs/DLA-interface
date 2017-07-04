#ifndef DLA_INTERFACE_TYPES_H
#define DLA_INTERFACE_TYPES_H

namespace dla_interface {
#ifdef DLA_HAVE_SCALAPACK
  using ScalapackDescriptor = int*;
  using BlacsContextType = int;
#endif
#ifdef HAVE_DPLASMA
  using DPlasmaDescriptor = two_dim_block_cyclic_t;
#endif

  enum Ordering { RowMajor = 'R', ColMajor = 'C' };
}

#endif  // DLA_INTERFACE_TYPES_H
