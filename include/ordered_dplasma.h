#ifndef DLA_INTERFACE_ORDERED_DPLASMA_H
#define DLA_INTERFACE_ORDERED_DPLASMA_H

// This files includes the dplasma header.
// It is required since some internal includes in dplasma are inside
// extern "C" regions.

#include <complex>
#include <tuple>
#include <mpi.h>

#define operator op
#include "dplasma.h"
#include "data_dist/matrix/two_dim_rectangle_cyclic.h"
#undef operator

// complex.h (included by core_blas.h included by dplasma.h)
// defines the macro I as the complex `i`.
// This macro breaks gtest.
#ifdef I
#undef I
#endif

#endif  // DLA_INTERFACE_ORDERED_DPLASMA_H
