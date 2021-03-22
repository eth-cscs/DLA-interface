#
# Distributed Linear Algebra Interface (DLAI)
#
# Copyright (c) 2018-2021, ETH Zurich
# All rights reserved.
#
# Please, refer to the LICENSE file in the root directory.
# SPDX-License-Identifier: BSD-3-Clause
#

add_library(DLAI
  communicator_grid.cpp
  communicator_manager.cpp

  # C/Fortran interface enabled only if ScaLAPACK is available
  $<$<BOOL:${SCALAPACK_FOUND}>:c_dla_interface.cpp>

  $<$<BOOL:${DLAI_WITH_FORTRAN}>:f_dla_interface.F90>
)

target_compile_options(DLAI PUBLIC
  $<$<COMPILE_LANGUAGE:Fortran>:-ffree-line-length-0>
)

target_compile_definitions(DLAI PUBLIC
  $<$<BOOL:${SCALAPACK_FOUND}>:DLAI_WITH_SCALAPACK>
  $<$<BOOL:${HWLOC_FOUND}>:DLAI_WITH_HWLOC>

  $<$<BOOL:${DLAI_WITH_DLAF}>:DLAI_WITH_DLAF>
  $<$<BOOL:${DLAI_WITH_DPLASMA}>:DLAI_WITH_DPLASMA>
  $<$<BOOL:${DLAI_WITH_ELPA}>:DLAI_WITH_ELPA>
)

target_include_directories(DLAI PUBLIC
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_SOURCE_DIR}/libs

  # *.mod files
  $<$<BOOL:${DLAI_WITH_FORTRAN}>:${CMAKE_CURRENT_BINARY_DIR}>
)

target_link_libraries(DLAI PUBLIC
  # backends
  $<$<BOOL:${DLAI_WITH_DLAF}>:DLAF::DLAF>
  $<$<BOOL:${DLAI_WITH_DPLASMA}>:DPlasma::dplasma>
  $<$<BOOL:${DLAI_WITH_ELPA}>:ELPA::ELPA>

  $<$<BOOL:${HWLOC_FOUND}>:HWLOC::HWLOC>

  MPI::MPI_CXX
  $<$<BOOL:${DLAI_WITH_FORTRAN}>:MPI::MPI_Fortran>

  LAPACK::LAPACK
  $<$<BOOL:${SCALAPACK_FOUND}>:SCALAPACK::SCALAPACK>

  OpenMP::OpenMP_CXX
  $<$<BOOL:${DLAI_WITH_FORTRAN}>:OpenMP::OpenMP_Fortran>
)
