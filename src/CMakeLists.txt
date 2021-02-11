#
# Distributed Linear Algebra Interface (DLAI)
#
# Copyright (c) 2018-2021, ETH Zurich
# All rights reserved.
#
# Please, refer to the LICENSE file in the root directory.
# SPDX-License-Identifier: BSD-3-Clause
#

set(DLA_SRC
    communicator_grid.cpp
    communicator_manager.cpp
)

# C/Fortran interface enabled only if ScaLAPACK is available
if(DLA_HAVE_SCALAPACK)
list(APPEND DLA_SRC
     c_dla_interface.cpp
)
endif()
if(DLA_WITH_FORTRAN)
list(APPEND DLA_SRC
     f_dla_interface.F90
)
endif()

add_library(dla_interface ${DLA_SRC})

target_compile_definitions(dla_interface PUBLIC
  $<$<BOOL:${DLAI_WITH_ELPA}>:DLAI_WITH_ELPA>
)

target_include_directories(dla_interface PUBLIC
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_SOURCE_DIR}/libs

  # *.mod files
  $<$<BOOL:DLAI_WITH_FORTRAN>:${PROJET_BINARY_DIR}/src>

  ${DLA_HPX_LINALG_INCLUDE_DIRS}
)

# Suppress warnings for DPLASMA headers.
target_include_directories(dla_interface PUBLIC SYSTEM ${DLA_DPLASMA_INCLUDE_DIRS})
# Suppress warnings for ELPA headers.
target_include_directories(dla_interface PUBLIC SYSTEM ${ELPA_INCLUDE_DIRS})

target_link_libraries(dla_interface PUBLIC
  # backends
  ${DLA_HPX_LINALG_LIBRARIES}
  ${DLA_ELPA_LIBRARY}
  ${DLA_DPLASMA_LIBS}

  # TODO HWLOC
  ${DLA_HWLOC_LIBS}

  MPI::MPI_CXX
  $<$<BOOL:${DLAI_WITH_FORTRAN}>:MPI::MPI_Fortran>
  ${LAPACK_TARGET}
  # TODO SCALAPACK_TARGET
  ${DLA_SCALAPACK_LIBRARY}
  OpenMP::OpenMP_CXX
  $<$<BOOL:${DLAI_WITH_FORTRAN}>:OpenMP::OpenMP_Fortran>
)
