#
# Distributed Linear Algebra Interface (DLAI)
#
# Copyright (c) 2018-2021, ETH Zurich
# All rights reserved.
#
# Please, refer to the LICENSE file in the root directory.
# SPDX-License-Identifier: BSD-3-Clause
#

cmake_minimum_required(VERSION 3.14)

project(DLAI)

# ---------------------------------------------------------------------------
# CMake configurations
# ---------------------------------------------------------------------------
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

include(CTest)
include(DLAI)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Build Type" FORCE)
endif()

option(DLAI_WITH_MKL "Enable MKL as provider for LAPACK" OFF)
option(DLAI_WITH_FORTRAN "Build with Fortran bindings" ON)
option(DLAI_BUILD_MINIAPPS "Build miniapps" ON)
option(DLAI_BUILD_TESTING "Build tests" ON)

# ----- BACKENDs
option(DLAI_WITH_DPLASMA "Enable DPLASMA backend" OFF)
option(DLAI_WITH_ELPA "Enable ELPA backend" OFF)
# TODO this is going to be replaced by DLAF
option(DLAI_WITH_HPXLINALG "Enable HPXLINALG backend" OFF)

# ----- others
option(DLAI_PRINT_DEBUG_INFO "Enables extra output for debugging purpose." OFF)
option(DLAI_PRINT_DEBUG_CALL_PARAM "Enables call parameters output for debugging purpose." OFF)

# ---------------------------------------------------------------------------
# Languages
# ---------------------------------------------------------------------------
enable_language(C)

enable_language(CXX)
set(CMAKE_CXX_STANDARD 14)

if (DLAI_WITH_FORTRAN)
  enable_language(Fortran)
endif()

# ---------------------------------------------------------------------------
# Libraries configuration
# ---------------------------------------------------------------------------

find_package(HWLOC)
find_package(OpenMP REQUIRED)

find_package(MPI REQUIRED)
mpi_preset_setup()

# ----- LAPACK/SCALAPACK
if (DLAI_WITH_MKL)
  find_package(MKL REQUIRED)

  set(MKL_LAPACK_TARGET "mkl::mkl_intel_32bit_seq_dyn"
    CACHE STRING "MKL LAPACK target (see FindMKL for details)")
  set(MKL_SCALAPACK_TARGET "mkl::scalapack_mpich_intel_32bit_seq_dyn"
    CACHE STRING "MKL ScaLAPACK target (see FindMKL for details)")

  # Note:
  # If specified targets are found, "standardize" them

  if (TARGET ${MKL_LAPACK_TARGET})
    set(LAPACK_FOUND TRUE)
    add_library(LAPACK::LAPACK INTERFACE IMPORTED GLOBAL)
    target_link_libraries(LAPACK::LAPACK INTERFACE ${MKL_LAPACK_TARGET})
  else()
    message(FATAL_ERROR "not available")
  endif()

  if (TARGET ${MKL_SCALAPACK_TARGET})
    set(SCALAPACK_FOUND TRUE)
    add_library(SCALAPACK::SCALAPACK INTERFACE IMPORTED GLOBAL)
    target_link_libraries(SCALAPACK::SCALAPACK INTERFACE ${MKL_SCALAPACK_TARGET})
  endif()
else()
  find_package(LAPACK REQUIRED)
  find_package(SCALAPACK REQUIRED)
endif()

# ----- BACKENDs
if (DLAI_WITH_DLAF)
  find_package(DLAF REQUIRED)
endif()

if (DLAI_WITH_DPLASMA)
  find_package(DPLASMA REQUIRED)
endif()

if (DLAI_WITH_ELPA)
  find_package(ELPA REQUIRED)
endif()

# ----- EXTERNAL
add_subdirectory(external)

# ---------------------------------------------------------------------------
# DLAI library
# ---------------------------------------------------------------------------
add_subdirectory(src)

# ---------------------------------------------------------------------------
# Miniapps
# ---------------------------------------------------------------------------
if (DLAI_BUILD_MINIAPPS)
  add_subdirectory(miniapp)
endif()

# ---------------------------------------------------------------------------
# Test
# ---------------------------------------------------------------------------
if (DLAI_BUILD_TESTING)
  add_subdirectory(test)
endif()
