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
option(DLAI_HAVE_DPLASMA "DPLASMA" OFF)
option(DLAI_HAVE_ELPA "ELPA" OFF)
option(DLAI_HAVE_HPXLINALG "HPXLINALG" OFF)

# ----- other
# TODO fix these
option(DLA_PRINT_DEBUG_INFO "Enables extra output for debugging purpose." OFF)
option(DLA_PRINT_DEBUG_CALL_PARAM "Enables call parameters output for debugging purpose." OFF)
option(DLA_COVERAGE_TEST "Enables coverage test." OFF)

option(DLA_ALL_TESTS_USE_RUNNER "Specify if the single rank tests has to use the runner as well" OFF)

# ---------------------------------------------------------------------------
# Languages
# ---------------------------------------------------------------------------
# ----- Fortran
if (DLAI_WITH_FORTRAN)
  enable_language(Fortran)
endif()

# ----- C
enable_language(C)

# ----- CXX
enable_language(CXX)
set(CMAKE_CXX_STANDARD 14)

# ---------------------------------------------------------------------------
# Libraries configuration
# ---------------------------------------------------------------------------
# ----- HWLOC
find_package(HWLOC)

# ----- MPI
find_package(MPI REQUIRED)

# ----- LAPACK/SCALAPACK
if (DLAI_WITH_MKL)
  find_package(MKL REQUIRED)

  set(LAPACK_TARGET "mkl::mkl_intel_32bit_seq_dyn")
  set(SCALAPACK_TARGET "mkl::scalapack_mpich_intel_32bit_seq_dyn")

  if (TARGET ${LAPACK_TARGET})
    add_library(LAPACK::LAPACK ALIAS ${LAPACK_TARGET})
  else()
    message(FATAL_ERROR "not available")
  endif()

  if (TARGET ${SCALAPACK_TARGET})
    set(SCALAPACK_FOUND TRUE)
    add_library(SCALAPACK::SCALAPACK ALIAS ${SCALAPACK_TARGET})
  endif()
else()
  if (LAPACK_TYPE OR LAPACK_INCLUDE_DIR OR LAPACK_LIBRARY)
    include(lapack)
  else()
    find_package(LAPACK REQUIRED)
  endif()

  find_package(SCALAPACK)
endif()

if (DLAI_WITH_FORTRAN)
  if (NOT SCALAPACK_FOUND)
    message(FATAL_ERROR "Cannot build Fortran bindings without Scalapack. Set DLAI_WITH_FORTRAN to OFF")
  endif()
endif()

# ----- OpenMP
find_package(OpenMP REQUIRED)

# ----- EXTERNAL
add_subdirectory(external)

# ----- BACKENDs
if (DLAI_WITH_DLAF)
  find_package(DLAF REQUIRED)
endif()

# TODO fix these
if (DLAI_WITH_DPLASMA)
  find_package(DPLASMA REQUIRED)
endif()

if (DLAI_WITH_ELPA)
  find_package(ELPA REQUIRED)
endif()

# ---- setup
# TODO fix these
if(DLA_PRINT_DEBUG_INFO)
  add_definitions("-DDLA_THREAD_DEBUG_INFO")
endif()
if(DLA_PRINT_DEBUG_CALL_PARAM)
  add_definitions("-DDLA_PRINT_DEBUG_CALL_PARAM")
endif()
if(DLA_COVERAGE_TEST)
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
    set(CMAKE_LD_FLAGS "${CMAKE_LD_FLAGS} --coverage")

    set(DLA_COVERAGE_SRC_DIRS
        ${CMAKE_SOURCE_DIR}/include/*
        ${CMAKE_SOURCE_DIR}/include/*/*
        ${CMAKE_SOURCE_DIR}/src/*
        ${CMAKE_SOURCE_DIR}/test/include/*
        ${CMAKE_SOURCE_DIR}/test/src/*
        ${CMAKE_SOURCE_DIR}/test/unit/*
        )

    add_custom_target("coverage"
      COMMAND lcov --capture --directory ${CMAKE_BINARY_DIR} --output-file coverage.info
      COMMAND lcov --extract coverage.info ${DLA_COVERAGE_SRC_DIRS} --output-file coverage-filt.info
      COMMAND genhtml coverage-filt.info --output-directory ${CMAKE_BINARY_DIR}/coverage
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

  else ()
    message(WARNING "Coverage test can be used only with GCC.")
    set(DLA_COVERAGE_TEST "OFF" FORCE)
  endif()
endif()

# ---------------------------------------------------------------------------
# DLAI library
# ---------------------------------------------------------------------------
add_subdirectory(src)

# ---------------------------------------------------------------------------
# mini Apps
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
