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
# TODO fix this
if (HWLOC_ROOT)
  set(CMAKE_PREFIX_PATH ${HWLOC_ROOT})
  find_package(PkgConfig REQUIRED)
  pkg_search_module(HWLOC REQUIRED hwloc)
  unset(CMAKE_PREFIX_PATH)
  set(DLA_HWLOC_LIBS "-L${HWLOC_LIBDIR};${HWLOC_LIBRARIES}" CACHE STRING "Libraries for HWLOC" FORCE)
  add_definitions(-DDLA_HAVE_HWLOC)
else()
  message(WARNING "Configurations without HWLOC may have bad performance due to bad thread binding. Please refers to DLAI documentation.")
endif()

# ----- MPI
find_package(MPI REQUIRED)

# ----- LAPACK/SCALAPACK
set(LAPACK_TARGET "")
set(DLAI_WITH_MKL OFF CACHE BOOL "Enable MKL as provider for LAPACK")
if (DLAI_WITH_MKL)
  set(MKL_CUSTOM_THREADING "Sequential")
  find_package(MKL REQUIRED COMPONENTS BLAS_32BIT_SEQ)

  set(LAPACK_TARGET "mkl::blas_32bit_seq")
else()
  set(LAPACK_TARGET "lapack::lapack")
  find_package(LAPACK REQUIRED)
endif()

if (DLAI_WITH_FORTRAN)
  if (NOT SCALAPACK_FOUND)
    message(FATAL_ERROR "Cannot build Fortran bindings without Scalapack. Set DLA_WITH_FORTRAN to OFF")
  endif()
endif()

# ----- OpenMP
find_package(OpenMP REQUIRED)

# ----- EXTERNAL
add_subdirectory(external)

# ----- MISC
# TODO fix these
if (DLAI_WITH_DPLASMA)
  find_package(DPLASMA REQUIRED)
endif()

if (DLAI_WITH_HPXLINALG)
  find_package(HPXLINALG REQUIRED)
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
# extra options
# ---------------------------------------------------------------------------
# TODO fix this
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-missing-braces")
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 7)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-dangling-else")
  endif()
endif()
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -ffree-line-length-0")


# ---------------------------------------------------------------------------
# DLAI library
# ---------------------------------------------------------------------------
add_subdirectory(src)
set(DLA_LIBS
    dla_interface
)

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
