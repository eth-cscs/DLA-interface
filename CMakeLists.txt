#
# Distributed Linear Algebra Interface (DLA)
#
# Copyright (c) 2018-2019, ETH Zurich
# All rights reserved.
#
# Please, refer to the LICENSE file in the root directory.
# SPDX-License-Identifier: BSD-3-Clause
#
cmake_minimum_required(VERSION 3.8)
project(DLA_interface CXX C)

set(CMAKE_CXX_STANDARD 14)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# ---------------------------------------------------------------------------
# CMake module path
# ---------------------------------------------------------------------------
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

include(dla_utils)
include(dla_lapack)
include(dla_scalapack)
include(dla_dplasma)
include(dla_hpx_linalg)
include(dla_hwloc)
include(dla_elpa)

setoption(DLA_WITH_FORTRAN BOOL ON "Build fortran binding")
if (${DLA_WITH_FORTRAN})
  enable_language(Fortran)
endif()

include(dla_openmp)
include(dla_mpi)

dla_find_lapack()

dla_find_scalapack()

dla_find_dplasma()

dla_find_hpx_linalg()

dla_find_hwloc()

dla_find_elpa()

if (${DLA_WITH_FORTRAN})
  if (NOT ${DLA_HAVE_SCALAPACK})
    message(FATAL_ERROR "Cannot build Fortran bindings without Scalapack. Set DLA_WITH_FORTRAN to OFF")
  endif()
endif()

setoption(DLA_PRINT_DEBUG_INFO BOOL "OFF" "Enables extra output for debugging purpose.")
if(DLA_PRINT_DEBUG_INFO)
  add_definitions("-DDLA_THREAD_DEBUG_INFO")
endif()
setoption(DLA_PRINT_DEBUG_CALL_PARAM BOOL "OFF" "Enables call parameters output for debugging purpose.")
if(DLA_PRINT_DEBUG_CALL_PARAM)
  add_definitions("-DDLA_PRINT_DEBUG_CALL_PARAM")
endif()
setoption(DLA_COVERAGE_TEST BOOL "OFF" "Enables coverage test.")
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

include(CTest)
enable_testing()
setoption(DLA_ALL_TESTS_USE_RUNNER BOOL "OFF" "Specify if the single rank tests has to use the runner as well")

include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/libs
    ${MPI_CXX_INCLUDE_PATH}
    ${DLA_LAPACK_INCLUDE_PATH}
    ${DLA_HPX_LINALG_INCLUDE_DIRS}
    ${HWLOC_INCLUDE_DIRS}
)
if (${DLA_WITH_FORTRAN})
  # .mod files
  include_directories(
    ${MPI_Fortran_INCLUDE_PATH}
    ${PROJECT_BINARY_DIR}/src
  )
endif()
# Suppress warnings for DPLASMA headers.
#include_directories(SYSTEM ${DLA_DPLASMA_INCLUDE_DIRS})
# Suppress warnings for ELPA headers.
#include_directories(SYSTEM ${DLA_ELPA_INCLUDE_DIRS})

# ---------------------------------------------------------------------------
# extra options
# ---------------------------------------------------------------------------
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-missing-braces")
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 7)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-dangling-else")
  endif()
endif()
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -ffree-line-length-0")

# ---------------------------------------------------------------------------
# setup libraries
# ---------------------------------------------------------------------------
add_subdirectory(src)

set(DLA_LIBS
  dla_interface
  ${DLA_HPX_LINALG_LIBRARIES}
  ${DLA_ELPA_LIBRARY}
  ${DLA_DPLASMA_LIBS}
  ${DLA_SCALAPACK_LIBRARY}
  ${DLA_LAPACK_LIBRARY}
  ${MPI_Fortran_LIBRARIES}
  ${MPI_CXX_LIBRARIES}
  ${DLA_HWLOC_LIBS}
)

# ---------------------------------------------------------------------------
# Add sub directories
# ---------------------------------------------------------------------------
add_subdirectory(miniapp)

# ---------------------------------------------------------------------------
# Test
# ---------------------------------------------------------------------------
set(DLA_WITH_TEST ON CACHE BOOL "Enables tests")
if (DLA_WITH_TEST)
  include(CTest)
  enable_testing()

  add_subdirectory(test)
endif()

# ---------------------------------------------------------------------------
# Internal Documetation - Doxygen
# ---------------------------------------------------------------------------
# ***!!! By default should be OFF, for testing it is set to ON !!!***
# ---------------------------------------------------------------------------
#
#option(DLA_INTERNAL_DOCUMENTATION "Build internal documentation" OFF)
option(DLA_INTERNAL_DOCUMENTATION "Build internal documentation" ON)

if (DLA_INTERNAL_DOCUMENTATION)
  add_subdirectory(doc)
endif()

