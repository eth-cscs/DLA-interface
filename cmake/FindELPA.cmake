#
# Distributed Linear Algebra Interface (DLAI)
#
# Copyright (c) 2018-2021, ETH Zurich
# All rights reserved.
#
# Please, refer to the LICENSE file in the root directory.
# SPDX-License-Identifier: BSD-3-Clause
#

# Find ELPA library
#
# This module finds the ELPA library specified with ELPA_PACKAGE_NAME variable with pkg-config.
# If you don't know how to set the ELPA_PACKAGE_NAME variable, please check the output of
#
# pkg-config --list-all | grep elpa
#
# and pick the one you want.
# If found, this module create the CMake target ELPA::ELPA

### Requirements
include(CheckLanguage)

check_language(Fortran)
if (CMAKE_Fortran_COMPILER)
  enable_language(Fortran)
else()
  message(FATAL_ERROR "Fortran is a requirement for ELPA")
endif()

find_package(SCALAPACK REQUIRED QUIET)
find_package(MPI REQUIRED QUIET COMPONENTS Fortran)

### Detect
include(FindPackageHandleStandardArgs)
find_package(PkgConfig)

if (NOT DEFINED ELPA_PACKAGE_NAME)
  message(SEND_ERROR "You should set ELPA_PACKAGE_NAME to pkg-config library name (see pkg-config --list-all | grep elpa)")
endif()

pkg_search_module(ELPA ${ELPA_PACKAGE_NAME})

### TEST
include(CheckFunctionExists)
include(CMakePushCheckState)

cmake_push_check_state(RESET)

set(CMAKE_REQUIRED_INCLUDES ${ELPA_INCLUDE_DIRS})
set(CMAKE_REQUIRED_LIBRARIES ${ELPA_LINK_LIBRARIES} MPI::MPI_Fortran SCALAPACK::SCALAPACK)

unset(ELPA_CHECK CACHE)
# Note:
# If the project does not enable the C language, this check_symbol_exists may fail because the compiler,
# by looking at the file extension of the test, may decide to build it as CXX and not as C.
# For this reason, here it is just checkde that the symbol is available at linking.
check_function_exists(elpa_allocate ELPA_CHECK)

cmake_pop_check_state()

### Package
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ELPA DEFAULT_MSG
  ELPA_LINK_LIBRARIES
  ELPA_INCLUDE_DIRS
  ELPA_CHECK
)

### CMake Target
if (ELPA_FOUND)
  if (NOT TARGET ELPA::ELPA)
    add_library(ELPA::ELPA INTERFACE IMPORTED GLOBAL)
  endif()

  target_include_directories(ELPA::ELPA
    INTERFACE ${ELPA_INCLUDE_DIRS})
  target_link_libraries(ELPA::ELPA
    INTERFACE ${ELPA_LINK_LIBRARIES})
endif()
