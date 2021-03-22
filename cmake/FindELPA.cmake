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
# This module finds the ELPA library specified with ELPA_VERSION variable with pkg-config.
# If you don't know how to set the ELPA_VERSION variable, please check the output of
#
# pkg-config --list-all | grep elpa
#
# and pick the one you want.
# If found, this module create the CMake target ELPA::ELPA

### REQUIREMENTS
if (NOT SCALAPACK_FOUND)
  message(FATAL_ERROR "SCALAPACK is a requirement")
endif()

include(FindPackageHandleStandardArgs)
find_package(PkgConfig)

if (NOT DEFINED ELPA_VERSION)
  message(SEND_ERROR "You should set ELPA_VERSION to pkg-config library name (see pkg-config --list-all | grep elpa)")
endif()

pkg_search_module(ELPA ${ELPA_VERSION})

### TEST
include(CheckSymbolExists)
include(CMakePushCheckState)

cmake_push_check_state(RESET)

set(CMAKE_REQUIRED_INCLUDES ${ELPA_INCLUDE_DIRS})
set(CMAKE_REQUIRED_LIBRARIES ${ELPA_LINK_LIBRARIES} SCALAPACK::SCALAPACK)

unset(ELPA_CHECK CACHE)
check_symbol_exists(elpa_allocate "elpa/elpa.h" ELPA_CHECK)

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
