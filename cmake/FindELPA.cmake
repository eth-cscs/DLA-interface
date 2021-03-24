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

### Detect
find_package(PkgConfig)

if (NOT DEFINED ELPA_PACKAGE_NAME)
  message(SEND_ERROR "You should set ELPA_PACKAGE_NAME to pkg-config library name (see pkg-config --list-all | grep elpa)")
endif()

pkg_search_module(PC_ELPA ${ELPA_PACKAGE_NAME})

find_path(ELPA_INCLUDE_DIR
  NAMES elpa/elpa.h
  PATHS ${PC_ELPA_INCLUDE_DIRS}
)

find_library(ELPA_LIBRARY
  NAMES elpa elpa_openmp
  PATHS ${PC_ELPA_LIBRARY_DIRS}
)

### TEST
include(CheckFunctionExists)
include(CMakePushCheckState)

cmake_push_check_state(RESET)

set(CMAKE_REQUIRED_LIBRARIES ${ELPA_LIBRARY} ${PC_ELPA_LDFLAGS} ${PC_ELPA_LDFLAGS_OTHER})

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
  ELPA_LIBRARY
  ELPA_INCLUDE_DIR
  ELPA_CHECK
)

### CMake Target
if (ELPA_FOUND)
  if (NOT TARGET ELPA::ELPA)
    add_library(ELPA::ELPA UNKNOWN IMPORTED)
  endif()

  set_target_properties(ELPA::ELPA PROPERTIES
    IMPORTED_LOCATION ${ELPA_LIBRARY}
    INTERFACE_COMPILE_OPTIONS ${PC_ELPA_CFLAGS_OTHER}
    INTERFACE_INCLUDE_DIRECTORIES ${ELPA_INCLUDE_DIR}
    INTERFACE_LINK_LIBRARIES "${PC_ELPA_LDFLAGS};${PC_ELPA_LDFLAGS_OTHER}"
  )
endif()
