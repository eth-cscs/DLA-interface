#
# Distributed Linear Algebra with Future (DLAF)
#
# Copyright (c) 2018-2021, ETH Zurich
# All rights reserved.
#
# Please, refer to the LICENSE file in the root directory.
# SPDX-License-Identifier: BSD-3-Clause
#

# Find LAPACK library
#
# This module finds an installed library that implements the LAPACK linear-algebra interface.
#
# This module sets the following variables:
#  LAPACK_FOUND - set to true if a library implementing the LAPACK interface is found
#
# Following options are allowed (for setting options from code see next section):
#   LAPACK_TYPE - it can be "Compiler" or "Custom"
#     - Compiler (Default): The compiler add the scalapack flag automatically therefore no
#                           extra link line has to be added.
#     - Custom: User can specify include folders and libraries through
#   LAPACK_INCLUDE_DIR - used if SCALAPACK_TYPE=Custom
#       ;-list of include paths
#   LAPACK_LIBRARY - used if SCALAPACK_TYPE=Custom
#       ;-list of {lib name, lib filepath, -Llibrary_folder}
#
# If you want to set options from CMake script you cannot use previous variables, instead you should use:
#   - LAPACK_CUSTOM_TYPE
#   - LAPACK_CUSTOM_INCLUDE_DIR
#   - LAPACK_CUSTOM_LIBRARY
#
# It creates target LAPACK::LAPACK

set(LAPACK_TYPE_OPTIONS "Compiler" "Custom")
set(LAPACK_TYPE "Compiler" CACHE STRING "BLAS/LAPACK type setting")
set_property(CACHE LAPACK_TYPE PROPERTY STRINGS ${LAPACK_TYPE_OPTIONS})

if(LAPACK_TYPE STREQUAL "Custom" OR LAPACK_INCLUDE_DIR OR LAPACK_LIBRARY)
  # user specifies values with
  # LAPACK_INCLUDE_DIR
  # LAPACK_LIBRARY
  if (NOT (DEFINED LAPACK_INCLUDE_DIR AND DEFINED LAPACK_LIBRARY))
    message(WARNING
      "You have not specified both variables for LAPACK. "
      "LAPACK_INCLUDE_DIR, LAPACK_LIBRARY")
  endif()
elseif(LAPACK_TYPE STREQUAL "Compiler")
  # using the compiler wrappers
else()
  message(FATAL_ERROR
    "Unknown LAPACK type: ${LAPACK_TYPE}. "
    "Available options: ${LAPACK_TYPE_OPTIONS}")
endif()

mark_as_advanced(
  LAPACK_TYPE
  LAPACK_INCLUDE_DIR
  LAPACK_LIBRARY
)

### Checks
include(CMakePushCheckState)
cmake_push_check_state(RESET)

include(CheckFunctionExists)

set(CMAKE_REQUIRED_INCLUDES ${LAPACK_INCLUDE_DIR})
set(CMAKE_REQUIRED_LIBRARIES ${LAPACK_LIBRARY})

unset(_LAPACK_CHECK_BLAS CACHE)
check_function_exists(dgemm_ _LAPACK_CHECK_BLAS)

unset(_LAPACK_CHECK CACHE)
check_function_exists(dpotrf_ _LAPACK_CHECK)

cmake_pop_check_state()

### Package
if (LAPACK_TYPE STREQUAL "Compiler")
  set(LAPACK_FOUND TRUE)
else()
  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(LAPACK DEFAULT_MSG
    LAPACK_LIBRARY
    _LAPACK_CHECK_BLAS
    _LAPACK_CHECK
  )
endif()

if (LAPACK_FOUND)
  set(LAPACK_INCLUDE_DIRS ${LAPACK_INCLUDE_DIR})
  set(LAPACK_LIBRARIES ${LAPACK_LIBRARY})

  if (NOT TARGET LAPACK::LAPACK)
    add_library(LAPACK::LAPACK INTERFACE IMPORTED GLOBAL)
  endif()

  set_target_properties(LAPACK::LAPACK
    PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${LAPACK_INCLUDE_DIRS}"
      INTERFACE_LINK_LIBRARIES "${LAPACK_LIBRARIES}"
  )
endif()
