#
# Distributed Linear Algebra Interface (DLAI)
#
# Copyright (c) 2018-2021, ETH Zurich
# All rights reserved.
#
# Please, refer to the LICENSE file in the root directory.
# SPDX-License-Identifier: BSD-3-Clause
#

# Find SCALAPACK library
#
# This module finds an installed library that implements the SCALAPACK linear-algebra interface.
#
# This module sets the following variables:
#  SCALAPACK_FOUND - set to true if a library implementing the SCALAPACK interface is found
#
# Following options are allowed (for setting options from code see next section):
#   SCALAPACK_TYPE - it can be "Compiler" or "Custom"
#     - Compiler (Default): The compiler add the scalapack flag automatically therefore no
#                           extra link line has to be added.
#     - Custom: User can specify include folders and libraries through
#   SCALAPACK_INCLUDE_DIR - used if SCALAPACK_TYPE=Custom
#       ;-list of include paths
#   SCALAPACK_LIBRARY - used if SCALAPACK_TYPE=Custom
#       ;-list of {lib name, lib filepath, -Llibrary_folder}
#
# If you want to set options from CMake script you cannot use previous variables, instead you should use:
#   - SCALAPACK_CUSTOM_TYPE
#   - SCALAPACK_CUSTOM_INCLUDE_DIR
#   - SCALAPACK_CUSTOM_LIBRARY
#
# It creates target SCALAPACK::SCALAPACK

# ===== Requirements
if (NOT MPI_FOUND OR NOT LAPACK_FOUND)
  message(FATAL_ERROR "MPI and LAPACK are requirements for ScaLAPACK")
endif()

# ===== Detection
set(SCALAPACK_TYPE_OPTIONS "Compiler" "Custom")
set(SCALAPACK_TYPE "Compiler" CACHE STRING "")
set_property(CACHE SCALAPACK_TYPE PROPERTY STRINGS ${SCALAPACK_TYPE_OPTIONS})

if(SCALAPACK_TYPE STREQUAL "Custom")
  # user specifies values with
  # SCALAPACK_INCLUDE_DIR
  # SCALAPACK_LIBRARY
  if (NOT (DEFINED SCALAPACK_INCLUDE_DIR AND DEFINED SCALAPACK_LIBRARY))
    message(WARNING
      "You have not specified both variables for ScaLAPACK. "
      "SCALAPACK_INCLUDE_DIR, SCALAPACK_LIBRARY")
  endif()
elseif(SCALAPACK_TYPE STREQUAL "Compiler")
  if (SCALAPACK_INCLUDE_DIR OR SCALAPACK_LIBRARY)
    message(FATAL_ERROR "You are not supposed to specify anything if SCALAPACK_TYPE=${SCALAPACK_TYPE}")
  endif()
else()
  message(FATAL_ERROR
    "Unknown ScaLAPACK type: ${SCALAPACK_TYPE}. "
    "Available options: ${SCALAPACK_TYPE_OPTIONS}")
endif()

mark_as_advanced(
  SCALAPACK_TYPE
  SCALAPACK_INCLUDE_DIR
  SCALAPACK_LIBRARY
)

# ===== TEST
include(CMakePushCheckState)
cmake_push_check_state(RESET)

include(CheckFunctionExists)

set(CMAKE_REQUIRED_LIBRARIES MPI::MPI_CXX LAPACK::LAPACK ${SCALAPACK_LIBRARY})

unset(_SCALAPACK_CHECK CACHE)
check_function_exists(pdpotrf_ _SCALAPACK_CHECK)

unset(_SCALAPACK_CHECK_BLACS CACHE)
check_function_exists(Cblacs_exit _SCALAPACK_CHECK_BLACS)

cmake_pop_check_state()

### Package
if (SCALAPACK_TYPE STREQUAL "Compiler")
  set(SCALAPACK_FOUND TRUE)
else()
  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(SCALAPACK DEFAULT_MSG
    SCALAPACK_LIBRARY
    _SCALAPACK_CHECK_BLACS
    _SCALAPACK_CHECK
  )
endif()

if (SCALAPACK_FOUND)
  if (NOT TARGET SCALAPACK::SCALAPACK)
    add_library(SCALAPACK::SCALAPACK INTERFACE IMPORTED GLOBAL)
  endif()

  target_include_directories(SCALAPACK::SCALAPACK
    INTERFACE ${SCALAPACK_INCLUDE_DIR})
  target_link_libraries(SCALAPACK::SCALAPACK
    INTERFACE ${SCALAPACK_LIBRARY})
endif()
