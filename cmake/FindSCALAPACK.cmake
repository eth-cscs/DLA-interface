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

cmake_minimum_required(VERSION 3.12)

### Detect

find_library(SCALAPACK_LIBRARY
  NAMES
    scalapack # netlib-scalapack
)

mark_as_advanced(
  SCALAPACK_LIBRARY
)

# ===== TEST
include(CMakePushCheckState)
cmake_push_check_state(RESET)

set(CMAKE_REQUIRED_LIBRARIES ${SCALAPACK_LIBRARY})

include(CheckFunctionExists)

unset(_SCALAPACK_CHECK CACHE)
check_function_exists(pdpotrf_ _SCALAPACK_CHECK)

unset(_SCALAPACK_CHECK_BLACS CACHE)
check_function_exists(Cblacs_exit _SCALAPACK_CHECK_BLACS)

cmake_pop_check_state()

### Package
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SCALAPACK DEFAULT_MSG
  _SCALAPACK_CHECK_BLACS
  _SCALAPACK_CHECK
)

if (SCALAPACK_FOUND)
  if (NOT TARGET SCALAPACK::SCALAPACK)
    add_library(SCALAPACK::SCALAPACK INTERFACE IMPORTED)
  endif()

  if (SCALAPACK_LIBRARY)
    set_target_properties(SCALAPACK::SCALAPACK PROPERTIES
      IMPORTED_LOCATION "${SCALAPACK_LIBRARY}"
    )
  endif()
endif()
