#
# CMake recipes
#
# Copyright (c) 2018-2021, ETH Zurich
# BSD 3-Clause License. All rights reserved.
#
# author: Alberto Invernizzi (a.invernizzi@cscs.ch)
#

# Find SCALAPACK library
#
# !!! WARNING !!!
# It is up to the user to ensure the compatibility between the ScaLAPACK library and any
# other BLAS/LAPACK and MPI implementations used in the project.
#
# This module checks if `netlib-scalapack` is available or if the compiler implicitly links
# to a ScaLAPACK implementation.
#
# Users can manually specify next variables (even by setting them empty to force use of
# the compiler implicit linking) to control which implementation they want to use:
#   SCALAPACK_LIBRARY
#       ;-list of {lib name, lib filepath, -Llibrary_folder}
#
# This module sets the following variables:
#   SCALAPACK_FOUND - set to true if a library implementing the SCALAPACK interface is found
#
# If ScaLAPACK symbols got found, it creates target SCALAPACK::SCALAPACK

cmake_minimum_required(VERSION 3.12)

macro(_scalapack_find_dependency dep)
  set(_scalapack_quiet_arg)
  if(SCALAPACK_FIND_QUIETLY)
    set(_scalapack_quiet_arg QUIET)
  endif()
  set(_scalapack_required_arg)
  if(SCALAPACK_FIND_REQUIRED)
    set(_scalapack_required_arg REQUIRED)
  endif()
  find_package(${dep} ${ARGN}
    ${_scalapack_quiet_arg}
    ${_scalapack_required_arg})
  if (NOT ${dep}_FOUND)
    set(SCALAPACK_NOT_FOUND_MESSAGE "SCALAPACK could not be found because dependency ${dep} could not be found.")
  endif()

  set(_scalapack_required_arg)
  set(_scalapack_quiet_arg)
endmacro()

macro(_scalapack_check_symbols _PREFIX)
  include(CMakePushCheckState)
  cmake_push_check_state(RESET)

  set(_libraries ${${_PREFIX}_LIBRARY})

  set(CMAKE_REQUIRED_QUIET TRUE)
  set(CMAKE_REQUIRED_LIBRARIES "${_libraries}")

  include(CheckFunctionExists)

  unset(_SCALAPACK_CHECK CACHE)
  check_function_exists(pdpotrf_ _SCALAPACK_CHECK)

  unset(_SCALAPACK_CHECK_BLACS CACHE)
  check_function_exists(Cblacs_exit _SCALAPACK_CHECK_BLACS)

  cmake_pop_check_state()

  if (_SCALAPACK_CHECK AND _SCALAPACK_CHECK_BLACS)
    if (_libraries)
      set(SCALAPACK_LIBRARY ${_libraries})
    else()
      set(SCALAPACK_LIBRARY "SCALAPACK-PLACEHOLDER-FOR-EMPTY-LIBRARIES")
    endif()
  else()
    set(SCALAPACK_LIBRARY FALSE)
  endif()

  unset(_libraries)
endmacro()


_scalapack_find_dependency(LAPACK)

if (NOT SCALAPACK_NOT_FOUND_MESSAGE)
  # check if compiler implicitly links to SCALAPACK...
  # ... or use what the users specified in SCALAPACK_LIBRARY
  _scalapack_check_symbols("SCALAPACK")

  # netlib-scalapack
  if (NOT SCALAPACK_LIBRARY)
    find_library(netlib_LIBRARY NAMES scalapack)
    _scalapack_check_symbols("netlib")
  endif()
endif()

if (SCALAPACK_NOT_FOUND_MESSAGE)
  set(SCALAPACK_NOT_FOUND_MESSAGE REASON_FAILURE_MESSAGE "${SCALAPACK_NOT_FOUND_MESSAGE}")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SCALAPACK
  REQUIRED_VARS SCALAPACK_LIBRARY
  ${SCALAPACK_NOT_FOUND_MESSAGE})

# In case compiler implicitly links, remove empty variable placeholder
if (SCALAPACK_LIBRARY MATCHES "SCALAPACK-PLACEHOLDER-FOR-EMPTY-LIBRARIES")
  set(SCALAPACK_LIBRARY)
endif()

# Create the CMake target
if (SCALAPACK_FOUND)
  if (NOT TARGET SCALAPACK::SCALAPACK)
    add_library(SCALAPACK::SCALAPACK INTERFACE IMPORTED)
  endif()

  # note: if the compiler implicitly link to SCALAPACK
  if (SCALAPACK_LIBRARY)
    set_target_properties(SCALAPACK::SCALAPACK PROPERTIES
      IMPORTED_LOCATION ${SCALAPACK_LIBRARY})
  endif()
endif()
