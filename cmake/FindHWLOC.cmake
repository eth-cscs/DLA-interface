#
# Distributed Linear Algebra Interface (DLAI)
#
# Copyright (c) 2018-2021, ETH Zurich
# All rights reserved.
#
# Please, refer to the LICENSE file in the root directory.
# SPDX-License-Identifier: BSD-3-Clause
#

# Find HWLOC library
#
# If HWLOC_ROOT is set, either in the environment or as CMake cache variable, it will
# be used for searching the library. Otherwise, pkg-config will be used.
#
# If found, this module creates the CMake target HWLOC::HWLOC

if (DEFINED HWLOC_ROOT OR DEFINED ENV{HWLOC_ROOT})
  find_path(HWLOC_INCLUDE_DIRS hwloc.h REQUIRED)
  find_library(HWLOC_LIBRARIES hwloc REQUIRED)
else()
  find_package(PkgConfig REQUIRED QUIET)
  pkg_search_module(HWLOC hwloc)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HWLOC DEFAULT_MSG
  HWLOC_INCLUDE_DIRS
  HWLOC_LIBRARIES
)

if (HWLOC_FOUND)
  if (NOT TARGET HWLOC::HWLOC)
    add_library(HWLOC::HWLOC INTERFACE IMPORTED GLOBAL)
  endif()

  target_include_directories(HWLOC::HWLOC INTERFACE ${HWLOC_INCLUDE_DIRS})
  target_link_libraries(HWLOC::HWLOC INTERFACE ${HWLOC_LIBRARIES})
endif()
