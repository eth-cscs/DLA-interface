#
# Distributed Linear Algebra with Future (DLAF)
#
# Copyright (c) 2018-2021, ETH Zurich
# All rights reserved.
#
# Please, refer to the LICENSE file in the root directory.
# SPDX-License-Identifier: BSD-3-Clause
#

add_library(test_dlai_main test_dlai_main.cpp gtest_mpi_listener.cpp)
target_include_directories(test_dlai_main PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>)
target_link_libraries(test_dlai_main PUBLIC gtest dla_interface)

if (DLA_WITH_FORTRAN)
  add_library(test_ftn_dlai_main test_ftn_dlai_main.cpp test_ftn_setup.f90)
  #add_dependencies(test_ftn_dlai_main dla_interface)
  target_include_directories(test_ftn_dlai_main PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>)
  target_link_libraries(test_ftn_dlai_main PUBLIC gtest dla_interface)
endif()

add_library(test_tools ref_scalapack_tools.cpp tile_matrix_tools.cpp)
target_include_directories(test_tools PUBLIC ../include)
target_link_libraries(test_tools PUBLIC dla_interface)
