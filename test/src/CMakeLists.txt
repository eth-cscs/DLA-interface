#
# Distributed Linear Algebra Interface (DLAI)
#
# Copyright (c) 2018-2021, ETH Zurich
# All rights reserved.
#
# Please, refer to the LICENSE file in the root directory.
# SPDX-License-Identifier: BSD-3-Clause
#

add_library(gtest_mpi_listener OBJECT gtest_mpi_listener.cpp)
target_include_directories(gtest_mpi_listener
  PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(gtest_mpi_listener
  PUBLIC gtest)

add_library(test_dlai_main test_dlai_main.cpp)
target_link_libraries(test_dlai_main
  PUBLIC
    gtest
    gtest_mpi_listener
    DLAI
    DLAI_test)

if (DLAI_WITH_FORTRAN)
  add_library(test_ftn_dlai_main
    test_ftn_dlai_main.cpp
    test_ftn_setup.f90)
  target_link_libraries(test_ftn_dlai_main
    PUBLIC
      gtest
      gtest_mpi_listener
      DLAI
      DLAI_test)
endif()

add_library(test_tools ref_scalapack_tools.cpp tile_matrix_tools.cpp)
target_include_directories(test_tools PUBLIC ../include)
target_link_libraries(test_tools PUBLIC DLAI)
