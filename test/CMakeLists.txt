#
# Distributed Linear Algebra Interface (DLAI)
#
# Copyright (c) 2018-2021, ETH Zurich
# All rights reserved.
#
# Please, refer to the LICENSE file in the root directory.
# SPDX-License-Identifier: BSD-3-Clause
#

add_library(DLAI_test INTERFACE)
target_include_directories(DLAI_test INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(DLAI_test INTERFACE gtest)

add_subdirectory(src)

add_subdirectory(unit)
