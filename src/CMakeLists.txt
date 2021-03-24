#
# Distributed Linear Algebra Interface (DLAI)
#
# Copyright (c) 2018-2021, ETH Zurich
# All rights reserved.
#
# Please, refer to the LICENSE file in the root directory.
# SPDX-License-Identifier: BSD-3-Clause
#

add_library(DLAI_CXX SHARED
  communicator_grid.cpp
  communicator_manager.cpp
  $<$<BOOL:${SCALAPACK_FOUND}>:c_dla_interface.cpp>
)

target_compile_definitions(DLAI_CXX PUBLIC
  $<$<BOOL:${SCALAPACK_FOUND}>:DLAI_WITH_SCALAPACK>
  $<$<BOOL:${HWLOC_FOUND}>:DLAI_WITH_HWLOC>

  $<$<BOOL:${DLAI_WITH_DLAF}>:DLAI_WITH_DLAF>
  $<$<BOOL:${DLAI_WITH_DPLASMA}>:DLAI_WITH_DPLASMA>
  $<$<BOOL:${DLAI_WITH_ELPA}>:DLAI_WITH_ELPA>

  $<$<BOOL:${DLAI_PRINT_DEBUG_INFO}>:DLAI_PRINT_DEBUG_INFO>
  $<$<BOOL:${DLAI_PRINT_DEBUG_CALL_PARAM}>:DLAI_PRINT_DEBUG_CALL_PARAM>
)

target_include_directories(DLAI_CXX PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>

  # TODO this folder contains just cxxopts includes, and it will be installed in "include" path
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/libs>
  $<INSTALL_INTERFACE:include>
)

target_link_libraries(DLAI_CXX
  PUBLIC
    # backends
    $<$<BOOL:${DLAI_WITH_DLAF}>:DLAF::DLAF>
    $<$<BOOL:${DLAI_WITH_DPLASMA}>:DPLASMA::dplasma>
    $<$<BOOL:${DLAI_WITH_ELPA}>:ELPA::ELPA>

    $<$<BOOL:${HWLOC_FOUND}>:HWLOC::HWLOC>

    LAPACK::LAPACK
    $<$<BOOL:${SCALAPACK_FOUND}>:SCALAPACK::SCALAPACK>

  PRIVATE
    MPI::MPI_CXX
    OpenMP::OpenMP_CXX
)

if (DLAI_WITH_FORTRAN)
  add_library(DLAI_Fortran SHARED
    f_dla_interface.F90
  )

  target_compile_options(DLAI_Fortran PUBLIC
    $<$<COMPILE_LANGUAGE:Fortran>:-ffree-line-length-0>
  )

  target_include_directories(DLAI_Fortran PUBLIC
    # *.mod files
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<INSTALL_INTERFACE:include>
  )

  target_link_libraries(DLAI_Fortran
    PUBLIC
      MPI::MPI_Fortran
    PRIVATE
      DLAI_CXX
      OpenMP::OpenMP_Fortran
  )
endif()

# ----- DEPLOY
include(GNUInstallDirs)

install(TARGETS
  DLAI_CXX DLAI_Fortran
  EXPORT DLAI-Targets
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# install includes
install(DIRECTORY ../include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(DIRECTORY ../libs/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# install custom FindModules
install(DIRECTORY ../cmake/
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${CMAKE_PROJECT_NAME}
  FILES_MATCHING PATTERN "Find*.cmake" # TODO problem with lapack.cmake
  PATTERN "template" EXCLUDE
)

# ----- CMake INTEGRATION
include(CMakePackageConfigHelpers)

# install targets configuration
install(EXPORT
  DLAI-Targets
  NAMESPACE DLAI::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${CMAKE_PROJECT_NAME}
)

# Config-file preparation and install
configure_package_config_file(
  ${CMAKE_CURRENT_LIST_DIR}/../cmake/template/DLAIConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/DLAIConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${CMAKE_PROJECT_NAME}
)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/DLAIConfig.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${CMAKE_PROJECT_NAME}
)
