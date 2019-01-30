add_library(gtest gtest/src/gtest-all.cpp)
target_include_directories(gtest PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/gtest>)
target_link_libraries(gtest ${MPI_CXX_LIBRARIES})

add_library(gtest_main gtest/src/gtest_main.cpp)
target_include_directories(gtest_main PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/gtest>)

add_library(test_dlai_main src/test_dlai_main.cpp)
target_include_directories(test_dlai_main PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/gtest>)
target_include_directories(test_dlai_main PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

if (${DLA_WITH_FORTRAN})
  add_library(test_ftn_dlai_main src/test_ftn_dlai_main.cpp src/test_ftn_setup.f90)
  add_dependencies(test_ftn_dlai_main dla_interface)
  target_include_directories(test_ftn_dlai_main PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/gtest>)
  target_include_directories(test_ftn_dlai_main PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
endif()

add_library(test_tools src/ref_scalapack_tools.cpp src/tile_matrix_tools.cpp)
target_include_directories(test_tools PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

include(dla_test)

add_subdirectory(unit)
add_subdirectory(extra)
