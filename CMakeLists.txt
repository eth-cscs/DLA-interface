cmake_minimum_required(VERSION 3.8)
project(DLA_interface CXX C)

set(CMAKE_CXX_STANDARD 14)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)
# ---------------------------------------------------------------------------
# CMake module path
# ---------------------------------------------------------------------------
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

include(dla_utils)
include(dla_lapack)
include(dla_scalapack)
include(dla_dplasma)
include(dla_hwloc)
# ---------------------------------------------------------------------------
# Boost
# ---------------------------------------------------------------------------
find_package(Boost COMPONENTS program_options REQUIRED)

include(dla_openmp)
include(dla_mpi)

dla_find_lapack()

dla_find_scalapack()

dla_find_dplasma()

dla_find_hwloc()

setoption(DLA_PRINT_DEBUG_INFO BOOL "OFF" "Enables extra output for debugging purpose.")
if(DLA_PRINT_DEBUG_INFO)
  add_definitions("-DDLA_THREAD_DEBUG_INFO")
endif()

include(CTest)
enable_testing()
setoption(DLA_ALL_TESTS_USE_RUNNER BOOL "OFF" "Specify if the single rank tests has to use the runner as well")

include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${MPI_CXX_INCLUDE_PATH}
    ${DLA_DPLASMA_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
    ${DLA_LAPACK_INCLUDE_PATH}
    ${HWLOC_INCLUDE_DIRS}
)

# ---------------------------------------------------------------------------
# extra options
# ---------------------------------------------------------------------------
add_compile_options(-Wall -Wno-missing-braces)

# ---------------------------------------------------------------------------
# setup libraries
# ---------------------------------------------------------------------------
add_subdirectory(src)
set(DLA_LIBS
    dla_interface
    ${Boost_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
    ${DLA_DPLASMA_LIBS}
    ${DLA_SCALAPACK_LIBRARY}
    ${DLA_LAPACK_LIBRARY}
    ${DLA_HWLOC_LIBS}
)

message(${DLA_LIBS})

# ---------------------------------------------------------------------------
# Add sub directories
# ---------------------------------------------------------------------------
add_subdirectory(miniapp)
add_subdirectory(test)

