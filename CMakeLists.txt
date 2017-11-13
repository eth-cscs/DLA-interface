cmake_minimum_required(VERSION 3.0)
project(cholesky CXX)

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
include(dla_mpi)
include(dla_lapack)
include(dla_scalapack)
# ---------------------------------------------------------------------------
# Boost
# ---------------------------------------------------------------------------
find_package(Boost COMPONENTS program_options REQUIRED)

dla_find_mpi()

dla_find_lapack()
dla_find_scalapack()

include(CTest)
enable_testing()

include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${MPI_CXX_INCLUDE_PATH}
)

# ---------------------------------------------------------------------------
# extra options
# ---------------------------------------------------------------------------
add_compile_options(-Wall -Wno-missing-braces)
if(DLA_HAVE_SCALAPACK)
  add_definitions(-DDLA_HAVE_SCALAPACK)
endif()

# ---------------------------------------------------------------------------
# setup libraries
# ---------------------------------------------------------------------------
add_subdirectory(src)
set(DLA_LA_LIBS
    ${DLA_SCALAPACK_LIBRARY}
    ${DLA_LAPACK_LIBRARY}
)
set(DLA_LIBS
    dla_interface
    ${Boost_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
    ${DLA_LA_LIBS}
)
message(${DLA_LIBS})

# ---------------------------------------------------------------------------
# Add sub directories
# ---------------------------------------------------------------------------
add_subdirectory(test)

