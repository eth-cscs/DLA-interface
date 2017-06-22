cmake_minimum_required(VERSION 3.0)
project(cholesky CXX)

set(CMAKE_CXX_STANDARD 14)

# ---------------------------------------------------------------------------
# CMake module path
# ---------------------------------------------------------------------------
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# ---------------------------------------------------------------------------
# Boost
# ---------------------------------------------------------------------------
find_package(Boost COMPONENTS program_options REQUIRED)

# ---------------------------------------------------------------------------
# MPI
# ---------------------------------------------------------------------------
# Check if a MPI compiler is used
include(CheckCXXSymbolExists)
CHECK_CXX_SYMBOL_EXISTS(MPI_Finalize mpi.h COMPILER_HAS_MPI)
# If the compiler is not an MPI compiler find MPI library
if(NOT COMPILER_HAS_MPI)
  find_package(MPI REQUIRED)
endif()

include(CTest)
enable_testing()

include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${MPI_C_INCLUDE_PATH}
)

# ---------------------------------------------------------------------------
# setup libraries
# ---------------------------------------------------------------------------
set(DLA_LIBS
    ${Boost_LIBRARIES}
    ${MPI_C_LIBRARIES}
)

# ---------------------------------------------------------------------------
# extra options
# ---------------------------------------------------------------------------
add_compile_options(-Wall -Wno-missing-braces)

# ---------------------------------------------------------------------------
# Add sub directories
# ---------------------------------------------------------------------------
add_subdirectory(test)

