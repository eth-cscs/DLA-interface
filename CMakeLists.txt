cmake_minimum_required(VERSION 3.0)
project(cholesky CXX)

set(CMAKE_CXX_STANDARD 14)

# ---------------------------------------------------------------------------
# CMake module path
# ---------------------------------------------------------------------------
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

include(dla_utils)

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
if(COMPILER_HAS_MPI)
  setoption(TEST_RUNNER FILEPATH "mpirun" "Executable for running MPI tests")
  setoption(TEST_RUNNER_NP_OPT STRING "-n" "Option passed to TEST_RUNNER to specify the number of ranks")
else()
  find_package(MPI REQUIRED)
  setoption(TEST_RUNNER FILEPATH "${MPIEXEC}" "Executable for running MPI tests")
  setoption(TEST_RUNNER_NP_OPT STRING "${MPIEXEC_NUMPROC_FLAG}" "Option passed to TEST_RUNNER to specify the number of ranks")
endif()

include(CTest)
enable_testing()

include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${MPI_CXX_INCLUDE_PATH}
)

# ---------------------------------------------------------------------------
# setup libraries
# ---------------------------------------------------------------------------
set(DLA_LIBS
    ${Boost_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
)

# ---------------------------------------------------------------------------
# extra options
# ---------------------------------------------------------------------------
add_compile_options(-Wall -Wno-missing-braces)

# ---------------------------------------------------------------------------
# Add sub directories
# ---------------------------------------------------------------------------
add_subdirectory(test)

