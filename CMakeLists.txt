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

include(dla_openmp)
include(dla_mpi)

dla_find_lapack()

dla_find_scalapack()

dla_find_dplasma()

dla_find_hwloc()

if (${DLA_HAVE_SCALAPACK})
  setoption(DLA_WITH_FORTRAN BOOL ON "Build fortran binding")
else()
  setoption(DLA_WITH_FORTRAN BOOL OFF "Build fortran binding")
endif()

if (${DLA_WITH_FORTRAN})
  if (${DLA_HAVE_SCALAPACK})
    enable_language(Fortran)
  else()
    message(WARNING "Cannot build Fortran bindings without Scalapack. Setting DLA_WITH_FORTRAN to OFF")
    set(${DLA_WITH_FORTRAN} OFF FORCE)
  endif()
endif()

setoption(DLA_PRINT_DEBUG_INFO BOOL "OFF" "Enables extra output for debugging purpose.")
if(DLA_PRINT_DEBUG_INFO)
  add_definitions("-DDLA_THREAD_DEBUG_INFO")
endif()
setoption(DLA_PRINT_DEBUG_CALL_PARAM BOOL "OFF" "Enables call parameters output for debugging purpose.")
if(DLA_PRINT_DEBUG_CALL_PARAM)
  add_definitions("-DDLA_PRINT_DEBUG_CALL_PARAM")
endif()
setoption(DLA_COVERAGE_TEST BOOL "OFF" "Enables coverage test.")
if(DLA_COVERAGE_TEST)
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
    set(CMAKE_LD_FLAGS "${CMAKE_LD_FLAGS} --coverage")

    set(DLA_COVERAGE_SRC_DIRS
        ${CMAKE_SOURCE_DIR}/include/*
        ${CMAKE_SOURCE_DIR}/include/*/*
        ${CMAKE_SOURCE_DIR}/src/*
        ${CMAKE_SOURCE_DIR}/test/include/*
        ${CMAKE_SOURCE_DIR}/test/src/*
        ${CMAKE_SOURCE_DIR}/test/unit/*
        )

    add_custom_target("coverage"
      COMMAND lcov --capture --directory ${CMAKE_BINARY_DIR} --output-file coverage.info
      COMMAND lcov --extract coverage.info ${DLA_COVERAGE_SRC_DIRS} --output-file coverage-filt.info
      COMMAND genhtml coverage-filt.info --output-directory ${CMAKE_BINARY_DIR}/coverage
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

  else ()
    message(WARNING "Coverage test can be used only with GCC.")
    set(DLA_COVERAGE_TEST "OFF" FORCE)
  endif()
endif()

include(CTest)
enable_testing()
setoption(DLA_ALL_TESTS_USE_RUNNER BOOL "OFF" "Specify if the single rank tests has to use the runner as well")

include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/libs
    ${MPI_CXX_INCLUDE_PATH}
    ${DLA_DPLASMA_INCLUDE_DIRS}
    ${DLA_LAPACK_INCLUDE_PATH}
    ${HWLOC_INCLUDE_DIRS}
)

# ---------------------------------------------------------------------------
# extra options
# ---------------------------------------------------------------------------
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-missing-braces")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -ffree-line-length-0")

# ---------------------------------------------------------------------------
# setup libraries
# ---------------------------------------------------------------------------
add_subdirectory(src)
set(DLA_LIBS
    dla_interface
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

