#
# Distributed Linear Algebra Interface (DLAI)
#
# Copyright (c) 2018-2021, ETH Zurich
# All rights reserved.
#
# Please, refer to the LICENSE file in the root directory.
# SPDX-License-Identifier: BSD-3-Clause
#

# MPI_PRESET_SETUP
#
# Helper macro enabling preset selection for MPI
macro(mpi_preset_setup)
  set(MPIEXEC_PRESET "custom" CACHE STRING "Select a preset to use")
  set_property(CACHE MPIEXEC_PRESET PROPERTY STRINGS "custom" "slurm")

  if (NOT MPIEXEC_PRESET OR MPIEXEC_PRESET STREQUAL "custom")
    # Classic, it uses default from MPI
  elseif (MPIEXEC_PRESET STREQUAL "slurm")
    execute_process(COMMAND which srun
      OUTPUT_VARIABLE SLURM_EXECUTABLE
      OUTPUT_STRIP_TRAILING_WHITESPACE)

    set(MPIEXEC_EXECUTABLE ${SLURM_EXECUTABLE} CACHE STRING "" FORCE)
    set(MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING "" FORCE)
    set(MPIEXEC_NUMCORE_FLAG "-c" CACHE STRING "" FORCE)
    set(MPIEXEC_NUMCORES "1" CACHE STRING "")

    message(STATUS "MPI preset: ${MPIEXEC_PRESET}")
  else()
    message(FATAL_ERROR "Preset ${MPIEXEC_PRESET} is not supported")
  endif()
endmacro()

# MPIEXEC_RUNLINE
#
# Returns in OUT_COMMAND the correct command line to run it with MPIEXEC_EXECUTABLE
# In particular, it will run EXECUTABLE_NAME on MPIRANKS, and in case MPIEXEC_NUMCORES
# is set (e.g. by the slurm MPI preset), it also computes the correct number of cores
# to assign to each rank (assuming single node allocation)
function(mpiexec_runline OUT_COMMAND EXECUTABLE_NAME MPIRANKS)
  if (MPIEXEC_NUMCORE_FLAG)
    if (MPIEXEC_NUMCORES)
      set(_NCORES ${MPIEXEC_NUMCORES})
    else()
      set(_NCORES 1)
    endif()

    math(EXPR CORE_PER_RANK "${_NCORES}/${MPIRANKS}")

    if (NOT CORE_PER_RANK)
      set(CORE_PER_RANK 1)
    endif()

    set(_MPI_CORE_ARGS ${MPIEXEC_NUMCORE_FLAG} ${CORE_PER_RANK})
  else()
    set(_MPI_CORE_ARGS "")
  endif()

  set(${OUT_COMMAND}
    ${MPIEXEC_EXECUTABLE}
    ${MPIEXEC_NUMPROC_FLAG} ${MPIRANKS}
    ${_MPI_CORE_ARGS}
    ${MPIEXEC_PREFLAGS} ${EXECUTABLE_NAME} ${MPIEXEC_POSTFLAGS}
    PARENT_SCOPE)
endfunction()

# DLAI_addTest(test_target_name
#   SOURCES <source1> [<source2> ...]
#   [COMPILE_DEFINITIONS <arguments for target_compile_definitions>]
#   [INCLUDE_DIRS <arguments for target_include_directories>]
#   [LIBRARIES <arguments for target_link_libraries>]
#   [MPIRANKS <number of rank>]
#   [USE_MAIN {PLAIN | MPI | FTN_MPI}]
# )
#
# At least one source file has to be specified, while other parameters are optional.
#
# COMPILE_DEFINITIONS, INCLUDE_DIRS and LIBRARIES are passed to respective cmake wrappers, so it is
# possible to specify PRIVATE/INTERFACE/PUBLIC modifiers.
#
# MPIRANKS specifies the number of ranks on which the test will be carried out and it implies a link with
# MPI library. At build time the constant NUM_MPI_RANKS=MPIRANKS is set.
#
# USE_MAIN links to an external main function, in particular:
#   - PLAIN: uses the classic gtest_main
#   - MPI: uses a main that initializes MPI
#   - FTN_MPI: uses a main that initializes MPI with Fortran
# If not specified, no external main is used and it should exist in the test source code.
#
# e.g.
#
# DLAI_addTest(example_test
#   SOURCE main.cpp testfixture.cpp
#   LIBRARIES
#     PRIVATE
#       boost::boost
# )
function(DLAI_addTest test_target_name)
  set(options "")
  set(oneValueArgs MPIRANKS USE_MAIN)
  set(multiValueArgs SOURCES COMPILE_DEFINITIONS INCLUDE_DIRS LIBRARIES ARGUMENTS)
  cmake_parse_arguments(DLAI_AT "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  option(DLAI_TEST_RUNALL_WITH_MPIEXEC "Run all tests using the workload manager." OFF)
  option(DLAI_CI_RUNNER_USES_MPIRUN "Remove mpiexec command for tests executed by ctest. This option is to be used if the CI runner executes the tests with <mpiexec + options> ctest -L RANK_<spawned MPI ranks>" OFF)
  option(DLAI_INSTALL_TESTS "If tests are built, it controls if they will be installed" OFF)

  ### Checks
  if (DLAI_AT_UNPARSED_ARGUMENTS)
    message(FATAL_ERROR "Unknown arguments ${DLAI_AT_UNPARSED_ARGUMENTS}")
  endif()

  if (NOT DLAI_AT_SOURCES)
    message(FATAL_ERROR "No sources specified for this test")
  endif()

  set(IS_AN_MPI_TEST FALSE)
  if (NOT DLAI_AT_USE_MAIN)
    set(_gtest_tgt gtest)
  elseif (DLAI_AT_USE_MAIN STREQUAL PLAIN)
    set(_gtest_tgt gtest_main)
  elseif (DLAI_AT_USE_MAIN STREQUAL MPI)
    set(_gtest_tgt test_dlai_main)
    set(IS_AN_MPI_TEST TRUE)
  elseif (DLAI_AT_USE_MAIN STREQUAL FTN_MPI)
    set(_gtest_tgt test_ftn_dlai_main)
    set(IS_AN_MPI_TEST TRUE)
  else()
    message(FATAL_ERROR "USE_MAIN=${DLAI_AT_USE_MAIN} is not a supported option")
  endif()

  if (IS_AN_MPI_TEST)
    if (NOT DLAI_AT_MPIRANKS)
      message(FATAL_ERROR "You are asking for an MPI external main without specifying MPIRANKS")
    endif()
    if (NOT DLAI_AT_MPIRANKS GREATER 0)
      message(FATAL_ERROR "Wrong MPIRANKS number ${DLAI_AT_MPIRANKS}")
    endif()
    if (DLAI_AT_MPIRANKS GREATER MPIEXEC_MAX_NUMPROCS)
      message(WARNING "\
      YOU ARE ASKING FOR ${DLAI_AT_MPIRANKS} RANKS, BUT THERE ARE JUST ${MPIEXEC_MAX_NUMPROCS} CORES.
      You can adjust MPIEXEC_MAX_NUMPROCS value to suppress this warning.
      Using OpenMPI may require to set the environment variable OMPI_MCA_rmaps_base_oversubscribe=1.")
    endif()
  else()
    if (DLAI_AT_MPIRANKS)
      message(FATAL_ERROR "You specified MPIRANKS and asked for an external main without MPI")
    else()
      set(DLAI_AT_MPIRANKS 1)
    endif()
  endif()

  ### Test target
  set(_TEST_ARGUMENTS ${DLAI_AT_ARGUMENTS})

  if (DLAI_TEST_RUNALL_WITH_MPIEXEC OR IS_AN_MPI_TEST)
    if(DLAI_CI_RUNNER_USES_MPIRUN)
      set(_TEST_COMMAND $<TARGET_FILE:${test_target_name}>)
    else()
      mpiexec_runline(_TEST_COMMAND $<TARGET_FILE:${test_target_name}> ${DLAI_AT_MPIRANKS})
    endif()
    set(_TEST_LABEL "RANK_${DLAI_AT_MPIRANKS}")

  # ----- Classic test
  else()
    set(_TEST_COMMAND $<TARGET_FILE:${test_target_name}>)
    set(_TEST_LABEL "RANK_1")
  endif()

  ### Test executable target
  add_executable(${test_target_name} ${DLAI_AT_SOURCES})
  target_link_libraries(${test_target_name}
    PRIVATE
      ${_gtest_tgt}
      DLAI_test
      ${DLAI_AT_LIBRARIES}
  )
  target_compile_definitions(${test_target_name}
    PRIVATE
      ${DLAI_AT_COMPILE_DEFINITIONS}
      $<$<BOOL:${IS_AN_MPI_TEST}>: NUM_MPI_RANKS=${DLAI_AT_MPIRANKS}>)
  target_include_directories(${test_target_name} PRIVATE ${DLAI_AT_INCLUDE_DIRS})
  target_add_warnings(${test_target_name})
  add_test(
    NAME ${test_target_name}
    COMMAND ${_TEST_COMMAND} ${_TEST_ARGUMENTS})
  set_tests_properties(${test_target_name} PROPERTIES LABELS "${_TEST_LABEL}")

  ### DEPLOY
  include(GNUInstallDirs)

  if (DLAI_INSTALL_TESTS)
    install(TARGETS
      ${test_target_name}
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  endif()
endfunction()
