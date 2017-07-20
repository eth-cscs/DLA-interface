include(dla_utils)
include(CheckCXXSymbolExists)

function(find_mpi)
  # Check if a MPI compiler is used, otherwise it looks for an MPI library.
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
endfunction()
