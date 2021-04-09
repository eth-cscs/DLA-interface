DLAI_addTest(test_util_math
  SOURCES test_util_math.cpp
  LIBRARIES DLAI_CXX
  USE_MAIN PLAIN)

DLAI_addTest(test_util_memory
  SOURCES test_util_memory.cpp
  LIBRARIES DLAI_CXX
  USE_MAIN PLAIN)

DLAI_addTest(test_memory_allocator
  SOURCES test_memory_allocator.cpp
  LIBRARIES DLAI_CXX
  USE_MAIN PLAIN)

DLAI_addTest(test_local_matrix
  SOURCES test_local_matrix.cpp
  LIBRARIES DLAI_CXX
  USE_MAIN PLAIN)

DLAI_addTest(test_types
  SOURCES test_types.cpp
  LIBRARIES DLAI_CXX
  USE_MAIN PLAIN)

DLAI_addTest(test_util_types
  SOURCES test_util_types.cpp
  LIBRARIES DLAI_CXX
  USE_MAIN PLAIN)

DLAI_addTest(test_util_check
  SOURCES test_util_check.cpp
  LIBRARIES DLAI_CXX
  USE_MAIN PLAIN)

DLAI_addTest(test_thread_blas
  SOURCES test_thread_blas.cpp
  LIBRARIES DLAI_CXX
  USE_MAIN PLAIN)

DLAI_addTest(test_thread_binding
  SOURCES test_thread_binding.cpp
  LIBRARIES DLAI_CXX
  USE_MAIN PLAIN)

DLAI_addTest(test_util_thread
  SOURCES test_util_thread.cpp
  LIBRARIES DLAI_CXX
  USE_MAIN MPI
  MPIRANKS 1)

DLAI_addTest(test_communicator_grid
  SOURCES test_communicator_grid.cpp
  LIBRARIES DLAI_CXX gtest_mpi_listener
  USE_MAIN MPI
  MPIRANKS 6)

DLAI_addTest(test_communicator_manager
  SOURCES test_communicator_manager.cpp
  LIBRARIES DLAI_CXX gtest_mpi_listener
  USE_MAIN MPI
  MPIRANKS 6)

DLAI_addTest(test_communicator_manager_init
  SOURCES test_communicator_manager.cpp
  LIBRARIES DLAI_CXX
  COMPILE_DEFINITIONS COMM_INITS_MPI
  USE_MAIN MPI
  MPIRANKS 6)

DLAI_addTest(test_distributed_matrix
  SOURCES test_distributed_matrix.cpp
  LIBRARIES DLAI_CXX test_tools
  USE_MAIN MPI
  MPIRANKS 6)

DLAI_addTest(test_dlai_cholesky_factorization
  SOURCES test_dlai_cholesky_factorization.cpp
  LIBRARIES DLAI_CXX test_tools
  USE_MAIN MPI
  MPIRANKS 6)

DLAI_addTest(test_dlai_cholesky_inverse
  SOURCES test_dlai_cholesky_inverse.cpp
  LIBRARIES DLAI_CXX test_tools
  USE_MAIN MPI
  MPIRANKS 6)

DLAI_addTest(test_dlai_matrix_multiplication
  SOURCES test_dlai_matrix_multiplication.cpp
  LIBRARIES DLAI_CXX test_tools
  USE_MAIN MPI
  MPIRANKS 6)

DLAI_addTest(test_dlai_lu_factorization
  SOURCES test_dlai_lu_factorization.cpp
  LIBRARIES DLAI_CXX test_tools
  USE_MAIN MPI
  MPIRANKS 6)

DLAI_addTest(test_dlai_hermitian_eigenvectors
  SOURCES test_dlai_hermitian_eigenvectors.cpp
  LIBRARIES DLAI_CXX test_tools
  USE_MAIN MPI
  MPIRANKS 6)

DLAI_addTest(test_dlai_triangular_inverse
  SOURCES test_dlai_triangular_inverse.cpp
  LIBRARIES DLAI_CXX test_tools
  USE_MAIN MPI
  MPIRANKS 6)

if (DLAI_WITH_FORTRAN)
  DLAI_addTest(test_ftn_dlai_cholesky_factorization
    SOURCES test_ftn_dlai_cholesky_factorization.cpp test_ftn_dlai_cholesky_factorization.F90
    LIBRARIES DLAI_Fortran test_tools
    USE_MAIN FTN_MPI
    MPIRANKS 6)

  DLAI_addTest(
    test_ftn_dlai_cholesky_inverse
    SOURCES test_ftn_dlai_cholesky_inverse.cpp test_ftn_dlai_cholesky_inverse.F90
    LIBRARIES DLAI_Fortran test_tools
    USE_MAIN FTN_MPI
    MPIRANKS 6)

  DLAI_addTest(test_ftn_dlai_matrix_multiplication
    SOURCES test_ftn_dlai_matrix_multiplication.cpp test_ftn_dlai_matrix_multiplication.F90
    LIBRARIES DLAI_Fortran test_tools
    USE_MAIN FTN_MPI
    MPIRANKS 6)

  DLAI_addTest(test_ftn_dlai_lu_factorization
    SOURCES test_ftn_dlai_lu_factorization.cpp test_ftn_dlai_lu_factorization.F90
    LIBRARIES DLAI_Fortran test_tools
    USE_MAIN FTN_MPI
    MPIRANKS 6)

  DLAI_addTest(test_ftn_dlai_hermitian_eigenvectors
    SOURCES test_ftn_dlai_hermitian_eigenvectors.cpp test_ftn_dlai_hermitian_eigenvectors.F90
    LIBRARIES DLAI_Fortran test_tools
    USE_MAIN FTN_MPI
    MPIRANKS 6)

  DLAI_addTest(test_ftn_dlai_triangular_inverse
    SOURCES test_ftn_dlai_triangular_inverse.cpp test_ftn_dlai_triangular_inverse.F90
    LIBRARIES DLAI_Fortran test_tools
    USE_MAIN FTN_MPI
    MPIRANKS 6)
endif()
