cmake_host_system_information(RESULT NCPU QUERY NUMBER_OF_PHYSICAL_CORES)

add_unit_test(gtest_unittest LIBS gtest_main)
add_unit_test(test_util_math LIBS gtest_main)
add_unit_test(test_util_memory LIBS gtest_main)
add_unit_test(test_memory_allocator LIBS gtest_main)
add_unit_test(test_local_matrix LIBS gtest_main)
add_unit_test(test_types LIBS gtest_main)
add_unit_test(test_util_types LIBS gtest_main)
add_unit_test(test_util_check LIBS gtest_main)
add_unit_test(test_thread_blas LIBS gtest_main ${DLA_LIBS})
add_unit_test(test_thread_binding LIBS gtest_main ${DLA_LIBS})
add_unit_test(test_util_thread MPIRANKS 1 LIBS ${DLA_LIBS})


if(NCPU GREATER_EQUAL 6)
    add_unit_test(test_communicator_grid MPIRANKS 6 LIBS ${DLA_LIBS})
    add_unit_test(test_communicator_manager MPIRANKS 6 LIBS ${DLA_LIBS})
    add_unit_test(test_communicator_manager_init SRC test_communicator_manager.cpp MPIRANKS 6 DEFS -DCOMM_INITS_MPI LIBS ${DLA_LIBS})
    add_unit_test(test_distributed_matrix MPIRANKS 6 LIBS test_tools ${DLA_LIBS})
    add_unit_test(test_dlai_cholesky_factorization MPIRANKS 6 LIBS test_tools test_dlai_main ${DLA_LIBS})
    add_unit_test(test_dlai_cholesky_inverse MPIRANKS 6 LIBS test_tools test_dlai_main ${DLA_LIBS})
    add_unit_test(test_dlai_matrix_multiplication MPIRANKS 6 LIBS test_tools test_dlai_main ${DLA_LIBS})
    add_unit_test(test_dlai_lu_factorization MPIRANKS 6 LIBS test_tools test_dlai_main ${DLA_LIBS})
    add_unit_test(test_dlai_hermitian_eigenvectors MPIRANKS 6 LIBS test_tools test_dlai_main ${DLA_LIBS})
    add_unit_test(test_dlai_triangular_inverse MPIRANKS 6 LIBS test_tools test_dlai_main ${DLA_LIBS})
  if (${DLA_WITH_FORTRAN})
    add_unit_test(test_ftn_dlai_cholesky_factorization SRC test_ftn_dlai_cholesky_factorization.cpp test_ftn_dlai_cholesky_factorization.F90 MPIRANKS 6 LIBS test_tools test_ftn_dlai_main ${DLA_LIBS})
    add_unit_test(test_ftn_dlai_cholesky_inverse SRC test_ftn_dlai_cholesky_inverse.cpp test_ftn_dlai_cholesky_inverse.F90 MPIRANKS 6 LIBS test_tools test_ftn_dlai_main ${DLA_LIBS})
    add_unit_test(test_ftn_dlai_matrix_multiplication SRC test_ftn_dlai_matrix_multiplication.cpp test_ftn_dlai_matrix_multiplication.F90 MPIRANKS 6 LIBS test_tools test_ftn_dlai_main ${DLA_LIBS})
    add_unit_test(test_ftn_dlai_lu_factorization SRC test_ftn_dlai_lu_factorization.cpp test_ftn_dlai_lu_factorization.F90 MPIRANKS 6 LIBS test_tools test_ftn_dlai_main ${DLA_LIBS})
    add_unit_test(test_ftn_dlai_hermitian_eigenvectors SRC test_ftn_dlai_hermitian_eigenvectors.cpp test_ftn_dlai_hermitian_eigenvectors.F90 MPIRANKS 6 LIBS test_tools test_ftn_dlai_main ${DLA_LIBS})
    add_unit_test(test_ftn_dlai_triangular_inverse SRC test_ftn_dlai_triangular_inverse.cpp test_ftn_dlai_triangular_inverse.F90 MPIRANKS 6 LIBS test_tools test_ftn_dlai_main ${DLA_LIBS})
  endif()
else()
    message(WARNING "Number of detected physical cores is ${NCPU}.\nDisabled tests which require 6 MPI slots:\n \
\ttest_communicator_grid\n \
\ttest_communicator_manager\n \
\ttest_communicator_manager_init\n \
\ttest_distributed_matrix\n \
\ttest_dlai_cholesky_factorization\n \
\ttest_dlai_cholesky_inverse\n \
\ttest_dlai_matrix_multiplication\n \
\ttest_dlai_lu_factorization\n \
\ttest_dlai_hermitian_eigenvectors\n \
\ttest_dlai_triangular_inverse\n \
\ttest_ftn_dlai_cholesky_factorization\n \
\ttest_ftn_dlai_cholesky_inverse\n \
\ttest_ftn_dlai_matrix_multiplication\n \
\ttest_ftn_dlai_lu_factorization\n \
\ttest_ftn_dlai_hermitian_eigenvectors\n \
\ttest_ftn_dlai_triangular_inverse")
endif()

