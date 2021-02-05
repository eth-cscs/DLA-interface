include(dla_miniapp)

add_miniapp(overlap_gaussian_orbitals LIBS ${DLA_LIBS})
add_miniapp(miniapp_matrix_multiplication EXE_NAME matrix_multiplication LIBS ${DLA_LIBS})
add_miniapp(small_test_communicator_grid LIBS ${DLA_LIBS})

if(DLA_HAVE_SCALAPACK)
add_miniapp(overlap_gaussian_orbitals_c LIBS ${DLA_LIBS})
endif()
