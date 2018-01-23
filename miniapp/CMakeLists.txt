include(dla_miniapp)

add_miniapp(overlap_gaussian_orbitals LIBS ${DLA_LIBS})

if(DLA_HAVE_SCALAPACK)
add_miniapp(overlap_gaussian_orbitals_c LIBS ${DLA_LIBS})
endif()
