add_executable(overlap_gaussian_orbitals overlap_gaussian_orbitals.cpp)
target_link_libraries(overlap_gaussian_orbitals PUBLIC DLAI)

add_executable(matrix_multiplication matrix_multiplication.cpp)
target_link_libraries(matrix_multiplication PUBLIC DLAI)

if (SCALAPACK_FOUND)
  add_executable(overlap_gaussian_orbitals_c overlap_gaussian_orbitals_c.cpp)
  target_link_libraries(overlap_gaussian_orbitals_c PUBLIC DLAI)
endif()
