module test_ftn_triangular_inverse
  use dla_interface
  use, intrinsic :: ISO_C_BINDING
  implicit none

  contains

#define TEST_FTN_TRIANGULAR_INVERSE(test_function_name, function_name, CType) \
    subroutine test_function_name(uplo, diag, n, a, ia, ja, desca, solver, info)         \
        bind(C, name="test_function_name")                                              ;\
      use, intrinsic :: ISO_C_BINDING                                                   ;\
      implicit none                                                                     ;\
      character(kind=C_CHAR), intent(in) :: uplo, diag                                  ;\
      integer(C_INT), intent(in) :: n, ia, ja                                           ;\
      integer(C_INT), dimension(*), intent(in) :: desca                                 ;\
      CType, dimension(desca(9),*), intent(inout) :: a                                  ;\
      character(len=1, kind=C_CHAR), dimension(*), intent(in) :: solver(*)              ;\
      integer(C_INT), intent(out) :: info                                               ;\
                                                                                        ;\
      integer lsolver                                                                   ;\
      character(len=50) :: f_solver                                                     ;\
                                                                                        ;\
      info = 0                                                                          ;\
      lsolver = 0                                                                       ;\
      do while (solver(lsolver + 1) .ne. C_NULL_CHAR .and. lsolver .lt. 50)             ;\
        lsolver = lsolver + 1                                                           ;\
        f_solver(lsolver:lsolver) = solver(lsolver)                                     ;\
      end do                                                                            ;\
                                                                                        ;\
      call function_name(uplo, diag, n, a, ia, ja, desca, f_solver(1:lsolver), info)    ;\
    end subroutine

    TEST_FTN_TRIANGULAR_INVERSE(test_ftn_s_triangular_inverse, dlai_s_triangular_inverse, real(C_FLOAT))
    TEST_FTN_TRIANGULAR_INVERSE(test_ftn_d_triangular_inverse, dlai_d_triangular_inverse, real(C_DOUBLE))
    TEST_FTN_TRIANGULAR_INVERSE(test_ftn_c_triangular_inverse, dlai_c_triangular_inverse, complex(C_FLOAT_COMPLEX))
    TEST_FTN_TRIANGULAR_INVERSE(test_ftn_z_triangular_inverse, dlai_z_triangular_inverse, complex(C_DOUBLE_COMPLEX))

end module
