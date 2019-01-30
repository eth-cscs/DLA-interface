module test_ftn_lu_factorization
  use dla_interface
  use, intrinsic :: ISO_C_BINDING
  implicit none

  contains

#define TEST_FTN_LU_FACTORIZATION(test_function_name, function_name, CType) \
    subroutine test_function_name(m, n, a, ia, ja, desca, ipiv, solver, info)            \
        bind(C, name="test_function_name")                                              ;\
      use, intrinsic :: ISO_C_BINDING                                                   ;\
      implicit none                                                                     ;\
      integer(C_INT), intent(in) :: m, n, ia, ja                                        ;\
      integer(C_INT), dimension(*), intent(in) :: desca                                 ;\
      CType, dimension(desca(9),*), intent(inout) :: a                                  ;\
      integer(C_INT), dimension(*), intent(out) :: ipiv                                 ;\
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
      call function_name(m, n, a, ia, ja, desca, ipiv, f_solver(1:lsolver), info)       ;\
    end subroutine

    TEST_FTN_LU_FACTORIZATION(test_ftn_s_lu_factorization, dlai_s_lu_factorization, real(C_FLOAT))
    TEST_FTN_LU_FACTORIZATION(test_ftn_d_lu_factorization, dlai_d_lu_factorization, real(C_DOUBLE))
    TEST_FTN_LU_FACTORIZATION(test_ftn_c_lu_factorization, dlai_c_lu_factorization, complex(C_FLOAT_COMPLEX))
    TEST_FTN_LU_FACTORIZATION(test_ftn_z_lu_factorization, dlai_z_lu_factorization, complex(C_DOUBLE_COMPLEX))

end module
