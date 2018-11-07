module test_ftn_matrix_multiplication
  use dla_interface
  use, intrinsic :: ISO_C_BINDING
  implicit none

  contains

#define TEST_FTN_MATRIX_MULTIPLICATION(test_function_name, function_name, CType) \
    subroutine test_function_name(trans_a, trans_b, m, n, k, alpha, a, ia, ja, desca,    \
        b, ib, jb, descb, beta, c, ic, jc, descc, solver, info)                          \
        bind(C, name="test_function_name")                                              ;\
      use, intrinsic :: ISO_C_BINDING                                                   ;\
      implicit none                                                                     ;\
      character(kind=C_CHAR), intent(in) :: trans_a, trans_b                            ;\
      integer(C_INT), intent(in) :: m, n, k, ia, ja, ib, jb, ic, jc                     ;\
      CType, intent(in) :: alpha, beta                                                  ;\
      integer(C_INT), dimension(*), intent(in) :: desca, descb, descc                   ;\
      CType, dimension(desca(9),*), intent(in) :: a                                     ;\
      CType, dimension(descb(9),*), intent(in) :: b                                     ;\
      CType, dimension(descc(9),*), intent(inout) :: c                                  ;\
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
      call function_name(trans_a, trans_b, m, n, k, alpha, a, ia, ja, desca, b, ib, jb,  \
                         descb, beta, c, ic, jc, descc, f_solver(1:lsolver), info)      ;\
    end subroutine

    TEST_FTN_MATRIX_MULTIPLICATION(test_ftn_s_matrix_multiplication, dlai_s_matrix_multiplication, real(C_FLOAT))
    TEST_FTN_MATRIX_MULTIPLICATION(test_ftn_d_matrix_multiplication, dlai_d_matrix_multiplication, real(C_DOUBLE))
    TEST_FTN_MATRIX_MULTIPLICATION(test_ftn_c_matrix_multiplication, dlai_c_matrix_multiplication, complex(C_FLOAT_COMPLEX))
    TEST_FTN_MATRIX_MULTIPLICATION(test_ftn_z_matrix_multiplication, dlai_z_matrix_multiplication, complex(C_DOUBLE_COMPLEX))

end module
