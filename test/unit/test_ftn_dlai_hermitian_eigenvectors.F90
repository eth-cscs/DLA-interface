module test_ftn_hermitian_eigenvectors
  use dla_interface
  use, intrinsic :: ISO_C_BINDING
  implicit none

  contains

#define TEST_FTN_HERMITIAN_EIGENVECTORS(test_function_name, function_name, CType, RealCType) \
    subroutine test_function_name(uplo, n, a, ia, ja, desca, evals, v, iv, jv, descv,    \
                                  solver, info) bind(C, name="test_function_name")      ;\
      use, intrinsic :: ISO_C_BINDING                                                   ;\
      implicit none                                                                     ;\
      character(kind=C_CHAR), intent(in) :: uplo                                        ;\
      integer(C_INT), intent(in) :: n, ia, ja, iv, jv                                   ;\
      integer(C_INT), dimension(*), intent(in) :: desca, descv                          ;\
      CType, dimension(desca(9),*), intent(inout) :: a                                  ;\
      RealCType, dimension(*), intent(out) :: evals                                     ;\
      CType, dimension(descv(9),*), intent(out) :: v                                    ;\
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
      call function_name(uplo, n, a, ia, ja, desca, evals, v, iv, jv, descv,             \
                         f_solver(1:lsolver), info)                                     ;\
    end subroutine

    TEST_FTN_HERMITIAN_EIGENVECTORS(test_ftn_s_hermitian_eigenvectors, dlai_s_hermitian_eigenvectors, real(C_FLOAT), real(C_FLOAT))
    TEST_FTN_HERMITIAN_EIGENVECTORS(test_ftn_d_hermitian_eigenvectors, dlai_d_hermitian_eigenvectors, real(C_DOUBLE), real(C_DOUBLE))
    TEST_FTN_HERMITIAN_EIGENVECTORS(test_ftn_c_hermitian_eigenvectors, dlai_c_hermitian_eigenvectors, complex(C_FLOAT_COMPLEX), real(C_FLOAT))
    TEST_FTN_HERMITIAN_EIGENVECTORS(test_ftn_z_hermitian_eigenvectors, dlai_z_hermitian_eigenvectors, complex(C_DOUBLE_COMPLEX), real(C_DOUBLE))

end module
