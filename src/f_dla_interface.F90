module dla_interface

  use, intrinsic :: ISO_C_BINDING

  interface
    subroutine dlai_initialize(nr_cores, init_mpi) &
        &bind(C, name="dlai_initialize")
      use, intrinsic :: ISO_C_BINDING
      implicit none
      integer(C_INT), intent(in) :: nr_cores
      integer(C_INT), intent(in) :: init_mpi
    end subroutine

    subroutine dlai_finalize() &
        &bind(C, name="dlai_finalize")
      use, intrinsic :: ISO_C_BINDING
      implicit none
    end subroutine

    integer(C_INT) function dlai_create_2d_grid(base_comm, row_size, col_size, ordering) &
        &bind(C, name="dlai_create_2d_grid")
      use, intrinsic :: ISO_C_BINDING
      implicit none
      integer(C_INT), intent(in) :: base_comm
      integer(C_INT), intent(in) :: row_size
      integer(C_INT), intent(in) :: col_size
      character(len=1, kind=C_CHAR), intent(in) :: ordering
    end function

    integer(C_INT) function dlai_create_2d_grid_blacs(blacs_handle, row_size, col_size, ordering) &
        &bind(C, name="dlai_create_2d_grid_blacs")
      use, intrinsic :: ISO_C_BINDING
      implicit none
      integer(C_INT), intent(in) :: blacs_handle
      integer(C_INT), intent(in) :: row_size
      integer(C_INT), intent(in) :: col_size
      character(len=1, kind=C_CHAR), intent(in) :: ordering
    end function

    subroutine dlai_free_2d_grid_blacs(blacs_context) &
        &bind(C, name="dlai_free_2d_grid_blacs")
      use, intrinsic :: ISO_C_BINDING
      implicit none
      integer(C_INT), intent(in) :: blacs_context
    end subroutine

    integer(C_INT) function dlai_get_print_timer_option() &
        &bind(C, name="dlai_get_print_timer_option")
      use, intrinsic :: ISO_C_BINDING
      implicit none
    end function

    subroutine dlai_set_print_timer_option(print_timer) &
        &bind(C, name="dlai_set_print_timer_option")
      use, intrinsic :: ISO_C_BINDING
      implicit none
      integer(C_INT), intent(in) :: print_timer
    end subroutine
  end interface

  contains

    function c_str(f_string) result(c_string)
        use, intrinsic :: ISO_C_BINDING
        implicit none
        character(len=*), intent(in)  :: f_string
        character(len=1, kind=C_CHAR) :: c_string(len_trim(f_string) + 1)
        integer i
        do i = 1, len_trim(f_string)
            c_string(i) = f_string(i:i)
        end do
        c_string(len_trim(f_string) + 1) = C_NULL_CHAR
    end function c_str

#define DLA_FTN_CHOLESKY_FACTORIZATION(function_name, function_name_aux, FtnType, CType) \
    subroutine function_name(uplo, n, a, ia, ja, desca, solver, info)                   ;\
      use, intrinsic :: ISO_C_BINDING                                                   ;\
      implicit none                                                                     ;\
      character(len=1), intent(in) :: uplo                                              ;\
      integer, intent(in) :: n                                                          ;\
      integer, dimension(*), intent(in) :: desca                                        ;\
      FtnType, dimension(desca(9),*), intent(inout) :: a                                ;\
      integer, intent(in) :: ia, ja                                                     ;\
      character(len=*), intent(in) :: solver                                            ;\
      integer, intent(out) :: info                                                      ;\
      interface                                                                         ;\
        integer(C_INT) function function_name_aux(uplo, n, a, ia, ja, desca, solver)     \
            bind(C, name="function_name")                                               ;\
          use, intrinsic :: ISO_C_BINDING                                               ;\
          implicit none                                                                 ;\
          character(len=1, kind=C_CHAR), intent(in) :: uplo                             ;\
          integer(C_INT), intent(in) :: n                                               ;\
          CType, intent(inout) :: a                                                     ;\
          integer(C_INT), intent(in) :: ia, ja                                          ;\
          integer(C_INT), intent(in) :: desca                                           ;\
          character(len=1, kind=C_CHAR), dimension(*), intent(in) :: solver             ;\
        end function                                                                    ;\
      end interface                                                                     ;\
      info = function_name_aux(uplo, n, a(1, 1), ia, ja, desca(1), c_str(solver))       ;\
    end subroutine

    DLA_FTN_CHOLESKY_FACTORIZATION(dlai_s_cholesky_factorization, dlai_s_cholesky_factorization_internal, real(4), real(C_FLOAT))
    DLA_FTN_CHOLESKY_FACTORIZATION(dlai_d_cholesky_factorization, dlai_d_cholesky_factorization_internal, real(8), real(C_DOUBLE))
    DLA_FTN_CHOLESKY_FACTORIZATION(dlai_c_cholesky_factorization, dlai_c_cholesky_factorization_internal, complex(4), complex(C_FLOAT_COMPLEX))
    DLA_FTN_CHOLESKY_FACTORIZATION(dlai_z_cholesky_factorization, dlai_z_cholesky_factorization_internal, complex(8), complex(C_DOUBLE_COMPLEX))

#define DLA_FTN_LU_FACTORIZATION(function_name, function_name_aux, FtnType, CType) \
    subroutine function_name(m, n, a, ia, ja, desca, solver, ipiv, info)                ;\
      use, intrinsic :: ISO_C_BINDING                                                   ;\
      implicit none                                                                     ;\
      integer, intent(in) :: m                                                          ;\
      integer, intent(in) :: n                                                          ;\
      integer, dimension(*), intent(in) :: desca                                        ;\
      FtnType, dimension(desca(9),*), intent(inout) :: a                                ;\
      integer, intent(in) :: ia, ja                                                     ;\
      integer, dimension(:), intent(in) :: ipiv                                         ;\
      character(len=*), intent(in) :: solver                                            ;\
      integer, intent(out) :: info                                                      ;\
      interface                                                                         ;\
        integer(C_INT) function function_name_aux(m, n, a, ia, ja, desca, ipiv, solver)  \
            bind(C, name="function_name")                                               ;\
          use, intrinsic :: ISO_C_BINDING                                               ;\
          implicit none                                                                 ;\
          integer(C_INT), intent(in) :: m                                               ;\
          integer(C_INT), intent(in) :: n                                               ;\
          CType, intent(inout) :: a                                                     ;\
          integer(C_INT), intent(in) :: ia, ja                                          ;\
          integer(C_INT), intent(in) :: desca                                           ;\
          integer(C_INT), intent(in) :: ipiv                                            ;\
          character(len=1, kind=C_CHAR), dimension(*), intent(in) :: solver             ;\
        end function                                                                    ;\
      end interface                                                                     ;\
      info = function_name_aux(m, n, a(1, 1), ia, ja, desca(1), ipiv(1), c_str(solver)) ;\
    end subroutine

    DLA_FTN_LU_FACTORIZATION(dlai_s_lu_factorization, dlai_s_lu_factorization_internal, real(4), real(C_FLOAT))
    DLA_FTN_LU_FACTORIZATION(dlai_d_lu_factorization, dlai_d_lu_factorization_internal, real(8), real(C_DOUBLE))
    DLA_FTN_LU_FACTORIZATION(dlai_c_lu_factorization, dlai_c_lu_factorization_internal, complex(4), complex(C_FLOAT_COMPLEX))
    DLA_FTN_LU_FACTORIZATION(dlai_z_lu_factorization, dlai_z_lu_factorization_internal, complex(8), complex(C_DOUBLE_COMPLEX))

#define DLA_FTN_MATRIX_MULTIPLICATION(function_name, function_name_aux, FtnType, CType)  \
    subroutine function_name(trans_a, trans_b, m, n, k, alpha, a, ia, ja, desca,         \
        b, ib, jb, descb, beta, c, ic, jc, descc, solver)                               ;\
      use, intrinsic :: ISO_C_BINDING                                                   ;\
      implicit none                                                                     ;\
      character(len=1), intent(in) :: trans_a, trans_b                                  ;\
      integer, intent(in) :: m, n, k                                                    ;\
      FtnType, intent(in) :: alpha, beta                                                ;\
      integer, dimension(*), intent(in) :: desca, descb, descc                          ;\
      FtnType, dimension(desca(9),*), intent(in) :: a                                   ;\
      FtnType, dimension(descb(9),*), intent(in) :: b                                   ;\
      FtnType, dimension(descc(9),*), intent(inout) :: c                                ;\
      integer, intent(in) :: ia, ja, ib, jb, ic, jc                                     ;\
      character(len=*), intent(in) :: solver                                            ;\
      integer :: i                                                                      ;\
      interface                                                                         ;\
        integer(C_INT) function function_name_aux(trans_a, trans_b, m, n, k, alpha,      \
          a, ia, ja, desca, b, ib, jb, descb, beta, c, ic, jc, descc, solver)            \
            bind(C, name="function_name")                                               ;\
          use, intrinsic :: ISO_C_BINDING                                               ;\
          implicit none                                                                 ;\
          character(len=1, kind=C_CHAR), intent(in) :: trans_a, trans_b                 ;\
          integer(C_INT), intent(in) :: m, n, k                                         ;\
          CType, intent(in) :: alpha, beta                                              ;\
          CType, intent(in) :: a, b                                                     ;\
          CType, intent(inout) :: c                                                     ;\
          integer(C_INT), intent(in) :: ia, ja, ib, jb, ic, jc                          ;\
          integer(C_INT), intent(in) :: desca, descb, descc                             ;\
          character(len=1, kind=C_CHAR), dimension(*), intent(in) :: solver             ;\
        end function                                                                    ;\
      end interface                                                                     ;\
      i = function_name_aux(trans_a, trans_b, m, n, k, alpha, a(1, 1), ia, ja, desca(1), \
     b(1, 1), ib, jb, descb(1), beta, c(1, 1), ic, jc, descc(1), c_str(solver))         ;\
    end subroutine

    DLA_FTN_MATRIX_MULTIPLICATION(dlai_s_matrix_multiplication, dlai_s_matrix_multiplication_internal, real(4), real(C_FLOAT))
    DLA_FTN_MATRIX_MULTIPLICATION(dlai_d_matrix_multiplication, dlai_d_matrix_multiplication_internal, real(8), real(C_DOUBLE))
    DLA_FTN_MATRIX_MULTIPLICATION(dlai_c_matrix_multiplication, dlai_c_matrix_multiplication_internal, complex(4), complex(C_FLOAT_COMPLEX))
    DLA_FTN_MATRIX_MULTIPLICATION(dlai_z_matrix_multiplication, dlai_z_matrix_multiplication_internal, complex(8), complex(C_DOUBLE_COMPLEX))

end module
