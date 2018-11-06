module test_ftn_dlai
  use dla_interface
  use mpi
  use, intrinsic :: ISO_C_BINDING
  implicit none

contains

  subroutine test_ftn_dlai_setup(icomms, set_comms) bind (C, name="test_ftn_dlai_setup")
    use, intrinsic :: ISO_C_BINDING
    implicit none
    integer(C_INT), dimension(8), intent(out) :: icomms
    integer(C_INT), intent(out) :: set_comms

    integer comm_size
    integer comm_rank
    integer ierr

    call dlai_initialize(2, 1);

    call MPI_Comm_size(MPI_COMM_WORLD, comm_size, ierr);
    if (comm_size .ne. 6) then
      print *, "This test need 6 MPI ranks (", comm_size,  " provided)!\n"
      stop 1
    endif

    call MPI_Comm_rank(MPI_COMM_WORLD, comm_rank, ierr);

    ! Create communicators used in the tests.
    icomms(1) = dlai_create_2d_grid(MPI_COMM_WORLD, 2, 3, "R");
    icomms(2) = dlai_create_2d_grid(MPI_COMM_WORLD, 3, 2, "R");
    icomms(3) = dlai_create_2d_grid(MPI_COMM_WORLD, 1, 6, "R");
    icomms(4) = dlai_create_2d_grid(MPI_COMM_WORLD, 2, 3, "C");
    icomms(5) = dlai_create_2d_grid(MPI_COMM_WORLD, 3, 2, "C");
    icomms(6) = dlai_create_2d_grid(MPI_COMM_WORLD, 1, 6, "C");
    set_comms = 6;
  end subroutine

  subroutine test_ftn_dlai_end() bind (C, name="test_ftn_dlai_end")
    call dlai_finalize();
  end subroutine
end module
