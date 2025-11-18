module check_module

#ifdef MPI  
  use mpi_module
#endif

  implicit none

  integer, parameter :: nchck = 10
  real(4) :: tchck(1:nchck) = 0.
  real(4) :: tchck0(1:nchck) = 0.
  real(4) :: t1, t2, tarr(2)
  integer :: icheck = 0

contains

  subroutine check_total

#ifdef MPI
    call mpi_reduce(tchck, tchck0, nchck, mpi_real, mpi_max, &
         0, mpi_comm_world, ier)
#endif
    tchck(1:nchck) = tchck0(1:nchck)

#ifdef MPI
    if(myrk == 0) then
#endif
       open(20, file = 'check.data', form = 'unformatted')
       write(20) tchck
       close(20)
#ifdef MPI
    endif
#endif

    return
  end subroutine check_total

end module check_module
