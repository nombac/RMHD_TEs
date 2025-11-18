
program main

  use bval_module, only : tcons
  use cfl_module, only : cfl
  use digit_module, only : digit3
  use dump_module, only : dump_write
  use floor_module, only : floor_v3
  use initial_module, only : initial
  use magnetic_module, only : magnetic
#ifdef MPI
  use mpi_module
#endif
  use radiation_module, only : fld
  use source_module, only : source
  use transport_module, only : transport
  use viscous_module, only : viscous

  implicit none

  real(8) :: t1, t2, elapsed

#ifdef MPI
  call mpi_begin
#endif

#ifdef MPI
  call mpi_barrier(mpi_comm_world, ier)
  if(myrk == 0) then
     t1 = mpi_wtime()
  end if
#else
  call cpu_time(t1)
#endif
  call initial

  call print_status(1)
  do 
     call fld

     call cfl
#ifdef SHEAR
     call tcons
#endif
     call source
     call viscous
     call transport
     call magnetic
     call floor_v3
#ifdef MPI
     call mpi_barrier(mpi_comm_world, ier)
     if(myrk == 0) then
        t2 = mpi_wtime()
     end if
     elapsed = t2 - t1
     call mpi_bcast(elapsed, 1, mpi_double, 0, mpi_comm_world, ier)
#else
     call cpu_time(t2)
     elapsed = t2 - t1
#endif
     if(intchk(elapsed)) exit
  enddo
  call print_status(0)

#ifdef MPI
  call mpi_end
#endif

  stop

contains




  logical function intchk(cputime)
    ! check stopping criteria
    use constant_module, only : PI
    use root_module, only : nhy, nlim, time, tlim, dt, dtmin, cpulim, &
                            dumpfile, histfile, dtdump, dthist, tdump, thist
    use param_module, only : omega
    use hist_module, only : hist_write
    use dump_module, only : dump_write
    real(8), intent(in) :: cputime
    integer, parameter :: NHY0 = 1000
    integer, parameter :: NHYA = 100
    integer :: iic
#ifdef VPP
    integer :: tvu
#endif
#ifdef MPI
    integer :: iicm
#endif  
    real(8) :: torbit
    integer :: iabort = 0

    torbit = (2d0 * PI / omega)

    intchk = .false.

    nhy = nhy + 1
    time = time + dt
#ifdef VPP
    call clockv(tvu, iic, 0, 0)
#else
    iic = int(cputime)
#endif
!!$#ifdef MPI
!!$    call mpi_allreduce(iic, iicm, 1, mpi_integer, mpi_max, mpi_comm_world, ier)
!!$    iic = iicm
!!$#endif

    if(time >= tdump + dtdump) then
       tdump = tdump + dtdump
       call dump_write(dumpfile)
    endif
    if(time >= thist + dthist) then
       thist = thist + dthist
       call hist_write(histfile)
    endif

    if(mod(nhy, NHY0) == 0) then
#ifdef MPI
       if(myrk == 0) then
#endif
       print *, 'nhy = ', nhy, &
                'time = ', real(time / torbit), &
                'dt = ', real(dt / torbit)
       open(90, file = 'cputime.data')
       write(90,*) 'cputime = ', iic, '[s] ', &
            real(iic)/real(cpulim)*100., ' %', &
            real(cpulim - iic)/60./60., ' hour left'
       close(90)
#ifdef MPI
       end if
#endif
    end if

    ! (1) elapse time
    if(iic > cpulim) then
       intchk = .true.
       print *, '--- elapsed time = ', iic, 'exceeded limit ', cpulim, ' ---'
    endif
    ! (2) dt
    if(dt < dtmin) then
       intchk = .true.
       print *, '### dt = ', real(dt / torbit), &
                ' < dtmin = ', real(dtmin / torbit), ' ###'
    endif
    ! (3) time
    if(time >= tlim) then
       intchk = .true.
       print *, '--- time = ', real(time / torbit), &
                ' >= tlim = ', real(tlim / torbit), ' ---'
    endif
    ! (4) nhy
    if(nhy >= nlim) then
       intchk = .true.
       print *, '--- nhy = ', nhy, ' >= nlim = ', nlim, ' ---'
    endif
    ! (5) abort
    if(mod(nhy, NHYA) == 0) then
#ifdef MPI
       if(myrk == 0) then
#endif
       open(99, file = 'abort.data', status = 'old')
       read(99,*) iabort
       close(99)
#ifdef MPI
       end if
       call mpi_bcast(iabort, 1, mpi_integer, 0, mpi_comm_world, ier)
#endif
       if(iabort == 1) then
          intchk = .true.
          print *, '*** program aborted by user ***'
       end if
    end if

    if(intchk) then
       call dump_write('restart.dat')
    endif

    return
  end function intchk




  subroutine print_status(istat)
    integer, intent(in) :: istat
#ifdef MPI
    if(myrk == 0) then
#endif
    open(99, file = 'status.data')
    write(99,*) istat
    close(99)
#ifdef MPI
    endif
#endif
    return
  end subroutine print_status



  
end program main
