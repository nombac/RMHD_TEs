
module dump_module

  implicit none

  private
  public :: dump_read, dump_write

contains




  subroutine dump_write(filename)
    use constant_module, only : PI
    use field_module, only : d, e, v1, v2, v3, b1, b2, b3, er, qmag, qkin, scol
    use flux_module, only : drhox, dadd
    use grid_module, only : in, jn, kn, x1b, x2b, x3b
#ifdef MPI
    use gather_module, only : gather
    use mpi_module, only : myrk
#endif
    use param_module, only : omega, dfloor
    use root_module, only : dtqq, time, nhy, tdump, thist, dumpfile, histfile
    use strtoi_module, only : strtoi
    character, intent(in) :: filename*11
    integer :: incr
    real(8) :: daddx

#ifdef MPI
    call gather(d,0)
    call gather(e,0)
    call gather(v1,1)
    call gather(v2,2)
    call gather(v3,3)
    call gather(b1,1)
    call gather(b2,0)
    call gather(b3,3)
    call gather(er,0)
#ifdef PARTETOTK
    call gather(qkin,4)
    call gather(scol,4)
#endif
#ifdef PARTETOT
    call gather(qmag,4)
#endif
#endif
    call cumul(daddx, drhox)
    dadd = dadd + daddx

    print *, 'restart dump written at time/orbit =', &
         time / (2d0 * PI / omega), 'cycle =', nhy, 'file =', filename

#ifdef MPI
    if(myrk == 0) then
#endif
    open(4, file = 'data/'//filename, status = 'unknown', form = 'unformatted')
    if(filename /= 'restart.dat') then
       incr = strtoi(dumpfile,4,6) + 1
       write(dumpfile, '(a3,i3.3,a5)') 'res', incr, '.data'
    endif
    write(4) in, jn, kn
    write(4) x1b, x2b, x3b
    write(4) d, e, v1, v2, v3, b1, b2, b3, er, qmag, qkin, scol
    write(4) dtqq, time, nhy, tdump, thist, dumpfile, histfile
    write(4) dfloor
    write(4) dadd
    close(4)

    open(4, file='data/resol.data')
    write(4,*) in, jn, kn
    close(4)
    open(4, file='data/grid-x.data', form='unformatted')
    write(4) real(x1b)
    close(4)
    open(4, file='data/grid-y.data', form='unformatted')
    write(4) real(x2b)
    close(4)
    open(4, file='data/grid-z.data', form='unformatted')
    write(4) real(x3b)
    close(4)
#ifdef MPI
    endif
#endif
    
    return
  end subroutine dump_write




  subroutine dump_read(filename, targ)
    use constant_module, only : PI
    use digit_module, only : digit3
    use field_module, only : d, e, v1, v2, v3, b1, b2, b3, er, qmag, qkin, scol
    use flux_module, only : dadd
    use grid_module, only : in, jn, kn, x1b, x2b, x3b
    use param_module, only : omega, dfloor
    use root_module, only : dtqq, time, nhy, tdump, thist, dumpfile, histfile
    use strtoi_module, only : strtoi
    character, intent(inout) :: filename*11
    character, intent(in) :: targ*5
    integer :: incr
    integer :: istat, ios, nfile
    integer, parameter :: nfile_max = 999

    open(99, file = 'status.data', status = 'old')
    read(99,*) istat
    close(99)

    if(istat /= 0 .and. filename == 'restart.dat') then
       do nfile = nfile_max, 0, -1
          filename = 'res'//digit3(nfile)//'.data'
          open(4, file = '../'//targ//'/data/'//filename, status = 'old', &
               & form = 'unformatted', iostat = ios)
          if(ios == 0) exit
          if(nfile == 0) then
             stop 'abnormal end at dump_read'
          endif
       enddo
    else
       open(4, file = '../'//targ//'/data/'//filename, status = 'old', &
            & form = 'unformatted')
    endif
    read(4) in, jn, kn
    read(4) x1b, x2b, x3b
    read(4) d, e, v1, v2, v3, b1, b2, b3, er, qmag, qkin, scol
    read(4) dtqq, time, nhy, tdump, thist, dumpfile, histfile
    read(4) dfloor
    read(4) dadd
    close(4)

    print *, 'restart dump read at time/orbit =', &
         time / (2d0 * PI / omega), 'cycle =', nhy, 'file =', filename

    incr = strtoi(histfile,4,6) + 1
    write(histfile, '(a3,i3.3,a5)') 'hst', incr, '.data'

    return
  end subroutine dump_read




  subroutine cumul(ax, a)
    use grid_module, only : is, ie, js, je, ks, ke, lx0, ly0, lz0
#ifdef MPI
    use mpi_module
#endif
    real(8), intent(out) :: ax
    real(8), intent(in) :: a(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    integer :: i, j, k
#ifdef MPI
    real(8) :: ax0
#endif

    ax = 0d0
    do k = ks, ke
    do j = js, je
    do i = is, ie
       ax = ax + a(i,j,k)
    enddo
    enddo
    enddo
#ifdef MPI
    call mpi_reduce(ax, ax0, 1, mpi_double_precision, mpi_sum, 0, &
         mpi_comm_world, ier)
    ax = ax0
#endif
    ax = ax / (lx0 * ly0 * lz0)

    return
  end subroutine cumul




end module dump_module
