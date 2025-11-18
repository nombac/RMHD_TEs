
module grid_module

  implicit none
  real(8), allocatable :: x1a(:), x1b(:), dx1a(:), dx1b(:)
  real(8), allocatable :: x2a(:), x2b(:), dx2a(:), dx2b(:)
  real(8), allocatable :: x3a(:), x3b(:), dx3a(:), dx3b(:)
  real(8) :: dx, dy, dz
  real(8) :: lx0, ly0, lz0
  integer :: in, jn, kn
  integer :: is, ie, is0, ie0
  integer :: js, je, js0, je0
  integer :: ks, ke, ks0, ke0
  integer :: nx1z, nx2z, nx3z

  public :: grid_initial, para_range

contains
  



  subroutine grid_initial
#ifdef MPI
    use mpi_module
#endif
    real(8) :: x1min, x1max, x2min, x2max, x3min, x3max
    real(8) :: h0
    integer :: in0, jn0, kn0
    integer :: i, j, k
    namelist /ggen1/ in0, x1min, x1max
    namelist /ggen2/ jn0, x2min, x2max
    namelist /ggen3/ kn0, x3min, x3max
    namelist /lnrmcon/ h0

    ! read simulation parameters
    open(1, file = 'z3dinput', status = 'old')
    read(1, ggen1) ! in0, x1min, x1max
    read(1, ggen2) ! jn0, x2min, x2max
    read(1, ggen3) ! kn0, x3min, x3max 
    close(1)

    ! read normalization length
    open(1,file = 'init/iparam.data', status = 'old')
    read(1, lnrmcon) ! h0
    close(1)

#ifdef SHEAR
    x1min = x1min * h0
    x1max = x1max * h0
    x2min = x2min * h0
    x2max = x2max * h0
    x3min = x3min * h0
    x3max = x3max * h0
#endif

    !-----------  GENERATE X1 GRID  ----------------------------------------
    in = in0 + 5
    allocate(x1a(in+1), x1b(in), dx1a(in), dx1b(in))
    
    is0 = 3
    ie0 = is0 + in0 - 1
    nx1z = ie0 - is0 + 1
    dx = (x1max - x1min) / float(in0)

    ! A-grid
    do i = is0, ie0
       dx1a(i) = dx
    enddo
    dx1a(is0-1) = dx1a(is0)
    dx1a(is0-2) = dx1a(is0)
    dx1a(ie0+1) = dx1a(ie0)
    dx1a(ie0+2) = dx1a(ie0)

    x1a(is0) = x1min
    do i = is0 + 1, ie0+1
       x1a(i) =  x1a(i-1) + dx1a(i-1)
    enddo
    x1a(is0-1) = x1a(is0  ) - dx1a(is0-1)
    x1a(is0-2) = x1a(is0-1) - dx1a(is0-2)
    x1a(ie0+2) = x1a(ie0+1) + dx1a(ie0+1)
    x1a(ie0+3) = x1a(ie0+2) + dx1a(ie0+2)

    ! B-grid
    x1b(is0-2) = x1a(is0-1) - 0.5 * dx1a(is0-2)
    do i = is0-1, ie0+2
       x1b(i) = x1a(i) + 0.5 * dx1a(i)
    enddo
    x1b(ie0+3) = x1a(ie0+3) + 0.5 * dx1a(ie0+2)

    dx1b(is0-2) = dx1a(is0-2)
    do i = is0-1, ie0+2
       dx1b(i) = x1b(i) - x1b(i-1)
    enddo

    lx0 = x1a(ie0+1) - x1a(is0)

    !-----------  GENERATE X2 GRID  ----------------------------------------
    jn = jn0 + 5
    allocate(x2a(jn+1), x2b(jn), dx2a(jn), dx2b(jn))

    js0 = 3
    je0 = js0 + jn0 - 1
    nx2z = je0 - js0 + 1
    dy = (x2max - x2min) / float(jn0)

    ! A-grid
    do j = js0, je0
       dx2a(j) = dy
    enddo
    dx2a(js0-1) = dx2a(js0)
    dx2a(js0-2) = dx2a(js0)
    dx2a(je0+1) = dx2a(je0)
    dx2a(je0+2) = dx2a(je0)

    x2a(js0) = x2min
    do j = js0 + 1, je0+1
       x2a(j) =  x2a(j-1) + dx2a(j-1)
    enddo
    x2a(js0-1) = x2a(js0  ) - dx2a(js0-1)
    x2a(js0-2) = x2a(js0-1) - dx2a(js0-2)
    x2a(je0+2) = x2a(je0+1) + dx2a(je0+1)
    x2a(je0+3) = x2a(je0+2) + dx2a(je0+2)

    ! B-grid
    x2b(js0-2) = x2a(js0-1) - 0.5 * dx2a(js0-2)
    do j = js0-1, je0+2
       x2b(j) = x2a(j) + 0.5 * dx2a(j)
    enddo
    x2b(je0+3) = x2a(je0+3) + 0.5 * dx2a(je0+2)

    dx2b(js0-2) = dx2a(js0-2)
    do j = js0-1, je0+2
       dx2b(j) = x2b(j) - x2b(j-1)
    enddo

    ly0 = x2a(je0+1) - x2a(js0)

    !-----------  GENERATE X3 GRID  ----------------------------------------
    kn = kn0 + 5
    allocate(x3a(kn+1), x3b(kn), dx3a(kn), dx3b(kn))

    ks0 = 3
    ke0 = ks0 + kn0 - 1
    nx3z = ke0 - ks0 + 1
    dz = (x3max - x3min) / float(kn0)

    ! A-grid
    do k = ks0, ke0
       dx3a(k) = dz
    enddo
    dx3a(ks0-1) = dx3a(ks0)
    dx3a(ks0-2) = dx3a(ks0)
    dx3a(ke0+1) = dx3a(ke0)
    dx3a(ke0+2) = dx3a(ke0)

    x3a(ks0) = x3min
    do k = ks0 + 1, ke0+1
       x3a(k) =  x3a(k-1) + dx3a(k-1)
    enddo
    x3a(ks0-1) = x3a(ks0  ) - dx3a(ks0-1)
    x3a(ks0-2) = x3a(ks0-1) - dx3a(ks0-2)
    x3a(ke0+2) = x3a(ke0+1) + dx3a(ke0+1)
    x3a(ke0+3) = x3a(ke0+2) + dx3a(ke0+2)

    ! B-grid
    x3b(ks0-2) = x3a(ks0-1) - 0.5 * dx3a(ks0-2)
    do k = ks0-1, ke0+2
       x3b(k) = x3a(k) + 0.5 * dx3a(k)
    enddo
    x3b(ke0+3) = x3a(ke0+3) + 0.5 * dx3a(ke0+2)

    dx3b(ks0-2) = dx3a(ks0-2)
    do k = ks0-1, ke0+2
       dx3b(k) = x3b(k) - x3b(k-1)
    enddo

    lz0 = x3a(ke0+1) - x3a(ks0)


    is = is0
    ie = ie0
    js = js0
    je = je0
    ks = ks0
    ke = ke0

#ifdef MPI
    call para_range(is0, ie0, iprc, myrki, is, ie)
    call para_range(js0, je0, jprc, myrkj, js, je)
    call para_range(ks0, ke0, kprc, myrkk, ks, ke)
#endif

    return

  end subroutine grid_initial




  subroutine para_range(n1, n2, nprocs, irank, ista, iend)

    integer, intent(in) :: n1
    integer, intent(in) :: n2
    integer, intent(in) :: nprocs
    integer, intent(in) :: irank
    integer, intent(out) :: ista
    integer, intent(out) :: iend

    integer :: iwork1, iwork2

    iwork1 = (n2 - n1 + 1) / nprocs
    iwork2 = mod(n2 - n1 + 1, nprocs)
    ista = irank * iwork1 + n1 + min(irank, iwork2)
    iend = ista+iwork1-1
    if(iwork2 > irank) then
       iend = iend + 1
    endif

    return

  end subroutine para_range




end module grid_module
