
#define MAXMIN min

#define ITERMAX 100
#define EPS 1d-6

#define ITERMAX_SLVSML 10000
#define EPS_SLVSML 1d-3

#define NCYCLE 1
#define NPRE   5
#define NPOST  5

module diff_module

#ifdef MPI
  use mpi_module
#endif

  implicit none

  private
  public :: diff

  real(8) :: lx, ly, lz, dt

  integer :: ng
#ifdef MPI
  integer :: ngx
#endif

  integer :: memlen, mem
  real(8), allocatable :: z(:)
  integer, allocatable :: ires(:), irho(:), irhs(:), iu(:)

  integer :: memlen0, mem0
  real(8), allocatable :: z0(:)
  integer, allocatable :: ido(:)

  integer :: memlen1, mem1
  real(8), allocatable :: z1(:)
  integer, allocatable :: ibo(:), ibi(:)

  integer, allocatable :: nx(:), ny(:), nz(:)
  integer, allocatable :: iq(:), jq(:), kq(:)

contains




  subroutine diff(u, d, i0, j0, k0, lx0, ly0, lz0, dt0, iter)

    use param_module, only : HUGE

    integer, intent(in) :: i0
    integer, intent(in) :: j0
    integer, intent(in) :: k0
    real(8), intent(inout) :: u(0:i0+1,0:j0+1,0:k0+1)
    real(8), intent(in) ::    d(0:i0+1,0:j0+1,0:k0+1)
    real(8), intent(in) :: lx0
    real(8), intent(in) :: ly0
    real(8), intent(in) :: lz0
    real(8), intent(in) :: dt0
    integer, intent(out) :: iter

    real(8) ::   v(0:i0+1,0:j0+1,0:k0+1)
    real(8) :: rhs(0:i0+1,0:j0+1,0:k0+1)
    real(8) :: res(0:i0+1,0:j0+1,0:k0+1)
    real(8) :: anorm
    
    lx = lx0
    ly = ly0
    lz = lz0
    dt = dt0

    call initial(i0, j0, k0)
    call coarsekb
    call coarsed(d, i0, j0, k0)
    call bndcz(u, i0, j0, k0, ng, 0)
    call bndcx(u, i0, j0, k0, ng)
    call bndcy(u, i0, j0, k0, ng)
    call src(rhs, u, i0, j0, k0)
    call resid(res, u, rhs, d, i0, j0, k0, ng)

    do iter = 1, ITERMAX
       call mglin(v, res, i0, j0, k0)
       call add(u, v, i0, j0, k0)
       call resid(res, u, rhs, d, i0, j0, k0, ng)
       call norm(anorm, res, rhs, i0, j0, k0, ng)
       if(anorm < EPS) then
          call final
          return ! normal exit
       end if
       if(anorm > HUGE) exit
    end do
    
    print *, 'diff: not converge: anorm =', anorm, ' iter =', iter
#ifdef MPI
    call mpi_finalize(ier)
#endif
    stop

    return

  end subroutine diff




  subroutine initial(i0, j0, k0)

    integer, intent(in) :: i0
    integer, intent(in) :: j0
    integer, intent(in) :: k0

    integer :: n
    integer :: array3d, array2d
    integer :: i, j, k
#ifdef MPI
    integer :: isw
#endif

    ! number of multi-grid levels
#ifdef MPI
    ng = log20(MAXMIN((i0 * iprc), (j0 * jprc), (k0 * kprc))) + 1
#else
    ng = log20(MAXMIN(i0, j0, k0)) + 1
#endif

    ! coarsening method
    allocate(nx(1:ng), ny(1:ng), nz(1:ng))
    allocate(iq(1:ng), jq(1:ng), kq(1:ng))
    nx(1:ng) = 1
    ny(1:ng) = 1
    nz(1:ng) = 1
    i = i0
    j = j0
    k = k0
#ifdef MPI
    isw = 0
    ngx = 1 ! default gathering/scattering level
#endif
    do n = ng, 1, -1
#ifdef MPI
       if((isw == 0) .and. ((i == 1) .or. (j == 1) .or. (k == 1))) then
          ngx = n ! redefine gathering/scattering level
          i = i * iprc
          j = j * jprc
          k = k * kprc
          isw = 1
       endif
#endif
       iq(n) = max(i, 1)
       jq(n) = max(j, 1)
       kq(n) = max(k, 1)
       if(i == 0) nx(n) = 0
       if(j == 0) ny(n) = 0
       if(k == 0) nz(n) = 0
       i = i / 2
       j = j / 2
       k = k / 2
    end do
#ifdef MPI
    if(ngx >= (ng-1)) then
       print *, 'ngx =', ngx,' >= (ng-1) =', ng-1
       call mpi_finalize(ier)
       stop
    end if
#endif

    ! memory allocation for 3d array (1)
    memlen = 0
    memlen0 = 0
    do n = ng, 1, -1
       i = iq(n)
       j = jq(n)
       k = kq(n)
       array3d = (i + 2) * (j + 2) * (k + 2)
       memlen = memlen + (array3d + 1) ! iu
       memlen = memlen + (array3d + 1) ! irhs
       if(n >= 2) then
          memlen = memlen + (array3d + 1) ! ires
       end if
       if(n <= ng-1) then
          memlen = memlen + (array3d + 1) ! irho
       end if
       memlen0 = memlen0 + (array3d + 1) ! id
    end do
    allocate(z(memlen))
    allocate(z0(memlen0))
    z(1:memlen) = 0d0
    z0(1:memlen0) = 0d0
    ! memory allocation for 3d array (2)
    allocate(iu(ng), irhs(ng), ires(ng), irho(ng))
    allocate(ido(ng))
    mem = 0
    mem0 = 0
    do n = ng, 1, -1
       i = iq(n)
       j = jq(n)
       k = kq(n)
       array3d = (i + 2) * (j + 2) * (k + 2)
       iu(n)   =  maloc(array3d)
       irhs(n) =  maloc(array3d)
       if(n >= 2) then
          ires(n) =  maloc(array3d)
       end if
       if(n <= ng - 1) then
          irho(n) =  maloc(array3d)
       end if
       ido(n)  = maloc0(array3d)
    end do

    ! memory allocation for 2d array (1)
    memlen1 = 0
    do n = ng, 1, -1
       i = iq(n)
       j = jq(n)
       array2d = i * j
       memlen1 = memlen1 + (array2d + 1) ! ibo
       memlen1 = memlen1 + (array2d + 1) ! ibi
    end do
    allocate(z1(memlen1))
    z1(1:memlen1) = 0d0
    ! memory allocation for 2d array (2)
    allocate(ibo(ng), ibi(ng))
    mem1 = 0
    do n = ng, 1, -1
       i = iq(n)
       j = jq(n)
       array2d = i * j
       ibo(n) = maloc1(array2d)
       ibi(n) = maloc1(array2d)
    end do

    return

  end subroutine initial




  subroutine final

    deallocate(z)
    deallocate(ires, irho, irhs, iu)
    
    deallocate(z0)
    deallocate(ido)
    
    deallocate(z1)
    deallocate(ibo, ibi)

    deallocate(nx, ny, nz)
    deallocate(iq, jq, kq)
    
    return

  end subroutine final




  subroutine coarsed(d, i0, j0, k0)

    real(8), intent(in) :: d(0:i0+1,0:j0+1,0:k0+1)
    integer, intent(in) :: i0
    integer, intent(in) :: j0
    integer, intent(in) :: k0

    integer :: n

    n = ng
    call copy(z0(ido(n)), d, iq(n), jq(n), kq(n))

    do n = ng - 1, 1, -1
       call rstrct(z0(ido(n)), z0(ido(n + 1)), iq(n), jq(n), kq(n), &
            nx(n), ny(n), nz(n), n)
       call bndcz(z0(ido(n)), iq(n), jq(n), kq(n), n, 1)
       call bndcx(z0(ido(n)), iq(n), jq(n), kq(n), n)
       call bndcy(z0(ido(n)), iq(n), jq(n), kq(n), n)
    end do

    return

  end subroutine coarsed




  subroutine src(rhs,u,i0,j0,k0)

    integer, intent(in) :: i0
    integer, intent(in) :: j0
    integer, intent(in) :: k0
    real(8), intent(out) :: rhs(0:i0+1,0:j0+1,0:k0+1)
    real(8), intent(in)  ::   u(0:i0+1,0:j0+1,0:k0+1)

    integer :: i, j, k
    integer :: l, lm, im, imjm, ista, iend, jsta, jend, ksta, kend

    ista = 0; iend = i0 + 1
    jsta = 0; jend = j0 + 1
    ksta = 0; kend = k0 + 1

    im =   (iend - ista + 1)
    imjm = (jend - jsta + 1) * im
    lm =   (kend - ksta + 1) * imjm

!CDIR LOOPCNT=LOOP_CNT
!OCL NOVREC
!OCL REPEAT(LOOP_CNT)
    do l = 0, lm - 1
       k = (l / imjm) + ksta
       j =((l - (k - ksta) * imjm) / im) + jsta
       i = (l - (k - ksta) * imjm - (j - jsta) * im) + ista
       rhs(i,j,k) = - u(i,j,k) / dt
    end do
    
    return

  end subroutine src




  subroutine add(u, v, i0, j0, k0)

    integer, intent(in) :: i0
    integer, intent(in) :: j0
    integer, intent(in) :: k0
    real(8), intent(inout) :: u(0:i0+1,0:j0+1,0:k0+1)
    real(8), intent(in)    :: v(0:i0+1,0:j0+1,0:k0+1)

    integer :: i, j, k
    integer :: l, lm, im, imjm, ista, iend, jsta, jend, ksta, kend

    ista = 0; iend = i0 + 1
    jsta = 0; jend = j0 + 1
    ksta = 0; kend = k0 + 1

    im =   (iend - ista + 1)
    imjm = (jend - jsta + 1) * im
    lm =   (kend - ksta + 1) * imjm

!CDIR LOOPCNT=LOOP_CNT
!OCL NOVREC
!OCL REPEAT(LOOP_CNT)
    do l = 0, lm - 1
       k = (l / imjm) + ksta
       j =((l - (k - ksta) * imjm) / im) + jsta
       i = (l - (k - ksta) * imjm - (j - jsta) * im) + ista
       u(i,j,k) = u(i,j,k) + v(i,j,k)
    end do
      
    return

  end subroutine add




  subroutine norm(anorm, res, rhs, i0, j0, k0, n)

    integer, intent(in) :: i0
    integer, intent(in) :: j0
    integer, intent(in) :: k0
    real(8), intent(out) :: anorm
    real(8), intent(in) :: res(0:i0+1,0:j0+1,0:k0+1)
    real(8), intent(in) :: rhs(0:i0+1,0:j0+1,0:k0+1)
    integer, intent(in) :: n

    integer :: i, j, k
    real(8) :: bunbo, bunsi
#ifdef MPI
    real(8) :: bunbo0, bunsi0
#endif
    integer :: l, lm, im, imjm, ista, iend, jsta, jend, ksta, kend

    ista = 1; iend = i0
    jsta = 1; jend = j0
    ksta = 1; kend = k0

    im =   (iend - ista + 1)
    imjm = (jend - jsta + 1) * im
    lm =   (kend - ksta + 1) * imjm

    bunbo = 0d0
    bunsi = 0d0
!CDIR LOOPCNT=LOOP_CNT
!OCL NOVREC
!OCL REPEAT(LOOP_CNT) 
   do l = 0, lm - 1
       k = (l / imjm) + ksta
       j =((l - (k - ksta) * imjm) / im) + jsta
       i = (l - (k - ksta) * imjm - (j - jsta) * im) + ista
       bunbo = bunbo + rhs(i,j,k)**2
       bunsi = bunsi + res(i,j,k)**2
    end do

#ifdef MPI
    if(n > ngx) then
       call mpi_allreduce(bunbo, bunbo0, 1, mpi_double_precision,&
            mpi_sum, mpi_comm_world, ier)
       call mpi_allreduce(bunsi, bunsi0, 1, mpi_double_precision,&
            mpi_sum, mpi_comm_world, ier)
       bunbo = bunbo0
       bunsi = bunsi0
    end if
#endif
    anorm = sqrt(bunsi / bunbo)

    return

  end subroutine norm




  subroutine mglin(u, rhs, i0, j0, k0)

    real(8), intent(inout) ::   u(0:i0+1,0:j0+1,0:k0+1)
    real(8), intent(in)    :: rhs(0:i0+1,0:j0+1,0:k0+1)
    integer, intent(in) :: i0
    integer, intent(in) :: j0
    integer, intent(in) :: k0

    integer :: l, m, n
    integer :: jcycle, jpost, jpre

    ! restriction down to the corasest level
    n = ng - 1
    call rstrct(z(irho(n)), rhs, iq(n), jq(n), kq(n), &
         nx(n), ny(n), nz(n), n)
    do n = ng - 2, 1, -1
       call rstrct(z(irho(n)), z(irho(n+1)), iq(n), jq(n), kq(n), &
            nx(n), ny(n), nz(n), n)
    end do

    ! solution at the coarsest level
    call slvsml(z(iu(1)), z(irho(1)), z0(ido(1)), iq(1), jq(1), kq(1))

    do l = 2, ng
       ! interpolation up to the next finer level
       call interp(z(iu(l)), z(iu(l-1)), iq(l), jq(l), kq(l), &
            nx(l-1), ny(l-1), nz(l-1), l-1)
       call bndcz(z(iu(l)), iq(l), jq(l), kq(l), l, 0)
       call bndcx(z(iu(l)), iq(l), jq(l), kq(l), l)
       call bndcy(z(iu(l)), iq(l), jq(l), kq(l), l)
       if(l /= ng) then
          call copy(z(irhs(l)), z(irho(l)), iq(l), jq(l), kq(l))
       else
          call copy(z(irhs(l)),        rhs, iq(l), jq(l), kq(l))
       end if

       do jcycle = 1, NCYCLE
          ! V-cycle downwards
          do m = l, 2, -1
             do jpre = 1, NPRE
                call relax(z(iu(m)), z(irhs(m)), z0(ido(m)), &
                     iq(m), jq(m), kq(m), m)
             end do
             call resid(z(ires(m)), z(iu(m)), z(irhs(m)), z0(ido(m)), &
                  iq(m), jq(m), kq(m), m)
             call rstrct(z(irhs(m-1)), z(ires(m)), iq(m-1), jq(m-1), kq(m-1), &
                  nx(m-1), ny(m-1), nz(m-1), m-1)
             call fill0(z(iu(m-1)), iq(m-1), jq(m-1), kq(m-1))
          end do

          ! solution at the coarsest level
          call slvsml(z(iu(1)), z(irhs(1)), z0(ido(1)), iq(1), jq(1), kq(1))


          ! V-cycle upwards
          do m = 2, l
             call interp(z(ires(m)), z(iu(m-1)), iq(m), jq(m), kq(m), &
                  nx(m-1), ny(m-1), nz(m-1), m-1)
             call add(z(iu(m)), z(ires(m)), iq(m), jq(m), kq(m))
             call bndcz(z(iu(m)), iq(m), jq(m), kq(m), m, 0)
             call bndcx(z(iu(m)), iq(m), jq(m), kq(m), m)
             call bndcy(z(iu(m)), iq(m), jq(m), kq(m), m)
             do jpost = 1, NPOST
                call relax(z(iu(m)), z(irhs(m)), z0(ido(m)), &
                     iq(m), jq(m), kq(m), m)
             end do
          end do
       end do
    end do
   
    call copy(u, z(iu(ng)), i0, j0, k0)

    return

  end subroutine mglin




  subroutine rstrct(uc, uf, ic0, jc0, kc0, nx, ny, nz, n)

    integer, intent(in) :: ic0
    integer, intent(in) :: jc0
    integer, intent(in) :: kc0
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    integer, intent(in) :: nz
    integer, intent(in) :: n
    real(8), intent(inout) :: uc(0:       ic0+1,0:       jc0+1,0:       kc0+1)
    real(8), intent(in)    :: uf(0:(nx+1)*ic0+1,0:(ny+1)*jc0+1,0:(nz+1)*kc0+1)

#ifdef MPI
    real(8), allocatable :: tmp(:,:,:)

    if(n == ngx) then
       allocate(tmp(0:iq(n)/iprc+1, 0:jq(n)/jprc+1, 0:kq(n)/kprc+1))
       call rstrct0(tmp, uf, iq(n)/iprc, jq(n)/jprc, kq(n)/kprc, nx, ny, nz)
       call gather(uc, tmp, iq(n), jq(n), kq(n))
       deallocate(tmp)
    else
       call rstrct0(uc, uf, ic0, jc0, kc0, nx, ny, nz)
    end if
#else
    call rstrct0(uc, uf, ic0, jc0, kc0, nx, ny, nz)
#endif

    return

  end subroutine rstrct


  subroutine rstrct0(uc, uf, ic0, jc0, kc0, nx, ny, nz)

    integer, intent(in) :: ic0
    integer, intent(in) :: jc0
    integer, intent(in) :: kc0
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    integer, intent(in) :: nz
    real(8), intent(inout) :: uc(0:       ic0+1,0:       jc0+1,0:       kc0+1)
    real(8), intent(in)    :: uf(0:(nx+1)*ic0+1,0:(ny+1)*jc0+1,0:(nz+1)*kc0+1)

    integer :: ic, if, jc, jf, kc, kf
    integer :: l, lm, im, imjm, ista, iend, jsta, jend, ksta, kend

    ista = 1; iend = ic0
    jsta = 1; jend = jc0
    ksta = 1; kend = kc0

    im =   (iend - ista + 1)
    imjm = (jend - jsta + 1) * im
    lm =   (kend - ksta + 1) * imjm

!CDIR LOOPCNT=LOOP_CNT
!OCL NOVREC
!OCL REPEAT(LOOP_CNT)
    do l = 0, lm - 1
       kc = (l / imjm) + ksta
       jc =((l - (kc - ksta) * imjm) / im) + jsta
       ic = (l - (kc - ksta) * imjm - (jc - jsta) * im) + ista
       if = (1 + nx) * ic - nx
       jf = (1 + ny) * jc - ny
       kf = (1 + nz) * kc - nz
       uc(ic,jc,kc) = .125d0 * (uf(if,jf   ,kf   ) + uf(if+nx,jf   ,kf   ) &
                              + uf(if,jf+ny,kf   ) + uf(if+nx,jf+ny,kf   ) &
                              + uf(if,jf   ,kf+nz) + uf(if+nx,jf   ,kf+nz) &
                              + uf(if,jf+ny,kf+nz) + uf(if+nx,jf+ny,kf+nz))
    end do

    return

  end subroutine rstrct0



  subroutine interp(uf, uc, if0, jf0, kf0, nx, ny, nz, n)
    integer, intent(in) :: if0
    integer, intent(in) :: jf0
    integer, intent(in) :: kf0
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    integer, intent(in) :: nz
    integer, intent(in) :: n
    real(8), intent(out) :: uf(0:if0       +1,0:jf0       +1,0:kf0       +1)
    real(8), intent(in)  :: uc(0:if0/(nx+1)+1,0:jf0/(ny+1)+1,0:kf0/(nz+1)+1)

#ifdef MPI
    real(8), allocatable :: tmp(:,:,:)

    if(n == ngx) then
       allocate(tmp(0:iq(n)/iprc+1,0:jq(n)/jprc+1,0:kq(n)/kprc+1))
       call scatter(tmp, uc, iq(n), jq(n), kq(n))
       call interp0(uf, tmp, iq(n+1), jq(n+1), kq(n+1), nx, ny, nz)
       deallocate(tmp)
    else
       call interp0(uf, uc, if0, jf0, kf0, nx, ny, nz)
    end if
#else
    call interp0(uf, uc, if0, jf0, kf0, nx, ny, nz)
#endif

    return
  end subroutine interp


  subroutine interp0(uf, uc, if0, jf0, kf0, nx, ny, nz)

    integer, intent(in) :: if0
    integer, intent(in) :: jf0
    integer, intent(in) :: kf0
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    integer, intent(in) :: nz
    real(8), intent(out) :: uf(0:if0       +1,0:jf0       +1,0:kf0       +1)
    real(8), intent(in)  :: uc(0:if0/(nx+1)+1,0:jf0/(ny+1)+1,0:kf0/(nz+1)+1)

    integer :: ic, if, jc, jf, kc, kf
    integer :: ic0, jc0, kc0
    integer :: l, lm, im, imjm, ista, iend, jsta, jend, ksta, kend

    ic0 = if0 / (nx + 1)
    jc0 = jf0 / (ny + 1)
    kc0 = kf0 / (nz + 1)

    ista = 1 - nx; iend = ic0
    jsta = 1 - ny; jend = jc0
    ksta = 1 - nz; kend = kc0

    im =   (iend - ista + 1)
    imjm = (jend - jsta + 1) * im
    lm =   (kend - ksta + 1) * imjm

!CDIR LOOPCNT=LOOP_CNT
!OCL NOVREC
!OCL REPEAT(LOOP_CNT)
    do l = 0, lm - 1
       kc = (l / imjm) + ksta
       jc =((l - (kc - ksta) * imjm) / im) + jsta
       ic = (l - (kc - ksta) * imjm - (jc - jsta) * im) + ista

       if = ic * (nx + 1)
       jf = jc * (ny + 1)
       kf = kc * (nz + 1)
       uf(if,jf,kf) = &
            (27d0 * uc(ic,jc   ,kc   ) + 9d0 * uc(ic+nx,jc   ,kc   )&
            + 9d0 * uc(ic,jc+ny,kc   ) + 3d0 * uc(ic+nx,jc+ny,kc   )&
            + 9d0 * uc(ic,jc   ,kc+nz) + 3d0 * uc(ic+nx,jc   ,kc+nz)&
            + 3d0 * uc(ic,jc+ny,kc+nz) + 1d0 * uc(ic+nx,jc+ny,kc+nz)&
            ) / 64d0

       if = ic * (nx + 1) + nx
       jf = jc * (ny + 1)
       kf = kc * (nz + 1)
       uf(if,jf,kf) = &
            ( 9d0 * uc(ic,jc   ,kc   ) +27d0 * uc(ic+nx,jc   ,kc   )&
            + 3d0 * uc(ic,jc+ny,kc   ) + 9d0 * uc(ic+nx,jc+ny,kc   )&
            + 3d0 * uc(ic,jc   ,kc+nz) + 9d0 * uc(ic+nx,jc   ,kc+nz)&
            + 1d0 * uc(ic,jc+ny,kc+nz) + 3d0 * uc(ic+nx,jc+ny,kc+nz)&
            ) / 64d0

       if = ic * (nx + 1)
       jf = jc * (ny + 1) + ny
       kf = kc * (nz + 1)
       uf(if,jf,kf) = &
            ( 9d0 * uc(ic,jc   ,kc   ) + 3d0 * uc(ic+nx,jc   ,kc   )&
            +27d0 * uc(ic,jc+ny,kc   ) + 9d0 * uc(ic+nx,jc+ny,kc   )&
            + 3d0 * uc(ic,jc   ,kc+nz) + 1d0 * uc(ic+nx,jc   ,kc+nz)&
            + 9d0 * uc(ic,jc+ny,kc+nz) + 3d0 * uc(ic+nx,jc+ny,kc+nz)&
            ) / 64d0

       if = ic * (nx + 1)
       jf = jc * (ny + 1)
       kf = kc * (nz + 1) + nz
       uf(if,jf,kf) = &
            ( 9d0 * uc(ic,jc   ,kc   ) + 3d0 * uc(ic+nx,jc   ,kc   )&
            + 3d0 * uc(ic,jc+ny,kc   ) + 1d0 * uc(ic+nx,jc+ny,kc   )&
            +27d0 * uc(ic,jc   ,kc+nz) + 9d0 * uc(ic+nx,jc   ,kc+nz)&
            + 9d0 * uc(ic,jc+ny,kc+nz) + 3d0 * uc(ic+nx,jc+ny,kc+nz)&
            ) / 64d0

       if = ic * (nx + 1)
       jf = jc * (ny + 1) + ny
       kf = kc * (nz + 1) + nz
       uf(if,jf,kf) =&
            ( 3d0 * uc(ic,jc   ,kc   ) + 1d0 * uc(ic+nx,jc   ,kc   )&
            + 9d0 * uc(ic,jc+ny,kc   ) + 3d0 * uc(ic+nx,jc+ny,kc   )&
            + 9d0 * uc(ic,jc   ,kc+nz) + 3d0 * uc(ic+nx,jc   ,kc+nz)&
            +27d0 * uc(ic,jc+ny,kc+nz) + 9d0 * uc(ic+nx,jc+ny,kc+nz)&
            ) / 64d0

       if = ic * (nx + 1) + nx
       jf = jc * (ny + 1)
       kf = kc * (nz + 1) + nz
       uf(if,jf,kf) =&
            ( 3d0 * uc(ic,jc   ,kc   ) + 9d0 * uc(ic+nx,jc   ,kc   )&
            + 1d0 * uc(ic,jc+ny,kc   ) + 3d0 * uc(ic+nx,jc+ny,kc   )&
            + 9d0 * uc(ic,jc   ,kc+nz) +27d0 * uc(ic+nx,jc   ,kc+nz)&
            + 3d0 * uc(ic,jc+ny,kc+nz) + 9d0 * uc(ic+nx,jc+ny,kc+nz)&
            ) / 64d0

       if = ic * (nx + 1) + nx
       jf = jc * (ny + 1) + ny
       kf = kc * (nz + 1)
       uf(if,jf,kf) =&
            ( 3d0 * uc(ic,jc   ,kc   ) + 9d0 * uc(ic+nx,jc   ,kc   )&
            + 9d0 * uc(ic,jc+ny,kc   ) +27d0 * uc(ic+nx,jc+ny,kc   )&
            + 1d0 * uc(ic,jc   ,kc+nz) + 3d0 * uc(ic+nx,jc   ,kc+nz)&
            + 3d0 * uc(ic,jc+ny,kc+nz) + 9d0 * uc(ic+nx,jc+ny,kc+nz)&
            ) / 64d0

       if = ic * (nx + 1) + nx
       jf = jc * (ny + 1) + ny
       kf = kc * (nz + 1) + nz
       uf(if,jf,kf) =&
            ( 1d0 * uc(ic,jc   ,kc   ) + 3d0 * uc(ic+nx,jc   ,kc   )&
            + 3d0 * uc(ic,jc+ny,kc   ) + 9d0 * uc(ic+nx,jc+ny,kc   )&
            + 3d0 * uc(ic,jc   ,kc+nz) + 9d0 * uc(ic+nx,jc   ,kc+nz)&
            + 9d0 * uc(ic,jc+ny,kc+nz) +27d0 * uc(ic+nx,jc+ny,kc+nz)&
            ) / 64d0
    end do

    return
    
  end subroutine interp0




  subroutine relax(u, rhs, d, i0, j0, k0, n)

    real(8), intent(inout) ::   u(0:i0+1,0:j0+1,0:k0+1)
    real(8), intent(in)    :: rhs(0:i0+1,0:j0+1,0:k0+1)
    real(8), intent(in)    ::   d(0:i0+1,0:j0+1,0:k0+1)
    integer, intent(in) :: i0
    integer, intent(in) :: j0
    integer, intent(in) :: k0
    integer, intent(in) :: n

    integer :: i, j, k
    integer :: ipass
    real(8) :: dim, dip, djm, djp, dkm, dkp
    real(8) :: dx2i, dy2i, dz2i
    real(8) :: uuf
    integer :: l, lm, im, imjm, ista, iend, jsta, jend, ksta, kend

    ista = 1; iend = i0
    jsta = 1; jend = j0
    ksta = 1; kend = k0

    im =   (iend - ista + 1)
    imjm = (jend - jsta + 1) * im
    lm =   (kend - ksta + 1) * imjm

    dx2i = 1d0 / (lx / dble(i0)) ** 2
    dy2i = 1d0 / (ly / dble(j0)) ** 2
    dz2i = 1d0 / (lz / dble(k0)) ** 2
#ifdef MPI
    if(n <= ngx) then
       dx2i = 1d0 / (lx / (dble(i0) / dble(iprc))) ** 2 
       dy2i = 1d0 / (ly / (dble(j0) / dble(jprc))) ** 2
       dz2i = 1d0 / (lz / (dble(k0) / dble(kprc))) ** 2
    end if
#endif

    do ipass = 1, 2
!CDIR NODEP
!CDIR LOOPCNT=LOOP_CNT
!OCL NOVREC
!OCL REPEAT(LOOP_CNT)
       do l = 0, lm - 1
          k = (l / imjm) + ksta
          j =((l - (k - ksta) * imjm) / im) + jsta
          i = (l - (k - ksta) * imjm - (j - jsta) * im) + ista
          if((ipass - mod(i + j + k + 1, 2)) == 1) then
          dim = .5d0 * (d(i-1,j,k) + d(i,j,k))
          dip = .5d0 * (d(i+1,j,k) + d(i,j,k))
          djm = .5d0 * (d(i,j-1,k) + d(i,j,k))
          djp = .5d0 * (d(i,j+1,k) + d(i,j,k))
          dkm = .5d0 * (d(i,j,k-1) + d(i,j,k))
          dkp = .5d0 * (d(i,j,k+1) + d(i,j,k))
          uuf = (dip + dim) * dx2i +(djp + djm) * dy2i +(dkp + dkm) * dz2i &
               + 1d0 / dt
          u(i,j,k) = ((dip * u(i+1,j,k) + dim * u(i-1,j,k)) * dx2i &
                     +(djp * u(i,j+1,k) + djm * u(i,j-1,k)) * dy2i &
                     +(dkp * u(i,j,k+1) + dkm * u(i,j,k-1)) * dz2i &
                     -rhs(i,j,k)) / uuf
          end if
       end do
       call bndcz(u, i0, j0, k0, n, 0)
       call bndcx(u, i0, j0, k0, n)
       call bndcy(u, i0, j0, k0, n)
    end do
    
    return
  end subroutine relax



     
  subroutine resid(res, u, rhs, d, i0, j0, k0, n)

    integer, intent(in) :: i0
    integer, intent(in) :: j0
    integer, intent(in) :: k0
    real(8), intent(out) :: res(0:i0+1,0:j0+1,0:k0+1)
    real(8), intent(in)  ::   u(0:i0+1,0:j0+1,0:k0+1)
    real(8), intent(in)  :: rhs(0:i0+1,0:j0+1,0:k0+1)
    real(8), intent(in)  ::   d(0:i0+1,0:j0+1,0:k0+1)
    integer, intent(in) :: n

    integer :: i, j, k
    real(8) :: dim, dip, djm, djp, dkm, dkp
    real(8) :: dx2i, dy2i, dz2i
    real(8) :: uu, uuf
    integer :: l, lm, im, imjm, ista, iend, jsta, jend, ksta, kend

    ista = 1; iend = i0
    jsta = 1; jend = j0
    ksta = 1; kend = k0

    im =   (iend - ista + 1)
    imjm = (jend - jsta + 1) * im
    lm =   (kend - ksta + 1) * imjm

    dx2i = 1d0 / (lx / dble(i0)) ** 2
    dy2i = 1d0 / (ly / dble(j0)) ** 2
    dz2i = 1d0 / (lz / dble(k0)) ** 2
#ifdef MPI
    if(n <= ngx) then
       dx2i = 1d0 / (lx / (dble(i0) / dble(iprc))) ** 2 
       dy2i = 1d0 / (ly / (dble(j0) / dble(jprc))) ** 2
       dz2i = 1d0 / (lz / (dble(k0) / dble(kprc))) ** 2
    end if
#endif

!CDIR LOOPCNT=LOOP_CNT
!OCL NOVREC
!OCL REPEAT(LOOP_CNT)
    do l = 0, lm - 1
       k = (l / imjm) + ksta
       j =((l - (k - ksta) * imjm) / im) + jsta
       i = (l - (k - ksta) * imjm - (j - jsta) * im) + ista
       dim = .5d0 * (d(i-1,j,k) + d(i,j,k))
       dip = .5d0 * (d(i+1,j,k) + d(i,j,k))
       djm = .5d0 * (d(i,j-1,k) + d(i,j,k))
       djp = .5d0 * (d(i,j+1,k) + d(i,j,k))
       dkm = .5d0 * (d(i,j,k-1) + d(i,j,k))
       dkp = .5d0 * (d(i,j,k+1) + d(i,j,k))
       uuf = (dip + dim) * dx2i + (djp + djm) * dy2i + (dkp + dkm) * dz2i &
            + 1d0 / dt
       uu = ((dip * u(i+1,j,k) + dim * u(i-1,j,k)) * dx2i &
            +(djp * u(i,j+1,k) + djm * u(i,j-1,k)) * dy2i &
            +(dkp * u(i,j,k+1) + dkm * u(i,j,k-1)) * dz2i - rhs(i,j,k)) / uuf
       res(i,j,k) = (u(i,j,k) - uu) * uuf
    end do

    return

  end subroutine resid




  subroutine copy(aout, ain, i0, j0, k0)

    integer, intent(in) :: i0
    integer, intent(in) :: j0
    integer, intent(in) :: k0
    real(8), intent(out) :: aout(0:i0+1,0:j0+1,0:k0+1)
    real(8), intent(in)  ::  ain(0:i0+1,0:j0+1,0:k0+1)

    integer :: i, j, k
    integer :: l, lm, im, imjm, ista, iend, jsta, jend, ksta, kend

    ista = 0; iend = i0 + 1
    jsta = 0; jend = j0 + 1
    ksta = 0; kend = k0 + 1

    im =   (iend - ista + 1)
    imjm = (jend - jsta + 1) * im
    lm =   (kend - ksta + 1) * imjm

!CDIR LOOPCNT=LOOP_CNT
!OCL NOVREC
!OCL REPEAT(LOOP_CNT)
    do l = 0, lm - 1
       k = (l / imjm) + ksta
       j =((l - (k - ksta) * imjm) / im) + jsta
       i = (l - (k - ksta) * imjm - (j - jsta) * im) + ista
       aout(i,j,k) = ain(i,j,k)
    end do

    return

  end subroutine copy




  subroutine fill0(u, i0, j0, k0)

    integer, intent(in) :: i0
    integer, intent(in) :: j0
    integer, intent(in) :: k0
    real(8), intent(out) :: u(0:i0+1,0:j0+1,0:k0+1)

    integer :: i, j, k
    integer :: l, lm, im, imjm, ista, iend, jsta, jend, ksta, kend

    ista = 0; iend = i0 + 1
    jsta = 0; jend = j0 + 1
    ksta = 0; kend = k0 + 1

    im =   (iend - ista + 1)
    imjm = (jend - jsta + 1) * im
    lm =   (kend - ksta + 1) * imjm

!CDIR LOOPCNT=LOOP_CNT
!OCL NOVREC
!OCL REPEAT(LOOP_CNT)
    do l = 0, lm - 1
       k = (l / imjm) + ksta
       j =((l - (k - ksta) * imjm) / im) + jsta
       i = (l - (k - ksta) * imjm - (j - jsta) * im) + ista
       u(i,j,k) = 0d0
    end do

    return

  end subroutine fill0




  subroutine slvsml(u, rhs, d, i0, j0, k0)

    integer, intent(in) :: i0
    integer, intent(in) :: j0
    integer, intent(in) :: k0
    real(8), intent(inout) ::   u(0:i0+1,0:j0+1,0:k0+1)
    real(8), intent(in)    :: rhs(0:i0+1,0:j0+1,0:k0+1)
    real(8), intent(in)    ::   d(0:i0+1,0:j0+1,0:k0+1)

    real(8) :: res(0:i0+1,0:j0+1,0:k0+1), anorm
    integer :: iter
    real(8) :: dz2

#ifdef MPI
    if(myrk /= 0) then
       return
    endif
#endif

    if((i0 == 1) .and. (j0 == 1) .and. (k0 ==1)) then
       ! exact solution
       dz2 = (lz) ** 2
#ifdef MPI
       dz2 = (lz / (1d0 / dble(kprc))) ** 2
#endif
       u(i0,j0,k0) = rhs(i0,j0,k0) &
            / (((d(i0,j0,k0+1) + d(i0,j0,k0)) * (z1(ibo(1)) - 1d0) &
              + (d(i0,j0,k0-1) + d(i0,j0,k0)) * (z1(ibi(1)) - 1d0)) &
              / (2d0 * dz2) - (1d0 / dt))
       u(  i0  ,  j0  ,  k0+1) = u(  i0  ,j0,  k0  ) * z1(ibo(1)) ! bndcz
       u(  i0  ,  j0  ,0     ) = u(  i0  ,j0,  k0  ) * z1(ibi(1)) ! bndcz
       u(  i0+1,  j0  ,0:k0+1) = u(  i0  ,j0,0:k0+1) ! bndcx
       u(0     ,  j0  ,0:k0+1) = u(  i0  ,j0,0:k0+1) ! bndcx
       u(0:i0+1,  j0+1,0:k0+1) = u(0:i0+1,j0,0:k0+1) ! bndcy
       u(0:i0+1,0     ,0:k0+1) = u(0:i0+1,j0,0:k0+1) ! bndcy
       return
    else
       ! approximate solution
       call fill0(u, i0, j0, k0)

       do iter = 1, ITERMAX_SLVSML
          call relax(u, rhs, d, i0, j0, k0, 1)
          call resid(res, u, rhs, d, i0, j0, k0, 1)
          call norm(anorm, res, rhs, i0, j0, k0, 1)
          if(anorm < EPS_SLVSML) then
             return ! normal exit
          end if
       end do

       print *, 'slvsml: did not converge: anorm=', anorm
#ifdef MPI
       call mpi_abort(ier)
#endif
       stop
    endif

    return

  end subroutine slvsml



      
  integer function maloc(len)

    integer, intent(in) :: len

    if(mem + len + 1 > memlen) then
       print *, 'insufficient memory in maloc', mem + len + 1, memlen
       stop
    end if
    z(mem + 1) = len
    maloc = mem + 2
    mem = mem + len + 1
    
    return

  end function maloc



  
  integer function maloc0(len)

    integer, intent(in) :: len

    if (mem0 + len + 1 > memlen0) then
       print *, 'insufficient memory in maloc0'
       stop
    end if
    z0(mem0 + 1) = len
    maloc0 = mem0 + 2
    mem0 = mem0 + len + 1
    
    return

  end function maloc0




  integer function maloc1(len)

    integer, intent(in) :: len
    
    if(mem1 + len + 1 > memlen1) then
       print *, 'insufficient memory in maloc1'
       stop
    end if
    z1(mem1 + 1) = len
    maloc1 = mem1 + 2
    mem1 = mem1 + len + 1
    
    return

  end function maloc1




  subroutine bndcx(u, i0, j0, k0, n)

    use bval_module, only : deltay

    real(8), intent(inout) :: u(0:i0+1,0:j0+1,0:k0+1)
    integer, intent(in) :: i0
    integer, intent(in) :: j0
    integer, intent(in) :: k0
    integer, intent(in) :: n

    integer :: j, k
    integer :: jremap, jplus
    real(8) :: eps, dq, delq2, epsi, epso, dy0
    real(8), allocatable :: q(:,:), f(:), delq(:)
    integer :: j00, jj, jshift
#ifdef MPI
    real(8), allocatable :: ui(:,:,:), uo(:,:,:)
#endif

    dy0 = ly / dble(j0)
#ifdef MPI
    if(n <= ngx) then
       dy0 = ly / (dble((j0)) / dble(jprc))
    end if
#endif

    jplus = int(deltay / dy0)
    epsi = dmod(deltay, dy0)
    epso = -epsi

    j00 = jq(n)
    jshift = 0
#ifdef MPI
    if(n > ngx) then
       j00 = j00 * jprc
       jshift = myrkj * j0
    end if
#endif

    allocate(q(-1:j00+2,0:k0+1), f(1:j00+1), delq( 0:j00+2))

#ifdef MPI
    allocate(ui(1:1,1:j00,0:k0+1), uo(i0:i0,1:j00,0:k0+1))
    call gatherx2(u, ui, uo, i0, j0, k0, 1, 1, i0, i0, 0, k0+1, j00)

    if(myrki == 0) then
#endif
       eps = epsi / dy0
       do k = 0, k0 + 1
!CDIR SHORT_LOOP
       do j = 1, j00
          jremap = j - jplus
          if(jremap <= 0) jremap = jremap + j00
#ifdef MPI       
          q(j,k) =uo(i0,jremap,k)
#else
          q(j,k) = u(i0,jremap,k)
#endif
       end do
       end do
       q( 0   ,0:k0+1) = q(j00  ,0:k0+1)
       q(-1   ,0:k0+1) = q(j00-1,0:k0+1)
       q(j00+1,0:k0+1) = q(    1,0:k0+1)
       do k = 0, k0 + 1
!CDIR SHORT_LOOP
          do j = 0, j0 + 1
             jj = j + jshift
             delq(jj) = q(jj,k) - q(jj - 1,k)
          end do
!CDIR SHORT_LOOP
          do j = 1, j0 + 1
             jj = j + jshift
             delq2 = delq(jj) * delq(jj - 1)
             dq = 0
             if(delq2 > 0) dq = delq2 / (delq(jj - 1) + delq(jj))
             f(jj) = eps * (q(jj - 1,k) + (1d0 - eps) * dq)
          end do
!CDIR SHORT_LOOP
          do j = 1, j0
             jj = j + jshift
             u(0,j,k) = q(jj,k) - (f(jj + 1) - f(jj))
          end do
       end do
#ifdef MPI
    end if
    if(myrki == iprc-1) then
#endif
       eps = epso / dy0
       do k = 0, k0 + 1
!CDIR SHORT_LOOP
       do j = 1, j00
          jremap = j + jplus
          if(jremap >= j00 + 1) jremap = jremap - j00
#ifdef MPI
          q(j,k) =ui(1,jremap,k)
#else
          q(j,k) = u(1,jremap,k)
#endif
       end do
       end do
       q(    0,0:k0+1) = q(j00,0:k0+1)
       q(j00+1,0:k0+1) = q(  1,0:k0+1)
       q(j00+2,0:k0+1) = q(  2,0:k0+1)
       do k = 0, k0 + 1
!CDIR SHORT_LOOP
          do j = 1, j0 + 2
             jj = j + jshift
             delq(jj) = q(jj,k) - q(jj - 1,k)
          end do
!CDIR SHORT_LOOP
          do j = 1, j0 + 1
             jj = j + jshift
             delq2 = delq(jj + 1) * delq(jj)
             dq = 0
             if(delq2 > 0) dq = delq2 / (delq(jj) + delq(jj + 1))
             f(jj) = eps * (q(jj,k) - (1d0 + eps) * dq)
          end do
!CDIR SHORT_LOOP
          do j = 1, j0
             jj = j + jshift
             u(i0+1,j,k) = q(jj,k) - (f(jj + 1) - f(jj))
          end do
       end do
#ifdef MPI
    end if

    deallocate(ui, uo)
#endif

    deallocate(q, f, delq)

    return

  end subroutine bndcx




  subroutine bndcy(u, i0, j0, k0, n)

    real(8), intent(inout) :: u(0:i0+1,0:j0+1,0:k0+1)
    integer, intent(in) :: i0
    integer, intent(in) :: j0
    integer, intent(in) :: k0
    integer, intent(in) :: n

    u(0:i0+1,   0,0:k0+1) = u(0:i0+1,j0,0:k0+1)
    u(0:i0+1,j0+1,0:k0+1) = u(0:i0+1, 1,0:k0+1)

#ifdef MPI
    if(n <= ngx) return
    call shiftx2(u, i0, j0, k0, 0, i0+1, j0+1, j0+1, 0, 0, 0, k0+1)
#endif

    return

  end subroutine bndcy




  subroutine bndcz(u, i0, j0, k0, n, id)

    real(8), intent(inout) :: u(0:i0+1,0:j0+1,0:k0+1)
    integer, intent(in) :: i0
    integer, intent(in) :: j0
    integer, intent(in) :: k0
    integer, intent(in) :: n
    integer, intent(in) :: id

    integer :: i, j
    
    if(id == 0) then
       do i = 1, i0
!CDIR SHORT_LOOP
       do j = 1, j0
          u(i,j,   0) = u(i,j, 1) * z1(ibi(n) + (i-1) + (j-1)*(i0))
          u(i,j,k0+1) = u(i,j,k0) * z1(ibo(n) + (i-1) + (j-1)*(i0))
       end do
       end do
    else
       do i = 1, i0
!CDIR SHORT_LOOP
       do j = 1, j0
          u(i,j,   0) = u(i,j,   2)
          u(i,j,k0+1) = u(i,j,k0-1)
       end do
       end do
    end if

#ifdef MPI
    if(n <= ngx) return
    call shiftx3(u,i0,j0,k0,1,i0,1,j0,k0+1,k0+1,0,0)
#endif

    return

  end subroutine bndcz




  integer function log20(n)

    integer, intent(in) :: n

    log20 = 0

    if(n ==   1) then
       log20 = 0
    else if(n ==   2) then
       log20 = 1
    else if(n ==   4) then
       log20 = 2
    else if(n ==   8) then
       log20 = 3
    else if(n ==  16) then
       log20 = 4
    else if(n ==  32) then
       log20 = 5
    else if(n ==  64) then
       log20 = 6
    else if(n == 128) then
       log20 = 7
    else if(n == 256) then
       log20 = 8
    else if(n == 512) then
       log20 = 9
    else
       print *,"error: log20, n = ", n
#ifdef MPi
       call mpi_finalize(ier)
#endif
       stop
    end if

    return

  end function log20




  subroutine coarsekb

    use field_module, only : erikb, erokb
    use grid_module, only : is, ie, js, je

    integer :: i, j, i0, j0
    integer n
    real(8) :: bko(ie-is+1,je-js+1), bki(ie-is+1,je-js+1)

    i0 = ie - is + 1
    j0 = je - js + 1

    do i = 1, i0
!CDIR SHORT_LOOP
    do j = 1, j0
       bko(i,j) = erokb(i - 1 + is, j - 1 + js)
       bki(i,j) = erikb(i - 1 + is, j - 1 + js)
    end do
    end do

    n = ng
    call copy2(z1(ibo(n)), bko, iq(n), jq(n))
    call copy2(z1(ibi(n)), bki, iq(n), jq(n))
    
    do n = ng - 1, 1, -1
#ifdef MPI
       call rstrct2(z1(ibo(n)), z1(ibo(n+1)), iq(n), jq(n), nx(n), ny(n), n, &
            top_root)
       call rstrct2(z1(ibi(n)), z1(ibi(n+1)), iq(n), jq(n), nx(n), ny(n), n, &
            bottom_root)
#else
       call rstrct20(z1(ibo(n)), z1(ibo(n+1)), iq(n), jq(n), nx(n), ny(n))
       call rstrct20(z1(ibi(n)), z1(ibi(n+1)), iq(n), jq(n), nx(n), ny(n))
#endif
    end do
    
    return

  end subroutine coarsekb




#ifdef MPI
  subroutine rstrct2(uc, uf, ic0, jc0, nx, ny, n, root)

    integer, intent(in) :: ic0
    integer, intent(in) :: jc0
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    integer, intent(in) :: n
    integer, intent(in) :: root
    real(8), intent(out) :: uc(1:       ic0,1:       jc0)
    real(8), intent(in)  :: uf(1:(nx+1)*ic0,1:(ny+1)*jc0)

    real(8), allocatable :: tmp(:,:)
    
    if(n == ngx) then
       allocate(tmp(iq(n)/iprc,jq(n)/jprc))
       call rstrct20(tmp, uf, iq(n)/iprc, jq(n)/jprc, nx, ny)
       call gather2(uc, tmp, iq(n), jq(n), root)
       deallocate(tmp)
    else
       call rstrct20(uc, uf, ic0, jc0, nx, ny)
    end if

    return

  end subroutine rstrct2
#endif


  subroutine rstrct20(uc, uf, ic0, jc0, nx, ny)

    integer, intent(in) :: ic0
    integer, intent(in) :: jc0
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    real(8), intent(out) :: uc(1:       ic0,1:       jc0)
    real(8), intent(in)  :: uf(1:(nx+1)*ic0,1:(ny+1)*jc0)

    integer :: ic, if, jc, jf

    do ic = 1, ic0
!CDIR SHORT_LOOP
    do jc = 1, jc0
       if = (1 + nx) * ic - nx
       jf = (1 + ny) * jc - ny
       uc(ic,jc) = .25d0 * (uf(if  ,jf   ) + uf(if+nx,jf   ) &
                          + uf(if  ,jf+ny) + uf(if+nx,jf+ny))
    end do
    end do

    return

  end subroutine rstrct20




  subroutine copy2(aout, ain, i0, j0)

    real(8), intent(out) :: aout(i0,j0)
    real(8), intent(in)  ::  ain(i0,j0)
    integer, intent(in) :: i0
    integer, intent(in) :: j0

    aout(1:i0, 1:j0) = ain(1:i0, 1:j0)

    return

  end subroutine copy2




#ifdef MPI
  subroutine scatter(u0, u, i0, j0, k0)

    use grid_module, only : para_range

    integer, intent(in) :: i0
    integer, intent(in) :: j0
    integer, intent(in) :: k0
    real(8), intent(out) :: u0(0:i0/iprc+1,0:j0/jprc+1,0:k0/kprc+1)
    real(8), intent(in) :: u(0:i0+1,0:j0+1,0:k0+1)

    real(8), allocatable :: tmp(:), tmp0(:)
    integer :: i, j, k, l
    integer :: is2, ie2, js2, je2, ks2, ke2
    integer :: io, jo, ko
    integer :: irk, jrk, krk

    io = i0 / iprc
    jo = j0 / jprc
    ko = k0 / kprc
    allocate(tmp0(0:((io + 2) * (jo + 2) * (ko + 2) - 1)))
    allocate( tmp(0:((io + 2) * (jo + 2) * (ko + 2) * nprc - 1)))

    if(myrk == 0) then
       l = 0
       irk = 0
       do jrk = 0, jprc - 1
       do krk = 0, kprc - 1
          call para_range(1, i0, iprc, irk, is2, ie2)
          call para_range(1, j0, jprc, jrk, js2, je2)
          call para_range(1, k0, kprc, krk, ks2, ke2)
          do k = ks2 - 1, ke2 + 1
          do j = js2 - 1, je2 + 1
          do i = is2 - 1, ie2 + 1
             tmp(l) = u(i,j,k)
             l = l + 1
          end do
          end do
          end do
       end do
       end do
    end if

    call mpi_scatter(tmp,  ((io + 2) * (jo + 2) * (ko + 2)), &
                     mpi_double_precision, &
                     tmp0, ((io + 2) * (jo + 2) * (ko + 2)), &
                     mpi_double_precision, 0, mpi_comm_world, ier)

    l = 0
    do k = 0, ko + 1
    do j = 0, jo + 1
    do i = 0, io + 1
       u0(i,j,k) = tmp0(l)
       l = l + 1
    end do
    end do
    end do

    deallocate(tmp0, tmp)

    return

  end subroutine scatter
  



  subroutine gather(u, u0, i0, j0, k0)

    use grid_module, only : para_range

    integer, intent(in) :: i0
    integer, intent(in) :: j0
    integer, intent(in) :: k0
    real(8), intent(in) :: u0(0:i0/iprc+1,0:j0/jprc+1,0:k0/kprc+1)
    real(8), intent(out) :: u(0:i0+1,0:j0+1,0:k0+1)

    real(8), allocatable :: tmp(:), tmp0(:)
    integer :: i, j, k, l
    integer :: is2, ie2, js2, je2, ks2, ke2
    integer :: io, jo, ko
    integer :: irk, jrk, krk

    io = i0 / iprc
    jo = j0 / jprc
    ko = k0 / kprc
    allocate(tmp0(0:(io * jo * ko - 1)))
    allocate( tmp(0:(io * jo * ko * nprc - 1)))

    l = 0
    do k = 1, ko
    do j = 1, jo
    do i = 1, io
       tmp0(l) = u0(i,j,k)
       l = l + 1
    end do
    end do
    end do

    call mpi_gather(tmp0, (io * jo * ko), mpi_double_precision, &
                    tmp,  (io * jo * ko), mpi_double_precision, &
                    0, mpi_comm_world, ier)

    if(myrk == 0) then
       l = 0
       irk = 0
       do jrk = 0, jprc - 1
       do krk = 0, kprc - 1
          call para_range(1, i0, iprc, irk, is2, ie2)
          call para_range(1, j0, jprc, jrk, js2, je2)
          call para_range(1, k0, kprc, krk, ks2, ke2)
          do k = ks2, ke2
          do j = js2, je2
          do i = is2, ie2
             u(i,j,k) = tmp(l)
             l = l + 1
          end do
          end do
          end do
       end do
       end do
    end if

    deallocate(tmp0, tmp)

    return

  end subroutine gather




  subroutine gather2(u, u0, i0, j0, iroot)

    use grid_module, only : para_range

    integer, intent(in) :: i0
    integer, intent(in) :: j0
    real(8), intent(in) :: u0(1:i0/iprc, 1:j0/jprc)
    real(8), intent(out) :: u(1:i0, 1:j0)
    integer, intent(in) :: iroot

    real(8), allocatable :: tmp(:), tmp0(:)
    integer :: i, j, l
    integer :: is2, ie2, js2, je2
    integer :: io, jo
    integer :: irk, jrk

    io = i0 / iprc
    jo = j0 / jprc
    allocate(tmp0(0:(io * jo - 1)))
    allocate( tmp(0:(io * jo * (iprc * jprc) - 1)))

    l = 0
    do j = 1, jo
    do i = 1, io
       tmp0(l) = u0(i,j)
       l = l + 1
    end do
    end do

    call mpi_gather(tmp0, (io * jo), mpi_double_precision, &
                    tmp,  (io * jo), mpi_double_precision, &
                    0, mpi_k_world, ier)

    if(myrk_k == 0) then
       l = 0
       irk = 0
       do jrk = 0, jprc - 1
          call para_range(1, i0, iprc, irk, is2, ie2)
          call para_range(1, j0, jprc, jrk, js2, je2)
          do j = js2, je2
          do i = is2, ie2
             u(i,j) = tmp(l)
             l = l + 1
          end do
          end do
       end do
    end if

    deallocate(tmp0, tmp)

    call mpi_bcast(u, (i0 * j0), mpi_double_precision, &
         iroot, mpi_comm_world, ier)

    return
    
  end subroutine gather2




  subroutine shiftx3(w, i0, j0, k0, ib, if, jb, jf, kb1, kf1, kb2, kf2)

    integer, intent(in) :: i0
    integer, intent(in) :: j0
    integer, intent(in) :: k0
    real(8), intent(inout) :: w(0:i0+1,0:j0+1,0:k0+1)
    integer, intent(in) :: ib
    integer, intent(in) :: if
    integer, intent(in) :: jb
    integer, intent(in) :: jf
    integer, intent(in) :: kb1
    integer, intent(in) :: kf1
    integer, intent(in) :: kb2
    integer, intent(in) :: kf2

    integer :: jl1, jl2, kds
    integer :: isd1, isd2, irv1, irv2
    integer :: i, j, k, l
    integer ist(mpi_status_size)
    real(8) :: wps(0:(if - ib + 1) * (jf - jb + 1) * (kf2 - kb2 + 1) - 1)
    real(8) :: wms(0:(if - ib + 1) * (jf - jb + 1) * (kf1 - kb1 + 1) - 1)
    real(8) :: wpr(0:(if - ib + 1) * (jf - jb + 1) * (kf1 - kb1 + 1) - 1)
    real(8) :: wmr(0:(if - ib + 1) * (jf - jb + 1) * (kf2 - kb2 + 1) - 1)
    integer :: im, imjm
    
    kds = k0

    im =   (if - ib + 1)
    imjm = (jf - jb + 1) * im
    jl1 = imjm * (kf1 - kb1 + 1)
    jl2 = imjm * (kf2 - kb2 + 1)

    ! send to k+ (1)
    if(myrkk /= kprc-1) then
       l = 0
       do k = kb2 + kds, kf2 + kds
       do j = jb, jf
!CDIR SHORT_LOOP
       do i = ib, if
          wps(l) = w(i,j,k)
          l = l + 1
       end do
       end do
       end do
    end if
    ! send to k- (2)
    if(myrkk /= 0) then
       l = 0
       do k = kb1 - kds, kf1 - kds
       do j = jb, jf
!CDIR SHORT_LOOP
       do i = ib, if
          wms(l) = w(i,j,k)
          l = l + 1
       end do
       end do
       end do
    end if
    call mpi_isend(wps,jl2,mpi_double_precision,ku,1,mpi_comm_world,isd2,ier)
    call mpi_irecv(wmr,jl2,mpi_double_precision,kd,1,mpi_comm_world,irv2,ier)
    call mpi_wait(isd2,ist,ier)
    call mpi_wait(irv2,ist,ier)
    call mpi_isend(wms,jl1,mpi_double_precision,kd,1,mpi_comm_world,isd1,ier)
    call mpi_irecv(wpr,jl1,mpi_double_precision,ku,1,mpi_comm_world,irv1,ier)
    call mpi_wait(isd1,ist,ier)
    call mpi_wait(irv1,ist,ier)
    ! receive from k+ (1)
    if(myrkk /= kprc-1) then
       l = 0
       do k = kb1, kf1
       do j = jb, jf
!CDIR SHORT_LOOP
       do i = ib, if
          w(i,j,k) = wpr(l)
          l = l + 1
       end do
       end do
       end do
    end if
    ! receive from k- (2)
    if(myrkk /= 0) then
       l = 0
       do k = kb2, kf2
       do j = jb, jf
!CDIR SHORT_LOOP
       do i = ib, if
          w(i,j,k) = wmr(l)
          l = l + 1
       end do
       end do
       end do
    end if
    
    return

  end subroutine shiftx3
#endif




#ifdef MPI
  subroutine shiftx2(w, i0, j0, k0, ib, if, jb1, jf1, jb2, jf2, kb, kf)

    integer, intent(in) :: i0
    integer, intent(in) :: j0
    integer, intent(in) :: k0
    real(8), intent(inout) :: w(0:i0+1,0:j0+1,0:k0+1)
    integer, intent(in) :: ib
    integer, intent(in) :: if
    integer, intent(in) :: jb1
    integer, intent(in) :: jf1
    integer, intent(in) :: jb2
    integer, intent(in) :: jf2
    integer, intent(in) :: kb
    integer, intent(in) :: kf

    integer :: jl1, jl2, jds
    integer :: isd1, isd2, irv1, irv2
    integer :: i, j, k, l
    integer ist(mpi_status_size)
    real(8) :: wps(0:(if - ib + 1) * (jf2 - jb2 + 1) * (kf - kb + 1) - 1)
    real(8) :: wmr(0:(if - ib + 1) * (jf2 - jb2 + 1) * (kf - kb + 1) - 1)
    real(8) :: wpr(0:(if - ib + 1) * (jf1 - jb1 + 1) * (kf - kb + 1) - 1)
    real(8) :: wms(0:(if - ib + 1) * (jf1 - jb1 + 1) * (kf - kb + 1) - 1)
    integer :: im, imjm
    
    jds = j0
    im =   (if - ib + 1)
    imjm = (kf - kb + 1) * im
    jl1 = imjm * (jf1 - jb1 + 1)
    jl2 = imjm * (jf2 - jb2 + 1)

    ! send to j+ (1)
    l = 0
    do k = kb, kf
    do j = jb2 + jds, jf2 + jds
!CDIR SHORT_LOOP
    do i = ib, if
       wps(l) = w(i,j,k)
       l = l + 1
    end do
    end do
    end do
    ! send to j- (2)
    l = 0
    do k = kb, kf
    do j = jb1 - jds, jf1 - jds
!CDIR SHORT_LOOP
    do i = ib, if
       wms(l) = w(i,j,k)
       l = l + 1
    end do
    end do
    end do
    if(jprc /= 1) then
       call mpi_isend(wps,jl2,mpi_double_precision,ju,1,mpi_comm_world,isd2,ier)
       call mpi_irecv(wmr,jl2,mpi_double_precision,jd,1,mpi_comm_world,irv2,ier)
       call mpi_wait(isd2,ist,ier)
       call mpi_wait(irv2,ist,ier)
       call mpi_isend(wms,jl1,mpi_double_precision,jd,1,mpi_comm_world,isd1,ier)
       call mpi_irecv(wpr,jl1,mpi_double_precision,ju,1,mpi_comm_world,irv1,ier)
       call mpi_wait(isd1,ist,ier)
       call mpi_wait(irv1,ist,ier)
    else
       wpr(0:jl1-1) = wms(0:jl1-1)
       wmr(0:jl2-1) = wps(0:jl2-1)
    end if
    ! receive from j+ (1)
    l = 0 
    do k = kb, kf
    do j = jb1, jf1
!CDIR SHORT_LOOP
    do i = ib, if
       w(i,j,k) = wpr(l)
       l = l + 1
    end do
    end do
    end do
    ! receive from j- (2)
    l = 0
    do k = kb, kf
    do j = jb2, jf2
!CDIR SHORT_LOOP
    do i = ib, if
       w(i,j,k) = wmr(l)
       l = l + 1
    end do
    end do
    end do
    
    return

  end subroutine shiftx2




  subroutine gatherx2(w, wi, wo, i0, j0, k0, ibi, ifi, ibo, ifo, kb, kf, j00)

    real(8), intent(in) :: w(0:i0+1,0:j0+1,0:k0+1)
    real(8), intent(out) :: wi(ibi:ifi, 1:j00, kb:kf)
    real(8), intent(out) :: wo(ibo:ifo, 1:j00, kb:kf)
    integer, intent(in) :: i0
    integer, intent(in) :: j0
    integer, intent(in) :: k0
    integer, intent(in) :: ibi
    integer, intent(in) :: ifi
    integer, intent(in) :: ibo
    integer, intent(in) :: ifo
    integer, intent(in) :: kb
    integer, intent(in) :: kf
    integer, intent(in) :: j00

    integer :: l
    integer :: li, li0
    integer :: lo, lo0
    integer :: i, j, k
    real(8), allocatable :: wis(:), wir(:)
    real(8), allocatable :: wos(:), wor(:)

    if(j00 == j0) then
       wi(ibi:ifi, 1:j00, kb:kf) = w(ibi:ifi, 1:j00, kb:kf)
       wo(ibo:ifo, 1:j00, kb:kf) = w(ibo:ifo, 1:j00, kb:kf)
       return
    end if

    li = (ifi - ibi + 1) * (kf - kb + 1) * (j0 -  1 + 1)
    li0= (ifi - ibi + 1) * (kf - kb + 1) * (j00-  1 + 1)
    lo = (ifo - ibo + 1) * (kf - kb + 1) * (j0 -  1 + 1)
    lo0= (ifo - ibo + 1) * (kf - kb + 1) * (j00-  1 + 1)

    allocate(wis(0:li - 1), wir(0:li0 - 1))
    allocate(wos(0:lo - 1), wor(0:lo0 - 1))

    if(myrki == 0) then
       l = 0
       do j = 1, j0
       do k = kb , kf 
       do i = ibi, ifi
          wis(l) = w(i,j,k)
          l = l + 1
       end do
       end do
       end do
       if(jprc /= 1) then
          call mpi_allgather(wis, li, mpi_double_precision, &
                             wir, li, mpi_double_precision, mpi_k_world, ier)
       else
          wir(0:li - 1) = wis(0:li - 1)
       end if
       l = 0
       do j = 1, j00
       do k = kb , kf 
       do i = ibi, ifi
          wi(i,j,k) = wir(l)
          l = l + 1
       end do
       end do
       end do
    end if

    if(myrki == iprc - 1) then
       l = 0
       do j = 1, j0
       do k = kb , kf 
       do i = ibo, ifo
          wos(l) = w(i,j,k)
          l = l + 1
       end do
       end do
       end do
       if(jprc /= 1) then
          call mpi_allgather(wos, lo, mpi_double_precision, &
                             wor, lo, mpi_double_precision, mpi_k_world, ier)
       else
          wor(0:lo - 1) = wos(0:lo - 1)
       end if
       l = 0
       do j = 1, j00
       do k = kb , kf 
       do i = ibo, ifo
          wo(i,j,k) = wor(l)
          l = l + 1
       end do
       end do
       end do
    end if

    deallocate(wis, wir)
    deallocate(wos, wor)

    return

  end subroutine gatherx2
#endif




#ifdef DEBUG
  subroutine dwrite(u, i0, j0, k0, n)

    use digit_module, only : digit2

    real(8), intent(in) :: u(0:i0+1,0:j0+1,0:k0+1)
    integer, intent(in) :: i0
    integer, intent(in) :: j0
    integer, intent(in) :: k0
    integer, intent(in) :: n

    real(8) :: x(0:i0+2), y(0:j0+2), z(0:k0+2)
    integer :: i, j, k
    real(8) :: dx, dy, dz

    dx = lx / dble(i0)
    dy = ly / dble(j0)
    dz = lz / dble(k0)

    do i = 0, i0 + 2
       x(i) = dx * (i - 1)
    end do
    do j = 0, j0 + 2
       y(j) = dy * (j - 1)
    end do
    do k = 0, k0 + 2
       z(k) = dz * (k - 1)
    end do

    open(88, file = 'd'//digit2(n)//'.data', form = 'unformatted')
    write(88) i0+3, j0+3, k0+3
    write(88) x, y, z
    write(88) u  
    close(88)

  end subroutine dwrite




  subroutine print2(a, i0, j0)

    real(8), intent(in) :: a(i0,j0)
    integer, intent(in) :: i0
    integer, intent(in) :: j0

    integer :: i, j

    do j = 1, j0
    do i = 1, i0
       write(10+myrk,*) i, j, a(i,j)
    end do
    end do

    return

  end subroutine print2




  subroutine print3(a, i0, j0, k0, im)

    real(8), intent(in) :: a(0:i0+1,0:j0+1,0:k0+1)
    integer, intent(in) :: i0
    integer, intent(in) :: j0
    integer, intent(in) :: k0
    integer, intent(in) :: im

    integer :: i, j, k

    if(im == 0) then
    do j = 0, j0+1
    do k = 0, k0+1
    do i = 0, i0+1
       write(10+myrk,*) i, j, k, a(i,j,k)
    end do
    end do
    end do
    else
    do j = 1, j0
    do k = 1, k0
    do i = 1, i0
       write(10+myrk,*) i, j, k, a(i,j,k)
    end do
    end do
    end do
    end if

    return

  end subroutine print3
#endif



end module diff_module
