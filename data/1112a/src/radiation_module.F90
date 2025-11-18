
module radiation_module

#ifdef MPI
  use mpi_module
#endif

  implicit none

  integer :: niter

  private
  public :: radiation_source, fld
  public :: radiation_transport
  public :: niter

contains




  subroutine radiation_source
    use bval_module, only : bvalx1, bvalx2, bvalx3
    use diff_module, only : diff
    use field_module, only : d, e, v1, v2, v3, er, dr, dr1, dr3, fr
    use flux_module, only : srad, sdvp, sena, frd1, frd3
    use grid_module, only : is, ie, js, je, ks, ke, dx, dz, &
                            x1a, x2a, x3a, dx1b, dx2b, dx3b
    use root_module, only : dt
!    use check_module, only : tchck

    integer :: i,j,k
    real(8) :: tma(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: tmb(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: tmc(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: ekin
    real(8) :: df
    integer :: i0,j0,k0
    real(8) :: lx,ly,lz
    real(8), allocatable :: er0(:,:,:)
    real(8), allocatable :: dr0(:,:,:)

!    REAL(8) :: T1, T2

    ! [3] acceleration by radiation force
    !     dv/dt = (chi d F/c)/d
    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       ekin=0.5*d(i,j,k)*0.5*&
           ( v1(i,j,k)**2 + v1(i+1,j,k)**2&
           + v2(i,j,k)**2 + v2(i,j+1,k)**2&
           + v3(i,j,k)**2 + v3(i,j,k+1)**2 )
       tma(i,j,k)=ekin+e(i,j,k)+er(i,j,k)
    enddo
    enddo
    enddo

    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       df = 0.5 * (d(i-1,j,k) + d(i,j,k))
       v1(i,j,k) = v1(i,j,k) + dt * ( &
            (0.5 * (fr(i-1,j,k) + fr(i,j,k))) * &
                   (er(i-1,j,k) - er(i,j,k)) / dx1b(i) / df)
    enddo
    enddo
    enddo
    call bvalx1(v1,          ks  ,ke  , 1)
    call bvalx2(v1,is-2,ie+3,ks  ,ke  , 0)
    call bvalx3(v1,is-2,ie+3,js-2,je+2, 0)

    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       df = 0.5 * (d(i,j-1,k) + d(i,j,k))
       v2(i,j,k) = v2(i,j,k) + dt * ( &
            (0.5 * (fr(i,j-1,k) + fr(i,j,k))) * &
                   (er(i,j-1,k) - er(i,j,k)) / dx2b(j) / df)
    enddo
    enddo
    enddo
    call bvalx1(v2,          ks  ,ke  , 0, 1)
    call bvalx2(v2,is-2,ie+2,ks  ,ke  , 1)
    call bvalx3(v2,is-2,ie+2,js-2,je+3, 0)
    do k = ks, ke+1
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       df = 0.5 * (d(i,j,k-1) + d(i,j,k))
       v3(i,j,k) = v3(i,j,k) + dt * ( &
            (0.5 * (fr(i,j,k-1) + fr(i,j,k))) * &
                   (er(i,j,k-1) - er(i,j,k)) / dx3b(k) / df)
    enddo
    enddo
    enddo
    call bvalx1(v3,          ks  ,ke+1, 0)
    call bvalx2(v3,is-2,ie+2,ks  ,ke+1, 0)
    call bvalx3(v3,is-2,ie+2,js-2,je+2, 1, 0)

! [4] radiation matter interaction
!     de/dt = -(4piB - cE)kappa d
!     dE/dt =  (4piB - cE)kappa d - (grad v):P
    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       tmb(i,j,k) = e(i,j,k) + er(i,j,k)
       tmc(i,j,k) = e(i,j,k)
    enddo
    enddo
    enddo

    call eerint

!   CALL CPU_TIME(T1)
    call COMPTON
!   CALL CPU_TIME(T2)
!   PRINT *, 'COMPTON=', T2 - T1

    call bvalx1(e,          ks-2,ke+2, 0)
    call bvalx2(e,is-2,ie+2,ks  ,ke  , 0)
    call bvalx3(e,is-2,ie+2,js-2,je+2, 0)
    
    call bvalx3(er,is  ,ie  ,js  ,je  , 0)
    call bvalx1(er,          ks-2,ke+2, 0)
    call bvalx2(er,is-2,ie+2,ks-2,ke+2, 0)

    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       ekin = 0.5 * d(i,j,k) * 0.5 * &
            ( v1(i,j,k)**2 + v1(i+1,j,k)**2&
            + v2(i,j,k)**2 + v2(i,j+1,k)**2&
            + v3(i,j,k)**2 + v3(i,j,k+1)**2 )
       srad(i,j,k) = (ekin + e(i,j,k) + er(i,j,k) - tma(i,j,k)) / dt
       sdvp(i,j,k) = (       e(i,j,k) + er(i,j,k) - tmb(i,j,k)) / dt
       sena(i,j,k) = (       e(i,j,k)              -tmc(i,j,k)) / dt
    enddo
    enddo
    enddo
   
! [5] radiation diffusion
!     dE/dt = -div(F)
    lx = x1a(ie+1) - x1a(is)
    ly = x2a(je+1) - x2a(js)
    lz = x3a(ke+1) - x3a(ks)
    i0 = ie - is + 1
    j0 = je - js + 1
    k0 = ke - ks + 1
    allocate(er0(0:i0+1,0:j0+1,0:k0+1))
    allocate(dr0(0:i0+1,0:j0+1,0:k0+1))
    do k = 0, k0+1
    do i = 0, i0+1
!CDIR SHORT_LOOP
    do j = 0, j0+1
       er0(i,j,k) = er(i+is-1,j+js-1,k+ks-1)
       dr0(i,j,k) = dr(i+is-1,j+js-1,k+ks-1)
    enddo
    enddo
    enddo

!   CALL CPU_TIME(T1)
    call diff(er0,dr0,i0,j0,k0,lx,ly,lz,dt,niter)
!   CALL CPU_TIME(T2)
!   TCHCK(3) = T2 - T1
!   PRINT *, 'DIFF=', TCHCK(3)

    do k = 1, k0
    do i = 1, i0
!CDIR SHORT_LOOP
    do j = 1, j0
       er(i+is-1,j+js-1,k+ks-1) = er0(i,j,k)
    enddo
    enddo
    enddo
    deallocate(er0)
    deallocate(dr0)
    call bvalx3(er,is  ,ie  ,js  ,je  , 0, 1)
    call bvalx1(er,          ks-2,ke+2, 0)
    call bvalx2(er,is-2,ie+2,ks-2,ke+2, 0)
    ! diffusioan flux
    do k = ks, ke
    do i = is, ie+1
!CDIR SHORT_LOOP
    do j = js, je
       frd1(i,j,k) = (dr1(i,j,k) * (er(i-1,j,k) - er(i,j,k)) / dx)
    enddo
    enddo
    enddo
    do k = ks, ke+1
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       frd3(i,j,k) = (dr3(i,j,k) * (er(i,j,k-1) - er(i,j,k)) / dz)
    enddo
    enddo
    enddo

    return
  end subroutine radiation_source




  subroutine fld
    use bval_module, only : bvalx3
    use constant_module, only : CLIGHT, SIGMAB, RGAS
    use field_module, only : d, e, er, &
         f11, f22, f12, f13, f23, dr, dr1, dr2, dr3, fr, rrr
    use grid_module, only : is, ie, js, je, ks, ke, dx1b, dx2b, dx3b
    use param_module, only : TINY
    integer :: i,j,k
    real(8) :: kap, sig, der1 ,der2, der3, dernorm, n1, n2, n3, erc, rr
    real(8) :: lmda, ef, omf, tfm1
    real(8) :: qa, qb
    real(8), parameter :: A3RD = 1d0 / 3d0, DIFFE = 1d-6

    do k = ks-1, ke+1
    do i = is-1, ie+1
!CDIR SHORT_LOOP
    do j = js-1, je+1
       kap = (1d52 * sqrt(d(i,j,k) ** 11) / sqrt(e(i,j,k) ** 7))
       sig = (0.33 * d(i,j,k))
       der1 = (er(i+1,j,k) - er(i-1,j,k)) / (dx1b(i) + dx1b(i+1))
       der2 = (er(i,j+1,k) - er(i,j-1,k)) / (dx2b(j) + dx2b(j+1))
       der3 = (er(i,j,k+1) - er(i,j,k-1)) / (dx3b(k) + dx3b(k+1))
       dernorm = sqrt(der1**2 + der2**2 + der3**2)
       n1 = der1 / (dernorm + TINY)
       n2 = der2 / (dernorm + TINY)
       n3 = der3 / (dernorm + TINY)
       erc = 0.125 * (er(i+1,j+1,k+1) + er(i-1,j+1,k+1) + &
                      er(i+1,j-1,k+1) + er(i-1,j-1,k+1) + &
                      er(i+1,j+1,k-1) + er(i-1,j+1,k-1) + &
                      er(i+1,j-1,k-1) + er(i-1,j-1,k-1) )
       rr = dernorm / ((kap + sig) * erc)

       lmda = A3RD
       if(rr > 0d0) then
          qa = 1d0 / tanh(rr)
          qb = 1d0 / rr
          if((qa-qb) > DIFFE * qa) lmda = (qa - qb) / rr
          if(lmda > A3RD) lmda = A3RD
          if(lmda < 0d0) lmda = 0d0
       endif

       dr(i,j,k) = CLIGHT * lmda / (kap + sig)
       fr(i,j,k) = lmda
       ef = lmda + (lmda * rr)**2
       omf = 1.0 - ef
       tfm1 = 3.0 * ef - 1.0
       f11(i,j,k) = 0.5 * (omf + tfm1 * n1 * n1)
       f22(i,j,k) = 0.5 * (omf + tfm1 * n2 * n2)
       f12(i,j,k) = 0.5 * (      tfm1 * n1 * n2)
       f13(i,j,k) = 0.5 * (      tfm1 * n1 * n3)
       f23(i,j,k) = 0.5 * (      tfm1 * n2 * n3)
       rrr(i,j,k) = rr
    enddo
    enddo
    enddo

    call bvalx3(dr,is-1,ie+1,js-1,je+1,0)

    do k = ks, ke
    do i = is, ie+1
!CDIR SHORT_LOOP
    do j = js, je
       dr1(i,j,k) = 0.5 * (dr(i-1,j,k) + dr(i,j,k))
    enddo
    enddo
    enddo
    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je+1
       dr2(i,j,k) = 0.5 * (dr(i,j-1,k) + dr(i,j,k))
    enddo
    enddo
    enddo
    do k = ks, ke+1
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       dr3(i,j,k) = 0.5 * (dr(i,j,k-1) + dr(i,j,k))
    enddo
    enddo
    enddo

    ! z boundary condition for radiation energy
    call bcerkb

    return
  end subroutine fld




  subroutine eerint
    use constant_module, only : CLIGHT, SIGMAB, RGAS
    use field_module, only : d, e, v1, v2, v3, er, f11, f22, f12, f13, f23
    use grid_module, only : is, ie, js, je, ks, ke, &
         dx1a, dx1b, dx2a, dx2b, dx3a, dx3b
    use param_module, only : GAMMA, MMW
    use root_module, only : dt

    real(8) :: gamm1, mult
    real(8) :: dx1, dx2, dx3!, divv
    real(8) :: gv11, gv22, gv33, gv12, gv21, gv13, gv31, gv23, gv32
    real(8) :: e00, er0, kappa
    real(8) :: a1, a2, a3, a4, qty1, qty2, enew, ernew
    real(8) :: c01, c02
    real(8) :: a, b, p, q, t
    integer :: i, j, k

    gamm1 = gamma - 1d0
    mult = 4d0 * SIGMAB * (mmw / RGAS)**4
      
    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       dx1 = dx1b(i+1) + dx1b(i)
       dx2 = dx2b(j+1) + dx2b(j)
       dx3 = dx3b(k+1) + dx3b(k)

       gv11 = (v1(i+1,j,k) - v1(i,j,k)) / dx1a(i)
       gv22 = (v2(i,j+1,k) - v2(i,j,k)) / dx2a(j)
       gv33 = (v3(i,j,k+1) - v3(i,j,k)) / dx3a(k)
!      divv = gv11 + gv22 + gv33
       gv12 = .5d0 * (v2(i+1,j  ,k) - v2(i-1,j  ,k) &
                    + v2(i+1,j+1,k) - v2(i-1,j+1,k)) / dx1
       gv21 = .5d0 * (v1(i  ,j+1,k) - v1(i  ,j-1,k) &
                    + v1(i+1,j+1,k) - v1(i+1,j-1,k)) / dx2
       gv13 = .5d0 * (v3(i+1,j,k  ) - v3(i-1,j,k  ) &
                    + v3(i+1,j,k+1) - v3(i-1,j,k+1)) / dx1
       gv31 = .5d0 * (v1(i  ,j,k+1) - v1(i  ,j,k-1) &
                    + v1(i+1,j,k+1) - v1(i+1,j,k-1)) / dx3
       gv23 = .5d0 * (v3(i,j+1,k  ) - v3(i,j-1,k  ) &
                    + v3(i,j+1,k+1) - v3(i,j-1,k+1)) / dx2
       gv32 = .5d0 * (v2(i,j  ,k+1) - v2(i,j  ,k-1) &
                    + v2(i,j+1,k+1) - v2(i,j+1,k-1)) / dx3

       e00 = e(i,j,k)
       er0 = er(i,j,k)

       kappa = (37d52 * sqrt(d(i,j,k) ** 11) / sqrt(e(i,j,k) ** 7))
       a1 = dt * mult * kappa * (gamm1 / d(i,j,k))**4
       a2 = dt * CLIGHT * kappa
       a3 = dt * (gv11 * f11(i,j,k) + gv22 * f22(i,j,k) &
             +  gv33 * (1d0 - f11(i,j,k) - f22(i,j,k)) &
             + (gv12 + gv21) * f12(i,j,k) &
             + (gv13 + gv31) * f13(i,j,k) &
             + (gv23 + gv32) * f23(i,j,k))
       a4 = 0d0
!      a4 = dt * gamm1 * divv

       qty1 = 1d0 + a2 + a3
       qty2 = a1 * (1d0 + a3)

       c01 = -(qty1 * e00 + a2 * er0) / qty2
       c02 = qty1 * (1d0 + a4) / qty2

!     solve (enew)^4 + c02 (enew) + c01 == 0 using Ferrari's formula
       p = -4d0 / 3d0 * c01
       q = -c02**2
       a = (sqrt(q**2 + 4d0 * p**3) - q) / 2d0
       b = (sqrt(q**2 + 4d0 * p**3) + q) / 2d0
       t = a**(1d0 / 3d0) - b**(1d0 / 3d0)
       enew = (-sqrt(t) + sqrt(-t + 2d0 * c02 / sqrt(t))) / 2d0

       ernew = (er0 + e00 - enew * (1d0 + a4)) / (1d0 + a3)
       e(i,j,k) = enew
       er(i,j,k) = ernew
    enddo
    enddo
    enddo

    return
  end subroutine eerint




  subroutine compton_sca

    use constant_module, only : CLIGHT, KBOLTZ, ARAD, MELE, MPRO, SIGMAT
    use field_module, only : d, e, er
    use grid_module, only : is, ie, js, je, ks, ke
    use param_module, only : GAMMA, MMW, MMWE
    use root_module, only : dt

    real(8), parameter :: AVEM = MMW * MPRO
    real(8), parameter :: MUE = MMWE * MPRO
    real(8), parameter :: EPS = 1d-6
    integer, parameter :: ITERMAX = 100
    real(8) :: s, t
    real(8) :: f0, dx, fprime
    real(8) :: x0, x02, x03, x04
    real(8) :: tcomb
    real(8) :: rho, erad, egas, erad1, egas1
    integer :: i, j, k
    integer :: iter

    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       erad = er(i,j,k)
       egas = e(i,j,k)
       rho = d(i,j,k)
       ! Two dimensionless ratios enter into the polynomial, S and T
       s = sqrt(sqrt(erad / ARAD)) * KBOLTZ * rho / ((GAMMA - 1d0) * AVEM)
       s = s / erad
       t = 0.25d0 * MELE * CLIGHT * MUE / AVEM
       t = t / ((GAMMA - 1d0) * SIGMAT)
       t = t / (dt * erad)
       tcomb = t - 1d0 - egas / erad
       
       ! Make the initial guess that the new Erad is equal to the old one.
       x0 = 1d0
       
       ! Evaluate the polynomial
       do iter = 1, ITERMAX
          x02 = x0 * x0
          x03 = x0 * x02
          x04 = x0 * x03
          f0 = x04 * (x04 + s * x0 + tcomb) - t
          
          ! Test for convergence and return if satisfied
          if(abs(f0) < EPS) then
             erad1 = x04 * erad
             egas1 = (egas + erad) - erad1
             er(i,j,k) = erad1
             e(i,j,k) = egas1
             exit
          endif
          
          ! If not converged, compute the next iteration
          fprime = x03 * (8d0 * x04 + 5d0 * s * x0 + 4d0 * tcomb)
          
          ! In some circumstances, there can be a minimum of f 
          ! at positive x whose value is negative.  When that occurs, 
          ! the physical root is at larger x, yet the Newton-Raphson 
          ! algorithm will move toward smaller x.
!!$!         if((fprime < 0d0) .and. (f0 < 0d0)) then
!!$!            x0 = 2d0 * x0
!!$!         else
!!$!            ! If the special test isn't satisfied, 
!!$!            ! the normal Newton-Raphson improvement should be correct.
!!$!            dx = -f0 / fprime
!!$!            x0 = x0 + dx
!!$!            ! By definition, negative x is unphysical.
!!$!            if(x0 < 0d0) then
!!$!               x0 = abs(x0)
!!$!            endif
!!$!         endif
          dx = -f0 / fprime
          x0 = x0 + dx * sign(1d0, sign(1d0, f0) + sign(1d0, fprime) + 1.5d0)
       enddo
       ! Don't permit more than itermax iterations
       if(abs(f0) >= eps) then
          print *, '*** compton: not converged ***'
#ifdef MPI
          call mpi_finalize(ier)
#endif
          stop
       endif
    enddo
    enddo
    enddo

    return
  end subroutine compton_sca




  subroutine compton_vec

    use constant_module, only : CLIGHT, KBOLTZ, ARAD, MELE, MPRO, SIGMAT
    use field_module, only : d, e, er
    use grid_module, only : is, ie, js, je, ks, ke
    use param_module, only : GAMMA, MMW, MMWE
    use root_module, only : dt
    
    real(8), parameter :: AVEM = MMW * MPRO
    real(8), parameter :: MUE = MMWE * MPRO
    real(8), parameter :: EPS = 1d-6
    integer, parameter :: ITERMAX = 100
    real(8) :: s, t
    real(8) :: f0, dx, fprime
    real(8) :: x0, x02, x03, x04
    real(8) :: tcomb
    real(8) :: rho, erad, egas, erad1, egas1
    integer :: i, j, k
    integer :: iter
    real(8) :: tcomb_w(js:je)
    real(8) :: x0_w(js:je)
    real(8) :: t_w(js:je), s_w(js:je)
    integer :: isw(js:je)
    integer :: aitermax = 0
  
    do k = ks, ke
    do i = is, ie
       do j = js, je
          erad = er(i,j,k)
          egas = e(i,j,k)
          rho = d(i,j,k)
          s = sqrt(sqrt(erad / arad)) * KBOLTZ * rho / ((GAMMA - 1d0) * AVEM)
          s = s / erad
          t = 0.25d0 * MELE * CLIGHT * MUE / AVEM
          t = t / ((GAMMA - 1d0) * SIGMAT)
          t = t / (dt * erad)
          tcomb_w(j) = t - 1d0 - egas / erad
          x0_w(j) = 1d0
          t_w(j) = t
          s_w(j) = s
       end do
     
       isw(js:je) = 0
     
       do iter = 1, ITERMAX

          do j = js, je
             x0    = x0_w(j)
             tcomb = tcomb_w(j)
             t = t_w(j)
             s = s_w(j)
             x02 = x0 * x0
             x03 = x0 * x02
             x04 = x0 * x03
             f0 = x04 * (x04 + s * x0 + tcomb) - t
             
             if((abs(f0) < EPS) .and. (isw(j) == 0)) then
                erad = er(i,j,k)
                egas = e(i,j,k)
                erad1 = x04 * erad
                egas1 = (erad + egas) - erad1
                er(i,j,k) = erad1
                e(i,j,k) = egas1
                isw(j) = 1
             endif
             
             fprime = x03 * (8d0 * x04 + 5d0 * s * x0 +  4d0 * tcomb)

!            if((fprime < 0d0) .and. (f0 < 0d0)) then
!               x0 = 2d0 * x0
!            else
!               dx = -f0 / fprime
!               x0 = x0 + dx
!               if(x0 < 0d0) then
!                  x0 = abs(x0)
!               endif
!            endif
             dx = -f0 / fprime
             x0 = x0 + dx * sign(1d0, sign(1d0, f0) + sign(1d0, fprime) + 1.5d0)

             x0_w(j) = x0
          end do

          if(sum(isw) == (je - js + 1)) exit

       enddo

       aitermax = max(aitermax, iter)

    enddo
    enddo

    if(aitermax >= ITERMAX) then
       print *, '*** compton: not converged ***'
#ifdef MPI
       call mpi_finalize(ier)
#endif
       stop
    endif

    return
  end subroutine compton_vec




   subroutine bcerkb
     ! boundary condition!
     use field_module, only : er, dr3, erikb, erokb
     use grid_module, only : is, ie, js, je, ks, ke, dz
     integer :: i, j
     real(8) :: xi, xo, a0, b0, c0

#ifdef MPI
     if(myrkk == 0) then
#endif
        do i = is, ie
!CDIR SHORT_LOOP
        do j = js, je
           a0 = ( dr3(i,j,ks  ) / dz)
           b0 = ( dr3(i,j,ks+1) / dz) - (-dr3(i,j,ks  ) / dz)
           c0 =                         (-dr3(i,j,ks+1) / dz)
           xi = b0 / a0 + c0 / a0 * er(i,j,ks+1) / er(i,j,ks)
           if(xi < 0) xi = 0
           if(xi > 1) xi = 1
           erikb(i,j) = xi
        enddo
        enddo
#ifdef MPI
     endif
     if(myrkk == kprc-1) then
#endif
        do i = is, ie
!CDIR SHORT_LOOP
        do j = js, je
           a0 = (-dr3(i,j,ke+1) / dz)
           b0 = (-dr3(i,j,ke  ) / dz) - ( dr3(i,j,ke+1) / dz)
           c0 =                         ( dr3(i,j,ke  ) / dz)
           xo = b0 / a0 + c0 / a0 * er(i,j,ke-1) / er(i,j,ke)
           if(xo < 0) xo = 0
           if(xo > 1) xo = 1
           erokb(i,j) = xo
        enddo
        enddo
#ifdef MPI
     endif
#endif
     
     return
   end subroutine bcerkb




   subroutine radiation_transport
     call tranr1
     call tranr2
     call tranr3
     return
   end subroutine radiation_transport




   subroutine tranr1
     use bval_module, only : bvalx1, bvalx2, bvalx3
     use field_module, only : v1, er
     use flux_module, only : fra1
     use grid_module, only : is, ie, js, je, ks, ke, dx
     use interp_module, only : x1intzc
     use root_module, only : dt
     integer :: i, j, k
     real(8) :: ertwid(is-2:ie+3,js-2:je+3,ks-2:ke+3)
     real(8) :: erflx(is-2:ie+3,js-2:je+3,ks-2:ke+3)

     call x1intzc(er,v1,ertwid)
     do k = ks, ke
     do i = is, ie+1
!CDIR SHORT_LOOP
     do j = js, je
        erflx(i,j,k) = v1(i,j,k) * dt * ertwid(i,j,k)
        fra1(i,j,k) = erflx(i,j,k) / dt
     enddo
     enddo
     enddo
     do k = ks, ke
     do i = is, ie
!CDIR SHORT_LOOP
     do j = js, je
        er(i,j,k) = er(i,j,k) - (erflx(i+1,j,k) - erflx(i,j,k)) / dx
     enddo
     enddo
     enddo
 
     call bvalx3(er,is  ,ie  ,js  ,je  , 0)
     call bvalx1(er,          ks-2,ke+2, 0)
     call bvalx2(er,is-2,ie+2,ks-2,ke+2, 0)
 
     return
   end subroutine tranr1




   subroutine tranr2 
     use bval_module, only : bvalx1, bvalx2, bvalx3
     use field_module, only : v2, er
     use grid_module, only : is, ie, js, je, ks, ke, dy
     use interp_module, only : x2intzc
     use root_module, only : dt
     integer :: i,j,k
     real(8) :: ertwid(is-2:ie+3,js-2:je+3,ks-2:ke+3)
     real(8) :: erflx(is-2:ie+3,js-2:je+3,ks-2:ke+3)
 
     call x2intzc(er,v2,ertwid)
     do k = ks, ke
     do i = is, ie
!CDIR SHORT_LOOP
     do j = js, je+1
        erflx(i,j,k) = v2(i,j,k) * dt * ertwid(i,j,k)
     enddo
     enddo
     enddo
     do k = ks, ke
     do i = is, ie
!CDIR SHORT_LOOP
     do j = js, je
        er(i,j,k) = er(i,j,k) - (erflx(i,j+1,k) - erflx(i,j,k)) / dy
     enddo
     enddo
     enddo
 
     call bvalx3(er,is  ,ie  ,js  ,je  , 0)
     call bvalx1(er,          ks-2,ke+2, 0)
     call bvalx2(er,is-2,ie+2,ks-2,ke+2, 0)
 
     return
   end subroutine tranr2
 
 
 
 
   subroutine tranr3
     use bval_module, only : bvalx1, bvalx2, bvalx3
     use field_module, only : v3, er
     use flux_module, only : fra3
     use grid_module, only : is, ie, js, je, ks, ke, dz
     use interp_module, only : x3intzc
     use root_module, only : dt
     integer :: i, j, k
     real(8) :: ertwid(is-2:ie+3,js-2:je+3,ks-2:ke+3)
     real(8) :: erflx(is-2:ie+3,js-2:je+3,ks-2:ke+3)
 
     call x3intzc(er,v3,ertwid)
     do k = ks, ke+1
     do i = is, ie
!CDIR SHORT_LOOP
     do j = js, je
        erflx(i,j,k) = v3(i,j,k) * dt * ertwid(i,j,k)
        fra3(i,j,k) = erflx(i,j,k) / dt
     enddo
     enddo
     enddo
     do k = ks, ke
     do i = is, ie
!CDIR SHORT_LOOP
     do j = js, je
        er(i,j,k) = er(i,j,k) - (erflx(i,j,k+1) - erflx(i,j,k)) / dz
     enddo
     enddo
     enddo
 
     call bvalx3(er,is  ,ie  ,js  ,je  , 0)
     call bvalx1(er,          ks-2,ke+2, 0)
     call bvalx2(er,is-2,ie+2,ks-2,ke+2, 0)
 
     return
   end subroutine tranr3




end module radiation_module
