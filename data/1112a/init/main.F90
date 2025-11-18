
!#define DFDTAU 1
!#define RADDOMINATED
!#define M     6.00d0
!#define MD    0.5d0
!#define ALPHA 0.05d0
!#define X     100d0
!#define ETA   0.1d0

!#define DFDTAU 0
!#define RADDOMINATED
!#define M     1d8
!#define MD    0.1d0
!#define ALPHA 0.5d0
!#define X     200d0
!#define ETA   0.1d0
!#define ZMAX 4.5d0

!#define DFDTAU 0
!#define RADDOMINATED
!#define M     6.62d0
!#define MD    0.1d0
!#define ALPHA 0.025d0
!#define X     30d0
!#define ETA   0.1d0
!#define ZMAX 6.0d0

!Hirose et al. (2006)
!#define DFDTAU 0
!#define GASDOMINATED
!#define M     6.62d0
!#define MD    0.1d0
!#define ALPHA 0.02d0
!#define X     300d0
!#define ETA   0.1d0
!#define ZMAX 7.6d0

!New Simulation 30
#define DFDTAU 1
#define RADDOMINATED
#define M     6.62d0
#define MD    0.1d0
#define ALPHA 0.025d0
#define X     30d0
#define ETA   0.1d0
#define ZMAX 5.0d0

!New Simulation 65
!#define DFDTAU 0
!#define RADDOMINATED
!#define M     6.62d0
!#define MD    0.1d0
!#define ALPHA 0.023d0
!#define X     30d0
!#define ETA   0.1d0
!#define ZMAX 4.0d0

!Krolik et al. (2007)
!#define DFDTAU 0
!#define GASDOMINATED
!#define M     6.62d0
!#define MD    0.1d0
!#define ALPHA 0.03d0
!#define X     150d0
!#define ETA   0.1d0
!#define ZMAX 5.6d0

! DFDTAU == 1 : Q / \rho = const / \sqrt{\tau}
! DFDTAU == 0 : Q / \rho = const

module vars
  implicit none
  real(8) rho0, tmp0, z0, dfloor
  real(8) f0, omega, tau0
end module vars




module consts
  implicit none
  real(8), parameter :: CLIGHT = 2.99792d+10
  real(8), parameter :: RGAS   = 8.31447d+07
  real(8), parameter :: ARAD   = 7.56760d-15
  real(8), parameter :: GUNIV  = 6.67300d-08
  real(8), parameter :: MSOL   = 1.98900d+33
  real(8), parameter :: PI     = 3.1415926535897932d0
! real(8), parameter :: PI     = acos(-1d0)
  real(8), parameter :: MMW   = 0.60d0   ! mean molecular weight
  real(8), parameter :: CHI   = 0.33d0   ! Thomson opacity
  real(8), parameter :: GAMMA = 5d0 / 3d0 ! specific heat ratio
end module consts




program main
  
  ! -----------------------------------------
  ! vertical profiles of SS disk -SH 01/25/05
  ! -----------------------------------------
  use vars
  use consts
  implicit none

  external :: derivs, erad

  integer, parameter :: nvar = NMAX
  real(8), parameter :: EPSILON = 1d-15
  real(8), parameter :: DMAX = 5d0
  real(8), parameter :: DMIN = 1d0 / DMAX
  real(8), parameter :: EPS = 1d-9
  real(8), parameter :: DH = 2d-5
  real(8), parameter :: HMIN0 = 1d-8
  real(8), parameter :: DFL = 1d-5
    
  real(8) :: p, er, eg, rtmp, tmp, flx, tau, rho
  real(8) :: sigma, h, hlocal
  real(8) :: dd, rr, ff
  real(8) :: z, dz, qdis
  real(8) :: er0, tau00, tauef
  real(8) :: dstart, d1, d2
  real(8) :: rhoold, egold
  integer :: itry, n, k, kount0
  
  integer :: nok, nbad, ierr
  real(8) :: ystart(nvar), x1, x2, h1, hmin
  integer :: kmax, kount
  real(8) :: dxsav, xp(KMAXX), yp(NMAX,KMAXX)
  common /path/ kmax, kount, dxsav, xp, yp
  
  real(8) :: dz1, va1, p1, rho1, lmax
  real(8), parameter :: BETA = 25d0
  real(8), parameter :: BPRATIO = 0.25d0
  integer, parameter :: K0 = 256
  real(8), parameter :: GMM = 2d0, GMG = 5d0/3d0, GMR = 4d0/3d0
  real(8) :: vm, vg, vr, vs, vmax, pm, pg, pr, lx0, cfl, sigmax
  
  ! orbital frequency and surface flux
  omega = 17.d00 * (M / 6.62d0)**(-1d0) &
       * (X / 150d0)**(-1.5d0)
  f0  = 3.7d19 * (M / 6.62d0)**(-1d0) &
       * (X / 150d0)**(-3d0) &
       * (MD / ETA)
  
  ! surface density and scale height (after Krolik 1999)
#ifdef RADDOMINATED
  sigma = 8.1d0 / ALPHA * (ETA / MD) * (X ** (1.5d0))
  h     = 2.2d5 * M * (MD / ETA)
#else
  sigma = 4.7d4 * (ALPHA / 0.02d0)**(-0.8d0) &
       * (MD / ETA)**(0.6d0) &
       * (M / 6.62d0)**(0.2d0) &
       * (X / 150d0)**(-0.6d0)
  h     = 3.1d6 * (ALPHA / 0.02d0)**(-0.1d0) &
       * (MD / ETA)**(0.2d0) &
       * (M / 6.62d0)**(0.9d0) &
       * (X / 150d0)**(1.05d0)
#endif

  print *, 'sigma = ', sigma, 'h =', h
  
  ! opacity at the midplane
  tau0  = sigma * chi + 1d0
  
  !========= solve tau > 1
  !---------
  ! initial guess for density at the midplane
  dstart = sigma / h
  d1 = dstart * DMIN
  d2 = dstart * DMAX
  
  ! integration range: x1 < x < x2
  x1 = tau0
  x2 = 1d0
  
  ! parameters for subroutine odeint
  h1 = (x1 - x2) * DH
  hmin = h1 * HMIN0
  
  !---------
  do ! *** main loop start ***
     
     ! set initial values at x = x1
     ystart(1) = dstart     ! d(tau0)
     ystart(2) = 0          ! z(tau0)
     
     ! set density floor
     dfloor = dstart * DFL
     
     ! integrate the system of differential equations
     call odeint(ystart,2,x1,x2,eps,h1,hmin,nok,nbad,ierr,derivs)
     
     ! if odeint returns error, raise dstart and integrate again
     if(ierr /= 0) then
        d1 = dstart
        dstart = (dstart + d2) / 2d0
        cycle
     endif
     
     ! if drho/dz(tau=1) > 0, lower dstart and integrate again
     if(yp(1,kount) > yp(1,kount-1)) then
        d2 = dstart
        dstart = (dstart + d1) / 2d0
        cycle
     endif
     
     ! lower dstart and integrate again to find smallest dstart
     d2 = dstart
     dstart = (dstart + d1) / 2d0
!     if(abs((dstart - d2) / dstart) < EPSILON) then
     if(abs((dstart - d1) / d1) < EPSILON) then
        exit
     endif
     
  enddo ! *** main loop end ***
  
  !---------
  ! open files for writing z, dens, Erad, and Egas
  open(2, file = 'init/initc.data')
  open(3, file = 'init/initd.data')
  
  tauef = 0d0
  vmax = 0d0
  sigmax = 0d0
  do k = 1, kount
     ! solutions from odeint (rho(tau), z(tau))
     rhoold = rho
     rho = yp(1,k)
     z = yp(2,k)
     tau = xp(k)

     if(k /= kount) then
        sigmax = sigmax + rho * (yp(2,k+1) - yp(2,k))
     endif
     
     ! dissipation rate
#if DFDTAU==0
     qdis = f0 * chi * rho / (tau0 - 1d0)
#elif DFDTAU==1
     qdis = f0 * chi * rho / (2d0 * (sqrt(tau0) - 1d0) * sqrt(tau))
#elif DFDTAU==2
     qdis = f0 * chi * rho / log(tau0) / tau
#endif
     
     ! compute erad, trad, flux from tau
     call erad(er, rtmp, flx, tau)
     
     ! compute other quantities
     tmp = rtmp ! Tgas
     p = rho * tmp * rgas / mmw ! Pgas
     egold = eg
     eg = p / (gamma - 1d0) ! Egas
     call coeff(rr, ff, dd, er, flx, rho)
     
     write(2,*) z, rho, er, eg
     write(3,*) z/h, rho, er
     write(3,*) eg, qdis, flx
     write(3,*) ff, tau, dd
     write(3,*) tmp, rtmp
     
     if(k == 1) then
        rho1 = rho
!       p1 = p
        p1 = p + er / 3d0
        cycle
     endif

     ! compute vmax
     vm = sqrt(gmm * pm / rho)
     vg = sqrt(gmg * pg / rho)
     vr = sqrt(gmr * pr / rho)
     vmax = max(vmax, vm + vg + vr)

     ! compute cumlative optical depth
     dz = (yp(2,k) - yp(2,k-1))
     tauef = tauef + (chid(rhoold,egold) + chid(rho,eg)) * dz / 2d0
  enddo
  close(2)
  close(3)
  
  rho0 = rho
  er0 = er
  z0 = z
  tmp0 = tmp
  tau00 = tau
  kount0 = kount
  
  write(6,*) 'tauef = ', tauef
  write(6,*) 'dstart = ', dstart, ' dfloor = ', rho0
  write(6,*) 'z0/h =', z0 / h
  print *, 'h/r =', 1.5d0 / X * (MD / ETA)

  dz1 = ZMAX * h / dble(K0)
  va1 = sqrt(2d0 * p1 / rho1 / BETA)
  hlocal = sqrt(2d0 * p1 /rho1) / omega
  lmax = 8d0 * PI / sqrt(15d0) * va1 / omega * BPRATIO
  print *, 'lmax/dz = ', lmax / dz1, va1, p1
  print *, 'h / hlocal =', h / hlocal
  
  !---------
  ! write disk parameters
  open(2,file='init/ipara.data')
  write(2,*) omega, h
  write(2,*) sigma, f0
  write(2,*) dstart, tau0
  write(2,*) M, MD, ALPHA, X, ETA
  write(2,*) (rho0 / dstart)
  close(2)
  open(2,file='init/iparam.data')
  write(2,*) '&lnrmcon   h0 = ', h, '   &end'
  write(2,*) '&tnrmcon   omega = ', omega, '   &end'
  write(2,*) '&dflcon   dfl = ', (rho0 / dstart), '   &end'
  close(2)

  print *, '+++ disk part done +++'
  
  !========= solve tau < 1
  !---------
  x1 = z0
  x2 = h * ZMAX
  dz = (x2 - x1) / (kount - 1)
  
  open(2,file='init/initc.data',position='append')
  open(3,file='init/initd.data',position='append')
  
  do k = 2, kount
     z = x1 + dz * (k - 1)
     rho = rho0
     flx = f0
     tmp = tmp0
!    er = er0 - (3d0 * chi * rho / clight) * flx * (z - z0)
     er = er0
     rtmp = (er / arad) ** 0.25d0
     eg = (rgas / mmw * rho * tmp) / (gamma - 1d0) ! Egas
     qdis = 0 ! Qdiss
     tau = tau00 - chi * rho * (z - z0)
     call coeff(rr, ff, dd, er, flx, rho)
     
     write(2,*) z, rho, er, eg
     write(3,*) z/h, rho, er
     write(3,*) eg, qdis, flx
     write(3,*) ff, tau, dd
     write(3,*) tmp, rtmp

     ! compute vmax
     pm = 0d0
     pg = (gmg - 1d0) * eg
     pr = (gmr - 1d0) * er
     vm = sqrt(gmm * pm / rho)
     vg = sqrt(gmg * pg / rho)
     vr = sqrt(gmr * pr / rho)
     vmax = max(vmax, vm + vg + vr)
  enddo
  
  close(3)
  close(2)

  lx0 = z0 / 2d0
  vs = 1.5d0 * omega * lx0
  vmax = max(vmax, vs)
  cfl = 0.5d0
  print *, 'dt/Torb =', (cfl * dz1 / vmax) / (2d0 * PI / omega)
  print *, 'sigmax / sigma =', sigmax / sigma
  
  !---------
  ! write resolution
  open(2,file='init/ireso.data')
  write(2,*) kount0 + kount - 1
  close(2)
  
  stop
  
contains
  
  
  
  
  real(8) function chid(rho, eg)
    real(8), intent(in) :: rho
    real(8), intent(in) :: eg
    real(8) :: absffp, absffr
    
    absffp = 3.7d53 * sqrt(rho ** 9 / eg ** 7)
    absffr = 1.0d52 * sqrt(rho ** 9 / eg ** 7)
    chid = sqrt(chi * absffp) * rho
    
    return
  end function chid
  
  
  
  
  subroutine coeff(rr, ff, dd, er, flx, rho)
    real(8), intent(out) :: rr
    real(8), intent(out) :: ff
    real(8), intent(out) :: dd
    real(8), intent(in) :: er
    real(8), intent(in) :: flx
    real(8), intent(in) :: rho
    real(8) :: lambda

    rr = 3d0 * flx / clight / er ! R
    lambda = (2d0 + rr) / (6d0 + 3d0 * rr + rr ** 2) ! lambda
    ff = lambda + (lambda * rr) ** 2 ! Eddington factor
    dd = clight * lambda / chi / rho ! diffusion coefficient
    
    return
  end subroutine coeff
  



end program main




! external programs

subroutine derivs(x, y, dydx)
  use consts, only : MMW, RGAS, CLIGHT, CHI
  use vars, only : omega, dfloor
  implicit none
  real(8), intent(in) :: x
  real(8), intent(in) :: y(4)
  real(8), intent(out) :: dydx(4)
  real(8) :: er, flx, tmp, rho, z
  
  call erad(er, tmp, flx, x)
  
  rho = y(1)
  z   = y(2)
  
  dydx(1) = MMW / RGAS / tmp &
       * (omega ** 2 / CHI * z - flx / CLIGHT) &
       -0.75d0 * rho / er * flx / CLIGHT
  dydx(2) = -1d0 / (chi * rho)
  
  if(rho < dfloor) then
     dydx(1) = 0
  endif
  
  return
end subroutine derivs




subroutine erad(er, rtmp, flx, tau)
  use consts, only : CLIGHT, ARAD
  use vars, only : f0, tau0
  implicit none
  real(8), intent(out) :: er
  real(8), intent(out) :: rtmp
  real(8), intent(out) :: flx
  real(8), intent(in) :: tau
  
#if DFDTAU==0
  er = 3d0 * f0 / CLIGHT * &
      (2d0 * tau0 * tau - tau**2 - 1d0) / (2d0 * (tau0 - 1d0))
  flx = f0 / (1d0 - tau0) * (tau - tau0)
#elif DFDTAU==1
  er = 3d0 * f0 / CLIGHT * &
       (sqrt(tau0) * tau - (1d0 + 2d0 * tau ** (1.5d0)) / 3d0) / &
       (sqrt(tau0) - 1d0)
  flx = f0 * (sqrt(tau0) - sqrt(tau)) / (sqrt(tau0) - 1d0)
#elif DFDTAU==2
  er  = 3d0 * f0 / CLIGHT * (tau - &
       (tau * log(tau) - tau + 1d0) / log(tau0))
  flx = f0 * (1d0 - log(tau) / log(tau0))
#endif
  rtmp = (er / ARAD) ** 0.25d0
  
  return
end subroutine erad
