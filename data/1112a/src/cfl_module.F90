
#undef DEBUG

module cfl_module

  implicit none

  private
  public :: cfl

contains




  subroutine cfl
    use field_module, only : d, e, v1, v2, v3, b1, b2, b3, er, &
                             f11, f22, f12, f13, f23
    use grid_module, only : is, ie, js, je, ks, ke, dx1a, dx2a, dx3a, lx0
#ifdef DEBUG
    use grid_module, only : x1b, x2b, x3b
#endif
#ifdef MPI
    use mpi_module
#endif
    use param_module, only : TINY, HUGE, GAMMA, omega, courno, GAMMAR
    use root_module, only : dt, dtqq

    integer :: i, j, k
    real(8) :: qa, qb, fm, dtnew, qsq
#ifdef MPI
    real(8) :: dt0
#endif
#ifdef DEBUG
    real(8) :: q1, q2, q3, q4, q5
    integer :: imin, jmin, kmin
#endif

    qa = GAMMA * (GAMMA - 1.0)
    qb = GAMMAR

    dtnew = HUGE
!   do k = ks, ke
    do k = ks-2, ke+2
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       qsq = sqrt(qa * e(i,j,k) / d(i,j,k))
       fm = max(f11(i,j,k), &
                f22(i,j,k), &
                1d0 - f11(i,j,k) - f22(i,j,k), &
                abs(f12(i,j,k)), &
                abs(f13(i,j,k)), &
                abs(f23(i,j,k)))
       qsq = qsq + sqrt(qb * fm * er(i,j,k) / d(i,j,k))
       qsq = qsq + sqrt(0.5d0 * (b1(i,j,k)**2 + b1(i+1,j,k)**2 &
                               + b2(i,j,k)**2 + b2(i,j+1,k)**2 &
                               + b3(i,j,k)**2 + b3(i,j,k+1)**2) / d(i,j,k))
       qsq = qsq + sqrt(0.5d0 * (v1(i,j,k)**2 + v1(i+1,j,k)**2 &
                               + v2(i,j,k)**2 + v2(i,j+1,k)**2 &
                               + v3(i,j,k)**2 + v3(i,j,k+1)**2))
#ifdef SHEAR
       qsq = qsq + 1.5d0 * omega * lx0
#endif
#ifdef DEBUG
       if(min(dx1a(i), dx2a(j), dx3a(k)) / (qsq + TINY) < dtnew) then
          imin = i
          jmin = j
          kmin = k
          q1 = sqrt(qa * e(i,j,k) / d(i,j,k))
          q2 = sqrt(qb * fm * er(i,j,k) / d(i,j,k))
          q3 = sqrt(0.5d0 * (b1(i,j,k)**2 + b1(i+1,j,k)**2 &
                               + b2(i,j,k)**2 + b2(i,j+1,k)**2 &
                               + b3(i,j,k)**2 + b3(i,j,k+1)**2) / d(i,j,k))
          q4 = sqrt(0.5d0 * (v1(i,j,k)**2 + v1(i+1,j,k)**2 &
                               + v2(i,j,k)**2 + v2(i,j+1,k)**2 &
                               + v3(i,j,k)**2 + v3(i,j,k+1)**2))
          q5 =  qsq + 1.5d0 * omega * lx0
          dtnew = min(dx1a(i), dx2a(j), dx3a(k)) / (qsq + TINY)
       endif
#else
       dtnew = min(dtnew, min(dx1a(i), dx2a(j), dx3a(k)) / (qsq + TINY))
#endif
    enddo
    enddo
    enddo

    dtnew = min(dtnew, dtqq)

    dt = courno * dtnew

#ifdef MPI
    call mpi_allreduce(dt, dt0, 1, mpi_double_precision, mpi_min, &
         mpi_comm_world, ier) 
    dt = dt0
#endif

#ifdef DEBUG
    print *, '*** CFL condition ***'
    print *, 'x =', x1b(imin)
    print *, 'y =', x2b(imin)
    print *, 'z =', x3b(imin)
    print *, 'gas sound speed', (q1 + TINY)
    print *, 'rad sound speed', (q2 + TINY)
    print *, 'mag sound speed', (q3 + TINY)
    print *, 'bulk velocity', (q4 + TINY)
    print *, 'shear velocity', (q5 + TINY)
    print *, '*********************'
#endif

    return
  end subroutine cfl




end module cfl_module
