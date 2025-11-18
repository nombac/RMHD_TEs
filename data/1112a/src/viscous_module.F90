
module viscous_module

  implicit none

  private
  public :: viscous

contains



  subroutine viscous
    use bval_module, only : bvalx1, bvalx2, bvalx3
    use convert_module, only : convert_et2e, convert_e2et
    use field_module, only : d, e, v1, v2, v3, et, scol, qq1, qq2, qq3, qkin
    use flux_module, only : fgw1, fgw3, sdek
    use grid_module, only : is, ie, js, je, ks, ke, dx, dy, dz
    use interp_module, only : x1intzc, x2intzc, x3intzc
#ifdef MPI
    use mpi_module
#endif
    use param_module, only : TINY, HUGE, qcon
    use root_module, only : dt, dtqq
    ! von Neumann-Richtmyer artificial viscosity
    ! [1] acceleration by pressure gradient
    !     dv/dt = -(grad q)/d
    ! [2] compression heating 
    !     de/dt = -q(div v)
    ! *** energy conservation ***
    ! [1]+[2]
    ! d(etot)/dt = -div(qv)
    ! ***************************
    integer :: i, j, k
    real(8) :: di,dvmin,qa
    real(8) :: dv1(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: q1(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: dv2(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: q2(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: dv3(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: q3(is-2:ie+3,js-2:je+3,ks-2:ke+3)
#ifdef MPI
    real(8) :: dvmin0
#endif
    real(8) :: qq1flx(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: qq2flx(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: qq3flx(is-2:ie+3,js-2:je+3,ks-2:ke+3)

    real(8) :: qq1twid(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: qq2twid(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: qq3twid(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: tma(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: tmb(is-2:ie+3,js-2:je+3,ks-2:ke+3)

    qa = dt * qcon

#ifdef PARTETOTK
    call convert_e2et(0)

    ! d(etot)/dt = -div(qv)

    do k = ks  , ke
    do i = is-2, ie+2
!CDIR SHORT_LOOP
    do j = js  , je
       dv1(i,j,k) = v1(i+1,j,k) - v1(i,j,k)
       if (dv1(i,j,k) > 0d0) dv1(i,j,k) = 0d0
       q1(i,j,k) = qa * d(i,j,k) * dv1(i,j,k)**2
       qq1(i,j,k) = q1(i,j,k) / dt
    enddo
    enddo
    enddo
    call x1intzc(qq1,v1,qq1twid)
    do k = ks, ke
    do i = is, ie+1
!CDIR SHORT_LOOP
    do j = js, je
       qq1flx(i,j,k) = v1(i,j,k) * dt * qq1twid(i,j,k)
       fgw1(i,j,k) = fgw1(i,j,k) + qq1flx(i,j,k) / dt
    enddo
    enddo
    enddo
    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       et(i,j,k) = et(i,j,k) - (qq1flx(i+1,j,k) - qq1flx(i,j,k)) / dx
    enddo
    enddo
    enddo

    do k = ks  , ke
    do i = is  , ie
!CDIR SHORT_LOOP
    do j = js-2, je+2
       dv2(i,j,k) = v2(i,j+1,k) - v2(i,j,k)
       if(dv2(i,j,k) > 0d0) dv2(i,j,k) = 0d0
       q2(i,j,k) = qa*d(i,j,k) * dv2(i,j,k)**2
       qq2(i,j,k) = q2(i,j,k) / dt
    enddo
    enddo
    enddo
    call x2intzc(qq2,v2,qq2twid)
    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je+1
       qq2flx(i,j,k) = v2(i,j,k) * dt * qq2twid(i,j,k)
    enddo
    enddo
    enddo
    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       et(i,j,k) = et(i,j,k) - (qq2flx(i,j+1,k) - qq2flx(i,j,k)) / dy
    enddo
    enddo
    enddo

    do k = ks-2, ke+2
    do i = is  , ie
!CDIR SHORT_LOOP
    do j = js  , je
       dv3(i,j,k) = v3(i,j,k+1) - v3(i,j,k)
       if(dv3(i,j,k) > 0d0) dv3(i,j,k) = 0d0
       q3(i,j,k) = qa*d(i,j,k) * dv3(i,j,k)**2
       qq3(i,j,k) = q3(i,j,k) / dt
    enddo
    enddo
    enddo
    call x3intzc(qq3,v3,qq3twid)
    do k = ks, ke+1
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       qq3flx(i,j,k) = v3(i,j,k) * dt * qq3twid(i,j,k)
       fgw3(i,j,k) = fgw3(i,j,k) + qq3flx(i,j,k) / dt
    enddo
    enddo
    enddo
    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       et(i,j,k) = et(i,j,k) - (qq3flx(i,j,k+1) - qq3flx(i,j,k)) / dz
    enddo
    enddo
    enddo
#endif

! dv/dt = -(grad q)/d
! de/dt = -q(div v)
    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       tma(i,j,k) = e(i,j,k)
    enddo
    enddo
    enddo

    dvmin = 0d0
!  Start von Neumann-Richtmyer artificial viscosity update.
!  Begin with v1
    do k = ks  , ke
    do i = is-2, ie+2
!CDIR SHORT_LOOP
    do j = js  , je
       dv1(i,j,k) = v1(i+1,j,k) - v1(i,j,k)
       if (dv1(i,j,k) > 0d0) dv1(i,j,k) = 0d0
       q1(i,j,k) = qa * d(i,j,k) * dv1(i,j,k)**2
    enddo
    enddo
    enddo
    do k = ks  , ke
    do i = is-1, ie+2
!CDIR SHORT_LOOP
    do j = js  , je
       di = 2d0 / (d(i-1,j,k) + d(i,j,k))
       e(i,j,k)= e(i,j,k)- q1(i,j,k) * dv1(i,j,k) / dx
       v1(i,j,k) = v1(i,j,k) - di * (q1(i,j,k) - q1(i-1,j,k)) / dx
       dvmin = min(dvmin,(dv1(i,j,k) / dx))
    enddo
    enddo
    enddo
    call bvalx1(v1,          ks  ,ke  , 1)
    call bvalx2(v1,is-2,ie+3,ks  ,ke  , 0)
    call bvalx3(v1,is-2,ie+3,js-2,je+2, 0)
#ifdef MPI
    call mpi_allreduce(dvmin,dvmin0,1,mpi_double_precision,mpi_min,&
         mpi_comm_world,ier)
    dvmin = dvmin0
#endif
!
!  Now do v2 
!
    do k = ks  , ke
    do i = is  , ie
!CDIR SHORT_LOOP
    do j = js-2, je+2
       dv2(i,j,k) = v2(i,j+1,k) - v2(i,j,k)
       if (dv2(i,j,k) > 0d0) dv2(i,j,k) = 0d0
       q2(i,j,k) = qa * d(i,j,k) * dv2(i,j,k)**2
    enddo
    enddo
    enddo
    do k = ks  , ke
    do i = is  , ie
!CDIR SHORT_LOOP
    do j = js-1, je+2
       di = 2d0 / (d(i,j-1,k) + d(i,j,k))
       e(i,j,k) = e(i,j,k) - q2(i,j,k) * dv2(i,j,k) / dy
       v2(i,j,k) = v2(i,j,k) - di*(q2(i,j,k) - q2(i,j-1,k)) / dy
       dvmin = min(dvmin, (dv2(i,j,k) / dy))
    enddo
    enddo
    enddo

    call bvalx1(v2,          ks  ,ke  , 0, 1)
    call bvalx2(v2,is-2,ie+2,ks  ,ke  , 1)
    call bvalx3(v2,is-2,ie+2,js-2,je+3, 0)
#ifdef MPI
    call mpi_allreduce(dvmin,dvmin0,1,mpi_double_precision,mpi_min,&
         mpi_comm_world,ier)
    dvmin = dvmin0
#endif
!
!  Now do v3 
!
    do k = ks-2, ke+2
    do i = is  , ie
!CDIR SHORT_LOOP
    do j = js  , je
       dv3(i,j,k) = v3(i,j,k+1) - v3(i,j,k)
       if(dv3(i,j,k) > 0d0) dv3(i,j,k) = 0d0
       q3(i,j,k) = qa * d(i,j,k) * dv3(i,j,k)**2
    enddo
    enddo
    enddo
    do k = ks-1, ke+2
    do i = is  , ie
!CDIR SHORT_LOOP
    do j = js  , je
       di = 2d0 / (d(i,j,k-1) + d(i,j,k))
       e(i,j,k) = e(i,j,k) - q3(i,j,k) * dv3(i,j,k) / dz
       v3(i,j,k) = v3(i,j,k) - di * (q3(i,j,k) - q3(i,j,k-1)) / dz
       dvmin = min(dvmin, (dv3(i,j,k) / dz))
    enddo
    enddo
    enddo
    call bvalx1(v3,          ks  ,ke+1, 0)
    call bvalx2(v3,is-2,ie+2,ks  ,ke+1, 0)
    call bvalx3(v3,is-2,ie+2,js-2,je+2, 1, 0)
#ifdef MPI
    call mpi_allreduce(dvmin,dvmin0,1,mpi_double_precision,mpi_min,&
         mpi_comm_world,ier)
    dvmin = dvmin0
#endif
!  viscus timestep.  Note we have found the minimum dv/dx above
!  since it is < 0 ; thus minimum dv/dx gives maximum absolute value.
    dtqq = min(1d0 / abs(4d0 * qcon * (dvmin - tiny)), huge)

    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       scol(i,j,k) =                 (e(i,j,k) - tma(i,j,k)) / dt
    enddo
    enddo
    enddo

    call bvalx1(e,          ks-2,ke+2, 0)
    call bvalx2(e,is-2,ie+2,ks  ,ke  , 0)
    call bvalx3(e,is-2,ie+2,js-2,je+2, 0)

#ifdef PARTETOTK
    call convert_et2e(0,tma,tmb)
    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       qkin(i,j,k) = qkin(i,j,k) + tma(i,j,k) / dt
       sdek(i,j,k) = sdek(i,j,k) + tmb(i,j,k) / dt
    enddo
    enddo
    enddo
#endif

    return
  end subroutine viscous



  
end module viscous_module
