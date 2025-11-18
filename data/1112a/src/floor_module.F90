
module floor_module

  implicit none

  private
  public :: floor_d, floor_e, floor_v3

contains




  subroutine floor_d
    use bval_module, only : bvalx1, bvalx2, bvalx3
    use field_module, only : d, v1, v2, v3
    use flux_module, only : sdrh, drhox
    use grid_module, only : is, ie, js, je, ks, ke
    use param_module, only : dfloor
    use root_module, only : dt
    integer :: i, j, k
    real(8) :: ekin
    real(8) :: tma(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: tmb(is-2:ie+3,js-2:je+3,ks-2:ke+3)

    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       ekin = 0.5 * d(i,j,k) * 0.5 * &
            ( v1(i,j,k)**2 + v1(i+1,j,k)**2&
            + v2(i,j,k)**2 + v2(i,j+1,k)**2&
            + v3(i,j,k)**2 + v3(i,j,k+1)**2 )
       tma(i,j,k) = ekin
       tmb(i,j,k) = d(i,j,k)
    enddo
    enddo
    enddo

    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       if(d(i,j,k) < dfloor) then
          d(i,j,k) = dfloor
       endif
    enddo
    enddo
    enddo
 
    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       ekin = 0.5 * d(i,j,k) * 0.5 * &
            ( v1(i,j,k)**2 + v1(i+1,j,k)**2&
            + v2(i,j,k)**2 + v2(i,j+1,k)**2&
            + v3(i,j,k)**2 + v3(i,j,k+1)**2 )
       sdrh(i,j,k)  =                (ekin - tma(i,j,k)) / dt
       drhox(i,j,k) = drhox(i,j,k) + (d(i,j,k) - tmb(i,j,k))
    enddo
    enddo
    enddo

    call bvalx1(d,          ks-2,ke+2, 0)
    call bvalx2(d,is-2,ie+2,ks  ,ke  , 0)
    call bvalx3(d,is-2,ie+2,js-2,je+2, 0)
    
    return
  end subroutine floor_d




  subroutine floor_e(q, qe, e0)
    use bval_module, only : bvalx1, bvalx2, bvalx3
    use field_module, only : e
    use grid_module, only : is, ie, js, je, ks, ke
    use param_module, only : fdamp, efac
    real(8), intent(out) :: q(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8), intent(out) :: qe(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8), intent(in) :: e0(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    integer :: i, j, k
    real(8) :: enew, eold

    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       q(i,j,k) = e0(i,j,k) - e(i,j,k)
       eold = e(i,j,k)
       enew = e0(i,j,k)
       qe(i,j,k) = 0d0
#ifdef FLOORE
       if(e0(i,j,k) < efac * e(i,j,k) .or. e0(i,j,k) < 0d0) then
          enew = e(i,j,k)
          qe(i,j,k) = e(i,j,k) - e0(i,j,k)
       endif
#endif
       e(i,j,k) = enew * (1d0 - fdamp(k)) + eold * fdamp(k)
       qe(i,j,k) = e(i,j,k) - e0(i,j,k)
    enddo
    enddo
    enddo

    call bvalx1(e,          ks-2,ke+2, 0)
    call bvalx2(e,is-2,ie+2,ks  ,ke  , 0)
    call bvalx3(e,is-2,ie+2,js-2,je+2, 0)
    
    return
  end subroutine floor_e




  subroutine floor_v3
    use bval_module, only : bvalx1, bvalx2, bvalx3
    use field_module, only : d, v1, v2, v3
    use flux_module, only : sdv3
    use grid_module, only : is, ie, js, je, ks, ke, lx0
    use param_module, only : fvcap, omega
    use root_module, only : dt
    integer :: i,j,k
    real(8) :: tma(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: ekin
    real(8) :: vcap

    vcap = 1.5 * omega * lx0 * fvcap

    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       ekin = 0.5 * d(i,j,k) * 0.5 * &
            ( v1(i,j,k)**2 + v1(i+1,j,k)**2 &
            + v2(i,j,k)**2 + v2(i,j+1,k)**2 &
            + v3(i,j,k)**2 + v3(i,j,k+1)**2 )
       tma(i,j,k) = ekin
    enddo
    enddo
    enddo

    do k = ks, ke
    do i = is, ie+1
!CDIR SHORT_LOOP
    do j = js, je
       if(v1(i,j,k) >  vcap) v1(i,j,k) =  vcap
       if(v1(i,j,k) < -vcap) v1(i,j,k) = -vcap
    enddo
    enddo
    enddo
    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je+1
       if(v2(i,j,k) >  vcap) v2(i,j,k) =  vcap
       if(v2(i,j,k) < -vcap) v2(i,j,k) = -vcap
    enddo
    enddo
    enddo
    do k = ks, ke+1
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       if(v3(i,j,k) >  vcap) v3(i,j,k) =  vcap
       if(v3(i,j,k) < -vcap) v3(i,j,k) = -vcap
    enddo
    enddo
    enddo

    call bvalx1(v1,          ks  ,ke  , 1)
    call bvalx2(v1,is-2,ie+3,ks  ,ke  , 0)
    call bvalx3(v1,is-2,ie+3,js-2,je+2, 0)
    call bvalx1(v2,          ks  ,ke  , 0, 1)
    call bvalx2(v2,is-2,ie+2,ks  ,ke  , 1)
    call bvalx3(v2,is-2,ie+2,js-2,je+3, 0)
    call bvalx1(v3,          ks  ,ke+1, 0)
    call bvalx2(v3,is-2,ie+2,ks  ,ke+1, 0)
    call bvalx3(v3,is-2,ie+2,js-2,je+2, 1, 0)

    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       ekin = 0.5 * d(i,j,k) * 0.5 * &
             ( v1(i,j,k)**2 + v1(i+1,j,k)**2 &
             + v2(i,j,k)**2 + v2(i,j+1,k)**2 &
             + v3(i,j,k)**2 + v3(i,j,k+1)**2 )
       sdv3(i,j,k) = (ekin - tma(i,j,k)) / dt
    enddo
    enddo
    enddo

    return
  end subroutine floor_v3
  



end module floor_module
