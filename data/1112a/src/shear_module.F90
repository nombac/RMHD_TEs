module shear_module

  implicit none

  private
  public :: shear_source

contains




  subroutine shear_source
    use bval_module, only : bvalx1, bvalx2, bvalx3
    use field_module, only : d, v1, v2, v3
    use flux_module, only : spot
    use grid_module, only : is, ie, js, je, ks, ke, x1a, x3a
    use root_module, only : dt
    use param_module, only : omega
    integer :: i, j, k
    real(8) :: tma(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: ekin
! [1] acceleration by Colioris force
!     dv/dt = - 2 omega x v
! [2] acceleration by tidal and gravitational force 
!     dv/dt = - grad(phi)

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
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       v1(i,j,k) = v1(i,j,k) + dt * ( 3.* omega**2 * x1a(i) + &
        0.5 * omega * (v2(i,j,k) + v2(i-1,j,k) + v2(i-1,j+1,k) + v2(i,j+1,k)))
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
       v2(i,j,k) = v2(i,j,k) + dt * ( &
       -0.5 * omega * (v1(i,j,k) + v1(i,j-1,k) + v1(i+1,j-1,k) + v1(i+1,j,k)))
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
       v3(i,j,k) = v3(i,j,k) + dt * (-omega**2 * x3a(k))
    enddo
    enddo
    enddo
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
       spot(i,j,k) = (ekin - tma(i,j,k)) / dt
    enddo
    enddo
    enddo         

    return
  end subroutine shear_source




end module shear_module
