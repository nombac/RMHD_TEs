
module convert_module
  
  implicit none

  private
  public :: convert_e2et, convert_et2e

contains




  subroutine convert_e2et(ieg)
    use field_module, only : d, e, v1, v2, v3, b1, b2, b3, et
    use grid_module, only : is, ie, js, je, ks, ke
    integer, intent(in) :: ieg
    integer :: i, j, k
    real(8) :: ekin, emag

    if(ieg == 0) then
       do k = ks, ke
       do i = is, ie
!CDIR SHORT_LOOP
       do j = js, je
          ekin = 0.5d0 * d(i,j,k) * 0.5d0 * &
               ( v1(i,j,k)**2 + v1(i+1,j,k)**2 &
               + v2(i,j,k)**2 + v2(i,j+1,k)**2 &
               + v3(i,j,k)**2 + v3(i,j,k+1)**2 )
          et(i,j,k) = e(i,j,k) + ekin
       enddo
       enddo
       enddo
    else if(ieg == 1) then
       do k = ks, ke
       do i = is, ie
!CDIR SHORT_LOOP
       do j = js, je
          ekin = 0.5d0 * d(i,j,k) * 0.5d0 * &
               ( v1(i,j,k)**2 + v1(i+1,j,k)**2 &
               + v2(i,j,k)**2 + v2(i,j+1,k)**2 &
               + v3(i,j,k)**2 + v3(i,j,k+1)**2 )
          emag = 0.5d0 * 0.5d0 * &
               ( b1(i,j,k)**2 + b1(i+1,j,k)**2 &
               + b2(i,j,k)**2 + b2(i,j+1,k)**2 &
               + b3(i,j,k)**2 + b3(i,j,k+1)**2 )
          et(i,j,k) = e(i,j,k) + ekin + emag
       enddo
       enddo
       enddo
    endif

    return
  end subroutine convert_e2et




  subroutine convert_et2e(ieg,tma,tmb)
    use field_module, only : d, v1, v2, v3, b1, b2, b3, et
    use floor_module, only : floor_e
    use grid_module, only : is, ie, js, je, ks, ke
    integer, intent(in) :: ieg
    real(8), intent(out) :: tma(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8), intent(out) :: tmb(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: e0(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    integer :: i, j, k
    real(8) :: ekin, emag

    if(ieg == 0) then
       do k = ks, ke
       do i = is, ie
!CDIR SHORT_LOOP
       do j = js, je
          ekin = 0.5d0 * d(i,j,k) * 0.5d0 * &
               ( v1(i,j,k)**2 + v1(i+1,j,k)**2 &
               + v2(i,j,k)**2 + v2(i,j+1,k)**2 &
               + v3(i,j,k)**2 + v3(i,j,k+1)**2 )
          e0(i,j,k) = et(i,j,k) - ekin
       enddo
       enddo
       enddo
    else if(ieg == 1) then
       do k = ks, ke
       do i = is, ie
!CDIR SHORT_LOOP
       do j = js, je
          ekin = 0.5d0 * d(i,j,k) * 0.5d0 * &
               ( v1(i,j,k)**2 + v1(i+1,j,k)**2 &
               + v2(i,j,k)**2 + v2(i,j+1,k)**2 &
               + v3(i,j,k)**2 + v3(i,j,k+1)**2 )
          emag = 0.5d0 * 0.5d0 * &
               ( b1(i,j,k)**2 + b1(i+1,j,k)**2 &
               + b2(i,j,k)**2 + b2(i,j+1,k)**2 &
               + b3(i,j,k)**2 + b3(i,j,k+1)**2 )
          e0(i,j,k) = et(i,j,k) - ekin - emag
       enddo
       enddo
       enddo
    endif

    call floor_e(tma,tmb,e0)

    return
  end subroutine convert_et2e




end module convert_module
