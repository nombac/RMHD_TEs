
module interp_module

  implicit none

  private
  public :: x1intfc, x1intzc, x2intfc, x2intzc, x3intfc, x3intzc

contains




  subroutine x1intfc(q, vel, qi)
    use grid_module, only : dx1a, dx1b, is, ie, js, je, ks, ke
    use root_module, only : dt
    real(8), intent(in) :: q(:,:,:)
    real(8), intent(in) :: vel(:,:,:)
    real(8), intent(out) :: qi(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    integer :: i, j, k
    real(8) :: deltq2, xi
    real(8) :: deltq(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: dq(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    
    do k = ks, ke
    do i = is-1, ie+2
!CDIR SHORT_LOOP
    do j = js, je
       deltq(i,j,k) = (q(i,j,k) - q(i-1,j,k)) / dx1a(i-1)
    enddo
    enddo
    enddo
    do k = ks, ke
    do i = is-1, ie+1
!CDIR SHORT_LOOP
    do j = js, je
       deltq2 = deltq(i,j,k) * deltq(i+1,j,k)
       dq(i,j,k) = 0
       if(deltq2 > 0) dq(i,j,k) = deltq2 / (deltq(i,j,k) + deltq(i+1,j,k))
    enddo
    enddo
    enddo
    do k = ks, ke
    do i = is-1, ie
!CDIR SHORT_LOOP
    do j = js, je
       xi = vel(i,j,k) * dt
       if(vel(i,j,k) >= 0) then
          qi(i,j,k) = q(i  ,j,k) + (dx1b(i  ) - xi) * dq(i  ,j,k)
       else 
          qi(i,j,k) = q(i+1,j,k) - (dx1b(i+1) + xi) * dq(i+1,j,k)
       endif
    enddo
    enddo
    enddo

    return
  end subroutine x1intfc




  subroutine x1intzc(q, vel, qi)
    use grid_module, only : dx1a, dx1b, is, ie, js, je, ks, ke
    use root_module, only : dt
    real(8), intent(in) :: q(:,:,:)
    real(8), intent(in) :: vel(:,:,:)
    real(8), intent(out) :: qi(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    integer :: i, j, k
    real(8) :: deltq2, xi
    real(8) :: deltq(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: dq(is-2:ie+3,js-2:je+3,ks-2:ke+3)

    do k = ks, ke+1
    do i = is-1, ie+2
!CDIR SHORT_LOOP
    do j = js, je
       deltq(i,j,k) = (q(i,j,k) - q(i-1,j,k)) / dx1b(i)
    enddo
    enddo
    enddo
    do k = ks, ke+1
    do i = is-1, ie+1
!CDIR SHORT_LOOP
    do j = js, je
       deltq2 = deltq(i,j,k) * deltq(i+1,j,k)
       dq(i,j,k) = 0
       if(deltq2 > 0) dq(i,j,k) = deltq2 / (deltq(i,j,k) + deltq(i+1,j,k))
    enddo
    enddo
    enddo
    do k = ks, ke+1
    do i = is, ie+1
!CDIR SHORT_LOOP
    do j = js, je
       xi = vel(i,j,k) * dt
       if(vel(i,j,k) >= 0) then
          qi(i,j,k) = q(i-1,j,k) + (dx1a(i-1) - xi) * dq(i-1,j,k)
       else
          qi(i,j,k) = q(i  ,j,k) - (dx1a(i  ) + xi) * dq(i  ,j,k)
       endif
    enddo
    enddo
    enddo
      
    return
  end subroutine x1intzc




  subroutine x2intfc(q, vel, qi)
    use grid_module, only : dx2a, dx2b, is, ie, js, je, ks, ke
    use root_module, only : dt
    real(8), intent(in) :: q(:,:,:)
    real(8), intent(in) :: vel(:,:,:)
    real(8), intent(out) :: qi(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    integer :: i, j, k
    real(8) :: deltq2, xi
    real(8) :: deltq(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: dq(is-2:ie+3,js-2:je+3,ks-2:ke+3)

    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js-1, je+2
       deltq(i,j,k) = (q(i,j,k) - q(i,j-1,k)) / dx2a(j-1)
    enddo
    enddo
    enddo
    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js-1, je+1
       deltq2 = deltq(i,j,k) * deltq(i,j+1,k)
       dq(i,j,k) = 0
       if(deltq2 > 0) dq(i,j,k) = deltq2 / (deltq(i,j,k) + deltq(i,j+1,k))
    enddo
    enddo
    enddo
    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js-1, je
       xi = vel(i,j,k) * dt
       if(vel(i,j,k) >= 0) then
          qi(i,j,k)= q(i,j  ,k) + (dx2b(j  ) - xi) * dq(i,j  ,k)
       else
          qi(i,j,k)= q(i,j+1,k) - (dx2b(j+1) + xi) * dq(i,j+1,k)
       endif
    enddo
    enddo
    enddo

    return
  end subroutine x2intfc




  subroutine x2intzc(q, vel, qi)
    use grid_module, only : dx2a, dx2b, is, ie, js, je, ks, ke
    use root_module, only : dt
    real(8), intent(in) :: q(:,:,:)
    real(8), intent(in) :: vel(:,:,:)
    real(8), intent(out) :: qi(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    integer :: i, j, k
    real(8) :: deltq2, xi
    real(8) :: deltq(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: dq(is-2:ie+3,js-2:je+3,ks-2:ke+3)

    do k = ks  , ke+1
    do i = is  , ie
!CDIR SHORT_LOOP
    do j = js-1, je+2
       deltq(i,j,k) = (q(i,j,k) - q(i,j-1,k)) / dx2b(j)
    enddo
    enddo
    enddo
    do k = ks  , ke+1
    do i = is  , ie
!CDIR SHORT_LOOP
    do j = js-1, je+1
       deltq2 = deltq(i,j,k) * deltq(i,j+1,k)
       dq(i,j,k) = 0d0
       if(deltq2 > 0d0) dq(i,j,k) = deltq2 / (deltq(i,j,k) + deltq(i,j+1,k))
    enddo
    enddo
    enddo
    do k = ks, ke+1
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je+1
       xi = vel(i,j,k) * dt
       if(vel(i,j,k) >= 0d0) then
          qi(i,j,k) = q(i,j-1,k)+ (dx2a(j-1) - xi) * dq(i,j-1,k)
       else
          qi(i,j,k) = q(i,j  ,k)- (dx2a(j  ) + xi) * dq(i,j  ,k)
       endif
    enddo
    enddo
    enddo
 
    return
  end subroutine x2intzc




  subroutine x3intfc(q, vel, qi)
    use grid_module, only : dx3a, dx3b, is, ie, js, je, ks, ke
    use root_module, only : dt
    real(8), intent(in) :: q(:,:,:)
    real(8), intent(in) :: vel(:,:,:)
    real(8), intent(out) :: qi(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    integer :: i, j, k
    real(8) :: deltq2, xi
    real(8) :: deltq(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: dq(is-2:ie+3,js-2:je+3,ks-2:ke+3)

    do k = ks-1, ke+3
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       deltq(i,j,k) = (q(i,j,k) - q(i,j,k-1)) / dx3a(k-1)
    enddo
    enddo
    enddo
    do k = ks-1, ke+2
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       deltq2 = deltq(i,j,k) * deltq(i,j,k+1)
       dq(i,j,k) = 0
       if (deltq2 > 0) dq(i,j,k) = deltq2 / (deltq(i,j,k) + deltq(i,j,k+1))
    enddo
    enddo
    enddo
    do k = ks-1,ke+1
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       xi = vel(i,j,k) * dt
       if(vel(i,j,k) > 0) then
          qi(i,j,k)= q(i,j,k  ) + (dx3b(k  ) - xi) * dq(i,j,k  )
       else
          qi(i,j,k)= q(i,j,k+1) - (dx3b(k+1) + xi) * dq(i,j,k+1)
       endif
    enddo
    enddo
    enddo

    return
  end subroutine x3intfc




  subroutine x3intzc(q, vel, qi)
    use grid_module, only : dx3a, dx3b, is, ie, js, je, ks, ke
    use root_module, only : dt
    real(8), intent(in) :: q(:,:,:)
    real(8), intent(in) :: vel(:,:,:)
    real(8), intent(out) :: qi(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    integer :: i, j, k
    real(8) :: deltq2, xi
    real(8) :: deltq(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: dq(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    
    do k = ks-1, ke+2
    do i = is  , ie
!CDIR SHORT_LOOP
    do j = js  , je
       deltq(i,j,k) = (q(i,j,k) - q(i,j,k-1)) / dx3b(k)
    enddo
    enddo
    enddo
    do k = ks-1, ke+1
    do i = is  , ie
!CDIR SHORT_LOOP
    do j = js  , je
       deltq2 = deltq(i,j,k) * deltq(i,j,k+1)
       dq(i,j,k) = 0d0
       if (deltq2 > 0d0) dq(i,j,k) = deltq2 / (deltq(i,j,k) + deltq(i,j,k+1))
    enddo
    enddo
    enddo
    do k = ks, ke+1
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       xi = vel(i,j,k) * dt
       if(vel(i,j,k) >= 0d0) then
          qi(i,j,k)= q(i,j,k-1) + (dx3a(k-1) - xi) * dq(i,j,k-1)
       else
          qi(i,j,k)= q(i,j,k  ) - (dx3a(k  ) + xi) * dq(i,j,k  )
       endif
   enddo
   enddo
   enddo
 
   return
 end subroutine x3intzc




end module interp_module
