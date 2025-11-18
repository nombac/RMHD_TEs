
module transport_module

  implicit none

  private
  public :: transport
  
contains




  subroutine transport
    use radiation_module, only : radiation_transport

    call radiation_transport
    call hydro_transport

    return
  end subroutine transport




  subroutine hydro_transport
    use convert_module, only : convert_e2et, convert_et2e
    use field_module, only : qkin
    use floor_module, only : floor_d
    use flux_module, only : sdek
    use grid_module, only : is, ie, js, je, ks, ke
    use root_module, only : dt

    integer i, j, k
#ifdef PARTETOTK
    real(8) :: tma(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: tmb(is-2:ie+3,js-2:je+3,ks-2:ke+3)
#endif

#ifdef PARTETOTK
    call convert_e2et(0)
#endif
    call vtos

    call tranx1
    call tranx2
    call tranx3

    call trans1
    call trans2
    call trans3

    call stov

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
#ifdef FLOORD
    call floor_d
#endif
    
    return
  end subroutine hydro_transport




  subroutine vtos
    use bval_module, only : bvalx1, bvalx2, bvalx3
    use field_module, only : d, v1, v2, v3, s1, s2, s3
    use grid_module, only : is, ie, js, je, ks, ke   
    integer :: i, j, k

    do k = ks, ke
    do i = is, ie+1
!CDIR SHORT_LOOP
    do j = js, je
       s1(i,j,k) = 0.5d0 * (d(i-1,j,k) + d(i,j,k)) * v1(i,j,k)
    enddo
    enddo
    enddo
    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je+1
       s2(i,j,k) = 0.5d0 * (d(i,j-1,k) + d(i,j,k)) * v2(i,j,k)
    enddo
    enddo
    enddo
    do k = ks, ke+1
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       s3(i,j,k) = 0.5d0 *(d(i,j,k-1) + d(i,j,k)) * v3(i,j,k)
    enddo
    enddo
    enddo

    call bvalx1(s1,          ks  ,ke  , 1)
    call bvalx2(s1,is-2,ie+2,ks  ,ke  , 0)
    call bvalx3(s1,is-2,ie+2,js-2,je+2, 0)
    call bvalx1(s2,          ks  ,ke  , 0, 0, d)
    call bvalx2(s2,is-2,ie+2,ks  ,ke  , 0)
    call bvalx3(s2,is-2,ie+2,js-2,je+2, 0)
    call bvalx1(s3,          ks  ,ke+1, 0)
    call bvalx2(s3,is-2,ie+2,ks  ,ke+1, 0)
    call bvalx3(s3,is-2,ie+2,js-2,je+2, 1, 0)

    return
  end subroutine vtos




  subroutine stov
    use bval_module, only : bvalx1, bvalx2, bvalx3
    use field_module, only : d, v1, v2, v3, s1, s2, s3
    use grid_module, only : is, ie, js, je, ks, ke
    integer :: i,j,k

    do k = ks  , ke
    do i = is-1, ie+2
!CDIR SHORT_LOOP
    do j = js  , je
       v1(i,j,k) = s1(i,j,k) / (0.5d0 * (d(i-1,j,k) + d(i,j,k)))
    enddo
    enddo
    enddo
    do k = ks  , ke
    do i = is-2, ie+2
!CDIR SHORT_LOOP
    do j = js  , je+1
       v2(i,j,k) = s2(i,j,k) / (0.5d0 * (d(i,j-1,k) + d(i,j,k)))
    enddo
    enddo
    enddo
    do k = ks  , ke+1
    do i = is-2, ie+2
!CDIR SHORT_LOOP
    do j = js  , je
       v3(i,j,k) = s3(i,j,k) / (0.5d0 * (d(i,j,k-1) + d(i,j,k)))
    enddo
    enddo
    enddo

    call bvalx1(v1,          ks  ,ke  , 1)
    call bvalx1(v2,          ks  ,ke  , 0, 1)
    call bvalx1(v3,          ks  ,ke+1, 0)
    call bvalx2(v1,is-2,ie+3,ks  ,ke  , 0)
    call bvalx2(v2,is-2,ie+2,ks  ,ke  , 1)
    call bvalx2(v3,is-2,ie+2,ks  ,ke+1, 0)
    call bvalx3(v1,is-2,ie+3,js-2,je+2, 0)
    call bvalx3(v2,is-2,ie+2,js-2,je+3, 0)
    call bvalx3(v3,is-2,ie+2,js-2,je+2, 1, 0)

    return
  end subroutine stov




  subroutine tranx1
    use bval_module, only : bvalx1, bvalx2, bvalx3
    use field_module, only : d, e, v1, et, v2
    use flux_module, only : fga1
    use grid_module, only : is, ie, js, je, ks, ke, dx1a
    use interp_module, only : x1intzc
    use root_module, only : dt
    integer :: i, j, k
    real(8) :: dtwid(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: etwid(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: dflx(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: eflx(is-2:ie+3,js-2:je+3,ks-2:ke+3)
#ifdef PARTETOTK
    real(8) :: ettwid(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: etflx(is-2:ie+3,js-2:je+3,ks-2:ke+3)
#endif

    call x1intzc(d ,v1, dtwid)
    call x1intzc(e ,v1, etwid)
#ifdef PARTETOTK
    call bvalx1(et,          ks  ,ke  , 0, 0, v2, d)
    call bvalx2(et,is-2,ie+2,ks  ,ke  , 0)
    call bvalx3(et,is-2,ie+2,js-2,je+2, 0)
    
    call x1intzc(et,v1,ettwid)
#endif
    do k = ks, ke
    do i = is, ie+1
!CDIR SHORT_LOOP
    do j = js, je
       dflx(i,j,k) = v1(i,j,k) * dt * dtwid(i,j,k)
       eflx(i,j,k) = v1(i,j,k) * dt * etwid(i,j,k)
    enddo
    enddo
    enddo
    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       d(i,j,k) = d(i,j,k) - (dflx(i+1,j,k) - dflx(i,j,k)) / dx1a(i)
       e(i,j,k) = e(i,j,k) - (eflx(i+1,j,k) - eflx(i,j,k)) / dx1a(i)
    enddo
    enddo
    enddo
#ifdef PARTETOTK
    do k = ks, ke
    do i = is, ie+1
!CDIR SHORT_LOOP
    do j = js, je
       etflx(i,j,k) = v1(i,j,k) * dt * ettwid(i,j,k)
       fga1(i,j,k) = etflx(i,j,k) / dt
    enddo
    enddo
    enddo
    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       et(i,j,k) = et(i,j,k) - (etflx(i+1,j,k) - etflx(i,j,k)) / dx1a(i)
    enddo
    enddo
    enddo
#endif

    call bvalx1(d,          ks-2,ke+2, 0)
    call bvalx2(d,is-2,ie+2,ks  ,ke  , 0)
    call bvalx3(d,is-2,ie+2,js-2,je+2, 0)
    
    call bvalx1(e,          ks-2,ke+2, 0)
    call bvalx2(e,is-2,ie+2,ks  ,ke  , 0)
    call bvalx3(e,is-2,ie+2,js-2,je+2, 0)
    
    return
  end subroutine tranx1




  subroutine trans1
    use bval_module, only : bvalx1, bvalx2, bvalx3
    use field_module, only : v1, s1, s2, s3, tv1, d
    use grid_module, only : is, ie, js, je, ks, ke, dx1a, dx1b
    use interp_module, only : x1intzc, x1intfc
    use root_module, only : dt
    integer :: i, j, k
    real(8) :: stwid(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: sflx(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    
    do k = ks  , ke
    do i = is-1, ie
!CDIR SHORT_LOOP
    do j = js  , je
       tv1(i,j,k) = 0.5 * (v1(i,j,k) + v1(i+1,j,k))
    enddo
    enddo
    enddo
    call x1intfc(s1,tv1,stwid)
    do k = ks  , ke
    do i = is-1, ie
!CDIR SHORT_LOOP
    do j = js  , je
       sflx(i,j,k) = tv1(i,j,k) * dt * stwid(i,j,k)
    enddo
    enddo
    enddo
    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       s1(i,j,k) = s1(i,j,k) - (sflx(i,j,k) - sflx(i-1,j,k)) / dx1b(i)
    enddo
    enddo
    enddo

    do k = ks, ke
    do i = is, ie+1
!CDIR SHORT_LOOP
    do j = js, je
       tv1(i,j,k) = 0.5 * (v1(i,j-1,k) + v1(i,j,k))
    enddo
    enddo
    enddo
    call x1intzc(s2,tv1,stwid)
    do k = ks, ke
    do i = is, ie+1
!CDIR SHORT_LOOP
    do j = js, je
       sflx(i,j,k) = tv1(i,j,k) * dt * stwid(i,j,k)
    enddo
    enddo
    enddo
    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       s2(i,j,k) = s2(i,j,k) - (sflx(i+1,j,k) - sflx(i,j,k)) / dx1a(i)
    enddo
    enddo
    enddo

    do k = ks, ke+1
    do i = is, ie+1
!CDIR SHORT_LOOP
    do j = js, je
       tv1(i,j,k) = 0.5 * (v1(i,j,k-1) + v1(i,j,k))
    enddo
    enddo
    enddo
    call x1intzc(s3,tv1,stwid)
    do k = ks, ke+1
    do i = is, ie+1
!CDIR SHORT_LOOP
    do j = js, je
       sflx(i,j,k) = tv1(i,j,k) * dt * stwid(i,j,k)
    enddo
    enddo
    enddo
    do k = ks, ke+1
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       s3(i,j,k) = s3(i,j,k) - (sflx(i+1,j,k) - sflx(i,j,k)) / dx1a(i)
    enddo
    enddo
    enddo

    call bvalx1(s1,          ks  ,ke  , 1)
    call bvalx2(s1,is-2,ie+2,ks  ,ke  , 0)
    call bvalx3(s1,is-2,ie+2,js-2,je+2, 0)
    
    call bvalx1(s2,          ks  ,ke  , 0, 0, d)
    call bvalx2(s2,is-2,ie+2,ks  ,ke  , 0)
    call bvalx3(s2,is-2,ie+2,js-2,je+2, 0)
    
    call bvalx1(s3,          ks  ,ke+1, 0)
    call bvalx2(s3,is-2,ie+2,ks  ,ke+1, 0)
    call bvalx3(s3,is-2,ie+2,js-2,je+2, 1, 0)

    return
  end subroutine trans1




  subroutine tranx2
    use bval_module, only : bvalx1, bvalx2, bvalx3
    use field_module, only : d, e, v2, et
    use grid_module, only : is, ie, js, je, ks, ke, dx2a
    use interp_module, only : x2intzc
    use root_module, only : dt

    integer i,j,k
    real(8) :: dtwid(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: dflx(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: etwid(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: eflx(is-2:ie+3,js-2:je+3,ks-2:ke+3)
#ifdef PARTETOTK
    real(8) :: ettwid(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: etflx(is-2:ie+3,js-2:je+3,ks-2:ke+3)
#endif

    call x2intzc( d,v2, dtwid)
    call x2intzc( e,v2, etwid)
#ifdef PARTETOTK
    call bvalx1(et,          ks  ,ke  , 0, 0, v2, d)
    call bvalx2(et,is-2,ie+2,ks  ,ke  , 0)
    call bvalx3(et,is-2,ie+2,js-2,je+2, 0)
    call x2intzc(et,v2,ettwid)
#endif
    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je+1
       dflx(i,j,k) = v2(i,j,k) * dt * dtwid(i,j,k)
       eflx(i,j,k) = v2(i,j,k) * dt * etwid(i,j,k)
#ifdef PARTETOTK
       etflx(i,j,k) = v2(i,j,k) * dt * ettwid(i,j,k)
#endif
    enddo
    enddo
    enddo
    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       d(i,j,k) = d(i,j,k) - (dflx(i,j+1,k) - dflx(i,j,k)) / dx2a(j)
       e(i,j,k) = e(i,j,k) - (eflx(i,j+1,k) - eflx(i,j,k)) / dx2a(j)
#ifdef PARTETOTK
       et(i,j,k) = et(i,j,k) - (etflx(i,j+1,k) - etflx(i,j,k)) / dx2a(j)
#endif
    enddo
    enddo
    enddo

    call bvalx1(d,          ks-2,ke+2, 0)
    call bvalx2(d,is-2,ie+2,ks  ,ke  , 0)
    call bvalx3(d,is-2,ie+2,js-2,je+2, 0)
    
    call bvalx1(e,          ks-2,ke+2, 0)
    call bvalx2(e,is-2,ie+2,ks  ,ke  , 0)
    call bvalx3(e,is-2,ie+2,js-2,je+2, 0)
    
    return
  end subroutine tranx2




  subroutine trans2
    use bval_module, only : bvalx1, bvalx2, bvalx3
    use field_module, only : v2, s1, s2, s3, tv2, d
    use grid_module, only : is, ie, js, je, ks, ke, dx2a, dx2b
    use interp_module, only : x2intzc, x2intfc
    use root_module, only : dt
    integer :: i, j, k
    real(8) :: stwid(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: sflx(is-2:ie+3,js-2:je+3,ks-2:ke+3)

    do k = ks  , ke
    do i = is  , ie
!CDIR SHORT_LOOP
    do j = js-1, je
       tv2(i,j,k) = 0.5 * (v2(i,j,k) + v2(i,j+1,k))
    enddo
    enddo
    enddo
    call x2intfc(s2,tv2,stwid)
    do k = ks  , ke
    do i = is  , ie
!CDIR SHORT_LOOP
    do j = js-1, je
       sflx(i,j,k) = tv2(i,j,k) * dt * stwid(i,j,k)
    enddo
    enddo
    enddo
    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       s2(i,j,k) = s2(i,j,k) - (sflx(i,j,k) - sflx(i,j-1,k)) / dx2b(j)
    enddo
    enddo
    enddo
    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je+1
       tv2(i,j,k) = 0.5 * (v2(i-1,j,k) + v2(i,j,k))
    enddo
    enddo
    enddo
    call x2intzc(s1,tv2,stwid)
    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je+1
       sflx(i,j,k) = tv2(i,j,k) * dt * stwid(i,j,k)
    enddo
    enddo
    enddo
    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       s1(i,j,k) = s1(i,j,k) - (sflx(i,j+1,k) - sflx(i,j,k)) / dx2a(j)
    enddo
    enddo
    enddo
    do k = ks, ke+1
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je+1
       tv2(i,j,k) = 0.5 * (v2(i,j,k-1) + v2(i,j,k))
    enddo
    enddo
    enddo
    call x2intzc(s3,tv2,stwid)
    do k = ks, ke+1
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je+1
       sflx(i,j,k) = tv2(i,j,k) * dt * stwid(i,j,k)
    enddo
    enddo
    enddo
    do k = ks, ke+1
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       s3(i,j,k) = s3(i,j,k) - (sflx(i,j+1,k) - sflx(i,j,k)) / dx2a(j)
    enddo
    enddo
    enddo

    call bvalx1(s1,          ks  ,ke  , 1)
    call bvalx2(s1,is-2,ie+2,ks  ,ke  , 0)
    call bvalx3(s1,is-2,ie+2,js-2,je+2, 0)
    
    call bvalx1(s2,          ks  ,ke  , 0, 0, d)
    call bvalx2(s2,is-2,ie+2,ks  ,ke  , 0)
    call bvalx3(s2,is-2,ie+2,js-2,je+2, 0)
    
    call bvalx1(s3,          ks  ,ke+1, 0)
    call bvalx2(s3,is-2,ie+2,ks  ,ke+1, 0)
    call bvalx3(s3,is-2,ie+2,js-2,je+2, 1, 0)

    return
  end subroutine trans2




  subroutine tranx3
    use bval_module, only : bvalx1, bvalx2, bvalx3
    use field_module, only : d, e, v3, et, v2
    use flux_module, only : fga3
    use grid_module, only : is, ie, js, je, ks, ke, dx3a
    use interp_module, only : x3intzc
    use root_module, only : dt
    integer :: i, j, k
    real(8) :: dtwid(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: dflx(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: etwid(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: eflx(is-2:ie+3,js-2:je+3,ks-2:ke+3)
#ifdef PARTETOTK
    real(8) :: ettwid(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: etflx(is-2:ie+3,js-2:je+3,ks-2:ke+3)
#endif

    call x3intzc( d,v3, dtwid)
    call x3intzc( e,v3, etwid)
#ifdef PARTETOTK
    call bvalx1(et,          ks  ,ke  , 0, 0, v2, d)
    call bvalx2(et,is-2,ie+2,ks  ,ke  , 0)
    call bvalx3(et,is-2,ie+2,js-2,je+2, 0)
    call x3intzc(et,v3,ettwid)
#endif
    do k = ks, ke+1
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       dflx(i,j,k) = v3(i,j,k) * dt * dtwid(i,j,k)
       eflx(i,j,k) = v3(i,j,k) * dt * etwid(i,j,k)
#ifdef PARTETOTK
       etflx(i,j,k) = v3(i,j,k) * dt * ettwid(i,j,k)
       fga3(i,j,k) = etflx(i,j,k) / dt
#endif
    enddo
    enddo
    enddo
    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       d(i,j,k) = d(i,j,k) - (dflx(i,j,k+1) - dflx(i,j,k)) / dx3a(k)
       e(i,j,k) = e(i,j,k) - (eflx(i,j,k+1) - eflx(i,j,k)) / dx3a(k)
#ifdef PARTETOTK
       et(i,j,k) = et(i,j,k) - (etflx(i,j,k+1) - etflx(i,j,k)) / dx3a(k)
#endif
    enddo
    enddo
    enddo

    call bvalx1(d,          ks-2,ke+2, 0)
    call bvalx2(d,is-2,ie+2,ks  ,ke  , 0)
    call bvalx3(d,is-2,ie+2,js-2,je+2, 0)
    
    call bvalx1(e,          ks-2,ke+2, 0)
    call bvalx2(e,is-2,ie+2,ks  ,ke  , 0)
    call bvalx3(e,is-2,ie+2,js-2,je+2, 0)
    
    return
  end subroutine tranx3




  subroutine trans3
    use bval_module, only : bvalx1, bvalx2, bvalx3
    use field_module, only : v3, s1, s2, s3, tv3, d
    use grid_module, only : is, ie, js, je, ks, ke, dx3a, dx3b
    use interp_module, only : x3intzc, x3intfc
    use root_module, only : dt
    integer :: i, j, k
    real(8) :: stwid(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: sflx(is-2:ie+3,js-2:je+3,ks-2:ke+3)

    do k = ks-2, ke+2
    do i = is  , ie
!CDIR SHORT_LOOP
    do j = js  , je
       tv3(i,j,k) = 0.5 * (v3(i,j,k) + v3(i,j,k+1))
    enddo
    enddo
    enddo
    call x3intfc(s3,tv3,stwid)
    do k = ks-1, ke+1
    do i = is  , ie
!CDIR SHORT_LOOP
    do j = js, je
       sflx(i,j,k) = tv3(i,j,k) * dt * stwid(i,j,k)
    enddo
    enddo
    enddo
    do k = ks, ke+1
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       s3(i,j,k) = s3(i,j,k) - (sflx(i,j,k) - sflx(i,j,k-1)) / dx3b(k)
    enddo
    enddo
    enddo
    do k = ks-2, ke+2
    do i = is  , ie
!CDIR SHORT_LOOP
    do j = js, je
       tv3(i,j,k) = 0.5 * (v3(i-1,j,k) + v3(i,j,k))
    enddo
    enddo
    enddo
    call x3intzc(s1,tv3,stwid)
    do k = ks, ke+1
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       sflx(i,j,k) = tv3(i,j,k) * dt * stwid(i,j,k)
    enddo
    enddo
    enddo
    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       s1(i,j,k) = s1(i,j,k) - (sflx(i,j,k+1) - sflx(i,j,k)) / dx3a(k)
    enddo
    enddo
    enddo
    do k = ks-2, ke+2
    do i = is  , ie
!CDIR SHORT_LOOP
    do j = js  , je
       tv3(i,j,k) = 0.5 * (v3(i,j-1,k) + v3(i,j,k))
    enddo
    enddo
    enddo
    call x3intzc(s2,tv3,stwid)
    do k = ks, ke+1
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       sflx(i,j,k) = tv3(i,j,k) * dt * stwid(i,j,k)
    enddo
    enddo
    enddo
    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       s2(i,j,k) = s2(i,j,k) - (sflx(i,j,k+1) - sflx(i,j,k)) / dx3a(k)
    enddo
    enddo
    enddo
    
    call bvalx1(s1,          ks  ,ke  , 1)
    call bvalx2(s1,is-2,ie+2,ks  ,ke  , 0)
    call bvalx3(s1,is-2,ie+2,js-2,je+2, 0)
    
    call bvalx1(s2,          ks  ,ke  , 0, 0, d)
    call bvalx2(s2,is-2,ie+2,ks  ,ke  , 0)
    call bvalx3(s2,is-2,ie+2,js-2,je+2, 0)
    
    call bvalx1(s3,          ks  ,ke+1, 0)
    call bvalx2(s3,is-2,ie+2,ks  ,ke+1, 0)
    call bvalx3(s3,is-2,ie+2,js-2,je+2, 1, 0)
    
  end subroutine trans3




end module transport_module
