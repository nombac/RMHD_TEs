
module source_module

  implicit none

  private
  public :: source

contains




  subroutine source
    use radiation_module, only : radiation_source
    use shear_module, only : shear_source

    call shear_source
    call radiation_source
    call hydro_source

    return
  end subroutine source



  
  subroutine hydro_source
    use bval_module, only : bvalx1, bvalx2, bvalx3
    use convert_module, only : convert_et2e, convert_e2et
    use field_module, only : d, v1, v2, v3, e, et, p, qkin
    use flux_module, only : sdek, fgw1, fgw3, divv
    use grid_module, only : is, ie, js, je, ks, ke, dz, &
         dx1a, dx1b, dx2a, dx2b, dx3a, dx3b
    use interp_module, only : x1intzc, x2intzc, x3intzc
    use param_module, only : GAMMA
    use root_module, only : dt
    integer :: i, j, k
    real(8) :: qa, df
#ifdef PARTETOTK
    real(8) :: qq1flx(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: qq2flx(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: qq3flx(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: qq1twid(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: qq2twid(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: qq3twid(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: tma(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: tmb(is-2:ie+3,js-2:je+3,ks-2:ke+3)
#endif
    real(8) :: div1, div2, div3, del
!-----------------------------------------------------------------------
#ifdef PARTETOTK
    call convert_e2et(0)
    !
    ! d(etot)/dt = -div(pv)
    !
    do k=ks-2,ke+2
    do i=is-2,ie+2
!CDIR SHORT_LOOP
    do j=js-2,je+2
       p(i,j,k)=(gamma-1.)*e(i,j,k)
    enddo
    enddo
    enddo
    call x1intzc(p,v1,qq1twid)
    do k=ks,ke
    do i=is,ie+1
!CDIR SHORT_LOOP
    do j=js,je
       qq1flx(i,j,k) = v1(i,j,k)*dt*qq1twid(i,j,k)
       fgw1 (i,j,k)=        +qq1flx(i,j,k)/dt
    enddo
    enddo
    enddo
    do k=ks,ke
    do i=is,ie
!CDIR SHORT_LOOP
    do j=js,je
       et(i,j,k) = et(i,j,k) - (qq1flx(i+1,j,k)-qq1flx(i,j,k))/dx1a(i)
    enddo
    enddo
    enddo
    call x2intzc(p,v2,qq2twid)
    do k=ks,ke
    do i=is,ie
!CDIR SHORT_LOOP
    do j=js,je+1
       qq2flx(i,j,k) = v2(i,j,k)*dt*qq2twid(i,j,k)
    enddo
    enddo
    enddo
    do k=ks,ke
    do i=is,ie
!CDIR SHORT_LOOP
    do j=js,je
       et(i,j,k) = et(i,j,k) - (qq2flx(i,j+1,k)-qq2flx(i,j,k))/dx2a(j)
    enddo
    enddo
    enddo
    call x3intzc(p,v3,qq3twid)
    do k=ks,ke+1
    do i=is,ie
!CDIR SHORT_LOOP
    do j=js,je
       qq3flx(i,j,k) = v3(i,j,k)*dt*qq3twid(i,j,k)
       fgw3 (i,j,k)=        +qq3flx(i,j,k)/dt
    enddo
    enddo
    enddo
    do k=ks,ke
    do i=is,ie
!CDIR SHORT_LOOP
    do j=js,je
       et(i,j,k) = et(i,j,k) - (qq3flx(i,j,k+1)-qq3flx(i,j,k))/dz
    enddo
    enddo
    enddo
#endif

    !
    ! de/dt = -p(div v)
    !
    qa = 0.5*dt*(gamma-1.0)
    divv(is:ie,js:je,ks:ke) = e(is:ie,js:je,ks:ke)
    do k=ks,ke
    do i=is,ie
!CDIR SHORT_LOOP
    do j=js,je
       div1 = (v1(i+1,j,k)-v1(i,j,k))/dx1a(i)
       div2 = (v2(i,j+1,k)-v2(i,j,k))/dx2a(j)
       div3 = (v3(i,j,k+1)-v3(i,j,k))/dx3a(k)
       del  = qa*(div1 + div2 + div3)
       e(i,j,k) = (1.0 - del)/(1.0 + del)*e(i,j,k)
    enddo
    enddo
    enddo
    divv(is:ie,js:je,ks:ke) = e(is:ie,js:je,ks:ke) - divv(is:ie,js:je,ks:ke)

    call bvalx1(e,          ks-2,ke+2, 0)
    call bvalx2(e,is-2,ie+2,ks  ,ke  , 0)
    call bvalx3(e,is-2,ie+2,js-2,je+2, 0)

    !
    ! dv/dt = -(grad p)/d
    !
    qa = gamma-1.0
    do k=ks,ke
    do i=is,ie
!CDIR SHORT_LOOP
    do j=js,je
       df=0.5*(d(i-1,j,k)+d(i,j,k))
       v1(i,j,k)=v1(i,j,k)+dt*(qa*(e(i-1,j,k)-e(i,j,k))/dx1b(i)/df)
    enddo
    enddo
    enddo
    call bvalx1(v1,          ks  ,ke  , 1)
    call bvalx2(v1,is-2,ie+3,ks  ,ke  , 0)
    call bvalx3(v1,is-2,ie+3,js-2,je+2, 0)
    do k=ks,ke
    do i=is,ie
!CDIR SHORT_LOOP
    do j=js,je
       df=0.5*(d(i,j-1,k)+d(i,j,k))
       v2(i,j,k)=v2(i,j,k)+dt*(qa*(e(i,j-1,k)-e(i,j,k))/dx2b(j)/df)
    enddo
    enddo
    enddo
    call bvalx1(v2,          ks  ,ke  , 0, 1)
    call bvalx2(v2,is-2,ie+2,ks  ,ke  , 1)
    call bvalx3(v2,is-2,ie+2,js-2,je+3, 0)
    do k=ks,ke+1
    do i=is,ie
!CDIR SHORT_LOOP
    do j=js,je
       df=0.5*(d(i,j,k-1)+d(i,j,k))
       v3(i,j,k)=v3(i,j,k)+dt*(qa*(e(i,j,k-1)-e(i,j,k))/dx3b(k)/df)
    enddo
    enddo
    enddo
    call bvalx1(v3,          ks  ,ke+1, 0)
    call bvalx2(v3,is-2,ie+2,ks  ,ke+1, 0)
    call bvalx3(v3,is-2,ie+2,js-2,je+2, 1, 0)

#ifdef PARTETOTK
    !
    ! e = (etot) - ekin
    !
    call convert_et2e(0,tma,tmb)
    do k=ks,ke
    do i=is,ie
!CDIR SHORT_LOOP
    do j=js,je
       qkin (i,j,k)=         tma(i,j,k)/dt
       sdek (i,j,k)=         tmb(i,j,k)/dt
    enddo
    enddo
    enddo
#endif

    return
  end subroutine hydro_source




end module source_module
