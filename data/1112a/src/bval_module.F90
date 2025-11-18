
module bval_module

#ifdef MPI
  use mpi_module
#endif

  implicit none

  private
  public :: bvalx1, bvalx2, bvalx3
  public :: tcons, epsi, epso, deltay, jplus
  
  real(8) :: epsi, epso
  real(8) :: deltay
  integer :: jplus

contains




  subroutine tcons
    use grid_module, only : lx0, ly0, dy
    use param_module, only : omega
    use root_module, only : time, dt

    deltay = 1.5d0 * omega * lx0 * (time + dt)
    deltay = dmod(deltay, ly0)
    jplus = int(deltay / dy)
    epsi = dmod(deltay, dy)
    epso = -epsi
    
    return
  end subroutine tcons




  subroutine bvalx1(q, kb, kf, isw, isw2, aa, bb)

    use grid_module, only : is, ie, js, je, js0, je0, lx0, dx2a, nx2z
    use param_module, only : omega

    real(8), intent(inout) :: q(:,:,:)
    integer, intent(in), optional :: isw2
    real(8), intent(in), optional :: aa(:,:,:)
    real(8), intent(in), optional :: bb(:,:,:)
    integer, intent(in) :: kb, kf, isw

    integer :: j, k, jremap
    real(8) :: eps, qa
    real(8) :: q1(js0-2:je0+3,kb:kf),q2(js0-2:je0+3,kb:kf),q3(js0-2:je0+3,kb:kf)
    real(8) :: a1(js0-2:je0+3,kb:kf),a2(js0-2:je0+3,kb:kf),a3(js0-2:je+3,kb:kf)
    real(8) :: b1(js0-2:je0+3,kb:kf),b2(js0-2:je0+3,kb:kf),b3(js0-2:je+3,kb:kf)
    real(8) :: delq1(js-2:je+3), delq2(js-2:je+3), delq3(js-2:je+3)
    real(8) ::    f1(js-2:je+3),    f2(js-2:je+3),    f3(js-2:je+3)
    real(8) :: delq12, delq22, delq32, dq1, dq2, dq3
#ifdef MPI
    real(8) :: qi(is:is+1+isw, js0:je0, kb:kf), qo(ie-1:ie, js0:je0, kb:kf)
    real(8) :: ai(is:is+1+isw, js0:je0, kb:kf), ao(ie-1:ie, js0:je0, kb:kf)
    real(8) :: bi(is:is+1+isw, js0:je0, kb:kf), bo(ie-1:ie, js0:je0, kb:kf)
#endif
    
    qa = 1.5 * omega * lx0

#ifdef MPI
    call gatherx2( q, qi, qo, is, is+1+isw, ie-1, ie, kb, kf)
    if(present(aa)) then
       call gatherx2(aa, ai, ao, is, is+1+isw, ie-1, ie, kb, kf)
    endif
    if(present(bb)) then
       call gatherx2(bb, bi, bo, is, is+1+isw, ie-1, ie, kb, kf)
    endif

    if(myrki == 0) then
#endif
    eps = epsi / dx2a(js0)
    if(present(aa)) then
       do k = kb  , kf
!CDIR SHORT_LOOP
          do j = js0+1, je0
#ifdef MPI
             a1(j,k) = 0.5 * (ao(ie  ,j,k) + ao(ie  ,j-1,k))
             a2(j,k) = 0.5 * (ao(ie-1,j,k) + ao(ie-1,j-1,k))
#else
             a1(j,k) = 0.5 * (aa(ie  ,j,k) + aa(ie  ,j-1,k))
             a2(j,k) = 0.5 * (aa(ie-1,j,k) + aa(ie-1,j-1,k))
#endif
          end do
#ifdef MPI
          a1(js0,k) = 0.5 * (ao(ie  ,js0,k) + ao(ie  ,je0,k))
          a2(js0,k) = 0.5 * (ao(ie-1,js0,k) + ao(ie-1,je0,k))
#else
          a1(js0,k) = 0.5 * (aa(ie  ,js0,k) + aa(ie  ,je0,k))
          a2(js0,k) = 0.5 * (aa(ie-1,js0,k) + aa(ie-1,je0,k))
#endif
       end do
    end if
    if(present(bb)) then
       do k = kb  , kf
!CDIR SHORT_LOOP
       do j = js0  , je0
#ifdef MPI
          b1(j,k) = bo(ie  ,j,k)
          b2(j,k) = bo(ie-1,j,k)
#else
          b1(j,k) = bb(ie  ,j,k)
          b2(j,k) = bb(ie-1,j,k)
#endif
       end do
       end do
    end if
    do k = kb, kf
!CDIR SHORT_LOOP
    do j = js0, je0
       jremap = j - jplus
       if(jremap <= js0-1) jremap = jremap + nx2z
#ifdef MPI
       q1(j,k) =qo(ie  ,jremap,k)
       q2(j,k) =qo(ie-1,jremap,k)
#else
       q1(j,k) = q(ie  ,jremap,k)
       q2(j,k) = q(ie-1,jremap,k)
#endif
       if(present(bb)) then
          q1(j,k) = q1(j,k) + (.5d0 * qa**2 + qa * a1(jremap,k)) * b1(jremap,k)
          q2(j,k) = q2(j,k) + (.5d0 * qa**2 + qa * a2(jremap,k)) * b2(jremap,k)
       else if(present(aa)) then
          q1(j,k) = q1(j,k) + qa * a1(jremap,k)
          q2(j,k) = q2(j,k) + qa * a2(jremap,k)
       else if(present(isw2)) then
          q1(j,k) = q1(j,k) + qa
          q2(j,k) = q2(j,k) + qa
       end if
    end do
    end do
    do k = kb, kf
       q1(js0-1,k) = q1(je0  ,k)
       q1(js0-2,k) = q1(je0-1,k)
       q1(je0+1,k) = q1(js0  ,k)
       q2(js0-1,k) = q2(je0  ,k)
       q2(js0-2,k) = q2(je0-1,k)
       q2(je0+1,k) = q2(js0  ,k)
    end do
    do k = kb, kf
!CDIR SHORT_LOOP
       do j = js-1, je+1
          delq1(j) = q1(j,k) - q1(j-1,k)
          delq2(j) = q2(j,k) - q2(j-1,k)
       end do
!CDIR SHORT_LOOP
       do j = js, je+1
          delq12 = delq1(j-1) * delq1(j)
          delq22 = delq2(j-1) * delq2(j)
          dq1 = 0
          dq2 = 0
          if(delq12 > 0) dq1 = delq12 / (delq1(j-1) + delq1(j))
          if(delq22 > 0) dq2 = delq22 / (delq2(j-1) + delq2(j))
          f1(j) = eps * (q1(j-1,k) + (1d0 - eps) * dq1)
          f2(j) = eps * (q2(j-1,k) + (1d0 - eps) * dq2)
       end do
!CDIR SHORT_LOOP
       do j = js, je
          q(is-1,j,k) = q1(j,k) - (f1(j+1)-f1(j))
          q(is-2,j,k) = q2(j,k) - (f2(j+1)-f2(j))
       end do
    end do
#ifdef MPI
    end if
    if(myrki == iprc-1) then
#endif
    eps = epso / dx2a(js0)
    if(present(aa)) then
       do k = kb ,kf 
!CDIR SHORT_LOOP
          do j = js0+1, je0
#ifdef MPI
             a1(j,k) = 0.5 * (ai(is  ,j,k) + ai(is  ,j-1,k))
             a2(j,k) = 0.5 * (ai(is+1,j,k) + ai(is+1,j-1,k))
             a3(j,k) = 0.5 * (ai(is+2,j,k) + ai(is+2,j-1,k))
#else
             a1(j,k) = 0.5 * (aa(is  ,j,k) + aa(is  ,j-1,k))
             a2(j,k) = 0.5 * (aa(is+1,j,k) + aa(is+1,j-1,k))
             a3(j,k) = 0.5 * (aa(is+2,j,k) + aa(is+2,j-1,k))
#endif
          end do
#ifdef MPI
          a1(js0,k) = 0.5 * (ai(is  ,js0,k) + ai(is  ,je0,k))
          a2(js0,k) = 0.5 * (ai(is+1,js0,k) + ai(is+1,je0,k))
          a3(js0,k) = 0.5 * (ai(is+2,js0,k) + ai(is+2,je0,k))
#else
          a1(js0,k) = 0.5 * (aa(is  ,js0,k) + aa(is  ,je0,k))
          a2(js0,k) = 0.5 * (aa(is+1,js0,k) + aa(is+1,je0,k))
          a3(js0,k) = 0.5 * (aa(is+2,js0,k) + aa(is+2,je0,k))
#endif
       end do
    end if
    if(present(bb)) then
       do k = kb  ,kf 
!CDIR SHORT_LOOP
       do j = js0  ,je0
#ifdef MPI
          b1(j,k) = bi(is  ,j,k)
          b2(j,k) = bi(is+1,j,k)
          b3(j,k) = bi(is+2,j,k)
#else
          b1(j,k) = bb(is  ,j,k)
          b2(j,k) = bb(is+1,j,k)
          b3(j,k) = bb(is+2,j,k)
#endif
       end do
       end do
    end if
    do k = kb, kf
!CDIR SHORT_LOOP
    do j = js0, je0
       jremap = j + jplus
       if(jremap >= je0+1) jremap = jremap - nx2z
#ifdef MPI
       q1(j,k) = qi(is  ,jremap,k)
       q2(j,k) = qi(is+1,jremap,k)
       q3(j,k) = qi(is+2,jremap,k)
#else
       q1(j,k) =  q(is  ,jremap,k)
       q2(j,k) =  q(is+1,jremap,k)
       q3(j,k) =  q(is+2,jremap,k)
#endif
       if(present(bb)) then
          q1(j,k) = q1(j,k) + (.5d0 * qa**2 - qa * a1(jremap,k)) * b1(jremap,k) 
          q2(j,k) = q2(j,k) + (.5d0 * qa**2 - qa * a2(jremap,k)) * b2(jremap,k) 
          q3(j,k) = q3(j,k) + (.5d0 * qa**2 - qa * a3(jremap,k)) * b3(jremap,k) 
       else if(present(aa)) then
          q1(j,k) = q1(j,k) - qa * a1(jremap,k)
          q2(j,k) = q2(j,k) - qa * a2(jremap,k)
          q3(j,k) = q3(j,k) - qa * a3(jremap,k)
       else if(present(isw2)) then
          q1(j,k) = q1(j,k) - qa
          q2(j,k) = q2(j,k) - qa
          q3(j,k) = q3(j,k) - qa
       end if
    end do
    end do
    do k = kb, kf
       q1(js0-1,k) = q1(je0  ,k)
       q1(je0+1,k) = q1(js0  ,k)
       q1(je0+2,k) = q1(js0+1,k)
       q2(js0-1,k) = q2(je0  ,k)
       q2(je0+1,k) = q2(js0  ,k)
       q2(je0+2,k) = q2(js0+1,k)
       q3(js0-1,k) = q3(je0  ,k)
       q3(je0+1,k) = q3(js0  ,k)
       q3(je0+2,k) = q3(js0+1,k)
    end do
    do k = kb, kf
!CDIR SHORT_LOOP
       do j = js, je+2
          delq1(j) = q1(j,k) - q1(j-1,k)
          delq2(j) = q2(j,k) - q2(j-1,k)
          delq3(j) = q3(j,k) - q3(j-1,k)
       end do
!CDIR SHORT_LOOP
       do j = js, je+1
          delq12 = delq1(j) * delq1(j+1)
          delq22 = delq2(j) * delq2(j+1)
          delq32 = delq3(j) * delq3(j+1)
          dq1 = 0
          dq2 = 0
          dq3 = 0
          if(delq12 > 0) dq1 = delq12 / (delq1(j) + delq1(j+1))
          if(delq22 > 0) dq2 = delq22 / (delq2(j) + delq2(j+1))
          if(delq32 > 0) dq3 = delq32 / (delq3(j) + delq3(j+1))
          f1(j) = eps*(q1(j,k) - (1d0 + eps) * dq1)
          f2(j) = eps*(q2(j,k) - (1d0 + eps) * dq2)
          f3(j) = eps*(q3(j,k) - (1d0 + eps) * dq3)
       end do
!CDIR SHORT_LOOP
       if(isw == 0) then
          do j = js, je
             q(ie+1,j,k) = q1(j,k) - (f1(j+1) - f1(j))
             q(ie+2,j,k) = q2(j,k) - (f2(j+1) - f2(j))
          end do
       else if(isw == 1) then
          do j = js, je
             q(ie+1,j,k) = q1(j,k) - (f1(j+1) - f1(j))
             q(ie+2,j,k) = q2(j,k) - (f2(j+1) - f2(j))
             q(ie+3,j,k) = q3(j,k) - (f3(j+1) - f3(j))
          end do
       else if(isw == 2) then
          do j = js, je
             q(ie+2,j,k) = q2(j,k) - (f2(j+1) - f2(j))
             q(ie+3,j,k) = q3(j,k) - (f3(j+1) - f3(j))
          end do
       end if
    end do
#ifdef MPI
    end if
#endif

    return

  end subroutine bvalx1




  subroutine bvalx2(q, ib, if, kb, kf, isw)

    use grid_module, only : js, je

    real(8), intent(inout) :: q(:,:,:)
    integer, intent(in) :: ib
    integer, intent(in) :: if
    integer, intent(in) :: kb
    integer, intent(in) :: kf
    integer, intent(in) :: isw

    integer :: i, k

    do k = kb, kf
!CDIR SHORT_LOOP
    do i = ib, if
       q(i,js-1,k) = q(i,je  ,k)
       q(i,js-2,k) = q(i,je-1,k)
       q(i,je+1,k) = q(i,js  ,k)
       q(i,je+2,k) = q(i,js+1,k)
       if(isw == 1) then
       q(i,je+3,k) = q(i,js+2,k)
       end if
    end do
    end do

#ifdef MPI
    call shiftx2(q, ib, if, je+1, je+2+isw, js-2, js-1, kb, kf)
#endif

    return

  end subroutine bvalx2




  subroutine bvalx3(q, ib, if, jb, jf, isw, isw2)

    use grid_module, only : ks, ke
    use field_module, only : erikb, erokb

    real(8), intent(inout) :: q(:,:,:)
    integer, intent(in) :: ib
    integer, intent(in) :: if
    integer, intent(in) :: jb
    integer, intent(in) :: jf
    integer, intent(in) :: isw
    integer, intent(in), optional :: isw2

    integer :: i, j

#ifdef MPI
     if(myrkk == 0) then
#endif
    do i = ib, if
!CDIR SHORT_LOOP
    do j = jb, jf
#ifdef FLOORV3
       if(present(isw2)) then
          if(isw2 == 0) then
             if(q(i,j,ks  ) > 0) then
                q(i,j,ks  ) = 0
             end if
          end if
       end if
#endif
       q(i,j,ks-1) = q(i,j,ks  )
       q(i,j,ks-2) = q(i,j,ks  )
       if(present(isw2)) then
          if(isw2 == 1) then
             q(i,j,ks-1) = q(i,j,ks  ) * erikb(i,j)
             q(i,j,ks-2) = q(i,j,ks  ) * erikb(i,j)
          end if
       end if
    end do
    end do
#ifdef MPI
     end if
     if(myrkk == kprc-1) then
#endif
    do i = ib, if
!CDIR SHORT_LOOP
    do j = jb, jf
#ifdef FLOORV3
       if(present(isw2)) then
          if(isw2 == 0) then
             if(q(i,j,ke+isw) < 0) then
                q(i,j,ke+isw) = 0
             end if
          end if
       end if
#endif
       q(i,j,ke+1+isw) = q(i,j,ke+isw)
       q(i,j,ke+2+isw) = q(i,j,ke+isw)
       if(present(isw2)) then
          if(isw2 == 1) then
             q(i,j,ke+1+isw) = q(i,j,ke+isw) * erokb(i,j)
             q(i,j,ke+2+isw) = q(i,j,ke+isw) * erokb(i,j)
          end if
       end if
    end do
    end do
#ifdef MPI
    end if
    call shiftx3(q,ib,if,jb,jf,ke+1+isw,ke+2+isw,ks-2,ks-1)
#endif
     
    return

  end subroutine bvalx3




#ifdef MPI
  subroutine shiftx3(w, ib, if, jb, jf, kb1, kf1, kb2, kf2)

    use grid_module, only : ks, ke

    real(8), intent(inout) :: w(:,:,:)
    integer, intent(in) :: ib
    integer, intent(in) :: if
    integer, intent(in) :: jb
    integer, intent(in) :: jf
    integer, intent(in) :: kb1
    integer, intent(in) :: kf1
    integer, intent(in) :: kb2
    integer, intent(in) :: kf2

    integer :: jl1, jl2, kds
    integer :: isd1, isd2, irv1, irv2
    integer :: i, j, k, l
    integer :: ist(mpi_status_size)
    real(8) :: wps(0:(if - ib + 1) * (jf - jb + 1) * (kf2 - kb2 + 1) - 1)
    real(8) :: wms(0:(if - ib + 1) * (jf - jb + 1) * (kf1 - kb1 + 1) - 1)
    real(8) :: wpr(0:(if - ib + 1) * (jf - jb + 1) * (kf1 - kb1 + 1) - 1)
    real(8) :: wmr(0:(if - ib + 1) * (jf - jb + 1) * (kf2 - kb2 + 1) - 1)
    integer :: im, imjm
     
    kds = ke - ks + 1
    im =   (if - ib + 1)
    imjm = (jf - jb + 1) * im
    jl1 = imjm * (kf1 - kb1 + 1)
    jl2 = imjm * (kf2 - kb2 + 1)

    ! send to k+ (1)
    if(myrkk /= kprc-1) then
       l = 0
       do k = kb2 + kds, kf2 + kds
       do j = jb, jf
!CDIR SHORT_LOOP
       do i = ib, if
          wps(l) = w(i,j,k)
          l = l + 1
       end do
       end do
       end do
    end if
    ! send to k- (2)
    if(myrkk /= 0) then
       l = 0
       do k = kb1 - kds, kf1 - kds
       do j = jb, jf
!CDIR SHORT_LOOP
       do i = ib, if
          wms(l) = w(i,j,k)
          l = l + 1
       end do
       end do
       end do
    end if
    call mpi_isend(wps,jl2,mpi_double_precision,ku,1,mpi_comm_world,isd2,ier)
    call mpi_irecv(wmr,jl2,mpi_double_precision,kd,1,mpi_comm_world,irv2,ier)
    call mpi_wait(isd2,ist,ier)
    call mpi_wait(irv2,ist,ier)
    call mpi_isend(wms,jl1,mpi_double_precision,kd,1,mpi_comm_world,isd1,ier)
    call mpi_irecv(wpr,jl1,mpi_double_precision,ku,1,mpi_comm_world,irv1,ier)
    call mpi_wait(isd1,ist,ier)
    call mpi_wait(irv1,ist,ier)
    ! receive from k+ (1)
    if(myrkk /= kprc-1) then
       l = 0
       do k = kb1, kf1
       do j = jb, jf
!CDIR SHORT_LOOP
       do i = ib, if
          w(i,j,k) = wpr(l)
          l = l + 1
       end do
       end do
       end do
    end if
    ! receive from k- (2)
    if(myrkk /= 0) then
       l = 0
       do k = kb2, kf2
       do j = jb, jf
!CDIR SHORT_LOOP
       do i = ib, if
          w(i,j,k) = wmr(l)
          l = l + 1
       end do
       end do
       end do
    end if
     
    return

  end subroutine shiftx3
#endif




#ifdef MPI
  subroutine shiftx2(w, ib, if, jb1, jf1, jb2, jf2, kb, kf)

    use grid_module, only : js, je

    real(8), intent(inout) :: w(:,:,:)
    integer, intent(in) :: ib
    integer, intent(in) :: if
    integer, intent(in) :: jb1
    integer, intent(in) :: jf1
    integer, intent(in) :: jb2
    integer, intent(in) :: jf2
    integer, intent(in) :: kb
    integer, intent(in) :: kf

    integer :: jl1, jl2, jds
    integer :: isd1, isd2, irv1, irv2
    integer :: i, j, k, l
    integer :: ist(mpi_status_size)
    real(8) :: wps(0:(if - ib + 1) * (jf2 - jb2 + 1) * (kf - kb + 1) - 1)
    real(8) :: wms(0:(if - ib + 1) * (jf1 - jb1 + 1) * (kf - kb + 1) - 1)
    real(8) :: wpr(0:(if - ib + 1) * (jf1 - jb1 + 1) * (kf - kb + 1) - 1)
    real(8) :: wmr(0:(if - ib + 1) * (jf2 - jb2 + 1) * (kf - kb + 1) - 1)
    integer :: im, imjm

    jds = je - js + 1
    im =   (if - ib + 1)
    imjm = (kf - kb + 1) * im
    jl1 =  imjm * (jf1 - jb1 + 1)
    jl2 =  imjm * (jf2 - jb2 + 1)

    ! send to j+ (1)
    l = 0
    do k = kb, kf
    do j = jb2 + jds, jf2 + jds
!CDIR SHORT_LOOP
    do i = ib, if
       wps(l) = w(i,j,k)
       l = l + 1
    end do
    end do
    end do
    ! send to j- (2)
    l = 0
    do k = kb, kf
    do j = jb1 - jds, jf1 - jds
!CDIR SHORT_LOOP
    do i = ib, if
       wms(l) = w(i,j,k)
       l = l + 1
    end do
    end do
    end do
    if(jprc /= 1) then
       call mpi_isend(wps,jl2,mpi_double_precision,ju,1,mpi_comm_world,isd2,ier)
       call mpi_irecv(wmr,jl2,mpi_double_precision,jd,1,mpi_comm_world,irv2,ier)
       call mpi_wait(isd2,ist,ier)
       call mpi_wait(irv2,ist,ier)
       call mpi_isend(wms,jl1,mpi_double_precision,jd,1,mpi_comm_world,isd1,ier)
       call mpi_irecv(wpr,jl1,mpi_double_precision,ju,1,mpi_comm_world,irv1,ier)
       call mpi_wait(isd1,ist,ier)
       call mpi_wait(irv1,ist,ier)
    else
       wpr(0:jl1-1) = wms(0:jl1-1)
       wmr(0:jl2-1) = wps(0:jl2-1)
    end if
    ! receive from j+ (1)
    l = 0
    do k = kb, kf
    do j = jb1, jf1
!CDIR SHORT_LOOP
    do i = ib, if
       w(i,j,k) = wpr(l)
       l = l + 1
    end do
    end do
    end do
    ! receive from j- (2)
    l = 0
    do k = kb, kf
    do j = jb2, jf2
!CDIR SHORT_LOOP
    do i = ib, if
       w(i,j,k) = wmr(l)
       l = l + 1
    end do
    end do
    end do
     
    return
    
  end subroutine shiftx2




  subroutine gatherx2(w, wi, wo, ibi, ifi, ibo, ifo, kb, kf)

    use grid_module, only : js, je, js0, je0

    real(8), intent(in) :: w(:,:,:)
    real(8), intent(out) :: wi(ibi:ifi, js0:je0, kb:kf)
    real(8), intent(out) :: wo(ibo:ifo, js0:je0, kb:kf)
    integer, intent(in) :: ibi
    integer, intent(in) :: ifi
    integer, intent(in) :: ibo
    integer, intent(in) :: ifo
    integer, intent(in) :: kb
    integer, intent(in) :: kf

    integer :: l
    integer :: li, li0
    integer :: lo, lo0
    integer :: i, j, k
    real(8), allocatable :: wis(:), wir(:)
    real(8), allocatable :: wos(:), wor(:)
     
    li = (ifi - ibi + 1) * (kf - kb + 1) * (je - js + 1)
    li0= (ifi - ibi + 1) * (kf - kb + 1) * (je0- js0+ 1)
    lo = (ifo - ibo + 1) * (kf - kb + 1) * (je - js + 1)
    lo0= (ifo - ibo + 1) * (kf - kb + 1) * (je0- js0+ 1)

    allocate(wis(0:li - 1), wir(0:li0 - 1))
    allocate(wos(0:lo - 1), wor(0:lo0 - 1))

    if(myrki == 0) then
       l = 0
       do j = js , je
       do k = kb , kf 
       do i = ibi, ifi
          wis(l) = w(i,j,k)
          l = l + 1
       end do
       end do
       end do
       if(jprc /= 1) then
          call mpi_allgather(wis, li, mpi_double_precision, &
                             wir, li, mpi_double_precision, mpi_k_world, ier)
       else
          wir(0:li - 1) = wis(0:li - 1)
       end if
       l = 0
       do j = js0, je0
       do k = kb , kf 
       do i = ibi, ifi
          wi(i,j,k) = wir(l)
          l = l + 1
       end do
       end do
       end do
    end if

    if(myrki == iprc - 1) then
       l = 0
       do j = js , je
       do k = kb , kf 
       do i = ibo, ifo
          wos(l) = w(i,j,k)
          l = l + 1
       end do
       end do
       end do
       if(jprc /= 1) then
          call mpi_allgather(wos, lo, mpi_double_precision, &
                             wor, lo, mpi_double_precision, mpi_k_world, ier)
       else
          wor(0:lo - 1) = wos(0:lo - 1)
       end if
       l = 0
       do j = js0, je0
       do k = kb , kf 
       do i = ibo, ifo
          wo(i,j,k) = wor(l)
          l = l + 1
       end do
       end do
       end do
    end if

    deallocate(wis, wir)
    deallocate(wos, wor)

    return

  end subroutine gatherx2
#endif




end module bval_module
