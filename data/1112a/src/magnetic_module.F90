
module magnetic_module

  implicit none

  private
  public :: magnetic

contains



  
  subroutine magnetic
    use convert_module, only : convert_e2et, convert_et2e
    use field_module, only : qmag
    use flux_module, only : sdem
    use grid_module, only : is, ie, js, je, ks, ke
    use root_module, only : dt

    integer :: i, j, k
    real(8) :: tma(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: tmb(is-2:ie+3,js-2:je+3,ks-2:ke+3)

#ifdef PARTETOT
    call convert_e2et(1)
#endif
    call mocct
    call moctf
    call magpres
#ifdef PARTETOT
    call convert_et2e(1,tma,tmb)
    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       qmag(i,j,k) =                 tma(i,j,k) / dt
       sdem(i,j,k) =                 tmb(i,j,k) / dt
    enddo
    enddo
    enddo
#endif

    return
  end subroutine magnetic
  



  subroutine mocct
    use bval_module, only : bvalx1, bvalx2, bvalx3
    use field_module, only : d, v1, v2, v3, b1, b2, b3, et, emfx, emfy, emfz
    use flux_module, only : fpf1, fpf3
    use grid_module, only : is, ie, js, je, ks, ke, dx, dy, dz, &
                            dx1a, dx2a, dx3a
    use root_module, only : dt
    use param_module, only : eta00, fdamp
    integer :: i, j, k
    real(8) :: delh(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: delv(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: dbzdx(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: dvzdx(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: dbxdz(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: dvxdz(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: dbydz(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: dvydz(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: dbzdy(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: dvzdy(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: dbxdy(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: dvxdy(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: dbydx(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: dvydx(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: rhop, rhon, bave, vave
    real(8) :: vcp, vcn, vpv, vnv, hpv, hnv
    real(8) :: qa, qb, qc, qd, sigma, sgn
    real(8) :: bxs, bys, bzs, vxs, vys, vzs
    real(8) :: bxsvz, vxsbz, bzsvx, vzsbx
    real(8) :: bzsvy, vzsby, bysvz, vysbz
    real(8) :: bxsvy, vxsby, bysvx, vysbx
#ifdef PARTETOT
    real(8) :: bxstar(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: bystar(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: bzstar(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: pf1(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: pf2(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: pf3(is-2:ie+3,js-2:je+3,ks-2:ke+3)
#endif
    real(8) :: eta0
    !=======================================================================

    eta0 = eta00 * (min(dx,dy,dz))**2 / dt

    ! *** begin to compute ey ***
    !     dvz/dx, dbz/dx
    do k = ks  , ke+1
    do i = is-2, ie+1
!CDIR SHORT_LOOP
    do j = js  , je
       delh(i,j,k) = (b3(i+1,j,k) - b3(i,j,k)) / dx
       delv(i,j,k) = (v3(i+1,j,k) - v3(i,j,k)) / dx
    enddo
    enddo
    enddo
    do k = ks  , ke+1
    do i = is-1, ie+1
!CDIR SHORT_LOOP
    do j = js  , je
       qa = delv(i,j,k) * delv(i-1,j,k)
       qb = delv(i,j,k) + delv(i-1,j,k)
       dvzdx(i,j,k) = 0
       if(qa > 0) dvzdx(i,j,k) = qa / qb
       qc = delh(i,j,k) * delh(i-1,j,k)
       qd = delh(i,j,k) + delh(i-1,j,k)
       dbzdx(i,j,k) = 0
       if(qc > 0) dbzdx(i,j,k) = qc / qd
    enddo
    enddo
    enddo

    !     dvx/dz, dbx/dz
    do k = ks-2, ke+1
    do i = is  , ie+1
!CDIR SHORT_LOOP
    do j = js  , je
       delh(i,j,k) = (b1(i,j,k+1) - b1(i,j,k)) / dz
       delv(i,j,k) = (v1(i,j,k+1) - v1(i,j,k)) / dz
    enddo
    enddo
    enddo
    do k = ks-1, ke+1
    do i = is  , ie+1
!CDIR SHORT_LOOP
    do j = js  , je
       qa = delv(i,j,k) * delv(i,j,k-1)
       qb = delv(i,j,k) + delv(i,j,k-1)
       dvxdz(i,j,k) = 0
       if(qa > 0) dvxdz(i,j,k) = qa / qb
       qc = delh(i,j,k) * delh(i,j,k-1)
       qd = delh(i,j,k) + delh(i,j,k-1)
       dbxdz(i,j,k) = 0
       if(qc > 0) dbxdz(i,j,k) = qc / qd
    enddo
    enddo
    enddo

    !     [ey]
    do k = ks, ke+1
    do i = is, ie+1
!CDIR SHORT_LOOP
    do j = js, je
       !     bzstar * vx, vzstar * bx
       sigma = .5d0 * (v3(i,j,k) + v3(i-1,j,k)) * dt
       if(sigma >= 0) then
          bave = b1(i,j,k-1) + (dz - sigma) * dbxdz(i,j,k-1)
          vave = v1(i,j,k-1) + (dz - sigma) * dvxdz(i,j,k-1)
       else
          bave = b1(i,j,k  ) - (dz + sigma) * dbxdz(i,j,k  )
          vave = v1(i,j,k  ) - (dz + sigma) * dvxdz(i,j,k  )
       endif
       rhop = sqrt(.5d0 * (d(i  ,j,k) + d(i  ,j,k-1)))
       rhon = sqrt(.5d0 * (d(i-1,j,k) + d(i-1,j,k-1)))
       vcp  = vave - abs(bave) / rhop
       vcn  = vave + abs(bave) / rhon
       sigma = vcp * dt
       if(sigma >= 0) then
          hpv = b3(i-1,j,k) + (dx - sigma) * dbzdx(i-1,j,k)
          vpv = v3(i-1,j,k) + (dx - sigma) * dvzdx(i-1,j,k)
       else
          hpv = b3(i  ,j,k) - (dx + sigma) * dbzdx(i  ,j,k)
          vpv = v3(i  ,j,k) - (dx + sigma) * dvzdx(i  ,j,k)
       endif
       sigma = vcn * dt
       if(sigma >= 0) then
          hnv = b3(i-1,j,k) + (dx - sigma) * dbzdx(i-1,j,k)
          vnv = v3(i-1,j,k) + (dx - sigma) * dvzdx(i-1,j,k)
       else
          hnv = b3(i  ,j,k) - (dx + sigma) * dbzdx(i  ,j,k)
          vnv = v3(i  ,j,k) - (dx + sigma) * dvzdx(i  ,j,k)
       endif
       sgn = sign(1d0,bave)
       bzs = sgn * ( (vpv - vnv) + sgn * (hpv / rhop + hnv / rhon)) &
            / (1d0 / rhop + 1d0 / rhon)
       vzs = .5d0 * ( (vpv + vnv) &
          + bzs * sgn * (1d0 / rhon - 1d0 / rhop) &
          + sgn * (hpv / rhop - hnv / rhon) )
#ifdef PARTETOT
       bzstar(i,j,k) = bzs
#endif
       vzsbx = vzs * bave
       bzsvx = bzs * vave

       !     bxstar * vz, vxstar * bz
       sigma = .5d0 * (v1(i,j,k) + v1(i,j,k-1)) * dt
       if(sigma >= 0.0) then
          bave = b3(i-1,j,k) + (dx - sigma) * dbzdx(i-1,j,k)
          vave = v3(i-1,j,k) + (dx - sigma) * dvzdx(i-1,j,k)
       else
          bave = b3(i  ,j,k) - (dx + sigma) * dbzdx(i,j,k  )
          vave = v3(i  ,j,k) - (dx + sigma) * dvzdx(i,j,k  )
       endif

       rhop = sqrt(.5d0 * (d(i,j,k  ) + d(i-1,j,k  )))
       rhon = sqrt(.5d0 * (d(i,j,k-1) + d(i-1,j,k-1)))
       vcn  = vave + abs(bave) / rhon
       vcp  = vave - abs(bave) / rhop
       sigma = vcp * dt
       if(sigma >= 0) then
          hpv = b1(i,j,k-1) + (dz - sigma) * dbxdz(i,j,k-1)
          vpv = v1(i,j,k-1) + (dz - sigma) * dvxdz(i,j,k-1)
       else
          hpv = b1(i,j,k  ) - (dz + sigma) * dbxdz(i,j,k  )
          vpv = v1(i,j,k  ) - (dz + sigma) * dvxdz(i,j,k  )
       endif
       sigma = vcn * dt
       if(sigma >= 0) then
          hnv = b1(i,j,k-1) + (dz - sigma) * dbxdz(i,j,k-1)
          vnv = v1(i,j,k-1) + (dz - sigma) * dvxdz(i,j,k-1)
       else
          hnv = b1(i,j,k  ) - (dz + sigma) * dbxdz(i,j,k  )
          vnv = v1(i,j,k  ) - (dz + sigma) * dvxdz(i,j,k  )
       endif
       sgn = sign(1d0,bave)
       bxs = sgn * ( (vpv - vnv) + sgn * (hpv / rhop + hnv / rhon) )&
            / (1d0 / rhop + 1d0 / rhon)
       vxs = .5d0 * ( (vpv + vnv) &
            + bxs * sgn * (1d0 / rhon - 1d0 / rhop)&
            + sgn * (hpv / rhop - hnv / rhon ) )
#ifdef PARTETOT
       bxstar(i,j,k) = bxs
#endif
       bxsvz = bxs * vave
       vxsbz = vxs * bave

       !     ey = vz bx - vx bz
       emfy(i,j,k) = .5d0 * ((vzsbx + bxsvz) - (vxsbz + bzsvx)) &
#ifdef DAMPING
            - eta0 * ((b1(i,j,k) - b1(i  ,j,k-1)) / dz &
                     -(b3(i,j,k) - b3(i-1,j,k  )) / dx) &
            * ((fdamp(k) + fdamp(k-1)) / 2d0)
#endif
    enddo
    enddo
    enddo

    call bvalx3(emfy,is  ,ie+1,js  ,je  , 1)
    call bvalx1(emfy,          ks-2,ke+3, 0)

#ifdef PARTETOT
    do k = ks, ke
    do i = is, ie+1
!CDIR SHORT_LOOP
    do j = js, je
       pf1(i,j,k) = -.5d0 * (emfy(i,j,k  ) * bzstar(i,j,k  ) &
                           + emfy(i,j,k+1) * bzstar(i,j,k+1))
    enddo
    enddo
    enddo
    do k = ks, ke+1
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       pf3(i,j,k) =  .5d0 * (emfy(i  ,j,k) * bxstar(i  ,j,k) &
                           + emfy(i+1,j,k) * bxstar(i+1,j,k))
    enddo
    enddo
    enddo
#endif


    ! *** begin to compute ez ***
    !     dvy/dx, dby/dx
    do k = ks  , ke
    do i = is-2, ie+1
!CDIR SHORT_LOOP
    do j = js  , je+1
       delh(i,j,k) = (b2(i+1,j,k) - b2(i,j,k)) / dx
       delv(i,j,k) = (v2(i+1,j,k) - v2(i,j,k)) / dx
    enddo
    enddo
    enddo
    do k = ks  , ke
    do i = is-1, ie+1
!CDIR SHORT_LOOP
    do j = js  , je+1
       qa = delv(i,j,k) * delv(i-1,j,k)
       qb = delv(i,j,k) + delv(i-1,j,k)
       dvydx(i,j,k) = 0
       if(qa > 0) dvydx(i,j,k) = qa / qb
       qc = delh(i,j,k) * delh(i-1,j,k)
       qd = delh(i,j,k) + delh(i-1,j,k)
       dbydx(i,j,k) = 0
       if(qc > 0) dbydx(i,j,k) = qc / qd
    enddo
    enddo
    enddo

    !     dvx/dy, dbx/dy
    do k = ks  , ke
    do i = is  , ie+1
!CDIR SHORT_LOOP
    do j = js-2, je+1
       delh(i,j,k) = (b1(i,j+1,k) - b1(i,j,k)) / dy
       delv(i,j,k) = (v1(i,j+1,k) - v1(i,j,k)) / dy
    enddo
    enddo
    enddo
    do k = ks  , ke
    do i = is  , ie+1
!CDIR SHORT_LOOP
    do j = js-1, je+1
       qa = delv(i,j,k) * delv(i,j-1,k)
       qb = delv(i,j,k) + delv(i,j-1,k)
       dvxdy(i,j,k) = 0
       if(qa > 0) dvxdy(i,j,k) = qa / qb
       qc = delh(i,j,k) * delh(i,j-1,k)
       qd = delh(i,j,k) + delh(i,j-1,k)
       dbxdy(i,j,k) = 0
       if(qc > 0) dbxdy(i,j,k) = qc / qd
    enddo
    enddo
    enddo

    !     [ez]
    do k = ks, ke
    do i = is, ie+1
!CDIR SHORT_LOOP
    do j = js, je+1
       !     bystar * vx, vystar * bx
       sigma = .5d0 * (v2(i,j,k) + v2(i-1,j,k)) * dt
       if(sigma >= 0) then
          bave = b1(i,j-1,k) + (dy - sigma) * dbxdy(i,j-1,k)
          vave = v1(i,j-1,k) + (dy - sigma) * dvxdy(i,j-1,k)
       else
          bave = b1(i,j  ,k) - (dy + sigma) * dbxdy(i,j  ,k)
          vave = v1(i,j  ,k) - (dy + sigma) * dvxdy(i,j  ,k)
       endif
       rhop = sqrt(.5d0 * (d(i  ,j,k) + d(i  ,j-1,k)))
       rhon = sqrt(.5d0 * (d(i-1,j,k) + d(i-1,j-1,k)))
       vcn  = vave + abs(bave) / rhon
       vcp  = vave - abs(bave) / rhop
       sigma = vcp * dt
       if(sigma >= 0) then
          hpv = b2(i-1,j,k) + (dx - sigma) * dbydx(i-1,j,k)
          vpv = v2(i-1,j,k) + (dx - sigma) * dvydx(i-1,j,k)
       else
          hpv = b2(i  ,j,k) - (dx + sigma) * dbydx(i  ,j,k)
          vpv = v2(i  ,j,k) - (dx + sigma) * dvydx(i  ,j,k)
       endif
       sigma = vcn * dt
       if(sigma >= 0) then
          hnv = b2(i-1,j,k) + (dx - sigma) * dbydx(i-1,j,k)
          vnv = v2(i-1,j,k) + (dx - sigma) * dvydx(i-1,j,k)
       else
          hnv = b2(i  ,j,k) - (dx + sigma) * dbydx(i  ,j,k)
          vnv = v2(i  ,j,k) - (dx + sigma) * dvydx(i  ,j,k)
       endif
       sgn = sign(1d0,bave)
       bys = sgn * ((vpv - vnv) + sgn * (hpv / rhop + hnv / rhon)) &
            / (1d0 / rhop + 1d0 / rhon)
       vys = .5d0 * ((vpv + vnv) &
            + bys * sgn * (1d0 / rhon - 1d0 / rhop)&
            + sgn * (hpv / rhop - hnv / rhon))
#ifdef PARTETOT
       bystar(i,j,k) = bys
#endif
       bysvx = bys * vave
       vysbx = vys * bave

       !     vxstar * by, bxstar * vy
       sigma = .5d0 * (v1(i,j-1,k) + v1(i,j,k)) * dt
       if(sigma >= 0) then
          bave = b2(i-1,j,k) + (dx - sigma) * dbydx(i-1,j,k)
          vave = v2(i-1,j,k) + (dx - sigma) * dvydx(i-1,j,k)
       else
          bave = b2(i  ,j,k) - (dx + sigma) * dbydx(i  ,j,k)
          vave = v2(i  ,j,k) - (dx + sigma) * dvydx(i  ,j,k)
       endif
       rhop = sqrt(.5d0 * (d(i,j  ,k) + d(i-1,j  ,k)))
       rhon = sqrt(.5d0 * (d(i,j-1,k) + d(i-1,j-1,k)))
       vcn  = vave + abs(bave) / rhon
       vcp  = vave - abs(bave) / rhop
       sigma = vcp * dt
       if(sigma >= 0) then
          hpv = b1(i,j-1,k) + (dy - sigma) * dbxdy(i,j-1,k)
          vpv = v1(i,j-1,k) + (dy - sigma) * dvxdy(i,j-1,k)
       else
          hpv = b1(i,j  ,k) - (dy + sigma) * dbxdy(i,j  ,k)
          vpv = v1(i,j  ,k) - (dy + sigma) * dvxdy(i,j  ,k)
       endif
       sigma = vcn * dt
       if(sigma >= 0) then
          hnv = b1(i,j-1,k) + (dy - sigma) * dbxdy(i,j-1,k)
          vnv = v1(i,j-1,k) + (dy - sigma) * dvxdy(i,j-1,k)
       else
          hnv = b1(i,j  ,k) - (dy + sigma) * dbxdy(i,j  ,k)
          vnv = v1(i,j  ,k) - (dy + sigma) * dvxdy(i,j  ,k)
       endif
       sgn = sign(1d0,bave)
       bxs = sgn * ((vpv - vnv) + sgn * (hpv / rhop + hnv / rhon))&
            / (1d0 / rhop + 1d0 / rhon)
       vxs = .5d0 * ((vpv + vnv) &
            + bxs * sgn * (1d0 / rhon - 1d0 / rhop)&
            + sgn * (hpv / rhop - hnv / rhon))
#ifdef PARTETOT
       bxstar(i,j,k) = bxs
#endif
       vxsby = vxs * bave
       bxsvy = bxs * vave
       
       !     ez = vx by - vy bx
       emfz(i,j,k) = .5d0 *((vxsby + bysvx) - (vysbx + bxsvy))&
#ifdef DAMPING
            - eta0 * ((b2(i,j,k) - b2(i-1,j  ,k)) / dx&
                     -(b1(i,j,k) - b1(i  ,j-1,k)) / dy)&
                     * fdamp(k)
#endif
    enddo
    enddo
    enddo

    call bvalx3(emfz,is  ,ie+1,js  ,je+1, 0)

#ifdef PARTETOT
    do k = ks, ke
    do i = is, ie+1
!CDIR SHORT_LOOP
    do j = js, je
       pf1(i,j,k) = pf1(i,j,k) + .5d0 * (emfz(i,j  ,k) * bystar(i,j  ,k)&
                                       + emfz(i,j+1,k) * bystar(i,j+1,k))
    enddo
    enddo
    enddo
    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je+1
       pf2(i,j,k) =            - .5d0 * (emfz(i  ,j,k) * bxstar(i  ,j,k)&
                                       + emfz(i+1,j,k) * bxstar(i+1,j,k))
    enddo 
    enddo
    enddo
#endif


    ! *** begin to compute ex ***
    !     dvy/dz, dby/dz
    do k = ks-2, ke+1
    do i = is  , ie
!CDIR SHORT_LOOP
    do j = js  , je+1
       delh(i,j,k) = (b2(i,j,k+1) - b2(i,j,k)) / dz
       delv(i,j,k) = (v2(i,j,k+1) - v2(i,j,k)) / dz
    enddo
    enddo
    enddo
    do k = ks-1, ke+1
    do i = is  , ie
!CDIR SHORT_LOOP
    do j = js  , je+1
       qa = delv(i,j,k) * delv(i,j,k-1)
       qb = delv(i,j,k) + delv(i,j,k-1)
       dvydz(i,j,k) = 0
       if(qa > 0) dvydz(i,j,k) = qa / qb
       qc = delh(i,j,k) * delh(i,j,k-1)
       qd = delh(i,j,k) + delh(i,j,k-1)
       dbydz(i,j,k) = 0
       if(qc > 0) dbydz(i,j,k) = qc / qd
    enddo
    enddo
    enddo

    !     dvz/dy, dbz/dy
    do k = ks  , ke+1
    do i = is  , ie
!CDIR SHORT_LOOP
    do j = js-2, je+1
       delh(i,j,k) = (b3(i,j+1,k) - b3(i,j,k)) / dy
       delv(i,j,k) = (v3(i,j+1,k) - v3(i,j,k)) / dy
    enddo
    enddo
    enddo
    do k = ks  , ke+1
    do i = is  , ie
!CDIR SHORT_LOOP
    do j = js-1, je+1
       qa = delv(i,j,k) * delv(i,j-1,k)
       qb = delv(i,j,k) + delv(i,j-1,k)
       dvzdy(i,j,k) = 0
       if(qa > 0) dvzdy(i,j,k) = qa / qb
       qc = delh(i,j,k) * delh(i,j-1,k)
       qd = delh(i,j,k) + delh(i,j-1,k)
       dbzdy(i,j,k) = 0
       if(qc > 0) dbzdy(i,j,k) = qc / qd
    enddo
    enddo
    enddo

    !     [ex]
    do k = ks, ke+1
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je+1
       !     bystar * vz, vystar * bz
       sigma = .5d0 * (v2(i,j,k) + v2(i,j,k-1)) * dt
       if(sigma >= 0) then
          bave = b3(i,j-1,k) + (dy - sigma) * dbzdy(i,j-1,k)
          vave = v3(i,j-1,k) + (dy - sigma) * dvzdy(i,j-1,k)
       else
          bave = b3(i,j  ,k) - (dy + sigma) * dbzdy(i,j  ,k)
          vave = v3(i,j  ,k) - (dy + sigma) * dvzdy(i,j  ,k)
       endif
       rhop = sqrt(.5d0 * (d(i,j,k  ) + d(i,j-1,k  )))
       rhon = sqrt(.5d0 * (d(i,j,k-1) + d(i,j-1,k-1)))
       vcn  = vave + abs(bave) / rhon
       vcp  = vave - abs(bave) / rhop
       sigma = vcp * dt
       if(sigma >= 0) then
          hpv = b2(i,j,k-1) + (dz - sigma) * dbydz(i,j,k-1)
          vpv = v2(i,j,k-1) + (dz - sigma) * dvydz(i,j,k-1)
       else
          hpv = b2(i,j,k  ) - (dz + sigma) * dbydz(i,j,k  )
          vpv = v2(i,j,k  ) - (dz + sigma) * dvydz(i,j,k  )
       endif
       sigma = vcn * dt
       if(sigma >= 0) then
          hnv = b2(i,j,k-1) + (dz - sigma) * dbydz(i,j,k-1)
          vnv = v2(i,j,k-1) + (dz - sigma) * dvydz(i,j,k-1)
       else
          hnv = b2(i,j,k  ) - (dz + sigma) * dbydz(i,j,k  )
          vnv = v2(i,j,k  ) - (dz + sigma) * dvydz(i,j,k  )
       endif
       sgn = sign(1d0,bave)
       bys = sgn * ((vpv - vnv) + sgn * (hpv / rhop + hnv / rhon))&
            / (1d0 / rhop + 1d0 / rhon)
       vys = .5d0 * ((vpv + vnv) &
            + bys * sgn *(1d0 / rhon - 1d0 / rhop)&
            + sgn * (hpv / rhop - hnv / rhon))
#ifdef PARTETOT
       bystar(i,j,k) = bys
#endif
       bysvz = bys * vave
       vysbz = vys * bave

       !     bzstar * vy, vzstar * by
       sigma = .5d0 * (v3(i,j,k) + v3(i,j-1,k)) * dt
       if(sigma >= 0) then
          bave = b2(i,j,k-1) + (dz - sigma) * dbydz(i,j,k-1)
          vave = v2(i,j,k-1) + (dz - sigma) * dvydz(i,j,k-1)
       else
          bave = b2(i,j,k  ) - (dz + sigma) * dbydz(i,j,k  )
          vave = v2(i,j,k  ) - (dz + sigma) * dvydz(i,j,k  )
       endif
       rhop = sqrt(.5d0 * (d(i,j  ,k) + d(i,j  ,k-1)))
       rhon = sqrt(.5d0 * (d(i,j-1,k) + d(i,j-1,k-1)))
       vcn  = vave + abs(bave) / rhon
       vcp  = vave - abs(bave) / rhop
       sigma = vcp * dt
       if(sigma >= 0) then
          hpv = b3(i,j-1,k) + (dy - sigma) * dbzdy(i,j-1,k)
          vpv = v3(i,j-1,k) + (dy - sigma) * dvzdy(i,j-1,k)
       else
          hpv = b3(i,j  ,k) - (dy + sigma) * dbzdy(i,j  ,k)
          vpv = v3(i,j  ,k) - (dy + sigma) * dvzdy(i,j  ,k)
       endif
       sigma = vcn * dt
       if(sigma >= 0) then
          hnv = b3(i,j-1,k) + (dy - sigma) * dbzdy(i,j-1,k)
          vnv = v3(i,j-1,k) + (dy - sigma) * dvzdy(i,j-1,k)
       else
          hnv = b3(i,j  ,k) - (dy + sigma) * dbzdy(i,j  ,k)
          vnv = v3(i,j  ,k) - (dy + sigma) * dvzdy(i,j  ,k)
       endif
       sgn = sign(1d0,bave)
       bzs = sgn * ((vpv - vnv) + sgn * (hpv / rhop + hnv / rhon))&
            / (1d0 / rhop + 1d0 / rhon)
       vzs = .5d0 * ((vpv + vnv) &
            + bzs * sgn *(1d0 / rhon - 1d0 / rhop)&
            + sgn * (hpv / rhop - hnv / rhon))
#ifdef PARTETOT
       bzstar(i,j,k) = bzs
#endif
       bzsvy = bzs * vave
       vzsby = vzs * bave

       !     ex = vy bz - vz by
       emfx(i,j,k) = .5d0 * ((vysbz + bzsvy) - (vzsby + bysvz))&
#ifdef DAMPING
            - eta0 * ((b3(i,j,k) - b3(i,j-1,k  )) / dy&
                     -(b2(i,j,k) - b2(i,j  ,k-1)) / dz)&
                     * ((fdamp(k) + fdamp(k-1)) / 2d0)
#endif
    enddo
    enddo
    enddo

    call bvalx3(emfx,is  ,ie  ,js  ,je+1, 1)

#ifdef PARTETOT
    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je+1
       pf2(i,j,k) = pf2(i,j,k) + .5d0 * (emfx(i,j,k  ) * bzstar(i,j,k  )&
                                       + emfx(i,j,k+1) * bzstar(i,j,k+1))
    enddo
    enddo
    enddo
    do k = ks, ke+1
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       pf3(i,j,k) = pf3(i,j,k) - .5d0 * (emfx(i,j  ,k) * bystar(i,j  ,k)&
                                       + emfx(i,j+1,k) * bystar(i,j+1,k))
    enddo
    enddo
    enddo
#endif


    ! CONSTRAINED TRANSPORT
    do k=ks-2,ke+2
    do i=is  ,ie+1
!CDIR SHORT_LOOP
    do j=js  ,je
       b1(i,j,k) = b1(i,j,k) - dt / dz * (emfy(i,j,k+1) - emfy(i,j,k))&
                             + dt / dy * (emfz(i,j+1,k) - emfz(i,j,k))
    enddo
    enddo
    enddo

    do k = ks-2,ke+2
    do i = is  ,ie
!CDIR SHORT_LOOP
    do j = js  ,je
       b2(i,j,k) = b2(i,j,k) + dt / dz * (emfx(i,j,k+1) - emfx(i,j,k))&
                             - dt / dx * (emfz(i+1,j,k) - emfz(i,j,k))
    enddo
    enddo
    enddo

    do k = ks-2,ke+3
    do i = is  ,ie
!CDIR SHORT_LOOP
    do j = js  ,je
       b3(i,j,k) = b3(i,j,k) + dt / dx * (emfy(i+1,j,k) - emfy(i,j,k))&
                             - dt / dy * (emfx(i,j+1,k) - emfx(i,j,k))
    enddo
    enddo
    enddo
    
    call bvalx1(b1,          ks-2,ke+2, 2)
    call bvalx2(b1,is-2,ie+2,ks-2,ke+2, 0)
    call bvalx1(b2,          ks-2,ke+2, 0)
    call bvalx2(b2,is-2,ie+2,ks-2,ke+2, 0)
    call bvalx1(b3,          ks-2,ke+3, 0)
    call bvalx2(b3,is-2,ie+2,ks-2,ke+3, 0)
!   call bvalx3(b1,is-2,ie+3,js-2,je+2, 0)
!   call bvalx3(b2,is-2,ie+2,js-2,je+3, 0)
!   call bvalx3(b3,is-2,ie+2,js-2,je+2, 0)

#ifdef PARTETOT
    do k = ks ,ke
    do i = is ,ie+1
!CDIR SHORT_LOOP
    do j = js ,je
       fpf1(i,j,k) = pf1(i,j,k)
    enddo
    enddo
    enddo
    do k = ks, ke+1
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       fpf3(i,j,k) = pf3(i,j,k)
    enddo
    enddo
    enddo
    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       et(i,j,k) = et(i,j,k) - dt * (pf1(i+1,j,k) - pf1(i,j,k)) / dx1a(i)
    enddo
    enddo
    enddo
    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       et(i,j,k) = et(i,j,k) - dt * (pf2(i,j+1,k) - pf2(i,j,k)) / dx2a(j)
    enddo
    enddo
    enddo
    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       et(i,j,k) = et(i,j,k) - dt * (pf3(i,j,k+1) - pf3(i,j,k)) / dx3a(k)
    enddo
    enddo
    enddo
#endif

    return
  end subroutine mocct




  subroutine moctf
    use field_module, only : d, v1, v2, v3, b1, b2, b3
    use grid_module, only : is, ie, js, je, ks, ke, dx, dy, dz
    use root_module, only : dt
    integer :: i, j, k
    real(8) :: v1st(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: v2st(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: v3st(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: q1b(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: q1v(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: dqb(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: q2b(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: q2v(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: dqv(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: bstar(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: vap(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: van(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: rhop(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: rhon(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: vp, vn, hp, hn
    real(8) :: valf, xi, sgn, hzed
    
    !=======================================================================
    !   s2 source term from waves in 3-direction  --------------------------
    !
    do k = ks, ke+1
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       rhop(i,j,k) = sqrt(.5d0 * (d(i,j,k  ) + d(i,j-1,k  )))
       rhon(i,j,k) = sqrt(.5d0 * (d(i,j,k-1) + d(i,j-1,k-1)))
       valf = abs(.5d0 * (b3(i,j-1,k) + b3(i,j,k)))
       vap(i,j,k) =  -1d0 * valf / rhop(i,j,k)
       van(i,j,k) =  valf / rhon(i,j,k)
    enddo
    enddo
    enddo

    do k = ks-1, ke+2
    do i = is  , ie
!CDIR SHORT_LOOP
    do j = js  , je
       q1b(i,j,k) = (b2(i,j,k) - b2(i,j,k-1)) / dz
       q1v(i,j,k) = (v2(i,j,k) - v2(i,j,k-1)) / dz
    enddo
    enddo
    enddo

    do k = ks-1, ke+1
    do i = is  , ie
!CDIR SHORT_LOOP
    do j = js  , je
       q2b(i,j,k) = q1b(i,j,k) * q1b(i,j,k+1)
       q2v(i,j,k) = q1v(i,j,k) * q1v(i,j,k+1)
       dqb(i,j,k) = 0
       dqv(i,j,k) = 0
       if(q2b(i,j,k) > 0)dqb(i,j,k) = q2b(i,j,k) / (q1b(i,j,k) + q1b(i,j,k+1))
       if(q2v(i,j,k) > 0)dqv(i,j,k) = q2v(i,j,k) / (q1v(i,j,k) + q1v(i,j,k+1))
    enddo
    enddo
    enddo

    do k = ks, ke+1
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       xi = vap(i,j,k) * dt
       hp = b2(i,j,k  ) - (dz + xi) * dqb(i,j,k  )
       vp = v2(i,j,k  ) - (dz + xi) * dqv(i,j,k  )
       xi = van(i,j,k) * dt
       hn = b2(i,j,k-1) + (dz - xi) * dqb(i,j,k-1)
       vn = v2(i,j,k-1) + (dz - xi) * dqv(i,j,k-1)

       sgn  = sign(1d0,(b3(i,j-1,k) + b3(i,j,k)))
       bstar(i,j,k) = sgn * ((vp - vn)&
          + sgn * (hp / rhop(i,j,k) + hn / rhon(i,j,k)))&
          / (1d0 / rhop(i,j,k) + 1d0 / rhon(i,j,k))
    enddo
    enddo
    enddo

    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       hzed = .25d0 * (b3(i,j,k) + b3(i,j-1,k)+ b3(i,j,k+1) + b3(i,j-1,k+1))
       v2st(i,j,k) = hzed * (bstar(i,j,k+1) - bstar(i,j,k)) / dz
    enddo
    enddo
    enddo

    !
    !   s3 source term from waves in 2-direction  --------------------------
    !
    do k = ks, ke+1
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je+1
       rhop(i,j,k) = sqrt(.5d0 * (d(i,j  ,k) + d(i,j  ,k-1)))
       rhon(i,j,k) = sqrt(.5d0 * (d(i,j-1,k) + d(i,j-1,k-1)))
       valf = abs(.5d0 * (b2(i,j,k-1) + b2(i,j,k)))
       vap(i,j,k) = -1d0 * valf / rhop(i,j,k)
       van(i,j,k) = valf / rhon(i,j,k)
    enddo
    enddo
    enddo

    do k = ks  , ke+1
    do i = is  , ie
!CDIR SHORT_LOOP
    do j = js-1, je+2
       q1b(i,j,k) = (b3(i,j,k) - b3(i,j-1,k)) / dy
       q1v(i,j,k) = (v3(i,j,k) - v3(i,j-1,k)) / dy
    enddo
    enddo
    enddo

    do k = ks  , ke+1
    do i = is  , ie
!CDIR SHORT_LOOP
    do j = js-1, je+1
       q2b(i,j,k) = q1b(i,j,k) * q1b(i,j+1,k)
       q2v(i,j,k) = q1v(i,j,k) * q1v(i,j+1,k)
       dqb(i,j,k) = 0
       dqv(i,j,k) = 0
       if(q2b(i,j,k) > 0)dqb(i,j,k) = q2b(i,j,k) / (q1b(i,j,k) + q1b(i,j+1,k))
       if(q2v(i,j,k) > 0)dqv(i,j,k) = q2v(i,j,k) / (q1v(i,j,k) + q1v(i,j+1,k))
    enddo
    enddo
    enddo

    do k = ks, ke+1
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je+1
       xi = vap(i,j,k) * dt
       hp = b3(i,j  ,k) - (dy + xi) * dqb(i,j  ,k)
       vp = v3(i,j  ,k) - (dy + xi) * dqv(i,j  ,k)
       xi = van(i,j,k) * dt
       hn = b3(i,j-1,k) + (dy - xi) * dqb(i,j-1,k)
       vn = v3(i,j-1,k) + (dy - xi) * dqv(i,j-1,k)
       sgn  = sign(1d0,(b2(i,j,k-1) + b2(i,j,k)))
       bstar(i,j,k) = sgn*((vp - vn)&
          + sgn * (hp / rhop(i,j,k) + hn / rhon(i,j,k)))&
          / (1d0 / rhop(i,j,k) + 1d0 / rhon(i,j,k))
    enddo
    enddo
    enddo

    do k = ks, ke+1
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       hzed = .25d0 * (b2(i,j,k) + b2(i,j,k-1)+ b2(i,j+1,k) + b2(i,j+1,k-1))
       v3st(i,j,k) = hzed * (bstar(i,j+1,k) - bstar(i,j,k)) / dy
    enddo
    enddo
    enddo


    !
    !   s1 source term from waves in 3-direction  --------------------------
    !
    do k = ks, ke+1
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       rhop(i,j,k) = sqrt(.5d0*(d(i,j,k  ) + d(i-1,j,k  )))
       rhon(i,j,k) = sqrt(.5d0*(d(i,j,k-1) + d(i-1,j,k-1)))
       valf = abs(.5d0 * (b3(i-1,j,k) + b3(i,j,k)))
       vap(i,j,k) = -1d0 * valf / rhop(i,j,k)
       van(i,j,k) =  valf / rhon(i,j,k)
    enddo
    enddo
    enddo

    do k = ks-1, ke+2
    do i = is  , ie
!CDIR SHORT_LOOP
    do j = js  , je
       q1b(i,j,k) = (b1(i,j,k) - b1(i,j,k-1)) / dz
       q1v(i,j,k) = (v1(i,j,k) - v1(i,j,k-1)) / dz
    enddo
    enddo
    enddo

    do k = ks-1, ke+1
    do i = is  , ie
!CDIR SHORT_LOOP
    do j = js  , je
       q2b(i,j,k) = q1b(i,j,k) * q1b(i,j,k+1)
       q2v(i,j,k) = q1v(i,j,k) * q1v(i,j,k+1)
       dqb(i,j,k) = 0
       dqv(i,j,k) = 0
       if(q2b(i,j,k) > 0)dqb(i,j,k) = q2b(i,j,k) / (q1b(i,j,k) + q1b(i,j,k+1))
       if(q2v(i,j,k) > 0)dqv(i,j,k) = q2v(i,j,k) / (q1v(i,j,k) + q1v(i,j,k+1))
    enddo
    enddo
    enddo

    do k = ks, ke+1
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       xi = vap(i,j,k) * dt
       hp = b1(i,j,k  ) - (dz + xi) * dqb(i,j,k  )
       vp = v1(i,j,k  ) - (dz + xi) * dqv(i,j,k  )
       xi = van(i,j,k) * dt
       hn = b1(i,j,k-1) + (dz - xi) * dqb(i,j,k-1)
       vn = v1(i,j,k-1) + (dz - xi) * dqv(i,j,k-1)
       sgn  = sign(1d0,(b3(i-1,j,k) + b3(i,j,k)))
       bstar(i,j,k) = sgn*((vp - vn)&
          + sgn * (hp / rhop(i,j,k) + hn / rhon(i,j,k)))&
          / (1d0 / rhop(i,j,k) + 1d0 / rhon(i,j,k))
    enddo
    enddo
    enddo

    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       hzed = .25d0 * (b3(i,j,k) + b3(i-1,j,k) + b3(i,j,k+1) + b3(i-1,j,k+1))
       v1st(i,j,k) = hzed * (bstar(i,j,k+1) - bstar(i,j,k)) / dz
    enddo
    enddo
    enddo


    !
    !   s3 source term from waves in 1-direction  --------------------------
    !
    do k = ks, ke+1
    do i = is, ie+1
!CDIR SHORT_LOOP
    do j = js, je
       rhop(i,j,k) = sqrt(.5d0 * (d(i  ,j,k) + d(i  ,j,k-1)))
       rhon(i,j,k) = sqrt(.5d0 * (d(i-1,j,k) + d(i-1,j,k-1)))
       valf = abs(.5d0 * (b1(i,j,k-1) + b1(i,j,k)))
       vap(i,j,k) = -1d0 * valf / rhop(i,j,k)
       van(i,j,k) = valf / rhon(i,j,k)
    enddo
    enddo
    enddo

    do k = ks  , ke+1
    do i = is-1, ie+2
!CDIR SHORT_LOOP
    do j = js  , je
       q1b(i,j,k) = (b3(i,j,k) - b3(i-1,j,k)) / dx
       q1v(i,j,k) = (v3(i,j,k) - v3(i-1,j,k)) / dx
    enddo
    enddo
    enddo

    do k = ks  , ke+1
    do i = is-1, ie+1
!CDIR SHORT_LOOP
    do j = js  , je
       q2b(i,j,k) = q1b(i,j,k) * q1b(i+1,j,k)
       q2v(i,j,k) = q1v(i,j,k) * q1v(i+1,j,k)
       dqb(i,j,k) = 0
       dqv(i,j,k) = 0
       if(q2b(i,j,k) > 0)dqb(i,j,k) = q2b(i,j,k) / (q1b(i,j,k) + q1b(i+1,j,k))
       if(q2v(i,j,k) > 0)dqv(i,j,k) = q2v(i,j,k) / (q1v(i,j,k) + q1v(i+1,j,k))
    enddo
    enddo
    enddo

    do k = ks, ke+1
    do i = is, ie+1
!CDIR SHORT_LOOP
    do j = js, je
       xi = vap(i,j,k) * dt
       hp = b3(i  ,j,k) - (dx + xi) * dqb(i,j  ,k)
       vp = v3(i  ,j,k) - (dx + xi) * dqv(i,j  ,k)
       xi = van(i,j,k) * dt
       hn = b3(i-1,j,k) + (dx - xi) * dqb(i-1,j,k)
       vn = v3(i-1,j,k) + (dx - xi) * dqv(i-1,j,k)
       sgn  = sign(1d0,(b1(i,j,k-1) + b1(i,j,k)))
       bstar(i,j,k) = sgn*((vp - vn)&
          + sgn * (hp / rhop(i,j,k) + hn / rhon(i,j,k)))&
          / (1d0 / rhop(i,j,k) + 1d0 / rhon(i,j,k))
    enddo
    enddo
    enddo

    do k = ks, ke+1
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       hzed = .25d0 * (b1(i,j,k) + b1(i,j,k-1) + b1(i+1,j,k) + b1(i+1,j,k-1))
       v3st(i,j,k) = v3st(i,j,k) + hzed * (bstar(i+1,j,k) - bstar(i,j,k)) / dx
    enddo
    enddo
    enddo

    !
    !   s1 source term from waves in 2-direction  --------------------------
    !
    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je+1
       rhop(i,j,k) = sqrt(.5d0 * (d(i,j  ,k) + d(i-1,j  ,k)))
       rhon(i,j,k) = sqrt(.5d0 * (d(i,j-1,k) + d(i-1,j-1,k)))
       valf = abs(.5d0 * (b2(i-1,j,k) + b2(i,j,k)))
       vap(i,j,k) = -1d0 * valf / rhop(i,j,k)
       van(i,j,k) = valf / rhon(i,j,k)
    enddo
    enddo
    enddo

    do k = ks  , ke
    do i = is  , ie
!CDIR SHORT_LOOP
    do j = js-1, je+2
       q1b(i,j,k) = (b1(i,j,k) - b1(i,j-1,k)) / dy
       q1v(i,j,k) = (v1(i,j,k) - v1(i,j-1,k)) / dy
    enddo
    enddo
    enddo

    do k = ks  , ke
    do i = is  , ie
!CDIR SHORT_LOOP
    do j = js-1, je+1
       q2b(i,j,k) = q1b(i,j,k) * q1b(i,j+1,k)
       q2v(i,j,k) = q1v(i,j,k) * q1v(i,j+1,k)
       dqb(i,j,k) = 0
       dqv(i,j,k) = 0
       if(q2b(i,j,k) > 0)dqb(i,j,k) = q2b(i,j,k) / (q1b(i,j,k) + q1b(i,j+1,k))
       if(q2v(i,j,k) > 0)dqv(i,j,k) = q2v(i,j,k) / (q1v(i,j,k) + q1v(i,j+1,k))
    enddo
    enddo
    enddo

    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je+1
       xi = vap(i,j,k) * dt
       hp = b1(i,j  ,k) - (dy + xi) * dqb(i,j  ,k)
       vp = v1(i,j  ,k) - (dy + xi) * dqv(i,j  ,k)
       xi = van(i,j,k) * dt
       hn = b1(i,j-1,k) + (dy - xi) * dqb(i,j-1,k)
       vn = v1(i,j-1,k) + (dy - xi) * dqv(i,j-1,k)
       sgn  = sign(1d0,(b2(i-1,j,k) + b2(i,j,k)))
       bstar(i,j,k) = sgn*((vp - vn)&
          + sgn * (hp / rhop(i,j,k) + hn / rhon(i,j,k)))&
          / (1d0 / rhop(i,j,k) + 1d0 / rhon(i,j,k))
    enddo
    enddo
    enddo

    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       hzed = .25d0 * (b2(i,j,k) + b2(i-1,j,k) + b2(i,j+1,k) + b2(i-1,j+1,k))
       v1st(i,j,k) = v1st(i,j,k) + hzed * (bstar(i,j+1,k) - bstar(i,j,k)) / dy
    enddo
    enddo
    enddo

    !
    !   s2 source term from waves in 1-direction  --------------------------
    !
    do k = ks, ke
    do i = is, ie+1
!CDIR SHORT_LOOP
    do j = js, je
       rhop(i,j,k) = sqrt(.5d0 * (d(i  ,j,k) + d(i  ,j-1,k)))
       rhon(i,j,k) = sqrt(.5d0 * (d(i-1,j,k) + d(i-1,j-1,k)))
       valf = abs(.5d0 * (b1(i,j-1,k) + b1(i,j,k)))
       vap(i,j,k) = -1d0 * valf / rhop(i,j,k)
       van(i,j,k) = valf / rhon(i,j,k)
    enddo
    enddo
    enddo

    do k = ks  , ke
    do i = is-1, ie+2
!CDIR SHORT_LOOP
    do j = js  , je
       q1b(i,j,k) = (b2(i,j,k) - b2(i-1,j,k)) / dx
       q1v(i,j,k) = (v2(i,j,k) - v2(i-1,j,k)) / dx
    enddo
    enddo
    enddo

    do k = ks  , ke
    do i = is-1, ie+1
!CDIR SHORT_LOOP
    do j = js  , je
       q2b(i,j,k) = q1b(i,j,k) * q1b(i+1,j,k)
       q2v(i,j,k) = q1v(i,j,k) * q1v(i+1,j,k)
       dqb(i,j,k) = 0
       dqv(i,j,k) = 0
       if(q2b(i,j,k) > 0)dqb(i,j,k) = q2b(i,j,k) / (q1b(i,j,k) + q1b(i+1,j,k))
       if(q2v(i,j,k) > 0)dqv(i,j,k) = q2v(i,j,k) / (q1v(i,j,k) + q1v(i+1,j,k))
    enddo
    enddo
    enddo

    do k = ks, ke
    do i = is, ie+1
!CDIR SHORT_LOOP
    do j = js, je
       xi = vap(i,j,k) * dt
       hp = b2(i  ,j,k) - (dx + xi) * dqb(i,j ,k )
       vp = v2(i  ,j,k) - (dx + xi) * dqv(i,j  ,k)
       xi = van(i,j,k) * dt
       hn = b2(i-1,j,k) + (dx - xi) * dqb(i-1,j,k)
       vn = v2(i-1,j,k) + (dx - xi) * dqv(i-1,j,k)
       sgn  = sign(1d0,(b1(i,j-1,k) + b1(i,j,k)))
       bstar(i,j,k) = sgn*((vp - vn)&
          + sgn * (hp / rhop(i,j,k) + hn / rhon(i,j,k)))&
          / (1d0 / rhop(i,j,k) + 1d0 / rhon(i,j,k))
    enddo
    enddo
    enddo

    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       hzed = .25d0 * (b1(i,j,k) + b1(i,j-1,k) + b1(i+1,j,k) + b1(i+1,j-1,k))
       v2st(i,j,k) = v2st(i,j,k) + hzed * (bstar(i+1,j,k) - bstar(i,j,k)) / dx
    enddo
    enddo
    enddo

!
!  Update velocities with source terms
!
    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       v1(i,j,k) = v1(i,j,k)+ dt * v1st(i,j,k) * 2d0 / (d(i,j,k) + d(i-1,j,k))
    enddo
    enddo
    enddo

    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       v2(i,j,k) = v2(i,j,k)+ dt * v2st(i,j,k) * 2d0 / (d(i,j,k) + d(i,j-1,k))
    enddo
    enddo
    enddo

    do k = ks, ke+1
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       v3(i,j,k) = v3(i,j,k)+ dt * v3st(i,j,k) * 2d0 / (d(i,j,k) + d(i,j,k-1))
    enddo
    enddo
    enddo

    return
  end subroutine moctf




  subroutine magpres
    use bval_module, only : bvalx1, bvalx2, bvalx3
    use field_module, only : d, v1, v2, v3, b1, b2, b3
    use grid_module, only : is, ie, js, je, ks, ke, dx, dy, dz
    use root_module, only : dt
    integer :: i, j, k
    real(8) :: db(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: av(is-2:ie+3,js-2:je+3,ks-2:ke+3)
    real(8) :: ri(is-2:ie+3,js-2:je+3,ks-2:ke+3)

!v1
    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je+1
       db(i,j,k) = (b2(i,j,k) - b2(i-1,j,k)) / dx
       av(i,j,k) = (b2(i,j,k) + b2(i-1,j,k))
       ri(i,j,k) = 2d0/(d(i,j,k) + d(i-1,j,k))
    enddo
    enddo
    enddo

    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       v1(i,j,k) = v1(i,j,k) - dt * 0.5d0 * (db(i,j,k) + db(i,j+1,k)) * &
            0.25d0 * (av(i,j,k) + av(i,j+1,k)) * ri(i,j,k)
    enddo
    enddo
    enddo

    do k = ks, ke+1
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       db(i,j,k) = (b3(i,j,k) - b3(i-1,j,k)) / dx
       av(i,j,k) = (b3(i,j,k) + b3(i-1,j,k))
    enddo
    enddo
    enddo

    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       v1(i,j,k) = v1(i,j,k) - dt * 0.5d0 * (db(i,j,k) + db(i,j,k+1)) * &
            0.25d0 * (av(i,j,k) + av(i,j,k+1)) * ri(i,j,k)
    enddo
    enddo
    enddo

    !v3
    do k = ks, ke+1
    do i = is, ie+1
!CDIR SHORT_LOOP
    do j = js, je
       db(i,j,k) = (b1(i,j,k) - b1(i,j,k-1)) / dz
       av(i,j,k) = (b1(i,j,k) + b1(i,j,k-1))
       ri(i,j,k) = 2d0 / (d(i,j,k) + d(i,j,k-1))
    enddo
    enddo
    enddo

    do k = ks, ke+1
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       v3(i,j,k) = v3(i,j,k) - dt * 0.5d0 * (db(i,j,k) + db(i+1,j,k)) * &
            0.25d0 * (av(i,j,k) + av(i+1,j,k)) * ri(i,j,k)
    enddo
    enddo
    enddo

    do k = ks, ke+1
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je+1
       db(i,j,k) = (b2(i,j,k) - b2(i,j,k-1)) / dz
       av(i,j,k) = (b2(i,j,k) + b2(i,j,k-1))
    enddo
    enddo
    enddo

    do k = ks, ke+1
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       v3(i,j,k) = v3(i,j,k) - dt * 0.5d0 * (db(i,j,k) + db(i,j+1,k)) * &
            0.25d0 * (av(i,j,k) + av(i,j+1,k)) * ri(i,j,k)
    enddo
    enddo
    enddo

    !v2
    do k = ks, ke
    do i = is, ie+1
!CDIR SHORT_LOOP
    do j = js, je
       db(i,j,k) = (b1(i,j,k) - b1(i,j-1,k)) / dy
       av(i,j,k) = (b1(i,j,k) + b1(i,j-1,k))
       ri(i,j,k) = 2d0 / (d(i,j,k) + d(i,j-1,k))
    enddo
    enddo
    enddo

    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       v2(i,j,k) = v2(i,j,k) - dt * 0.5d0 * (db(i,j,k) + db(i+1,j,k)) * &
            0.25d0 * (av(i,j,k) + av(i+1,j,k)) * ri(i,j,k)
    enddo
    enddo
    enddo

    do k = ks, ke+1
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       db(i,j,k) = (b3(i,j,k) - b3(i,j-1,k)) / dy
       av(i,j,k) = (b3(i,j,k) + b3(i,j-1,k))
    enddo
    enddo
    enddo

    do k = ks, ke
    do i = is, ie
!CDIR SHORT_LOOP
    do j = js, je
       v2(i,j,k) = v2(i,j,k) - dt * 0.5d0 * (db(i,j,k) + db(i,j,k+1)) * &
            0.25d0 * (av(i,j,k) + av(i,j,k+1)) * ri(i,j,k)
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
  end subroutine magpres




end module magnetic_module
