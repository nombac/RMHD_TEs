
!       vave        have
!
!  01    t/orbit    emag
!  02   dt/orbit    erad
!  03   d           ekin
!  04   e           e
!  05   ekin1       d
!  06   ekin2       chi d dz
!  07   ekin3       lm/dz
!  08   emag1       T
!  09   emag2       spot
!  10   emag3       vz
!  11   erad        qmag
!  12   -b1b2       qkin
!  13   dv1dv2      srad
!  14   P12         sdek+sdem
!  15   -           sdrh
!  16   ekin2d      R
!  17   dadd        D
!  18   -           -b1b2
!  19   -           dv1dv2 
!  20   -           P12
!  21   -           scol
!  22   -           fga3
!  23   -           fpf3
!  24   -           ekin2d
!  25   -           sdv3
!  26   -           gradv:P
!  27   -           frd3
!  28   -           -(4piB-cE)kd
!  29   -           fra3
!  30   -           vA**2
!  31   -           v1**2+dv2**2+v3**2
!  32   -           j**2
!  33   -           emag3 !divv !vr**2
!  34   -    
!  35   niter
!  36   -    
!  37   -    
!  38   -    
!  39   -    
!  40   sdek
!  41   sdem
!  42   sdrh
!  43   sdv3
!  44   srad
!  45   spot
!  46   scol
!  47   qkin
!  48   qmag
!  49   fga1
!  50   fga3
!  51   fgw1
!  52   fgw3
!  53   fpf1
!  54   fpf3
!  55   fra1
!  56   fra3
!  57   frd1
!  58   frd3
!  59   gradv:P
!  60   -(4piB-cE)kd
!  61   dt/dt_D

module hist_module

  implicit none

contains




  subroutine hist_write(filename)
    use constant_module, only : PI, RGAS
    use digit_module, only : digit3
    use field_module, only : d, e, v1, v2, v3, b1, b2, b3, er, &
         qmag, qkin, scol, f12, dr, rrr
    use flux_module, only : drhox, frd1, frd3, fra1, fra3, fgw1, fgw3, &
         fga1, fga3, fpf1, fpf3, srad, spot, sdek, sdem, sdrh, sdv3, &
         sdvp, sena, divv
    use grid_module, only : is, ie, js, je, ks, ke, dx, dy, dz, &
         lx0, ly0, lz0, kn, x1b, dx1a, dx2a, dx3a
#ifdef MPI
    use mpi_module
#endif
    use param_module, only : MMW, GAMMA, omega, GAMMAR
    use radiation_module, only : niter
    use root_module, only : dt, time, histfile, nhy
    use strtoi_module, only : strtoi

    character, intent(in) :: filename*11
    integer, parameter :: NVAVE = 61, NHAVE = 33
    real(8) :: vave(NVAVE)
    real(8), allocatable :: have(:,:)
    real(8) :: drmax
#ifdef MPI
    real(8) :: vave0(NVAVE)
    real(8), allocatable :: have0(:,:)
    real(8) :: drmax0
#endif
    integer :: i, j, k
    integer :: incr
    real(8) :: v1c , v2c , dv2c , v3c , b1c , b2c , b3c
    real(8) :: v1sc, v2sc, dv2sc, v3sc, b1sc, b2sc, b3sc
    real(8) :: cj1, cj2, cj3
    

    allocate(have(kn,NHAVE))
#ifdef MPI
    allocate(have0(kn,NHAVE))
#endif


    ! ---------------------------------------------------------------
    ! volume and horizontally averaging
    vave(1:NVAVE) = 0d0
    have(1:kn,1:NHAVE) = 0d0
    do k = ks, ke
    do j = js, je
    do i = is, ie
       ! density
       have(k, 5) = have(k, 5) + d(i,j,k)
       vave(   3) = vave(   3) + d(i,j,k)
       vave(  17) = vave(  17) + drhox(i,j,k)
       
       ! energy
       ! internal
       have(k, 4) = have(k, 4) + e(i,j,k)
       vave(   4) = vave(   4) + e(i,j,k)
       ! kinetic
       v1c  = 0.5d0 * (v1(i,j,k)    + v1(i+1,j,k)   )
       v2c  = 0.5d0 * (v2(i,j,k)    + v2(i,j+1,k)   )
       v3c  = 0.5d0 * (v3(i,j,k)    + v3(i,j,k+1)   )
       v1sc = 0.5d0 * (v1(i,j,k)**2 + v1(i+1,j,k)**2)
       v2sc = 0.5d0 * (v2(i,j,k)**2 + v2(i,j+1,k)**2)
       v3sc = 0.5d0 * (v3(i,j,k)**2 + v3(i,j,k+1)**2)
       dv2c = v2c + (1.5d0 * omega * x1b(i))
       dv2sc = 0.5*((v2(i,j  ,k) + (1.5d0 * omega * x1b(i)))**2 &
                  + (v2(i,j+1,k) + (1.5d0 * omega * x1b(i)))**2)
       vave(   5) = vave(   5) + (d(i,j,k) * 0.5d0 * v1sc)
       vave(   6) = vave(   6) + (d(i,j,k) * 0.5d0 * v2sc)
       vave(   7) = vave(   7) + (d(i,j,k) * 0.5d0 * v3sc)
       vave(  16) = vave(  16) + (d(i,j,k) * 0.5d0 * dv2sc)
       have(k, 3) = have(k, 3) + (d(i,j,k) * 0.5d0 * (v1sc + v2sc + v3sc))
       have(k,24) = have(k,24) + (d(i,j,k) * 0.5d0 * (v1sc + dv2sc + v3sc))
       ! magnetic
       b1c  = 0.5d0 * (b1(i,j,k)    + b1(i+1,j,k)   )
       b2c  = 0.5d0 * (b2(i,j,k)    + b2(i,j+1,k)   )
       b3c  = 0.5d0 * (b3(i,j,k)    + b3(i,j,k+1)   )
       b1sc = 0.5d0 * (b1(i,j,k)**2 + b1(i+1,j,k)**2)
       b2sc = 0.5d0 * (b2(i,j,k)**2 + b2(i,j+1,k)**2)
       b3sc = 0.5d0 * (b3(i,j,k)**2 + b3(i,j,k+1)**2)
       vave(   8) = vave(   8) + (0.5d0 * b1sc)
       vave(   9) = vave(   9) + (0.5d0 * b2sc)
       vave(  10) = vave(  10) + (0.5d0 * b3sc)
       have(k, 1) = have(k, 1) + (0.5d0 * (b1sc + b2sc + b3sc))
       ! radiation energy
       vave(  11) = vave(  11) + er(i,j,k)
       have(k, 2) = have(k, 2) + er(i,j,k)

       ! stress
       ! magnetic
       vave(  12) = vave(  12) + (-b1c * b2c)
       have(k,18) = have(k,18) + (-b1c * b2c)
       ! Reynolds
       vave(  13) = vave(  13) + (d(i,j,k) * (v1c * dv2c))
       have(k,19) = have(k,19) + (d(i,j,k) * (v1c * dv2c))
       ! radiation
       vave(  14) = vave(  14) + (f12(i,j,k) * er(i,j,k))
       have(k,20) = have(k,20) + (f12(i,j,k) * er(i,j,k))
       
       ! artificial energy injection rate
       ! efloor- related (kinetic)
       vave(  40) = vave(  40) +  sdek(i,j,k)
       ! efloor- related (magnetic)
       vave(  41) = vave(  41) +  sdem(i,j,k)
       have(k,14) = have(k,14) + (sdek(i,j,k) + sdem(i,j,k))
       ! dfloor- related
       vave(  42) = vave(  42) +  sdrh(i,j,k)
       have(k,15) = have(k,15) +  sdrh(i,j,k)
       ! vcap- related
       vave(  43) = vave(  43) +  sdv3(i,j,k)
       have(k,25) = have(k,25) +  sdv3(i,j,k)
       
       ! energy source
       ! radiation force
       vave(  44) = vave(  44) +  srad(i,j,k)
       have(k,13) = have(k,13) +  srad(i,j,k)
#ifdef SHEAR
       ! potential gradient force
       vave(  45) = vave(  45) +  spot(i,j,k)
       have(k, 9) = have(k, 9) +  spot(i,j,k)
#endif
       
       ! energy dissipation rate
       ! kinetic(numerical)
       vave(  47) = vave(  47) +  qkin(i,j,k)
       have(k,12) = have(k,12) +  qkin(i,j,k)
       ! magnetic(numerical)
       vave(  48) = vave(  48) +  qmag(i,j,k)
       have(k,11) = have(k,11) +  qmag(i,j,k)
       ! via artificial viscosity
       vave(  46) = vave(  46) +  scol(i,j,k)
       have(k,21) = have(k,21) +  scol(i,j,k)
       
       ! energy flux
       ! NOTE - SIGN DEFINITION: (come in thru x1 > 0, go out thru x3 >0)
       ! gas advection (x1)
       vave(  49) = vave(  49) - (         fga1(i+1,j,k) - fga1(i,j,k) ) / dx
       ! gas advection (x3)
       vave(  50) = vave(  50) + (         fga3(i,j,k+1) - fga3(i,j,k) ) / dz
       have(k,22) = have(k,22) + (0.5d0 * (fga3(i,j,k+1) + fga3(i,j,k)))   
       ! gas surface-work (x1)
       vave(  51) = vave(  51) - (         fgw1(i+1,j,k) - fgw1(i,j,k) ) / dx
       ! gas surface-work (x3)
       vave(  52) = vave(  52) + (         fgw3(i,j,k+1) - fgw3(i,j,k) ) / dz
       ! Poynting flux (x1)
       vave(  53) = vave(  53) - (         fpf1(i+1,j,k) - fpf1(i,j,k) ) / dx
       ! Poynting flux (x3)
       vave(  54) = vave(  54) + (         fpf3(i,j,k+1) - fpf3(i,j,k) ) / dz
       have(k,23) = have(k,23) + (0.5d0 * (fpf3(i,j,k+1) + fpf3(i,j,k)))   
       ! rad advection (x1)
       vave(  55) = vave(  55) - (         fra1(i+1,j,k) - fra1(i,j,k) ) / dx
       ! rad advection (x3)
       vave(  56) = vave(  56) + (         fra3(i,j,k+1) - fra3(i,j,k) ) / dz
       have(k,29) = have(k,29) + (0.5d0 * (fra3(i,j,k+1) + fra3(i,j,k)))   
       ! rad diffusion (x1)
       vave(  57) = vave(  57) - (         frd1(i+1,j,k) - frd1(i,j,k) ) / dx
       ! rad diffusion (x3)
       vave(  58) = vave(  58) + (         frd3(i,j,k+1) - frd3(i,j,k) ) / dz
       have(k,27) = have(k,27) + (0.5d0 * (frd3(i,j,k+1) + frd3(i,j,k)))   
       
       ! energy exchange rate
       ! gradv:P
       vave(  59) = vave(  59) + sdvp(i,j,k)
       have(k,26) = have(k,26) + sdvp(i,j,k)
       ! -(4piB-cE)kd
       vave(  60) = vave(  60) + sena(i,j,k)
       have(k,28) = have(k,28) + sena(i,j,k)
       
       ! miscellaneous
       ! squared current density
       cj1 = (b3(i,j+1,k) - b3(i,j,k)) / dx2a(j) - &
             (b2(i,j,k+1) - b2(i,j,k)) / dx3a(k)
       cj2 = (b1(i,j,k+1) - b1(i,j,k)) / dx3a(k) - &
             (b3(i+1,j,k) - b3(i,j,k)) / dx1a(i)
       cj3 = (b2(i+1,j,k) - b2(i,j,k)) / dx1a(i) - &
             (b1(i,j+1,k) - b1(i,j,k)) / dx2a(j)
       have(k,32) = have(k,32) + (cj1**2 + cj2**2 + cj3**2)
       ! wavelength of maxium growth over dz
       have(k, 7) = have(k, 7) + &
            ((6.49d0 / omega * sqrt(b3c**2 / (d(i,j,k)))) / dz)
       ! squared Alfven speed
       have(k,30) = have(k,30) + ((b1sc + b2sc + b3sc) / d(i,j,k))
       ! temperature
       have(k, 8) = have(k, 8) + &
            (MMW * e(i,j,k) / RGAS / d(i,j,k) * (GAMMA - 1d0))
       ! vz
       have(k,10) = have(k,10) + v3c
       ! velocity dispersion
       have(k,31) = have(k,31) + (v1sc + dv2sc + v3sc)
       ! diffusion coefficient
       have(k,17) = have(k,17) + dr(i,j,k)
       ! opacity parameter R
       have(k,16) = have(k,16) + rrr(i,j,k)
       ! total opacity per zone
       have(k, 6) = have(k, 6) + &
            (((1d52 * sqrt(d(i,j,k)**11) / sqrt(e(i,j,k)**7)) + &
            (0.33*d(i,j,k))) * dz)
!       ! squared radiation sound speed
!       have(k,33) = have(k,33)+ (GAMMAR * er(i,j,k) / d(i,j,k))
!       ! - divv * p
!       have(k,33) = have(k,33)+ divv(i,j,k)
       have(k,33) = have(k,33)+ (0.5d0 * b3sc)
    enddo
    enddo
    enddo 
#ifdef MPI
    call mpi_reduce(vave, vave0, NVAVE, mpi_double_precision, mpi_sum, 0, &
         mpi_comm_world, ier)
    vave(1:NVAVE) = vave0(1:NVAVE)

    call mpi_reduce(have, have0, NHAVE*kn, mpi_double_precision, mpi_sum, 0, &
         mpi_comm_world, ier)
    have(1:kn,1:NHAVE) = have0(1:kn,1:NHAVE)
#endif
    ! normalization
    vave(1:NVAVE) = vave(1:NVAVE) * (dx * dy * dz) / (lx0 * ly0 * lz0)
    have(1:kn,1:NHAVE) = have(1:kn,1:NHAVE) * (dx * dy) / (lx0 * ly0)
    ! ---------------------------------------------------------------

    ! non averaging quantities
    vave( 1) = time * omega / (2d0 * PI)
    vave( 2) = dt   * omega / (2d0 * PI)
    vave(35) = niter
    drmax = maxval(dr)
#ifdef MPI
    call mpi_allreduce(drmax, drmax0, 1, mpi_double_precision, mpi_max, &
         mpi_comm_world, ier)
    drmax = drmax0
#endif
    vave(61) = dt / (min(dx**2,dy**2,dz**2) / drmax)


    ! write out
#ifdef MPI
    if(myrk == 0) then
#endif
    incr = strtoi(filename,4,6)
    open(2, file = 'hist/hr'//digit3(incr)//'.data')
    write(2,*) NVAVE, NHAVE
    close(2)
    open(3,file = 'hist/'//histfile, form = 'unformatted', position = 'append')
    write(3) vave
    write(3) have
    close(3)
#ifdef MPI
    endif
#endif

    print *, 'history dump written at time/orbit =', &
         time / (2d0 * PI / omega), 'cycle =', nhy, 'file =', filename
    
    deallocate(have)
#ifdef MPI
    deallocate(have0)
#endif


    return
  end subroutine hist_write
  



end module hist_module
