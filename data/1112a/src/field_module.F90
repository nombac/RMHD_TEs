module field_module

  implicit none

  real(8), allocatable :: d(:,:,:) ! density
  real(8), allocatable :: e(:,:,:) ! internal energy
  real(8), allocatable :: v1(:,:,:) ! velocity
  real(8), allocatable :: v2(:,:,:) ! velocity
  real(8), allocatable :: v3(:,:,:) ! velocity
  real(8), allocatable :: b1(:,:,:) ! magnetic field
  real(8), allocatable :: b2(:,:,:) ! magnetic field
  real(8), allocatable :: b3(:,:,:) ! magnetic field
  real(8), allocatable :: er(:,:,:) ! radiation energy

  real(8), allocatable :: qmag(:,:,:) ! numerical energy dissipation
  real(8), allocatable :: qkin(:,:,:) ! numerical energy dissipation
  real(8), allocatable :: scol(:,:,:) ! energy diss. via artificial viscosity

  real(8), allocatable :: f11(:,:,:) ! eddington tensor
  real(8), allocatable :: f22(:,:,:) ! eddington tensor
  real(8), allocatable :: f12(:,:,:) ! eddington tensor
  real(8), allocatable :: f13(:,:,:) ! eddington tensor
  real(8), allocatable :: f23(:,:,:) ! eddington tensor
  real(8), allocatable :: dr(:,:,:)
  real(8), allocatable :: fr(:,:,:)
  real(8), allocatable :: dr1(:,:,:)
  real(8), allocatable :: dr2(:,:,:)
  real(8), allocatable :: dr3(:,:,:)
  real(8), allocatable :: rrr(:,:,:)

  real(8), allocatable :: emfx(:,:,:)
  real(8), allocatable :: emfy(:,:,:)
  real(8), allocatable :: emfz(:,:,:)
  real(8), allocatable :: et(:,:,:) ! total energy (ekin + eint + emag)
  real(8), allocatable :: s1(:,:,:) ! momentum
  real(8), allocatable :: s2(:,:,:) ! momentum
  real(8), allocatable :: s3(:,:,:) ! momentum

  real(8), allocatable :: erikb(:,:)
  real(8), allocatable :: erokb(:,:)

  real(8), allocatable :: tv1(:,:,:)
  real(8), allocatable :: tv2(:,:,:)
  real(8), allocatable :: tv3(:,:,:)

  real(8), allocatable :: p(:,:,:)
  real(8), allocatable :: qq1(:,:,:)
  real(8), allocatable :: qq2(:,:,:)
  real(8), allocatable :: qq3(:,:,:)

contains



  subroutine field_initial
    use grid_module, only : in, jn, kn, is, ie, js, je, ks, ke
    use root_module, only : irest
    
    allocate(d(in,jn,kn)) ! density
    allocate(e(in,jn,kn)) ! internal energy
    allocate(v1(in,jn,kn)) ! velocity
    allocate(v2(in,jn,kn)) ! velocity
    allocate(v3(in,jn,kn)) ! velocity
    allocate(b1(in,jn,kn)) ! magnetic field
    allocate(b2(in,jn,kn)) ! magnetic field
    allocate(b3(in,jn,kn)) ! magnetic field
    allocate(er(in,jn,kn)) ! radiation energy
    allocate(dr(in,jn,kn))
    
    allocate(qmag(in,jn,kn)) ! numerical energy dissipation
    allocate(qkin(in,jn,kn)) ! numerical energy dissipation
    allocate(scol(in,jn,kn)) ! energy dissipation via artificial viscosity
    
    allocate(f11(is-2:ie+3,js-2:je+3,ks-2:ke+3)) ! eddington tensor
    allocate(f22(is-2:ie+3,js-2:je+3,ks-2:ke+3)) ! eddington tensor
    allocate(f12(is-2:ie+3,js-2:je+3,ks-2:ke+3)) ! eddington tensor
    allocate(f13(is-2:ie+3,js-2:je+3,ks-2:ke+3)) ! eddington tensor
    allocate(f23(is-2:ie+3,js-2:je+3,ks-2:ke+3)) ! eddington tensor
    allocate(fr(is-2:ie+3,js-2:je+3,ks-2:ke+3))
    allocate(dr1(is-2:ie+3,js-2:je+3,ks-2:ke+3))
    allocate(dr2(is-2:ie+3,js-2:je+3,ks-2:ke+3))
    allocate(dr3(is-2:ie+3,js-2:je+3,ks-2:ke+3))
    allocate(rrr(is-2:ie+3,js-2:je+3,ks-2:ke+3))

    allocate(emfx(in,jn,kn))
    allocate(emfy(in,jn,kn))
    allocate(emfz(in,jn,kn))
    allocate(et(in,jn,kn)) ! total energy (ekin + eint + emag)
    allocate(s1(in,jn,kn)) ! momentum
    allocate(s2(in,jn,kn)) ! momentum
    allocate(s3(in,jn,kn)) ! momentum

    allocate(p(in,jn,kn)) ! total energy (ekin + eint + emag)
    allocate(qq1(in,jn,kn)) ! momentum
    allocate(qq2(in,jn,kn)) ! momentum
    allocate(qq3(in,jn,kn)) ! momentum
    allocate(tv1(in,jn,kn)) ! momentum
    allocate(tv2(in,jn,kn)) ! momentum
    allocate(tv3(in,jn,kn)) ! momentum

    allocate(erikb(in,jn))
    allocate(erokb(in,jn))

    if(irest == 0) then
       call problem
    endif

    return
  end subroutine field_initial



  subroutine problem

    use constant_module, only : PI
    use grid_module, only : in, jn, kn, dx, dz, x1a, x1b, x3a, x3b
    use param_module, only : omega, GAMMA, dfloor, GAMMAR
    use m_rng_uniform_mt, c_rng_uniform => c_rng_uniform_mt

    integer :: i, j, k
    integer :: l, lmax, l0
    real(8), allocatable :: zl(:), dl(:), erl(:), el(:)
    real(8), allocatable :: ay(:,:)
    real(8) :: dc, erc, ec
    real(8) :: a0
    real(8) :: x_loop, z_loop, r_loop, ptc, b0, b0_poloidal, r
    real(8) :: v_sound
    real(8), parameter :: RLOOP = 0.75d0, BETA = 10.d0, BPRATIO = 0.25d0
    real(8), parameter :: AMPLITUDE = 1d-2
    real(8) :: dfl
    type(c_rng_uniform) :: rng    ! for using m_rng_uniform_mt
    real(8), allocatable :: ra(:) ! for using m_rng_uniform_mt
    integer :: nrand              ! for using m_rng_uniform_mt
    integer :: seed = 4357 ! seed for random number
    integer :: ios
    namelist /dflcon/ dfl
    real(8) :: x000, z000

    ! read seed parameter
    open(1,file = 'seed.data', status = 'old', iostat = ios)
    if(ios == 0) then
       read(1,*) seed
    endif
    close(1)

    ! read disk parameters
    open(1,file = 'init/iparam.data', status = 'old')
    read(1, dflcon) !dfl
    close(1)

    ! read initial profiles in the z direction
    open(8, file = 'init/ireso.data', status = 'old')
    read(8,*) lmax
    close(8)

    allocate(zl(lmax), dl(lmax), erl(lmax), el(lmax))

    open(8, file = 'init/initc.data', status = 'old')
    do l = 1, lmax
       read(8,*) zl(l), dl(l), erl(l), el(l)
    enddo
    close(8)

    ! values at midplane
    dc = dl(1)
    erc = erl(1)
    ec = el(1)

    ! set density floor
    dfloor = dc * dfl

    do k = 1, kn
       l0 = 1
       do l = 1, lmax-1
          if(zl(l) < abs(x3b(k))) then
             l0 = l
          endif
       enddo
       do j = 1, jn
       do i = 1, in
          d(i,j,k) = ((zl(l0+1) - abs(x3b(k))) * dl(l0  ) + &
               (abs(x3b(k)) - zl(l0  )) * dl(l0+1)) / (zl(l0+1) - zl(l0))
          e(i,j,k) = ((zl(l0+1) - abs(x3b(k))) * el(l0  ) + &
               (abs(x3b(k)) - zl(l0  )) * el(l0+1)) / (zl(l0+1) - zl(l0))
          er(i,j,k) = ((zl(l0+1) - abs(x3b(k))) * erl(l0  ) + &
               (abs(x3b(k)) - zl(l0  )) * erl(l0+1)) / (zl(l0+1) - zl(l0))
       enddo
       enddo
    enddo

    deallocate(zl, dl, erl, el)

    ! search disk surface
    do k = kn, 1, -1
       if(d(in/2,jn/2,k) > dfloor) exit
    end do
    z000 = (x3b(k) * 2d0) * 0.8d0
    x000 = (x1a(in-2) - x1a(3))
    
    !*** magnetic field ***
    x_loop = 0d0
    z_loop = 0d0
    r_loop = RLOOP * (x000 * 0.5d0)

    ptc = ec * (GAMMA - 1d0) + erc * (GAMMAR - 1d0) ! total pressure
    b0 = sqrt(8d0 * PI * ptc / BETA)
    b0_poloidal = b0 * BPRATIO
    a0 = r_loop * b0_poloidal / PI
    allocate(ay(in,kn))
    do k = 1, kn
    do i = 1, in
       r = sqrt((x1a(i) - x_loop) ** 2 + ((x3a(k) - z_loop) / (z000 / x000)) ** 2)
       if(r < r_loop) then
          ay(i,k) = a0 * (1d0 + cos((r / r_loop) * PI))
       else
          ay(i,k) = 0d0
       endif
    enddo
    enddo
    b1(1:in,1:jn,1:kn) = 0d0
    b2(1:in,1:jn,1:kn) = 0d0
    b3(1:in,1:jn,1:kn) = 0d0


    ! zero net field
    do k = 2, kn - 1
    do j = 2, jn - 1
    do i = 2, in - 2
       b1(i,j,k) =  (ay(i  ,k+1) - ay(i  ,k-1)) / (2d0 * dz)
       b3(i,j,k) = -(ay(i+1,k  ) - ay(i-1,k  )) / (2d0 * dx)
       if(b1(i,j,k) ** 2 + b3(i,j,k) ** 2 /= 0d0) then
          b2(i,j,k) = sqrt(b0 ** 2 - (b1(i,j,k) ** 2 + b3(i,j,k) ** 2))
       else
          b2(i,j,k) = 0d0
       endif
    enddo
    enddo
    enddo


    ! unit conversion
    b1(1:in,1:jn,1:kn) = b1(1:in,1:jn,1:kn) / sqrt(4d0 * PI)
    b2(1:in,1:jn,1:kn) = b2(1:in,1:jn,1:kn) / sqrt(4d0 * PI)
    b3(1:in,1:jn,1:kn) = b3(1:in,1:jn,1:kn) / sqrt(4d0 * PI)
    deallocate(ay)

    
    !*** velocity filed
    allocate(ra(1:2*(in*jn*kn)))
    call init(rng, SEED)
    call gen_rand_array(rng, ra(1:2*(in*jn*kn))) ! using m_rng_uniform_mt
    nrand = 0
    do k = 1, kn
    do j = 1, jn
    do i = 1, in
       v_sound = sqrt( 4d0 / 9d0             * er(i,j,k) / d(i,j,k) + &
                       GAMMA * (GAMMA - 1d0) *  e(i,j,k) / d(i,j,k)) 
       nrand = nrand + 1
       v1(i,j,k) = -AMPLITUDE * v_sound * 2d0 * (ra(nrand) - 0.5d0)
       nrand = nrand + 1
       v3(i,j,k) = -AMPLITUDE * v_sound * 2d0 * (ra(nrand) - 0.5d0)
       v2(i,j,k) = -1.5d0 * omega * x1b(i)
    enddo
    enddo
    enddo
    call final(rng) ! using m_rng_uniform_mt

    return
  end subroutine problem



end module field_module
