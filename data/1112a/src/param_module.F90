module param_module

  implicit none

  real(8), parameter :: TINY = 1d-30
  real(8), parameter :: HUGE = 1d+30

  real(8), parameter :: MMW = 0.61d0 ! mean molecular weight
  real(8), parameter :: MMWE = 1.2d0
  real(8), parameter :: GAMMA = 5.0d0 / 3.0d0 ! specific heat ratio
  real(8), parameter :: GAMMAR = 4.0d0 / 3.0d0

  real(8) :: eta00
  real(8) :: dfloor
  real(8) :: fvcap
  real(8) :: efac

  real(8) :: courno
  real(8) :: qcon

  real(8) :: omega

  real(8), allocatable :: fdamp(:)

contains




  subroutine param_initial
    use constant_module, only : PI
    use grid_module, only : x3a, x3b, kn
    real(8) :: zdmin, zdmax, zabs
    integer :: k, kdamp
    namelist /dampcon/ kdamp
    namelist /hycon/ qcon, courno
    namelist /etacon/ eta00
    namelist /flocon/ efac, fvcap
    namelist /tnrmcon/ omega


    ! read disk parameters
    open(2,file = 'init/iparam.data', status = 'old')
    read(2, tnrmcon) ! omega
    close(2)

    ! read simulation parameters
    open(1, file = 'z3dinput', status = 'old')
    read(1, hycon) ! qcon, courno
    read(1, etacon) ! eta00
    read(1, flocon) ! efac, fvcap
    read(1, dampcon) ! kdamp
    close(1)
    
    allocate(fdamp(kn))

! damping function
    zdmin = x3a(kn-2 - kdamp)
    zdmax = x3a(kn-2)
    do k = 1, kn
       zabs = abs(x3b(k))
       if(zabs < zdmin) then
          fdamp(k) = 0d0
       else if(zabs < zdmax) then
          fdamp(k) = (sin((PI / 2d0) * (((zabs - zdmin) - (zdmax - zabs)) / &
                     (zdmax - zdmin))) + 1d0) / 2d0
       else
          fdamp(k) = 1d0
       endif
    enddo

    return
  end subroutine param_initial




end module param_module
