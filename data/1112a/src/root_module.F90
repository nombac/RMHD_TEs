
module root_module

  implicit none

  ! necessary for restart
  real(8) :: time
  integer :: nhy
  real(8) :: dtqq

  real(8) :: thist
  real(8) :: tdump

  character :: histfile*11
  character :: dumpfile*11

  ! not necessary for restart
  real(8) :: dt
  real(8) :: dtmin
  real(8) :: tlim
  integer :: nlim
  integer :: cpulim

  real(8) :: dtdump
  real(8) :: dthist

  integer :: irest
  character :: targ*5
  character :: restfile*11

contains




  subroutine root_initial
    use constant_module, only : PI
    use param_module, only : omega, HUGE
#ifdef SHEAR
    real(8) :: torbit
#endif
    namelist /rescon/ restfile, targ
    namelist /pcon/ nlim, tlim, dtmin, cpulim
    namelist /iocon/ dthist, dtdump
    
    ! read simulation parameters
    open(1, file = 'z3dinput' , status = 'old')
    read(1, rescon) ! irestart, restfile, targ
    read(1, pcon) ! nlim, tlim, dtmin, cpulim
    read(1, iocon) ! dthist, dtdump
    close(1)

    ! read irest
    open(1, file = 'restart.data', status = 'old')
    read(1, *) irest
    close(1)

#ifdef SHEAR
    torbit = (2d0 * PI / omega)

    tlim = tlim * torbit
    dtmin = dtmin * torbit

    dthist = dthist * torbit
    dtdump = dtdump * torbit
#endif
    
    if(irest == 0) then
       time = 0d0
       nhy = 0
       dtqq = HUGE
       
       tdump = 0d0
       thist = 0d0

       write(histfile,'(a3, i3.3, a5)') 'hst', 0, '.data'
       write(dumpfile,'(a3, i3.3, a5)') 'res', 0, '.data'
    endif
    
    return
  end subroutine root_initial




end module root_module
