
module initial_module

  implicit none

  private
  public :: initial

contains




  subroutine initial
    use dump_module, only : dump_read, dump_write
    use field_module, only : field_initial
    use flux_module, only : flux_initial
#ifdef MPI
    use gather_module, only : gather_initial
#endif
    use grid_module, only : grid_initial
    use hist_module, only : hist_write
    use param_module, only : param_initial
    use root_module, only : irest, restfile, targ, &
                            dumpfile, histfile, root_initial
    use radiation_module, only : fld

    call grid_initial
    call param_initial
    call root_initial
#ifdef MPI
    call gather_initial
#endif
    call field_initial
    call flux_initial

    if(irest == 0) then
       call fld
       call dump_write(dumpfile)
       call hist_write(histfile)
    else
       call dump_read(restfile,targ)
    endif
    
    irest = 1
    open(1, file = 'restart.data', status = 'old')
    write(1, *) irest
    close(1)

    return
  end subroutine initial




end module initial_module
