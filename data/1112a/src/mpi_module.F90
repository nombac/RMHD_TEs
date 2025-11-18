
module mpi_module
#ifdef MPI
  implicit none
  include 'mpif.h'
  integer :: ju, jd, ku, kd
  integer :: myrki, myrkj, myrkk
  integer :: iprc, jprc, kprc
  integer :: nprc, myrk
  integer :: ier
  integer :: mpi_k_world
  integer :: nprc_k, myrk_k, top_root, bottom_root

contains




  subroutine mpi_begin

    use digit_module

    integer :: irnk
    integer :: j, k
    integer, allocatable :: itab(:,:)
    integer :: icolor, ikey
    
    ! mpi initialization
    call mpi_init(ier)
    call mpi_comm_size(mpi_comm_world, nprc, ier)
    call mpi_comm_rank(mpi_comm_world, myrk, ier)
    
    ! process allocation
    call alloc_proc(nprc, jprc, kprc)
    iprc = 1
    
    ! communication table
    allocate(itab(-1:jprc, -1:kprc))
    do j = -1, jprc
    do k = -1, kprc
       itab(j,k) = mpi_proc_null
    enddo
    enddo
    irnk = 0
    do j = 0, jprc-1
    do k = 0, kprc-1
       itab(j,k) = irnk
       if(myrk == irnk) then
          myrkj = j
          myrkk = k
       endif
       irnk = irnk + 1
    enddo
    enddo
    myrki = 0
    ku = itab(myrkj    , myrkk + 1)
    kd = itab(myrkj    , myrkk - 1)
    ju = itab(myrkj + 1, myrkk    )
    jd = itab(myrkj - 1, myrkk    )
    if(myrkj == 0) then
       jd = itab(jprc-1, myrkk)
    endif
    if(myrkj == jprc-1) then
       ju = itab(     0, myrkk)
    endif
    deallocate(itab)

    top_root = kprc-1
    bottom_root = 0
    icolor = myrkk
    ikey = myrkj
    call mpi_comm_split(mpi_comm_world, icolor, ikey, mpi_k_world, ier)
    call mpi_comm_size(mpi_k_world, nprc_k, ier)
    call mpi_comm_rank(mpi_k_world, myrk_k, ier)

    return

  end subroutine mpi_begin




  subroutine mpi_end

    call mpi_finalize(ier)

    return

  end subroutine mpi_end
  



  subroutine alloc_proc(nprc, jprc, kprc)
    integer, intent(in) :: nprc
    integer, intent(out) :: jprc
    integer, intent(out) :: kprc

    if     (nprc == 1) then
       jprc = 1
       kprc = 1
    else if(nprc == 2) then
       jprc = 1
       kprc = 2
    else if(nprc == 4) then
       jprc = 1
       kprc = 4
    else if(nprc == 8) then
       jprc = 1
       kprc = 8
    else if(nprc == 16) then
       jprc = 1
       kprc = 16
    else if(nprc == 32) then
       jprc = 1
       kprc = 32
    else if(nprc == 56) then
       jprc = 1
       kprc = 56
    else if(nprc == 60) then
       jprc = 1
       kprc = 60
    else if(nprc == 64) then
       jprc = 1
       kprc = 64
    else if(nprc == 80) then
       jprc = 1
       kprc = 80
    else if(nprc == 96) then
       jprc = 1
       kprc = 96
    else if(nprc == 120) then
       jprc = 2
       kprc = 60
    else if(nprc == 128) then
       jprc = 2
       kprc = 64
    else if(nprc == 192) then
       jprc = 2
       kprc = 96
    else if(nprc == 224) then
       jprc = 2
       kprc = 112
    else if(nprc == 256) then
       jprc = 2
       kprc = 128
    else if(nprc == 512) then
       jprc = 4
       kprc = 128
    else
       print *, 'process number error: nprc = ', nprc
       stop
    endif

    return

  end subroutine alloc_proc



#endif
end module mpi_module
