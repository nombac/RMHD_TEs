
module gather_module
#ifdef MPI
  implicit none

  private
  public :: gather, gather_initial

  integer, allocatable :: ids(:,:)
  integer, allocatable :: irc(:,:)
  integer, allocatable :: kks(:), kke(:)
  integer, allocatable :: jjs(:), jje(:)
  integer, allocatable :: iis(:), iie(:)

contains




  subroutine gather_initial
    use grid_module, only : para_range, is0, ie0, js0, je0, ks0, ke0
    use mpi_module
    integer, parameter :: NP = 5
    integer :: j, k, m
    integer :: nrk, is1, ie1, js1, je1, ks1, ke1, inum, jnum, knum

    allocate(ids(0:nprc,0:NP-1))
    allocate(irc(0:nprc-1,0:NP-1))
    allocate(kks(0:NP-1),kke(0:NP-1))
    allocate(jjs(0:NP-1),jje(0:NP-1))
    allocate(iis(0:NP-1),iie(0:NP-1))

    ! pattern 0
    kks(0) = 2; kke(0) = 2
    jjs(0) = 2; jje(0) = 2
    iis(0) = 2; iie(0) = 2
    ! pattern 1
    kks(1) = 2; kke(1) = 2
    jjs(1) = 2; jje(1) = 2
    iis(1) = 2; iie(1) = 3
    ! pattern 2
    kks(2) = 2; kke(2) = 2
    jjs(2) = 2; jje(2) = 3
    iis(2) = 2; iie(2) = 2
    ! pattern 3
    kks(3) = 2; kke(3) = 3
    jjs(3) = 2; jje(3) = 2
    iis(3) = 2; iie(3) = 2
    ! pattern 4
    kks(4) = 0; kke(4) = 0
    jjs(4) = 0; jje(4) = 0
    iis(4) = 0; iie(4) = 0

    ids(0,0:NP-1) = 0
    do j = 0, jprc-1
    do k = 0, kprc-1
       nrk = k + j * kprc
       call para_range(is0, ie0,    1, 0, is1, ie1)
       call para_range(js0, je0, jprc, j, js1, je1)
       call para_range(ks0, ke0, kprc, k, ks1, ke1)
       do m = 0, NP-1
          inum = (ie1 + iie(m)) - (is1 - iis(m)) + 1
          jnum = (je1 + jje(m)) - (js1 - jjs(m)) + 1
          knum = (ke1 + kke(m)) - (ks1 - kks(m)) + 1
          irc(nrk,m) = inum * jnum * knum
          ids(nrk+1,m) = ids(nrk,m) + irc(nrk,m)
       enddo
    enddo
    enddo

    return
  end subroutine gather_initial

 


  subroutine gather(q,m)
    use grid_module, only : para_range, is, ie, js, je, ks, ke, &
                            is0, ie0, js0, je0, ks0, ke0, in, jn, kn
    use mpi_module
    real(8), intent(inout) :: q(in,jn,kn)
    integer, intent(in) :: m
    real(8), allocatable :: tmp0(:)
    real(8), allocatable :: tmp(:)
    integer is2, ie2, js2, je2, ks2, ke2
    integer i, j, k, l
    integer jrk, krk

    allocate(tmp0(0:irc(myrk,m)-1))
    allocate(tmp(0:ids(nprc,m)-1))

    l = 0
    do k = ks-kks(m), ke+kke(m)
    do j = js-jjs(m), je+jje(m)
    do i = is-iis(m), ie+iie(m)
       tmp0(l) = q(i,j,k)
       l = l + 1
    enddo
    enddo
    enddo
    call mpi_gatherv(tmp0, irc(myrk,m), mpi_double_precision, &
         tmp, irc(0:nprc-1,m), ids(0:nprc,m), mpi_double_precision, &
         0, mpi_comm_world, ier)

    if(myrk == 0) then
       l = 0
       do jrk = 0, jprc-1
       do krk = 0, kprc-1
          call para_range(is0, ie0,    1,   0, is2, ie2)
          call para_range(js0, je0, jprc, jrk, js2, je2)
          call para_range(ks0, ke0, kprc, krk, ks2, ke2)
          do k = ks2-kks(m), ke2+kke(m)
          do j = js2-jjs(m), je2+jje(m)
          do i = is2-iis(m), ie2+iie(m)
             q(i,j,k) = tmp(l)
             l = l + 1
          enddo
          enddo
          enddo
       enddo
       enddo
    endif

    deallocate(tmp0)
    deallocate(tmp)

    return
  end subroutine gather



#endif
end module gather_module
