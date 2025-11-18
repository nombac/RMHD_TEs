module flux_module

  implicit none

  real(8), allocatable :: frd1(:,:,:)
  real(8), allocatable :: frd3(:,:,:)
  real(8), allocatable :: fra1(:,:,:)
  real(8), allocatable :: fra3(:,:,:)
  real(8), allocatable :: fgw1(:,:,:)
  real(8), allocatable :: fgw3(:,:,:)
  real(8), allocatable :: fga1(:,:,:)
  real(8), allocatable :: fga3(:,:,:)
  real(8), allocatable :: fpf1(:,:,:)
  real(8), allocatable :: fpf3(:,:,:)
  real(8), allocatable :: srad(:,:,:)
  real(8), allocatable :: spot(:,:,:)
  real(8), allocatable :: sdek(:,:,:)
  real(8), allocatable :: sdem(:,:,:)
  real(8), allocatable :: sdrh(:,:,:)
  real(8), allocatable :: sdv3(:,:,:)
  real(8), allocatable :: sdvp(:,:,:)
  real(8), allocatable :: sena(:,:,:)
  real(8), allocatable :: drhox(:,:,:)
  real(8), allocatable :: divv(:,:,:)

  real(8) :: dadd

contains




  subroutine flux_initial
    ! initialize variavles
    use grid_module, only : is, ie, js, je, ks, ke

    ! initialize variavles
    allocate(frd1(is-2:ie+3,js-2:je+3,ks-2:ke+3))
    allocate(frd3(is-2:ie+3,js-2:je+3,ks-2:ke+3))
    allocate(fra1(is-2:ie+3,js-2:je+3,ks-2:ke+3))
    allocate(fra3(is-2:ie+3,js-2:je+3,ks-2:ke+3))
    allocate(fgw1(is-2:ie+3,js-2:je+3,ks-2:ke+3))
    allocate(fgw3(is-2:ie+3,js-2:je+3,ks-2:ke+3))
    allocate(fga1(is-2:ie+3,js-2:je+3,ks-2:ke+3))
    allocate(fga3(is-2:ie+3,js-2:je+3,ks-2:ke+3))
    allocate(fpf1(is-2:ie+3,js-2:je+3,ks-2:ke+3))
    allocate(fpf3(is-2:ie+3,js-2:je+3,ks-2:ke+3))
    allocate(srad(is-2:ie+3,js-2:je+3,ks-2:ke+3))
    allocate(spot(is-2:ie+3,js-2:je+3,ks-2:ke+3))
    allocate(sdek(is-2:ie+3,js-2:je+3,ks-2:ke+3))
    allocate(sdem(is-2:ie+3,js-2:je+3,ks-2:ke+3))
    allocate(sdrh(is-2:ie+3,js-2:je+3,ks-2:ke+3))
    allocate(sdv3(is-2:ie+3,js-2:je+3,ks-2:ke+3))
    allocate(sdvp(is-2:ie+3,js-2:je+3,ks-2:ke+3))
    allocate(sena(is-2:ie+3,js-2:je+3,ks-2:ke+3))
    allocate(drhox(is-2:ie+3,js-2:je+3,ks-2:ke+3))
    allocate(divv(is-2:ie+3,js-2:je+3,ks-2:ke+3))

    drhox(is-2:ie+3,js-2:je+3,ks-2:ke+3) = 0d0

    dadd = 0d0

    return
  end subroutine flux_initial




end module flux_module
