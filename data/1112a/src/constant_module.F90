module constant_module

  implicit none

  real(8), parameter :: PI = 3.1415926535898d+00

  real(8), parameter :: CLIGHT = 2.99792d+10
  real(8), parameter :: SIGMAB = 5.67040d-05
  real(8), parameter :: RGAS   = 8.31447d+07
  real(8), parameter :: ARAD   = SIGMAB / CLIGHT * 4d0
  real(8), parameter :: KBOLTZ = 1.38065d-16
  real(8), parameter :: MELE   = 9.10938d-28
  real(8), parameter :: MPRO   = 1.67262d-24
  real(8), parameter :: SIGMAT = 6.65245d-25

end module constant_module
