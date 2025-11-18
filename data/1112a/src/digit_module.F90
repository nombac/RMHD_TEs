
module digit_module
  !  convert integer to digits

  implicit none

  private
  public :: digit2, digit3, digit4, digit5

contains
  
  
  function digit2(n)
    implicit none
    integer, intent(in) :: n
    character :: digit2*2
    character :: ch1, ch2
    ch1 = char(int(mod(n,100)/10) + 48)
    ch2 = char(int(mod(n, 10)/ 1) + 48)
    digit2 = ch1//ch2
    return
  end function digit2
  
  
  function digit3(n)
    implicit none
    integer, intent(in) :: n
    character :: digit3*3
    character :: ch1, ch2, ch3
    ch1 = char(int(mod(n,1000)/100) + 48)
    ch2 = char(int(mod(n, 100)/ 10) + 48)
    ch3 = char(int(mod(n,  10)/  1) + 48)
    digit3 = ch1//ch2//ch3
    return
  end function digit3
  
  
  function digit4(n)
    implicit none
    integer, intent(in) :: n
    character :: digit4*4
    character :: ch1, ch2, ch3, ch4
    ch1 = char(int(mod(n,10000)/1000) + 48)
    ch2 = char(int(mod(n, 1000)/ 100) + 48)
    ch3 = char(int(mod(n,  100)/  10) + 48)
    ch4 = char(int(mod(n,   10)/   1) + 48)
    digit4 = ch1//ch2//ch3//ch4
    return
  end function digit4
  
  
  function digit5(n)
    implicit none
    integer, intent(in) :: n
    character :: digit5*5
    character :: ch1, ch2, ch3, ch4, ch5
    ch1 = char(int(mod(n,100000)/10000) + 48)
    ch2 = char(int(mod(n, 10000)/ 1000) + 48)
    ch3 = char(int(mod(n,  1000)/  100) + 48)
    ch4 = char(int(mod(n,   100)/   10) + 48)
    ch5 = char(int(mod(n,    10)/    1) + 48)
    digit5 = ch1//ch2//ch3//ch4//ch5
    return
  end function digit5
  
  
end module digit_module
