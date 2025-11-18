
module strtoi_module

  implicit none
  private
  public :: strtoi

contains




  integer function strtoi(string, istrt, iend)
    character, intent(in) :: string*8
    integer, intent(in) :: istrt
    integer, intent(in) :: iend
    integer :: asciic(48:57), ishift, ival, i
    data asciic / 0,1,2,3,4,5,6,7,8,9 /
    save asciic

    ishift = 1
    ival = 0
    do i = iend, istrt, -1
       ival = ival + ishift * asciic(ichar(string(i:i)))
       ishift = ishift * 10
    enddo
    strtoi = ival

    return
  end function strtoi


  

end module strtoi_module
