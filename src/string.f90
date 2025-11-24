! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module string

implicit none

public

interface str
  module procedure :: str_int
end interface str

contains

!********************************************************************************
!
! This function converts uppercase to lowercase
!
!********************************************************************************
function lowercase(str) result(str_res)
  character(len=*), intent(in) :: str
  character(len=len(str))      :: str_res
  !local variables
  character(len=26)            ::  up="ABCDEFGHIJKLMNOPQRSTUVWXYZ"
  character(len=26)            :: low="abcdefghijklmnopqrstuvwxyz"
  integer                      :: i, n, pos
  !
  n=len(str)
  do i = 1, n
    pos = index(up, str(i:i))
    if (pos > 0) then
      str_res(i:i) = low(pos:pos)
    else
      str_res(i:i) = str(i:i)
    end if
  end do
end function lowercase


!********************************************************************************
!
! Check if the system is little endian
!
!********************************************************************************
logical function little_endian()
  use constants, only: i4
  integer(i4)      :: i
  character(len=4) :: c
  character        :: first_byte

  i = 1
  c = transfer(i, c) 
  first_byte = c(1:1)

  little_endian = iachar(first_byte) == 1
end function little_endian


!******************************************************************************** 
!
! Convert numerical vlaue into string
!
!******************************************************************************** 
pure function str_int(val) result(str)
  use constants, only: charlen
  integer, intent(in)           :: val
  character(len=:), allocatable :: str
  !local
  character(len=charlen)        :: temp_str
  
  write(temp_str,"(i10)") val
  str = trim(adjustl(temp_str))
end function str_int


!******************************************************************************** 
!       Routine to separate key and value from read line in INPUT file
!               INPUT:
!                       key     -       contains hole line from the INPUT file
!                       del     -       option to set delimiter                    
!               OUTPUT:
!                       key     -       contains only key string before delimiter
!                       value   -       extracted value from the input line
!                       info    -       1 if delimiter found; 0 if not
!******************************************************************************** 
subroutine split_tag(key, value_, info, del)
  character(len=*), intent(inout)        :: key
  character(len=*), intent(out)          :: value_
  integer, optional, intent(out)         :: info
  character(len=1), optional, intent(in) :: del
  !local variables
  integer                                :: pos
  character(len=1)                       :: del_

  if (present(del)) then
    del_ = del
  else 
    del_ = "="
  end if

  pos = scan(key, "=")

  if (pos > 0) then
    value_ = trim(key(pos+1:))
    key = trim(key(1:pos-1))
    key = adjustl(key)
  else
    pos = scan(key, " ")
    if (pos > 0) then
      value_ = trim(key(pos+1:))
      key = trim(key(1:pos-1))
    end if
  end if
  
  if (present(info)) info = min(1, pos)
end subroutine split_tag


!******************************************************************************** 
!
!       in operator for two strings - is left string present in the right string
!
!******************************************************************************** 
logical function isin(left, right)
  character(len=*), intent(in)  :: left
  character(len=*), intent(in)  :: right
  !local variables
  character(len=:), allocatable :: left_, right_, check
  integer                       :: pos
  !
  isin = .false.
  left_ = trim(lowercase(left))
  right_ = trim(lowercase(right))
  !
  do  
    pos = scan(right_, " ")
    if (pos > 0) then
      check = trim(right_(1:pos-1))
      right_ = trim(right_(pos+1:))
      if (trim(left_) == trim(check)) isin = .true.
    else
      if (trim(left_) == trim(right_)) isin = .true.
      exit
    end if
  end do
end function isin


!******************************************************************************** 
!
!       Extract the rest of the string if the left is found in right
!
!******************************************************************************** 
 function rest_str(left, right, sep) result(str)
  character(len=*), intent(in)           :: left
  character(len=*), intent(in)           :: right
  character(len=:), allocatable          :: str
  character(len=1), optional, intent(in) :: sep    
  !local variables
  character(len=:), allocatable          :: left_, right_, pattern
  character(len=1)                       :: sep_    
  integer                                :: pos
  !
  left_ = trim(lowercase(left))
  right_ = trim(lowercase(right))
  !
  if (present(sep)) then
    sep_ = sep
  else
    sep_ = " "
  end if
  !
  do  
    pos = scan(right_, sep_)
    if (pos > 0) then
      right_ = adjustl(right_(pos+1:))
      pattern = right_(1:len(left_))
      if (trim(left_) == trim(pattern)) then
        str = trim(right_(len(left_)+1:))
        exit
      end if
    else
      str = "0"
      exit
    end if
  end do
end function rest_str


!******************************************************************************** 
!
!       Extract first string from the larger string
!
!******************************************************************************** 
function extract_word(str) result(str_res)
  character(len=*), intent(in)  :: str
  character(len=:), allocatable :: str_res
  !local variables
  character(len=:), allocatable :: temp
  integer                       :: pos
  
  temp = trim(lowercase(str))
  pos = scan(temp, " ")
  
  if (pos > 0) then
    str_res = trim(temp(1:pos-1))
  else
    str_res = temp 
  end if
end function extract_word


!******************************************************************************** 
!
!       Read integer array from the string
!
!       Delimiter del is used to separate integer elements
!
!******************************************************************************** 
function int_arr_from_string(str, del) result(x)
  character(len=*), intent(in)           :: str
  character(len=1), optional, intent(in) :: del
  integer, allocatable                   :: x(:)
  !local variables
  integer, parameter                     :: maxlen = 1000
  integer                                :: n, posl, posr
  integer                                :: arr(maxlen)
  character(len=1)                       :: del_
  character(len=:), allocatable          :: str_

  if (present(del)) then
    del_ = del
  else 
    del_ = " "
  end if

  n = 0
  posl = 0

  str_ = trim(str)

  do 
    posr = posl + scan(str_(posl+1:), del_)
    if (posr-posl > 1) then
      n = n + 1
      read(str_(posl+1:posr-1), *) arr(n)
      posl = posr
    else if (posr-posl == 1) then
      posl = posr
      cycle
    else 
      n = n + 1
      read(str_(posl+1:), *) arr(n)
      exit
    end if
  end do

  if (allocated(x)) deallocate(x)
  allocate(x(n))
  x  = arr(1:n)
end function int_arr_from_string

end module string
