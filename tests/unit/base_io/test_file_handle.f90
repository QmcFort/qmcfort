! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module test_file_handle

use file_handle, only: FileHandle
use fruit, only: assert_equals, run_test_case, test_set_summary

implicit none

private
public :: test_file_handle_set

contains

subroutine test_file_handle_set
  call run_test_case(test_file_handle_init_exist, "test_file_handle_init_exist should be false if used for non existing file")
  call run_test_case(test_file_handle_open_close_opened, &
                     "fopen should be .true. once the file is opened, and .false. when closed")
  call run_test_case(test_file_handle_open_close_exist,  &
                     "fexist should be .true. once the file is opened, and .false. when deleted")
  call run_test_case(test_file_handle_units_two_files, "Difference between file units of two subsequently opened files is 1")
  call run_test_case(test_file_handle_read_integer_bin, "write and read integer value from the binary file")
  call run_test_case(test_file_handle_read_reals_bin, &
                     "write and read real value from the binary file")
  call run_test_case(test_file_handle_read_wrong_reclen,  &
                     "write direct acc bin file with actual reclen and read with assumed reclen")

  call test_set_summary("src/base_io/file_handle.f90")
end subroutine test_file_handle_set

subroutine test_file_handle_init_exist
  type(FileHandle)  :: temp_file
  logical           :: expected
  
  expected = .false.
  temp_file = FileHandle("file")
  
  call assert_equals(expected, temp_file%exists())
end subroutine test_file_handle_init_exist

subroutine test_file_handle_open_close_opened
  type(FileHandle)  :: fh
  logical           :: expected_1, expected_2
  logical           :: result_1, result_2
  
  expected_1 = .true.
  expected_2 = .false. 
  
  fh = FileHandle("test_file")

  call fh%open(status="new", action="write")
  result_1 = fh%is_open()

  call fh%close(status="delete")
  result_2 = fh%is_open()
  
  call assert_equals(expected_1, result_1)
  call assert_equals(expected_2, result_2)
end subroutine test_file_handle_open_close_opened

subroutine test_file_handle_open_close_exist
  type(FileHandle)  :: fh
  logical           :: expected_1, expected_2
  logical           :: result_1, result_2
  
  expected_1 = .true.
  expected_2 = .false. 

  fh = FileHandle("test_file")

  call fh%open(status="new", action="write")
  result_1 = fh%exists()

  call fh%close(status="delete")
  result_2 = fh%exists()
  
  call assert_equals(expected_1, result_1)
  call assert_equals(expected_2, result_2)
end subroutine test_file_handle_open_close_exist

subroutine test_file_handle_units_two_files
  type(FileHandle) :: fh1, fh2
  logical          :: expected, result_
  
  expected = .true.

  fh1 = FileHandle("test_file")
  fh2 = FileHandle("test_file2")
  call fh1%open()
  call fh2%open()

  result_ = fh2%funit /= fh1%funit

  call fh1%close(status="delete")
  call fh2%close(status="delete")
  
  call assert_equals(expected, result_)
end subroutine test_file_handle_units_two_files

subroutine test_file_handle_read_integer_bin
  type(FileHandle)  :: fh
  integer           :: var_read, var_write
  
  var_read = 1234

  fh = FileHandle("file")
  call fh%open(status="new", access="direct", form="unformatted", recl=4)
  write(unit=fh%funit, rec=1) var_read
  read(unit=fh%funit, rec=1) var_write
  call fh%close(status="delete")
  
  call assert_equals(var_read, var_write)
end subroutine test_file_handle_read_integer_bin

subroutine test_file_handle_read_reals_bin
  use constants, only: wp
  type(FileHandle)  :: temp_file
  real(wp)          :: expected, result_
  integer           :: i, reclen, num_records=10, irec=7
  
  expected = real(irec, wp)
  temp_file = FileHandle("no_name")
  inquire(iolength=reclen) expected
  call temp_file%open(status="new", access="direct", form="unformatted", recl=reclen)

  do i = 1, num_records
    write(unit=temp_file%funit, rec=i) real(i, wp)
  end do
  read(unit=temp_file%funit, rec=irec) result_

  call temp_file%close(status="delete")
  
  call assert_equals(expected, result_)
end subroutine test_file_handle_read_reals_bin

subroutine test_file_handle_read_wrong_reclen
  use constants, only: wp
  type(FileHandle)  :: temp_file
  integer, parameter :: n=10, m=100
  integer            :: reclen
  integer            :: i
  real(wp)           :: expected(n), result_(n), array(m)
  
  !initialize data
  do i = 1, n
    expected(i) = real(i, wp)
  end do
  do i = 1, m
    array(i) = real(i, wp)
  end do

  !write data - larger reclen will be used
  inquire(iolength=reclen) array
  temp_file = FileHandle("no_name")
  call temp_file%open(status="new", access="direct", form="unformatted", recl=reclen)    
  write(unit=temp_file%funit, rec=1) expected
  write(unit=temp_file%funit, rec=2) array
  call temp_file%close()

  !read data with wrong reclen - smaller reclen of the smaller array is used
  inquire(iolength=reclen) result_
  call temp_file%open(status="old", access="direct", form="unformatted", recl=reclen)
  read(unit=temp_file%funit, rec=1) result_
  call temp_file%close(status="delete")
  
  call assert_equals(expected, result_, n)
end subroutine test_file_handle_read_wrong_reclen

end module test_file_handle