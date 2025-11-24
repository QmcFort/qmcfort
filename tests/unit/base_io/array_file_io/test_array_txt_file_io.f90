! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module test_array_txt_file_io

use array_txt_file_io, only: ArrayTxtFileIO
use constants, only: wp
use file_handle, only: FileHandle
use fruit, only: assert_equals, run_test_case, test_set_summary

implicit none

private
public :: test_array_txt_file_io_set

contains

subroutine test_array_txt_file_io_set
  call run_test_case(test_array_txt_file_io_1d_read_write, "write and read 1d array in formatted txt file")

  call test_set_summary("src/base_io/array_file_io/array_txt_file_io.f90")
end subroutine test_array_txt_file_io_set

subroutine test_array_txt_file_io_1d_read_write
  integer, parameter    :: n = 100
  real(wp), parameter   :: tol = 1.0e-06_wp
  real(wp), allocatable :: expected(:), result_(:)
  type(FileHandle)      :: fh
  type(ArrayTxtFileIO)  :: array_io
  
  allocate(expected(n))
  call random_number(expected)
  
  fh = FileHandle("txt_file")

  call fh%open(status="replace", form="formatted", action="write")
  call array_io%write(fh%funit, expected)
  call fh%close()

  call fh%open(status="old", form="formatted", action="read")
  call array_io%read(fh%funit, result_)
  call fh%close(status="delete")
  
  call assert_equals(expected, result_, n, tol)
end subroutine test_array_txt_file_io_1d_read_write

end module   test_array_txt_file_io