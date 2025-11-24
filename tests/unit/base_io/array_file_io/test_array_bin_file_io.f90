! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module test_array_bin_file_io

use array_bin_file_io, only: ArrayBinFileIO
use constants, only: wp
use file_handle, only: FileHandle
use fruit, only: assert_equals, run_test_case, test_set_summary

implicit none

private
public :: test_array_bin_file_io_set

contains

subroutine test_array_bin_file_io_set
  call run_test_case(test_array_bin_file_io_2d_save_load, "save/load 2d array to/from binary file")

  call test_set_summary("src/base_io/array_file_io/array_bin_file_io.f90")
end subroutine test_array_bin_file_io_set

subroutine test_array_bin_file_io_2d_save_load
  integer, parameter          :: n=100, m=10, small_fac=10
  real(wp), parameter         :: tol = small_fac * tiny(1.0_wp)
  character(len=*), parameter :: fname = "test_file_bin"
  real(wp), allocatable       :: expected(:,:), result_(:,:)
  type(FileHandle)            :: fh
  type(ArrayBinFileIO)        :: array_io
  
  allocate(expected(n,m))
  call random_number(expected)
  
  call array_io%save(fname, expected)
  call array_io%load(fname, result_)
  
  fh = FileHandle(fname)
  call fh%delete()

  call assert_equals(expected, result_, n, m, tol)
end subroutine test_array_bin_file_io_2d_save_load

end module test_array_bin_file_io
