! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module test_array_file_io_factories

use array_file_io, only: ArrayFileIO
use array_io_factories, only: array_io_factory
use constants, only: wp
use file_handle, only: FileHandle
use fruit, only: assert_equals, run_test_case, test_set_summary

implicit none

private
public :: test_array_file_io_factories_set

contains

subroutine test_array_file_io_factories_set
  call run_test_case(test_array_io_factory, "array_io factory: save/load 3d numpy file")

  call test_set_summary("src/base_io/array_file_io/array_file_io_factories.f90")
end subroutine test_array_file_io_factories_set

subroutine test_array_io_factory
  integer, parameter              :: n=5, m=4, k=3, small_fac=10
  real(wp), parameter             :: tol = small_fac * tiny(1.0_wp)
  character(len=*), parameter     :: fname = "test_file_np"
  real(wp), allocatable           :: expected(:,:,:), result_(:,:,:)
  type(FileHandle)                :: fh
  class(ArrayFileIO), allocatable :: array_io
  
  allocate(expected(n,m,k))
  call random_number(expected)
  
  call array_io_factory("numpy", array_io)
  call array_io%save(fname, expected)
  call array_io%load(fname, result_)

  fh = FileHandle(fname)
  call fh%delete()
  
  call assert_equals(reshape(expected, shape=[n*m*k]), reshape(result_, shape=[n*m*k]), n*m*k, tol)
end subroutine test_array_io_factory

end module test_array_file_io_factories