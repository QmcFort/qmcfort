! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module test_array_numpy_file_io

use array_numpy_file_io, only: ArrayNumpyFileIO
use constants, only: wp
use file_handle, only: FileHandle
use fruit, only: assert_equals, run_test_case, test_set_summary

implicit none

private
public :: test_array_numpy_file_io_set

contains

subroutine test_array_numpy_file_io_set
  call run_test_case(test_array_np_file_io_3d_save_load, "save/load 3d array to/from numpy file")

  call test_set_summary("src/base_io/array_file_io/array_numpy_file_io.f90")
end subroutine test_array_numpy_file_io_set

subroutine test_array_np_file_io_3d_save_load
  integer, parameter          :: n=20, m=10, k=5, small_fac=10
  real(wp), parameter         :: tol = small_fac * tiny(1.0_wp)
  character(len=*), parameter :: fname = "test_file_np"
  complex(wp), allocatable    :: expected(:,:,:), result_(:,:,:)
  real(wp), allocatable       :: expected_r(:,:,:), expected_i(:,:,:)
  type(FileHandle)            :: fh
  type(ArrayNumpyFileIO)      :: array_io
  
  allocate(expected(n,m,k), expected_r(n,m,k), expected_i(n,m,k))
  call random_number(expected_r)
  call random_number(expected_i)
  expected = cmplx(expected_r, expected_i, kind=wp)
  
  call array_io%save(fname, expected)
  call array_io%load(fname, result_)

  fh = FileHandle(fname)
  call fh%delete()
  
  call assert_equals(reshape(expected, shape=[n*m*k]), reshape(result_, shape=[n*m*k]), n*m*k, tol)
end subroutine test_array_np_file_io_3d_save_load

end module test_array_numpy_file_io
