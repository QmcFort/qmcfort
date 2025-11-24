! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module test_qmcfort_in

use mpi
use qmcfort_in

use constants, only: wp
use file_handle, only: FileHandle
use fruit, only: assert_equals, run_test_case, test_set_summary
use qmcfort_io, only: io

implicit none

private
public :: test_qmcfort_in_set

contains

subroutine test_qmcfort_in_set
  call run_test_case(test_readin_int, "Test qmcfort_in read integer value")
  call run_test_case(test_readin_real_1d, "Test qmcfort_in read of the 1d real array")
  call run_test_case(test_readin_real_1d_alloc, "Test qmcfort_in read of the 1d real allocatable array")

  call test_set_summary("src/base_io/qmcfort_in.f90")
end subroutine test_qmcfort_in_set

subroutine test_readin_int
  integer, parameter  :: size_ = 4
  integer             :: expected
  integer             :: result_

  if (.not. comm_world%initialized) call comm_world%init()

  io%qmcfort_in = FileHandle("qmcfort_in")
  call io%qmcfort_in%open(status="replace", action="write")
  write(io%qmcfort_in%funit,*) "hf = .true.            "
  write(io%qmcfort_in%funit,*) "test = 5               "
  write(io%qmcfort_in%funit,*) "tau = 10               "
  write(io%qmcfort_in%funit,*) "nblocks = -100         "
  write(io%qmcfort_in%funit,*) "diis = .false.         "
  call io%qmcfort_in%close()
    
  expected = 5
  call add_input("test", result_)

  call assert_equals(expected, result_)
  
  call io%qmcfort_in%delete()
end subroutine test_readin_int

subroutine test_readin_real_1d
  integer, parameter    :: size_ = 4
  real(wp), parameter   :: tol = 1.0e-06_wp
  real(wp)              :: expected(size_)
  real(wp), allocatable :: result_(:)
  
  if (.not. comm_world%initialized) call comm_world%init

  io%qmcfort_in = FileHandle("qmcfort_in")
  call io%qmcfort_in%open(status="replace", action="write")
  write(io%qmcfort_in%funit,*) "hf = .true.            "
  write(io%qmcfort_in%funit,*) "test = 0.1 0.2 .3 0.4  "
  write(io%qmcfort_in%funit,*) "tau = 10               "
  write(io%qmcfort_in%funit,*) "nblocks = -100         "
  write(io%qmcfort_in%funit,*) "diis = .false.         "
  call io%qmcfort_in%close()
  
  expected = [0.1_wp, 0.2_wp, 0.3_wp, 0.4_wp]
  call add_input("test", result_)

  call assert_equals(expected, result_, size_, tol)
  
  call io%qmcfort_in%delete()
end subroutine test_readin_real_1d

subroutine test_readin_real_1d_alloc
  integer, parameter    :: size_ = 4
  real(wp), parameter   :: tol = 1.0e-06_wp
  real(wp)              :: expected(size_)
  real(wp), allocatable :: result_(:)
  
  if (.not. comm_world%initialized) call comm_world%init
  
  io%qmcfort_in = FileHandle("qmcfort_in")
  call io%qmcfort_in%open(status="replace", action="write")
  write(io%qmcfort_in%funit,*) "hf = .true.            "
  write(io%qmcfort_in%funit,*) "test = 0.1 0.2 .3 0.4  "
  write(io%qmcfort_in%funit,*) "tau = 10               "
  write(io%qmcfort_in%funit,*) "nblocks = -100         "
  write(io%qmcfort_in%funit,*) "diis = .false.         "
  call io%qmcfort_in%close()
  
  expected = [0.1_wp, 0.2_wp, 0.3_wp, 0.4_wp]
  call add_input("test", result_)

  call assert_equals(expected, result_, size_, tol)
  
  call io%qmcfort_in%delete()
end subroutine test_readin_real_1d_alloc

end module test_qmcfort_in