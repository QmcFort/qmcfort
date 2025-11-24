! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module test_mc_descriptor

use mc_descriptor

use constants, only: wp
use fruit, only: assert_equals, run_test_case, test_set_summary

implicit none
 
private 
public :: test_mc_descriptor_set

contains

subroutine test_mc_descriptor_set
  call run_test_case(test_mcdes_increment, "test step and block variables in mc_descriptor")
  call run_test_case(test_mcdes_increment_after_eq, "test equilibration status after equilibration in mc_descriptor")
  call run_test_case(test_mcdes_tau, "test tau and sqmtau in mc_descriptor")
  call run_test_case(test_mcdes_nsteps, "test nsteps and time functions in mc_descriptor")

  call test_set_summary("src/mc_descriptor.f90")
end subroutine test_mc_descriptor_set

subroutine test_mcdes_increment
  type(McDescriptor) :: mcdes
  integer, parameter :: n = 2410
  integer            :: i
  integer            :: expected_step, expected_block
  logical            :: expected_eq
  
  mcdes = McDescriptor(nwalkers=1000, tau=0.01_wp, steps_per_block=100, nblocks=100, eqblocks=100, &
                       reorth=1, sample=1, samplex=1)

  do i = 1, n
    call mcdes%increment()
  end do

  expected_step = n
  expected_block = 25
  expected_eq = .false.

  call assert_equals(expected_step, mcdes%step)
  call assert_equals(expected_block, mcdes%block)
  call assert_equals(expected_eq, mcdes%is_eq())
end subroutine test_mcdes_increment

subroutine test_mcdes_increment_after_eq
  type(McDescriptor) :: mcdes
  integer, parameter :: n = 312, bblock=32
  integer            :: i
  integer            :: expected_step, expected_block
  logical            :: expected_eq
  
  mcdes = McDescriptor(nwalkers=1000, tau=0.01_wp, steps_per_block=10, nblocks=80, eqblocks=20, &
                       reorth=1, sample=1, samplex=1)

  do i = 1, n
    call mcdes%increment()
  end do

  expected_step = n
  expected_block = bblock
  expected_eq = .true.


  call assert_equals(expected_step, mcdes%step)
  call assert_equals(expected_block, mcdes%block)
  call assert_equals(expected_eq, mcdes%is_eq())
end subroutine test_mcdes_increment_after_eq

subroutine test_mcdes_tau
  type(McDescriptor)  :: mcdes
  real(wp), parameter :: tol=1.0e-06_wp
  real(wp)            :: expected_tau
  complex(wp)         :: expected_isqtau
  
  expected_tau = 0.09_wp
  expected_isqtau = (0.00_wp, 0.3_wp)

  mcdes = McDescriptor(nwalkers=1000, tau=expected_tau, steps_per_block=100, nblocks=100, eqblocks=100, &
                       reorth=1, sample=1, samplex=1)

  call assert_equals(expected_tau, mcdes%tau, tol)
  call assert_equals(expected_isqtau, mcdes%isqtau, tol)
end subroutine test_mcdes_tau

subroutine test_mcdes_nsteps
  type(McDescriptor)  :: mcdes
  real(wp), parameter :: tol=1.0e-06_wp
  integer             :: expected_nsteps
  real(wp)            :: expected_time
  
  mcdes = McDescriptor(nwalkers=1000, tau=0.01_wp, steps_per_block=10, nblocks=80, eqblocks=20, &
                       reorth=1, sample=1, samplex=1)
  
  expected_nsteps = 1000
  expected_time = 10.0_wp

  call assert_equals(expected_nsteps, mcdes%steps_tot())
  call assert_equals(expected_time, mcdes%time_tot(), tol)
end subroutine test_mcdes_nsteps
end module test_mc_descriptor