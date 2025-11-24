! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module test_diis_mixer

use diis_mixer
use fruit, only: assert_equals, run_test_case, test_set_summary

implicit none

private
public :: test_diis_mixer_set

contains

subroutine test_diis_mixer_set
  call run_test_case(test_diis_simple, "Simple test for diis scheme")

  call test_set_summary("src/diis/diis_mixer.f90")
end subroutine test_diis_mixer_set

subroutine test_diis_simple
  use constants, only: wp
  integer, parameter    :: n=3, m=2
  real(wp), parameter   :: tol = 1.0e-08_wp
  integer               :: i 
  type(DiisConfig)      :: diis_cfg
  type(DiisMixer)       :: diis
  real(wp), allocatable :: x(:,:), e(:,:)
  real(wp), allocatable :: expected(:), result_(:)


  diis_cfg = DiisConfig(start_sampling_step=1, start_solving_step=1)
  diis = DiisMixer(n, diis_cfg)
  
  allocate(x(n,m), e(n,m), expected(n), result_(n))
  x(:,1) = [1.0_wp, 0.0_wp, 0.0_wp]
  x(:,2) = [0.0_wp, 1.0_wp, 0.0_wp]
  e(:,1) = [1.0_wp, 0.0_wp, 0.0_wp]
  e(:,2) = [0.5_wp, 1.0_wp, 0.0_wp]
  expected = [0.6_wp, 0.4_wp, 0.0_wp]
  
  do i = 1, m
    call diis%update(x(:,i), e(:,i))
  end do

  call diis%solve(result_)

  call assert_equals(expected, result_, n, tol)
end subroutine test_diis_simple

end module test_diis_mixer
