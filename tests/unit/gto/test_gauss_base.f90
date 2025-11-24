! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module test_gauss_base

use gauss_base

use constants, only: wp
use fruit, only: assert_equals, run_test_case, test_set_summary

implicit none

private
public :: test_gauss_base_set

contains

subroutine test_gauss_base_set
  call run_test_case(test_cart_indx, "Test ordering of catresian GTOs")
  call run_test_case(test_harm_cart_trafo_d_shell, "Test harmonic-cartesian GTO transformation for l=2")
  call run_test_case(test_harm_cart_trafo_f_shell, "Test harmonic-cartesian GTO transformation for l=3")

  call test_set_summary("src/gto/gauss_base.f90")
end subroutine test_gauss_base_set

subroutine test_cart_indx
  integer, parameter  :: l=2
  integer             :: expected(2)
  integer             :: result_(2)

  expected = [1, 4]
  
  result_(1) = get_cart_indx(2, 0, 0)
  result_(2) = get_cart_indx(1, 0, 1)

  call assert_equals(expected, result_, 2)
end subroutine test_cart_indx

subroutine test_harm_cart_trafo_d_shell
  integer, parameter    :: l=2, size1=6, size2=5
  real(wp), allocatable :: expected(:,:)
  real(wp), allocatable :: result_(:,:)
  real(wp), parameter   :: tol = 1.0e-06_wp

  allocate(expected(size1, size2))
  expected = 0.0_wp
  expected(2,5) = sqrt(3.0_wp)
  expected(5,4) = sqrt(3.0_wp)
  expected(1,3) = -0.5_wp
  expected(3,3) = -0.5_wp
  expected(6,3) = 1.0_wp
  expected(4,2) = sqrt(3.0_wp)
  expected(1,1) = sqrt(3.0_wp) / 2.0_wp
  expected(3,1) = -sqrt(3.0_wp) / 2.0_wp

  result_ = harm_cart_trafo(l)

  call assert_equals(size1, size(result_, 1))
  call assert_equals(size2, size(result_, 2))
  call assert_equals(expected, result_, size1, size2, tol)
end subroutine test_harm_cart_trafo_d_shell

subroutine test_harm_cart_trafo_f_shell
  integer, parameter    :: l=3, size1=10, size2=7
  real(wp), allocatable :: expected(:,:)
  real(wp), allocatable :: result_(:,:)
  real(wp), parameter   :: tol = 1.0e-06_wp

  allocate(expected(size1, size2))
  expected = 0.0_wp
  expected(2,7) = 3.0_wp * sqrt(10.0_wp) / 4.0_wp
  expected(4,7) = -sqrt(10.0_wp) / 4.0_wp
  expected(6,6) = sqrt(15.0_wp)
  expected(2,5) = -sqrt(6.0_wp) / 4.0_wp
  expected(4,5) = -sqrt(6.0_wp) / 4.0_wp
  expected(9,5) = sqrt(6.0_wp)
  expected(5,4) = -3.0_wp / 2.0_wp
  expected(7,4) = -3.0_wp / 2.0_wp
  expected(10,4) = 1.0_wp
  expected(1,3) = -sqrt(6.0_wp) / 4.0_wp
  expected(3,3) = -sqrt(6.0_wp) / 4.0_wp
  expected(8,3) = sqrt(6.0_wp)
  expected(5,2) = sqrt(15.0_wp) / 2.0_wp
  expected(7,2) = -sqrt(15.0_wp) / 2.0_wp
  expected(1,1) = sqrt(10.0_wp) / 4.0_wp 
  expected(3,1) = -3.0_wp * sqrt(10.0_wp) / 4.0_wp 

  result_ = harm_cart_trafo(l)

  call assert_equals(size1, size(result_, 1))
  call assert_equals(size2, size(result_, 2))
  call assert_equals(expected, result_, size1, size2, tol)
end subroutine test_harm_cart_trafo_f_shell

end module test_gauss_base