! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module test_standalone

use standalone

use constants, only: sp, dp, wp
use fruit, only: assert_equals, run_test_case, test_set_summary

implicit none

private
public :: test_standalone_set

contains

subroutine test_standalone_set
  call run_test_case(test_trace_zero_i, "test_trace_zero_i should give 0 trace")
  call run_test_case(test_trace_from_1_to_n_i, "test_trace_from_1_to_n should give n*(n+1)/2")
  call run_test_case(test_factorial_arg5, "factorial with argument 5")
  call run_test_case(test_factorial_arg0, "factorial with argument 0")
  call run_test_case(test_factorial_rec_arg5, "recursive factorial with argument 5")
  call run_test_case(test_factorial_rec_arg0, "recursive factorial with argument 0")
  call run_test_case(test_factorial_consistency, "check consistency of the factorials") 
  call run_test_case(test_dfactorial_arg5, "double factorial with argument 5")
  call run_test_case(test_dfactorial_arg6, "double factorial with argument 6")
  call run_test_case(test_dfactorial_arg0, "double factorial with argument 0")
  call run_test_case(test_dfactorial_rec_arg5, "recursive double factorial with argument 5")
  call run_test_case(test_dfactorial_rec_arg6, "recursive double factorial with argument 6")
  call run_test_case(test_dfactorial_rec_arg0, "recursive double factorial with argument 0")
  call run_test_case(test_dfactorial_consistency, "check consistency of the double factorials") 
  call run_test_case(test_trans_symm_i_trivial, "transpose lower part of the symmetric integer matrix")
  call run_test_case(test_trans_symm_i_up_trivial, "transpose upper part of the symmetric integer matrix")
  call run_test_case(test_trans_symm_c_trivial, "transpose lower part of the hermitian complex matrix")
  call run_test_case(test_extract_matrix_cols, "test extracting eigenvectors and eigenrows using mask")
  call run_test_case(test_phase_cmplx, "test phases obtained from complex numbers")
  call run_test_case(test_phase_real, "test 'phases' obtained from real numbers")
  call run_test_case(test_near_zero_real_dprec, "test near_zero for dprec real numbers")
  call run_test_case(test_near_zero_cmplx_sprec, "test near_zero for sprec complex numbers")

  call test_set_summary("src/standalone.f90")
end subroutine test_standalone_set

subroutine test_trace_zero_i
  integer, parameter :: n = 3
  integer            :: matrix(n,n)
  integer            :: expected
  integer            :: result_
  !
  expected = 0
  matrix = 0
  !
  result_ = trace(matrix)
  call assert_equals(expected, result_)
end subroutine test_trace_zero_i

subroutine test_trace_from_1_to_n_i
  integer, parameter :: n = 10
  integer            :: matrix(n,n)
  integer            :: expected
  integer            :: result_
  integer            :: i
  !
  expected = n * (n+1) / 2
  matrix = 0
  do i = 1, n
    matrix(i,i) = i
  end do
  !
  result_ = trace(matrix)
  call assert_equals(expected, result_)
end subroutine test_trace_from_1_to_n_i

subroutine test_factorial_arg5
  integer :: expected
  integer :: result_
  !
  expected = 120
  result_  = factorial(5)
  !
  call assert_equals(expected, result_)
end subroutine test_factorial_arg5

subroutine test_factorial_arg0
  integer :: expected
  integer :: result_
  !
  expected = 1
  result_  = factorial(0)
  !
  call assert_equals(expected, result_)
end subroutine test_factorial_arg0

subroutine test_factorial_rec_arg5
  integer :: expected
  integer :: result_
  !
  expected = 120
  result_  = factorial_rec(5)
  !
  call assert_equals(expected, result_)
end subroutine test_factorial_rec_arg5

subroutine test_factorial_rec_arg0
  integer :: expected
  integer :: result_
  !
  expected = 1
  result_  = factorial_rec(0)
  !
  call assert_equals(expected, result_)
end subroutine test_factorial_rec_arg0

subroutine test_factorial_consistency
  integer :: n
  integer :: expected
  integer :: result_
  !
  n = 10
  expected = factorial(n)
  result_ = n * factorial_rec(n-1)
  !
  call assert_equals(expected, result_)
end subroutine test_factorial_consistency

subroutine test_dfactorial_arg5
  integer :: expected
  integer :: result_
  !
  expected = 15
  result_  = double_factorial(5)
  !
  call assert_equals(expected, result_)
end subroutine test_dfactorial_arg5

subroutine test_dfactorial_arg6
  integer :: expected
  integer :: result_
  !
  expected = 48
  result_  = double_factorial(6)
  !
  call assert_equals(expected, result_)
end subroutine test_dfactorial_arg6

subroutine test_dfactorial_arg0
  integer :: expected
  integer :: result_
  !
  expected = 1
  result_  = double_factorial(0)
  !
  call assert_equals(expected, result_)
end subroutine test_dfactorial_arg0

subroutine test_dfactorial_rec_arg5
  integer :: expected
  integer :: result_
  !
  expected = 15
  result_  = double_factorial_rec(5)
  !
  call assert_equals(expected, result_)
end subroutine test_dfactorial_rec_arg5

subroutine test_dfactorial_rec_arg6
  integer :: expected
  integer :: result_
  !
  expected = 48
  result_  = double_factorial_rec(6)
  !
  call assert_equals(expected, result_)
end subroutine test_dfactorial_rec_arg6

subroutine test_dfactorial_rec_arg0
  integer :: expected
  integer :: result_
  !
  expected = 1
  result_  = double_factorial_rec(0)
  !
  call assert_equals(expected, result_)
end subroutine test_dfactorial_rec_arg0

subroutine test_dfactorial_consistency
  integer :: n
  integer :: expected
  integer :: result_
  !
  n = 10
  expected = double_factorial(n)
  result_ = n * double_factorial_rec(n-2)
  !
  call assert_equals(expected, result_)
end subroutine test_dfactorial_consistency

subroutine test_trans_symm_i_trivial
  integer, parameter :: n = 3
  integer            :: expected(n,n)
  integer            :: result_(n,n)   
  !
  expected = reshape([1,2,3,2,4,5,3,5,6], shape=[n, n])
  result_  = reshape([1,2,3,0,4,5,0,0,6], shape=[n, n])
  !
  call trans_symm(result_)
  call assert_equals(expected, result_, n, n)
end subroutine test_trans_symm_i_trivial

subroutine test_trans_symm_i_up_trivial
  integer, parameter :: n = 3
  integer            :: expected(n,n)
  integer            :: result_(n,n)  
  !
  expected  = reshape([1,2,4,2,3,5,4,5,6], shape=[n, n])
  result_   = reshape([1,0,0,2,3,0,4,5,6], shape=[n, n])
  !
  call trans_symm(result_, "up")
  call assert_equals(expected, result_, n, n)
end subroutine test_trans_symm_i_up_trivial

subroutine test_trans_symm_c_trivial
  integer, parameter :: n = 2
  complex(dp)        :: expected(n,n)
  complex(dp)        :: result_(n,n)
  !
  expected = reshape([(1.0, 1.0),(1.0, -1.0),(1.0, 1.0),(2.0, 1.0)], shape=[n, n])
  result_  = reshape([(1.0, 1.0),(1.0, -1.0),(0.0, 0.0),(2.0, 1.0)], shape=[n, n])
  !
  call trans_symm(result_)
  call assert_equals(expected, result_, n, n)
end subroutine test_trans_symm_c_trivial

subroutine test_extract_matrix_cols
  integer, parameter    :: n=5, m=3
  real(sp), parameter   :: tol=1.0e-06_sp
  real(sp)              :: u(n,n)
  real(sp)              :: e(n)
  real(sp)              :: expected_1(n,m)
  real(sp)              :: expected_2(m)
  real(sp), allocatable :: result_1(:,:)
  real(sp), allocatable :: result_2(:)
  logical               :: mask(n)
  
  !
  e = [1.0, 0.0, 0.0, 3.0, 1.0]
  mask = [.true., .false., .false., .true., .true.]
  u(:,1) = [1.0, 2.0, 3.0, 4.0, 5.0]
  u(:,2) = [3.0, 5.0, 6.0, 8.0, 9.0]
  u(:,3) = [0.0, 1.0, 9.0, 4.0, 3.0]
  u(:,4) = [2.0, 2.0, 2.0, 1.0, 1.0]
  u(:,5) = [1.0, 2.0, 1.0, 1.0, 1.0]
  !
  expected_1(:,1) = u(:,1)
  expected_1(:,2) = u(:,4)
  expected_1(:,3) = u(:,5)
  expected_2 = [1.0, 3.0, 1.0]
  !
  call extract_matrix_cols(u, e, mask, result_1, result_2)
  !
  call assert_equals(expected_1, result_1, n, m, tol)
  call assert_equals(expected_2, result_2, m, tol)
end subroutine test_extract_matrix_cols

subroutine test_phase_cmplx
  integer, parameter     :: n = 8
  real(wp), parameter    :: tol = 1.0e-04_wp
  complex(wp), parameter :: imag = (0.0_wp, 1.0_wp)
  real(wp)               :: phases(n)
  real(wp)               :: result_(n)
  real(wp)               :: expected(n)
  !
  phases = [0.0_wp, pi/4.0_wp, pi, 2.0_wp*pi, 3.0_wp*pi, -pi/2.0_wp, -pi, -2.0_wp*pi]
  expected = [0.0_wp, pi/4.0_wp, pi, 0.0_wp, pi, -pi/2.0_wp, -pi, 0.0_wp]
  result_ = phase(exp(imag * phases))
  !
  call assert_equals(expected, result_, n, tol)
end subroutine test_phase_cmplx

subroutine test_phase_real
  integer, parameter  :: n = 8
  real(wp), parameter :: tol = 1.0e-04_wp
  real(wp)            :: phases(n)
  real(wp)            :: result_(n)
  real(wp)            :: expected(n)
  !
  phases = [0.0_wp, pi/4.0_wp, pi, 2.0_wp*pi, 3.0_wp*pi, -pi/2.0_wp, -pi, -2.0_wp*pi]
  expected = [pi, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, pi, pi, pi]
  result_ = phase(phases)
  !
  call assert_equals(expected, result_, n, tol)
end subroutine test_phase_real

subroutine test_near_zero_real_dprec
  integer, parameter  :: n = 3
  real(dp), parameter :: tol = 0.000001_dp
  real(dp)            :: x(n), y(n)
  logical             :: expected(n), result_(n)

  x = [1.0_dp, -3.0_dp, 10.0_dp]
  y = [1.0000001_dp, -3.001_dp, 10.0000000000000000000001_dp]

  expected = [.true., .false., .true.]
  result_ = near_zero(x, y , tol=tol)

  call assert_equals(expected, result_, n)
end subroutine test_near_zero_real_dprec

subroutine test_near_zero_cmplx_sprec
  integer, parameter :: n = 3
  complex(sp)        :: x(n), y(n)
  logical            :: expected(n), result_(n)

  x = [(1.0_sp, 0.0_sp), (-3.0_sp, 1.0_sp), (10.0_sp, -2.0_sp)]
  y = [(1.0000000001_sp, 0.0_sp), (-3.001_sp, 1.0_sp), (10.1_sp, -2.1_sp)]

  expected = [.true., .false., .false.]
  result_ = near_zero(x, y)

  call assert_equals(expected, result_, n)
end subroutine test_near_zero_cmplx_sprec
end module test_standalone

