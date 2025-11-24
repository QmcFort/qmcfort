! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module test_regression

use regression

use constants, only: wp
use fruit, only: assert_equals, run_test_case, test_set_summary

implicit none

private 
public :: test_regression_set

contains

subroutine test_regression_set
  call run_test_case(test_quadratic_fit_simple, "Check quadratic fit function")
  call run_test_case(test_linear_fit_w, "Check linear fit with uncertainties in y")
  call run_test_case(test_polyfit_quad, "Check polyfit (n=2) using linear fit routine")

  call test_set_summary("src/regression.f90")
end subroutine test_regression_set

subroutine test_quadratic_fit_simple
    integer, parameter  :: n = 4
    real(wp), parameter :: tol = 0.0001_wp
    real(wp)            :: x(n), y(n)
    real(wp)            :: coeff(2), dcoeff(2), rsq
    real(wp)            :: coeff_(2), dcoeff_(2), rsq_
    
    x = [0.0_wp, 1.0_wp, 2.0_wp, 3.0_wp]
    y = [0.95_wp, 4.10_wp, 12.80_wp, 27.70_wp]
    
    coeff = [2.9627551_wp, 1.0178571_wp]
    dcoeff = 0.0_wp
    rsq = 0.9999444_wp

    call linear_fit(x**2, y, coeff_, dcoeff_, rsq_)

    call assert_equals(coeff, coeff_, 2, tol)
    call assert_equals(rsq, rsq_, tol)
end subroutine test_quadratic_fit_simple

subroutine test_linear_fit_w
    integer, parameter  :: n = 6
    real(wp), parameter :: tol = 0.2_wp
    real(wp)            :: x(n), y(n), dy(n)
    real(wp)            :: coeff(2), rsq
    real(wp)            :: coeff_(2), dcoeff_(2), rsq_
    
    x = [1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp, 6.0_wp]
    y = [1.04_wp, 3.01_wp, 4.97_wp, 6.99_wp, 9.00_wp, 11.04_wp]
    dy = [0.06_wp, 0.04_wp, 0.05_wp, 0.03_wp, 0.02_wp, 0.05_wp]
    
    coeff = [2.0_wp, -1.0_wp]

    call linear_fit(x, y, coeff_, dcoeff_, rsq_, dy)

    call assert_equals(coeff, coeff_, 2, tol)
end subroutine test_linear_fit_w

subroutine test_polyfit_quad
    integer, parameter  :: n = 4
    integer             :: order
    real(wp), parameter :: tol = 0.001_wp
    real(wp)            :: x(n), y(n)
    real(wp)            :: coeff(2), dcoeff(2), rsq
    real(wp)            :: coeff_(2), dcoeff_(2), rsq_
    
    order = 2
    x = [1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp]
    y = [1.98_wp, 12.03_wp, 26.03_wp, 48.05_wp]
    
    call linear_fit(x**2, y, coeff, dcoeff, rsq)
    call polyfit(x, y, order, coeff_, dcoeff_, rsq_)

    call assert_equals(coeff, coeff_, 2, tol)
    call assert_equals(dcoeff, dcoeff_, 2, tol)
    call assert_equals(rsq, rsq_, tol)
end subroutine test_polyfit_quad

end module test_regression