! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module test_boys_functions

use boys_functions

use constants, only: wp
use fruit, only: assert_equals, run_test_case, test_set_summary

implicit none

private
public :: test_boys_functions_set

contains

subroutine test_boys_functions_set
  call run_test_case(test_boys_functions_small_an, "compare small and analytic boysfunctions in range [5 15]")
  call run_test_case(test_boys_functions_up_an, "compare upward rec and analytic boysfunctions in range [15 35]")
  call run_test_case(test_boys_functions_large_an, "compare upward rec and analytic boysfunctions in range [35 50]")
  call run_test_case(test_boys_functions_rec_large, "compare boysfun with boysfun_up for large argument")
  call run_test_case(test_boys_functions_rec_med, "compare boysfun with boysfun_up for medium size argument")
  call run_test_case(test_boys_functions_table_up, "compare boys functions table with upward recursion with analytic one")
  call run_test_case(test_boys_functions_table_down, &
                     "compare boys functions table with downward recursion with analytic one")

  call test_set_summary("src/gto/test_boys_functions.f90")
end subroutine test_boys_functions_set

subroutine test_boys_functions_small_an
  integer, parameter  :: n=50, m=11
  real(wp), parameter :: tol = 1.0e-04_wp
  integer             :: i, j
  real(wp)            :: x0, x, dx
  real(wp)            :: boys1(n,m)
  real(wp)            :: boys2(n,m)
  ! 
  dx = 1.0_wp
  x0 = 5.0_wp
  !range [5.0, 15.0]
  do j = 1, m
    x = x0 + (j-1) * dx
    do i = 1, n
      boys1(i,j) = boysfun_small(x, i)
      boys2(i,j) = boysfun_an(x, i)  
    end do
  end do
  !
  call assert_equals(boys1, boys2, n, m, tol)
end subroutine test_boys_functions_small_an

subroutine test_boys_functions_up_an
  integer, parameter  :: n=50, m=11
  real(wp), parameter :: tol = 1.0e-04_wp
  integer             :: i, j
  real(wp)            :: x0, x, dx
  real(wp)            :: boys1(n,m)
  real(wp)            :: boys2(n,m)
  ! 
  dx = 2.0_wp
  x0 = 15.0_wp
  !range [15.0, 35.0]
  do j = 1, m
    x = x0 + (j-1) * dx
    do i = 1, n
      boys1(i,j) = boysfun_up(x, i)
      boys2(i,j) = boysfun_an(x, i)  
    end do
  end do
  !
  call assert_equals(boys1, boys2, n, m, tol)
end subroutine test_boys_functions_up_an

subroutine test_boys_functions_large_an
  integer, parameter  :: n=50, m=16
  real(wp), parameter :: tol = 1.0e-06_wp
  integer             :: i, j
  real(wp)            :: x0, x, dx
  real(wp)            :: boys1(n,m)
  real(wp)            :: boys2(n,m)
  ! 
  dx = 2.0_wp
  x0 = 35.0_wp
  !range [35.0, 65.0]
  do j = 1, m
    x = x0 + (j-1) * dx
    do i = 1, n
      boys1(i,j) = boysfun_large(x, i)
      boys2(i,j) = boysfun_an(x, i)  
    end do
  end do
  !
  call assert_equals(boys1, boys2, n, m, tol)
end subroutine test_boys_functions_large_an

subroutine test_boys_functions_rec_large
  integer, parameter   :: n = 50
  real(wp), parameter  :: tol = 1.0e-06_wp
  integer              :: i
  real(wp)             :: x
  real(wp)             :: boys1(n), boys2(n)
  !
  x = 125.0_wp
  do i = 1, n
    boys1(i) = boysfun_large(x, i)
    boys2(i) = boysfun_up(x, i)
  end do
  !
  call assert_equals(boys1, boys2, n, tol)
end subroutine test_boys_functions_rec_large

subroutine test_boys_functions_rec_med
  integer, parameter   :: n = 50
  real(wp), parameter  :: tol = 1.0e-06_wp
  integer              :: i
  real(wp)             :: x
  real(wp)             :: boys1(n), boys2(n)
  !
  x = 18.0_wp
  do i = 1, n
    boys1(i) = boysfun(x, i)
    boys2(i) = boysfun_up(x, i)
  end do
  !
  call assert_equals(boys1, boys2, n, tol)
end subroutine test_boys_functions_rec_med

subroutine test_boys_functions_table_down
  integer, parameter    :: max_=50
  real(wp), parameter   :: tol = 1.0e-06_wp
  integer               :: n
  real(wp)              :: x
  real(wp), allocatable :: boys1(:), boys2(:)
  !
  allocate(boys1(0:max_))
  x = 0.8_wp
  do n = 0, max_
    boys1(n) = boysfun(x, n)
  end do  
  call boysfun_table(max_, x, boys2)
  call assert_equals(boys1, boys2, max_+1, tol)
end subroutine test_boys_functions_table_down

subroutine test_boys_functions_table_up
  integer, parameter    :: max_=50
  real(wp), parameter   :: tol = 1.0e-06_wp
  integer               :: n
  real(wp)              :: x
  real(wp), allocatable :: boys1(:), boys2(:)
  !
  allocate(boys1(0:max_))
  x = 22.9_wp
  do n = 0, max_
    boys1(n) = boysfun(x, n)
  end do  
  call boysfun_table(max_, x, boys2)
  call assert_equals(boys1, boys2, max_+1, tol)
end subroutine test_boys_functions_table_up

end module test_boys_functions
