! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module test_statistics

use mpi

use constants, only: wp
use dist_rng, only: DistRng
use fruit, only: assert_equals, run_test_case, test_set_summary
use lcg48_basic_rng, only: Lcg48BasicRng

implicit none

private
public :: test_statistics_set

contains

subroutine test_statistics_set
  call run_test_case(test_mean_simple, "test mean value of the simple array ")
  call run_test_case(test_mean_w_simple, "test weighted mean value of the simple array ")
  call run_test_case(test_mean_ax, "test mean value and ax argument of the random array ")
  call run_test_case(test_cov_matrix, "test covariance matrix of real weighted 3x3 sample")
  call run_test_case(test_cov_matrix_rc, "test covariance matrix of complex weighted 3x3 sample")
  call run_test_case(test_var_cov, "test variance of sum from the covariance")
  call run_test_case(test_histogram_simple, "test simple histogram")
  call run_test_case(test_histogram_auto, "test simple histogram without giving minmax values")

  call test_set_summary("stc/statistic.f90")
end subroutine test_statistics_set

subroutine test_mean_simple
  use statistics, only: mean
  integer, parameter  :: n=5
  real(wp), parameter :: tol=1.0e-06_wp
  real(wp)            :: x(n)
  real(wp)            :: expected, result_
  
  x = [1.0_wp, 0.0_wp, -1.0_wp, 2.0_wp, 1.0_wp]
 
  expected = 3.0_wp/5.0_wp
  result_ = mean(x)
  
  call assert_equals(expected, result_, tol)
end subroutine test_mean_simple

subroutine test_mean_w_simple
  use statistics, only: mean
  integer, parameter  :: n=5
  real(wp), parameter :: tol=1.0e-06_wp
  complex(wp)         :: w(n), x(n)
  complex(wp)         :: expected, result_
  
  w = [1.0_wp, 0.0_wp, -1.0_wp, 2.0_wp, 1.0_wp]
  x = [(2.0_wp, 0.0_wp), (0.0_wp, 1.0_wp), (0.0_wp, 0.0_wp), (1.0_wp, 1.0_wp), (-1.0_wp, 1.0_wp)]
  
  expected = (1.0_wp, 1.0_wp)
  result_ = mean(w, x)
  
  call assert_equals(expected, result_, tol)
end subroutine test_mean_w_simple

subroutine test_mean_ax
  use statistics, only: mean
  integer, parameter    :: n=10, m=4
  real(wp), parameter   :: tol=1.0e-06_wp
  real(wp), allocatable :: a(:,:), b(:,:)
  real(wp), allocatable :: expected(:), result_(:)
  type(Lcg48BasicRng)   :: brng

  if (.not. comm_world%initialized) call comm_world%init()
  brng = Lcg48BasicRng()

  allocate(a(n, m), b(m,n))
  call brng%rand(a)
  b = transpose(a)
 
  expected = mean(a, ax=1)
  result_ = mean(b, ax=2)
  
  call assert_equals(expected, result_, m, tol)
end subroutine test_mean_ax

subroutine test_cov_matrix
  use statistics, only: cov
  integer, parameter    :: n=3
  real(wp), parameter   :: tol=1.0e-04_wp
  real(wp)              :: w(n)
  real(wp)              :: x(n,n)
  real(wp)              :: expected(n,n)
  real(wp), allocatable :: result_(:,:)
  
  w = [1.0_wp, 1.5_wp, 0.5_wp]
  x(:,1) = [1.0_wp, -1.0_wp, 1.0_wp]
  x(:,2) = [2.1_wp,  2.0_wp, 1.9_wp]
  x(:,3) = [-0.5_wp, 1.0_wp, 2.0_wp]
  
  expected(:,1) = [ 1.636363_wp,  0.027272_wp, -0.545454_wp]
  expected(:,2) = [ 0.027272_wp,  0.007727_wp, -0.100000_wp]
  expected(:,3) = [-0.545454_wp, -0.100000_wp,  1.318181_wp]

  result_ = cov(w, x, 1)

  call assert_equals(expected, result_, n, n, tol)
end subroutine test_cov_matrix

subroutine test_cov_matrix_rc
  use statistics, only: cov
  integer, parameter       :: n=3
  real(wp), parameter      :: tol=1.0e-04_wp
  real(wp)                 :: w(n)
  complex(wp)              :: x(n,n)
  complex(wp)              :: expected(n,n)
  complex(wp), allocatable :: result_(:,:)
  
  w = [1.0_wp, 1.5_wp, 0.5_wp]
  x(:,1) = [(1.0_wp,1.0_wp), (-1.0_wp,-1.0_wp), (1.0_wp,0.0_wp)]
  x(:,2) = [(2.1_wp,0.2_wp), (2.0_wp,-0.1_wp), (1.9_wp,0.1_wp)]
  x(:,3) = [(-0.5_wp,0.0_wp), (1.0_wp,1.0_wp), (2.0_wp,-1.0_wp)]
  
  expected(:,1) = [(2.954545_wp, 0.000000_wp), (0.227272_wp, -0.1590909_wp), (-1.272727_wp, 0.181818_wp)]
  expected(:,2) = [(0.227272_wp, 0.1590909_wp), (0.038636_wp, 0.000000_wp), (-0.227272_wp, -0.136363_wp)]
  expected(:,3) = [(-1.272727_wp, -0.181818_wp), (-0.227272_wp, 0.136363_wp), (2.227272_wp, 0.000000_wp)]

  allocate(result_(n,n))
  result_ = cov(w, x, 1)

  call assert_equals(expected, result_, n, n, tol)
end subroutine test_cov_matrix_rc

subroutine test_var_cov
  use statistics, only: std, cov
  use mpi
  integer, parameter    :: n = 100
  real(wp), parameter   :: tol = 1.0e-06_wp
  real(wp)              :: w(n)
  complex(wp)           :: x(n), y(n), z(n)
  complex(wp)           :: expected, result_
  type(Lcg48BasicRng)   :: brng
  type(DistRng)         :: drng

  if (.not. comm_world%initialized) call comm_world%init()
  brng = Lcg48BasicRng()
  drng = DistRng(brng)

  call drng%normal(w, 1.0_wp, 0.05_wp)
  call drng%normal(x)
  call drng%normal(y)
  z = x + y
 
  expected = std(w, z)**2
  result_ = std(w,x)**2 + std(w,y)**2 + cov(w,x,y) + cov(w,y,x)

  call assert_equals(expected, result_, tol)
end subroutine test_var_cov

subroutine test_histogram_simple
  use statistics, only: histogram
  integer, parameter    :: n=6, nbins=11
  real(wp), parameter   :: a=1.0_wp/6.0_wp, b=2.0_wp/6.0_wp, tol=1.0e-06_wp
  real(wp), allocatable :: xhist(:), xhist_(:)
  real(wp), allocatable :: hist(:), hist_(:)
  real(wp)              :: x(n), xmin, xmax
  
  allocate(xhist(nbins), hist(nbins))
  xmin = 0.0_wp
  xmax = 11.0_wp
  x = [1.3_wp, 2.2_wp, 4.1_wp, 1.3_wp, 8.5_wp, 3.2_wp]
  xhist = [0.5_wp, 1.5_wp, 2.5_wp, 3.5_wp, 4.5_wp, 5.5_wp, 6.5_wp, 7.5_wp, 8.5_wp, 9.5_wp, 10.5_wp]
  hist = [0.0_wp, b, a, a, a, 0.0_wp, 0.0_wp, 0.0_wp, a, 0.0_wp, 0.0_wp]
  
  call histogram(x, hist_, xhist_, xmin, xmax, nbins)
  
  call assert_equals(xhist, xhist_, nbins, tol)
  call assert_equals(hist, hist_, nbins, tol)
end subroutine test_histogram_simple

subroutine test_histogram_auto
  use statistics, only: histogram
  integer, parameter    :: n=6, nbins=6
  real(wp), parameter   :: a=1.0_wp/6.0_wp, b=3.0_wp/6.0_wp, tol=1.0e-06_wp
  real(wp), allocatable :: hist(:), hist_(:)
  real(wp)              :: x(n)
  
  allocate(hist(nbins))
  x = [1.3_wp, 2.2_wp, 4.1_wp, 1.3_wp, 8.5_wp, 3.2_wp]
  hist = [b, a, a, 0.0_wp, 0.0_wp, a]
  
  call histogram(x, hist_, nbins=nbins)
  
  call assert_equals(hist, hist_, nbins, tol)
end subroutine test_histogram_auto

end module test_statistics
