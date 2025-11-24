! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module test_expm_mod

use expm_mod
use mpi

use constants, only: wp
use fruit, only: assert_equals, run_test_case, test_set_summary
use lcg48_basic_rng, only: Lcg48BasicRng
use standalone, only: imag

implicit none

private
public :: test_expm_mod_set

contains

subroutine test_expm_mod_set
  call run_test_case(test_expm_cheb, "test matrix exponentiation via cheb. polynomials")
  call run_test_case(test_expm_cheb_act, "test action of matrix exponentiation via cheb. polynomials")

  call test_set_summary("src/expm/expm_mod.f90")
end subroutine test_expm_mod_set

subroutine test_expm_cheb
  integer, parameter    :: n = 10
  real(wp), parameter   :: tol = 0.001_wp
  real(wp), allocatable :: A(:,:), X(:,:), B(:,:), B_(:,:)
  type(Lcg48BasicRng)   :: brng

  if (.not. comm_world%initialized) call comm_world%init()
  brng = Lcg48BasicRng()

  allocate(X(n,n), A(n,n), B(n,n), B_(n,n))

  call brng%rand(X)
  A = 0.005_wp * (X + transpose(X))

  call expm(1.0_wp, A, B, expm_mode=EXPM_EXACT_MODE)
  call expm(1.0_wp, A, B_, kmax=10, expm_mode=EXPM_CHEB_MODE)

  call assert_equals(B, B_, n, n, tol)
end subroutine test_expm_cheb

subroutine test_expm_cheb_act
  integer, parameter       :: n=10, m=5
  real(wp), parameter      :: tol = 0.001_wp
  complex(wp), allocatable :: A(:,:), C_(:,:), C(:,:)
  real(wp), allocatable    :: X(:,:), Y(:,:), XX(:,:), YY(:,:)
  type(Lcg48BasicRng)      :: brng

  if (.not. comm_world%initialized) call comm_world%init()
  brng = Lcg48BasicRng()

  allocate(X(n,n), Y(n,n),  XX(n,m), YY(n,m))
  allocate(A(n,n), C(n,m))

  call brng%rand(X)
  call brng%rand(Y)
  X = 0.005_wp * (X + transpose(X))
  Y = 0.005_wp * (Y - transpose(Y))
  A = cmplx(X, Y, kind=wp)

  call brng%rand(XX)
  call brng%rand(YY)
  C = 0.05 * cmplx(XX, YY, wp)
  allocate(C_, source=C)

  call expm_act((0.1_wp, 0.1_wp), A, C, expm_mode=EXPM_EXACT_MODE)
  call expm_act((0.1_wp, 0.1_wp), A, C_, kmax=10, expm_mode=EXPM_CHEB_MODE)

  call assert_equals(C_, C, n, m, tol)
end subroutine test_expm_cheb_act

end module test_expm_mod