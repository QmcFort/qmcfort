! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module test_expm_cheb_mod

use expm_cheb_mod
use expm_exact_mod
use mpi

use constants, only: wp
use fruit, only: assert_equals, run_test_case, test_set_summary
use lcg48_basic_rng, only: Lcg48BasicRng

implicit none

private
public ::  test_expm_cheb_mod_set

contains

subroutine test_expm_cheb_mod_set
  call run_test_case(test_expm_cheb_tol, "test two versions of chebyshev expansion with tolerance argument")
  call run_test_case(test_expm_cheb_r_tol, "test Chebyshev expansion against exact diagonalization")

  call test_set_summary("src/expm/expm_cheb_mod.f90")
end subroutine test_expm_cheb_mod_set

subroutine test_expm_cheb_tol
  integer, parameter    :: n=10, m=5
  real(wp), parameter   :: tol=1.0e-04_wp, tol2=1.0e-06_wp
  real(wp), allocatable :: A(:,:), C(:,:), C_(:,:)
  type(Lcg48BasicRng)   :: brng

  if (.not. comm_world%initialized) call comm_world%init()
  brng = Lcg48BasicRng()

  allocate(A(n,n))
  call brng%rand(A)
  A = 0.005_wp * 0.5_wp * (A + transpose(A))

  allocate(C(n,m), C_(n,m))
  call brng%rand(C)
  C = 0.05_wp * C
  C_ = C

  call expm_cheb_act(1.0_wp, A, C, tol=tol2)
  call expm_cheb_clen_act(1.0_wp, A, C_, kmax=25)

  call assert_equals(C, C_, n, m, tol)
end subroutine test_expm_cheb_tol

subroutine test_expm_cheb_r_tol
  integer, parameter    :: n=10, m=5
  real(wp), parameter   :: tol=1.0e-04_wp, tol2=1.0e-06_wp
  real(wp), allocatable :: A(:,:), C(:,:), C_(:,:)
  type(Lcg48BasicRng)   :: brng

  if (.not. comm_world%initialized) call comm_world%init()
  brng = Lcg48BasicRng()

  allocate(A(n,n))
  call brng%rand(A)
  A = 0.1_wp * (A + transpose(A))

  allocate(C(n,m), C_(n,m))
  call brng%rand(C)
  C_ = C

  call expm_cheb_act(1.0_wp, A, C, kmax=45)
  call expm_exact_act(1.0_wp, A, C_)

  call assert_equals(C, C_, n, m, tol)
end subroutine test_expm_cheb_r_tol

end module test_expm_cheb_mod