! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module test_expm_krylov_mod

#include "../../../src/preproc.inc"

use expm_exact_mod
use expm_krylov_mod
use expm_taylor_mod
use mpi

use constants, only: wp
use fruit, only: assert_equals, run_test_case, test_set_summary
use lcg48_basic_rng, only: Lcg48BasicRng

implicit none

private
public ::  test_expm_krylov_mod_set

contains

subroutine test_expm_krylov_mod_set
  call run_test_case(test_expm_krylov, "compare krylov_exp with expm_taylor_act")
  call run_test_case(test_expm_krylov_rc, "compare krylov_exp with expm_taylor_act for rc arrays and sprec")
  call run_test_case(test_expm_block_krylov_rc, "compare block_krylov_exp with expm_taylor_act for rc and sprec")
  call run_test_case(test_expm_block_krylov_c, "compare block_krylov_exp with expm_taylor for complex matrices")

  call test_set_summary("src/expm/expm_krylov_mod.f90")
end subroutine test_expm_krylov_mod_set

subroutine test_expm_krylov
  integer, parameter         :: n=50, m=5, k=10, kk=11, prec=wp
  real(prec), parameter      :: tol=1.0e-04_prec
  real(prec), allocatable    :: Ar(:,:), Ai(:,:), vr(:,:), vi(:,:)
  complex(prec), allocatable :: A(:,:), B(:,:), H(:,:), v(:,:), v1(:,:), v2(:,:)
  type(Lcg48BasicRng)        :: brng

  if (.not. comm_world%initialized) call comm_world%init()
  brng = Lcg48BasicRng()

  allocate(Ar(n,n), Ai(n,n), A(n,n), vr(n,m), vi(n,m), v(n,m), v1(n,m), v2(n,m))
  call brng%rand(Ar)
  call brng%rand(Ai)
  call brng%rand(vr)
  call brng%rand(vi)
  A = 0.01_prec * cmplx(Ar, Ai, kind=prec)
  v = cmplx(vr, vi, kind=prec)
  v1 = v
  v2 = v

  call expm_taylor_act((1.0_prec, 0.0_prec), A, v, kmax=k)
  call expm_exact_act((1.0_prec, 0.0_prec), A, v2)
  call expm_krylov((1.0_prec, 0.0_prec), A, v1, k)

  call assert_equals(v-v2, v1-v2, n, m, tol)
end subroutine test_expm_krylov

subroutine test_expm_krylov_rc
  integer, parameter         :: n=50, k=8, kk=9, prec=wp
  real(prec), parameter      :: tol=1.0e-04_prec
  real(prec), allocatable    :: A(:,:), vr(:), vi(:)
  complex(prec), allocatable :: B(:,:), H(:,:), v(:), v1(:,:)
  type(Lcg48BasicRng)        :: brng

  if (.not. comm_world%initialized) call comm_world%init()
  brng = Lcg48BasicRng()

  allocate(A(n,n), vr(n), vi(n), v(n), v1(n,1))
  call brng%rand(A)
  call brng%rand(vr)
  call brng%rand(vi)
  A = 0.01_prec * A
  v = cmplx(vr, vi, kind=wp)
  v1(:,1) = v

  call expm_taylor_act(onec, A, v1, kmax=k)
  call expm_krylov(onec, A, v, k)

  call assert_equals(v, v1(:,1), n, tol)
end subroutine test_expm_krylov_rc

subroutine test_expm_block_krylov_rc
  integer, parameter         :: n=50, m=5, k=4, kt=20, prec=wp
  integer                    :: i
  real(prec), parameter      :: tol=1.0e-04_prec
  real(prec), allocatable    :: A(:,:), Vr(:,:), Vi(:,:)
  complex(prec), allocatable :: B(:,:), H(:,:), V(:,:), V_(:,:)
  type(Lcg48BasicRng)        :: brng

  if (.not. comm_world%initialized) call comm_world%init()
  brng = Lcg48BasicRng()

  allocate(A(n,n), Vr(n,m), Vi(n,m), V(n,m), V_(n,m))
  call brng%rand(A) 
  call brng%rand(Vr) 
  call brng%rand(Vi) 
  A = 0.1_prec * A
  V = cmplx(0.1_prec*Vr, 0.05_prec*Vi, kind=prec)
  do i = 1, m
    V(i,i) = 1.0_prec
  end do
  V_ = V

  call expm_taylor_act(onec, A, V, kmax=kt)
  call expm_block_krylov(onec, A, V_, k=5)

  call assert_equals(V, V_, n, m, tol)
end subroutine test_expm_block_krylov_rc

subroutine test_expm_block_krylov_c
  integer, parameter         :: n=100, m=10, prec=wp
  integer                    :: i
  real(prec), parameter      :: tol=1.0e-04_prec
  real(prec), allocatable    :: Vr(:,:), Vi(:,:), Ar(:,:), Ai(:,:)
  complex(prec), allocatable :: A(:,:), B(:,:), H(:,:), V(:,:), V_(:,:)
  type(Lcg48BasicRng)        :: brng

  if (.not. comm_world%initialized) call comm_world%init()
  brng = Lcg48BasicRng()

  allocate(A(n,n), Ar(n,n), Ai(n,n), Vr(n,m), Vi(n,m), V(n,m), V_(n,m))
  call brng%rand(Ar)
  call brng%rand(Ai) 
  call brng%rand(Vr) 
  call brng%rand(Vi) 
  A = cmplx(0.1_prec*Ar, 0.1_prec*Ai, kind=prec)
  V = cmplx(0.1_prec*Vr, 0.1_prec*Vi, kind=prec)
  do i = 1, m
    V(i,i) = 1.0_prec
  end do
  V_ = V

  call expm_taylor_act(onec, A, V, kmax=30)
  call expm_block_krylov(onec, A, V_, k=6)

  call assert_equals(V, V_, n, m, tol)
end subroutine test_expm_block_krylov_c

end module test_expm_krylov_mod