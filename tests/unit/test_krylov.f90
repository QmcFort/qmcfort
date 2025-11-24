! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module test_krylov

#include "../../src/preproc.inc"

use krylov
use mpi

use constants, only: wp
use fruit, only: assert_equals, run_test_case, test_set_summary
use lcg48_basic_rng, only: Lcg48BasicRng
use standalone, only: eye

implicit none 

private
public :: test_krylov_set

contains

subroutine test_krylov_set
  call run_test_case(test_arnoldi_orth, "test whether arnoldi process gives independent basis")
  call run_test_case(test_arnoldi_proj, "test projection of the matrix on the subspace")
  call run_test_case(test_block_arnoldi_proj, "test projection matrix in block_arnoldi")
  call run_test_case(test_arnoldi_proj_c, "test projection of the matrix on the subspace")

  call test_set_summary("src/krylov.f90")
end subroutine test_krylov_set

subroutine test_arnoldi_orth
  integer, parameter       :: n=20, k=5, kk=6
  real(wp), parameter      :: tol=1.0e-04_wp
  real(wp), allocatable    :: Ar(:,:), Ai(:,:), vr(:), vi(:)
  complex(wp), allocatable :: A(:,:), B(:,:), H(:,:), v(:)
  complex(wp), allocatable :: id(:,:), id2(:,:)
  type(Lcg48BasicRng)      :: brng

  if (.not. comm_world%initialized) call comm_world%init()
  brng = Lcg48BasicRng()

  allocate(Ar(n,n), Ai(n,n), A(n,n), vr(n), vi(n), v(n))
  allocate(id(kk,kk), id2(kk,kk))

  call brng%rand(Ar)
  call brng%rand(Ai)
  call brng%rand(vr)
  call brng%rand(vi)
  A = 0.005_wp * cmplx(Ar, Ai, kind=wp)
  v = cmplx(vr, vi, kind=wp)

  call arnoldi(A, v, k, B, H)
  call gemm("c","n",kk,kk,n,onec,B,n,B,n,zeroc,id,kk)
  call eye(id2)

  call assert_equals(id, id2, kk, kk, tol)
end subroutine test_arnoldi_orth

subroutine test_arnoldi_proj_c
  integer, parameter       :: n=100, k=8, kk=9
  real(wp), parameter      :: tol=1.0e-04_wp
  real(wp), allocatable    :: Ar(:,:), Ai(:,:), vr(:), vi(:)
  complex(wp), allocatable :: A(:,:), v(:), B(:,:), H(:,:), H_(:,:)
  type(Lcg48BasicRng)      :: brng

  if (.not. comm_world%initialized) call comm_world%init()
  brng = Lcg48BasicRng()

  allocate(Ar(n,n), Ai(n,n), vr(n), vi(n))
  allocate(A(n,n), v(n), H_(kk,kk))

  call brng%rand(Ar)
  call brng%rand(Ai)
  call brng%rand(vr)
  call brng%rand(vi)
  A = 0.01_wp * cmplx(Ar, Ai, kind=wp)
  v = cmplx(vr, vi, kind=wp)

  call arnoldi(A, v, k, B, H)
  H_ = matmul(transpose(conjg(B)), matmul(A,B))
 
  call assert_equals(H, H_, kk, kk, tol)
end subroutine test_arnoldi_proj_c

subroutine test_arnoldi_proj
  integer, parameter    :: n=100, k=8, kk=9
  real(wp), parameter   :: tol=1.0e-04_wp
  real(wp), allocatable :: A(:,:), v(:), B(:,:), H(:,:), H_(:,:)
  type(Lcg48BasicRng)   :: brng

  if (.not. comm_world%initialized) call comm_world%init()
  brng = Lcg48BasicRng()

  allocate(A(n,n), v(n), H_(kk,kk))

  call brng%rand(A)
  call brng%rand(v)
  A = 0.01_wp * A

  call arnoldi(A, v, k, B, H)
  H_ = matmul(transpose(B), matmul(A,B))
 
  call assert_equals(H, H_, kk, kk, tol)
end subroutine test_arnoldi_proj

subroutine test_block_arnoldi_proj
  integer, parameter    :: n=50, m=5, k=2, kk=3
  integer               :: mk, i
  real(wp), parameter   :: tol=1.0e-04_wp
  real(wp), allocatable :: A(:,:), V(:,:), B(:,:), H(:,:), H_(:,:)
  type(Lcg48BasicRng)   :: brng

  if (.not. comm_world%initialized) call comm_world%init()
  brng = Lcg48BasicRng()

  mk = m * kk
  allocate(A(n,n), V(n,m), H_(mk,mk))

  call brng%rand(A)
  call brng%rand(V)
  A = 0.01_wp * A
  V = 0.005_wp * V
  do i = 1, m
    V(i,i) = 1.0_wp
  end do

  call block_arnoldi(A, V, k, B, H)
  H_ = matmul(transpose(B), matmul(A,B))
 
  call assert_equals(H, H_, mk, mk, tol)
end subroutine test_block_arnoldi_proj

end module test_krylov