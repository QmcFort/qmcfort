! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module expm_krylov_mod

#include "../preproc.inc"

use constants
use expm_utils
use krylov
use lapack

use expm_taylor_mod, only: expm_taylor_act

implicit none

private
public :: expm_krylov, expm_block_krylov

interface expm_krylov
  module procedure expm_krylov_r_1d, expm_krylov_c_1d, expm_krylov_cr_1d
  module procedure expm_krylov_r_2d, expm_krylov_c_2d, expm_krylov_cr_2d
end interface expm_krylov

interface expm_block_krylov
  module procedure expm_block_krylov_r, expm_block_krylov_c, expm_block_krylov_cr
end interface expm_block_krylov

real(wp), parameter :: etol = 1.0e-09_wp

contains

subroutine expm_krylov_r_1d(f, A, v, k)
  real(wp), intent(in)          :: f
  real(wp), intent(in)          :: A(:,:)
  real(wp), intent(inout)       :: v(:)
  integer, optional, intent(in) :: k
  !local
  integer                       :: k_, kk, i, n
  real(wp)                      :: norma
  real(wp), allocatable         :: B(:,:), H(:,:), e(:,:)

  if (present(k)) then
    k_ = k
  else 
    k_ = default_expm_order
  end if

  n = size(v)
  kk = min(n, k_+1)

  allocate(e(kk,1))
  e = zeror
  e(1,1) = oner

  call arnoldi(A, v, k_, B, H)
  call expm_taylor_act(f, H, e, tol=etol)
  norma = norm2(v)
  call gemv("n", n, kk, oner, B, n, e(:,1), 1, zeror, v, 1)
  v = norma * v
end subroutine expm_krylov_r_1d

subroutine expm_krylov_c_1d(f, A, v, k)
  complex(wp), intent(in)       :: f
  complex(wp), intent(in)       :: A(:,:)
  complex(wp), intent(inout)    :: v(:)
  integer, optional, intent(in) :: k
  !local
  integer                       :: k_, kk, i, n
  complex(wp)                   :: norma
  complex(wp), allocatable      :: B(:,:), H(:,:), e(:,:)

  if (present(k)) then
    k_ = k
  else 
    k_ = default_expm_order
  end if

  n = size(v)
  kk = min(n, k_+1)

  allocate(e(kk,1))
  e = zeroc
  e(1,1) = onec

  call arnoldi(A, v, k_, B, H)
  call expm_taylor_act(f, H, e, tol=etol)
  norma = sqrt(real(dot_product(v, v), kind=wp))
  call gemv("n", n, kk, onec, B, n, e(:,1), 1, zeroc, v, 1)
  v = norma * v
end subroutine expm_krylov_c_1d

subroutine expm_krylov_cr_1d(f, A, v, k)
  complex(wp), intent(in)       :: f
  real(wp), intent(in)          :: A(:,:)
  complex(wp), intent(inout)    :: v(:)
  integer, optional, intent(in) :: k
  !local
  integer                       :: k_, kk, i, n
  complex(wp)                   :: norma
  complex(wp), allocatable      :: B(:,:), H(:,:), e(:,:)

  if (present(k)) then
    k_ = k
  else 
    k_ = default_expm_order
  end if

  n = size(v)
  kk = min(n, k_+1)

  allocate(e(kk,1))
  e = zeroc
  e(1,1) = onec

  call arnoldi(A, v, k_, B, H)
  call expm_taylor_act(f, H, e, tol=etol)
  norma = sqrt(real(dot_product(v, v), kind=wp))
  call gemv("n", n, kk, onec, B, n, e(:,1), 1, zeroc, v, 1)
  v = norma * v
end subroutine expm_krylov_cr_1d

subroutine expm_krylov_r_2d(f, A, V, k)
  real(wp), intent(in)          :: f
  real(wp), intent(in)          :: A(:,:)
  real(wp), intent(inout)       :: V(:,:)
  integer, optional, intent(in) :: k
  !local
  integer                       :: k_, kk, p, n, m
  real(wp)                      :: norma
  real(wp), allocatable         :: B(:,:,:), H(:,:,:), e(:,:), e1(:,:)

  if (present(k)) then
    k_ = k
  else 
    k_ = default_expm_order
  end if

  n = size(V, 1)
  m = size(V, 2)
  kk = min(n, k_+1)

  allocate(e1(kk,1), e(kk,1))
  e1 = zeroc
  e1(1,1) = onec

  call arnoldi(A, V, k_, B, H)

  do p = 1, m
    e = e1
    call expm_taylor_act(f, H(p,:,:), e, tol=etol)
    norma = norm2(V(:,p))
    call gemv("n", n, kk, oner, B(:,p,:), n, e(:,1), 1, zeror, V(:,p), 1)
    V(:,p) = norma * V(:,p)
  end do
end subroutine expm_krylov_r_2d

subroutine expm_krylov_c_2d(f, A, V, k)
  complex(wp), intent(in)       :: f
  complex(wp), intent(in)       :: A(:,:)
  complex(wp), intent(inout)    :: V(:,:)
  integer, optional, intent(in) :: k
  !local
  integer                       :: k_, kk, p, n, m
  complex(wp)                   :: norma
  complex(wp), allocatable      :: B(:,:,:), H(:,:,:), e(:,:), e1(:,:)

  if (present(k)) then
    k_ = k
  else 
    k_ = default_expm_order
  end if

  n = size(V, 1)
  m = size(V, 2)
  kk = min(n, k_+1)

  allocate(e1(kk,1), e(kk,1))
  e1 = zeroc
  e1(1,1) = onec

  call arnoldi(A, V, k_, B, H)

  do p = 1, m
    e = e1
    call expm_taylor_act(f, H(p,:,:), e, tol=etol)
    norma = sqrt(real(dot_product(V(:,p), V(:,p)), kind=wp))
    call gemv("n", n, kk, onec, B(:,p,:), n, e(:,1), 1, zeroc, V(:,p), 1)
    V(:,p) = norma * V(:,p)
  end do
end subroutine expm_krylov_c_2d

subroutine expm_krylov_cr_2d(f, A, V, k)
  complex(wp), intent(in)       :: f
  real(wp), intent(in)          :: A(:,:)
  complex(wp), intent(inout)    :: V(:,:)
  integer, optional, intent(in) :: k
  !local
  integer                       :: k_, kk, p, n, m
  complex(wp)                   :: norma
  complex(wp), allocatable      :: B(:,:,:), H(:,:,:), e(:,:), e1(:,:)

  if (present(k)) then
    k_ = k
  else 
    k_ = default_expm_order
  end if

  n = size(V, 1)
  m = size(V, 2)
  kk = min(n, k_+1)

  allocate(e1(kk,1), e(kk,1))
  e1 = zeroc
  e1(1,1) = onec

  call arnoldi(A, V, k_, B, H)

  do p = 1, m
    e = e1
    call expm_taylor_act(f, H(p,:,:), e, tol=etol)
    norma = sqrt(real(dot_product(V(:,p), V(:,p)), kind=wp))
    call gemv("n", n, kk, onec, B(:,p,:), n, e(:,1), 1, zeroc, V(:,p), 1)
    V(:,p) = norma * V(:,p)
  end do
end subroutine expm_krylov_cr_2d

subroutine expm_block_krylov_r(f, A, V, k)
  real(wp), intent(in)          :: f
  real(wp), intent(in)          :: A(:,:)
  real(wp), intent(inout)       :: V(:,:)
  integer, optional, intent(in) :: k
  !local
  integer                       :: k_, kk, p, n, m, mk
  real(wp), allocatable         :: B(:,:), H(:,:), e(:,:), V_(:,:), R(:,:)

  if (present(k)) then
    k_ = k
  else 
    k_ = defualt_expm_block_order
  end if

  n = size(V, 1)
  m = size(V, 2)
  kk = min(n/m, k_+1)
  mk = m*kk

  allocate(V_(n,m), R(m,m), e(mk,m))

  V_ = V
  call getqr(V_, R)
  e  = zeror
  e(1:m,1:m) = R

  call block_arnoldi(A, V, k_, B, H)
  call expm_taylor_act(f, H, e, tol=etol)
  call gemm("n", "n", n, m, mk, oner, B, n, e, mk, zeror, V, n)
end subroutine expm_block_krylov_r

subroutine expm_block_krylov_c(f, A, V, k)
  complex(wp), intent(in)       :: f
  complex(wp), intent(in)       :: A(:,:)
  complex(wp), intent(inout)    :: V(:,:)
  integer, optional, intent(in) :: k
  !local
  integer                       :: k_, kk, p, n, m, mk
  complex(wp), allocatable      :: B(:,:), H(:,:), e(:,:), V_(:,:), R(:,:)

  if (present(k)) then
    k_ = k
  else 
    k_ = defualt_expm_block_order
  end if

  n = size(V, 1)
  m = size(V, 2)
  kk = min(n/m, k_+1)
  mk = m*kk

  allocate(V_(n,m), R(m,m), e(mk,m))

  V_ = V
  call getqr(V_, R)
  e  = zeroc
  e(1:m,1:m) = R

  call block_arnoldi(A, V, k_, B, H)
  call expm_taylor_act(f, H, e, tol=etol)
  call gemm("n", "n", n, m, mk, onec, B, n, e, mk, zeroc, V, n)
end subroutine expm_block_krylov_c

subroutine expm_block_krylov_cr(f, A, V, k)
  complex(wp), intent(in)       :: f
  real(wp), intent(in)          :: A(:,:)
  complex(wp), intent(inout)    :: V(:,:)
  integer, optional, intent(in) :: k
  !local
  integer                       :: k_, kk, p, n, m, mk
  complex(wp), allocatable      :: B(:,:), H(:,:), e(:,:), V_(:,:), R(:,:)

  if (present(k)) then
    k_ = k
  else 
    k_ = defualt_expm_block_order
  end if

  n = size(V, 1)
  m = size(V, 2)
  kk = min(n/m, k_+1)
  mk = m*kk

  allocate(V_(n,m), R(m,m), e(mk,m))

  V_ = V
  call getqr(V_, R)
  e  = zeroc
  e(1:m,1:m) = R

  call block_arnoldi(A, V, k_, B, H)
  call expm_taylor_act(f, H, e, tol=etol)
  call gemm("n", "n", n, m, mk, onec, B, n, e, mk, zeroc, V, n)
end subroutine expm_block_krylov_cr

end module expm_krylov_mod