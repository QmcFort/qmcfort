! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module expm_exact_mod

#include "../preproc.inc"

use constants
use expm_utils
use lapack
use profiling
use standalone

implicit none

private
public :: expm_exact, expm_exact_act

interface expm_exact
  module procedure expm_exact_r, expm_exact_c, expm_exact_cr
end interface expm_exact

interface expm_exact_act
  module procedure expm_exact_act_r, expm_exact_act_c, expm_exact_act_cr
end interface expm_exact_act

contains

!******************************************************************************** 
!
! Matrix exponential exp(f*A) = U exp(f*ew) * U^{-1}
!
! INPUT:
!         f - prefactor in exponential
!         A - Input matrix A(NxN)
!
! OUTPUT:
!         B - Output matrix B = exp(f*A)
!
!******************************************************************************** 
subroutine expm_exact_r(f, A, B)
  real(wp), intent(in)     :: f
  real(wp), intent(in)     :: A(:,:)
  real(wp), intent(out)    :: B(:,:)
  !local
  integer                  :: i, n
  real(wp)                 :: residue
  real(wp), parameter      :: tol = 1.0e-06_wp
  complex(wp), allocatable :: U(:,:), U_(:,:), B_(:,:), ew(:)

  n = size(A, 1)
  allocate(U(n,n), U_(n,n), B_(n,n), ew(n))

  B = zeror
  B_ = zeroc

  call geev(A, ew, U)

  !calculate U^{-1}
  U_ = U
  call inverse(U_)

  !calculate exp(f*ew)
  ew = exp(f*ew)
  
  !calculate U*exp(f*ew)
  do i = 1, n
    U(:,i) = U(:,i) * ew(i)
  end do

  !calculate U*exp(f*ew)*U^{-1}
  call gemm("n", "n", n, n, n, onec, U, n, U_, n, zeroc, B_, n)

  B = real(B_, kind=wp)

  residue = sum(abs(B_ - B))
  if (residue > tol) write(*,*) "Non-vanishing imaginary part of the real matrix exponential"
end subroutine expm_exact_r

subroutine expm_exact_c(f, A, B)
  complex(wp), intent(in)  :: f
  complex(wp), intent(in)  :: A(:,:)
  complex(wp), intent(out) :: B(:,:)
  !local
  integer                  :: i, n 
  complex(wp), allocatable :: U(:,:), U_(:,:), ew(:)

  n = size(A, 1)
  allocate(U(n,n), U_(n,n), ew(n))

  B = zeroc

  call geev(A, ew, U)

  !calculate U^{-1}
  U_ = U
  call inverse(U_)

  !calculate exp(f*ew)
  ew = exp(f*ew)

  !calculate U*exp(f*ew)
  do i = 1, n
    U(:,i) = U(:,i) * ew(i)
  end do

  !calculate U*exp(f*ew)*U^{-1}
  call gemm("n", "n", n, n, n, onec, U, n, U_, n, zeroc, B, n)
end subroutine expm_exact_c

subroutine expm_exact_cr(f, A, B)
  complex(wp), intent(in)  :: f
  real(wp), intent(in)     :: A(:,:)
  complex(wp), intent(out) :: B(:,:)
  !local
  integer                  :: i, n 
  complex(wp), allocatable :: U(:,:), U_(:,:), ew(:)

  n = size(A, 1)
  allocate(U(n,n), U_(n,n), ew(n))

  B = zeroc

  call geev(A, ew, U)

  !calculate U^{-1}
  U_ = U
  call inverse(U_)

  !calculate exp(f*ew)
  ew = exp(f*ew)

  !calculate U*exp(f*ew)
  do i = 1, n
    U(:,i) = U(:,i) * ew(i)
  end do

  !calculate U*exp(f*ew)*U^{-1}
  call gemm("n", "n", n, n, n, onec, U, n, U_, n, zeroc, B, n)
end subroutine expm_exact_cr


!******************************************************************************** 
!
! Action of the matrix expoential on the matrix 
!
! exp(f*A) = U exp(f*ew) * U^{-1}
!
! INPUT:
!         f - prefactor in exponential
!         A - Input matrix A(NxN)
!
! OUTPUT:
!         B - Output matrix B = exp(f*A)
!
!******************************************************************************** 
subroutine expm_exact_act_r(f, A, C)
  real(wp), intent(in)     :: f
  real(wp), intent(in)     :: A(:,:)
  real(wp), intent(inout)  :: C(:,:)
  !local
  real(wp), parameter      ::  tol = 1.0e-06_wp
  real(wp)                 :: residue
  integer                  :: i, n, m 
  complex(wp), allocatable :: U(:,:), U_(:,:), C_(:,:), C__(:,:), ew(:)

  n = size(C, 1)
  m = size(C, 2)
  allocate(U(n,n), U_(n,n), C_(n,m), C__(n,m), ew(n))

  call geev(A, ew, U)

  !calculate U^{-1}
  U_ = U
  call inverse(U_)

  !calculate exp(f*ew)
  ew = exp(f*ew)

  !calculate U^{-1}*phi
  call gemm("n", "n", n, m, n, onec, U_, n, C, n, zeroc, C_, n)

  !calculate exp(f*ew)*U^{-1}*phi
  do i = 1, n
    C_(i,:) = ew(i) * C_(i,:)
  end do

  !calculate U*exp(f*ew)*U^{-1}*phi
  call gemm("n", "n", n, m, n, onec, U, n, C_, n, zeroc, C__, n)

  C = real(C__, kind=wp)

  residue = sum(abs(C__ - C))
  if (residue > tol) write(*,*) "Non-vanishing imaginary part of the real matrix exponential"
end subroutine expm_exact_act_r

subroutine expm_exact_act_c(f, A, C)
  complex(wp), intent(in)    :: f
  complex(wp), intent(in)    :: A(:,:)
  complex(wp), intent(inout) :: C(:,:)
  !local
  integer                    :: i, n, m 
  complex(wp), allocatable   :: U(:,:), U_(:,:), C_(:,:), ew(:)

  n = size(C, 1)
  m = size(C, 2)
  allocate(U(n,n), U_(n,n), C_(n, m), ew(n))
  
  call geev(A, ew, U)

  !calculate U^{-1}
  U_ = U
  call inverse(U_)

  !calculate exp(f*ew)
  ew = exp(f * ew)

  !calculate U^{-1}*phi
  call gemm("n", "n", n, m, n, onec, U_, n, C, n, zeroc, C_, n)

  !calculate exp(f*ew)*U^{-1}*phi
  do i = 1, n
    C_(i,:) = ew(i) * C_(i,:)
  end do

  !calculate U*exp(f*ew)*U^{-1}*phi
  call gemm("n", "n", n, m, n, onec, U, n, C_, n, zeroc, C, n)
end subroutine expm_exact_act_c

subroutine expm_exact_act_cr(f, A, C)
  complex(wp), intent(in)    :: f
  real(wp), intent(in)       :: A(:,:)
  complex(wp), intent(inout) :: C(:,:)
  !local
  integer                    :: i, n, m 
  complex(wp), allocatable   :: U(:,:), U_(:,:), C_(:,:), ew(:)

  n = size(C, 1)
  m = size(C, 2)
  allocate(U(n,n), U_(n,n), C_(n, m), ew(n))
  
  call geev(A, ew, U)

  !calculate U^{-1}
  U_ = U
  call inverse(U_)

  !calculate exp(f*ew)
  ew = exp(f * ew)

  !calculate U^{-1}*phi
  call gemm("n", "n", n, m, n, onec, U_, n, C, n, zeroc, C_, n)

  !calculate exp(f*ew)*U^{-1}*phi
  do i = 1, n
    C_(i,:) = ew(i) * C_(i,:)
  end do

  !calculate U*exp(f*ew)*U^{-1}*phi
  call gemm("n", "n", n, m, n, onec, U, n, C_, n, zeroc, C, n)
end subroutine expm_exact_act_cr

end module expm_exact_mod