! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module propagators
!******************************************************************************** 
!              
!       Module that impelements some basic propagators
!
!******************************************************************************** 

#include "preproc.inc"
use constants
use profiling
use lapack
 
implicit none

public

interface propagate_euler
  module procedure propagate_euler_r, propagate_euler_rc, propagate_euler_c
end interface propagate_euler

interface propagate_cn
  module procedure propagate_cn_r, propagate_cn_c
end interface propagate_cn

character(len=*), parameter :: gemm_prop = "gemm_prop"

contains

!******************************************************************************** 
!
!       Euler forward propagator
!               Input:
!                   h - Hamiltonian
!                   coeff - orbitals to propagate
!               Output:
!                   coeff - updated orbitals - inplace update
!
!                       E(h) coeff = e^{h} coeff  = (1 + h) U
!
!******************************************************************************** 
subroutine propagate_euler_r(coeff, h)
  real(wp), intent(inout) :: coeff(:,:)
  real(wp), intent(in)    :: h(:,:)
  !local variables
  real(wp), allocatable   :: coeff_(:,:)
  integer                 :: n, ne
  
  allocate(coeff_, source=coeff)
  n = size(coeff, 1)
  ne = size(coeff, 2)
  
  if (profile_code) call start_profiling(gemm_prop)
  call gemm("n", "n", n, ne, n, oner, h, n, coeff_, n, oner, coeff, n)
  if (profile_code) call end_profiling(gemm_prop, n, ne, n, "r")
end subroutine propagate_euler_r

subroutine propagate_euler_rc(coeff, h)
  complex(wp), intent(inout) :: coeff(:,:)
  real(wp), intent(in)       :: h(:,:)
  !local variables
  complex(wp), allocatable   :: coeff_(:,:)
  integer                    :: n, ne
  
  allocate(coeff_, source=coeff)
  n = size(coeff, 1)
  ne = size(coeff, 2)
  
  if (profile_code) call start_profiling(gemm_prop)
  call gemm("n", "n", n, ne, n, onec, h, n, coeff_, n, onec, coeff, n)
  if (profile_code) call end_profiling(gemm_prop, n, ne, n, "rc")
end subroutine propagate_euler_rc

subroutine propagate_euler_c(coeff, h)
  complex(wp), intent(inout) :: coeff(:,:)
  complex(wp), intent(in)    :: h(:,:)
  !local variables
  complex(wp), allocatable   :: coeff_(:,:)
  integer                    :: n, ne
  
  allocate(coeff_, source=coeff)
  n = size(coeff, 1)
  ne = size(coeff, 2)
  
  if (profile_code) call start_profiling(gemm_prop)
  call gemm("n", "n", n, ne, n, onec, h, n, coeff_, n, onec, coeff, n)
  if (profile_code) call end_profiling(gemm_prop, n, ne, n, "c")
end subroutine propagate_euler_c


!******************************************************************************** 
!
!       Crnak-Nicolson propagator
!               Input:
!                   h - Hamiltonian matrix
!                   coeff - orbitals to propagate
!               Output:
!                   coeff - updated orbitals - inplace update
!
!       CN(h) coeff = (1 - h/2) / (1 + h/2) coeff
!
!******************************************************************************** 
subroutine propagate_cn_r(coeff, h)
  use standalone, only: eyef_r
  real(wp), intent(inout) :: coeff(:,:)
  real(wp), intent(in)    :: h(:,:)
  !local variables
  real(wp), allocatable   :: h_(:,:), coeff_(:,:)
  integer                 :: n, ne
  
  allocate(coeff_, source=coeff)
  allocate(h_, source=h)
  n = size(coeff, 1)
  ne = size(coeff, 2)
  
  if (profile_code) call start_profiling(gemm_prop)
  call gemm("n", "n", n, ne, n, oner/2.0_wp, h, n, coeff_, n, oner, coeff, n)
  if (profile_code)call end_profiling(gemm_prop, n, ne, n, "r")

  h_ = eyef_r(n) - h_ / 2.0_wp
  call solve_ls(h_, coeff)
end subroutine propagate_cn_r

subroutine propagate_cn_c(coeff, h)
  use standalone, only: eyef_c
  complex(wp), intent(inout) :: coeff(:,:)
  complex(wp), intent(in)    :: h(:,:)
  !local variables
  complex(wp), allocatable   :: h_(:,:), coeff_(:,:)
  integer                    :: n, ne
  
  allocate(coeff_, source=coeff)
  allocate(h_, source=h)
  n = size(coeff, 1)
  ne = size(coeff, 2)
  
  if (profile_code) call start_profiling(gemm_prop)
  call gemm("n", "n", n, ne, n, onec/2.0_wp, h, n, coeff_, n, onec, coeff, n)
  if (profile_code) call end_profiling(gemm_prop, n, ne, n, "c")

  h_ = eyef_c(n) - h_ / 2.0_wp
  call solve_ls(h_, coeff)
end subroutine propagate_cn_c


!******************************************************************************** 
!
!       s4 propagator:
!
!               s_4(x) = s_2(fx) s_2((1-2f)x) s_2(fx),
!       where
!               s_2(x) = e^{h1x/2} e^{h2x} e^{h1x/2}         with
!
!               f = 1/(2-sqrt{3}{2})
!
!******************************************************************************** 
end module propagators
