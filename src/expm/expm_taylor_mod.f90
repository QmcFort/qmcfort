! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module expm_taylor_mod

#include "../preproc.inc"

use constants
use expm_utils
use lapack
use profiling

use standalone, only: eyef_r, eyef_c

implicit none

private
public :: expm_taylor, expm_taylor_act

interface expm_taylor
  module procedure expm_taylor_r, expm_taylor_c, expm_taylor_cr
end interface expm_taylor

interface expm_taylor_act
   module procedure expm_taylor_act_r, expm_taylor_act_c, expm_taylor_act_cr
end interface expm_taylor_act

contains

!******************************************************************************** 
!
! Taylor expansion to approximate matrix exponential
!
! INPUT:  
!         f - prefactor in exponential
!         A - matrix to exponentiate
!         kmax - truncation order of the Taylor series
!         tol - stopping criterion
!
! OUTPUT:  
!         expm_A - exp(f*A)
!
!******************************************************************************** 
subroutine expm_taylor_r(f, A, expm_A, kmax, tol)
  real(wp), intent(in)           :: f
  real(wp), intent(in)           :: A(:,:)
  real(wp), intent(out)          :: expm_A(:,:)
  integer, optional, intent(in)  :: kmax
  real(wp), optional, intent(in) :: tol
  !local
  integer                        :: i, n, kmax_
  real(wp)                       :: res, fact
  real(wp), allocatable          :: residual(:,:), temp(:,:)

  if (present(kmax)) then
    kmax_ = kmax
  else
    kmax_ = default_expm_order
  end if

  if (present(tol)) kmax_ = max_expm_order

  n = size(A, 1) 
  allocate(residual(n,n), temp(n,n))

  fact = f
  expm_A = eyef_r(n) + f*A

  residual = A
  do i = 2, kmax_
    fact = fact * f / real(i, wp)
    temp = residual
    call gemm("n", "n", n, n, n, oner, temp, n, A, n, zeror, residual, n)
    expm_A = expm_A + fact*residual
    if (present(tol)) then
      res = maxval(abs(fact*residual))
      if (res < tol) return
    end if
  end do  

  if (present(tol)) write(*,"(1x,a,2es12.4)") "expm_taylor_r did not converge: res tol = ", res, tol
end subroutine expm_taylor_r

subroutine expm_taylor_c(f, A, expm_A, kmax, tol)
  complex(wp), intent(in)        :: f
  complex(wp), intent(in)        :: A(:,:)
  complex(wp), intent(out)       :: expm_A(:,:)
  integer, optional, intent(in)  :: kmax
  real(wp), optional, intent(in) :: tol 
  !local 
  integer                        :: i, n, kmax_
  real(wp)                       :: res
  complex(wp)                    :: fact 
  complex(wp), allocatable       :: residual(:,:), temp(:,:)

  if (present(kmax)) then
    kmax_ = kmax
  else
    kmax_ = default_expm_order
  end if

  if (present(tol)) kmax_ = max_expm_order

  n = size(A, 1) 
  allocate(residual(n,n), temp(n,n))

  fact = f
  expm_A = eyef_c(n) + f*A

  residual = A
  do i = 2, kmax_
    fact = fact * f / real(i, wp)
    temp = residual
    call gemm("n", "n", n, n, n, onec, temp, n, A, n, zeroc, residual, n)
    expm_A = expm_A + fact*residual
    if (present(tol)) then
      res = maxval(abs(fact*residual))
      if (res < tol) return
    end if
  end do  

  if (present(tol)) write(*,"(1x,a,2es12.4)") "expm_taylor_c did not converge: res tol = ", res, tol
end subroutine expm_taylor_c 

subroutine expm_taylor_cr(f, A, expm_A, kmax, tol)
  complex(wp), intent(in)        :: f
  real(wp), intent(in)           :: A(:,:)
  complex(wp), intent(out)       :: expm_A(:,:)
  integer, optional, intent(in)  :: kmax
  real(wp), optional, intent(in) :: tol
  !local
  integer                        :: i, n, kmax_
  real(wp)                       :: res
  complex(wp)                    :: fact 
  real(wp), allocatable          :: residual(:,:), temp(:,:)

  if (present(kmax)) then
    kmax_ = kmax
  else
    kmax_ = default_expm_order
  end if

  if (present(tol)) kmax_ = max_expm_order

  n = size(A, 1) 
  allocate(residual(n,n), temp(n,n))

  fact = f
  expm_A = eyef_c(n) + f*A

  residual = A
  do i = 2, kmax_
    fact = fact * f / real(i, wp)
    temp = residual
    call gemm("n", "n", n, n, n, oner, temp, n, A, n, zeror, residual, n)
    expm_A = expm_A + fact*residual
    if (present(tol)) then
      res = maxval(abs(fact*residual))
      if (res < tol) return
    end if
  end do  

  if (present(tol)) write(*,"(1x,a,2es12.4)") "expm_taylor_cr did not converge: res tol = ", res, tol
end subroutine expm_taylor_cr


!******************************************************************************** 
!
! Taylor expansion to approximate action of the matrix exponetial
!
! INPUT:
!         f - prefactor in exponential
!         A - matrix to exponentiate
!         C - matrix to act on
!         kmax - truncation order of the Taylor series
!         tol - stopping criterion
!
! OUTPUT:
!         C - exp(f*A) C 
!
!******************************************************************************** 
subroutine expm_taylor_act_r(f, A, C, kmax, tol)
  real(wp), intent(in)           :: f
  real(wp), intent(in)           :: A(:,:)
  real(wp), intent(inout)        :: C(:,:)
  integer, optional, intent(in)  :: kmax
  real(wp), optional, intent(in) :: tol
  !local
  integer                        :: i, n, m, kmax_
  real(wp)                       :: res, fact 
  real(wp), allocatable          :: residual(:,:), temp(:,:)

  if (present(kmax)) then
    kmax_ = kmax
  else
    kmax_ = default_expm_order
  end if

  if (present(tol)) kmax_ = max_expm_order

  n = size(C, 1)
  m = size(C, 2)
  allocate(residual(n,m), temp(n,m))

  fact = oner
  residual = C

  do i = 1, kmax_
    fact = fact * f / real(i, wp)
    temp = residual
    if (profile_code) call start_profiling(gemm_expm_timer)
    call gemm("n", "n", n, m, n, oner, A, n, temp, n, zeror, residual, n)
    if (profile_code) call end_profiling(gemm_expm_timer, n, m, n, "r")
    C = C + fact*residual
    if (present(tol)) then
      res = maxval(abs(fact*residual))
      if (res < tol) return
    end if
  end do  

  if (present(tol)) write(*,"(1x,a,2es12.4)") "expm_taylor_act_r did not converge: res tol = ", res, tol
end subroutine expm_taylor_act_r

subroutine expm_taylor_act_c(f, A, C, kmax, tol)
  complex(wp), intent(in)        :: f
  complex(wp), intent(in)        :: A(:,:)
  complex(wp), intent(inout)     :: C(:,:)
  integer, optional, intent(in)  :: kmax
  real(wp), optional, intent(in) :: tol
  !local
  integer                        :: i, n, m, kmax_
  real(wp)                       :: res
  complex(wp)                    :: fact 
  complex(wp), allocatable       :: residual(:,:), temp(:,:)

  if (present(kmax)) then
    kmax_ = kmax
  else
    kmax_ = default_expm_order
  end if

  if (present(tol)) kmax_ = max_expm_order

  n = size(C, 1)
  m = size(C, 2)
  allocate(residual(n,m), temp(n,m))

  fact = onec
  residual = C

  do i = 1, kmax_
    fact = fact * f / real(i, wp)
    temp = residual
    if (profile_code) call start_profiling(gemm_expm_timer)
    call gemm("n", "n", n, m, n, onec, A, n, temp, n, zeroc, residual, n)
    if (profile_code) call end_profiling(gemm_expm_timer, n, m, n, "c")
    C = C + fact*residual
    if (present(tol)) then
      res = maxval(abs(fact*residual))
      if (res < tol) return
    end if
  end do  

  if (present(tol)) write(*,"(1x,a,2es12.4)") "expm_taylor_act_c did not converge: res tol = ", res, tol
end subroutine expm_taylor_act_c

subroutine expm_taylor_act_cr(f, A, C, kmax, tol)
  complex(wp), intent(in)        :: f
  real(wp), intent(in)           :: A(:,:)
  complex(wp), intent(inout)     :: C(:,:)
  integer, optional, intent(in)  :: kmax
  real(wp), optional, intent(in) :: tol
  !local
  integer                        :: i, n, m, kmax_
  real(wp)                       :: res
  complex(wp)                    :: fact 
  complex(wp), allocatable       :: residual(:,:), temp(:,:)

  if (present(kmax)) then
    kmax_ = kmax
  else
    kmax_ = default_expm_order
  end if

  if (present(tol)) kmax_ = max_expm_order

  n = size(C, 1)
  m = size(C, 2)
  allocate(residual(n,m), temp(n,m))

  fact = onec
  residual = C

  do i = 1, kmax_
    fact = fact * f / real(i, wp)
    temp = residual
    if (profile_code) call start_profiling(gemm_expm_timer)
    call gemm("n", "n", n, m, n, onec, A, n, temp, n, zeroc, residual, n)
    if (profile_code) call end_profiling(gemm_expm_timer, n, m, n, "rc")
    C = C + fact*residual
    if (present(tol)) then
      res = maxval(abs(fact*residual))
      if (res < tol) return
    end if
  end do  

  if (present(tol)) write(*,"(1x,a,2es12.4)") "expm_taylor_act_cr did not converge: res tol = ", res, tol
end subroutine expm_taylor_act_cr

end module expm_taylor_mod