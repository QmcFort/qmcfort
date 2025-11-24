! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module expm_mod

#include "../preproc.inc"

use constants
use expm_cheb_mod
use expm_exact_mod
use expm_krylov_mod
use expm_taylor_mod
use expm_utils

use log_manager_mod, only: logging

implicit none 

public

interface expm
  module procedure expm_r, expm_c, expm_cr
end interface expm

interface expm_act
  module procedure expm_act_r, expm_act_c, expm_act_cr
end interface expm_act

contains

!******************************************************************************** 
!
! Driver for the matrix exponentiation 
!
! INPUT:
!     f - prefactor in exponential
!     A - matrix to exponentiate
!     kmax - exponentiation order
!     tol - stopping criterion
!
! OUTPUT:
!     B - matrix exponential B = exp(f*A)
!
!******************************************************************************** 
subroutine expm_r(f, A, B, kmax, tol, expm_mode)
  real(wp), intent(in)           :: f
  real(wp), intent(in)           :: A(:,:)
  real(wp), intent(out)          :: B(:,:)
  integer, optional, intent(in)  :: kmax
  real(wp), optional, intent(in) :: tol
  integer, optional, intent(in)  :: expm_mode
  !local
  integer                        :: expm_mode_
  
  if (present(expm_mode)) then
    expm_mode_ = expm_mode
  else 
    expm_mode_ = default_expm_mode
  end if 

  select case(expm_mode_)
    case (EXPM_EXACT_MODE)
      call expm_exact(f, A, B)
    case (EXPM_TAYLOR_MODE)
      call expm_taylor(f, A, B, kmax=kmax, tol=tol)
    case (EXPM_CHEB_MODE)
      call expm_cheb(f, A, B, kmax=kmax, tol=tol)
    case default
      call logging%error(.true., "invalid value of the named integer constant expm_mode", __FILE__, __LINE__)
  end select
end subroutine expm_r

subroutine expm_c(f, A, B, kmax, tol, expm_mode)
  complex(wp), intent(in)        :: f
  complex(wp), intent(in)        :: A(:,:)
  complex(wp), intent(out)       :: B(:,:)
  integer, optional, intent(in)  :: kmax
  real(wp), optional, intent(in) :: tol
  integer, optional, intent(in)  :: expm_mode
  !local
  integer                        :: expm_mode_

  if (present(expm_mode)) then
    expm_mode_ = expm_mode
  else 
    expm_mode_ = default_expm_mode
  end if 
  
  select case(expm_mode_)
    case (EXPM_EXACT_MODE)
      call expm_exact(f, A, B)
    case (EXPM_TAYLOR_MODE)
      call expm_taylor(f, A, B, kmax=kmax, tol=tol)
    case (EXPM_CHEB_MODE)
      call expm_cheb(f, A, B, kmax=kmax, tol=tol)
    case default
      call logging%error(.true., "invalid value of the named integer constant expm_mode", __FILE__, __LINE__)
  end select
end subroutine expm_c

subroutine expm_cr(f, A, B, kmax, tol, expm_mode)
  complex(wp), intent(in)        :: f
  real(wp), intent(in)           :: A(:,:)
  complex(wp), intent(out)       :: B(:,:)
  integer, optional, intent(in)  :: kmax
  real(wp), optional, intent(in) :: tol
  integer, optional, intent(in)  :: expm_mode
  !local
  integer                        :: expm_mode_

  if (present(expm_mode)) then
    expm_mode_ = expm_mode
  else 
    expm_mode_ = default_expm_mode
  end if 
  
  select case(expm_mode_)
    case (EXPM_EXACT_MODE)
      call expm_exact(f, A, B)
    case (EXPM_TAYLOR_MODE)
      call expm_taylor(f, A, B, kmax=kmax, tol=tol)
    case (EXPM_CHEB_MODE)
      call expm_cheb(f, A, B, kmax=kmax, tol=tol)
    case default
      call logging%error(.true., "invalid value of the named integer constant expm_mode", __FILE__, __LINE__)
  end select
end subroutine expm_cr


!******************************************************************************** 
!
! Driver for the acting of exponential on the matrix
!
! INPUT:
!     f - prefactor in exponential
!     A - matrix to exponentiate
!     C - matrix to act on
!     kmax - exponentiation order
!     tol - stopping criterion
!
! OUTPUT:
!     C - exp(fA) * C
!
!******************************************************************************** 
subroutine expm_act_r(f, A, C, kmax, tol, expm_mode)
  real(wp), intent(in)           :: f
  real(wp), intent(in)           :: A(:,:)
  real(wp), intent(inout)        :: C(:,:)
  integer, optional, intent(in)  :: kmax
  real(wp), optional, intent(in) :: tol
  integer, optional, intent(in)  :: expm_mode
  !local
  integer                        :: expm_mode_
  
  if (present(expm_mode)) then
    expm_mode_ = expm_mode
  else 
    expm_mode_ = default_expm_mode
  end if 

  select case(expm_mode_)
    case (EXPM_EXACT_MODE)
      call expm_exact_act(f, A, C)
    case (EXPM_TAYLOR_MODE)
      call expm_taylor_act(f, A, C, kmax=kmax, tol=tol)
    case (EXPM_CHEB_MODE)
      call expm_cheb_act(f, A, C, kmax=kmax, tol=tol)
    case (EXPM_KRYLOV_MODE)
      call expm_krylov(f, A, C, k=kmax)
    case (EXPM_BLOCK_KRYLOV_MODE)
      call expm_block_krylov(f, A, C, k=kmax)
    case default
      call logging%error(.true., "invalid value of the named integer constant expm_mode", __FILE__, __LINE__)
  end select
end subroutine expm_act_r

subroutine expm_act_c(f, A, C, kmax, tol, expm_mode)
  complex(wp), intent(in)        :: f
  complex(wp), intent(in)        :: A(:,:)
  complex(wp), intent(inout)     :: C(:,:)
  integer, optional, intent(in)  :: kmax
  real(wp), optional, intent(in) :: tol
  integer, optional, intent(in)  :: expm_mode
  !local
  integer                        :: expm_mode_
  
  if (present(expm_mode)) then
    expm_mode_ = expm_mode
  else 
    expm_mode_ = default_expm_mode
  end if 

  select case(expm_mode_)
    case (EXPM_EXACT_MODE)
      call expm_exact_act(f, A, C)
    case (EXPM_TAYLOR_MODE)
      call expm_taylor_act(f, A, C, kmax=kmax, tol=tol)
    case (EXPM_CHEB_MODE)
      call expm_cheb_act(f, A, C, kmax=kmax, tol=tol)
    case (EXPM_KRYLOV_MODE)
      call expm_krylov(f, A, C, k=kmax)
    case (EXPM_BLOCK_KRYLOV_MODE)
      call expm_block_krylov(f, A, C, k=kmax)
    case default
      call logging%error(.true., "invalid value of the named integer constant expm_mode", __FILE__, __LINE__)
  end select
end subroutine expm_act_c

subroutine expm_act_cr(f, A, C, kmax, tol, expm_mode)
  complex(wp), intent(in)        :: f
  real(wp), intent(in)           :: A(:,:)
  complex(wp), intent(inout)     :: C(:,:)
  integer, optional, intent(in)  :: kmax
  real(wp), optional, intent(in) :: tol
  integer, optional, intent(in)  :: expm_mode
  !local
  integer                        :: expm_mode_
  
  if (present(expm_mode)) then
    expm_mode_ = expm_mode
  else 
    expm_mode_ = default_expm_mode
  end if 

  select case(expm_mode_)
    case (EXPM_EXACT_MODE)
      call expm_exact_act(f, A, C)
    case (EXPM_TAYLOR_MODE)
      call expm_taylor_act(f, A, C, kmax=kmax, tol=tol)
    case (EXPM_CHEB_MODE)
      call expm_cheb_act(f, A, C, kmax=kmax, tol=tol)
    case (EXPM_KRYLOV_MODE)
      call expm_krylov(f, A, C, k=kmax)
    case (EXPM_BLOCK_KRYLOV_MODE)
      call expm_block_krylov(f, A, C, k=kmax)
    case default
      call logging%error(.true., "invalid value of the named integer constant expm_mode", __FILE__, __LINE__)
  end select
end subroutine expm_act_cr

end module expm_mod