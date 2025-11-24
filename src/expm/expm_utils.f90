! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0


module expm_utils

use constants

use string, only: lowercase

implicit none

public

interface get_expm_mode
  module procedure get_expm_mode_int, get_expm_mode_char
end interface get_expm_mode

!Expm types
integer, parameter :: EXPM_EXACT_MODE = 1
integer, parameter :: EXPM_TAYLOR_MODE = 2 
integer, parameter :: EXPM_CHEB_MODE = 3
integer, parameter :: EXPM_KRYLOV_MODE = 4
integer, parameter :: EXPM_BLOCK_KRYLOV_MODE = 5

integer                :: default_expm_mode = EXPM_TAYLOR_MODE
integer                :: default_expm_order = 6 
integer                :: defualt_expm_block_order = 4
integer                :: max_expm_order = 50

character(len=*), parameter :: gemm_expm_timer = "gemm_expm"

contains

!******************************************************************************** 
!
! Return named integer constant from a character 
!
!******************************************************************************** 
function get_expm_mode_int(expm_mode_char) result(expm_mode_int)
  character(len=*), intent(in) :: expm_mode_char
  integer                      :: expm_mode_int

  select case (lowercase(expm_mode_char))
    case ("exact")
      expm_mode_int = EXPM_EXACT_MODE
    case ("taylor")
      expm_mode_int = EXPM_TAYLOR_MODE
    case ("cheb", "chebyshev") 
      expm_mode_int = EXPM_CHEB_MODE
    case ("krylov")
      expm_mode_int = EXPM_KRYLOV_MODE
    case ("block_krylov")
      expm_mode_int = EXPM_BLOCK_KRYLOV_MODE
    case default 
      expm_mode_int = EXPM_TAYLOR_MODE
  end select
end function get_expm_mode_int


!******************************************************************************** 
!
! Return character for a given named constant
!
!******************************************************************************** 
function get_expm_mode_char(expm_mode_int) result(expm_mode_char)
  integer, intent(in)           :: expm_mode_int
  character(len=:), allocatable :: expm_mode_char

  select case (expm_mode_int)
    case (EXPM_EXACT_MODE)
      expm_mode_char = "exact"
    case (EXPM_TAYLOR_MODE)
      expm_mode_char = "taylor"
    case (EXPM_CHEB_MODE) 
      expm_mode_char = "cheb"
    case (EXPM_KRYLOV_MODE)
      expm_mode_char = "krylov"
    case (EXPM_BLOCK_KRYLOV_MODE)
      expm_mode_char = "block_krylov"
    case default 
      expm_mode_char = "taylor"
  end select
end function get_expm_mode_char

end module expm_utils