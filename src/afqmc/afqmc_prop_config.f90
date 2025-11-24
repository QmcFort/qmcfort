! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module afqmc_prop_config

use constants
use expm_utils

use qmcfort_in, only: add_input

implicit none

private
public :: AfqmcPropConfig, get_prop_mode, PROP_S0_MODE, PROP_S1_MODE, PROP_S2_MODE

!prop types
integer, parameter :: PROP_S0_MODE = 1
integer, parameter :: PROP_S1_MODE = 2
integer, parameter :: PROP_S2_MODE = 3

type AfqmcPropConfig
  integer :: prop_mode = PROP_S2_MODE
  integer :: expm_mode = EXPM_TAYLOR_MODE
  integer :: expm_order = 6
contains
  procedure :: read_tags => read_afqmc_prop_config_tags
end type AfqmcPropConfig

interface AfqmcPropConfig
  module procedure finit_afqmc_prop_config
end interface AfqmcProPConfig

interface get_prop_mode
  module procedure get_prop_mode_int, get_prop_mode_char
end interface get_prop_mode

contains

!******************************************************************************** 
!
! Initialization of the AfqmcPropConfig object
!
!******************************************************************************** 
subroutine init_afqmc_prop_config(self, prop_mode, expm_mode, expm_order)
  type(AfqmcPropConfig), intent(out) :: self
  integer, optional, intent(in)      :: prop_mode
  integer, optional, intent(in)      :: expm_mode
  integer, optional, intent(in)      :: expm_order

  if (present(prop_mode)) self%prop_mode = prop_mode
  if (present(expm_mode)) self%expm_mode = expm_mode
  if (present(expm_order)) self%expm_order = expm_order
end subroutine init_afqmc_prop_config


!******************************************************************************** 
!
! AfqmcPropConfig constructor
!
!******************************************************************************** 
function finit_afqmc_prop_config(prop_mode, expm_mode, expm_order) result(self)
  integer, optional, intent(in) :: prop_mode
  integer, optional, intent(in) :: expm_mode
  integer, optional, intent(in) :: expm_order
  type(AfqmcPropConfig)         :: self

  call init_afqmc_prop_config(self, prop_mode, expm_mode, expm_order)
end function finit_afqmc_prop_config


!******************************************************************************** 
!
! Read AfqmcProPConfig entries from input tags
!
!******************************************************************************** 
subroutine read_afqmc_prop_config_tags(self)
  class(AfqmcProPConfig), intent(inout) :: self
  !local
  logical                               :: is_found
  character(len=charlen)                :: prop_mode, expm_mode

  call add_input("prop", prop_mode, found=is_found)
  if (is_found) then
    self%prop_mode = get_prop_mode(prop_mode)
  end if

  call add_input("expm", expm_mode, found=is_found)
  if (is_found) then
    self%expm_mode = get_expm_mode(expm_mode)
  end if

  call add_input("expm_order", self%expm_order)
end subroutine read_afqmc_prop_config_tags


!******************************************************************************** 
!
! Return named integer constant from a character 
!
!******************************************************************************** 
function get_prop_mode_int(prop_mode_char) result(prop_mode_int)
  character(len=*), intent(in) :: prop_mode_char
  integer                      :: prop_mode_int

  select case (lowercase(prop_mode_char))
    case ("s0")
      prop_mode_int = PROP_S0_MODE
    case ("s1")
      prop_mode_int = PROP_S1_MODE
    case ("s2")
      prop_mode_int = PROP_S2_MODE
    case default 
      prop_mode_int = PROP_S2_MODE
  end select
end function get_prop_mode_int


!******************************************************************************** 
!
! Return character for a given named constant
!
!******************************************************************************** 
function get_prop_mode_char(prop_mode_int) result(prop_mode_char)
  integer, intent(in)           :: prop_mode_int
  character(len=:), allocatable :: prop_mode_char

  select case (prop_mode_int)
    case (PROP_S0_MODE)
      prop_mode_char = "s0"
    case (PROP_S1_MODE)
      prop_mode_char = "s1"
    case (PROP_S2_MODE)
      prop_mode_char = "s2"
    case default 
      prop_mode_char = "s2"
  end select
end function get_prop_mode_char

end module afqmc_prop_config