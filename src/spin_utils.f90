! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module spin_utils

#include "preproc.inc"
use constants

implicit none

public

contains

!******************************************************************************** 
!
! Get spin channel for wave functions
!
!******************************************************************************** 
function get_spin_channel(spin, ispin) result(sp)
  integer, intent(in) :: spin
  integer, intent(in) :: ispin
  integer             :: sp

  if (ispin == 1) then
    sp = 1
  else
    sp = spin
  end if
end function get_spin_channel


!******************************************************************************** 
!
! Get spin channel for one-body operators
!
!******************************************************************************** 
function get_spin1_channel(spin, ispin1) result(sp)
  integer, intent(in) :: spin
  integer, intent(in) :: ispin1
  integer             :: sp

  if (ispin1 == 1) then
    sp = 1
  else
    sp = spin
  end if
end function get_spin1_channel


!******************************************************************************** 
!
! Get spin channel for two-body operators
!
!******************************************************************************** 
function get_spin2_channel(spin1, spin2, ispin2) result(sp)
  integer, intent(in) :: spin1
  integer, intent(in) :: spin2
  integer, intent(in) :: ispin2
  integer             :: sp

  if (ispin2 == 1) then
    sp = 1
  else if (spin1 == spin2) then
    sp = spin1
  else
    sp = 3
  end if
end function get_spin2_channel

end module spin_utils