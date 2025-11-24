! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module noci_selector_des

use constants

implicit none

private 
public :: NOCISelectorDes

type NOCISelectorDes
  integer  :: ndet_max
  real(wp) :: ovlp_thresh
  real(wp) :: energy_thresh
end type NOCISelectorDes

interface NOCISelectorDes
  module procedure finit_noci_selector_des
end interface NOCISelectorDes

contains

!******************************************************************************** 
!
! Initialization of the NOCISelectorDes
!
!******************************************************************************** 
subroutine init_noci_selector_des(ndet_max, ovlp_thresh, energy_thresh, self)
  integer, intent(in)                :: ndet_max
  real(wp), intent(in)               :: ovlp_thresh
  real(wp), intent(in)               :: energy_thresh
  type(NOCISelectorDes), intent(out) :: self

  self%ndet_max = ndet_max
  self%ovlp_thresh = ovlp_thresh
  self%energy_thresh = energy_thresh
end subroutine init_noci_selector_des


!******************************************************************************** 
!
! NOCISelectorDes constructor
!
!******************************************************************************** 
function finit_noci_selector_des(ndet_max, ovlp_thresh, energy_thresh) result(self)
  integer, intent(in)   :: ndet_max
  real(wp), intent(in)  :: ovlp_thresh
  real(wp), intent(in)  :: energy_thresh
  type(NOCISelectorDes) :: self

  call init_noci_selector_des(ndet_max, ovlp_thresh, energy_thresh, self)
end function finit_noci_selector_des

end module noci_selector_des