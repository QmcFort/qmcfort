! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module noci_intags

use constants
use mpi, only: comm_world
use qmcfort_in, only: add_input

implicit none

private
public :: NociTags

type NociTags
  !NOCISelector related variables
  integer  :: ndet_max
  real(wp) :: ovlp_thresh
  real(wp) :: energy_thresh

  !NociAfqmcSelection related variables
  integer  :: nepochs
  integer  :: steps_per_epoch
  integer  :: select_from_step
  real(wp) :: energy_thresh_min
  real(wp) :: energy_thresh_max
  real(wp) :: sigma_eloc
  logical  :: update_afqmc_trial
contains
  procedure :: reader => read_noci_tags
end type NociTags

interface NociTags
  module procedure finit_noci_tags
end interface NociTags

contains 

!******************************************************************************** 
!
! Initialization of the NociTags object
!
!******************************************************************************** 
subroutine init_noci_tags(self)
  type(NociTags), intent(out) :: self

  call self%reader()
end subroutine init_noci_tags


!******************************************************************************** 
!
! NociTags constructor 
!
!******************************************************************************** 
function finit_noci_tags() result(self)
  type(NociTags) :: self

  call init_noci_tags(self)
end function finit_noci_tags


!******************************************************************************** 
!
! Read NOCI realated flags
!
!******************************************************************************** 
subroutine read_noci_tags(self)
  class(NociTags), intent(inout) :: self
  !local variables
  character(len=:), allocatable  :: tag
  logical                        :: found

  tag = "ndet_max"
  self%ndet_max = 100
  call add_input(tag, self%ndet_max, found)

  tag = "ovlp_thresh"
  self%ovlp_thresh = 1.0e-04_wp
  call add_input(tag, self%ovlp_thresh, found)

  tag = "energy_thresh"
  self%energy_thresh = 1.0e-04_wp
  call add_input(tag, self%energy_thresh, found)


  !NociAfqmcSelection variables
  tag = "nepochs"
  self%nepochs = 1
  call add_input(tag, self%nepochs, found)

  tag = "steps_per_epoch"
  self%steps_per_epoch = 1
  call add_input(tag, self%steps_per_epoch, found)

  tag = "select_from_step"
  self%select_from_step = 0
  call add_input(tag, self%select_from_step, found)

  tag = "energy_thresh_max"
  self%energy_thresh_max = 1.0e-04_wp
  call add_input(tag, self%energy_thresh_max, found)

  tag = "energy_thresh_min"
  self%energy_thresh_min = 1.0e-05_wp
  call add_input(tag, self%energy_thresh_min, found)

  tag = "sigma_eloc"
  self%sigma_eloc = 2.0_wp
  call add_input(tag, self%sigma_eloc, found)

  tag = "update_afqmc_trial"
  self%update_afqmc_trial = .true.
  call add_input(tag, self%update_afqmc_trial, found)
end subroutine read_noci_tags

end module noci_intags