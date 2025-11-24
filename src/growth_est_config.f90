! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module growth_est_config

use constants
use qmcfort_in, only: add_input

implicit none

private
public GrowthEstConfig

type GrowthEstConfig
  logical  :: active = .false.
  integer  :: period = 1
  real(wp) :: damp = 0.1_wp
contains
  procedure :: read_tags => read_growth_est_config_tags
end type GrowthEstConfig

interface GrowthEstConfig
  module procedure finit_growth_est_config
end interface GrowthEstConfig

contains

!******************************************************************************** 
!
! Initialization of the GrowthEstConfig object
!
!******************************************************************************** 
subroutine init_growth_est_config(self, active, period, damp)
  type(GrowthEstConfig), intent(out) :: self
  logical, optional, intent(in)      :: active
  integer, optional, intent(in)      :: period
  real(wp), optional, intent(in)     :: damp

  if (present(active)) self%active = active
  if (present(period)) self%period = period
  if (present(damp)) self%damp = damp
end subroutine init_growth_est_config


!******************************************************************************** 
!
! GrowthEstConfig object constructor
!
!******************************************************************************** 
function finit_growth_est_config(active, period, damp) result(self)
  logical, optional, intent(in)      :: active
  integer, optional, intent(in)      :: period
  real(wp), optional, intent(in)     :: damp
  type(GrowthEstConfig)              :: self

  call init_growth_est_config(self, active, period, damp)
end function finit_growth_est_config


!******************************************************************************** 
!
! GrowthEstConfig reader
!
!******************************************************************************** 
subroutine read_growth_est_config_tags(self)
  class(GrowthEstConfig), intent(inout) :: self

  call add_input("growth_est_active", self%active)
  call add_input("growth_est_period", self%period)
  call add_input("growth_est_damp", self%damp)
end subroutine read_growth_est_config_tags

end module growth_est_config