! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module afqmc_proj_config

use constants

use qmcfort_in, only: add_input

implicit none 

private
public AfqmcProjConfig

type AfqmcProjConfig
  character(len=charlen) :: projection = "phaseless"
  logical                :: auto_hybrid = .false.
  real(wp)               :: hybrid = 1.0_wp
  real(wp)               :: mix = 1.0_wp
  real(wp)               :: min_cos_hp = -0.5_wp
  logical                :: accept_reject = .false.
contains
  procedure :: read_tags => read_afqmc_proj_config_tags
end type AfqmcProjConfig

interface AfqmcProjConfig
  module procedure finit_afqmc_proj_config
end interface AfqmcProjConfig

contains

!******************************************************************************** 
!
! Initialization of the AfqmcProjConfig object
!
!******************************************************************************** 
subroutine init_afqmc_proj_config(self, projection, auto_hybrid, hybrid, mix, min_cos_hp, accept_reject)
  type(AfqmcProjConfig), intent(out)     :: self
  character(len=*), optional, intent(in) :: projection
  logical, optional, intent(in)          :: auto_hybrid
  real(wp), optional, intent(in)         :: hybrid
  real(wp), optional, intent(in)         :: mix
  real(wp), optional, intent(in)         :: min_cos_hp
  logical, optional, intent(in)          :: accept_reject

  if (present(projection)) self%projection = projection
  if (present(auto_hybrid)) self%auto_hybrid = auto_hybrid
  if (present(hybrid)) self%hybrid = hybrid
  if (present(mix)) self%mix = mix
  if (present(min_cos_hp)) self%min_cos_hp = min_cos_hp
  if (present(accept_reject)) self%accept_reject = accept_reject
end subroutine init_afqmc_proj_config


!******************************************************************************** 
!
! AfqmcProjConfig constructor
!
!******************************************************************************** 
function finit_afqmc_proj_config(projection, auto_hybrid, hybrid, mix, min_cos_hp, accept_reject) result(self)
  character(len=*), optional, intent(in) :: projection
  logical, optional, intent(in)          :: auto_hybrid
  real(wp), optional, intent(in)         :: hybrid
  real(wp), optional, intent(in)         :: mix
  real(wp), optional, intent(in)         :: min_cos_hp
  logical, optional, intent(in)          :: accept_reject
  type(AfqmcProjConfig)                  :: self

  call init_afqmc_proj_config(self, projection, auto_hybrid, hybrid, mix, min_cos_hp, accept_reject)
end function finit_afqmc_proj_config


!******************************************************************************** 
!
! Read AfqmcProjConfig entries from input tags
!
!******************************************************************************** 
subroutine read_afqmc_proj_config_tags(self)
  class(AfqmcProjConfig), intent(inout) :: self

  call add_input("projection", self%projection)
  call add_input("auto_hybrid", self%auto_hybrid)
  call add_input("hybrid", self%hybrid)
  call add_input("afqmc_mix", self%mix)
  call add_input("min_cos_hp", self%min_cos_hp)
  call add_input("accept_reject", self%accept_reject)
end subroutine read_afqmc_proj_config_tags

end module afqmc_proj_config