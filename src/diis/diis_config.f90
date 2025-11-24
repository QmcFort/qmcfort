! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module diis_config

use constants
use qmcfort_in, only: add_input

implicit none

private
public DiisConfig

type DiisConfig
  logical :: enabled = .true.
  integer :: max_diis_size = 12
  integer :: start_sampling_step = 2
  integer :: start_solving_step = 5
contains
  procedure :: read_tags => read_diis_config_tags
  procedure :: validate => validate_diis_config
end type DiisConfig

interface DiisConfig
  module procedure finit_diis_config
end interface DiisConfig

contains

!******************************************************************************** 
!
! Initialization of the DiisConfig object
!
!******************************************************************************** 
subroutine init_diis_config(self, enabled, max_diis_size, start_sampling_step, start_solving_step)
  type(DiisConfig), intent(out) :: self
  logical, optional, intent(in) :: enabled
  integer, optional, intent(in) :: max_diis_size
  integer, optional, intent(in) :: start_sampling_step
  integer, optional, intent(in) :: start_solving_step

  if (present(enabled)) self%enabled = enabled
  if (present(max_diis_size)) self%max_diis_size = max_diis_size
  if (present(start_sampling_step)) self%start_sampling_step = start_sampling_step
  if (present(start_solving_step)) self%start_solving_step = start_solving_step
end subroutine init_diis_config


!******************************************************************************** 
!
! DiisConfig constructor
!
!******************************************************************************** 
function finit_diis_config(enabled, max_diis_size, start_sampling_step, start_solving_step) result(self)
  logical, optional, intent(in) :: enabled
  integer, optional, intent(in) :: max_diis_size
  integer, optional, intent(in) :: start_sampling_step
  integer, optional, intent(in) :: start_solving_step
  type(DiisConfig)              :: self

  call init_diis_config(self, enabled, max_diis_size, start_sampling_step, start_solving_step)
end function finit_diis_config


!******************************************************************************** 
!
! Read DiisConfig entries from input tags
!
!******************************************************************************** 
subroutine read_diis_config_tags(self)
  class(DiisConfig), intent(inout) :: self

  call add_input("use_diis", self%enabled)
  call add_input("max_diis_size", self%max_diis_size)
  call add_input("diis_start_sampling_step", self%start_sampling_step)
  call add_input("diis_start_solving_step", self%start_solving_step)
end subroutine read_diis_config_tags


!******************************************************************************** 
!
! Validate DiisConfig entries
!
!******************************************************************************** 
subroutine validate_diis_config(self)
  class(DiisConfig), intent(inout) :: self

  if (self%max_diis_size < 2) self%max_diis_size = 2
end subroutine validate_diis_config

end module diis_config