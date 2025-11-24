! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module afqmc_intags
!******************************************************************************** 
!
! Type to store and validate AFQMC input tags
!
! It also defines default 
!
!******************************************************************************** 

use constants
use mpi, only: comm_world
use qmcfort_in, only: add_input

implicit none

private
public :: AfqmcTags

type AfqmcTags
  !AfqmcManager related variables
  logical                :: g_resolved
  logical                :: sp_proj
  logical                :: dump_walker
  integer                :: local_energy_test

  !McDescriptor related variables
  real(wp)               :: tau
  integer                :: nwalkers
  integer                :: steps_per_block
  integer                :: nblocks
  integer                :: eqblocks
  integer                :: reorth
  integer                :: sample
  integer                :: samplex

  !ImportanceSampler related variables
  character(len=charlen) :: is_shift
  character(len=charlen) :: is_scale
  logical                :: subtract_mean_field
  character(len=charlen) :: cut_shift
  character(len=charlen) :: shift_cutoff

  !AuxFieldGen related variables
  character(len=charlen) :: aux_field_dist
  integer                :: nfields_max

  !WaveTrialDes related variables
  character(len=charlen) :: wave_trial_request
  character(len=charlen) :: exchange_mode
  integer                :: h_block_size
  integer                :: x_block_size
  logical                :: compress_x
  real(wp)               :: compress_x_tol

  !growth_estimator related variables
  logical                :: growth_estimator
contains
  procedure :: reader => read_afqmc_tags
end type AfqmcTags

interface AfqmcTags
  module procedure finit_afqmc_tags
end interface AfqmcTags

contains 

!******************************************************************************** 
!
! Initialization of the AfqmcTags object
!
!******************************************************************************** 
subroutine init_afqmc_tags(self)
  type(AfqmcTags), intent(out) :: self

  call self%reader()
end subroutine init_afqmc_tags


!******************************************************************************** 
!
! afqmc_tags constructor 
!
!******************************************************************************** 
function finit_afqmc_tags() result(self)
  type(AfqmcTags) :: self

  call init_afqmc_tags(self)
end function finit_afqmc_tags


!******************************************************************************** 
!
! Read and validate AFQMC associated tags 
!
!    1. values specified in the input file are read to the dedicated typeM
!    2. before read, default values are initialized
!
! Implementation follows 3 classes of the tags:
!
!    1. mandatory tags : program stops if the tags are omitted
!    2. important tags : warnings raised if the tags are omitted
!    3. normal tags    : default value will be used silently
!
!******************************************************************************** 
subroutine read_afqmc_tags(self)
  class(AfqmcTags), intent(inout) :: self
  !local variables
  character(len=:), allocatable   :: tag
  logical                         :: found

  tag = "tau"
  call add_input(tag, self%tau, found)
  if (.not. found) then
    if(comm_world%mpirank == 0) write(*,100) tag
    call comm_world%abort()
    stop
  end if

  tag = "nwalkers"
  call add_input(tag, self%nwalkers, found)
  if (.not. found) then
    if(comm_world%mpirank == 0) write(*,100) tag
    call comm_world%abort()
    stop
  end if

  tag = "steps_per_block"
  call add_input(tag, self%steps_per_block, found)
  if (.not. found) then
    if(comm_world%mpirank == 0) write(*,100) tag
    call comm_world%abort()
    stop
  end if

  tag = "nblocks"
  call add_input(tag, self%nblocks, found)
  if (.not. found) then
    if(comm_world%mpirank == 0) write(*,100) tag
    call comm_world%abort()
    stop
  end if

  tag = "eqblocks"
  self%eqblocks = 0
  call add_input(tag, self%eqblocks, found)
  if (.not. found) then
    if(comm_world%mpirank == 0) write(*,101) tag, "0"
  end if

  tag = "reorth"
  self%reorth = 1
  call add_input(tag, self%reorth, found)
  if (.not. found) then
    if(comm_world%mpirank == 0) write(*,101) tag, "1"
  end if

  tag = "subtract_mean_field"
  self%subtract_mean_field = .true.
  call add_input(tag, self%subtract_mean_field, found)

  tag = "importance_sampling_shift"
  self%is_shift = "mix"
  call add_input(tag, self%is_shift, found)

  tag = "importance_sampling_scale"
  self%is_scale = "none"
  call add_input(tag, self%is_scale, found)

  if (self%is_scale == "mix") then
    self%g_resolved = .true.
  else
    self%g_resolved = .false.
  end if

  tag = "cut_shift"
  self%cut_shift = "remove"
  call add_input(tag, self%cut_shift, found)

  tag = "shift_cutoff"
  self%shift_cutoff = "MZ"
  call add_input(tag, self%shift_cutoff, found)

  tag = "growth_estimator"
  self%growth_estimator = .false.
  call add_input(tag, self%growth_estimator, found)

  tag = "local_energy_test"
  self%local_energy_test = 0
  call add_input(tag, self%local_energy_test, found)

  tag = "dump_walker"
  self%dump_walker = .false.
  call add_input(tag, self%dump_walker, found)

  tag = "spin_projected_walkers"
  self%sp_proj = .false.
  call add_input(tag, self%sp_proj, found)

  tag = "sample"
  self%sample = 1
  call add_input(tag, self%sample, found)
  if (self%g_resolved) then
    self%sample = 1
    if (comm_world%mpirank == 0) write(102,*) "tags: sample = 1 required with g_resolved = True"
  else 
    if (.not. found) then
      if(comm_world%mpirank == 0) write(*,101) tag, "1"
    end if
  end if

  tag = "samplex"
  self%samplex = 1
  call add_input(tag, self%samplex, found)
  if (self%g_resolved) then
    self%sample = 1
    if (comm_world%mpirank == 0) write(102,*) "tags: samplex = 1 required with g_resolved = True"
  else 
    if (.not. found) then
      if(comm_world%mpirank == 0) write(*,101) tag, "1"
    end if
  end if

  tag = "wave_trial"
  self%wave_trial_request = "none"
  call add_input(tag, self%wave_trial_request, found)

  tag = "exchange_mode"
  self%exchange_mode = "cholesky"
  call add_input(tag, self%exchange_mode, found)
  if (self%g_resolved) then
    self%exchange_mode = "cholesky"
    if (comm_world%mpirank == 0) write(102,*) "tags: exchange_mode = 'cholesky' required with g_resolved = True"
  end if

  tag = "compress_x"
  self%compress_x = .false.
  call add_input(tag, self%compress_x, found)
  if (self%compress_x .and. self%exchange_mode == "eri") then
    self%compress_x = .false.
    if (comm_world%mpirank == 0) write(102,*) "tags: compress_x = .true. can not be used with exchange_mode = eri"
  end if

  tag = "compress_x_tol"
  self%compress_x_tol = 1.0e-06_wp
  call add_input(tag, self%compress_x_tol, found)

  tag = "h_block_size"
  self%h_block_size = self%nwalkers / comm_world%totsize
  call add_input(tag, self%h_block_size, found)

  tag = "x_block_size"
  select case (trim(self%exchange_mode))
    case ("cholesky")
      self%x_block_size = 5
    case ("eri")
      self%x_block_size = self%h_block_size
    case default
      if (comm_world%mpirank == 0) write(*,*) "exchange_mode not recognized: program stops now"
      call exit
  end select
  call add_input(tag, self%x_block_size, found)

  tag = "aux_field_dist"
  self%aux_field_dist = "normal"
  call add_input(tag, self%aux_field_dist, found)

  tag = "nfields_max"
  self%nfields_max = unset   
  call add_input(tag, self%nfields_max, found)

  100 format(1x,"tags: mandatory tag ", a, " not specified. Please add it to the qmcfort_in file!")
  101 format(1x,"tags: tag ", a, " not specified. Default value ", a, " will be used!")
  102 format(1x,a)
end subroutine read_afqmc_tags

end module afqmc_intags