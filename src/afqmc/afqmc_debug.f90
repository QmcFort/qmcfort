! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module afqmc_debug
!******************************************************************************** 
!
!       Impementation of the debugging files produced for the postprocessing    
!   
!******************************************************************************** 

use constants
use mpi
use file_handle, only: FileHandle
use hamilton_vars
use hamilton_type
use mc_descriptor, only: McDescriptor
use afqmc_proj, only: AfqmcProj
use afqmc_walker
use afqmc_rebalance
use array_numpy_file_io, only: ArrayNumpyFileIO

implicit none

private
public afqmc_debugger

type afqmc_debugger
  logical                         :: active 
  integer                         :: npop
  integer                         :: step
  integer                         :: step_pop

  type(McDescriptor), pointer     :: mc_des
  type(AfqmcProj), pointer        :: af_proj
  type(pop_controll), pointer     :: pop_des

  real(wp), allocatable           :: weight(:,:)                     ! W_k^w
  real(wp), allocatable           :: phase(:,:)                      ! theta_k^w
  complex(wp), allocatable        :: el(:,:)                         ! EL_k^w
  complex(wp), allocatable        :: eh(:,:)                         ! EH_k^w
  complex(wp), allocatable        :: overlap(:,:)                    ! S_k^w
  complex(wp), allocatable        :: isf(:,:)                        ! Ik^w = p_target(x) / p_proposal(x)
  integer, allocatable            :: rebalance(:,:)                  ! rebalancing of walkers
contains
  procedure                :: npops
  procedure                :: init => init_debugger
  procedure                :: reader => debugger_reader
  procedure                :: update => update_debugger
  procedure                :: update_rebalance => update_debugger_rebalance
  procedure                :: write_files => write_debug_files
end type afqmc_debugger

contains

!******************************************************************************** 
!
! Number of rebalancing calls
!
!******************************************************************************** 
integer function npops(self)
  class(afqmc_debugger), intent(in) :: self

  npops = self%mc_des%steps_tot() / self%pop_des%freq
end function npops


!******************************************************************************** 
!
! Initialize afqmc debugger
! 
!    nsteps   - total number of steps including the equilibration phase
!    nwalkers - number of walkers per MPI process
!
!******************************************************************************** 
subroutine init_debugger(self, mc_des, af_proj, pop_des, ham)
  class(afqmc_debugger), intent(inout)   :: self
  type(McDescriptor), target, intent(in) :: mc_des
  type(AfqmcProj), target, intent(in)    :: af_proj
  type(pop_controll), target, intent(in) :: pop_des
  type(Hamilton), intent(in)             :: ham
  !local variables
  integer                                :: nsteps
  integer                                :: nwalkers
  integer                                :: g

  self%active = .false.
  self%step = 0
  self%step_pop = 0

  self%mc_des => mc_des
  self%af_proj => af_proj
  self%pop_des => pop_des

  call self%reader

  if (.not. self%active) return

  nsteps = self%mc_des%steps_tot()
  nwalkers = self%mc_des%walkers_per_rank

  allocate(self%weight(nsteps, nwalkers))
  allocate(self%phase(nsteps, nwalkers))
  allocate(self%el(nsteps, nwalkers))
  allocate(self%eh(nsteps, nwalkers))
  allocate(self%overlap(nsteps, nwalkers))
  allocate(self%isf(nsteps, nwalkers))
end subroutine init_debugger


!******************************************************************************** 
!
! Read input data for afqmc debugger
!
!******************************************************************************** 
subroutine debugger_reader(self)
  use qmcfort_in, only: add_input
  class(afqmc_debugger), intent(inout) :: self

  call add_input("afqmc_debug", self%active)
end subroutine debugger_reader


!******************************************************************************** 
!
! Update debugger arrays
!
!******************************************************************************** 
subroutine update_debugger(self, af_walk)
  class(afqmc_debugger), intent(inout) :: self
  type(AfqmcWalker), intent(in)        :: af_walk(:)
  !local varaibles
  integer                              :: step, g

  self%step = self%step + 1
 
  self%weight(self%step,:) = af_walk(:)%weight
  self%phase(self%step,:) = af_walk(:)%phase
  self%el(self%step,:) = af_walk(:)%energy%electron_energy()
  self%eh(self%step,:) = af_walk(:)%energyh
  self%overlap(self%step,:) = af_walk(:)%doverlap
  self%isf(self%step,:) = af_walk(:)%isf
end subroutine update_debugger


!******************************************************************************** 
!
! Update walker rebalance arrays
!
!******************************************************************************** 
subroutine update_debugger_rebalance(self, map)
  class(afqmc_debugger), intent(inout) :: self
  integer, intent(in)                  :: map(:,:)

  self%step_pop = self%step_pop + 1

  if (.not. allocated(self%rebalance)) allocate(self%rebalance(self%npops(), size(map)))
  self%rebalance(self%step_pop,:) = reshape(map, shape=[size(map)])
end subroutine update_debugger_rebalance


!******************************************************************************** 
!
! Gather data to the master node and write them to the files
!
!******************************************************************************** 
subroutine write_debug_files(self, comm, enuc, h0)
  class(afqmc_debugger), intent(inout) :: self
  type(mpi_communicator), intent(in)   :: comm
  real(wp), intent(in)                 :: enuc
  real(wp), intent(in)                 :: h0
  !local
  character(len=*), parameter          :: dir_name = "afqmc_debugger/"
  real(wp), allocatable                :: real_data(:,:)
  complex(wp), allocatable             :: cmplx_data(:,:)
  type(ArrayNumpyFileIO)               :: numpy_io

  if (comm%mpirank == 0) call system("mkdir -p " // dir_name)

  !write config yaml file
  if (comm%mpirank == 0) call dump_yaml_file(dir_name//"input.yaml", self, enuc, h0)

  !gather weights and write them to the file
  call comm%gather(self%weight, real_data, root=0)
  if (comm%mpirank == 0) call numpy_io%save(dir_name // "w_file", real_data)

  !gather phases and write them to the file
  call comm%gather(self%phase, real_data, root=0)
  if (comm%mpirank == 0) call numpy_io%save(dir_name // "p_file", real_data)
  
  !gather local energies and write them to the file
  call comm%gather(self%el, cmplx_data, root=0)
  if (comm%mpirank == 0) call numpy_io%save(dir_name // "el_file", cmplx_data)
  
  !gather hybrid energies and write them to the file
  call comm%gather(self%eh, cmplx_data, root=0)
  if (comm%mpirank == 0) call numpy_io%save(dir_name // "eh_file", cmplx_data)
  
  !gather overlaps and write them to the file
  call comm%gather(self%overlap, cmplx_data, root=0)
  if (comm%mpirank == 0) call numpy_io%save(dir_name // "s_file", cmplx_data)
  
  !gather isf and write them to the file
  call comm%gather(self%isf, cmplx_data, root=0)
  if (comm%mpirank == 0) call numpy_io%save(dir_name // "isf_file", cmplx_data)
  
  !write pop to the file
  if (comm%mpirank == 0) call numpy_io%save(dir_name // "pop_file", self%rebalance)
end subroutine write_debug_files


!******************************************************************************** 
!
! Write yaml file for python/en_
!
!debug: yaml - temporary solution until proper yaml implementation
!******************************************************************************** 
subroutine dump_yaml_file(fname, af_debug, enuc, h0)
  character(len=*), intent(in)     :: fname
  type(afqmc_debugger), intent(in) :: af_debug
  real(wp), intent(in)             :: enuc
  real(wp), intent(in)             :: h0
  !local variables
  type(FileHandle)                 :: yaml_file

  yaml_file = FileHandle(fname)
  call yaml_file%open(status="replace", action="write")

  write(yaml_file%funit,102) "tau: ", af_debug%mc_des%tau
  write(yaml_file%funit,101) "nblocks: ", af_debug%mc_des%nblocks
  write(yaml_file%funit,101) "eqblocks: ", af_debug%mc_des%eqblocks
  write(yaml_file%funit,101) "steps_per_block: ", af_debug%mc_des%steps_per_block
  write(yaml_file%funit,100) "projection: ", af_debug%af_proj%cfg%projection
  write(yaml_file%funit,101) "pop_freq: ", af_debug%pop_des%freq
  write(yaml_file%funit,102) "enuc: ", enuc
  write(yaml_file%funit,102) "ecore: ", h0

  call yaml_file%close()

  100 format(1x, a, a)
  101 format(1x, a, i8)
  102 format(1x, a, f12.6)
end subroutine dump_yaml_file

end module afqmc_debug