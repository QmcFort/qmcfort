! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module afqmc_descriptor
!******************************************************************************** 
!
! Contains all derived types needed in afqmc calculation
!
!******************************************************************************** 

use constants
use qmcfort_in, only: add_input
use mpi, only: mpi_communicator, comm_world
use hamilton_vars, only: hamil_descriptor
use mc_descriptor, only: McDescriptor
use afqmc_proj, only: AfqmcProj

implicit none

type AfqmcDescriptor
  type(McDescriptor), pointer :: mc_des
  type(mpi_communicator)      :: comm

  logical                     :: g_resolved
  logical                     :: sp_proj
  logical                     :: dump_walker
  integer                     :: local_energy_test
end type AfqmcDescriptor

interface AfqmcDescriptor
  procedure :: init_afqmc_descriptor
end interface AfqmcDescriptor

contains

!******************************************************************************** 
!
! AfqmcDescriptor constructor
!
!******************************************************************************** 
function init_afqmc_descriptor(hdes, mc_des, g_resolved, sp_proj, dump_walker, local_energy_test) result(self)
!debug:
!shouldnt have side effects on hdes
  type(hamil_descriptor), intent(inout)  :: hdes
!debug: these two should probably be POINTERS NOT TARGETS
  type(McDescriptor), target, intent(in) :: mc_des
  logical, intent(in)                    :: g_resolved
  logical, intent(in)                    :: sp_proj
  logical, intent(in)                    :: dump_walker
  integer, intent(in)                    :: local_energy_test
  type(AfqmcDescriptor)                  :: self
  
  self%mc_des => mc_des
  self%comm = comm_world

  call afqmc_reader(self, hdes)

  self%g_resolved = g_resolved
  self%sp_proj = sp_proj
  self%dump_walker = dump_walker
  self%local_energy_test = local_energy_test
end function init_afqmc_descriptor


!******************************************************************************** 
!
! AfqmcDescriptor reader 
!
!debug: should be removed on the long run
!******************************************************************************** 
subroutine afqmc_reader(self, ham_des)
  use qmcfort_in, only: add_input
  type(AfqmcDescriptor), intent(inout) :: self
  type(hamil_descriptor), intent(inout) :: ham_des
  
  call add_input("threshold_spars",     ham_des%tol_spars)
  call add_input("tol_spars",           ham_des%tol_spars)
end subroutine afqmc_reader

end module afqmc_descriptor