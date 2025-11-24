! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module afqmc_chkpt
!******************************************************************************** 
!       
! Implements routines for AFQMC checkpoints
!
!    write/read AfqmcWalker array
!    write/read AfqmcEnergy object
!
!******************************************************************************** 


#include "../preproc.inc"
use constants
use file_handle, only: FileHandle
use profiling
use energy_types
use mc_descriptor, only: McDescriptor
use hamilton_vars, only: hamil_descriptor
use afqmc_walker
use afqmc_energy

implicit none

contains

!******************************************************************************** 
!
! Read afqmc checkpoint file:
!      
!    AfqmcWalker array and AfqmcEnergy type are read from the file
!
!******************************************************************************** 
subroutine read_afqmc_chkpt(hdes, mc_des, af_en, walkers, comm)
  type(hamil_descriptor), intent(in)           :: hdes
  type(McDescriptor), intent(inout)            :: mc_des
  type(AfqmcEnergy), intent(inout)             :: af_en
  type(AfqmcWalker), allocatable, intent(out)  :: walkers(:)
  type(mpi_communicator), intent(inout)        :: comm
  !local variables
  integer                                      :: n, nel, ng, w
  type(FileHandle)                             :: fh
  type(AfqmcWalker), allocatable               :: walkers_tot(:)
  character(len=*), parameter                  :: proc_name = "read_afqmc_chkpt"
  
  if (use_profiling) call start_profiling(proc_name)

  if (comm%mpirank == 0) then
    fh = FileHandle("afqmc_chkpt")
    call fh%open(status="old", access="stream", form="unformatted", action="read")
    
    read(fh%funit) mc_des%nwalkers, n, nel, ng
    call mc_des%adjust_walkers()
    call alloc_walker(walkers_tot, hdes, mc_des%nwalkers)
    do w = 1, mc_des%nwalkers
      call walkers_tot(w)%read_(fh)
    end do

    call af_en%from_file(fh)
    
    call fh%close()
  end if
  
  call dist_afqmc_data_r(hdes, mc_des, af_en, walkers_tot, walkers, comm)
  call dealloc_walker(walkers_tot)

  if (use_profiling) call end_profiling(proc_name)
end subroutine read_afqmc_chkpt


!******************************************************************************** 
!
! Write afqmc checkpoint file:
!      
!    AfqmcWalker array and AfqmcEnergy object are written to the file
!
!******************************************************************************** 
subroutine write_afqmc_chkpt(hdes, mc_des, af_en, walkers, comm)
  type(hamil_descriptor), intent(in)    :: hdes
  type(McDescriptor), intent(in)        :: mc_des
  type(AfqmcEnergy), intent(inout)      :: af_en
  type(AfqmcWalker), intent(in)         :: walkers(:)
  type(mpi_communicator), intent(inout) :: comm
  !local
  integer                               :: w
  logical                               :: lreturn
  integer, allocatable                  :: frac(:)
  type(FileHandle)                      :: fh
  type(AfqmcWalker), allocatable        :: walkers_tot(:)
  character(len=*), parameter           :: proc_name = "write_afqmc_chkpt"
  
  allocate(frac(4))
  frac = [25, 50, 75, 100]
  frac = frac * (mc_des%steps_tot()-mc_des%steps_eq()) / 100
  where (frac == 0) frac = unset

  lreturn = .not. mc_des%is_eq() .or. all(mc_des%step-mc_des%steps_eq() /= frac)
  if (lreturn) return

  if (use_profiling) call start_profiling(proc_name)

  call dist_afqmc_data_w(hdes, mc_des, af_en, walkers, walkers_tot, comm)
  
  if (comm%mpirank == 0) then
    fh = FileHandle("afqmc_chkpt")
    call fh%open(status="replace", access="stream", form="unformatted", action="write")
    
    write(fh%funit) mc_des%nwalkers, size(walkers_tot(1)%coeff, 1), size(walkers_tot(1)%coeff, 2), size(walkers_tot(1)%lgmean)
    
    do w = 1, mc_des%nwalkers
      call walkers_tot(w)%write_(fh)
    end do
    
    call af_en%to_file(fh)

    call fh%close()
  end if
  
  call dealloc_walker(walkers_tot)

  if (use_profiling) call end_profiling(proc_name)
end subroutine write_afqmc_chkpt


!******************************************************************************** 
!
! Distribute data over the nodes on the input/read
!
!******************************************************************************** 
subroutine dist_afqmc_data_r(hdes, mc_des, af_en, walkers_tot, walkers, comm)
  type(hamil_descriptor), intent(in)            :: hdes
  type(McDescriptor), intent(inout)             :: mc_des
  type(AfqmcEnergy), intent(inout)              :: af_en
  type(AfqmcWalker), allocatable, intent(inout) :: walkers_tot(:)
  type(AfqmcWalker), allocatable, intent(inout) :: walkers(:)
  type(mpi_communicator), intent(inout)         :: comm

  call mc_des%mpi_bcast(comm, root=0)

  !note mc_des must be broadcasted, since it is a pointer within af_en
  call af_en%mpi_bcast(comm, root=0)

  call alloc_walker(walkers, hdes, mc_des%walkers_per_rank)
  call scatter_walkers(comm, walkers_tot, walkers)
end subroutine dist_afqmc_data_r


!******************************************************************************** 
!
! Distribute data over the nodes on the output/write
!
!******************************************************************************** 
subroutine dist_afqmc_data_w(hdes, mc_des, af_en, walkers, walkers_tot, comm)
  type(hamil_descriptor), intent(in)          :: hdes
  type(McDescriptor), intent(in)              :: mc_des
  type(AfqmcEnergy), intent(inout)           :: af_en
  type(AfqmcWalker), intent(in)               :: walkers(:)
  type(AfqmcWalker), allocatable, intent(out) :: walkers_tot(:)
  type(mpi_communicator), intent(inout)       :: comm

  call gather_walkers(comm, hdes, walkers, walkers_tot, root=0)
end subroutine dist_afqmc_data_w

end module afqmc_chkpt