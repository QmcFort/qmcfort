! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module mc_descriptor
!******************************************************************************** 
!       
! Type to store / control Monte Carlo variables
!
!******************************************************************************** 

#include "preproc.inc"
use constants
use mpi

implicit none

private
public McDescriptor

type McDescriptor
  real(wp)    :: tau                        !Monte Carlo step size
  integer     :: nwalkers                   !Number of walkers
  integer     :: steps_per_block            !Number of steps per block
  integer     :: nblocks                    !Number of sampling blocks
  integer     :: eqblocks                   !Number of equilibration blocks
  integer     :: reorth                     !Frequency of the reorthogonalization frequency
  integer     :: sample                     !Frequency of the energy sampling
  integer     :: samplex                    !Frequency of the exchange energy sampling

  integer     :: walkers_per_rank           !Number of walkers per MPI rank
  integer     :: print_screen
  complex(wp) :: isqtau                     !sqrt(-tau)

  integer     :: step                       !Current Step
  integer     :: block                      !Current block
contains
  procedure   :: increment => increment_mc_descriptor
  procedure   :: is_running => is_mc_running

  procedure   :: steps_eq => mc_steps_eq
  procedure   :: steps_samp => mc_steps_samp
  procedure   :: steps_tot => mc_steps_total
  procedure   :: steps_sample => mc_steps_sample
  procedure   :: steps_samplex => mc_steps_samplex

  procedure   :: time_eq => get_equil_time 
  procedure   :: time_samp => get_sampling_time
  procedure   :: time_tot => get_total_time 
  procedure   :: time_current => get_current_time

  procedure   :: is_sample => is_mc_sample_step
  procedure   :: is_samplex => is_mc_samplex_step
  procedure   :: is_block_end => is_mc_block_end
  procedure   :: block_range => mc_block_range

  procedure   :: is_eq => is_equilibrated
  procedure   :: eq_finish => does_equilibration_finish

  procedure   :: set_walkers => set_mc_walkers
  procedure   :: adjust_walkers => adjust_mc_walkers

  procedure   :: set_print_screen => set_print_screen

  procedure   :: mpi_bcast => mpi_bcast_mc_descriptor
end type McDescriptor

interface McDescriptor
  module procedure init_mc_descriptor
end interface McDescriptor

contains

!******************************************************************************** 
! 
! McDescriptor constructor
!
!******************************************************************************** 
function init_mc_descriptor(nwalkers, tau, steps_per_block, nblocks, eqblocks, reorth, sample, samplex) result(self)
  integer, intent(in)  :: nwalkers
  real(wp), intent(in) :: tau
  integer, intent(in)  :: steps_per_block
  integer, intent(in)  :: nblocks
  integer, intent(in)  :: eqblocks
  integer, intent(in)  :: reorth
  integer, intent(in)  :: sample
  integer, intent(in)  :: samplex
  type(McDescriptor)  :: self

  self%tau = tau
  self%nwalkers = nwalkers
  self%steps_per_block = steps_per_block
  self%nblocks = nblocks
  self%eqblocks = eqblocks
  self%reorth = reorth
  self%sample = sample
  self%samplex = samplex

  self%step = 0
  self%block = 0

  self%isqtau = - self%tau !* onec
  self%isqtau = sqrt(self%isqtau)

  call self%adjust_walkers()

  call self%set_print_screen()
end function init_mc_descriptor


!******************************************************************************** 
!
! Increment state of the MC descriptor
!
!******************************************************************************** 
subroutine increment_mc_descriptor(self)
  class(McDescriptor), intent(inout) :: self

  self%step = self%step + 1
  if (mod(self%step, self%steps_per_block) == 1) self%block = self%block + 1
end subroutine increment_mc_descriptor


!******************************************************************************** 
!
! Is MC procedure still running
!
!******************************************************************************** 
logical function is_mc_running(self)
  class(McDescriptor), intent(in) :: self

  is_mc_running = self%step < self%steps_tot()
end function is_mc_running


!******************************************************************************** 
!
! Number of equilibration steps
!
!******************************************************************************** 
integer function mc_steps_eq(self)
  class(McDescriptor), intent(in) :: self

  mc_steps_eq = self%steps_per_block * self%eqblocks
end function mc_steps_eq


!******************************************************************************** 
!
! Number of sampling steps (number of steps after equilibration)
!
!******************************************************************************** 
integer function mc_steps_samp(self)
  class(McDescriptor), intent(in) :: self

  mc_steps_samp = self%steps_per_block * self%nblocks
end function mc_steps_samp


!******************************************************************************** 
!
! Determine total number of the MC steps
!
!******************************************************************************** 
integer function mc_steps_total(self)
  class(McDescriptor), intent(in) :: self

  mc_steps_total = self%steps_per_block * (self%eqblocks + self%nblocks)
end function mc_steps_total


!******************************************************************************** 
!
! Determine number of the MC steps where energy is sampled
!
!******************************************************************************** 
integer function mc_steps_sample(self)
  class(McDescriptor), intent(in) :: self

  mc_steps_sample = self%steps_tot() / self%sample
end function mc_steps_sample


!******************************************************************************** 
!
! Determine number of the MC steps where exchange energy is sampled
!
!******************************************************************************** 
integer function mc_steps_samplex(self)
  class(McDescriptor), intent(in) :: self

  mc_steps_samplex = self%steps_sample() / self%samplex
end function mc_steps_samplex


!******************************************************************************** 
!
! Determine equilibration time of the MC procedure
!
!******************************************************************************** 
real(wp) function get_equil_time(self)
  class(McDescriptor), intent(in) :: self

  get_equil_time = self%tau * self%steps_eq()
end function get_equil_time


!******************************************************************************** 
!
! Determine equilibration time of the MC procedure
!
!******************************************************************************** 
real(wp) function get_sampling_time(self)
  class(McDescriptor), intent(in) :: self

  get_sampling_time = self%tau * self%steps_samp()
end function get_sampling_time


!******************************************************************************** 
!
! Determine total simulation time of the MC procedure
!
!******************************************************************************** 
real(wp) function get_total_time(self)
  class(McDescriptor), intent(in) :: self

  get_total_time = self%tau * self%steps_tot()
end function get_total_time


!******************************************************************************** 
!
! Determine current simulation time of the MC procedure
!
!******************************************************************************** 
real(wp) function get_current_time(self)
  class(McDescriptor), intent(in) :: self

  get_current_time = self%tau * self%step
end function get_current_time


!******************************************************************************** 
!
! Is the energy sample step
!
!******************************************************************************** 
logical function is_mc_sample_step(self)
  use standalone, only: modul
  class(McDescriptor), intent(in) :: self

  is_mc_sample_step = modul(self%step, self%sample)
end function is_mc_sample_step


!******************************************************************************** 
!
! Is the exchange energy sample step
!
!******************************************************************************** 
logical function is_mc_samplex_step(self)
  use standalone, only: modul
  class(McDescriptor), intent(in) :: self

  is_mc_samplex_step = modul(self%step, self%sample*self%samplex)
end function is_mc_samplex_step

!******************************************************************************** 
!
! Is the end of the MC block
!
!******************************************************************************** 
logical function is_mc_block_end(self)
  use standalone, only: modul
  class(McDescriptor), intent(in) :: self

  is_mc_block_end = modul(self%step, self%steps_per_block)
end function is_mc_block_end


!******************************************************************************** 
!
! Calculate first and the last MC step in the given MC block
!
!******************************************************************************** 
function mc_block_range(self, bblock) result(block_range)
  class(McDescriptor), intent(in) :: self
  integer, intent(in)             :: bblock
  integer                         :: block_range(2)

  if (bblock == 0) then
    block_range = 0
  else
    block_range(1) = (bblock-1) * self%steps_per_block + 1
    block_range(2) = bblock * self%steps_per_block
  end if
end function mc_block_range


!******************************************************************************** 
!
! Check whether sampling is already equilibrated
!
!******************************************************************************** 
logical function is_equilibrated(self)
  class(McDescriptor), intent(in) :: self

  is_equilibrated = self%block > self%eqblocks
end function is_equilibrated


!******************************************************************************** 
!
! Check if the equilibration finishes at a current block
!
!******************************************************************************** 
logical function does_equilibration_finish(self)
  class(McDescriptor), intent(in) :: self

  does_equilibration_finish = self%block == self%eqblocks
end function does_equilibration_finish


!******************************************************************************** 
!
! Set walkers variables
!
!******************************************************************************** 
subroutine set_mc_walkers(self, nwalkers)
  class(McDescriptor), intent(inout) :: self
  integer, intent(in)                 :: nwalkers

  self%nwalkers = nwalkers
  call self%adjust_walkers()
end subroutine set_mc_walkers


!******************************************************************************** 
!
! Adjust number of walkers to be multiple of the MPI ranks
!
!******************************************************************************** 
subroutine adjust_mc_walkers(self)
  class(McDescriptor), intent(inout) :: self
  !local variables
  integer                             :: nwalkers

  nwalkers = (self%nwalkers / comm_world%mpisize) * comm_world%mpisize

  if (nwalkers/=self%nwalkers .and. comm_world%mpirank==0) then
      write(*,*) "mc_des: Number of walkers is modified from ", self%nwalkers, "to ", nwalkers
  end if

  self%nwalkers = nwalkers
  self%walkers_per_rank  = self%nwalkers / comm_world%mpisize
end subroutine adjust_mc_walkers


!******************************************************************************** 
!
! Determine the frequency to write output to std output 
!
!******************************************************************************** 
subroutine set_print_screen(self)
  class(McDescriptor), intent(inout) :: self
  
  if (self%nblocks <= 20) then
    self%print_screen = 1
  else if (self%nblocks <= 50) then
    self%print_screen = 2
  else if (self%nblocks <= 100) then
    self%print_screen = 5
  else if (self%nblocks <= 250) then
    self%print_screen = 10
  else if (self%nblocks <= 500) then
    self%print_screen = 20
  else
    self%print_screen = self%nblocks / 20 
  end if
end subroutine set_print_screen


!******************************************************************************** 
!
! MPI Bcast McDescriptor object
!
!******************************************************************************** 
subroutine mpi_bcast_mc_descriptor(self, comm, root)
  class(McDescriptor), intent(inout)  :: self
  class(mpi_communicator), intent(in) :: comm
  integer, intent(in)                 :: root

  call comm%bcast(self%step, root)
  call comm%bcast(self%block, root)

  call comm%bcast(self%tau, root)
  call comm%bcast(self%nwalkers, root)
  call comm%bcast(self%walkers_per_rank, root)
  call comm%bcast(self%steps_per_block, root)
  call comm%bcast(self%eqblocks, root)
  call comm%bcast(self%nblocks, root)
  call comm%bcast(self%reorth, root)
  call comm%bcast(self%sample, root)
  call comm%bcast(self%samplex, root)

  call comm%bcast(self%isqtau, root)
end subroutine mpi_bcast_mc_descriptor

end module mc_descriptor