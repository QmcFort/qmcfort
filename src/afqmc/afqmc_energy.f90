! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module afqmc_energy
!******************************************************************************** 
!
!    The AfqmcEnergy type manages the calculation of various energy averages
! in AFQMC simulations. It stores time-averaged weights, local and hybrid
! energies, as well as intermediate quantities required for computing 
! standard deviations.
!
!    At each time step, after computing walker energies, the routine
! AfqmcEnergy%update(walkers, af_projection) is called to accumulate
! walker-averaged quantities and store them in internal buffers.
!
!    The type supports multiple methods for computing time averages,
! including block and running averages, for quantities such as local
! and hybrid energies, and the average phase/sign.
!
!******************************************************************************** 

#include "../preproc.inc"
use constants
use mpi
use profiling
use file_handle, only: FileHandle
use statistics
use energy_types
use mc_descriptor, only: McDescriptor
use afqmc_walker, only: AfqmcWalker
use afqmc_proj, only: AfqmcProj

implicit none

private
public AfqmcEnergy

type AfqmcEnergy
  type(McDescriptor), pointer :: mc_des
  type(mpi_communicator)      :: comm
  type(CEnergy)               :: e_init

  real(wp), allocatable       :: weight(:)
  complex(wp), allocatable    :: sample_weight(:)
  type(CEnergy), allocatable  :: el(:)
  complex(wp), allocatable    :: eh(:)
  type(CEnergy), allocatable  :: el_nw(:)
  complex(wp), allocatable    :: eh_nw(:)
  complex(wp), allocatable    :: el_var(:)
  complex(wp), allocatable    :: eh_var(:)

  !internal counter to know what is already populated
  integer                     :: step__

  !fast cumulative estiamtes
  complex(wp)                 :: sample_weight__
  type(CEnergy)               :: el__
  complex(wp)                 :: eh__
contains
  procedure :: sizeof => sizeof_afqmc_energy
  procedure :: alloc_arrays => allocate_afqmc_energy_arrays
  procedure :: setup_cumulatives => setup_afqmc_energy_cumulatives
  procedure :: update => update_afqmc_energy

  procedure :: start_indx => afqmc_energy_start_indx

  procedure :: last_sample_weight => afqmc_last_sample_weight
  procedure :: block_sample_weight => afqmc_block_sample_weight

  procedure :: block_sign => afqmc_block_sign
  procedure :: mean_sign => afqmc_mean_sign

  procedure :: block_el => afqmc_block_el
  procedure :: block_eh => afqmc_block_eh

  procedure :: block_el_sd => afqmc_block_el_sd
  procedure :: block_eh_sd => afqmc_block_eh_sd

  procedure :: mean_el => afqmc_mean_el
  procedure :: mean_eh => afqmc_mean_eh

  procedure :: mean_el_sd => afqmc_mean_el_sd
  procedure :: mean_eh_sd => afqmc_mean_eh_sd

  procedure :: mean_el_sdm => afqmc_mean_el_sdm
  procedure :: mean_eh_sdm => afqmc_mean_eh_sdm

  procedure :: mean_el_nw => afqmc_mean_el_nw
  procedure :: mean_eh_nw => afqmc_mean_eh_nw

  procedure :: mean_el_sdm_nw => afqmc_mean_el_sdm_nw
  procedure :: mean_eh_sdm_nw => afqmc_mean_eh_sdm_nw

  procedure :: mean_el_tnw => afqmc_mean_el_tnw
  procedure :: mean_eh_tnw => afqmc_mean_eh_tnw

  procedure :: mean_el_sdm_tnw => afqmc_mean_el_sdm_tnw
  procedure :: mean_eh_sdm_tnw => afqmc_mean_eh_sdm_tnw

  procedure :: fast_mean_el => afqmc_fast_mean_el
  procedure :: fast_mean_eh => afqmc_fast_mean_eh

  procedure :: from_file => read_afqmc_energy_from_file
  procedure :: to_file => write_afqmc_energy_to_file
  procedure :: mpi_bcast => mpi_bcast_afqmc_energy
end type AfqmcEnergy

interface AfqmcEnergy
  module procedure init_afqmc_energy
end interface AfqmcEnergy

interface downsample_array
  module procedure downsample_array_cmplx, downsample_array_energy
end interface downsample_array

contains

!******************************************************************************** 
!
! AfqmcEnergy constructor
!
!******************************************************************************** 
function init_afqmc_energy(mc_des, comm, e_init) result(self)
  type(McDescriptor), pointer, intent(in) :: mc_des
  type(mpi_communicator), intent(in)      :: comm
  type(CEnergy), intent(in)               :: e_init
  type(AfqmcEnergy)                       :: self
  !local 
  integer                                 :: i

  self%mc_des => mc_des
  self%comm = comm
  self%e_init = e_init

  self%step__ = 0
  
  call self%alloc_arrays()
  call self%setup_cumulatives()
end function init_afqmc_energy


!******************************************************************************** 
!
! Memory size of the AfqmcEnergy instance
!
!******************************************************************************** 
function sizeof_afqmc_energy(self) result(mem)
  class(AfqmcEnergy), intent(in) :: self
  real(wp)                       :: mem

  mem = sizeof(self) + sizeof(self%weight) + sizeof(self%sample_weight) 
  mem = mem + size(self%el) * sizeof(self%el(1)) + size(self%el_nw) * sizeof(self%el_nw(1))
  mem = mem + sizeof(self%eh) + sizeof(self%eh_nw) + sizeof(self%el_var) + sizeof(self%eh_var)
end function sizeof_afqmc_energy


!******************************************************************************** 
!
! Allocate AfqmcEnergy arrays
!
! 0-th element is calculated using initial set of walkers
!
!******************************************************************************** 
subroutine allocate_afqmc_energy_arrays(self)
  class(AfqmcEnergy), intent(inout) :: self
  !local
  integer                           :: steps_tot

  steps_tot = self%mc_des%steps_tot()

  if (allocated(self%sample_weight)) deallocate(self%sample_weight)
  allocate(self%sample_weight(0:steps_tot))

  if (allocated(self%weight)) deallocate(self%weight)
  allocate(self%weight(0:steps_tot))

  if (allocated(self%el)) deallocate(self%el)
  allocate(self%el(0:steps_tot))

  if (allocated(self%eh)) deallocate(self%eh)
  allocate(self%eh(0:steps_tot))

  if (allocated(self%el_nw)) deallocate(self%el_nw)
  allocate(self%el_nw(0:steps_tot))

  if (allocated(self%eh_nw)) deallocate(self%eh_nw)
  allocate(self%eh_nw(0:steps_tot))

  if (allocated(self%el_var)) deallocate(self%el_var)
  allocate(self%el_var(0:steps_tot))

  if (allocated(self%eh_var)) deallocate(self%eh_var)
  allocate(self%eh_var(0:steps_tot))
end subroutine allocate_afqmc_energy_arrays


!******************************************************************************** 
!
! Recalculate cumulative quantites that are used for fast mean energies
!
!******************************************************************************** 
subroutine setup_afqmc_energy_cumulatives(self)
  class(AfqmcEnergy), intent(inout) :: self
  !local
  integer                           :: i

  self%sample_weight__ = zeroc
  self%eh__ = zeroc
  self%el__ = CEnergy(enuc=self%e_init%enuc, e0=self%e_init%e0)

  do i = 1, self%mc_des%step
    self%sample_weight__ = self%sample_weight__ + self%sample_weight(i)
    self%eh__ = self%eh__ + self%sample_weight(i) * self%eh(i)
    self%el__ = self%el__ + self%sample_weight(i) * self%el(i)
  end do
end subroutine setup_afqmc_energy_cumulatives


!******************************************************************************** 
!
! Update AfqmcEnergy arrays after walker propagation
!
!    1. calculate weighted and non-weighted walker averages for diff. observables
!    2. MPI synchronization after each time step
!
!******************************************************************************** 
subroutine update_afqmc_energy(self, walkers, af_proj)
  class(AfqmcEnergy), intent(inout) :: self
  type(AfqmcWalker), intent(in)     :: walkers(:)
  type(AfqmcProj), intent(in)       :: af_proj
  !local                  
  integer                           :: w, step
  real(wp)                          :: weight_sum
  complex(wp)                       :: sample_weight, sample_weight_sum
  complex(wp)                       :: sample_eh, eh, eh_nw, el_var, eh_var, mu
  type(CEnergy)                     :: sample_el, el, el_nw
  character(len=*), parameter       :: proc_name = "update_afqmc_energy"

  if (use_profiling) call start_profiling(proc_name)

  mu = self%e_init%electron_energy()

  weight_sum = zeror
  sample_weight_sum = zeroc

  el = CEnergy(enuc=self%e_init%enuc, e0=self%e_init%e0)
  eh = zeroc

  el_nw = CEnergy(enuc=self%e_init%enuc, e0=self%e_init%e0)
  eh_nw = zeroc

  el_var = zeroc
  eh_var = zeroc

  do w = 1, size(walkers)
    sample_weight = af_proj%sample_weight(walkers(w))
    weight_sum = weight_sum + walkers(w)%weight
    sample_weight_sum = sample_weight_sum + sample_weight

    sample_eh = af_proj%sample_energy_val(walkers(w)%energyh)
    eh = eh + sample_weight * sample_eh
    eh_nw = eh_nw + sample_eh
    eh_var = eh_var + sample_weight * abs(sample_eh - mu)**2

    if (self%mc_des%is_sample()) then
      sample_el = af_proj%sample_energy_type(walkers(w)%energy)
      el = el + sample_weight * sample_el
      el_nw = el_nw + sample_el
      el_var = el_var + sample_weight * abs(sample_el%electron_energy() - mu)**2
    end if
  end do

  call self%comm%mpisum(weight_sum)
  call self%comm%mpisum(sample_weight_sum)
  call self%comm%mpisum(eh)
  call self%comm%mpisum(eh_nw)
  call el%mpisum(self%comm)
  call el_nw%mpisum(self%comm)
  call self%comm%mpisum(el_var)
  call self%comm%mpisum(eh_var)

  if (self%mc_des%is_samplex()) then
    self%sample_weight__ = self%sample_weight__ + sample_weight_sum
    self%el__ = self%el__ + el
    self%eh__ = self%eh__ + eh
  end if

  el = el / sample_weight_sum
  el_nw = el_nw / real(self%mc_des%nwalkers, kind=wp)
  eh = eh / sample_weight_sum
  eh_nw = eh_nw / self%mc_des%nwalkers
  el_var = el_var / sample_weight_sum
  eh_var = eh_var / sample_weight_sum

  !sync internal AfqmcEnergy counter with McDescriptor
  self%step__ = self%mc_des%step
  step = self%step__

  self%weight(step) = weight_sum
  self%sample_weight(step) = sample_weight_sum
  self%el(step) = el
  self%el_nw(step) = el_nw
  self%eh(step) = eh
  self%eh_nw(step) = eh_nw
  self%el_var(step) = el_var
  self%eh_var(step) = eh_var

  if (use_profiling) call end_profiling(proc_name)
end subroutine update_afqmc_energy


!******************************************************************************** 
!
! Determine the starting index for evaluating time averages
!
!******************************************************************************** 
integer function afqmc_energy_start_indx(self) result(start_indx)
  class(AfqmcEnergy), intent(in) :: self

  if (self%mc_des%is_eq()) then
    start_indx = self%mc_des%steps_eq() + 1
  else
    start_indx = 1
  end if
end function afqmc_energy_start_indx


!******************************************************************************** 
!
! Calculate sample weight (summed over walkers) for the last MC step
!
!    \sum{w} W_{kw} / nwalkers
!
!******************************************************************************** 
function afqmc_last_sample_weight(self) result(last_weight)
  class(AfqmcEnergy), intent(in) :: self
  complex(wp)                    :: last_weight
  !local
  integer                        :: step

  step = self%step__
  last_weight = self%sample_weight(step) / self%mc_des%nwalkers
end function afqmc_last_sample_weight


!******************************************************************************** 
!
! Calculate average sample weight within the selected MC block
!
!    \sum_{k in B} W_k = \sum{k in B, w} W_{kw} / (nwalkers * steps_per_block) 
!
!******************************************************************************** 
function afqmc_block_sample_weight(self) result(block_weight)
  class(AfqmcEnergy), intent(in) :: self
  complex(wp)                    :: block_weight
  !local
  integer                        :: bblock, block_range(2)
  
  bblock = self%step__ / self%mc_des%steps_per_block
  block_range = self%mc_des%block_range(bblock)

  if (bblock == 0) then
    block_weight = self%sample_weight(0) / self%mc_des%nwalkers
  else
    block_weight = sum(self%sample_weight(block_range(1):block_range(2)))
    block_weight = block_weight / self%mc_des%steps_per_block / self%mc_des%nwalkers
  end if
end function afqmc_block_sample_weight


!******************************************************************************** 
!
! Estimate block local energy for the current MC block
!
!    E_B = \sum_{k in B, w} W_{kw} E_{kw} / \sum_{k in B, w} W_{kw}
!
!******************************************************************************** 
function afqmc_block_el(self) result(block_el)
  class(AfqmcEnergy), intent(in) :: self
  type(CEnergy)                  :: block_el
  !local
  integer                        :: bblock, block_range(2)
  complex(wp), allocatable       :: weight(:), weight_x(:)
  type(CEnergy), allocatable     :: en(:), en_x(:)

  bblock = self%step__ / self%mc_des%steps_per_block
  block_range = self%mc_des%block_range(bblock)

  if (bblock == 0) then
    block_el = self%el(0)
  else
    weight = downsample_array(self%sample_weight(block_range(1):block_range(2)), self%mc_des%sample)
    weight_x = downsample_array(weight, self%mc_des%samplex)

    en = downsample_array(self%el(block_range(1):block_range(2)), self%mc_des%sample)
    en_x = downsample_array(en, self%mc_des%samplex)

    block_el = CEnergy(enuc=self%e_init%enuc, e0=self%e_init%e0)
    block_el%e1 = mean(weight, en%e1)
    block_el%eh = mean(weight, en%eh)
    block_el%ex = mean(weight_x, en_x%ex)
  end if
end function afqmc_block_el


!******************************************************************************** 
!
! Estimate block hybrid energy for the current MC block
!
!    E_B = \sum_{k in B, w} W_{kw} E_{kw} / \sum_{k in B, w} W_{kw}
!
!******************************************************************************** 
function afqmc_block_eh(self) result(block_eh)
  class(AfqmcEnergy), intent(in) :: self
  complex(wp)                    :: block_eh
  !local
  integer                        :: bblock, block_range(2)

  bblock = self%step__ / self%mc_des%steps_per_block
  block_range = self%mc_des%block_range(bblock)

  if (bblock == 0) then
    block_eh = self%eh(0)
  else
    associate(weight => self%sample_weight(block_range(1):block_range(2)), &
              eh => self%eh(block_range(1):block_range(2)))
      block_eh = mean(weight, eh)
    end associate
  end if
end function afqmc_block_eh


!******************************************************************************** 
!
! Estimate standard deviation of the local energy for the current MC block
!
!    sigma_B = \sum_{k in B, w} W_{kw} (E_{kw} - <E>)^2 / \sum_{k in B, w} W_{kw}
!
! Note: for stability we use: <E^2> - <E>^2 = <(E-mu)^2> - (<E>-mu)^2,
!       where mu is the variational energy of the trial wave function
!
!******************************************************************************** 
function afqmc_block_el_sd(self) result(block_el_sd)
  class(AfqmcEnergy), intent(in) :: self
  real(wp)                       :: block_el_sd
  !local
  integer                        :: bblock, block_range(2)
  complex(wp)                    :: variance, mu, block_el
  complex(wp), allocatable       :: weight(:), el_var(:)
  type(CEnergy)                  :: en

  bblock = self%step__ / self%mc_des%steps_per_block
  block_range = self%mc_des%block_range(bblock)

  if (bblock == 0) then
    variance = self%el_var(0)
  else
    weight = downsample_array(self%sample_weight(block_range(1):block_range(2)), self%mc_des%sample*self%mc_des%samplex)
    el_var = downsample_array(self%el_var(block_range(1):block_range(2)), self%mc_des%sample*self%mc_des%samplex)

    variance = mean(weight, el_var)
  end if

  mu = self%e_init%electron_energy()
  en = self%block_el()
  block_el = en%electron_energy()

  block_el_sd = sqrt(abs(variance - abs(mu-block_el)**2))
end function afqmc_block_el_sd


!******************************************************************************** 
!
! Estimate standard deviation of the hybrid energy for the current MC block
!
!    sigma_B = \sum_{k in B, w} W_{kw} (E_{kw} - <E>)^2 / \sum_{k in B, w} W_{kw}
!
! Note: for stability we used: <E^2> - <E>^2 = <(E-mu)^2> - (<E>-mu)^2,
!       where mu is the variational energy of the trial wave function
!
!******************************************************************************** 
function afqmc_block_eh_sd(self) result(block_eh_sd)
  class(AfqmcEnergy), intent(in) :: self
  real(wp)                       :: block_eh_sd
  !local
  integer                        :: bblock, block_range(2)
  complex(wp)                    :: variance, mu, block_eh

  bblock = self%step__ / self%mc_des%steps_per_block
  block_range = self%mc_des%block_range(bblock)

  if (bblock == 0) then
    variance = self%eh_var(0)
  else
    associate(weight => self%sample_weight(block_range(1):block_range(2)), &
              eh_var => self%eh_var(block_range(1):block_range(2)))
      variance = mean(weight, eh_var)
    end associate
  end if

  mu = self%e_init%electron_energy()
  block_eh = self%block_eh()

  block_eh_sd = sqrt(abs(variance - abs(mu-block_eh)**2))
end function afqmc_block_eh_sd


!******************************************************************************** 
!
! Estimate running avarage of the local energy
!
!    \sum_{kw} W_{kw} E_{kw} / \sum_{kw} W_{kw}
!
! Note: exchange contribution is calculated on the downsampled grid according to 
!       the samplex flag
!
!******************************************************************************** 
function afqmc_mean_el(self) result(mean_el)
  class(AfqmcEnergy), intent(in) :: self
  type(CEnergy)                  :: mean_el
  !local
  integer                        :: step, start_step
  complex(wp), allocatable       :: weight(:), weight_x(:)
  type(CEnergy), allocatable     :: en(:), en_x(:)
  character(len=*), parameter    :: proc_name = "afqmc_mean_el"

  if (use_profiling) call start_profiling(proc_name)

  step = self%step__
  start_step = self%start_indx()

  weight = downsample_array(self%sample_weight(start_step:step), self%mc_des%sample)
  weight_x = downsample_array(weight, self%mc_des%samplex)

  en = downsample_array(self%el(start_step:step), self%mc_des%sample)
  en_x = downsample_array(en, self%mc_des%samplex)

  if (size(en_x) == 0) then
    mean_el = self%el(0)
  else 
    mean_el = CEnergy(enuc=self%e_init%enuc, e0=self%e_init%e0)
    mean_el%e1 = mean(weight, en%e1)
    mean_el%eh = mean(weight, en%eh)
    mean_el%ex = mean(weight_x, en_x%ex)
  end if

  if (use_profiling) call end_profiling(proc_name)
end function afqmc_mean_el


!******************************************************************************** 
!
! Estimate running avarage of the hybrid energy
!
!    \sum_{kw} W_{kw} E_{kw} / \sum_{kw} W_{kw}
!
!******************************************************************************** 
function afqmc_mean_eh(self) result(mean_eh)
  class(AfqmcEnergy), intent(in) :: self
  complex(wp)                    :: mean_eh
  !local
  integer                        :: step, start_step

  step = self%step__
  start_step = self%start_indx()

  if (step < start_step) then
    mean_eh = self%eh(step)
  else
    associate (weight => self%sample_weight(start_step:step), &
               eh => self%eh(start_step:step))
      mean_eh = mean(weight, eh)
    end associate
  end if
end function afqmc_mean_eh


!******************************************************************************** 
!
! Estimate running avarage of the local energy:
!    use weighted average over walkers, but no weights for time averages
!
! Note: exchange contribution is calculated on the downsampled grid according to 
!       the samplex flag
!
!******************************************************************************** 
function afqmc_mean_el_tnw(self) result(mean_el_tnw)
  class(AfqmcEnergy), intent(in) :: self
  type(CEnergy)                  :: mean_el_tnw
  !local
  integer                        :: step, start_step
  type(CEnergy), allocatable     :: en(:), en_x(:)
  character(len=*), parameter    :: proc_name = "afqmc_mean_el"

  if (use_profiling) call start_profiling(proc_name)

  step = self%step__
  start_step = self%start_indx()

  en = downsample_array(self%el(start_step:step), self%mc_des%sample)
  en_x = downsample_array(en, self%mc_des%samplex)

  if (size(en_x) == 0) then
    mean_el_tnw = self%el(0)
  else 
    mean_el_tnw = CEnergy(enuc=self%e_init%enuc, e0=self%e_init%e0)
    mean_el_tnw%e1 = mean(en%e1)
    mean_el_tnw%eh = mean(en%eh)
    mean_el_tnw%ex = mean(en_x%ex)
  end if

  if (use_profiling) call end_profiling(proc_name)
end function afqmc_mean_el_tnw


!******************************************************************************** 
!
! Estimate running avarage of the hybrid energy:
!    use weighted average over walkers, but no weights for time averages
!
!******************************************************************************** 
function afqmc_mean_eh_tnw(self) result(mean_eh_tnw)
  class(AfqmcEnergy), intent(in) :: self
  complex(wp)                    :: mean_eh_tnw
  !local
  integer                        :: step, start_step

  step = self%step__
  start_step = self%start_indx()

  if (step < start_step) then
    mean_eh_tnw = self%eh(step)
  else
    associate (eh => self%eh(start_step:step))
      mean_eh_tnw = mean(eh)
    end associate
  end if
end function afqmc_mean_eh_tnw


!******************************************************************************** 
!
! Estimate non_weighted running avarage of the local energy
!
! Note: exchange contribution is calculated on the downsampled grid according to 
!       the samplex flag
!
!******************************************************************************** 
function afqmc_mean_el_nw(self) result(mean_el_nw)
  class(AfqmcEnergy), intent(in) :: self
  type(CEnergy)                  :: mean_el_nw
  !local
  integer                        :: step, start_step
  type(CEnergy), allocatable     :: en(:), en_x(:)
  character(len=*), parameter    :: proc_name = "afqmc_mean_el"

  if (use_profiling) call start_profiling(proc_name)

  step = self%step__
  start_step = self%start_indx()

  en = downsample_array(self%el_nw(start_step:step), self%mc_des%sample)
  en_x = downsample_array(en, self%mc_des%samplex)

  if (size(en_x) == 0) then
    mean_el_nw = self%el_nw(0)
  else 
    mean_el_nw = CEnergy(enuc=self%e_init%enuc, e0=self%e_init%e0)
    mean_el_nw%e1 = mean(en%e1)
    mean_el_nw%eh = mean(en%eh)
    mean_el_nw%ex = mean(en_x%ex)
  end if

  if (use_profiling) call end_profiling(proc_name)
end function afqmc_mean_el_nw


!******************************************************************************** 
!
! Estimate non-weighted running avarage of the hybrid energy
!
!******************************************************************************** 
function afqmc_mean_eh_nw(self) result(mean_eh_nw)
  class(AfqmcEnergy), intent(in) :: self
  complex(wp)                    :: mean_eh_nw
  !local
  integer                        :: step, start_step

  step = self%step__
  start_step = self%start_indx()

  if (step < start_step) then
    mean_eh_nw = self%eh_nw(step)
  else
    associate (eh_nw => self%eh_nw(start_step:step))
      mean_eh_nw = mean(eh_nw)
    end associate
  end if
end function afqmc_mean_eh_nw


!******************************************************************************** 
!
! Estimate standard deviation of the local energy (only for the total energy)
!
!    sigma = \sum_{kw} W_{kw} (E_{kw} - <E>)^2 / \sum_{k in B, w} W_{kw}
!
! Note: for stability we use: <E^2> - <E>^2 = <(E-mu)^2> - (<E>-mu)^2,
!       where mu is the variational energy of the trial wave function
!
!******************************************************************************** 
function afqmc_mean_el_sd(self) result(el_sd)
  class(AfqmcEnergy), intent(in) :: self
  real(wp)                       :: el_sd
  !local
  integer                        :: step, start_step
  complex(wp)                    :: variance, mu, mean_el
  complex(wp), allocatable       :: weight(:), el_var(:)
  type(CEnergy)                  :: en
  character(len=*), parameter    :: proc_name = "afqmc_mean_el_sd"

  if (use_profiling) call start_profiling(proc_name)

  step = self%step__
  start_step = self%start_indx()

  weight = downsample_array(self%sample_weight(start_step:step), self%mc_des%sample*self%mc_des%samplex)
  el_var = downsample_array(self%el_var(start_step:step), self%mc_des%sample*self%mc_des%samplex)

  if (size(el_var) == 0) then
    variance = self%el_var(step)
  else 
    variance = mean(weight, el_var)
  end if

  mu = self%e_init%electron_energy()
  en = self%mean_el()
  mean_el = en%electron_energy()

  el_sd = sqrt(abs(variance - abs(mu-mean_el)**2))

  if (use_profiling) call end_profiling(proc_name)
end function afqmc_mean_el_sd


!******************************************************************************** 
!
! Estimate standard deviation of the hybrid energy
!
!    sigma = \sum_{kw} W_{kw} (E_{kw} - <E>)^2 / \sum_{k in B, w} W_{kw}
!
! Note: for stability we use: <E^2> - <E>^2 = <(E-mu)^2> - (<E>-mu)^2,
!       where mu is the variational energy of the trial wave function
!
!******************************************************************************** 
function afqmc_mean_eh_sd(self) result(eh_sd)
  class(AfqmcEnergy), intent(in) :: self
  real(wp)                       :: eh_sd
  !local
  integer                        :: step, start_step
  complex(wp)                    :: mu, mean_eh, variance

  step = self%step__
  start_step = self%start_indx()

  if (step < start_step) then
    variance = self%eh_var(step)
  else
    associate (weight => self%sample_weight(start_step:step), &
               var => self%eh_var(start_step:step))
      variance = mean(weight, var)
    end associate
  end if

  mu = self%e_init%electron_energy()
  mean_eh = self%mean_eh()

  eh_sd = sqrt(abs(variance - abs(mu-mean_eh)**2))
end function afqmc_mean_eh_sd


!******************************************************************************** 
!
! Estimate standard deviation of the mean (SDM) of the local energy
!
! if block present:
!    use fix block size to perform block averaging
! else:
!    perform full block averaging, to estimate uncorrelated SDM
!
! Note: exchange contribution is calculated on the downsampled grid according to 
!       the samplex flag
!
!******************************************************************************** 
function afqmc_mean_el_sdm(self, block) result(el_sdm)
  class(AfqmcEnergy), intent(in) :: self
  integer, optional, intent(in)  :: block
  type(EnergyStd)                :: el_sdm
  !local
  integer                        :: step, start_step
  complex(wp), allocatable       :: weight(:), weight_x(:)
  type(CEnergy), allocatable     :: en(:), en_x(:)
  character(len=*), parameter    :: proc_name = "afqmc_mean_el_sdm"

  if (use_profiling) call start_profiling(proc_name)

  step = self%step__
  start_step = self%start_indx()

  weight = downsample_array(self%sample_weight(start_step:step), self%mc_des%sample)
  weight_x = downsample_array(weight, self%mc_des%samplex)

  en = downsample_array(self%el(start_step:step), self%mc_des%sample)
  en_x = downsample_array(en, self%mc_des%samplex)

  if (size(en_x) == 0) then
    el_sdm = EnergyStd(0.0_wp)
  else 
    if (present(block)) then
      call blocking(weight, en%e1, el_sdm%e1, el_sdm%e1_corr_len, bsize=block/self%mc_des%sample)
      call blocking(weight, en%eh, el_sdm%eh, el_sdm%eh_corr_len, bsize=block/self%mc_des%sample)
      call blocking(weight_x, en_x%ex, el_sdm%ex, el_sdm%ex_corr_len, bsize=block/self%mc_des%sample/self%mc_des%samplex)
      call blocking(weight_x, en_x%electron_energy(), el_sdm%e, el_sdm%e_corr_len, bsize=block/self%mc_des%sample/self%mc_des%samplex)
    else
      call blocking(weight, en%e1, el_sdm%e1, el_sdm%e1_corr_len)
      call blocking(weight, en%eh, el_sdm%eh, el_sdm%eh_corr_len)
      call blocking(weight_x, en_x%ex, el_sdm%ex, el_sdm%ex_corr_len)
      call blocking(weight_x, en_x%electron_energy(), el_sdm%e, el_sdm%e_corr_len)
    end if
  end if

  if (use_profiling) call end_profiling(proc_name)
end function afqmc_mean_el_sdm


!******************************************************************************** 
!
! Estimate standard deviation of the mean (SDM) of the local energy
!
! if block present:
!    use fix block size to perform block averaging
! else:
!    perform full block averaging, to estimate uncorrelated SDM
!
!******************************************************************************** 
function afqmc_mean_eh_sdm(self, block) result(eh_sdm)
  class(AfqmcEnergy), intent(in) :: self
  integer, optional, intent(in)  :: block
  type(CorrStddev)               :: eh_sdm
  !local
  integer                        :: step, start_step

  step = self%step__
  start_step = self%start_indx()

  if (step < start_step) then
    eh_sdm = CorrStddev(0.0_wp, 0.0_wp)
  else
    associate (weight => self%sample_weight(start_step:step), &
               eh => self%eh(start_step:step))

      if (present(block)) then
        call blocking(weight, eh, eh_sdm%stddev, eh_sdm%corr_len, bsize=block)
      else
        call blocking(weight, eh, eh_sdm%stddev, eh_sdm%corr_len)
      end if
    end associate
  end if
end function afqmc_mean_eh_sdm


!******************************************************************************** 
!
! Estimate standard deviation of the mean (SDM) of the non-weighted local energy
!
! if block present:
!    use fix block size to perform block averaging
! else:
!    perform full block averaging, to estimate uncorrelated SDM
!
! Note: exchange contribution is calculated on the downsampled grid according to 
!       the samplex flag
!
!******************************************************************************** 
function afqmc_mean_el_sdm_nw(self, block) result(el_sdm_nw)
  class(AfqmcEnergy), intent(in) :: self
  integer, optional, intent(in)  :: block
  type(EnergyStd)                :: el_sdm_nw
  !local
  integer                        :: step, start_step
  type(CEnergy), allocatable     :: en(:), en_x(:)
  character(len=*), parameter    :: proc_name = "afqmc_mean_el_sdm_nw"

  if (use_profiling) call start_profiling(proc_name)

  step = self%step__
  start_step = self%start_indx()

  en = downsample_array(self%el_nw(start_step:step), self%mc_des%sample)
  en_x = downsample_array(en, self%mc_des%samplex)

  if (size(en_x) == 0) then
    el_sdm_nw = EnergyStd(0.0_wp)
  else 
    if (present(block)) then
      call blocking(en%e1, el_sdm_nw%e1, el_sdm_nw%e1_corr_len, bsize=block/self%mc_des%sample)
      call blocking(en%eh, el_sdm_nw%eh, el_sdm_nw%eh_corr_len, bsize=block/self%mc_des%sample)
      call blocking(en_x%ex, el_sdm_nw%ex, el_sdm_nw%ex_corr_len, bsize=block/self%mc_des%sample/self%mc_des%samplex)
      call blocking(en_x%electron_energy(), el_sdm_nw%e, el_sdm_nw%e_corr_len, bsize=block/self%mc_des%sample/self%mc_des%samplex)
    else
      call blocking(en%e1, el_sdm_nw%e1, el_sdm_nw%e1_corr_len)
      call blocking(en%eh, el_sdm_nw%eh, el_sdm_nw%eh_corr_len)
      call blocking(en_x%ex, el_sdm_nw%ex, el_sdm_nw%ex_corr_len)
      call blocking(en_x%electron_energy(), el_sdm_nw%e, el_sdm_nw%e_corr_len)
    end if
  end if

  if (use_profiling) call end_profiling(proc_name)
end function afqmc_mean_el_sdm_nw


!******************************************************************************** 
!
! Estimate standard deviation of the mean (SDM) of the non-weighted hybrid energy
!
! if block present:
!    use fix block size to perform block averaging
! else:
!    perform full block averaging, to estimate uncorrelated SDM
!
!******************************************************************************** 
function afqmc_mean_eh_sdm_nw(self, block) result(eh_sdm_nw)
  class(AfqmcEnergy), intent(in) :: self
  integer, optional, intent(in)  :: block
  type(CorrStddev)               :: eh_sdm_nw
  !local
  integer                        :: step, start_step

  step = self%step__
  start_step = self%start_indx()

  if (step < start_step) then
    eh_sdm_nw = CorrStddev(0.0_wp, 0.0_wp)
  else
    associate (eh_nw => self%eh_nw(start_step:step))
      if (present(block)) then
        call blocking(eh_nw, eh_sdm_nw%stddev, eh_sdm_nw%corr_len, bsize=block)
      else
        call blocking(eh_nw, eh_sdm_nw%stddev, eh_sdm_nw%corr_len)
      end if
    end associate
  end if
end function afqmc_mean_eh_sdm_nw


!******************************************************************************** 
!
! Estimate standard deviation of the mean (SDM) of the local energy:
!    use weighted average over walkers, but no weights for time averages
!
! if block present:
!    use fix block size to perform block averaging
! else:
!    perform full block averaging, to estimate uncorrelated SDM
!
! Note: exchange contribution is calculated on the downsampled grid according to 
!       the samplex flag
!
!******************************************************************************** 
function afqmc_mean_el_sdm_tnw(self, block) result(el_sdm_tnw)
  class(AfqmcEnergy), intent(in) :: self
  integer, optional, intent(in)  :: block
  type(EnergyStd)                :: el_sdm_tnw
  !local
  integer                        :: step, start_step
  type(CEnergy), allocatable     :: en(:), en_x(:)
  character(len=*), parameter    :: proc_name = "afqmc_mean_el_sdm_tnw"

  if (use_profiling) call start_profiling(proc_name)

  step = self%step__
  start_step = self%start_indx()

  en = downsample_array(self%el(start_step:step), self%mc_des%sample)
  en_x = downsample_array(en, self%mc_des%samplex)

  if (size(en_x) == 0) then
    el_sdm_tnw = EnergyStd(0.0_wp)
  else 
    if (present(block)) then
      call blocking(en%e1, el_sdm_tnw%e1, el_sdm_tnw%e1_corr_len, bsize=block/self%mc_des%sample)
      call blocking(en%eh, el_sdm_tnw%eh, el_sdm_tnw%eh_corr_len, bsize=block/self%mc_des%sample)
      call blocking(en_x%ex, el_sdm_tnw%ex, el_sdm_tnw%ex_corr_len, bsize=block/self%mc_des%sample/self%mc_des%samplex)
      call blocking(en_x%electron_energy(), el_sdm_tnw%e, el_sdm_tnw%e_corr_len, bsize=block/self%mc_des%sample/self%mc_des%samplex)
    else
      call blocking(en%e1, el_sdm_tnw%e1, el_sdm_tnw%e1_corr_len)
      call blocking(en%eh, el_sdm_tnw%eh, el_sdm_tnw%eh_corr_len)
      call blocking(en_x%ex, el_sdm_tnw%ex, el_sdm_tnw%ex_corr_len)
      call blocking(en_x%electron_energy(), el_sdm_tnw%e, el_sdm_tnw%e_corr_len)
    end if
  end if

  if (use_profiling) call end_profiling(proc_name)
end function afqmc_mean_el_sdm_tnw


!******************************************************************************** 
!
! Estimate standard deviation of the mean (SDM) of the hybrid energy:
!    use weighted average over walkers, but no weights for time averages
!
! if block present:
!    use fix block size to perform block averaging
! else:
!    perform full block averaging, to estimate uncorrelated SDM
!
!******************************************************************************** 
function afqmc_mean_eh_sdm_tnw(self, block) result(eh_sdm_tnw)
  class(AfqmcEnergy), intent(in) :: self
  integer, optional, intent(in)  :: block
  type(CorrStddev)               :: eh_sdm_tnw
  !local
  integer                        :: step, start_step

  step = self%step__
  start_step = self%start_indx()

  if (step < start_step) then
    eh_sdm_tnw = CorrStddev(0.0_wp, 0.0_wp)
  else
    associate (eh => self%eh(start_step:step))
      if (present(block)) then
        call blocking(eh, eh_sdm_tnw%stddev, eh_sdm_tnw%corr_len, bsize=block)
      else
        call blocking(eh, eh_sdm_tnw%stddev, eh_sdm_tnw%corr_len)
      end if
    end associate
  end if
end function afqmc_mean_eh_sdm_tnw


!******************************************************************************** 
!
! Calculate current mean energy using cumulative sums of w and w*E
!
!******************************************************************************** 
function afqmc_fast_mean_el(self) result(mean_el)
  class(AfqmcEnergy), intent(in) :: self
  type(CEnergy)                  :: mean_el

  if (abs(self%sample_weight__) < 0.1_wp) then
    mean_el = self%e_init
  else 
    mean_el = self%el__ / self%sample_weight__
  end if
end function afqmc_fast_mean_el


!******************************************************************************** 
!
! Calculate current mean hybrid energy using cumulative sums of w and w*E
!
!******************************************************************************** 
function afqmc_fast_mean_eh(self) result(mean_eh)
  class(AfqmcEnergy), intent(in) :: self
  complex(wp)                    :: mean_eh

  if (abs(self%sample_weight__) < 0.1_wp) then
    mean_eh = self%e_init%electron_energy()
  else 
    mean_eh = self%eh__ / self%sample_weight__
  end if
end function afqmc_fast_mean_eh


!******************************************************************************** 
!
! Estimate block average sign for a given MC block
!
!    s_B =  |\sum_{k in B} W_k| / \sum_{k in B} |W_k|
!
!    W_k = \sum_w W_{kw} e^{i \theta_{kw}}
!
!******************************************************************************** 
function afqmc_block_sign(self) result(block_sign)
  class(AfqmcEnergy), intent(in) :: self
  real(wp)                       :: block_sign
  !local
  integer                        :: bblock, block_range(2)

  bblock = self%step__ / self%mc_des%steps_per_block
  block_range = self%mc_des%block_range(bblock)

  if (bblock == 0) then
    block_sign = 1.0_wp
  else
    associate(weight => self%weight(block_range(1):block_range(2)), &
              sample_weight => self%sample_weight(block_range(1):block_range(2)))

    block_sign = abs(sum(sample_weight)) / sum(weight)
    end associate
  end if
end function afqmc_block_sign


!******************************************************************************** 
!
! Estimate average sign in the sampling phase
!
!    s =  |\sum_{k} W_k| / \sum_{k} |W_k|
!
!    W_k = \sum_w W_{kw} e^{i \theta_{kw}}
!
!******************************************************************************** 
function afqmc_mean_sign(self) result(mean_sign)
  class(AfqmcEnergy), intent(in) :: self
  real(wp)                       :: mean_sign
  !local
  integer                        :: step, start_step

  step = self%step__
  start_step = self%start_indx()

  if (step < start_step) then
    mean_sign = 1.0_wp
  else
    associate (weight => self%weight(start_step:step), &
               sample_weight => self%sample_weight(start_step:step))

    mean_sign = abs(sum(sample_weight)) / sum(weight)
    end associate
  end if
end function afqmc_mean_sign


!******************************************************************************** 
!
! Read AfqmcEnergy object from the file (for afqmc_chkpt.f)
!
!******************************************************************************** 
subroutine read_afqmc_energy_from_file(self, fh)
  class(AfqmcEnergy), intent(inout) :: self
  type(FileHandle), intent(in)      :: fh
  !local
  integer                           :: i, step

  read(fh%funit) step
  self%step__ = step

  !adjust McDescriptor
  self%mc_des%step = step
  self%mc_des%block = step / self%mc_des%steps_per_block
  !if the total number of blocks is smaller than the current block,
  !adjust nblocks to aviod allocation errors
  if (self%mc_des%block > self%mc_des%eqblocks+self%mc_des%nblocks) then
    self%mc_des%nblocks = self%mc_des%block - self%mc_des%eqblocks
  end if

  !note: af_en arrays reallocate since mc_des changes
  call self%alloc_arrays()

  read(fh%funit) self%weight(0:step)
  read(fh%funit) self%sample_weight(0:step)
  read(fh%funit) self%eh(0:step)
  read(fh%funit) self%eh_nw(0:step)
  read(fh%funit) self%el_var(0:step)
  read(fh%funit) self%eh_var(0:step)
  call read_c_energies(self%el(0:step), fh%funit)
  call read_c_energies(self%el_nw(0:step), fh%funit)

  call self%setup_cumulatives()
end subroutine read_afqmc_energy_from_file


!******************************************************************************** 
!
! Write AfqmcEnergy object to the file (for afqmc_chkpt.f)
!
!******************************************************************************** 
subroutine write_afqmc_energy_to_file(self, fh)
  class(AfqmcEnergy), intent(in) :: self
  type(FileHandle), intent(in)   :: fh
  !local 
  integer                        :: step

  step = self%step__

  write(fh%funit) step
  write(fh%funit) self%weight(0:step)
  write(fh%funit) self%sample_weight(0:step)
  write(fh%funit) self%eh(0:step)
  write(fh%funit) self%eh_nw(0:step)
  write(fh%funit) self%el_var(0:step)
  write(fh%funit) self%eh_var(0:step)
  call write_c_energies(self%el(0:step), fh%funit)
  call write_c_energies(self%el_nw(0:step), fh%funit)
end subroutine write_afqmc_energy_to_file
                 
                 
!******************************************************************************** 
!
! MPI Bcast AfqmcEnergy object
!
!******************************************************************************** 
subroutine mpi_bcast_afqmc_energy(self, comm, root)
  class(AfqmcEnergy), intent(inout)   :: self
  type(mpi_communicator), intent(in)  :: comm
  integer, intent(in)                 :: root

  if (comm%mpirank /= root) call self%alloc_arrays()

  call comm%bcast(self%step__, root)
  call comm%bcast(self%weight, root)
  call comm%bcast(self%sample_weight, root)
  call comm%bcast(self%eh, root)
  call comm%bcast(self%eh_nw, root)
  call comm%bcast(self%el_var, root)
  call comm%bcast(self%eh_var, root)
  call comm%bcast(self%el%enuc, root)
  call comm%bcast(self%el%e0, root)
  call comm%bcast(self%el%e1, root)
  call comm%bcast(self%el%eh, root)
  call comm%bcast(self%el%ex, root)
  call comm%bcast(self%el_nw%enuc, root)
  call comm%bcast(self%el_nw%e0, root)
  call comm%bcast(self%el_nw%e1, root)
  call comm%bcast(self%el_nw%eh, root)
  call comm%bcast(self%el_nw%ex, root)

  call self%setup_cumulatives()
end subroutine mpi_bcast_afqmc_energy


!******************************************************************************** 
!                
! Downsample array to the coarser exchange grid
!
!******************************************************************************** 
function downsample_array_cmplx(array, subsample) result(array_sub)
  complex(wp), intent(in)  :: array(:)
  integer, intent(in)      :: subsample
  complex(wp), allocatable :: array_sub(:)
  !local
  integer                  :: i, n, ii

  n = size(array) / subsample
  
  if (allocated(array_sub)) deallocate(array_sub)
  allocate(array_sub(n))

  do i = 1, n
    ii = i * subsample
    array_sub(i) = array(ii)
  end do
end function downsample_array_cmplx


!******************************************************************************** 
!
! Downsample CEnergy array to the coarser exchange grid
!
!******************************************************************************** 
function downsample_array_energy(array, subsample) result(array_sub)
  type(CEnergy), intent(in)  :: array(:)
  integer, intent(in)        :: subsample
  type(CEnergy), allocatable :: array_sub(:)
  !local
  integer                    :: i, n, ii

  n = size(array) / subsample
  
  if (allocated(array_sub)) deallocate(array_sub)
  allocate(array_sub(n))

  do i = 1, n
    ii = i * subsample
    array_sub(i) = array(ii)
  end do
end function downsample_array_energy

end module afqmc_energy