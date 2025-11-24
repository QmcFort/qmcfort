! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module afqmc_walker_factories

#include "../preproc.inc"
use constants
use mpi
use profiling
use afqmc_walker
use wave_trial, only: WaveTrial
use wave_trial_sd, only: WaveTrialSD
use wave_trial_cas, only: WaveTrialCAS
use wave_trial_noci, only: WaveTrialNOCI
use noci_selector_des, only: NOCISelectorDes
use noci_selector, only: NOCISelector

use standalone, only: dump_array

implicit none

private
public :: afqmc_walker_factory_coeff, afqmc_walker_factory_trial

contains


!******************************************************************************** 
!
! AfqmcWalker Factory
!
!    all walkers initialized according to the given orbitals
!
!******************************************************************************** 
subroutine afqmc_walker_factory_coeff(phi_t, coeff, walkers, nwalkers, sp_proj)
  class(WaveTrial), intent(inout)             :: phi_t
  complex(wp), intent(in)                     :: coeff(:,:,:)
  type(AfqmcWalker), allocatable, intent(out) :: walkers(:)
  integer, intent(in)                         :: nwalkers
  logical, optional, intent(in)               :: sp_proj
  !local
  integer                                     :: w
  type(AfqmcWalker)                           :: walker
  character(len=*), parameter                 :: proc_name = "afqmc_walker_factory_coeff"

  if (use_profiling) call start_profiling(proc_name)

  allocate(walkers(nwalkers))

  walker = AfqmcWalker(phi_t, coeff, sp_proj)

  do w = 1, nwalkers
    walkers(w) =  walker
  end do

  if (use_profiling) call end_profiling(proc_name)
end subroutine afqmc_walker_factory_coeff


!******************************************************************************** 
!
! AfqmcWalker Factory
!
!    walkers initialized according to the trial wave function
!
!******************************************************************************** 
subroutine afqmc_walker_factory_trial(phi_t, walkers, nwalkers, sp_proj)
  class(WaveTrial), intent(inout)             :: phi_t
  type(AfqmcWalker), allocatable, intent(out) :: walkers(:)
  integer, intent(in)                         :: nwalkers
  logical, optional, intent(in)               :: sp_proj
  !local
  integer                                     :: i, w, ww
  integer, allocatable                        :: ncopies(:)
  complex(wp), allocatable                    :: coeff(:,:,:,:)
  type(AfqmcWalker), allocatable              :: walkers_det(:)

  call get_nwalkers_per_det(phi_t, coeff, ncopies, nwalkers)

  allocate(walkers_det(size(ncopies)))
  do i = 1, size(ncopies)
    if (ncopies(i) == 0) cycle
    walkers_det(i) = AfqmcWalker(phi_t, coeff(:,:,:,i), sp_proj)
    walkers_det(i)%id = i
  end do

  allocate(walkers(nwalkers))

  ww = 0
  do i = 1, size(ncopies)
    do w = 1, ncopies(i)
      ww = ww + 1
      walkers(ww) = walkers_det(i)
    end do
  end do
end subroutine afqmc_walker_factory_trial


!******************************************************************************** 
!
! Extract orbitals from determinants in WaveTrial
! and number of walker copies per determinant
!
!******************************************************************************** 
subroutine get_nwalkers_per_det(phi_t, coeff, ncopies, nwalkers)
  class(WaveTrial), intent(in)          :: phi_t
  complex(wp), allocatable, intent(out) :: coeff(:,:,:,:)
  integer, allocatable, intent(out)     :: ncopies(:)
  integer, intent(in)                   :: nwalkers
  !local
  integer                               :: diff
  real(wp), allocatable                 :: weights(:)
  type(NOCISelectorDes)                 :: noci_des
  type(NOCISelector)                    :: noci

  select type (phi_t)
    type is (WaveTrialSD)
      allocate(ncopies(1))
      ncopies = nwalkers

      allocate(coeff(size(phi_t%coeff,1), size(phi_t%coeff,2), size(phi_t%coeff,3), 1))
      coeff(:,:,:,1) = phi_t%coeff
    type is (WaveTrialNOCI)
      noci_des = NOCISelectorDes(1000, 0.1_wp, 0.00001_wp)
      noci = NOCISelector(noci_des, phi_t)

      allocate(coeff, source=phi_t%coeff)

      weights = noci%get_weights()
      allocate(ncopies(size(weights)))
      ncopies = nint(weights * nwalkers)
      diff = nwalkers - sum(ncopies)
      ncopies(1) = ncopies(1) + diff
    type is (WaveTrialCAS)
      if (comm_world%mpirank == 0) write(*,*) "can not create walker ensemble from WaveTrialCAS type"
      stop
    class default
      if (comm_world%mpirank == 0) write(*,*) "Unknown WaveTrial type - can not create walker ensemble from it"
      stop
  end select
end subroutine get_nwalkers_per_det

end module afqmc_walker_factories