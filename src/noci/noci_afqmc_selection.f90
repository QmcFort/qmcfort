! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module noci_afqmc_selection

#include "../preproc.inc"
!$ use omp_lib
use constants
use mpi
use profiling
use qmcfort_io
use qmcfort_pos
use standalone, only: concatenate
use method_base
use hamilton_vars, only: hamil_descriptor
use hamilton_type
use wave_trial_noci, only: WaveTrialNOCI
use noci_intags
use noci_selector_des, only: NOCISelectorDes
use noci_selector, only: NOCISelector
use wave_trial_factories, only: wave_trial_factory_afqmc, wave_trial_factory_noci
use afqmc_walker
use afqmc_walker_factories, only: afqmc_walker_factory_coeff, afqmc_walker_factory_trial
use afqmc, only: setup_afqmc_hamiltonian, &
                 calculate_observables, reorthogonalize_walkers, afqmc_propagator

implicit none

private
public :: do_noci_afqmc_selection

type NociAfqmcDes
  integer  :: nepochs
  integer  :: steps_per_epoch
  integer  :: select_from_step
  real(wp) :: energy_thresh_max
  real(wp) :: energy_thresh_min
  real(wp) :: sigma_eloc
  logical  :: update_afqmc_trial

  real(wp) :: delta_e
end type NociAfqmcDes

type NociAfqmcManager
  type(NociAfqmcDes)     :: des
  type(afqmc_propagator) :: afqmc_mng
  type(NOCISelector)     :: noci_selector
contains
  procedure :: report_setup => report_setup_noci_afqmc
end type NociAfqmcManager

interface NociAfqmcDes
  procedure :: finit_noci_afqmc_des
end interface NociAfqmcDes

contains 

!******************************************************************************** 
!
! Initialization of the NociAfqmcDes object
!
!******************************************************************************** 
subroutine init_noci_afqmc_des(nepochs, steps_per_epoch, select_from_step, energy_thresh_max, energy_thresh_min, sigma_eloc, update_afqmc_trial, self)
  integer, intent(in)             :: nepochs
  integer, intent(in)             :: steps_per_epoch
  integer, intent(in)             :: select_from_step
  real(wp), intent(in)            :: energy_thresh_max
  real(wp), intent(in)            :: energy_thresh_min
  real(wp), intent(in)            :: sigma_eloc
  logical, intent(in)             :: update_afqmc_trial
  type(NociAfqmcDes), intent(out) :: self

  self%nepochs = nepochs
  self%steps_per_epoch = steps_per_epoch
  self%select_from_step = select_from_step
  self%energy_thresh_max = energy_thresh_max
  self%energy_thresh_min = energy_thresh_min
  self%sigma_eloc = sigma_eloc
  self%update_afqmc_trial = update_afqmc_trial

  if (self%nepochs == 1) then
    self%delta_e = 1.0_wp
  else 
    self%delta_e = (energy_thresh_min/energy_thresh_max)**(1.0_wp/(self%nepochs-1))
  end if
end subroutine init_noci_afqmc_des


!******************************************************************************** 
!
! NociAfqmcDes constructor
!
!******************************************************************************** 
function finit_noci_afqmc_des(nepochs, steps_per_epoch, select_from_step, energy_thresh_max, energy_thresh_min, sigma_eloc, update_afqmc_trial) result(self)
  integer, intent(in)  :: nepochs
  integer, intent(in)  :: steps_per_epoch
  integer, intent(in)  :: select_from_step
  real(wp), intent(in) :: energy_thresh_max
  real(wp), intent(in) :: energy_thresh_min
  real(wp), intent(in) :: sigma_eloc
  logical, intent(in)  :: update_afqmc_trial
  type(NociAfqmcDes)   :: self
  
  call init_noci_afqmc_des(nepochs, steps_per_epoch, select_from_step, energy_thresh_max, energy_thresh_min, sigma_eloc, update_afqmc_trial, self)
end function finit_noci_afqmc_des


!******************************************************************************** 
!
! Main routine in NOCI AFQMC Selection module
!
!******************************************************************************** 
subroutine do_noci_afqmc_selection(ham, struc)
  type(Hamilton), target, intent(inout) :: ham
  type(structure), intent(in)           :: struc
  !local
  integer                               :: epoch, step, nwalkers, walkers_per_rank
  real(wp)                              :: tau
  complex(wp), allocatable              :: coeff(:,:,:)
  type(NociAfqmcManager)                :: mng
  type(AfqmcWalker), allocatable        :: walkers(:), sel_walkers(:)
  type(QmcFortMethod)                   :: method

  method = QmcFortMethod(method="noci_afqmc", def_active=.false., basis="mo ao_orth", integral_mode="cholesky")
  if (.not. method%is_active()) return
  call ham%meet_method_req(method)

  call setup_noci_afqmc_manager(mng, ham, struc)
  call mng%report_setup()

  allocate(coeff(size(ham%phi,1), size(ham%phi,2), size(ham%phi,3)))
  coeff = ham%phi

  tau = mng%afqmc_mng%af_des%mc_des%tau
  nwalkers = mng%afqmc_mng%af_des%mc_des%nwalkers
  walkers_per_rank = mng%afqmc_mng%af_des%mc_des%walkers_per_rank

  do epoch = 1, mng%des%nepochs
    !update AFQMC trial wave function to enhance the random walk
    if (mng%des%update_afqmc_trial .and. epoch>1) then
      call wave_trial_factory_afqmc(mng%afqmc_mng%tags%wave_trial_request, mng%afqmc_mng%wtdes, ham, mng%afqmc_mng%phi_t)
      call mng%afqmc_mng%phi_t%setup(lself=.true.)
    end if

    call set_noci_energy_thresh(mng%noci_selector%des, mng%des, epoch)
    !call afqmc_walker_factory_coeff(mng%afqmc_mng%phi_t, coeff, walkers, walkers_per_rank)
    call afqmc_walker_factory_trial(mng%afqmc_mng%phi_t, walkers, walkers_per_rank)

    call report_noci_afqmc_epoch(mng%noci_selector%des, epoch)

    do step = 1, mng%des%steps_per_epoch
      call mng%afqmc_mng%af_prop%propagate_walkers(mng%afqmc_mng%afgen, mng%afqmc_mng%imp_samp, walkers)
      call calculate_observables(mng%afqmc_mng, walkers)
      call reorthogonalize_walkers(mng%afqmc_mng, walkers)
      if (step <= mng%des%select_from_step) cycle
      call preselect_walkers(mng, walkers, sel_walkers)
      call report_noci_afqmc_step(step, tau, size(sel_walkers), nwalkers)
      call noci_walker_selection(mng%noci_selector, sel_walkers)
    end do

    if (mng%des%update_afqmc_trial .and. comm_world%mpirank==0) call mng%noci_selector%phi%write_to_file()
    call comm_world%barrier()
  end do

!debug: nocq
  !call mng%noci_selector%compress(0.001_wp)

  call finalize_noci_afqmc_selection(mng)
end subroutine do_noci_afqmc_selection


!******************************************************************************** 
!
! Setup manager variable for the NOCI AFQC Selection
!
!******************************************************************************** 
subroutine setup_noci_afqmc_manager(mng, ham, struc)
  type(NociAfqmcManager), intent(out)         :: mng
  type(Hamilton), target, intent(inout)       :: ham
  type(structure), intent(in)                 :: struc
  !local
  type(NociTags)                              :: noci_tags
  type(NOCISelectorDes)                       :: noci_des
  type(WaveTrialNOCI)                         :: phi_t_noci

  noci_tags = NociTags()

  mng%des = NociAfqmcDes(noci_tags%nepochs, noci_tags%steps_per_epoch, noci_tags%select_from_step, noci_tags%energy_thresh_max, &
                         noci_tags%energy_thresh_min, noci_tags%sigma_eloc, noci_tags%update_afqmc_trial)

  call setup_afqmc_hamiltonian(mng%afqmc_mng, ham, struc)

  !create NOCISelectroDes for NOCI selection
  noci_des = NOCISelectorDes(noci_tags%ndet_max, noci_tags%ovlp_thresh, noci_tags%energy_thresh)

  !NOCISelector initialization
  call wave_trial_factory_noci(mng%afqmc_mng%phi_t%des, ham, phi_t_noci)
  call phi_t_noci%setup(lself=.false.)
  mng%noci_selector = NOCISelector(noci_des, phi_t_noci)
end subroutine setup_noci_afqmc_manager


subroutine set_noci_energy_thresh(noci_selector_des, noci_afqmc_des, epoch)
  type(NOCISelectorDes), intent(inout) :: noci_selector_des
  type(NociAfqmcDes), intent(in)       :: noci_afqmc_des
  integer, intent(in)                  :: epoch
  !local
  real(wp)                             :: emax, delta_e

  emax = noci_afqmc_des%energy_thresh_max
  delta_e = noci_afqmc_des%delta_e

  noci_selector_des%energy_thresh = emax * delta_e**(epoch-1)
end subroutine set_noci_energy_thresh


!******************************************************************************** 
!
! Pre-Select walkers for the NOCI selection
!
!    Walkers with E_L < \bar E_L - n sigma[E_L] are accepted
!
!******************************************************************************** 
subroutine preselect_walkers(mng, walkers, sel_walkers)
  type(NociAfqmcManager), intent(in)          :: mng
  type(AfqmcWalker), intent(in)               :: walkers(:)
  type(AfqmcWalker), allocatable, intent(out) :: sel_walkers(:)
  !local
  integer                                     :: w, nw, nw_sel, indx
  integer, allocatable                        :: sel_indices(:)
  real(wp)                                    :: mean_locen, std_locen
  real(wp), allocatable                       :: locen(:)
  type(AfqmcWalker), allocatable              :: temp_walkers(:)
  character(len=*), parameter                 :: proc_name = "preselect_walkers"

  if (use_profiling) call start_profiling(proc_name)

  nw = size(walkers)
  allocate(locen(nw))
  locen = real(walkers%energy%total_energy(), wp)

  mean_locen = mean(locen, comm=comm_world)
  std_locen = std(locen, comm=comm_world)

  nw_sel = count(locen<=mean_locen-mng%des%sigma_eloc*std_locen)
  allocate(sel_indices(nw_sel))
  sel_indices = pack([(w, w=1,nw)], locen<=mean_locen-mng%des%sigma_eloc*std_locen)

  allocate(temp_walkers(nw_sel))
  do w = 1, nw_sel
    indx = sel_indices(w)
    temp_walkers(w) = walkers(indx)
  end do

  call gather_walkers(comm_world, mng%afqmc_mng%ham%des, temp_walkers, sel_walkers, root=0) 

  if (use_profiling) call end_profiling(proc_name)
end subroutine preselect_walkers


!******************************************************************************** 
!
! Test AFQMC walkers for the NOCI selection
!
!******************************************************************************** 
subroutine noci_walker_selection(noci_selector, walkers)
  type(NOCISelector), intent(inout) :: noci_selector
  type(AfqmcWalker), intent(in)     :: walkers(:)
  !local
  integer                           :: w
  complex(wp), allocatable          :: coeff(:,:,:)
  character(len=*), parameter       :: proc_name = "noci_walker_selection"

  if (size(noci_selector%phi%ci_coeff) >= noci_selector%des%ndet_max) return  

  if (use_profiling) call start_profiling(proc_name)
  
  if (comm_world%mpirank == 0) then
    do w = 1, size(walkers)
      coeff = prepare_coeffs(walkers(w)%coeff, noci_selector%phi%ham%des)
      call noci_selector%update(coeff)
    end do
  end if

  if (use_profiling) call end_profiling(proc_name)
contains
  !small subroutine to transform 2d coeff into 3d coeff
  function prepare_coeffs(coeff_in, hdes) result(coeff_out)
    complex(wp), intent(in)            :: coeff_in(:,:)
    type(hamil_descriptor), intent(in) :: hdes
    complex(wp), allocatable           :: coeff_out(:,:,:)
    !local
    integer                            :: spin, nocc, nel, w1, w2
    character(len=*), parameter        :: proc_name = "prepare_coeffs"

    if (profile_code) call start_profiling(proc_name)

    nocc = maxval(hdes%nel)
    allocate(coeff_out(hdes%n, nocc, hdes%ispin))
    coeff_out = zeroc

    do spin = 1, hdes%ispin
      nel = hdes%nel(spin)
      w1 = hdes%nell(spin) + 1
      w2 = hdes%nell(spin) + nel
      coeff_out(:,1:nel,spin) = coeff_in(:,w1:w2)
    end do

    if (profile_code) call end_profiling(proc_name)
  end function prepare_coeffs
end subroutine noci_walker_selection


!******************************************************************************** 
!
! Report and write out about NOCI AFQMC Selection
!
!******************************************************************************** 
subroutine finalize_noci_afqmc_selection(mng)
  type(NociAfqmcManager), intent(in) :: mng

  !dump out selected determinants
  if (comm_world%mpirank==0) call mng%noci_selector%phi%write_to_file()

  !write out collected determinants
  call mng%noci_selector%report(io%screen)
  call mng%noci_selector%report(io%qmcfort_log)
  call mng%noci_selector%report(io%qmcfort_out)
end subroutine finalize_noci_afqmc_selection


!******************************************************************************** 
!
! Report NOCI AFQMC Selection setup
!
!******************************************************************************** 
subroutine report_setup_noci_afqmc(self)
  class(NociAfqmcManager), intent(inout) :: self
  !local
  character(len=100)                     :: intfor, logfor, realfor, expfor, charfor
  
  intfor  = '(1x,t5,a,t50,"= ",10i8)'
  logfor  = '(1x,t5,a,t50,"= ",10l6)'
  realfor = '(1x,t5,a,t50,"= ",10f12.6)'
  expfor = '(1x,t5,a,t50,"= ",10es12.4)'
  charfor = '(1x,t5,a,t50,"= ",a)'

  if (comm_world%mpirank == 0) then
    write(io%qmcfort_out%funit,*)
    write(io%qmcfort_out%funit,*)      "NOCI-AFQMC Setup"
    write(io%qmcfort_out%funit,*)      starshort
    write(io%qmcfort_out%funit,*)      "NOCI selection descriptor:"
    write(io%qmcfort_out%funit,intfor) "number of epochs",                 self%des%nepochs
    write(io%qmcfort_out%funit,intfor) "steps per epoch",                  self%des%steps_per_epoch
    write(io%qmcfort_out%funit,realfor)"local energy preselection width",  self%des%sigma_eloc
    write(io%qmcfort_out%funit,expfor) "energy threshold max",             self%des%energy_thresh_max
    write(io%qmcfort_out%funit,expfor) "energy threshold min",             self%des%energy_thresh_min
    write(io%qmcfort_out%funit,expfor) "overlap threshold",                self%noci_selector%des%ovlp_thresh
    write(io%qmcfort_out%funit,intfor) "maximal number of determinants",   self%noci_selector%des%ndet_max
    write(io%qmcfort_out%funit,*)
    write(io%qmcfort_out%funit,*)      "AFQMC descriptor:"
    write(io%qmcfort_out%funit,realfor)"timestep",                         self%afqmc_mng%af_des%mc_des%tau
    write(io%qmcfort_out%funit,intfor) "number of afqmc walkers",          self%afqmc_mng%af_des%mc_des%nwalkers
    write(io%qmcfort_out%funit,intfor) "number of afqmc walkers per proc", self%afqmc_mng%af_des%mc_des%walkers_per_rank
  end if
    
  call self%afqmc_mng%imp_samp%report_setup(comm_world, io%qmcfort_out)
  if (comm_world%mpirank == 0) call self%afqmc_mng%af_proj%report(io%qmcfort_out%funit)
  call self%afqmc_mng%afgen%report(comm_world, io%qmcfort_out)
  call self%afqmc_mng%af_prop%report(comm_world, io%qmcfort_out%funit)
  call self%afqmc_mng%phi_t%report(comm_world, io%qmcfort_out)
  call self%afqmc_mng%phi_t%report(comm_world, io%screen)

  if (comm_world%mpirank == 0) then
    write(io%qmcfort_out%funit,*)
    write(io%qmcfort_out%funit,*)      "NOCI-AFQMC"
    write(io%qmcfort_out%funit,*)      starshort

    write(io%qmcfort_log%funit,*)
    write(io%qmcfort_log%funit,*)      "NOCI-AFQMC"
    write(io%qmcfort_log%funit,*)      starshort

    write(io%screen%funit,*)
    write(io%screen%funit,*)      "NOCI-AFQMC"
    write(io%screen%funit,*)      starshort
  end if
end subroutine report_setup_noci_afqmc


subroutine report_noci_afqmc_epoch(noci_des, epoch)
  type(NOCISelectorDes), intent(in) :: noci_des
  integer, intent(in)               :: epoch

  if (comm_world%mpirank == 0) then
    write(io%qmcfort_out%funit,100) "  Epoch ", epoch, " :  energy_thresh = ", noci_des%energy_thresh, ",  ovlp_thresh = ", noci_des%ovlp_thresh
    write(io%qmcfort_log%funit,100) "  Epoch ", epoch, " :  energy_thresh = ", noci_des%energy_thresh, ",  ovlp_thresh = ", noci_des%ovlp_thresh
    write(io%screen%funit,100)      "  Epoch ", epoch, " :  energy_thresh = ", noci_des%energy_thresh, ",  ovlp_thresh = ", noci_des%ovlp_thresh
  end if
  100 format (1x,a,i3,a,es10.2,a,es10.2)
end subroutine report_noci_afqmc_epoch

subroutine report_noci_afqmc_step(step, tau, sel_walkers, tot_walkers)
  integer, intent(in)  :: step
  real(wp), intent(in) :: tau
  integer, intent(in)  :: sel_walkers
  integer, intent(in)  :: tot_walkers
  !local
  real(wp)             :: ratio

  if (comm_world%mpirank == 0) then
    ratio = real(sel_walkers, wp) / tot_walkers

    write(io%qmcfort_out%funit,100) "    Step ", step, " :  tau = ", step*tau, ",  Nw for selection = ", sel_walkers, ",  ratio = ", ratio
    write(io%qmcfort_log%funit,100) "    Step ", step, " :  tau = ", step*tau, ",  Nw for selection = ", sel_walkers, ",  ratio = ", ratio
    write(io%screen%funit,100)      "    Step ", step, " :  tau = ", step*tau, ",  Nw for selection = ", sel_walkers, ",  ratio = ", ratio

    write(io%qmcfort_out%funit,101) " NOCI size ", " old NOCI energy ", " new NOCI energy ", " energy residual ", " ovlp residual "
    write(io%qmcfort_out%funit,101) " ========= ", " =============== ", " =============== ", " =============== ", " ============= "
    write(io%qmcfort_log%funit,101) " NOCI size ", " old NOCI energy ", " new NOCI energy ", " energy residual ", " ovlp residual "
    write(io%qmcfort_log%funit,101) " ========= ", " =============== ", " =============== ", " =============== ", " ============= "
    write(io%screen%funit,101)      " NOCI size ", " old NOCI energy ", " new NOCI energy ", " energy residual ", " ovlp residual "
    write(io%screen%funit,101)      " ========= ", " =============== ", " =============== ", " =============== ", " ============= "
  end if
  100 format (1x,a,i4,a,f8.4,a,i6,a,f8.4)
  101 format (1x,t18,a,t32,a,t51,a,t70,a,t87,a)
end subroutine report_noci_afqmc_step
end module noci_afqmc_selection