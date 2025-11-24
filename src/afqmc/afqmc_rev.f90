! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module afqmc_rev

#include "../preproc.inc"

use afqmc_energy
use afqmc_walker
use constants
use energy_types
use profiling
use qmcfort_io

use afqmc_descriptor, only: AfqmcDescriptor
use afqmc_proj, only: AfqmcProj
use file_handle, only: FileHandle
use growth_est, only: GrowthEst
use hamilton_vars, only: hamil_descriptor
use mc_descriptor, only: McDescriptor
use mpi, only: mpi_communicator, comm_world
use standalone, only: cap_value, modul

implicit none 

type :: rare_events
  logical                :: active, el_capping, eh_capping, r_capping, t_capping, i_capping, w_capping, p_capping
  integer                :: nsteps, nsteps_eq, sample, nwalkers, ng
  integer                :: number_elu, number_ell, number_eluu, number_ellu
  integer                :: number_ehu, number_ehl, number_ehuu, number_ehlu
  integer                :: number_ru, number_rl, number_ruu, number_rlu
  integer                :: number_tu, number_tl, number_tuu, number_tlu
  integer                :: number_iu, number_il, number_iuu, number_ilu
  integer                :: number_p, number_pu
  integer                :: number_w, number_wu
  integer                :: number_fb
  integer                :: number_acc
  integer                :: rev_lifemax, ndeath_rev, ndeath_revu
  integer                :: pers_lifemax, ndeath_pers, ndeath_persu
  integer                :: el_capping_mode
  real(wp)               :: ediff_l, ediff_h
  real(wp)               :: width_el, width_eh
  real(wp)               :: width_r, width_t
  real(wp)               :: el_ratio, eh_ratio
  real(wp)               :: r_ratio, t_ratio
  real(wp)               :: maxweight, impwmax, impwmin
  real(wp)               :: mean_el, mean_el2, mean_el_eq, mean_el2_eq
  real(wp)               :: mean_eh, mean_eh2, mean_eh_eq, mean_eh2_eq
  real(wp)               :: mean_r, mean_r2, mean_r_eq, mean_r2_eq
  real(wp)               :: mean_t, mean_t2, mean_t_eq, mean_t2_eq
  complex(wp)            :: sumw_l, sumw_h
contains
  procedure              :: init => init_rare_events
  procedure              :: reader => rare_events_reader
  generic                :: detect => detect_rare_events, detect_rare_events_1d
  procedure              :: adjust => adjust_rev_energy
  procedure              :: get_ediff_l => get_ediff_l
  procedure              :: finalize => finalize_rare_events
  procedure              :: finalize_eq => finalize_rare_events_eq
  procedure              :: report_setup => report_rare_events_setup
  procedure              :: report => report_rare_events
  procedure, private     :: detect_rare_events, detect_rare_events_1d
end type rare_events

contains

!******************************************************************************** 
!
! Initialize the structure for rare events
!
!******************************************************************************** 
subroutine init_rare_events(self, hdes, mc_des, tau)
  class(rare_events), intent(out)    :: self
  type(hamil_descriptor), intent(in) :: hdes
  type(McDescriptor), intent(in)     :: mc_des
  real(wp), intent(in)               :: tau
  
  self%active = .true.
  self%p_capping = .true.
  self%r_capping = .false.
  self%t_capping = .false.
  self%el_capping = .true.
  self%eh_capping = .false.
  self%i_capping = .true.
  self%w_capping = .true.
  self%el_capping_mode = 1

  self%r_ratio = 5.0_wp
  self%t_ratio = 5.0_wp
  self%el_ratio = 1.0_wp
  self%eh_ratio = 1.0_wp
  self%width_r = 0.0_wp
  self%width_t = 0.0_wp
  self%width_el = 0.0_wp
  self%width_eh = 0.0_wp

  self%impwmax = 10.0_wp
  self%impwmin = 0.005_wp
  self%rev_lifemax = 3
  self%pers_lifemax = 50
  self%maxweight = min(100.0_wp, 0.1_wp*real(mc_des%nwalkers, wp))
  if (self%maxweight < 10.0_wp) self%maxweight = 10.0_wp
  self%sumw_l = zeroc
  self%sumw_h = zeroc

  self%nwalkers =  mc_des%nwalkers
  self%sample = mc_des%sample
  self%ng = hdes%ng

  call self%reader

  self%number_p = 0; self%number_pu = 0
  self%number_ru = 0; self%number_ruu = 0
  self%number_rl = 0; self%number_rlu = 0
  self%number_tu = 0; self%number_tuu = 0
  self%number_tl = 0; self%number_tlu = 0
  self%number_elu = 0; self%number_eluu = 0
  self%number_ell = 0; self%number_ellu = 0
  self%number_ehu = 0; self%number_ehuu = 0
  self%number_ehl = 0; self%number_ehlu = 0
  self%number_iu = 0; self%number_iuu = 0
  self%number_il = 0; self%number_ilu = 0
  self%number_w  = 0; self%number_wu = 0
  self%ndeath_rev = 0; self%ndeath_revu = 0
  self%ndeath_pers = 0; self%ndeath_persu = 0
  self%number_acc = 0
  
  self%mean_r = 0.0_wp; self%mean_r2 = 0.0_wp; self%mean_r_eq = 0.0_wp; self%mean_r2_eq = 0.0_wp
  self%mean_t = 0.0_wp; self%mean_t2 = 0.0_wp; self%mean_t_eq = 0.0_wp; self%mean_t2_eq = 0.0_wp
  self%ediff_l = 0.0_wp; self%ediff_h = 0.0_wp
  self%mean_el = 0.0_wp; self%mean_el2 = 0.0_wp; self%mean_el_eq = 0.0_wp; self%mean_el2_eq = 0.0_wp
  self%mean_eh = 0.0_wp; self%mean_eh2 = 0.0_wp; self%mean_eh_eq = 0.0_wp; self%mean_eh2_eq = 0.0_wp
  
  if (abs(self%width_r) < 1.0e-06_wp) self%width_r = sqrt(sum(hdes%nel) * tau) / 4.0_wp
  if (abs(self%width_t) < 1.0e-06_wp) self%width_t = sqrt(tau) / 4.0_wp
  if (abs(self%width_el) < 1.0e-06_wp) then
    if (self%el_capping_mode == 2) then
      self%width_el = sqrt(2.0_wp/tau)
    else
      self%width_el = 0.5_wp*sqrt(sum(hdes%nel)/tau) + sqrt(sum(hdes%nel)*tau)
    end if
  end if
  if (abs(self%width_eh) < 1.0e-06_wp) self%width_eh = self%width_el
end subroutine init_rare_events


!******************************************************************************** 
!
! Reader for the rare_events structure
!
!******************************************************************************** 
subroutine rare_events_reader(self)
  use qmcfort_in, only: add_input
  class(rare_events), intent(inout) :: self
  
  call add_input("rev",             self%active)
  call add_input("lrev",            self%active)
  call add_input("rare_events",     self%active)
  call add_input("p_capping",       self%p_capping)
  call add_input("r_capping",       self%r_capping)
  call add_input("t_capping",       self%t_capping)
  call add_input("el_capping",      self%el_capping)
  call add_input("eh_capping",      self%eh_capping)
  call add_input("w_capping",       self%w_capping)
  call add_input("i_capping",       self%i_capping)
  call add_input("r_capping_ratio", self%r_ratio)
  call add_input("t_capping_ratio", self%t_ratio)
  call add_input("el_capping_ratio",self%el_ratio)
  call add_input("eh_capping_ratio",self%eh_ratio)
  call add_input("el_capping_mode", self%el_capping_mode)
  call add_input("cap_width_el",    self%width_el)
  call add_input("cap_width_eh",    self%width_eh)
  call add_input("rev_lifemax",     self%rev_lifemax)
  call add_input("pers_lifemax",    self%pers_lifemax)
end subroutine rare_events_reader


!******************************************************************************** 
!
! Detect rare events
!
! Chierarchy of rare events:
!    1. phase rare events                  p_capping 
!    2. weight rare events                 w_capping
!    3. local energy rare events           el_capping
!    4. hybrid energy rare events          eh_capping
!    5. reweighting factor rare events     i_capping
!    6. overlap ratio rare events          r_capping
!    7. phase change rare events           t_capping
!
!******************************************************************************** 
subroutine detect_rare_events(self, af_des, af_proj, mean_energy, grow, tau, walker)
  class(rare_events), intent(inout) :: self
  type(AfqmcDescriptor), intent(in) :: af_des
  type(AfqmcProj), intent(in)       :: af_proj
  type(CEnergy), intent(in)         :: mean_energy
  type(GrowthEst), intent(in)       :: grow
  real(wp), intent(in)              :: tau
  type(AfqmcWalker), intent(inout)  :: walker
  !local_variables
  logical                           :: rare_event,  p_event, rti_event, death_rev, lsample, lsamplex, equil
  real(wp)                          :: mean_r, mean_t, mean_el, mean_eh
  real(wp)                          :: width_r, width_t, width_el, width_eh
  real(wp)                          :: diff, val
  real(wp)                          :: mu

  if (.not. self%active) return

  rare_event = .false.
  p_event = .false.
  rti_event = .false.
  death_rev = .false.

  lsamplex = modul(af_des%mc_des%step, af_des%mc_des%sample*af_des%mc_des%samplex)
  lsample = modul(af_des%mc_des%step, af_des%mc_des%sample)

  mean_r = 1.0_wp
  width_r = self%r_ratio * self%width_r

  mean_t = 0.0_wp
  width_t = self%t_ratio * self%width_t

  mean_eh = real(mean_energy%electron_energy(), wp)
  width_eh = self%eh_ratio * self%width_eh

  equil = af_des%mc_des%is_eq()

  if (equil) self%sumw_h = self%sumw_h + af_proj%sample_weight(walker)
  if (equil .and. lsample) self%sumw_l = self%sumw_l + af_proj%sample_weight(walker)

  if (lsample) then
    if (lsamplex) then
      mean_el = real(mean_energy%electron_energy(), wp)
      width_el = self%el_ratio * self%width_el
    else
      mean_el = real(mean_energy%electron_energy()-mean_energy%ex, wp)
      width_el = 1.5_wp * self%el_ratio * self%width_el
    end if
  end if

  !1. phase capping
  if (self%p_capping) then
    if (abs(walker%dphase) >= pi/2.0_wp) then 
      walker%weight = 0.0_wp
      if (equil) then
        self%number_p = self%number_p + 1
        if (.not. rare_event) then
          self%number_pu = self%number_pu + 1
          rare_event = .true.
          p_event = .true.
        end if
      end if
    end if
  end if

  !2. weight capping
  if (self%w_capping) then
    if (walker%weight > self%maxweight) then
      walker%weight = self%maxweight
      if (equil) then
        self%number_w = self%number_w + 1 
        if (.not. rare_event) then
          self%number_wu = self%number_wu + 1
          rare_event = .true.
        end if
      end if
    end if
  end if

  !3. local energy capping
  if (lsample .and. self%el_capping) then
    if (real(walker%energyl,wp) > mean_el + width_el) then
      diff = real(walker%energyl,wp) - (mean_el+width_el)
      walker%energyl = cmplx(mean_el+width_el, aimag(walker%energyl))
      call self%adjust(walker)
      if (equil) then 
        self%ediff_l = self%ediff_l + af_proj%sample_weight(walker)*diff
        self%number_elu = self%number_elu + 1
        if (.not. rare_event) then
          self%number_eluu = self%number_eluu + 1
          rare_event = .true.
        end if
      end if
    else if (real(walker%energyl,wp) < mean_el - width_el) then
      diff = real(walker%energyl,wp) - (mean_el-width_el)
      walker%energyl = cmplx(mean_el-width_el, aimag(walker%energyl))
      call self%adjust(walker)
      if (equil) then 
        self%ediff_l = self%ediff_l + af_proj%sample_weight(walker)*diff
        self%number_ell = self%number_ell + 1
        if (.not. rare_event) then
          self%number_ellu = self%number_ellu + 1
          rare_event = .true.
        end if
      end if
    end if
  end if

  !4. hybrid energy capping
  if (self%eh_capping) then
    if (real(walker%energyh,wp) > mean_eh+width_eh) then
      diff = real(walker%energyh,wp) - mean_eh - width_eh
      !walker%energyh = cmplx(mean_eh+width_eh, aimag(walker%energyh))
      if (equil) then 
        self%ediff_h = self%ediff_h + af_proj%sample_weight(walker)*diff
        self%number_ehu = self%number_ehu + 1
        if (.not. rare_event) then
          self%number_ehuu = self%number_ehuu + 1
          rare_event = .true.
        end if
      end if
    else if (real(walker%energyh,wp) < mean_eh-width_eh) then
      diff = real(walker%energyh,wp) - mean_eh + width_eh
      walker%energyh = cmplx(mean_eh-width_eh, aimag(walker%energyh))
      if (equil) then 
        self%ediff_h = self%ediff_h + af_proj%sample_weight(walker)*diff
        self%number_ehl = self%number_ehl + 1
        if (.not. rare_event) then
          self%number_ehlu = self%number_ehlu + 1
          rare_event = .true.
        end if
      end if
    end if
  end if

  mu = real(mean_energy%electron_energy(), wp) + grow%energy
  call af_proj%reweighting_factor(walker, tau, mu)

  !5. importance function capping
  if (self%i_capping) then
    if (walker%rew_w > self%impwmax) then
      walker%rew_w = 0.0_wp
      if (equil) then
        self%number_iu = self%number_iu + 1
        if (.not. rare_event) then
          self%number_iuu = self%number_iuu + 1
          rare_event = .true.
        end if
      end if
    else if (walker%rew_w < self%impwmin) then
      walker%rew_w = 0.0_wp
      if (equil) then
        self%number_il = self%number_il + 1
        if (.not. rare_event) then
          self%number_ilu = self%number_ilu + 1
          rare_event = .true.
        end if
      end if
    end if
  end if
  
  !6. overlap ratio capping
  if (self%r_capping) then
    if (walker%rew_w > mean_r+width_r) then
      if (.not. p_event) rti_event = .true.
      if (equil) then
        self%number_ru = self%number_ru + 1
        if (.not. rare_event) then
          self%number_ruu = self%number_ruu + 1
          rare_event = .true.
        end if
      end if
    else if (walker%rew_w < mean_r-width_r) then
      if (.not. p_event)  rti_event = .true.
      if (equil) then
        self%number_rl = self%number_rl + 1
        if (.not. rare_event) then
          self%number_rlu = self%number_rlu + 1
          rare_event = .true.
        end if
      end if
    end if
  end if

  !7. phase change capping
  if (self%t_capping) then
    if (walker%rew_p > mean_t+width_t) then
      if (.not. p_event)  rti_event = .true.
      if (equil) then
        self%number_tu = self%number_tu + 1
        if (.not. rare_event) then
          self%number_tuu = self%number_tuu + 1
          rare_event = .true.
        end if
      end if
    else if (walker%rew_p < mean_t-width_t) then
      if (.not. p_event)  rti_event = .true.
      if (equil) then
        self%number_tl = self%number_tl + 1
        if (.not. rare_event) then
          self%number_tlu = self%number_tlu + 1
          rare_event = .true.
        end if
      end if
    end if
  end if

  !update life variable and decide whether the walker dies
  if (rti_event) then 
    walker%revlife = walker%revlife + 1
    if (walker%revlife == self%rev_lifemax) then 
      walker%weight = 0.0_wp
      self%ndeath_rev = self%ndeath_rev + 1
      if (.not. p_event) self%ndeath_revu = self%ndeath_revu + 1
      death_rev = .true.
    end if
  else 
    walker%revlife = max(0, walker%revlife-1)
  end if

  !8. acceptance rate
  if (.not. walker%move) then
    if (equil) self%number_acc = self%number_acc + 1
    walker%life = walker%life + 1
    if (walker%life == self%pers_lifemax) then
      walker%weight = 0.0_wp
      if (equil) then
        self%ndeath_pers = self%ndeath_pers + 1
        if (.not. (p_event .or. death_rev)) self%ndeath_persu = self%ndeath_persu + 1
      end if
    end if
  else 
    walker%life = 0
  end if

  !sample self values
  if (equil) then
    val = cap_value(walker%rew_w, self%mean_r_eq, self%mean_r2_eq, self%r_ratio)
    self%mean_r = self%mean_r + val
    self%mean_r2 = self%mean_r2 + val**2
    val = cap_value(abs(walker%rew_p), self%mean_t_eq, self%mean_t2_eq, self%t_ratio)
    self%mean_t = self%mean_t + val
    self%mean_t2 = self%mean_t2 + val**2
    val = cap_value(real(walker%energyh,wp), self%mean_eh_eq, self%mean_eh2_eq, 5.0_wp)
    self%mean_eh  = self%mean_eh + val
    self%mean_eh2 = self%mean_eh2 + val**2
    if (lsamplex) then
      val = cap_value(real(walker%energyl,wp), self%mean_el_eq, self%mean_el2_eq, 5.0_wp)
      self%mean_el  = self%mean_el + val
      self%mean_el2 = self%mean_el2 + val**2
    end if
  else 
    self%mean_r_eq = self%mean_r_eq + walker%rew_w
    self%mean_r2_eq = self%mean_r2_eq + walker%rew_w**2
    self%mean_t_eq = self%mean_t_eq + walker%rew_p
    self%mean_t2_eq = self%mean_t2_eq + walker%rew_p**2
    self%mean_eh_eq = self%mean_eh_eq + real(walker%energyh,wp)
    self%mean_eh2_eq = self%mean_eh2_eq + real(walker%energyh,wp)**2
    if (lsamplex) then
      self%mean_el_eq  = self%mean_el_eq + real(walker%energyl,wp)
      self%mean_el2_eq = self%mean_el2_eq + real(walker%energyl,wp)**2
    end if
  end if
end subroutine detect_rare_events


!******************************************************************************** 
!
! Detect rare events (working with an array of walkers)
!
!******************************************************************************** 
subroutine detect_rare_events_1d(self, af_des, af_proj, mean_energy, grow, tau, walkers)
  class(rare_events), intent(inout) :: self
  type(AfqmcDescriptor), intent(in) :: af_des
  type(AfqmcProj), intent(in)       :: af_proj
  type(CEnergy), intent(in)         :: mean_energy
  type(GrowthEst), intent(in)       :: grow
  real(wp), intent(in)              :: tau
  type(AfqmcWalker), intent(inout)  :: walkers(:)
  !local_variables
  integer                           :: w
  character(len=*), parameter       :: proc_name = "detect_rare_events"

  if (use_profiling) call start_profiling(proc_name)

  do w = 1, size(walkers)
    call self%detect(af_des, af_proj, mean_energy, grow, tau, walkers(w))
  end do

  if (use_profiling) call end_profiling(proc_name)
end subroutine detect_rare_events_1d


!******************************************************************************** 
!
! Adjust walker energy if it is the rare event
!
!******************************************************************************** 
subroutine adjust_rev_energy(self, walker)
  class(rare_events), intent(inout) :: self
  type(AfqmcWalker), intent(inout)  :: walker

!debug: energy_types
  !call walker%energy%set_ediff(real((walker%energyl - walker%energy%electron_energy()), wp))
  !walker%energy%e = walker%energyl
  !call walker%energy%set_toten_()
end subroutine adjust_rev_energy


!******************************************************************************** 
!
! Obtain ediff_l during the sampling (before finalize_rare_events is called)
!
!******************************************************************************** 
function get_ediff_l(self, comm) result(ediff_l)
  class(rare_events), intent(in)     :: self
  type(mpi_communicator), intent(in) :: comm
  real(wp)                           :: ediff_l
  !local variables
  complex(wp)                        :: sumw_l

  sumw_l = self%sumw_l
  ediff_l = self%ediff_l

  call comm%mpisum(ediff_l)
  call comm%mpisum(sumw_l)

  ediff_l = ediff_l / sumw_l
end function get_ediff_l


!******************************************************************************** 
!
! Collect rare event across all processes
!
!******************************************************************************** 
subroutine finalize_rare_events(self, mc_des)
  class(rare_events), intent(inout) :: self
  type(McDescriptor), intent(in)    :: mc_des
  !local variables
  integer                           :: nsamples, nsamples_

  call comm_world%mpisum(self%number_p); call comm_world%mpisum(self%number_pu)
  call comm_world%mpisum(self%number_ru); call comm_world%mpisum(self%number_ruu)
  call comm_world%mpisum(self%number_rl); call comm_world%mpisum(self%number_rlu)
  call comm_world%mpisum(self%number_tu); call comm_world%mpisum(self%number_tuu)
  call comm_world%mpisum(self%number_tl); call comm_world%mpisum(self%number_tlu)
  call comm_world%mpisum(self%number_elu); call comm_world%mpisum(self%number_eluu)
  call comm_world%mpisum(self%number_ell); call comm_world%mpisum(self%number_ellu)
  call comm_world%mpisum(self%number_ehu); call comm_world%mpisum(self%number_ehuu)
  call comm_world%mpisum(self%number_ehl); call comm_world%mpisum(self%number_ehlu)
  call comm_world%mpisum(self%number_iu); call comm_world%mpisum(self%number_iuu)
  call comm_world%mpisum(self%number_il); call comm_world%mpisum(self%number_ilu)
  call comm_world%mpisum(self%number_w); call comm_world%mpisum(self%number_wu)
  call comm_world%mpisum(self%ndeath_rev); call comm_world%mpisum(self%ndeath_revu)
  call comm_world%mpisum(self%ndeath_pers); call comm_world%mpisum(self%ndeath_persu)
  call comm_world%mpisum(self%number_acc)

  call comm_world%mpisum(self%mean_r)
  call comm_world%mpisum(self%mean_r2)
  call comm_world%mpisum(self%mean_t)
  call comm_world%mpisum(self%mean_t2)
  call comm_world%mpisum(self%mean_el)
  call comm_world%mpisum(self%mean_el2)
  call comm_world%mpisum(self%mean_eh)
  call comm_world%mpisum(self%mean_eh2)
  call comm_world%mpisum(self%ediff_l)
  call comm_world%mpisum(self%ediff_h)
  call comm_world%mpisum(self%sumw_l)
  call comm_world%mpisum(self%sumw_h)

  self%nsteps = mc_des%step
  nsamples = mc_des%steps_samp() * mc_des%nwalkers
  nsamples_ = nsamples / mc_des%sample / mc_des%samplex

  self%mean_r = self%mean_r / nsamples
  self%mean_t = self%mean_t / nsamples
  self%mean_el = self%mean_el / nsamples_
  self%mean_eh = self%mean_eh / nsamples
  self%ediff_l = self%ediff_l / real(self%sumw_l, wp)
  self%ediff_h = self%ediff_h / real(self%sumw_h, wp)

  self%mean_r2 = sqrt(self%mean_r2/nsamples - self%mean_r**2)
  self%mean_t2 = sqrt(self%mean_t2/nsamples - self%mean_t**2)
  self%mean_el2 = sqrt(self%mean_el2/nsamples_ - self%mean_el**2)
  self%mean_eh2 = sqrt(self%mean_eh2/nsamples - self%mean_eh**2)
end subroutine finalize_rare_events


!******************************************************************************** 
!
! Collect rare event averages after equilibration      
!
!******************************************************************************** 
subroutine finalize_rare_events_eq(self, mc_des)
  class(rare_events), intent(inout) :: self
  type(McDescriptor), intent(in)    :: mc_des
  !local variables
  integer                           :: nsamples, nsamples_

  call comm_world%mpisum(self%mean_r_eq)
  call comm_world%mpisum(self%mean_r2_eq)
  call comm_world%mpisum(self%mean_t_eq)
  call comm_world%mpisum(self%mean_t2_eq)
  call comm_world%mpisum(self%mean_el_eq)
  call comm_world%mpisum(self%mean_el2_eq)
  call comm_world%mpisum(self%mean_eh_eq)
  call comm_world%mpisum(self%mean_eh2_eq)

  self%nsteps_eq = mc_des%steps_eq()
  nsamples = mc_des%steps_eq() * mc_des%nwalkers
  nsamples_ = nsamples / mc_des%sample / mc_des%samplex

  self%mean_r_eq = self%mean_r_eq / nsamples
  self%mean_t_eq = self%mean_t_eq / nsamples
  self%mean_el_eq = self%mean_el_eq / nsamples_
  self%mean_eh_eq = self%mean_eh_eq / nsamples

  self%mean_r2_eq = sqrt(self%mean_r2_eq/nsamples - self%mean_r_eq**2)
  self%mean_t2_eq = sqrt(self%mean_t2_eq/nsamples - self%mean_t_eq**2)
  self%mean_el2_eq = sqrt(self%mean_el2_eq/nsamples_ - self%mean_el_eq**2)
  self%mean_eh2_eq = sqrt(self%mean_eh2_eq/nsamples - self%mean_eh_eq**2)
end subroutine finalize_rare_events_eq


!******************************************************************************** 
!
! Final rare events io 
!
!******************************************************************************** 
subroutine report_rare_events(self, fh)
  class(rare_events), intent(inout) :: self
  type(FileHandle), intent(in)      :: fh
  !local variables
  real(wp)                          :: nsamp_w, nsamp_we
  integer                           :: g, i

  nsamp_w = max(real(self%nsteps * self%nwalkers), 1.0_wp)
  nsamp_we = max(nsamp_w/self%sample, 1.0_wp)
  
  write(fh%funit,*)
  write(fh%funit,*)   "Rare events statistics:"
  write(fh%funit,*)   "-----------------------"

  !1. phaseless capping
  write(fh%funit,100) "phaseless rare events             = ", self%p_capping
  write(fh%funit,101) "  # of rare events                = ", self%number_p, self%number_pu, &
                      write_prom(self%number_pu/nsamp_w)
  write(fh%funit,*)

  !2. overlap ratio capping
  write(fh%funit,100) "overlap ratio capping             = ", self%r_capping
  write(fh%funit,103) "  sampled mean                    = ", self%mean_r, self%mean_r_eq
  write(fh%funit,103) "  sampled stddev                  = ", self%mean_r2, self%mean_r2_eq
  write(fh%funit,103) "  capping width                   = ", self%width_r
  write(fh%funit,103) "  capping ratio                   = ", self%r_ratio
  write(fh%funit,101) "  # of rare events upper          = ", self%number_ru, self%number_ruu, &
                      write_prom(self%number_ruu/nsamp_w)
  write(fh%funit,101) "  # of rare events lower          = ", self%number_rl, self%number_rlu, &
                      write_prom(self%number_rlu/nsamp_w)
  write(fh%funit,*)

  !3. phase change capping
  write(fh%funit,100) "phase change capping              = ", self%t_capping
  write(fh%funit,103) "  sampled mean                    = ", self%mean_t, self%mean_t_eq
  write(fh%funit,103) "  sampled stddev                  = ", self%mean_t2, self%mean_t2_eq
  write(fh%funit,103) "  capping width                   = ", self%width_t
  write(fh%funit,103) "  capping ratio                   = ", self%t_ratio
  write(fh%funit,101) "  # of rare events upper          = ", self%number_tu, self%number_tuu, &
                      write_prom(self%number_tuu/nsamp_w)
  write(fh%funit,101) "  # of rare events lower          = ", self%number_tl, self%number_tlu, &
                      write_prom(self%number_tlu/nsamp_w)
  write(fh%funit,*)

  !persistent events
  write(fh%funit,100) "persistent death events"
  write(fh%funit,105) "  maximal rare event lifetime     = ", self%rev_lifemax
  write(fh%funit,101) "  # of death events               = ", self%ndeath_rev, self%ndeath_revu, &
                      write_prom(self%ndeath_revu/nsamp_w)
  write(fh%funit,*)

  !4. weight capping
  write(fh%funit,100) "weight rare events                = ", self%w_capping
  write(fh%funit,103) "  upper value                     = ", self%maxweight
  write(fh%funit,101) "  # of rare events                = ", self%number_w, self%number_wu, &
                      write_prom(self%number_wu/nsamp_w)
  write(fh%funit,*)

  !5. local energy capping
  write(fh%funit,100) "local energy capping              = ", self%el_capping
  write(fh%funit,103) "  sampled mean                    = ", self%mean_el, self%mean_el_eq
  write(fh%funit,103) "  sampled stddev                  = ", self%mean_el2, self%mean_el2_eq
  write(fh%funit,103) "  capping width                   = ", self%width_el
  write(fh%funit,103) "  capping ratio                   = ", self%el_ratio
  write(fh%funit,104) "  bare-capped local energy        = ", self%ediff_l
  write(fh%funit,101) "  # of rare events upper          = ", self%number_elu, self%number_eluu, &
                      write_prom(self%number_eluu/nsamp_we)
  write(fh%funit,101) "  # of rare events lower          = ", self%number_ell, self%number_ellu, &
                      write_prom(self%number_ellu/nsamp_we)
  write(fh%funit,*)

  !6. hybrid energy capping
  write(fh%funit,100) "hybrid energy capping             = ", self%eh_capping
  write(fh%funit,103) "  sampled mean                    = ", self%mean_eh, self%mean_eh_eq
  write(fh%funit,103) "  sampled stddev                  = ", self%mean_eh2, self%mean_eh2_eq
  write(fh%funit,103) "  capping width                   = ", self%width_eh
  write(fh%funit,103) "  capping ratio                   = ", self%eh_ratio
  write(fh%funit,104) "  bare-capped hybrid energy       = ", self%ediff_h
  write(fh%funit,101) "  # of rare events upper          = ", self%number_ehu, self%number_ehuu, &
                      write_prom(self%number_ehuu/nsamp_w)
  write(fh%funit,101) "  # of rare events lower          = ", self%number_ehl, self%number_ehlu, &
                      write_prom(self%number_ehlu/nsamp_w)
  write(fh%funit,*)

  !7. importance function capping
  write(fh%funit,100) "imp function rare events          = ", self%i_capping
  write(fh%funit,103) "  upper value                     = ", self%impwmax
  write(fh%funit,103) "  lower value                     = ", self%impwmin
  write(fh%funit,101) "  # of rare events  upper         = ", self%number_iu, self%number_iuu, &
                      write_prom(self%number_iuu/nsamp_w)
  write(fh%funit,101) "  # od rare events lower          = ", self%number_il, self%number_ilu, &
                      write_prom(self%number_ilu/nsamp_w)
  write(fh%funit,*)

  !8. Acceptance rate 
  write(fh%funit,100) "Accept/Reject step"
  write(fh%funit,102) "  acceptance rate                 = ", self%number_acc, write_prom(self%number_acc/nsamp_w) 
  write(fh%funit,105) "  maximal rare event lifetime     = ", self%pers_lifemax
  write(fh%funit,101) "  # of death events               = ", self%ndeath_pers, self%ndeath_persu, &
                      write_prom(self%ndeath_persu/nsamp_w)
  write(fh%funit,*)

  100 format (1x,a,l6)
  101 format (1x,a,2i10," ",a)
  102 format (1x,a,i10," ",a)
  103 format (1x,a,2f14.6)
  104 format (1x,a,es14.6)
  105 format (1x,a,i6)
end subroutine report_rare_events


!******************************************************************************** 
!
! Rare events setup io 
!
!******************************************************************************** 
subroutine report_rare_events_setup(self, fh)
  class(rare_events), intent(in) :: self
  type(FileHandle), intent(in)   :: fh
  
  write(fh%funit,*)   "Rare Events setup"
  write(fh%funit,*)   "-----------------"
  write(fh%funit,100) "  rare events activated           ", self%active
  write(fh%funit,100) "  phase capping                   ", self%p_capping
  write(fh%funit,100) "  overlap ratio capping           ", self%r_capping
  write(fh%funit,100) "  phase change capping            ", self%t_capping
  write(fh%funit,100) "  local energy capping            ", self%el_capping
  write(fh%funit,100) "  hybrid energy capping           ", self%eh_capping
  write(fh%funit,100) "  importance weight capping       ", self%i_capping
  write(fh%funit,100) "  weight capping                  ", self%w_capping
  
  100 format (1x,t5,a,t50,"= ",10l6)
  101 format (1x,t5,a,t50,"= ",a)
end subroutine report_rare_events_setup

end module afqmc_rev