! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module afqmc_proj
!******************************************************************************** 
!
! Implementation of the AfqmcProj
!
!    Behavior-oriented class used to:
!        - calculate reweighting factor
!        - update walker weight
!        - calculate the sampling weight and energies used for averages
!
!******************************************************************************** 

#include "../preproc.inc"

use constants
use profiling

use afqmc_proj_config, only: AfqmcProjConfig
use afqmc_walker, only: AfqmcWalker
use energy_types, only: CEnergy
use standalone, only: phase

implicit none

private
public AfqmcProjConfig, AfqmcProj

type AfqmcProj
  type(AfqmcProjConfig)                   :: cfg
  real(wp)                                :: tau

  procedure(Iupdate_weight), pointer      :: update_weight
  procedure(Isample_weight), pointer      :: sample_weight
  procedure(Isample_energy_val), pointer  :: sample_energy_val
  procedure(Isample_energy_type), pointer :: sample_energy_type
contains
  generic            :: assignment(=) => assign_afqmc_proj
  procedure, private :: assign_afqmc_proj

  procedure          :: report => report_afqmc_proj

  procedure          :: reweighting_factor => get_reweighting_factor
  procedure          :: update_weights => update_walker_weights

  procedure, private :: get_hybrid

  procedure, private :: update_weight_factory
  procedure, private :: update_weight_free, update_weight_phaseless 

  procedure, private :: sample_weight_factory
  procedure, private :: sample_weight_free, sample_weight_free_hp, sample_weight_phaseless

  procedure, private :: sample_energy_val_factory
  procedure, private :: sample_energy_val_free, sample_energy_val_phaseless

  procedure, private :: sample_energy_type_factory
  procedure, private :: sample_energy_type_free, sample_energy_type_phaseless
end type AfqmcProj

interface AfqmcProj
  module procedure finit_afqmc_proj
end interface AfqmcProj

abstract interface
  subroutine Iupdate_weight(self, walker)
    import AfqmcProj, AfqmcWalker
    class(AfqmcProj), intent(in)     :: self
    type(AFqmcWalker), intent(inout) :: walker
  end subroutine Iupdate_weight

  function Isample_weight(self, walker) result(sample_weight)
    import AfqmcProj, AfqmcWalker, wp
    class(AfqmcProj), intent(in)  :: self
    type(AFqmcWalker), intent(in) :: walker
    complex(wp)                   :: sample_weight
  end function Isample_weight

  function Isample_energy_val(self, energy) result(sample_energy)
    import AfqmcProj, wp
    class(AfqmcProj), intent(in) :: self
    complex(wp), intent(in)      :: energy
    complex(wp)                  :: sample_energy
  end function Isample_energy_val

  function Isample_energy_type(self, energy) result(sample_energy)
    import AfqmcProj, CEnergy
    class(AfqmcProj), intent(in) :: self
    type(CEnergy), intent(in)    :: energy
    type(CEnergy)                :: sample_energy
  end function Isample_energy_type
end interface

contains

!******************************************************************************** 
! 
! Initialization of the AfqmcProj object
!
!******************************************************************************** 
subroutine init_afqmc_proj(self, cfg, tau)
  type(AfqmcProj), intent(out)      :: self
  type(AfqmcProjConfig), intent(in) :: cfg
  real(wp), intent(in)              :: tau

  self%tau = tau
  self%cfg = cfg

  !setup procedure pointers explicitly
  call self%update_weight_factory()
  call self%sample_weight_factory()
  call self%sample_energy_val_factory()
  call self%sample_energy_type_factory()
end subroutine init_afqmc_proj


!******************************************************************************** 
! 
! AfqmcProj constructor
!
!******************************************************************************** 
function finit_afqmc_proj(cfg, tau) result(self)
  type(AfqmcProjConfig), intent(in) :: cfg
  type(AfqmcProj)                   :: self
  real(wp), intent(in)              :: tau

  call init_afqmc_proj(self, cfg, tau)
end function finit_afqmc_proj


!******************************************************************************** 
!
! Assignment Operator for the AfqmcProj type
!
! User-defined assignment defined to ensure procedure pointers are setup correctly.
!
!******************************************************************************** 
subroutine assign_afqmc_proj(to, from)
  class(AfqmcProj), intent(inout) :: to
  type(AfqmcProj), intent(in)     :: from
  
  to%tau = from%tau
  to%cfg = from%cfg

  !setup procedure pointers explicitly
  call to%update_weight_factory()
  call to%sample_weight_factory()
  call to%sample_energy_val_factory()
  call to%sample_energy_type_factory()
end subroutine assign_afqmc_proj


!******************************************************************************** 
!
! Report AfqmcProj object
!
!******************************************************************************** 
subroutine report_afqmc_proj(self, funit)
  class(AfqmcProj), intent(in) :: self
  integer, intent(in)          :: funit

  write(funit,*)
  write(funit,100) "Afqmc Projection (AfqmcProj):"
  write(funit,100) "-----------------------------"
  write(funit,101) "  projection", self%cfg%projection
  write(funit,103) "  auto_hybrid", self%cfg%auto_hybrid 
  write(funit,102) "  hybrid", self%get_hybrid()
  write(funit,102) "  mix", self%cfg%mix
  write(funit,102) "  min_cos_hp", self%cfg%min_cos_hp
  write(funit,103) "  accept reject", self%cfg%accept_reject

  100 format (1x,a)
  101 format (1x,a,t50,"= ",a)
  102 format (1x,a,t50,"= ",f12.4)
  103 format (1x,a,t50,"= ",l)
end subroutine report_afqmc_proj


!******************************************************************************** 
! 
! Depending on self%cfg%auto_hybrid calculates the portion of hybrid energy to 
! be used in AFQMC reweighting
! 
!******************************************************************************** 
function get_hybrid(self) result(hybrid)
  class(AfqmcProj), intent(in) :: self
  real(wp)                     :: hybrid

  if (self%cfg%auto_hybrid) then
    hybrid = 1.0_wp / sqrt(1.0_wp + 2.0_wp*self%tau)
  else
    hybrid = self%cfg%hybrid
  end if
end function get_hybrid


!******************************************************************************** 
!
! Calculate AFQMC reweighting factor
!
!    E = h (mix*Eh(i) + (1-mix)Eh(i-1)) + (1-h) (mix*El(i) + (1-mix)El(i-1))
!    impw = exp(-dt * (E-mu))
!
!******************************************************************************** 
subroutine get_reweighting_factor(self, walker, tau, mu)
  class(AfqmcProj), intent(in)     :: self
  type(AfqmcWalker), intent(inout) :: walker
  real(wp), intent(in)             :: tau
  real(wp), intent(in)             :: mu
  !local
  real(wp)                         :: hybrid
  complex(wp)                      :: rew

  hybrid = self%get_hybrid()

  associate(el=>walker%energyl, el_old=>walker%energyl_old, &
            eh=>walker%energyh, eh_old=>walker%energyh_old)

  rew = (1.0_wp - hybrid) * (self%cfg%mix*el + (1.0_wp-self%cfg%mix)*el_old) &
                + hybrid  * (self%cfg%mix*eh + (1.0_wp-self%cfg%mix)*eh_old) 

  rew = exp(-tau * (rew-mu))
  end associate

  walker%rew_w = abs(rew)
  walker%rew_p = phase(rew)
!debug: vmc try
!walker%rew_w = 1.0_wp
!walker%rew_p = 0.0_wp
end subroutine get_reweighting_factor


!******************************************************************************** 
!
! Factory for the update_weight method
!
!******************************************************************************** 
subroutine update_weight_factory(self)
  class(AfqmcProj), intent(inout) :: self

  select case (trim(self%cfg%projection))
    case ("phaseless")
      self%update_weight => update_weight_phaseless
    case("free")
      self%update_weight => update_weight_free
    case ("free_hp")
      self%update_weight => update_weight_free
    case default
      self%update_weight => update_weight_phaseless
  end select
end subroutine update_weight_factory


!******************************************************************************** 
!
! Update weight of the AfqmcWalker (Free projection)
!
!    w *= rew_w
!    phase += rew_p
!
!******************************************************************************** 
subroutine update_weight_free(self, walker)
  class(AfqmcProj), intent(in)     :: self
  type(AfqmcWalker), intent(inout) :: walker

    walker%weight = walker%weight * walker%rew_w
    walker%bare_weight = walker%bare_weight * walker%rew_w
    walker%phase = walker%phase + walker%rew_p
end subroutine update_weight_free


!******************************************************************************** 
!
! Update weight of the AfqmcWalker (phaseless projection)
!
!    w *= rew_w * max(0, cos(dphase))
!    phase = 0
!
!******************************************************************************** 
subroutine update_weight_phaseless(self, walker)
  class(AfqmcProj), intent(in)     :: self
  type(AfqmcWalker), intent(inout) :: walker

  call self%update_weight_free(walker)
  walker%weight = walker%weight * max(0.0_wp, cos(walker%dphase))
  walker%phase = 0.0_wp
end subroutine update_weight_phaseless


!******************************************************************************** 
!
! Update weight for the array of walkers 
!
!******************************************************************************** 
subroutine update_walker_weights(self, walkers)
  class(AfqmcProj), intent(in)     :: self
  type(AfqmcWalker), intent(inout) :: walkers(:)
  !local
  integer                          :: w
  character(len=*), parameter      :: proc_name = "update_walker_weights"

  if (use_profiling) call start_profiling(proc_name)

  do w = 1, size(walkers)
    call self%update_weight(walkers(w))
  end do

  if (use_profiling) call end_profiling(proc_name)
end subroutine update_walker_weights


!******************************************************************************** 
!
! Factory for the sample_weight method
!
!******************************************************************************** 
subroutine sample_weight_factory(self)
  class(AfqmcProj), intent(inout) :: self

  select case (trim(self%cfg%projection))
    case ("phaseless")
      self%sample_weight => sample_weight_phaseless
    case ("free")
      self%sample_weight => sample_weight_free
    case ("free_hp")
      self%sample_weight => sample_weight_free_hp
    case default
      self%sample_weight => sample_weight_phaseless
  end select
end subroutine sample_weight_factory


!******************************************************************************** 
!
! Sample weight for the free-projection AFQMC:
!
!    w_sample = w * e^{i*phase}
!
!******************************************************************************** 
function sample_weight_free(self, walker) result(sample_weight)
  class(AfqmcProj), intent(in)  :: self
  type(AfqmcWalker), intent(in) :: walker
  complex(wp)                   :: sample_weight

  sample_weight = walker%weight * exp(im * walker%phase)
end function sample_weight_free


!******************************************************************************** 
!
! Sample weight for the free-projection AFQMC in a halfplane:
!
!    w_sample = w * e^{i*phase}, if phase in the right halfplane
! 
!    w_sample = 0, otherwise
!
!******************************************************************************** 
function sample_weight_free_hp(self, walker) result(sample_weight)
  class(AfqmcProj), intent(in)  :: self
  type(AfqmcWalker), intent(in) :: walker
  complex(wp)                   :: sample_weight

  if (cos(walker%phase) >= self%cfg%min_cos_hp) then
    sample_weight = walker%weight * exp(im * walker%phase)
  else
    sample_weight = zeroc
  end if
end function sample_weight_free_hp


!******************************************************************************** 
!
! Sample weight for the free-projection AFQMC:
!
!    w_sample = w
!
!******************************************************************************** 
function sample_weight_phaseless(self, walker) result(sample_weight)
  class(AfqmcProj), intent(in)  :: self
  type(AfqmcWalker), intent(in) :: walker
  complex(wp)                   :: sample_weight

  sample_weight = walker%weight
end function sample_weight_phaseless


!******************************************************************************** 
!
! Factory for the sample_energy_val method
!
!******************************************************************************** 
subroutine sample_energy_val_factory(self)
  class(AfqmcProj), intent(inout) :: self

  select case (trim(self%cfg%projection))
    case ("phaseless")
      self%sample_energy_val => sample_energy_val_phaseless
    case("free")
      self%sample_energy_val => sample_energy_val_free
    case ("free_hp")
      self%sample_energy_val => sample_energy_val_free
    case default
      self%sample_energy_val => sample_energy_val_phaseless
  end select
end subroutine sample_energy_val_factory


!******************************************************************************** 
!
! Sample energy for the free-projection AFQMC:
!
!    sample_energy = energy
!
!******************************************************************************** 
function sample_energy_val_free(self, energy) result(sample_energy)
  class(AfqmcProj), intent(in) :: self
  complex(wp), intent(in)      :: energy
  complex(wp)                  :: sample_energy

  sample_energy = energy
end function sample_energy_val_free


!******************************************************************************** 
!
! Sample energy for the phaseless AFQMC:
!
!    sample_energy = real(energy)
!
!******************************************************************************** 
function sample_energy_val_phaseless(self, energy) result(sample_energy)
  class(AfqmcProj), intent(in) :: self
  complex(wp), intent(in)      :: energy
  complex(wp)                  :: sample_energy

  sample_energy = real(energy, kind=wp)
end function sample_energy_val_phaseless


!******************************************************************************** 
!
! Factory for the sample_energy_type method
!
!******************************************************************************** 
subroutine sample_energy_type_factory(self)
  class(AfqmcProj), intent(inout) :: self

  select case (trim(self%cfg%projection))
    case ("phaseless")
      self%sample_energy_type => sample_energy_type_phaseless
    case("free")
      self%sample_energy_type => sample_energy_type_free
    case ("free_hp")
      self%sample_energy_type => sample_energy_type_free
    case default
      self%sample_energy_type => sample_energy_type_phaseless
  end select
end subroutine sample_energy_type_factory


!******************************************************************************** 
!
! Sample energy (CEnergy type) for the free-projection AFQMC:
!
!    sample_energy = energy
!
!******************************************************************************** 
function sample_energy_type_free(self, energy) result(sample_energy)
  class(AfqmcProj), intent(in) :: self
  type(CEnergy), intent(in)    :: energy
  type(CEnergy)                :: sample_energy

  sample_energy = energy
end function sample_energy_type_free


!******************************************************************************** 
!
! Sample energy (CEnergy) for the phaseless AFQMC:
!
!    sample_energy = real(energy)
!
!******************************************************************************** 
function sample_energy_type_phaseless(self, energy) result(sample_energy)
  class(AfqmcProj), intent(in) :: self
  type(CEnergy), intent(in)    :: energy
  type(CEnergy)                :: sample_energy

  sample_energy = energy
  sample_energy%e1 = real(energy%e1, kind=wp)
  sample_energy%eh = real(energy%eh, kind=wp)
  sample_energy%ex = real(energy%ex, kind=wp)
end function sample_energy_type_phaseless

end module afqmc_proj