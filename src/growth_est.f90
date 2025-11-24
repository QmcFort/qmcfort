! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module growth_est

!******************************************************************************** 
!
! GrowthEst calculates the chemical potential for the weight update
!
! Eg_{k} = Eg_{k-1} + damp * \Delta Eg
!   
! \Delta Eg = - log(W_{kp} / W_{(k-1)*p}) / tau
!
! References:
!   G.H. Booth, J. Chem. Phys. 131, 054106 (2009)
!   C. Umrigar, J. Chem. Phys. 99, 2865–2890 (1993)
!   M. Motta, WIREs Comput. Mol. Sci. 8, e1364 (2018) 
!
!******************************************************************************** 

#include "preproc.inc"
use constants
use growth_est_config, only: GrowthEstConfig

type GrowthEst
  type(GrowthEstConfig) :: cfg
  real(wp)              :: tau

  integer               :: step
  real(wp)              :: dweight
  real(wp)              :: energy
contains
  procedure :: update => update_growth_est
end type GrowthEst

interface GrowthEst
  module procedure :: finit_growth_est
end interface GrowthEst

contains

!******************************************************************************** 
! 
! Initialization of the GrowthEst object
!
!******************************************************************************** 
subroutine init_growth_est(self, cfg, tau)
  type(GrowthEst), intent(out)      :: self
  type(GrowthEstConfig), intent(in) :: cfg
  real(wp), intent(in)              :: tau

  self%cfg = cfg
  self%tau = tau

  self%step = 0
  self%dweight = 1.0_wp
  self%energy = 0.0_wp
end subroutine init_growth_est


!******************************************************************************** 
! 
! GrowthEst constructor
!
!******************************************************************************** 
function finit_growth_est(cfg, tau) result(self)
  type(GrowthEstConfig), intent(in) :: cfg
  real(wp), intent(in)              :: tau
  type(GrowthEst)                   :: self

  call init_growth_est(self, cfg, tau)
end function finit_growth_est


!******************************************************************************** 
! 
! Update GrowthEst state and update the growth energy
!
!******************************************************************************** 
subroutine update_growth_est(self, dweight)
  class(GrowthEst), intent(inout) :: self
  real(wp), intent(in)            :: dweight
  !local
  integer                         :: step
  real(wp)                        :: energy_change

  if (.not. self%cfg%active) return

  self%step = self%step + 1 

  !update the weight ratio
  self%dweight = self%dweight * dweight

  if (mod(self%step, self%cfg%period) /= 0) return

  !claculate change in chemical potential
  energy_change = - log(self%dweight) / (self%cfg%period * self%tau)

  !reset the weight ratio after evaluating energy change
  self%dweight = 1.0_wp

  !update the growth energy
  self%energy = self%energy + self%cfg%damp * energy_change
end subroutine update_growth_est

end module growth_est