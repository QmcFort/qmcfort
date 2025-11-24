! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module test_afqmc_proj

use afqmc_walker

use afqmc_proj, only: AfqmcProjConfig, AfqmcProj
use constants, only: wp
use fruit, only: assert_equals, run_test_case, test_set_summary

implicit none

private
public :: test_afqmc_proj_set

contains

subroutine test_afqmc_proj_set
  call run_test_case(test_sample_weight_free, "get sample_weight for the AfqmcWalker in free regime")
  call run_test_case(test_update_weight_ph, "update weight for the phaseless projection")
  call run_test_case(test_reweighting_factor, "test reweighting_factor + update_weight together")

  call test_set_summary("src/afqmc/afqmc_proj.f90")
end subroutine test_afqmc_proj_set

subroutine test_update_weight_ph
  real(wp), parameter   :: tol = 1.0e-04_wp
  type(AfqmcProjConfig) :: af_proj_cfg
  type(AfqmcProj)       :: af_proj
  type(AfqmcWalker)     :: walker
  real(wp)              :: expected, result_

  af_proj_cfg = AfqmcProjConfig()
  af_proj = AFqmcProj(af_proj_cfg, tau=0.01_wp)
                        
  walker%weight = 1.0_wp
  walker%dphase = pi/6.0_wp
  walker%rew_w = 0.8_wp
  walker%rew_p = pi/12.0_wp
  walker%move = .true.

  expected = walker%weight * walker%rew_w * cos(walker%dphase)
  call af_proj%update_weight(walker)
  result_ = walker%weight

  call assert_equals(expected, result_, tol)
end subroutine test_update_weight_ph

subroutine test_sample_weight_free
  real(wp), parameter   :: tol = 1.0e-04_wp
  type(AfqmcProjConfig) :: af_proj_cfg
  type(AfqmcProj)       :: af_proj
  type(AfqmcWalker)     :: walker
  complex(wp)           :: expected, result_

  af_proj_cfg = AfqmcProjConfig(projection="free")
  af_proj = AFqmcProj(af_proj_cfg, tau=0.01_wp)
                        
  walker%weight = 1.0_wp
  walker%phase = pi/6.0_wp

  expected = walker%weight * exp((0.0_wp,1.0_wp)*walker%phase)
  result_ = af_proj%sample_weight(walker)

  call assert_equals(expected, result_, tol)
end subroutine test_sample_weight_free

subroutine test_reweighting_factor
  real(wp), parameter   :: tol = 1.0e-04_wp, hybrid=0.5_wp
  type(AfqmcProjConfig) :: af_proj_cfg
  type(AfqmcProj)       :: af_proj
  type(AfqmcWalker)     :: walker
  real(wp)              :: tau, mu
  real(wp)              :: expected_weight, result_weight
  real(wp)              :: expected_phase, result_phase

  tau = 0.01_wp
  mu = -100.0_wp

  walker%weight = 1.0_wp
  walker%acc = 1.0_wp
  walker%phase = 0.0_wp
  walker%energyh = (-100.0_wp, -0.1_wp)
  walker%energyl = -101.0_wp
  walker%move = .true.

  af_proj_cfg = AfqmcProjConfig(projection="free", hybrid=hybrid)
  af_proj = AfqmcProj(af_proj_cfg, tau)

  expected_weight = walker%weight * exp(-tau*(hybrid*real(walker%energyh+walker%energyl,wp)-mu))
  expected_phase = walker%phase - tau*aimag(hybrid*walker%energyh)

  call af_proj%reweighting_factor(walker, tau, mu)
  call af_proj%update_weight(walker)

  result_weight = walker%weight
  result_phase = walker%phase

  call assert_equals(expected_weight, result_weight, tol)
  call assert_equals(expected_phase, result_phase, tol)
end subroutine test_reweighting_factor

end module test_afqmc_proj