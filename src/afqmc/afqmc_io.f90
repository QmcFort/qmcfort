! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module afqmc_io
!******************************************************************************** 
!       
! IO routines for afqmc module
!
!******************************************************************************** 

use afqmc_importance_sampling
use afqmc_energy
use afqmc_rebalance
use afqmc_rev
use constants
use energy_types
use profiling
use standalone
use statistics
use qmcfort_io

use afqmc_descriptor, only: AfqmcDescriptor
use afqmc_proj, only: AfqmcProj
use afqmc_prop, only: AfqmcProp
use aux_field_gen, only: AuxFieldGen
use growth_est, only: GrowthEst
use mpi, only: comm_world
use wave_trial, only: WaveTrial

implicit none

public

contains

!******************************************************************************** 
!
! Write AFQMC log labels to the files
!
!******************************************************************************** 
subroutine print_afqmc_log_labels()
  write(io%screen%funit,100)      "  tag   ", "  step  ", "   B Re El  ", " B SD El", "   T Re El  ", " T SD El", "   SDM El  ", " B <S>  ", " T <S>  " 
  write(io%screen%funit,100)      "  ---   ", "  ----  ", "   -------  ", " -------", "   -------  ", " -------", "   ------  ", " -----  ", " -----  "

  write(io%qmcfort_log%funit,100) "  tag   ", "  step  ", "   B Re El  ", " B SD El", "   T Re El  ", " T SD El", "   SDM El  ", " B <S>  ", " T <S>  " 
  write(io%qmcfort_log%funit,100) "  ---   ", "  ----  ", "   -------  ", " -------", "   -------  ", " -------", "   ------  ", " -----  ", " -----  "

  write(io%qmcfort_out%funit,101) "   tag  ", "  step  ", "   time   ", " B Re W ", " B Im W ", &
                                  "   B Re El  ", "B Im El ", " B SD El", "   T Re El  ", " T SD El", "   SDM El  ", &
                                  "   B Re Eh  ", "B Im Eh ", " B SD Eh", "   T Re Eh  ", " T SD Eh", "   SDM Eh  ", &
                                  "    Eg   ", " B <S>  ", " T <S>  "
  write(io%qmcfort_out%funit,101) "   ---  ", "  ----  ", "   ----   ", " ------ ", " ------ ", &
                                  "   -------  ", "------- ", " -------", "   -------  ", " -------", "   ------  ", &
                                  "   -------  ", "------- ", " -------", "   -------  ", " -------", "   ------  ", &
                                  "    --   ", " -----  ", " -----  "
  100 format (10(1x,a))
  101 format (20(1x,a))
end subroutine print_afqmc_log_labels


!******************************************************************************** 
!
! Write a table row of the AFQMC log
!
!******************************************************************************** 
subroutine print_afqmc_log_row(af_en, grow)
  type(AfqmcEnergy), intent(in)       :: af_en
  type(GrowthEst), intent(in)         :: grow
  !local
  logical                             :: equil
  real(wp)                            :: tau, block_sign, mean_sign
  real(wp)                            :: block_el_sd, block_eh_sd, mean_el_sd, mean_eh_sd
  complex(wp)                         :: block_weight, block_eh, mean_eh
  character(len=:), allocatable       :: eq_stat
  type(CEnergy)                       :: block_el, mean_el
  type(EnergyStd)                     :: mean_el_sdm
  type(CorrStddev)                    :: mean_eh_sdm

  associate(mc_des => af_en%mc_des)

  tau = mc_des%time_current()
  equil = mc_des%is_eq()

  block_weight = af_en%block_sample_weight()

  block_sign = af_en%block_sign()
  mean_sign = af_en%mean_sign()

  block_el = af_en%block_el()
  block_eh = af_en%block_eh()

  block_el_sd = af_en%block_el_sd()
  block_eh_sd = af_en%block_eh_sd()

  mean_el = af_en%mean_el()
  mean_eh = af_en%mean_eh()

  mean_el_sd = af_en%mean_el_sd()
  mean_eh_sd = af_en%mean_eh_sd()

  mean_el_sdm = af_en%mean_el_sdm(block=mc_des%steps_per_block)
  mean_eh_sdm = af_en%mean_eh_sdm(block=mc_des%steps_per_block)

  block_eh = block_eh + mean_el%enuc
  mean_eh = mean_eh + mean_el%enuc

  if (equil) then
    eq_stat = "afqmc:s:"
  else
    eq_stat = "afqmc:e:"
    mean_el_sdm%e = 0.0_wp
    mean_eh_sdm%stddev = 0.0_wp
  end if

  if (modul(mc_des%block, mc_des%print_screen)) then
    write(io%screen%funit,100)      eq_stat, mc_des%step, real(block_el%total_energy(), wp), block_el_sd, &
                                    real(mean_el%total_energy(), wp), mean_el_sd, mean_el_sdm%e, block_sign, mean_sign

    write(io%qmcfort_log%funit,100) eq_stat, mc_des%step, real(block_el%total_energy(), wp), block_el_sd, &
                                    real(mean_el%total_energy(), wp), mean_el_sd, mean_el_sdm%e, block_sign, mean_sign
  end if

  write(io%qmcfort_out%funit,101) eq_stat, mc_des%step, tau, block_weight, &
                                  block_el%total_energy(), block_el_sd, real(mean_el%total_energy(), wp), mean_el_sd, mean_el_sdm%e, &
                                  block_eh, block_eh_sd, real(mean_eh, wp), mean_eh_sd, mean_eh_sdm%stddev, grow%energy, block_sign, mean_sign

  100 format (1x,a8,1x,i8,1x,f12.6,1x,f8.4,1x,f12.6,1x,f8.4,es11.4,2(f8.4,1x))
  101 format (1x,a8,1x,i8,1x,f10.4,1x,f8.4,1x,f8.4,1x,2(f12.6,1x,f8.4,1x,f8.4,1x,f12.6,1x,f8.4,1x,es11.4,1x),f9.4,2(f8.4,1x))

  end associate
end subroutine print_afqmc_log_row


!******************************************************************************** 
!
! Print info before the main AFQMC loop
!
!******************************************************************************** 
subroutine print_pre_afqmc(af_des, pop, af_en, grow, imp_samp, af_proj, rev, afgen, af_prop, phi_t)
  type(AfqmcDescriptor), intent(in)   :: af_des
  type(pop_controll), intent(in)      :: pop
  type(AfqmcEnergy), intent(in)       :: af_en
  type(GrowthEst), intent(in)         :: grow
  type(ImportanceSampler), intent(in) :: imp_samp
  type(AfqmcProj), intent(in)         :: af_proj
  type(rare_events), intent(in)       :: rev
  type(AuxFieldGen), intent(in)       :: afgen
  type(AfqmcProp), intent(in)         :: af_prop
  class(WaveTrial), intent(inout)     :: phi_t
  !local
  character(len=100)                  :: intfor, logfor, realfor, expfor, charfor

  intfor  = '(1x,t5,a,t50,"= ",10i8)'
  logfor  = '(1x,t5,a,t50,"= ",10l6)'
  realfor = '(1x,t5,a,t50,"= ",10f12.6)'
  expfor = '(1x,t5,a,t50,"= ",10es12.4)'
  charfor = '(1x,t5,a,t50,"= ",a)'

  if (comm_world%mpirank == 0) then
  associate(mc_des => af_des%mc_des)
    write(io%qmcfort_out%funit,*)
    write(io%qmcfort_out%funit,*)      "AFQMC Setup"
    write(io%qmcfort_out%funit,*)      starshort
    write(io%qmcfort_out%funit,*)      "AFQMC descriptor:"
    write(io%qmcfort_out%funit,*)      "afqmc mc controll variables:"
    write(io%qmcfort_out%funit,realfor)"timestep",                         mc_des%tau
    write(io%qmcfort_out%funit,realfor)"total proapgation time",           mc_des%time_tot()
    write(io%qmcfort_out%funit,intfor) "number of blocks",                 mc_des%nblocks
    write(io%qmcfort_out%funit,intfor) "steps per block",                  mc_des%steps_per_block
    write(io%qmcfort_out%funit,intfor) "total number of steps",            mc_des%steps_tot()
    write(io%qmcfort_out%funit,intfor) "print to the screen",              mc_des%print_screen
    write(io%qmcfort_out%funit,intfor) "energy Sampling",                  mc_des%sample
    write(io%qmcfort_out%funit,intfor) "exchange energy Sampling",         mc_des%samplex
    write(io%qmcfort_out%funit,intfor) "reorthogonalization",              mc_des%reorth
    write(io%qmcfort_out%funit,intfor) "population control",               pop%freq
    write(io%qmcfort_out%funit,logfor) "preserved weight",                 pop%lnorm
    write(io%qmcfort_out%funit,realfor)"minimal weight allowed",           pop%min_weight
    write(io%qmcfort_out%funit,realfor)"maximal weight allowed",           pop%max_weight
    write(io%qmcfort_out%funit,intfor) "number of afqmc walkers",          mc_des%nwalkers
    write(io%qmcfort_out%funit,intfor) "number of afqmc walkers per proc", mc_des%walkers_per_rank
  end associate
  end if
    
  call imp_samp%report_setup(comm_world, io%qmcfort_out)
  call afgen%report(comm_world, io%qmcfort_out)
  call af_prop%report(comm_world, io%qmcfort_out%funit)
  if (comm_world%mpirank == 0) call af_proj%report(io%qmcfort_out%funit)
  call phi_t%report(comm_world, io%qmcfort_out)
  call phi_t%report(comm_world, io%screen)

  if (comm_world%mpirank == 0) then
    call rev%report_setup(io%qmcfort_out)

    write(io%screen%funit,*)
    write(io%screen%funit,*)   "AFQMC CALCULATION"
    write(io%screen%funit,*)   starshort
    write(io%screen%funit,*)   "Start AFQMC propagation"

    write(io%qmcfort_log%funit,*)
    write(io%qmcfort_log%funit,*)   "AFQMC CALCULATION"
    write(io%qmcfort_log%funit,*)   starshort
    write(io%qmcfort_log%funit,*)   "Start AFQMC propagation"
    
    write(io%qmcfort_out%funit,*)
    write(io%qmcfort_out%funit,*)   "Start AFQMC propagation"

    call print_afqmc_log_labels()
    call print_afqmc_log_row(af_en, grow)
  end if

end subroutine print_pre_afqmc


!******************************************************************************** 
!
! Print info after the main AFQMC loop
!
!******************************************************************************** 
subroutine print_post_afqmc(pop, af_en, imp_samp, rev)
  type(pop_controll), intent(in)      :: pop
  type(AfqmcEnergy), intent(inout)    :: af_en
  type(ImportanceSampler), intent(in) :: imp_samp
  type(rare_events), intent(inout)    :: rev
  !local variables
  real(wp)                            :: mean_el_sd, mean_eh_sd
  complex(wp)                         :: mean_eh, mean_eh_nw, mean_eh_tnw
  type(CEnergy)                       :: mean_el, mean_el_nw, mean_el_tnw
  type(EnergyStd)                     :: mean_el_sdm, mean_el_sdm_nw, mean_el_sdm_tnw
  type(CorrStddev)                    :: mean_eh_sdm, mean_eh_sdm_nw, mean_eh_sdm_tnw

  mean_el = af_en%mean_el()
  mean_eh = af_en%mean_eh()

  mean_el_nw = af_en%mean_el_nw()
  mean_eh_nw = af_en%mean_eh_nw()

  mean_el_tnw = af_en%mean_el_tnw()
  mean_eh_tnw = af_en%mean_eh_tnw()

  mean_el_sd = af_en%mean_el_sd()
  mean_eh_sd = af_en%mean_eh_sd()

  mean_el_sdm = af_en%mean_el_sdm()
  mean_eh_sdm = af_en%mean_eh_sdm()

  mean_el_sdm_tnw = af_en%mean_el_sdm_tnw()
  mean_eh_sdm_tnw = af_en%mean_eh_sdm_tnw()

  mean_el_sdm_nw = af_en%mean_el_sdm_nw()
  mean_eh_sdm_nw = af_en%mean_eh_sdm_nw()

  mean_eh = mean_eh + mean_el%enuc
  mean_eh_nw = mean_eh_nw + mean_el%enuc

  if (comm_world%mpirank == 0) then
    call mean_el%report(mean_el_sdm, io%screen, "afqmc")
    call mean_el%report(mean_el_sdm, io%qmcfort_log, "afqmc")
    call mean_el%report(mean_el_sdm, io%qmcfort_out, "afqmc")
    call mean_el_tnw%report(mean_el_sdm_tnw, io%qmcfort_out, "afqmc_tnw")
    call mean_el_nw%report(mean_el_sdm_nw, io%qmcfort_out, "afqmc_nw")

    call report_energy_value(comm_world, io%screen, "afqmc hybrid energy", mean_eh, mean_eh_sdm%stddev, mean_eh_sdm%corr_len, "hyb")
    call report_energy_value(comm_world, io%qmcfort_log, "afqmc hybrid energy", mean_eh, mean_eh_sdm%stddev, mean_eh_sdm%corr_len, "hyb")
    call report_energy_value(comm_world, io%qmcfort_out, "afqmc hybrid energy", mean_eh, mean_eh_sdm%stddev, mean_eh_sdm%corr_len, "hyb")
    call report_energy_value(comm_world, io%qmcfort_out, "afqmc time-non-weighted hybrid energy", mean_eh_tnw, mean_eh_sdm_tnw%stddev, mean_eh_sdm_tnw%corr_len, "hyb_tnw")
    call report_energy_value(comm_world, io%qmcfort_out, "afqmc non-weighted hybrid energy", mean_eh_nw, mean_eh_sdm_nw%stddev, mean_eh_sdm_nw%corr_len, "hyb_nw")

    write(io%qmcfort_out%funit,*) 
    write(io%qmcfort_out%funit,*) "*************** Local and Hybrid energy distributions **************"
    write(io%qmcfort_out%funit,106) "SD Eloc           = ", mean_el_sd
    write(io%qmcfort_out%funit,106) "SD Ehyb           = ", mean_eh_sd
    write(io%qmcfort_out%funit,106) "SD Eloc / SD Ehyb = ", mean_el_sd / mean_eh_sd
    write(io%qmcfort_out%funit,107) "Var(Eloc)/<Eloc>  = ", mean_el_sd**2/abs(real(mean_el%e1+mean_el%eh+mean_el%ex, kind=wp))
    write(io%qmcfort_out%funit,107) "Var(Ehyb)/<Ehyb>  = ", mean_eh_sd**2/abs(real(mean_eh-mean_el%enuc-mean_el%e0, kind=wp))

    call rev%report(io%qmcfort_out)
  end if

  call imp_samp%report(comm_world, io%qmcfort_out)
    
  if (comm_world%mpirank == 0) then
    write(io%qmcfort_out%funit,*)
    write(io%qmcfort_out%funit,*) "*************** Population control statistics **************"
    write(io%qmcfort_out%funit,*)
    write(io%qmcfort_out%funit, 101) "number of population controll calls         = ", pop%count_
    write(io%qmcfort_out%funit, 101) "number of walkers sent around               = ", pop%sent_walkers_tot, pop%average()
    write(io%qmcfort_out%funit, 101) "maximal number of clones                    = ", pop%max_walkers_tot

    write(io%screen%funit,*) starlong 
    write(io%qmcfort_log%funit,*) starlong 
    write(io%qmcfort_out%funit,*) starlong 
  end if

  101 format (1x,a,10i8)
  102 format (1x,a,4es14.6)
  103 format (1x,a)
  104 format (1x,a,2f12.6)
  105 format (1x,a,f10.4,a)
  106 format (1x,a,f12.6)
  107 format (1x,a,es11.4)
end subroutine print_post_afqmc


!******************************************************************************** 
!
! Printing at the end of the equilibration phase
!
!******************************************************************************** 
subroutine print_statistics_equilibration_end()
  if (comm_world%mpirank == 0) then
    write(io%screen%funit,*)        "Start AFQMC propagation - Sampling phase:"
    write(io%qmcfort_log%funit,*)   "Start AFQMC propagation - Sampling phase:"
    write(io%qmcfort_out%funit,*)   "Start AFQMC propagation - Sampling phase:"
    call print_afqmc_log_labels()
  end if
end subroutine print_statistics_equilibration_end

end module afqmc_io