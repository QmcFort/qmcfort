! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module afqmc
!******************************************************************************** 
!
!       This Module implements Auxiliary Field Quantum Monte Carlo
!               Method used with Gaussian basis sets
!
!********************************************************************************   

#include "../preproc.inc"
!$ use omp_lib
use afqmc_chkpt
use afqmc_debug
use afqmc_descriptor
use afqmc_energy
use afqmc_io
use afqmc_rebalance
use afqmc_rev
use afqmc_walker
use constants
use energy_types
use hamilton_type
use lapack
use method_base
use mpi
use profiling
use qmcfort_io
use qmcfort_pos
use standalone

use afqmc_importance_sampling, only: ImportanceSampler, ISParameters
use afqmc_intags, only: AFQMCTags
use afqmc_proj, only: AfqmcProjConfig, AfqmcProj
use afqmc_prop, only: AfqmcPropConfig, AfqmcProp
use afqmc_walker_factories, only: afqmc_walker_factory_coeff, afqmc_walker_factory_trial
use aux_field_gen, only: AuxFieldGen, AuxFieldGenDes
use file_handle, only: FileHandle
use growth_est, only: GrowthEst
use growth_est_config, only: GrowthEstConfig
use nosd, only: get_coeff_occ
use wave_trial, only: WaveTrial
use wave_trial_cas, only: WaveTrialCAS
use wave_trial_des, only: WaveTrialDes
use wave_trial_factories, only: wave_trial_factory_afqmc
use wave_trial_noci, only: WaveTrialNOCI
use wave_trial_sd, only: WaveTrialSD

implicit none

private
public :: do_afqmc, setup_afqmc_hamiltonian, &
          calculate_observables, reorthogonalize_walkers, afqmc_propagator

type afqmc_propagator
  type(AfqmcDescriptor)         :: af_des
  type(AfqmcTags)               :: tags
  type(pop_controll)            :: pop
  type(ImportanceSampler)       :: imp_samp
  type(AuxFieldGen)             :: afgen
  type(GrowthEst)               :: grow
  type(rare_events)             :: rev
  type(AfqmcProj)               :: af_proj
  type(AfqmcProp)               :: af_prop
  type(AfqmcEnergy)             :: af_en
  type(afqmc_debugger)          :: af_debug
  type(WaveTrialDes)            :: wtdes
  class(WaveTrial), allocatable :: phi_t

  !Hamiltonian dependent values
  type(structure), pointer      :: struc
  type(Hamilton), pointer       :: ham
  real(wp), allocatable         :: h1tau(:,:,:)
end type afqmc_propagator

contains

!********************************************************************************   
!
! Main routine of the AFQMC module
!
!********************************************************************************   
subroutine do_afqmc(ham, struc)
  type(Hamilton), target, intent(inout) :: ham
  type(structure), intent(in)           :: struc
  !local variables
  integer                               :: i
  logical                               :: lpop
  type(afqmc_propagator)                :: prop
  type(AfqmcWalker), allocatable        :: walkers(:), walkers_(:)
  type(QmcFortMethod)                   :: method
  character(len=*), parameter           :: proc_name = "do_afqmc"
  
  method = QmcFortMethod(method="afqmc", def_active=.false., basis="mo ao_orth", integral_mode="any")
  if (.not. method%is_active()) return
  call ham%meet_method_req(method)

  if (use_profiling) call start_profiling(proc_name)

  call setup_afqmc_hamiltonian(prop, ham, struc)
  call init_afqmc_population(prop, walkers)
  call prop%af_en%update(walkers, prop%af_proj)

  call afqmc_memory_req(prop, walkers)
  call afqmc_runtime(prop%af_des, ham%des)
  call print_pre_afqmc(prop%af_des, prop%pop, prop%af_en, prop%grow, prop%imp_samp, &
                       prop%af_proj, prop%rev, prop%afgen, prop%af_prop, prop%phi_t)
  
  associate(mc_des => prop%af_des%mc_des)
  do while (mc_des%is_running())
    call mc_des%increment()
    call copy_afqmc_walkers(walkers, walkers_)
    call prop%af_prop%propagate_walkers(prop%afgen, prop%imp_samp, walkers)
    call calculate_observables(prop, walkers)
    call accept_reject_walkers(prop, walkers, walkers_)
    call update_afqmc_state(prop, walkers)
    if (modul(mc_des%step, mc_des%reorth)) call reorthogonalize_walkers(prop, walkers)

    !rebalance walkers according to their weights
    call prop%pop%rebalance_walkers(walkers, mc_des%step, lpop)
    if (lpop) call prop%af_en%update(walkers, prop%af_proj)
    if (prop%af_debug%active .and. lpop .and. prop%af_des%comm%mpirank==0) call prop%af_debug%update_rebalance(prop%pop%map)

    call write_afqmc_chkpt(prop%ham%des, mc_des, prop%af_en, walkers, prop%af_des%comm)
  end do
  end associate
 
  if (prop%af_des%dump_walker) call dump_afqmc_wlaker(walkers(1))
  if (prop%af_des%local_energy_test > 0) call test_local_energy_sampling(prop, walkers(1), prop%af_des%local_energy_test)

  call finalize_afqmc(prop)

  if (use_profiling) call end_profiling(proc_name)
end subroutine do_afqmc


!******************************************************************************** 
!
! Print memory requirements for AFQMC
!
!******************************************************************************** 
subroutine afqmc_memory_req(prop, walkers)
  type(afqmc_propagator), intent(in) :: prop
  type(AfqmcWalker), intent(in)      :: walkers(:)
  !local variables
  integer                            :: nwalkers
  real(dp)                           :: temp, total, temp2

  if (prop%af_des%comm%mpirank /= 0) return

  write(io%screen%funit,*)     "Memory reqirements in AFQMC:"
  write(io%screen%funit,*)     "---------------------------------------------"
  write(io%qmcfort_out%funit,*) "Memory reqirements in AFQMC:"
  write(io%qmcfort_out%funit,*) "---------------------------------------------"

  nwalkers = prop%af_des%mc_des%nwalkers
  total = 0.0_dp

  !walkers:
  temp = nwalkers * walkers(1)%sizeof()
  total = total + temp
  write(io%screen%funit,*)      "    Walkers              : ", trim(write_bytes(temp))
  write(io%qmcfort_out%funit,*) "    Walkers              : ", trim(write_bytes(temp))

  !hamiltonian:
  temp = prop%af_des%comm%mpisize * prop%ham%sizeof()
  total = total + temp
  write(io%screen%funit,*)      "    Hamiltonian          : ", trim(write_bytes(temp))
  write(io%qmcfort_out%funit,*) "    Hamiltonian          : ", trim(write_bytes(temp))

  !wave_trial: 
  temp = prop%af_des%comm%mpisize * prop%phi_t%sizeof()
  total = total + temp
  write(io%screen%funit,*)      "    Trial wave function  : ", trim(write_bytes(temp))
  write(io%qmcfort_out%funit,*) "    Trial wave function  : ", trim(write_bytes(temp))

  !AfqmcEnergy:
  temp = prop%af_en%sizeof()
  total = total + temp
  write(io%screen%funit,*)      "    AfqmcEnergy          : ", trim(write_bytes(temp))
  write(io%qmcfort_out%funit,*) "    AfqmcEnergy          : ", trim(write_bytes(temp))

  !propagators:
  temp2 = 0.0_dp
  temp = prop%ham%des%n**2 * nwalkers * prop%ham%des%ispin1 * 16.0_dp !eff_ham
  temp2 = temp2 + temp
  temp = prop%ham%des%n**2 * 2.0_dp * 8.0_dp !one-el hamiltonian
  temp2 = temp2 + temp
  temp = prop%ham%des%ng * nwalkers * (16.0_dp + 8.0_dp) !random fields
  temp2 = temp2 + temp
  total = total + temp2
  write(io%screen%funit,*)      "    Propagators          : ", trim(write_bytes(temp2))
  write(io%qmcfort_out%funit,*) "    Propagators          : ", trim(write_bytes(temp2))


  write(io%screen%funit,*)     "    -----------------------------------------"
  write(io%screen%funit,*)     "    Total                : ", trim(write_bytes(total))
  write(io%screen%funit,*)
  write(io%screen%funit,*)     starshort

  write(io%qmcfort_out%funit,*) "    -----------------------------------------"
  write(io%qmcfort_out%funit,*) "    Total                : ", trim(write_bytes(total))
  write(io%qmcfort_out%funit,*) 
  write(io%qmcfort_out%funit,*) starshort
end subroutine afqmc_memory_req


!******************************************************************************** 
!
! Estimate simulation time in AFQMC
!
!******************************************************************************** 
subroutine afqmc_runtime(af_des, hdes)
  type(AfqmcDescriptor), intent(in)  :: af_des
  type(hamil_descriptor), intent(in) :: hdes
  !local variables
  integer                            :: nsteps, nwalkers
  real(dp)                           :: temp, total_p, total_e, tim

  if (af_des%comm%mpirank /= 0) return 

  total_p = 0.0_dp
  total_e = 0.0_dp

  nsteps = af_des%mc_des%steps_tot()
  nwalkers = af_des%mc_des%nwalkers

  !propagation
  temp = 32.0_dp*hdes%n**2*hdes%nocc*nsteps*nwalkers
  total_p = total_p + temp

  !auxiliary fields
  temp = 4.0_dp*hdes%n**2*hdes%ng*hdes%ispin1*nsteps*nwalkers
  total_p = total_p + temp

  !hartree energy
  temp = 4.0_dp*hdes%ng*hdes%n*hdes%nocc*nsteps*nwalkers
  total_e = total_e + temp

  !exchange energy
  temp = 4.0_dp*hdes%n*hdes%ng*hdes%nocc**2*nsteps*nwalkers / & 
                   (af_des%mc_des%sample*af_des%mc_des%samplex)
  total_e = total_e + temp

  tim = (total_e + total_p) / (10.0_dp**(10) * af_des%comm%totsize)

  write(io%screen%funit,*) "AFQMC runtime estimation:"
  write(io%screen%funit,*) "---------------------------------------------"
  write(io%screen%funit,*) "    FLOP estimate for AFQMC propagation : ", trim(write_flop(total_p))
  write(io%screen%funit,*) "    FLOP estimate for AFQMC energy eval : ", trim(write_flop(total_e))
  write(io%screen%funit,*) "    -----------------------------------------"
  write(io%screen%funit,*) "    Total FLOP estimate in AFQMC        : ", trim(write_flop(total_e+total_p))
  write(io%screen%funit,*)
  write(io%screen%funit,*) "    Average performance assumed         : ", "10 GFLOPS per core" 
  write(io%screen%funit,*) "    Total number of cores               : ", af_des%comm%totsize
  write(io%screen%funit,*) "    Total runtime                       : ", trim(write_time(tim))
  write(io%screen%funit,*) 
  write(io%screen%funit,*) starshort

  write(io%qmcfort_out%funit,*) "AFQMC runtime estimation:"
  write(io%qmcfort_out%funit,*) "---------------------------------------------"
  write(io%qmcfort_out%funit,*) "    FLOP estimate for AFQMC propagation : ", trim(write_flop(total_p))
  write(io%qmcfort_out%funit,*) "    FLOP estimate for AFQMC energy eval : ", trim(write_flop(total_e))
  write(io%qmcfort_out%funit,*) "    -----------------------------------------"
  write(io%qmcfort_out%funit,*) "    Total FLOP estimate in AFQMC        : ", trim(write_flop(total_e+total_p))
  write(io%qmcfort_out%funit,*)
  write(io%qmcfort_out%funit,*) "    average performance assumed         : ", "10 GFLOPS per core" 
  write(io%qmcfort_out%funit,*) "    Total runtime                       : ", trim(write_time(tim))
  write(io%qmcfort_out%funit,*) 
  write(io%qmcfort_out%funit,*) starshort
end subroutine afqmc_runtime


!******************************************************************************** 
!
! Update AFQMC state after propagation
!
!******************************************************************************** 
subroutine update_afqmc_state(prop, walkers)
  type(afqmc_propagator), intent(inout) :: prop
  type(AfqmcWalker), intent(inout)      :: walkers(:)
  !local
  logical                               :: lblock, equil
  real(wp)                              :: last_weight_old, last_weight_new
  type(CEnergy)                         :: mean_energy

  lblock = modul(prop%af_des%mc_des%step, prop%af_des%mc_des%steps_per_block)
  mean_energy = prop%af_en%fast_mean_el()

  last_weight_old = real(prop%af_en%last_sample_weight(), wp)

  call prop%rev%detect(prop%af_des, prop%af_proj, mean_energy, prop%grow, prop%af_des%mc_des%tau, walkers)
  call prop%af_proj%update_weights(walkers)
  call prop%af_en%update(walkers, prop%af_proj)

  last_weight_new = real(prop%af_en%last_sample_weight(), wp)
  call prop%grow%update(last_weight_new/last_weight_old)

  if (prop%af_debug%active) call prop%af_debug%update(walkers)

  if (lblock) then
    if (comm_world%mpirank == 0) call print_afqmc_log_row(prop%af_en, prop%grow)

    if (prop%af_des%mc_des%eq_finish()) then
      call prop%rev%finalize_eq(prop%af_des%mc_des)
      call print_statistics_equilibration_end()
    end if
  end if
end subroutine update_afqmc_state


!******************************************************************************** 
!
! Finalize AFQMC run
!
!******************************************************************************** 
subroutine finalize_afqmc(prop)
  type(afqmc_propagator), intent(inout) :: prop

  call prop%rev%finalize(prop%af_des%mc_des)
  call print_post_afqmc(prop%pop, prop%af_en, prop%imp_samp, prop%rev)

  if (prop%af_debug%active) call prop%af_debug%write_files(comm_world, prop%ham%enuc, prop%ham%h0)
end subroutine finalize_afqmc


!******************************************************************************** 
!
! Setup of the Hamiltonian for afqmc procedure
!
!******************************************************************************** 
subroutine setup_afqmc_hamiltonian(prop, ham, struc)
  type(afqmc_propagator), target, intent(inout) :: prop
  type(Hamilton), target, intent(inout)         :: ham
  type(structure), target, intent(in)           :: struc
  !local variables
  complex(wp)                                   :: e1_trial, e1_init
  complex(wp), allocatable                      :: lgmean_trial(:), eh_trial(:), ex_trial(:), eself_trial(:)
  complex(wp), allocatable                      :: lgmean_init(:), eh_init(:), ex_init(:), eself_init(:)
  type(CEnergy)                                 :: e_trial, e_init
  type(AuxFieldGenDes)                          :: afgen_des
  type(ISParameters)                            :: is_param
  type(AfqmcProjConfig)                         :: af_proj_cfg
  type(AfqmcPropConfig)                         :: af_prop_cfg
  type(GrowthEstConfig)                         :: grow_cfg
  type(AfqmcTags), pointer                      :: tags
  type(McDescriptor), pointer                   :: mc_des
  character(len=*), parameter                   :: proc_name = "setup_afqmc_hamiltonian"

  if (use_profiling) call start_profiling(proc_name)

  !read afqmc related tags
  prop%tags = AfqmcTags()
  tags => prop%tags

  !setup McDescriptor for AFQMC Monte Carlo 
  allocate(mc_des)
  mc_des = McDescriptor(tags%nwalkers, tags%tau, tags%steps_per_block, tags%nblocks, tags%eqblocks, tags%reorth, tags%sample, tags%samplex)

  !setup AfqmcProjConfig
  af_proj_cfg = AfqmcProjConfig()
  call af_proj_cfg%read_tags()

  !setup AfqmcProj
  prop%af_proj = AfqmcProj(af_proj_cfg, tags%tau)

  !setup afqmc_descriptor
  prop%af_des = AfqmcDescriptor(ham%des, mc_des, tags%g_resolved, tags%sp_proj, tags%dump_walker, tags%local_energy_test)

  prop%struc => struc
  prop%ham => ham

  !setup AuxFieldGenDes 
  afgen_des = AuxFieldGenDes(tags%aux_field_dist, tags%nfields_max)

  !setup AuxFieldGen - for creation of the random fields and auxiliary potential
  prop%afgen = AuxFieldGen(afgen_des, ham)

  !setup WaveTrialDes
  prop%wtdes = WaveTrialDes(ham%des, tags%exchange_mode, tags%h_block_size, tags%x_block_size, &
                       tags%compress_x, tags%compress_x_tol)

  !setup trial wave function
  call wave_trial_factory_afqmc(tags%wave_trial_request, prop%wtdes, ham, prop%phi_t)
  call prop%phi_t%setup(lself=.true.)

  !debug: wave_trial
  call prop%phi_t%trial_energy(e1_trial, lgmean_trial, eh_trial, ex_trial, eself_trial)
  call prop%phi_t%initial_energy(e1_init, lgmean_init, eh_init, ex_init, eself_init)
  e_trial = prop%phi_t%pack_energy_gres(e1_trial, eh_trial, ex_trial)
  e_init = prop%phi_t%pack_energy_gres(e1_init, eh_init, ex_init)

  !subtract mean_field
  if (tags%subtract_mean_field) call ham%subtract_mean_field(lgmean_trial)

  !setup importance sampling parameters 
  is_param = ISParameters(lgmean_trial, ex_trial, eself_trial)

  !setup ImportanceSampler object
  prop%imp_samp = ImportanceSampler(tags%is_shift, tags%is_scale, tags%subtract_mean_field, mc_des%isqtau, tags%cut_shift, tags%shift_cutoff, is_param)
  
  !calculate self energy if not calculated in Hamiltonian
  if (.not. allocated(ham%hself)) call ham%self_potential(ham%hself)

  !Setup AfqmcPropConfig object
  af_prop_cfg = AfqmcPropConfig()
  call af_prop_cfg%read_tags()

  !Setup AfqmcProp
  prop%af_prop = AfqmcProp(af_prop_cfg, tags%tau, ham)

  !setup GrowthEstConfig object
  grow_cfg = GrowthEstConfig()
  call grow_cfg%read_tags()

  !setup GrowthEst object
  prop%grow = GrowthEst(grow_cfg, tags%tau)

  call prop%rev%init(prop%ham%des, prop%af_des%mc_des, prop%af_des%mc_des%tau)
  
  !setup AfqmcEnergy object
  prop%af_en = AfqmcEnergy(mc_des, prop%af_des%comm, e_init)

  call prop%pop%init(prop%af_des%mc_des%nwalkers)

  call prop%af_debug%init(prop%af_des%mc_des, prop%af_proj, prop%pop, prop%ham)

  if (use_profiling) call end_profiling(proc_name)
end subroutine setup_afqmc_hamiltonian


!******************************************************************************** 
!
! Initialize walkers and arrays for energy statistics
!
!    If afqmc_chkpt file is present, everything read from the file
!
!    Otherwise population initialized according to the trial state
!
!******************************************************************************** 
subroutine init_afqmc_population(prop, walkers)
  type(afqmc_propagator), intent(inout)       :: prop
  type(AfqmcWalker), allocatable, intent(out) :: walkers(:)
  !local variables
  integer                                     :: nwalkers
  logical                                     :: lexist
  complex(wp), allocatable                    :: coeff(:,:,:)
  
  inquire(file="afqmc_chkpt", exist=lexist)
 
  if (lexist) then
    call read_afqmc_chkpt(prop%ham%des, prop%af_des%mc_des, prop%af_en, walkers, prop%af_des%comm)
  else
    nwalkers = prop%af_des%mc_des%walkers_per_rank
    allocate(coeff(size(prop%ham%phi,1), size(prop%ham%phi,2), size(prop%ham%phi,3)))
    coeff = prop%ham%phi
    call afqmc_walker_factory_coeff(prop%phi_t, coeff, walkers, nwalkers, prop%af_des%sp_proj)
  end if 
end subroutine init_afqmc_population


!******************************************************************************** 
!
!  accept/reject walkers
!
!******************************************************************************** 
subroutine accept_reject_walkers(prop, walkers, walkers_old)
  type(afqmc_propagator), intent(inout) :: prop
  type(AfqmcWalker), intent(inout)      :: walkers(:)
  type(AfqmcWalker), intent(in)         :: walkers_old(:)
  !local variables
  integer                               :: w
  complex(wp)                           :: I1, I2
  real(wp), allocatable                 :: eta(:)
  complex(wp), allocatable              :: shift(:)
  character(len=*), parameter           :: proc_name = "accept_reject_walkers"

  if (use_profiling) call start_profiling(proc_name)

  allocate(eta(size(walkers)))
  call random_number(eta)
  
  do w = 1, size(walkers)
    walkers(w)%acc = 1.0_wp

    if (prop%af_proj%cfg%accept_reject) then
      call move_alloc(walkers(w)%shift, shift)
      call prop%imp_samp%shift(walkers(w))

      I1 = exp(-sum(real(walkers(w)%shift - shift + walkers(w)%random_field, kind=wp)**2/2.0_wp))
      I2 = exp(-sum(real(walkers(w)%random_field, kind=wp)**2/2.0_wp))

      walkers(w)%acc = min(1.0_wp, abs(I1/I2))
    end if

    if (eta(w) <= abs(walkers(w)%acc)) then
      walkers(w)%move = .true.
    else 
      walkers(w) = walkers_old(w)
      walkers(w)%move = .false.
    end if
  end do

  if (use_profiling) call end_profiling(proc_name)
end subroutine accept_reject_walkers


!********************************************************************************
!
! Estimate local energy of the walker using small time-step propagation
!
!    E = lim_{t--->0}  -1/t \frac{<Phi|1 - exp(-t H)|Phi>}{<Phi|Phi>}
!
!********************************************************************************
subroutine local_energy_sampling(prop, walker, walkers_per_rank, energy, stddev)
  type(afqmc_propagator), intent(inout) :: prop
  type(AfqmcWalker), intent(in)         :: walker
  integer, intent(inout)                :: walkers_per_rank
  complex(wp), intent(out)              :: energy
  real(wp), intent(out)                 :: stddev
  !local variables
  integer                               :: nwalkers, walkers_per_block, nblocks
  integer                               :: w, b
  complex(wp)                           :: val, ovlp, ovlp_new, imp
  type(AfqmcWalker), allocatable        :: walkers(:)
  
  walkers_per_block = min(500, walkers_per_rank)
  nblocks = walkers_per_rank / walkers_per_block
  walkers_per_rank = nblocks * walkers_per_block
  nwalkers = walkers_per_rank * comm_world%mpisize

  energy   = zeroc
  stddev = 0.0_wp
  
  call prop%phi_t%get_ovlp(walker%coeff, ovlp)
  allocate(walkers(walkers_per_block))

  do b = 1, nblocks
    do w = 1, walkers_per_block
      walkers(w) = walker
    end do

    call prop%af_prop%propagate_walkers(prop%afgen, prop%imp_samp, walkers)

    do w = 1, walkers_per_block
      call prop%phi_t%get_ovlp(walkers(w)%coeff, ovlp_new)
      val = -(log(ovlp_new / ovlp) + walkers(w)%isf_mf + walkers(w)%isf) / prop%af_des%mc_des%tau
      energy = energy + val
      stddev = stddev + real(val, wp)**2
    end do
  end do

  call comm_world%mpisum(energy, 0)
  call comm_world%mpisum(stddev, 0)

  if (comm_world%mpirank == 0) then
    energy  = energy / real(nwalkers, wp)
    stddev = stddev/real(nwalkers, wp) - real(energy, wp)**2
    stddev = sqrt(stddev / real(nwalkers, wp))
  end if
end subroutine local_energy_sampling


!********************************************************************************
!
! Test local energy evaluation using small time step propagation
!
!********************************************************************************
subroutine test_local_energy_sampling(prop, walker_in, nwalkers)
  type(afqmc_propagator), intent(inout) :: prop
  type(AfqmcWalker), intent(in)         :: walker_in
  integer, intent(inout)                :: nwalkers
  !local variables
  logical                               :: is_walker_file
  real(wp)                              :: stddev
  complex(wp)                           :: energy
  type(AfqmcWalker)                     :: walker
  type(FileHandle)                      :: fh

  inquire(file="afqmc_walker", exist=is_walker_file)

  if (is_walker_file) then
    fh = FileHandle("afqmc_walker")
    call fh%open(status="old", access="stream", form="unformatted", action="read")

    call walker%alloc(prop%ham%des)
    if (comm_world%mpirank == 0) then
      call walker%read_(fh)
      write(*,*) "WALKER READ FROM THE FILE"
    end if

    call fh%close()
  else
    walker = walker_in
  end if
  
  !synchronize walker over all nodes
  call broadcast_afqmc_walker(walker, comm_world, 0)

  call local_energy_sampling(prop, walker, nwalkers, energy, stddev)
  energy = energy + prop%af_en%e_init%enuc

  if (comm_world%mpirank == 0) then
    call print_small_afqmc(io%screen, walker, prop%af_des%mc_des%tau, energy, stddev)
    call print_small_afqmc(io%qmcfort_out, walker, prop%af_des%mc_des%tau, energy, stddev)
  end if

contains 
  subroutine print_small_afqmc(fh, walker, tau,  mean, stddev)
    type(FileHandle), intent(in)  :: fh
    type(AfqmcWalker), intent(in) :: walker
    real(wp), intent(in)          :: tau
    complex(wp), intent(in)       :: mean
    real(wp), intent(in)          :: stddev

    write(fh%funit,*)
    write(fh%funit,100) "Small time-step energy evaluation"
    write(fh%funit,105) "walker energy = ", walker%energy%total_energy()
    write(fh%funit,102) "Number of walkers per MPI rank nwalkers = ", nwalkers
    write(fh%funit,100) "***********************************************"
    write(fh%funit,103) "safqmc:", tau, real(mean, wp), sqrt(tau)*aimag(mean), stddev
    write(fh%funit,*) "*******************************************************************"

    100 format(1x,a)
    101 format(1x,a,f14.8,es14.6)
    102 format(1x,a,i10)
    103 format (1x,a,t10,es12.4,t25,f12.6,t40,f12.6,t55,f12.6)
    104 format(1x,a,es14.6,es14.6)
    105 format(1x,a,f12.6,f12.6)
  end subroutine print_small_afqmc
end subroutine test_local_energy_sampling


!********************************************************************************
!
! Dump one AFQMC walker to the file afqmc_walker
!
! Useful for local_energy_sampling test
!
!********************************************************************************
subroutine dump_afqmc_wlaker(walker)
  type(AfqmcWalker), intent(in) :: walker
  !local variables
  type(FileHandle)              :: fh

  if (comm_world%mpirank == 0) then 
    fh = FileHandle("afqmc_walker")
    call fh%open(status="replace", access="stream", form="unformatted", action="write")
    call walker%write_(fh)
    call fh%close()
  end if
end subroutine dump_afqmc_wlaker


!******************************************************************************** 
!
! Observables calculation in AFQMC
!
!******************************************************************************** 
subroutine calculate_observables(prop, walkers)
  type(afqmc_propagator), intent(inout) :: prop
  type(AfqmcWalker), intent(inout)      :: walkers(:)
  !local variables
  integer                               :: i, w
  logical                               :: is_sample, is_samplex, is_self
  complex(wp), allocatable              :: e1(:), eh(:), ex(:), eself(:)
  complex(wp), allocatable              :: lgmean(:,:), eh_(:,:), ex_(:,:), eself_(:,:)
  complex(wp), allocatable              :: coeff(:,:), ovlp(:)
  type(CEnergy)                         :: en
  type(CEnergy), allocatable            :: walker_energy(:)
  character(len=*), parameter           :: proc_name = "calculate_observables"
  
  if (use_profiling) call start_profiling(proc_name)

  is_sample = prop%af_des%mc_des%is_sample()
  is_samplex = prop%af_des%mc_des%is_samplex()
  is_self = .true.

  call batch_walkers(prop%ham%des, walkers, coeff)

  call prop%phi_t%get_ovlp(coeff, ovlp)

  if (is_sample) then 
    if (prop%af_des%g_resolved) then
      call prop%phi_t%local_energy_gres(coeff, e1, lgmean, eh_, ex_, eself_, is_samplex, is_self, walker_energy)
    else
      call prop%phi_t%local_energy(coeff, e1, lgmean, eh, ex, eself, is_samplex, is_self, walker_energy)
    end if
  else 
    call prop%phi_t%h2_gres_mean(coeff, lgmean, eh_)
  end if

  if (prop%af_des%g_resolved) then
    call unbatch_walkers_gres(walkers, ovlp, lgmean, eh_, ex_, eself_, walker_energy)
  else 
    call unbatch_walkers(walkers, ovlp, lgmean, walker_energy)
  end if

  en = prop%af_en%fast_mean_el()

  !$omp  parallel do private(w) 
  do w  = 1, size(walkers)
    walkers(w)%energyl_old = walkers(w)%energyl 
    walkers(w)%energyl = walkers(w)%energy%electron_energy()
    walkers(w)%energyh_old = walkers(w)%energyh
    walkers(w)%energyh = walkers(w)%hybrid_energy(prop%af_des%mc_des%tau)
    walkers(w)%doverlap = walkers(w)%overlap / walkers(w)%overlap_old * exp(walkers(w)%isf_mf + prop%af_des%mc_des%tau*(real(en%electron_energy(),wp) + prop%grow%energy))
    walkers(w)%dphase = phase(walkers(w)%doverlap)
  end do
  !$omp end parallel do

  if (use_profiling) call end_profiling(proc_name)
end subroutine calculate_observables


!******************************************************************************** 
!
! Dispatch routine for walker reorthogonalization
!
!******************************************************************************** 
subroutine reorthogonalize_walkers(prop, walkers)
  type(afqmc_propagator), intent(in) :: prop
  type(AfqmcWalker), intent(inout)   :: walkers(:)
  !local variables
  integer                            :: w
  character(len=*), parameter        :: proc_name = "reorthogonalize_walkers"

  if (use_profiling) call start_profiling(proc_name)

  !$omp parallel do private(w) 
  do w = 1, size(walkers)
    !call reort_walker_qr(prop%ham%des, walkers(w))
    call reort_walker_chol(prop%ham%des, walkers(w))
  end do
  !$omp end parallel do

  if (use_profiling) call end_profiling(proc_name)
end subroutine reorthogonalize_walkers


!******************************************************************************** 
!
! Reorthogonalize walkers via QR factorization
!
! Actual set of coefficients C = walker%coeff
!
!    C = Q x R
!
!    Q    -    new set of orbitals
!    R    -    det(R) used to update walker overlap       
!
!******************************************************************************** 
subroutine reort_walker_qr(ham_des, walker)
  type(hamil_descriptor), intent(in) :: ham_des
  type(AfqmcWalker), intent(inout)   :: walker
  !local variables
  integer                            :: spin, w1, w2, nocc
  complex(wp)                        :: det
  
  do spin = 1, ham_des%ispin
    nocc = ham_des%nel(spin)
    if (nocc == 0) cycle
    w1 = ham_des%nell(spin) + 1
    w2 = ham_des%nell(spin) + nocc
    call getqr_det(walker%coeff(:,w1:w2), det)
    if (real(det,wp) < 0.0_wp) then
      det = -det
      walker%coeff(:,w1:w2) = walker%coeff(:,w1:w2) * exp(pi*(0.0_wp,1.0_wp)/real(nocc,wp))
    end if
    walker%overlap = walker%overlap / det**ham_des%rspin
  end do
end subroutine reort_walker_qr


!******************************************************************************** 
!
! Reorthogonalize walkers via Cholesky decomposition
!
! Actual set of coefficients C = walker%coeff
!
!    S = C^{H} C 
!    S = L L^{H}
!
!    C --> L^{-1} C
!
!******************************************************************************** 
subroutine reort_walker_chol(ham_des, walker)
  type(hamil_descriptor), intent(in) :: ham_des
  type(AfqmcWalker), intent(inout)   :: walker
  !local variables
  integer                            :: i, spin, w1, w2, nocc, n 
  complex(wp)                        :: det
  complex(wp), allocatable           :: ovlp(:,:)
  
  n = size(walker%coeff, 1)

  do spin = 1, ham_des%ispin
    nocc = ham_des%nel(spin)
    if (nocc == 0) cycle
    w1 = ham_des%nell(spin) + 1
    w2 = ham_des%nell(spin) + nocc
    call get_ovlp_matrix(walker%coeff(:,w1:w2), walker%coeff(:,w1:w2), ovlp)
    call potrf(ovlp, det=det)
    call inverse_tri(ovlp)
    call trmm("r", "u", "n", "n", n, nocc, (1.0_wp,0.0_wp), ovlp, nocc, walker%coeff(:,w1:w2), n)
    walker%overlap = walker%overlap / det**(ham_des%rspin)
  end do
end subroutine reort_walker_chol

end module afqmc