! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module afqmc_prop

#include "../preproc.inc"

use afqmc_prop_config
use constants
use expm_mod
use lapack
use profiling

use afqmc_importance_sampling, only: ImportanceSampler
use afqmc_walker, only: AfqmcWalker
use aux_field_gen, only: AuxFieldGen
use hamilton_type, only: Hamilton
use log_manager_mod, only: logging
use mpi, only: mpi_communicator, comm_world

implicit none

private
public :: AfqmcProp, AfqmcPropConfig, PROP_S0_MODE, PROP_S1_MODE, PROP_S2_MODE

type AfqmcProp
  type(AfqmcPropConfig)           :: cfg
  real(wp)                        :: tau
  type(Hamilton), pointer         :: ham

  real(wp), allocatable           :: prop(:,:,:)
  
  procedure(Isetup_prop), pointer :: setup_prop
  procedure(Ipropagate), pointer  :: propagate
contains
  procedure :: report => report_afqmc_prop
  procedure :: afqmc_prop_factory

  procedure :: propagate_walkers
  procedure :: apply_propagator
  procedure :: propagate_h0
  procedure :: propagate_aux_field

  procedure, private :: setup_prop_s0, setup_prop_s1, setup_prop_s2
  procedure, private :: propagate_s0, propagate_s1, propagate_s2
end type AfqmcProp

interface AfqmcProp
  module procedure finit_afqmc_prop
end interface AfqmcProp

abstract interface
  subroutine Isetup_prop(self)
    import AfqmcProp
    class(AfqmcProp), intent(inout) :: self
  end subroutine Isetup_prop

  subroutine Ipropagate(self, af_gen, imp_samp, walkers)
    import AfqmcProp, ImportanceSampler, AuxFieldGen, AfqmcWalker
    class(AfqmcProp), intent(in)           :: self
    type(AuxFieldGen), intent(in)          :: af_gen
    type(ImportanceSampler), intent(inout) :: imp_samp
    type(AfqmcWalker), intent(inout)       :: walkers(:)
  end subroutine Ipropagate
end interface 

!weight threshold
real(wp), parameter :: weight_tol = 1.0E-06_wp

contains

!******************************************************************************** 
!
! Initialization of the AfqmcProp object
!
!******************************************************************************** 
subroutine init_afqmc_prop(self, cfg, tau, ham)
  type(AfqmcProp), intent(out)       :: self
  type(AfqmcPropConfig), intent(in)  :: cfg
  real(wp), intent(in)               :: tau
  type(Hamilton), target, intent(in) :: ham

  self%cfg = cfg
  self%tau = tau
  self%ham => ham

  call self%afqmc_prop_factory()
  call self%setup_prop()
end subroutine init_afqmc_prop


!******************************************************************************** 
!
! AfqmcProp constructor
!
!******************************************************************************** 
function finit_afqmc_prop(cfg, tau, ham) result(self)
  type(AfqmcPropConfig), intent(in)  :: cfg
  real(wp), intent(in)               :: tau
  type(Hamilton), target, intent(in) :: ham
  type(AfqmcProp)                    :: self

  call init_afqmc_prop(self, cfg, tau, ham)
end function finit_afqmc_prop


!******************************************************************************** 
!
! Report AfqmcProp object
!
!******************************************************************************** 
subroutine report_afqmc_prop(self, comm, funit)
  class(AfqmcProp), intent(in)       :: self
  type(mpi_communicator), intent(in) :: comm
  integer, intent(in)                :: funit
  !local
  character(len=:), allocatable      :: prop_mode_char, expm_mode_char

  prop_mode_char = get_prop_mode(self%cfg%prop_mode)
  expm_mode_char = get_expm_mode(self%cfg%expm_mode)

  if (comm%mpirank == 0) then
    write(funit,*)
    write(funit,100) "AFQMC Propagator (AfqmcProp):"
    write(funit,100) "-----------------------------"
    write(funit,101) "  propagator", prop_mode_char
    write(funit,101) "  expm", expm_mode_char
    write(funit,102) "  expm_order", self%cfg%expm_order
  end if

  100 format (1x,a)
  101 format (1x,a,t50,"= ",a)
  102 format (1x,a,t50,"= ",i6)
end subroutine report_afqmc_prop


!******************************************************************************** 
!
! AfqmcProp Factory:
!
!    Seetup the setup_prop and propagate pointers
!    according to the prop_mode flag
!
!******************************************************************************** 
subroutine afqmc_prop_factory(self)
  class(AfqmcProp), intent(inout) :: self

  select case (self%cfg%prop_mode)
    case (PROP_S0_MODE)
      self%setup_prop => setup_prop_s0
      self%propagate => propagate_s0
    case (PROP_S1_MODE)
      self%setup_prop => setup_prop_s1
      self%propagate => propagate_s1
    case (PROP_S2_MODE)
      self%setup_prop => setup_prop_s2
      self%propagate => propagate_s2
    case default
      call logging%error(.true., "invalid value of the named integer constant prop_mode", __FILE__, __LINE__)
  end select
end subroutine afqmc_prop_factory


!******************************************************************************** 
!
! Precompute propagator parts for the S0 propagator:
!
!    prop = - tau (H1 + Hself)
!
!******************************************************************************** 
subroutine setup_prop_s0(self)
  class(AfqmcProp), intent(inout) :: self

  if (.not. allocated(self%prop)) allocate(self%prop, mold=self%ham%h1)
  self%prop = - self%tau * (self%ham%h1 + self%ham%hself)
end subroutine setup_prop_s0


!******************************************************************************** 
!
! Precompute propagator parts for the S1 propagator:
!
!    prop = exp{-tau (H1 + Hself)}
!
!******************************************************************************** 
subroutine setup_prop_s1(self)
  class(AfqmcProp), intent(inout) :: self
  !local
  integer                         :: spin
  real(wp), allocatable           :: h1(:,:,:)

  allocate(h1, mold=self%ham%h1)
  if (.not. allocated(self%prop)) allocate(self%prop, mold=self%ham%h1)

  h1 = - self%tau * (self%ham%h1 + self%ham%hself)
  do spin = 1, size(h1, 3)
    call expm(oner, h1(:,:,spin) , self%prop(:,:,spin), expm_mode=EXPM_EXACT_MODE)
  end do
end subroutine setup_prop_s1


!******************************************************************************** 
!
! Precompute propagator parts for the S2 propagator:
!
!    prop = exp{-tau (H1 + Hself) / 2}
!
!******************************************************************************** 
subroutine setup_prop_s2(self)
  class(AfqmcProp), intent(inout) :: self
  !local
  integer                         :: spin
  real(wp), allocatable           :: h1(:,:,:)

  allocate(h1, mold=self%ham%h1)
  if (.not. allocated(self%prop)) allocate(self%prop, mold=self%ham%h1)

  h1 = - self%tau * (self%ham%h1 + self%ham%hself) / 2.0_wp
  do spin = 1, size(h1, 3)
    call expm(oner, h1(:,:,spin) , self%prop(:,:,spin), expm_mode=EXPM_EXACT_MODE)
  end do
end subroutine setup_prop_s2


!******************************************************************************** 
!
! S0 propagator
!  
!    e^{-tau H} = e^{-tau (H1 + Haux)}
!
!******************************************************************************** 
subroutine propagate_s0(self, af_gen, imp_samp, walkers)
  class(AfqmcProp), intent(in)           :: self
  type(AuxFieldGen), intent(in)          :: af_gen
  type(ImportanceSampler), intent(inout) :: imp_samp
  type(AfqmcWalker), intent(inout)       :: walkers(:)
  !local
  integer                                :: w
  complex(wp), allocatable               :: aux_pot(:,:,:,:)

  call af_gen%setup_auxiliary_field(imp_samp, walkers, aux_pot)

  !add H1 to Haux
  do w = 1, size(aux_pot, 3)
    aux_pot(:,:,w,:) = aux_pot(:,:,w,:) + self%prop
  end do

  call self%propagate_aux_field(aux_pot, walkers)
  call self%propagate_h0(self%ham%h0+imp_samp%h0mean, walkers)
end subroutine propagate_s0


!******************************************************************************** 
!
! S1 propagator
!  
!    e^{-tau H} = e^{-tau H1} e^{-tau Haux)}
!
!******************************************************************************** 
subroutine propagate_s1(self, af_gen, imp_samp, walkers)
  class(AfqmcProp), intent(in)           :: self
  type(AuxFieldGen), intent(in)          :: af_gen
  type(ImportanceSampler), intent(inout) :: imp_samp
  type(AfqmcWalker), intent(inout)       :: walkers(:)
  !local
  complex(wp), allocatable               :: aux_pot(:,:,:,:)

  call self%apply_propagator(walkers)

  call af_gen%setup_auxiliary_field(imp_samp, walkers, aux_pot)
  call self%propagate_aux_field(aux_pot, walkers)

  call self%propagate_h0(self%ham%h0+imp_samp%h0mean, walkers)
end subroutine propagate_s1


!******************************************************************************** 
!
! S2 propagator
!  
!    e^{-tau H} = e^{-tau H1/2} e^{-tau Haux)} e^{-tau H1/2}
!
!******************************************************************************** 
subroutine propagate_s2(self, af_gen, imp_samp, walkers)
  class(AfqmcProp), intent(in)           :: self
  type(AuxFieldGen), intent(in)          :: af_gen
  type(ImportanceSampler), intent(inout) :: imp_samp
  type(AfqmcWalker), intent(inout)       :: walkers(:)
  !local
  complex(wp), allocatable               :: aux_pot(:,:,:,:)

  call self%apply_propagator(walkers)

  call af_gen%setup_auxiliary_field(imp_samp, walkers, aux_pot)
  call self%propagate_aux_field(aux_pot, walkers)

  call self%apply_propagator(walkers)

  call self%propagate_h0(self%ham%h0+imp_samp%h0mean, walkers)
end subroutine propagate_s2


!******************************************************************************** 
!
! Propagate all walkers 
! 
!    OpenMP parallelization over walker batches
!
! debug: consider moving this to AfqmcState and pack it alltogether so that one 
!        can simply use it as af_state%propagate(walkers)
!
!******************************************************************************** 
subroutine propagate_walkers(self, af_gen, imp_samp, walkers)
  class(AfqmcProp), intent(in)           :: self
  type(AuxFieldGen), intent(in)          :: af_gen
  type(ImportanceSampler), intent(inout) :: imp_samp
  type(AfqmcWalker), intent(inout)       :: walkers(:)
  !local
  integer                                :: walkers_per_thread, walker_groups
  integer                                :: w, w1, w2
  character(len=*), parameter            :: proc_name = "propagate_walkers"
  
  if (use_profiling) call start_profiling(proc_name)

  walkers_per_thread = size(walkers) / comm_world%ompsize
  walker_groups = ceiling(real(size(walkers),wp) / walkers_per_thread)

  !$omp parallel do private(w, w1, w2)
  do w = 1, walker_groups
    w1 = (w-1) * walkers_per_thread + 1
    w2 = min(w*walkers_per_thread, size(walkers))
    call self%propagate(af_gen, imp_samp, walkers(w1:w2))
  end do
  !$omp end parallel do

  if (use_profiling) call end_profiling(proc_name)
end subroutine propagate_walkers


!******************************************************************************** 
!
! Apply precomputed propagator P on the walker orbitals
!
!******************************************************************************** 
subroutine apply_propagator(self, walkers)
  class(AfqmcProp), intent(in)     :: self
  type(AfqmcWalker), intent(inout) :: walkers(:)
  !local
  integer                          :: w, n, nel, i1, i2, spin, sp1
  complex(wp), allocatable         :: coeff(:,:)

  n = self%ham%des%n

  allocate(coeff(n, maxval(self%ham%des%nel)))

  do spin = 1, self%ham%des%ispin
    sp1 = self%ham%des%get_spin1(spin)
    nel = self%ham%des%nel(spin)
    i1 = self%ham%des%nell(spin) + 1
    i2 = self%ham%des%nell(spin) + nel

    do w = 1, size(walkers)
      if (walkers(w)%weight <= weight_tol) cycle

      coeff(:,1:nel) = walkers(w)%coeff(:,i1:i2)
      call gemm("n", "n", n, nel, n, onec, self%prop(:,:,sp1), n, & 
                coeff(:,1:nel), n, zeroc, walkers(w)%coeff(:,i1:i2), n)
    end do
  end do
end subroutine apply_propagator


!******************************************************************************** 
!
! Apply e^{-tau H0 / Nel} on the walker orbitals
!
!******************************************************************************** 
subroutine propagate_h0(self, h0, walkers)
  class(AfqmcProp), intent(in)     :: self
  real(wp), intent(in)             :: h0
  type(AfqmcWalker), intent(inout) :: walkers(:)
  !local
  integer                          :: w
  complex(wp)                      :: mu

  mu = -self%tau * h0 / sum(self%ham%des%nel)
  
  do w = 1, size(walkers)
    if (walkers(w)%weight <= weight_tol) cycle
    walkers(w)%coeff = exp(mu) * walkers(w)%coeff
  end do
end subroutine propagate_h0


!******************************************************************************** 
!
! Propagate random fields
!
!******************************************************************************** 
subroutine propagate_aux_field(self, aux_field, walkers)
  class(AfqmcProp), intent(in)     :: self
  complex(wp), intent(in)          :: aux_field(:,:,:,:)
  type(AfqmcWalker), intent(inout) :: walkers(:)
  !local
  integer                          :: w, i1, i2, spin, sp1

  do spin = 1, self%ham%des%ispin
    sp1 = self%ham%des%get_spin1(spin)
    i1 = self%ham%des%nell(spin) + 1
    i2 = self%ham%des%nell(spin) + self%ham%des%nel(spin)

    do w = 1, size(walkers)
      if (walkers(w)%weight <= weight_tol) cycle
      call expm_act(onec, aux_field(:,:,w,sp1), walkers(w)%coeff(:,i1:i2), & 
                    kmax=self%cfg%expm_order, expm_mode=self%cfg%expm_mode)
    end do
  end do
end subroutine propagate_aux_field

end module afqmc_prop