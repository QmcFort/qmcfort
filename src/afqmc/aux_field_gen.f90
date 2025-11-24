! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module aux_field_gen

#include "../preproc.inc"
use constants
use lapack
use mpi
use profiling
use afqmc_walker
use afqmc_importance_sampling,  only: ImportanceSampler
use drng_factories, only: drng_glb
use file_handle, only: FileHandle
use hamilton_type, only: Hamilton

implicit none

private
public :: AuxFieldGen, AuxFieldGenDes

type AuxFieldGenDes
  character(len=:), allocatable :: aux_field_dist
  integer                       :: nfields_max
end type AuxFieldGenDes

type AuxFieldGen
  type(AuxFieldGenDes)                    :: des
  type(Hamilton), pointer                 :: ham
  procedure(Isetup_random_field), pointer :: setup_random_field => null()
contains
  procedure          :: report => report_aux_field_gen
  procedure          :: setup_aux_pot
  procedure          :: setup_auxiliary_field
  procedure, private :: setup_random_field_factory
end type AuxFieldGen

interface AuxFieldGenDes
  module procedure finit_aux_field_gen_des
end interface AuxFieldGenDes

interface AuxFieldGen
  module procedure finit_aux_field_gen
end interface AuxFieldGen

abstract interface
  subroutine Isetup_random_field(self, random_field, nw)
    import AuxFieldGen, wp
    class(AuxFieldGen), intent(in)     :: self
    real(wp), allocatable, intent(out) :: random_field(:,:)
    integer, intent(in)                :: nw
  end subroutine Isetup_random_field
end interface
  
contains

!******************************************************************************** 
!
! Initialization of the AuxFieldGenDes object
!
!******************************************************************************** 
subroutine init_aux_field_gen_des(aux_field_dist, nfields_max, self)
  character(len=*), intent(in)      :: aux_field_dist
  integer, intent(in)               :: nfields_max
  type(AuxFieldGenDes), intent(out) :: self

  self%aux_field_dist = aux_field_dist
  self%nfields_max = nfields_max
end subroutine init_aux_field_gen_des


!******************************************************************************** 
!
! AuxFieldGenDes constructor
!
!******************************************************************************** 
function finit_aux_field_gen_des(aux_field_dist, nfields_max) result(self)
  character(len=*), intent(in) :: aux_field_dist
  integer, intent(in)          :: nfields_max
  type(AuxFieldGenDes)         :: self

  call init_aux_field_gen_des(aux_field_dist, nfields_max, self)
end function finit_aux_field_gen_des


!******************************************************************************** 
!
! Initialization of the AuxFieldGen object
!
!******************************************************************************** 
subroutine init_aux_field_gen(des, ham, self)
  type(AuxFieldGenDes), intent(in)   :: des
  type(Hamilton), target, intent(in) :: ham
  type(AuxFieldGen), intent(out)     :: self

  self%des = des
  self%ham => ham

  !check whether nfields_max is setup properly
  if (self%des%nfields_max < 0) self%des%nfields_max = self%ham%des%ng

  call self%setup_random_field_factory()
end subroutine init_aux_field_gen


!******************************************************************************** 
!
! AuxFieldGen constructor
!
!******************************************************************************** 
function finit_aux_field_gen(des, ham) result(self)
  type(AuxFieldGenDes), intent(in)   :: des
  type(Hamilton), target, intent(in) :: ham
  type(AuxFieldGen)                  :: self

  call init_aux_field_gen(des, ham, self)
end function finit_aux_field_gen


!******************************************************************************** 
!
! Report AuxFieldGen object
!
!******************************************************************************** 
subroutine report_aux_field_gen(self, comm, fh)
  class(AuxFieldGen), intent(in)     :: self
  type(mpi_communicator), intent(in) :: comm
  type(FileHandle), intent(in)       :: fh

  if (comm%mpirank == 0) then
    write(fh%funit,*)
    write(fh%funit,100) "Random fields / Auxiliary potential setup (AuxFieldGen):"
    write(fh%funit,100) "--------------------------------------------------------"
    write(fh%funit,101) "  distribution of the random fields", self%des%aux_field_dist
    write(fh%funit,102) "  number of fields for the propagation", self%des%nfields_max
  end if

  100 format (1x,a)
  101 format (1x,a,t50,"= ",a)
  102 format (1x,a,t50,"= ",i8)
end subroutine report_aux_field_gen


!******************************************************************************** 
!
! Associate setup_random_field routine
!
!******************************************************************************** 
subroutine setup_random_field_factory(self)
  class(AuxFieldGen), intent(inout) :: self

  select case (self%des%aux_field_dist)
    case ('normal')
      self%setup_random_field => setup_random_field_normal
    case ('uniform')
      self%setup_random_field => setup_random_field_uniform
    case ('none')
      self%setup_random_field => setup_random_field_none
    case default
      self%des%aux_field_dist = "normal"
      self%setup_random_field => setup_random_field_normal
  end select
end subroutine setup_random_field_factory


!******************************************************************************** 
!
! Setup normally distributed random field
!
!******************************************************************************** 
subroutine setup_random_field_normal(self, random_field, nw)
  class(AuxFieldGen), intent(in)     :: self
  real(wp), allocatable, intent(out) :: random_field(:,:)
  integer, intent(in)                :: nw
  !local 
  integer                            :: ng
  
  ng = self%des%nfields_max
  allocate(random_field(ng, nw))
  call drng_glb%normal(random_field)
end subroutine setup_random_field_normal


!******************************************************************************** 
!
! Setup uniformly distributed random field
!
!******************************************************************************** 
subroutine setup_random_field_uniform(self, random_field, nw)
  class(AuxFieldGen), intent(in)     :: self
  real(wp), allocatable, intent(out) :: random_field(:,:)
  integer, intent(in)                :: nw
  !local
  integer                            :: ng

  ng = self%des%nfields_max
  allocate(random_field(ng, nw))
  call drng_glb%uniform(random_field, -sqrt(3.0_wp), sqrt(3.0_wp))
end subroutine setup_random_field_uniform


!******************************************************************************** 
!
! Zero out random fields
!
!******************************************************************************** 
subroutine setup_random_field_none(self, random_field, nw)
  class(AuxFieldGen), intent(in)     :: self
  real(wp), allocatable, intent(out) :: random_field(:,:)
  integer, intent(in)                :: nw
  !local
  integer                            :: ng

  ng = self%des%nfields_max
  allocate(random_field(ng, nw))
  random_field = 0.0_wp
end subroutine setup_random_field_none


!******************************************************************************** 
!
! Setup random field and calculate auxiliary potential
!
!******************************************************************************** 
subroutine setup_auxiliary_field(self, imp_samp, walkers, aux_pot)
  class(AuxFieldGen), intent(in)         :: self
  type(ImportanceSampler), intent(inout) :: imp_samp
  type(AfqmcWalker), intent(inout)       :: walkers(:)
  complex(wp), allocatable, intent(out)  :: aux_pot(:,:,:,:)
  !local
  integer                                :: w
  real(wp), allocatable                  :: random_field(:,:)
  complex(wp), allocatable               :: shifted_field(:,:) 
  character(len=*), parameter            :: proc_name = "setup_auxiliary_field"

  if (use_profiling) call start_profiling(proc_name)

  call self%setup_random_field(random_field, size(walkers))
  allocate(shifted_field(size(random_field,1), size(random_field,2)))

  do w = 1, size(walkers)
    if (.not. allocated(walkers(w)%random_field)) allocate(walkers(w)%random_field(size(random_field,1)))
    walkers(w)%random_field = random_field(:,w)
    call imp_samp%move(walkers(w)) !walker%shifted_field setup here
    walkers(w)%isf = imp_samp%weight(walkers(w))
    walkers(w)%isf_mf = imp_samp%mean_field_weight(walkers(w))
    shifted_field(:,w) = imp_samp%isqtau * walkers(w)%shifted_field
  end do

  call self%setup_aux_pot(shifted_field, aux_pot) 

  if (use_profiling) call end_profiling(proc_name)
end subroutine setup_auxiliary_field




!******************************************************************************** 
!
! Create auxiliary potential: convolute random fields with Cholesky vectors
!
! Calculate:
!
!    H_{pq}^{w} = \sum_g \bar{x}_{g}^{w} L_{g,pq}
!    \bar{x}_{g}^{w} = i sqrt(dt) (x_{g}^{w} + f_{g}^{w})
!    f_{g}^{w} = i sqrt(dt) <L_g - l_g>
!
!******************************************************************************** 
subroutine setup_aux_pot(self, shifted_field, aux_pot)
  class(AuxFieldGen), intent(in)                :: self
  complex(wp), intent(in)                       :: shifted_field(:,:)
  complex(wp), allocatable, target, intent(out) :: aux_pot(:,:,:,:)
  !local variables
  integer                                       :: n, nn, ng, nw, ispin, spin
  complex(wp), contiguous, pointer              :: aux_pot_spin(:,:)
  character(len=*), parameter                   :: proc_name = "setup_aux_pot"
  character(len=*), parameter                   :: timer_name = "gemm_fields"
  
  if (use_profiling) call start_profiling(proc_name)

  n = self%ham%des%n
  nn = n * n
  ng = size(shifted_field, 1)
  nw = size(shifted_field, 2)
  ispin = size(self%ham%h2_gres, 3)
  
  allocate(aux_pot(n, n, nw, ispin))

  do spin = 1, ispin
    aux_pot_spin(1:nn, 1:nw) => aux_pot(:,:,:,spin)

    if (profile_code) call start_profiling(timer_name)
    call gemm("n", "n", nn, nw, ng, onec, self%ham%h2_gres(:,1:ng,spin), nn, shifted_field, ng, zeroc, &
              aux_pot_spin, nn)
    if (profile_code) call end_profiling(timer_name, nn, nw, ng, "rc")
  end do

  if (use_profiling) call end_profiling(proc_name)
end subroutine setup_aux_pot

end module aux_field_gen