! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module afqmc_importance_sampling

#include "../preproc.inc"
use constants
use mpi
use file_handle, only: FileHandle
use afqmc_walker

implicit none

private
public ImportanceSampler, ISParameters

type ISParameters
  complex(wp), allocatable :: lgmean(:)
  complex(wp), allocatable :: ex(:)
  complex(wp), allocatable :: eself(:)
end type ISParameters

type ImportanceSampler
  character(len=:), allocatable  :: is_shift
  character(len=:), allocatable  :: is_scale
  logical                        :: subtract_mean_field
  complex(wp)                    :: isqtau
  character(len=:), allocatable  :: i_cut_shift
  character(len=:), allocatable  :: i_shift_cutoff
  type(ISParameters)             :: is_param

  complex(wp), allocatable       :: mean_field_shift(:)
  real(wp)                       :: h0mean

  real(wp)                       :: shift_cutoff
  integer, allocatable           :: rare_events(:)
  complex(wp)                    :: tau

  procedure(Ishift), pointer     :: shift
  procedure(Iscale), pointer     :: scale
  procedure(Icut_shift), pointer :: cut_shift
contains
  procedure          :: move => is_move
  procedure          :: weight => is_weight
  procedure          :: mean_field_weight => is_mean_field_weight
  procedure          :: report => is_report
  procedure          :: report_setup => is_report_setup

  procedure, private :: shift_factory => is_shift_factory
  procedure, private :: scale_factory => is_scale_factory

  procedure, private :: is_shift_none, is_shift_trial, is_shift_mix
  procedure, private :: is_scale_none, is_scale_trial, is_scale_mix

  procedure, private :: shift_cutoff_factory => is_shift_cutoff_factory
  procedure, private :: cut_shift_factory => is_cut_shift_factory

  procedure, private :: is_cut_shift_none, is_cut_shift_cap, is_cut_shift_rescale
  procedure, private :: is_cut_shift_remove, is_cut_shift_remove_field, is_cut_shift_unr
end type ImportanceSampler

!debug:
real(wp) :: max_scale = 0.0_wp
real(wp) :: min_scale = 1.0_wp
real(wp) :: max_shift = 0.0_wp
real(wp) :: min_shift = 1.0_wp

interface ISParameters
  module procedure init_is_parameters
end interface ISParameters

interface ImportanceSampler
  module procedure init_importance_sampler
end interface ImportanceSampler

abstract interface
  subroutine Ishift(self, walker)
    import ImportanceSampler, AfqmcWalker
    class(ImportanceSampler), intent(inout) :: self
    type(AfqmcWalker), intent(inout)        :: walker 
  end subroutine Ishift

  subroutine Iscale(self, walker)
    import ImportanceSampler, AfqmcWalker
    class(ImportanceSampler), intent(inout) :: self
    type(AfqmcWalker), intent(inout)        :: walker
  end subroutine Iscale

  subroutine Icut_shift(self, walker)
    import ImportanceSampler, AfqmcWalker
    class(ImportanceSampler), intent(inout) :: self
    type(AfqmcWalker), intent(inout)        :: walker
  end subroutine Icut_shift
end interface

contains 

!******************************************************************************** 
!
! ISParameters Constructor
!
! Store trial/initial values of lgmean, ex, eself 
!
!******************************************************************************** 
function init_is_parameters(lgmean, ex, eself) result(self)
  complex(wp), allocatable, intent(in) :: lgmean(:)
  complex(wp), allocatable, intent(in) :: ex(:)
  complex(wp), allocatable, intent(in) :: eself(:)
  type(ISParameters)                   :: self

  allocate(self%lgmean, source=lgmean)
  allocate(self%ex, source=ex)
  allocate(self%eself, source=eself)
end function init_is_parameters


!******************************************************************************** 
!
! ImportanceSampler Constructor
!
! Stores information and provides routines for mean-field subtraction and 
! importance sampling in AFQMC
!
!******************************************************************************** 
function init_importance_sampler(is_shift, is_scale, subtract_mean_field, isqtau, i_cut_shift, i_shift_cutoff, is_param) result(self)
  character(len=*), intent(in)   :: is_shift
  character(len=*), intent(in)   :: is_scale
  logical, intent(in)            :: subtract_mean_field
  complex(wp), intent(in)        :: isqtau
  character(len=*), intent(in)   :: i_cut_shift
  character(len=*), intent(in)   :: i_shift_cutoff
  type(ISParameters), intent(in) :: is_param
  type(ImportanceSampler)        :: self

  self%is_shift = is_shift
  self%is_scale = is_scale
  self%subtract_mean_field = subtract_mean_field
  self%isqtau = isqtau
  self%i_cut_shift = i_cut_shift
  self%i_shift_cutoff = i_shift_cutoff
  self%is_param = is_param
  
  allocate(self%mean_field_shift, source=is_param%lgmean)
  if (.not. subtract_mean_field) self%mean_field_shift = zeroc
  self%h0mean = -0.5_wp * sum(abs(self%mean_field_shift)**2)

  self%tau = - self%isqtau**2

  allocate(self%rare_events(size(self%mean_field_shift)))
  self%rare_events = 0

  call self%shift_factory()
  call self%scale_factory()
  call self%shift_cutoff_factory()
  call self%cut_shift_factory()
end function init_importance_sampler


!******************************************************************************** 
!
! Monte Carlo move in AFQMC
!
!    q_g --> q_g + shift_g + scale_g x_g 
!
!    shift_g = i sqrt(tau) (l_g - l_{g, mf})
!    scale_g = 1 / sqrt(1 + 2tau(Ex_g - Eself_g))
!
!******************************************************************************** 
subroutine is_move(self, walker)
  class(ImportanceSampler), intent(inout) :: self
  type(AfqmcWalker), intent(inout)        :: walker

  call self%shift(walker)
  call self%scale(walker)
  if (.not. allocated(walker%shifted_field)) allocate(walker%shifted_field(size(walker%random_field)))
  walker%shifted_field = walker%random_field / sqrt(onec + walker%scale) + walker%shift

!debug:
  max_scale = max(max_scale, maxval(abs(onec/sqrt(onec+walker%scale))))
  min_scale = min(min_scale, minval(abs(onec/sqrt(onec+walker%scale))))

  max_shift = max(max_shift, maxval(abs(walker%shift)))
  min_shift = min(min_shift, minval(abs(walker%shift)))
end subroutine is_move


!******************************************************************************** 
!
! Importance sampling weight:
!
!    w(x) = p_target(x) / p_prop(x)
!
!******************************************************************************** 
function is_weight(self, walker) result(weight)
  class(ImportanceSampler), intent(in) :: self
  type(AfqmcWalker), intent(inout)     :: walker
  complex(wp)                          :: weight

  weight = 0.5_wp * sum((onec + walker%scale)*(walker%shifted_field-walker%shift)**2) &
         - 0.5_wp * sum(walker%shifted_field**2) &
         - 0.5_wp * sum(log(onec + walker%scale))
!FOR DETBAL
!  weight = 0.5_wp * sum((onec + walker%scale)*abs(walker%shifted_field-walker%shift)**2) &
!         - 0.5_wp * sum(abs(walker%shifted_field)**2) &
!         - 0.5_wp * sum(log(onec + walker%scale))
end function is_weight


!******************************************************************************** 
!
! Importane sampling weight from mean-field subtraction
!
!    w(x) = - i sqrt(tau) * sum((scale_g x_g + shift_g) * l_{g, mf})
!
!******************************************************************************** 
function is_mean_field_weight(self, walker) result(weight)
  class(ImportanceSampler), intent(in) :: self
  type(AfqmcWalker), intent(inout)     :: walker
  complex(wp)                          :: weight
  !local
  integer                              :: ng

  ng = size(walker%shifted_field)

  weight = - self%isqtau * sum(walker%shifted_field * self%mean_field_shift(1:ng))
end function is_mean_field_weight


!******************************************************************************** 
!
! Report setup of the ImportanceSampler object
!
!******************************************************************************** 
subroutine is_report_setup(self, comm, fh)
  class(ImportanceSampler), intent(in) :: self
  type(mpi_communicator), intent(in)   :: comm
  type(FileHandle), intent(in)         :: fh

  if (comm%mpirank == 0) then
    write(fh%funit,*)
    write(fh%funit,100) "Importance sampling setup:"
    write(fh%funit,100) "--------------------------"
    write(fh%funit,101) "  subtract mean field", self%subtract_mean_field
    write(fh%funit,102) "  importance sampling shift", self%is_shift
    write(fh%funit,102) "  importance sampling scale", self%is_scale
    write(fh%funit,102) "  importance sampling cut shift ", self%i_cut_shift
    write(fh%funit,103) "  importance sampling shift cutoff ", self%shift_cutoff
  end if

  100 format (1x,a)
  101 format (1x,a,t50,"= ",l)
  102 format (1x,a,t50,"= ",a)
  103 format (1x,a,t50,"= ",f12.4)
end subroutine is_report_setup


!******************************************************************************** 
!
! Report ImportanceSampler object
!
!******************************************************************************** 
subroutine is_report(self, comm, fh)
  class(ImportanceSampler), intent(in) :: self
  type(mpi_communicator), intent(in)   :: comm
  type(FileHandle), intent(in)         :: fh
  !local variables
  integer                              :: i, g
  integer, allocatable                 :: rare_events(:)

  allocate(rare_events, source=self%rare_events)
  call comm%mpisum(rare_events)

!debug:
call comm%reduce(max_scale, op=mpi_max)
call comm%reduce(min_scale, op=mpi_min)
call comm%reduce(max_shift, op=mpi_max)
call comm%reduce(min_shift, op=mpi_min)

  if (comm%mpirank == 0) then
    write(fh%funit,100) "Cut shift mode                    = ", self%i_cut_shift  
    write(fh%funit,103) "Shift cutoff                      = ", self%shift_cutoff
    write(fh%funit,101) "  total number of rare events     = ", sum(rare_events)
    do i = 1, min(size(rare_events), 10)
      g = maxloc(rare_events, 1)
      if (rare_events(g) == 0) exit
      write(fh%funit, 102) "    ", i, " g index = ", g, " nrev = ", rare_events(g)
      rare_events(g) = 0
    end do
!debug:
write(*,*) "max_scale = ", max_scale
write(*,*) "min_scale = ", min_scale
write(*,*)
write(*,*) "max_shift = ", max_shift
write(*,*) "min_shift = ", min_shift
  end if

  100 format (1x,a,a)
  101 format (1x,a,i10)
  102 format (1x,a,i6,a,i10,a,i10)
  103 format (1x,a,f12.4)
end subroutine is_report


!******************************************************************************** 
!
! Factory to select routine to calculatte importance sampling shift
!
!******************************************************************************** 
subroutine is_shift_factory(self)
  class(ImportanceSampler), intent(inout) :: self

  select case (trim(self%is_shift))
    case ("none")
      self%shift => is_shift_none
    case("trial")
      self%shift => is_shift_trial
    case("mix")
      self%shift => is_shift_mix
    case default
      self%shift => is_shift_mix
  end select
end subroutine is_shift_factory


!******************************************************************************** 
!
! Importance sampling shift for is_shift = "none"
!
!******************************************************************************** 
subroutine is_shift_none(self, walker)
  class(ImportanceSampler), intent(inout) :: self
  type(AfqmcWalker), intent(inout)        :: walker
  !local
  integer                                 :: ng

  ng = size(walker%random_field)

  if (.not. allocated(walker%shift)) allocate(walker%shift(ng))
  walker%shift = zeroc
end subroutine is_shift_none


!******************************************************************************** 
!
! Importance sampling shift for is_shift = "trial"
!
!******************************************************************************** 
subroutine is_shift_trial(self, walker)
  class(ImportanceSampler), intent(inout) :: self
  type(AfqmcWalker), intent(inout)        :: walker
  !local
  integer                                 :: ng

  ng = size(walker%random_field)

  if (.not. allocated(walker%shift)) allocate(walker%shift(ng))
  walker%shift = self%isqtau * (self%is_param%lgmean(1:ng) - self%mean_field_shift(1:ng))
end subroutine is_shift_trial


!******************************************************************************** 
!
! Importance sampling shift for is_shift = "mix"
!
!******************************************************************************** 
subroutine is_shift_mix(self, walker)
  class(ImportanceSampler), intent(inout) :: self
  type(AfqmcWalker), intent(inout)        :: walker
  !local
  integer                                 :: ng

  ng = size(walker%random_field)

  if (.not. allocated(walker%shift)) allocate(walker%shift(ng))
  walker%shift = self%isqtau * (walker%lgmean(1:ng) - self%mean_field_shift(1:ng))
  call self%cut_shift(walker)
end subroutine is_shift_mix


!******************************************************************************** 
!
! Factory to select routine to calculatte importance sampling scale
!
!******************************************************************************** 
subroutine is_scale_factory(self)
  class(ImportanceSampler), intent(inout) :: self

  select case (trim(self%is_scale))
    case("none")
      self%scale => is_scale_none
    case("trial")
      self%scale => is_scale_trial
    case("mix")
      self%scale => is_scale_mix
    case default
      self%scale => is_scale_none
  end select
end subroutine is_scale_factory


!******************************************************************************** 
!
! Importance sampling scale for is_scale = "none"
!
!******************************************************************************** 
subroutine is_scale_none(self, walker)
  class(ImportanceSampler), intent(inout) :: self
  type(AfqmcWalker), intent(inout)        :: walker
  !local
  integer                                 :: ng
  
  ng = size(walker%random_field)

  if (.not. allocated(walker%scale)) allocate(walker%scale(ng))
  walker%scale = zeroc
end subroutine is_scale_none


!******************************************************************************** 
!
! Importance sampling scale for is_scale = "trial"
!
!******************************************************************************** 
subroutine is_scale_trial(self, walker)
  class(ImportanceSampler), intent(inout) :: self
  type(AfqmcWalker), intent(inout)        :: walker
  !local
  integer                                 :: ng
  
  ng = size(walker%random_field)

  if (.not. allocated(walker%scale)) allocate(walker%scale(ng))
  walker%scale = 2.0_wp * self%tau * (self%is_param%ex(1:ng) - self%is_param%eself(1:ng))
end subroutine is_scale_trial


!******************************************************************************** 
!
! Importance sampling scale for is_scale = "mix"
!
!******************************************************************************** 
subroutine is_scale_mix(self, walker)
  class(ImportanceSampler), intent(inout) :: self
  type(AfqmcWalker), intent(inout)        :: walker
  !local
  integer                                 :: ng
  
  ng = size(walker%random_field)

  if (.not. allocated(walker%scale)) allocate(walker%scale(ng))
  walker%scale = 2.0_wp * self%tau * real(walker%ex(1:ng) - walker%eself(1:ng), kind=wp)
end subroutine is_scale_mix


!******************************************************************************** 
!
! Factory to select shift cutoff value:
!
!    MZ:
!        M. Motta, S. Zhang, WIRES Comp. Mol. Science 8, 1364 (2018)
!    DRV:
!        M.F. DePasquale, S.M. Rothstein, J. Vrbik, J. Chem. Phys. 89, 3629 (1988)
!
!******************************************************************************** 
subroutine is_shift_cutoff_factory(self)
  class(ImportanceSampler), intent(inout) :: self

  select case (trim(self%i_shift_cutoff))
    case("mz")
      self%shift_cutoff = 1.0_wp
    case("drv")
      self%shift_cutoff = 1.0_wp / abs(self%isqtau)
    case default
      self%shift_cutoff = 1.0_wp
  end select
end subroutine is_shift_cutoff_factory


!******************************************************************************** 
!
! Factory to select cut shift mode
!
!******************************************************************************** 
subroutine is_cut_shift_factory(self)
  class(ImportanceSampler), intent(inout) :: self

  select case (trim(self%i_cut_shift))
    case("rescale")
      self%cut_shift => is_cut_shift_rescale
    case("cap")
      self%cut_shift => is_cut_shift_cap
    case("remove")
      self%cut_shift => is_cut_shift_remove
    case("remove_field")
      self%cut_shift => is_cut_shift_remove_field
    case("unr")
      self%cut_shift => is_cut_shift_unr
    case ("none")
      self%cut_shift => is_cut_shift_none 
    case default
      self%cut_shift => is_cut_shift_remove
  end select
end subroutine is_cut_shift_factory


!******************************************************************************** 
!
! Shift cutoff mode: rescale
!
!******************************************************************************** 
subroutine is_cut_shift_rescale(self, walker)
  class(ImportanceSampler), intent(inout) :: self
  type(AfqmcWalker), intent(inout)        :: walker

  where (abs(walker%shift) > self%shift_cutoff)
    walker%shift = walker%shift / abs(walker%shift)
    self%rare_events = self%rare_events + 1
  end where
end subroutine is_cut_shift_rescale


!******************************************************************************** 
!
! Shift cutoff mode: cap
!
!******************************************************************************** 
subroutine is_cut_shift_cap(self, walker)
  class(ImportanceSampler), intent(inout) :: self
  type(AfqmcWalker), intent(inout)        :: walker

  where (abs(walker%shift) > self%shift_cutoff)
    walker%shift = onec
    self%rare_events = self%rare_events + 1
  end where
end subroutine is_cut_shift_cap


!******************************************************************************** 
!
! Shift cutoff mode: remove
!
!******************************************************************************** 
subroutine is_cut_shift_remove(self, walker)
  class(ImportanceSampler), intent(inout) :: self
  type(AfqmcWalker), intent(inout)        :: walker

  where (abs(walker%shift) > self%shift_cutoff)
    walker%shift = zeroc
    self%rare_events = self%rare_events + 1
  end where
end subroutine is_cut_shift_remove


!******************************************************************************** 
!
! Shift cutoff mode: remove_field
!
!******************************************************************************** 
subroutine is_cut_shift_remove_field(self, walker)
  class(ImportanceSampler), intent(inout) :: self
  type(AfqmcWalker), intent(inout)        :: walker

  where (abs(walker%shift) > self%shift_cutoff)
    walker%shift = zeroc
    walker%random_field = zeror
    self%rare_events = self%rare_events + 1
  end where
end subroutine is_cut_shift_remove_field


!******************************************************************************** 
!
! Shift cutoff mode: none
!
!******************************************************************************** 
subroutine is_cut_shift_none(self, walker)
  class(ImportanceSampler), intent(inout) :: self
  type(AfqmcWalker), intent(inout)        :: walker

  where (abs(walker%shift) > self%shift_cutoff)
    self%rare_events = self%rare_events + 1
  end where
end subroutine is_cut_shift_none


!******************************************************************************** 
!
! Shift cutoff mode: unr
!
!    C. Umrigar, J. Chem. Phys. 99, 2865–2890 (1993)
!    
!******************************************************************************** 
subroutine is_cut_shift_unr(self, walker)
  class(ImportanceSampler), intent(inout) :: self
  type(AfqmcWalker), intent(inout)        :: walker
  !local variables
  real(wp), parameter                     :: a = 1.0_wp

  where(abs(walker%shift) > 0.0001_wp)
    walker%shift = (sqrt(1.0_wp + 2.0_wp*a*abs(walker%shift)**2) - 1.0_wp) / (a*abs(walker%shift)**2) * walker%shift
  end where
end subroutine is_cut_shift_unr

end module afqmc_importance_sampling