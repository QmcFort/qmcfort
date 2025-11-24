! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

#ifdef MKL
include "mkl_vsl.f90"
#endif

module vsl_basic_rng

use constants
use mpi, only: mpi_communicator
use rng_utils, only: sp32_resolution, dp32_resolution, def_brng_seed
use basic_rng, only: BasicRng

#ifdef MKL
use mkl_vsl
use mkl_vsl_type

implicit none

private
public :: VslBasicRng, init_vsl_brng

type, extends(BasicRng) :: VslBasicRng
  integer                :: seed
  integer                :: rng_type
  type(vsl_stream_state) :: stream
contains
  procedure, private :: set_vsl_brng_type

  !deferred routines
  procedure :: clone => clone_vsl
  procedure :: get_seed => get_seed_vsl
  procedure :: set_seed => set_seed_vsl
  procedure :: reset => reset_vsl
  procedure :: skip => skip_vsl
  procedure :: rand_raw_real_1d_sp
  procedure :: rand_raw_real_1d_dp
end type VslBasicRng

interface VslBasicRng
  module procedure :: finit_vsl_brng
end interface VslBasicRng

contains

!******************************************************************************** 
!
! Initialization of the VslBasicRng object
!
!******************************************************************************** 
subroutine init_vsl_brng(name, self, seed, comm, mpi_distinct, omp_distinct)
  character(len=*), intent(in)                 :: name
  type(VslBasicRng), intent(out)               :: self
  integer(i8), optional, intent(in)            :: seed
  type(mpi_communicator), optional, intent(in) :: comm
  logical, optional, intent(in)                :: mpi_distinct
  logical, optional, intent(in)                :: omp_distinct
  !local
  integer(i8)                                  :: seed_

  call self%base_init(name, comm, mpi_distinct, omp_distinct)

  call self%set_vsl_brng_type(name)

  seed_ = def_brng_seed
  if (present(seed)) seed_ = seed

  call  self%set_seed(seed_)
end subroutine init_vsl_brng


!******************************************************************************** 
!
! VslRng constructor
!
!******************************************************************************** 
function finit_vsl_brng(name, seed, comm, mpi_distinct, omp_distinct) result(self)
  character(len=*), intent(in)                 :: name
  integer(i8), optional, intent(in)            :: seed
  type(mpi_communicator), optional, intent(in) :: comm
  logical, optional, intent(in)                :: mpi_distinct
  logical, optional, intent(in)                :: omp_distinct
  type(VslBasicRng)                            :: self

  call init_vsl_brng(name, self, seed, comm, mpi_distinct, omp_distinct)
end function finit_vsl_brng


!******************************************************************************** 
!
! Setup RNG type from VSL library 
!
!******************************************************************************** 
subroutine set_vsl_brng_type(self, name)
  use string, only: lowercase
  class(VslBasicRng), intent(inout) :: self
  character(len=*), intent(in)      :: name

  select case (lowercase(name))
    case ("vsl_brng_mcg31")
      self%rng_type = VSL_BRNG_MCG31
    case ("vsl_brng_r250")
      self%rng_type = VSL_BRNG_R250 !fast skipping does not work for this one
    case ("vsl_brng_mt19937")
      self%rng_type = VSL_BRNG_MT19937
    case default
!debug: rng: warning and select default rng
      self%rng_type = VSL_BRNG_MCG31
  end select
end subroutine set_vsl_brng_type


!******************************************************************************** 
!
! Clone an VslRng object
!
! Needed for safe drawing of random numbers within parallel OpenMP regions
!
!******************************************************************************** 
subroutine clone_vsl(self, clone)
  class(VslBasicRng), intent(in)            :: self
  class(BasicRng), allocatable, intent(out) :: clone
  !local
  integer                                   :: ierr
  type(VslBasicRng), allocatable            :: clone_
  type(vsl_stream_state)                    :: stream

  ierr = VslCopyStream(stream, self%stream)

  allocate(clone_, source=self)
  clone_%stream = stream

  call move_alloc(clone_, clone)
end subroutine clone_vsl


!******************************************************************************** 
!
! Retrieve a current seed
!
!******************************************************************************** 
function get_seed_vsl(self) result(seed)
  class(VslBasicRng), intent(in) :: self
  integer(i8)                    :: seed

  seed = self%seed
end function get_seed_vsl


!******************************************************************************** 
!
! Setup seed and initialize the rng state
!
! Performs a safe move of the bits between int32 and int64.
!
!******************************************************************************** 
subroutine set_seed_vsl(self, seed) 
  class(VslBasicRng), intent(inout) :: self
  integer(i8), intent(in)           :: seed
  !local
  integer                           :: err

  self%seed = int(iand(seed, maskr(32, kind=i8)), kind=i4)
  err = VslNewStream(self%stream, self%rng_type, self%seed)
end subroutine set_seed_vsl


!******************************************************************************** 
!
! State is returned to the saved seed
!
!******************************************************************************** 
subroutine reset_vsl(self)
  class(VslBasicRng), intent(inout) :: self
  !local
  integer                           :: err

  err = VslNewStream(self%stream, self%rng_type, self%seed)
end subroutine reset_vsl


!******************************************************************************** 
!
! Skip nskips states x_n --> x_{n+nskips}
!
!    calls fast skipping if implemented
!
!    otherwise draw nskips random numbers
!
!******************************************************************************** 
subroutine skip_vsl(self, nskips)
  class(VslBasicRng), intent(inout) :: self
  integer, intent(in)               :: nskips
  !local
  integer                           :: err
  integer(i8)                       :: nskips_
  real(sp), allocatable             :: sample(:)

  nskips_ = nskips
  err = VslSkipAheadStream(self%stream, nskips_)  

  if (err == VSL_RNG_ERROR_SKIPAHEAD_UNSUPPORTED) then
    allocate(sample(nskips))  
    call self%rand_raw(sample)
  end if
end subroutine skip_vsl


!******************************************************************************** 
!
! Draw uniform 1d random array in signle precision
!
!******************************************************************************** 
subroutine rand_raw_real_1d_sp(self, sample)
  class(VslBasicRng), intent(inout) :: self
  real(sp), intent(out)             :: sample(:)
  !local 
  integer                           :: err
  real(sp), parameter               :: a=0.0_sp, b=1.0_sp

  err = VsRngUniform(VSL_RNG_METHOD_UNIFORM_STD, self%stream, size(sample), sample, a, b)

  !Ensure that random numbers are in the range (0,1].
  !vsl rng yields random numbers in range [0,1).
  sample = 1.0_sp - sample
end subroutine rand_raw_real_1d_sp


!******************************************************************************** 
!
! Draw uniform 1d random array in double precision
!
!******************************************************************************** 
subroutine rand_raw_real_1d_dp(self, sample)
  class(VslBasicRng), intent(inout) :: self
  real(dp), intent(out)             :: sample(:)
  !local 
  integer                           :: err
  real(dp), parameter               :: a=0.0_dp, b=1.0_dp

  err = VdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, self%stream, size(sample), sample, a, b)

  !Ensure that random numbers are in the range (0,1].
  !vsl rng yields random numbers in range [0,1).
  sample = 1.0_dp - sample
end subroutine rand_raw_real_1d_dp
#endif
end module vsl_basic_rng