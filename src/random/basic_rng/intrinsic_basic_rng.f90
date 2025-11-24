! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module intrinsic_basic_rng

use constants
use mpi, only: mpi_communicator
use rng_utils, only: sp32_resolution, dp32_resolution, def_brng_seed
use basic_rng, only: BasicRng

implicit none

private
public :: IntrinsicBasicRng, init_intrinsic_brng

type, extends(BasicRng) :: IntrinsicBasicRng
contains
  !deferred routines
  procedure :: clone => clone_intrinsic
  procedure :: get_seed => get_seed_intrinsic
  procedure :: set_seed => set_seed_intrinsic
  procedure :: reset => reset_intrinsic
  procedure :: skip => skip_intrinsic
  procedure :: rand_raw_real_1d_sp
  procedure :: rand_raw_real_1d_dp
end type IntrinsicBasicRng

interface IntrinsicBasicRng
  module procedure :: finit_intrinsic_brng
end interface IntrinsicBasicRng

contains

!******************************************************************************** 
!
! Initialization of the IntrinsicRng object
!
!******************************************************************************** 
subroutine init_intrinsic_brng(self, seed, comm, mpi_distinct, omp_distinct)
  type(IntrinsicBasicRng), intent(out)         :: self
  integer(i8), optional, intent(in)            :: seed
  type(mpi_communicator), optional, intent(in) :: comm
  logical, optional, intent(in)                :: mpi_distinct
  logical, optional, intent(in)                :: omp_distinct
  !local
  integer(i8)                                  :: seed_

  call self%base_init("intrinsic", comm, mpi_distinct, omp_distinct)

  seed_ = def_brng_seed
  if (present(seed)) seed_ = seed

  call  self%set_seed(seed_)
end subroutine init_intrinsic_brng


!******************************************************************************** 
!
! IntrinsicRng constructor
!
!******************************************************************************** 
function finit_intrinsic_brng(seed, comm, mpi_distinct, omp_distinct) result(self)
  integer(i8), optional, intent(in)            :: seed
  type(mpi_communicator), optional, intent(in) :: comm
  logical, optional, intent(in)                :: mpi_distinct
  logical, optional, intent(in)                :: omp_distinct
  type(IntrinsicBasicRng)                      :: self

  call init_intrinsic_brng(self, seed, comm, mpi_distinct, omp_distinct)
end function finit_intrinsic_brng


!******************************************************************************** 
!
! Clone an IntrinsicRng object
!
! Needed for safe drawing of random numbers within parallel OpenMP regions
!
!******************************************************************************** 
subroutine clone_intrinsic(self, clone)
  class(IntrinsicBasicRng), intent(in)      :: self
  class(BasicRng), allocatable, intent(out) :: clone
  !local
  type(IntrinsicBasicRng), allocatable      :: clone_

  allocate(clone_, source=self)
  call move_alloc(clone_, clone)
end subroutine clone_intrinsic


!******************************************************************************** 
!
! Retrieve a current seed
!
! Note: seed in BasicRng is implemented as int64, while 
!       random_seed() uses default integer size
!
!******************************************************************************** 
function get_seed_intrinsic(self) result(seed)
  class(IntrinsicBasicRng), intent(in) :: self
  integer(i8)                          :: seed
  !local 
  integer                              :: seed_size
  integer, allocatable                 :: seed_(:)

  call random_seed(size=seed_size)

  allocate(seed_(seed_size))
  call random_seed(get=seed_)

  !copy at most 64 bits from the intrinsic seed to the qmcfort seed
  if (bit_size(seed_size) == 32) then
    if (seed_size == 1) then
      seed = seed_(1)
    else
      seed = ior(shiftl(int(seed_(2), kind=i8), 32), int(seed_(1), kind=i8))
    end if
  else
    seed = seed_(1)
  end if
end function get_seed_intrinsic


!******************************************************************************** 
!
! Setup seed and initialize the rng state
!
! Performs a safe move of the bits between int32 and int64.
!
!******************************************************************************** 
subroutine set_seed_intrinsic(self, seed) 
  class(IntrinsicBasicRng), intent(inout) :: self
  integer(i8), intent(in)                 :: seed
  !local
  integer                                 :: seed_size, i
  integer(i8)                             :: mask
  integer, allocatable                    :: seed_(:)

  call random_seed(size=seed_size)
  allocate(seed_(seed_size))
  seed_ = 0

  !copy at most 64 bits to the intrinsic seed
  if (bit_size(seed_size) == 32) then
    seed_(1) = int(iand(seed, maskr(32, kind=i8)), kind=i4)
    if (seed_size >= 1) then
      seed_(2) = int(shiftr(seed, 32), kind=i4)
    end if
  else
    seed_(1) = seed
  end if

  if (self%comm%mpirank == 0) call random_seed(put=seed_)
end subroutine set_seed_intrinsic


!******************************************************************************** 
!
! State is returned to the saved seed
!
!******************************************************************************** 
subroutine reset_intrinsic(self)
  class(IntrinsicBasicRng), intent(inout) :: self
  !local
  integer                                 :: seed_size
  integer, allocatable                    :: seed(:)

  if (self%comm%mpirank == 0) then
    call random_seed(size=seed_size)
    allocate(seed(seed_size))

    call random_seed(get=seed)
    call random_seed(put=seed)
  end if
end subroutine reset_intrinsic


!******************************************************************************** 
!
! Skip n states without iterating state by state using 
!
!    x_{n} = a^{n} x_0 + c \frac{a^{n}-1}{a-1}  mod m
!
! Scales as O(log2(n)) instead of O(n)
!
!******************************************************************************** 
subroutine skip_intrinsic(self, nskips)
  class(IntrinsicBasicRng), intent(inout) :: self
  integer, intent(in)                     :: nskips
  !local
  real(sp)                                :: sample(nskips)

  call random_number(sample)
end subroutine skip_intrinsic


!******************************************************************************** 
!
! Draw uniform 1d random array in signle precision
!
!******************************************************************************** 
subroutine rand_raw_real_1d_sp(self, sample)
  class(IntrinsicBasicRng), intent(inout) :: self
  real(sp), intent(out)                   :: sample(:)

  call random_number(sample)

  !Ensure that random numbers are in the range (0,1].
  !Intrinsic rng yields random numbers in range [0,1).
  sample = 1.0_sp - sample
end subroutine rand_raw_real_1d_sp


!******************************************************************************** 
!
! Draw uniform 1d random array in double precision
!
!******************************************************************************** 
subroutine rand_raw_real_1d_dp(self, sample)
  class(IntrinsicBasicRng), intent(inout) :: self
  real(dp), intent(out)                   :: sample(:)

  call random_number(sample)

  !Ensure that random numbers are in the range (0,1].
  !Intrinsic rng yields random numbers in range [0,1).
  sample = 1.0_dp - sample
end subroutine rand_raw_real_1d_dp

end module intrinsic_basic_rng