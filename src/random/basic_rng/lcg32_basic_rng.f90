! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module lcg32_basic_rng

use constants
use mpi, only: mpi_communicator
use rng_utils, only: sp32_resolution, dp32_resolution, def_brng_seed
use basic_rng, only: BasicRng

implicit none

private
public :: Lcg32BasicRng, init_lcg32_brng

integer(i4), parameter :: def_lcg32_mult = int(mod(2891336453_i8, 2**32_i8), kind=i4) !2891336453 can't be represented using signed integers
integer(i4), parameter :: def_lcg32_inc = 1043493431_i4

type, extends(BasicRng) :: Lcg32BasicRng
  integer(i4) :: state
  integer(i4) :: seed
  integer(i4) :: mult
  integer(i4) :: inc
contains
  procedure, private :: next_state => next_state_lcg32

  !deferred routines
  procedure :: clone => clone_lcg32
  procedure :: get_seed => get_seed_lcg32
  procedure :: set_seed => set_seed_lcg32
  procedure :: reset => reset_lcg32
  procedure :: skip => skip_lcg32
  procedure :: rand_raw_real_1d_sp
  procedure :: rand_raw_real_1d_dp
end type Lcg32BasicRng

interface Lcg32BasicRng
  module procedure :: finit_lcg32_brng
end interface Lcg32BasicRng

contains

!******************************************************************************** 
!
! Initialization of the Lcg32Rng object
!
!******************************************************************************** 
subroutine init_lcg32_brng(self, seed, comm, mpi_distinct, omp_distinct)
  type(Lcg32BasicRng), intent(out)             :: self
  integer(i8), optional, intent(in)            :: seed
  type(mpi_communicator), optional, intent(in) :: comm
  logical, optional, intent(in)                :: mpi_distinct
  logical, optional, intent(in)                :: omp_distinct
  !local
  integer(i8)                                  :: seed_

  call self%base_init("LCG32", comm, mpi_distinct, omp_distinct)

  self%mult = def_lcg32_mult
  self%inc = def_lcg32_inc

  seed_ = def_brng_seed
  if (present(seed)) seed_ = seed

  call  self%set_seed(seed_)
end subroutine init_lcg32_brng


!******************************************************************************** 
!
! Lcg32Rng constructor
!
!******************************************************************************** 
function finit_lcg32_brng(seed, comm, mpi_distinct, omp_distinct) result(self)
  integer(i8), optional, intent(in)            :: seed
  type(mpi_communicator), optional, intent(in) :: comm
  logical, optional, intent(in)                :: mpi_distinct
  logical, optional, intent(in)                :: omp_distinct
  type(Lcg32BasicRng)                          :: self

  call init_lcg32_brng(self, seed, comm, mpi_distinct, omp_distinct)
end function finit_lcg32_brng


!******************************************************************************** 
!
! Update Lcg32Rng state:
!
!    x_{n+1} = a x_{n} + c mod m
!
! Relies on integer overflow, so x_n \in [-2^{31}, ..., 0, 1, ...m 2^{32}-1]
!
!******************************************************************************** 
subroutine next_state_lcg32(self)
  class(Lcg32BasicRng), intent(inout) :: self

  !mod arithmetic done automatically with integer overflow mechanism
  !gnu compiler may fail with integer overflow: try to use -fwrapv compiler option
  self%state = self%mult * self%state + self%inc
end subroutine next_state_lcg32


!******************************************************************************** 
!
! Clone an Lcg32Rng object
!
! Needed for safe drawing of random numbers within parallel OpenMP regions
!
!******************************************************************************** 
subroutine clone_lcg32(self, clone)
  class(Lcg32BasicRng), intent(in)          :: self
  class(BasicRng), allocatable, intent(out) :: clone
  !local
  type(Lcg32BasicRng), allocatable          :: clone_

  allocate(clone_, source=self)
  call move_alloc(clone_, clone)
end subroutine clone_lcg32


!******************************************************************************** 
!
! Retrieve a current seed
!
! Performs a safe move of the bits between int32 and int64.
!
!******************************************************************************** 
function get_seed_lcg32(self) result(seed)
  class(Lcg32BasicRng), intent(in) :: self
  integer(i8)                      :: seed

  seed = int(self%seed, kind=i8)
end function get_seed_lcg32


!******************************************************************************** 
!
! Setup seed and initialize the rng state
!
! Performs a safe move of the bits between int32 and int64.
!
!******************************************************************************** 
subroutine set_seed_lcg32(self, seed) 
  class(Lcg32BasicRng), intent(inout) :: self
  integer(i8), intent(in)             :: seed

  self%seed = int(iand(seed, maskr(32, kind=i8)), kind=i4)
  self%state = self%seed
end subroutine set_seed_lcg32


!******************************************************************************** 
!
! State is returned to the saved seed
!
!******************************************************************************** 
subroutine reset_lcg32(self)
  class(Lcg32BasicRng), intent(inout) :: self

  self%state = self%seed
end subroutine reset_lcg32


!******************************************************************************** 
!
! Skip n states without iterating state by state using 
!
!    x_{n} = a^{n} x_0 + c \frac{a^{n}-1}{a-1}  mod m
!
! Scales as O(log2(n)) instead of O(n)
!
!******************************************************************************** 
subroutine skip_lcg32(self, nskips)
  class(Lcg32BasicRng), intent(inout) :: self
  integer, intent(in)                 :: nskips
  !local
  integer(i4)                         :: acc_mult, acc_inc
  integer(i4)                         :: cur_mult, cur_inc
  integer(i4)                         :: n

  acc_mult = 1_i4
  acc_inc = 0_i4

  cur_mult = self%mult
  cur_inc = self%inc

  n = int(nskips, kind=i4)

  do while (n > 0)
    if (iand(n, 1_i4) == 1_i4) then
      acc_inc = acc_inc * cur_mult + cur_inc
      acc_mult = acc_mult * cur_mult
    end if

    cur_inc = cur_inc * (cur_mult + 1_i4)
    cur_mult = cur_mult * cur_mult

    n = ishft(n, -1)
  end do
  
  self%state = acc_mult * self%state + acc_inc
end subroutine skip_lcg32


!******************************************************************************** 
!
! Draw uniform 1d random array in signle precision
!
!******************************************************************************** 
subroutine rand_raw_real_1d_sp(self, sample)
  class(Lcg32BasicRng), intent(inout) :: self
  real(sp), intent(out)               :: sample(:)
  !local
  integer                             :: i

  do i = 1, size(sample)
    call self%next_state()
    sample(i) = sp32_resolution * real(self%state, sp)
  end do

  !Ensure that random numbers are in the range (0,1].
  !Due to the integer overflow arithmetic in Fortran,
  !state can have values from -2^{31} ... 2^{31} - 1 
  sample = 0.5_sp - sample
end subroutine rand_raw_real_1d_sp


!******************************************************************************** 
!
! Draw uniform 1d random array in double precision
!
!******************************************************************************** 
subroutine rand_raw_real_1d_dp(self, sample)
  class(Lcg32BasicRng), intent(inout) :: self
  real(dp), intent(out)               :: sample(:)
  !local  
  integer                             :: i

  do i = 1, size(sample)
    call self%next_state()
    sample(i) = dp32_resolution * real(self%state, dp)
  end do

  !Ensure that random numbers are in the range (0,1].
  !Due to the integer overflow arithmetic in Fortran,
  !state can have values from -2^{31} ... 2^{31} - 1 
  sample = 0.5_dp - sample
end subroutine rand_raw_real_1d_dp

end module lcg32_basic_rng