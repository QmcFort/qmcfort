! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module lcg64_basic_rng

use constants
use mpi, only: mpi_communicator
use rng_utils, only: sp64_resolution, dp64_resolution, def_brng_seed
use basic_rng, only: BasicRng

implicit none

private
public :: Lcg64BasicRng, init_lcg64_brng

integer(i8), parameter :: def_lcg64_mult = 2862933555777941757_i8
integer(i8), parameter :: def_lcg64_inc = 1848390349349385475_i8

type, extends(BasicRng) :: Lcg64BasicRng
  integer(i8) :: state
  integer(i8) :: seed
  integer(i8) :: mult
  integer(i8) :: inc
contains
  procedure, private :: next_state => next_state_lcg64

  !deferred routines
  procedure :: clone => clone_lcg64
  procedure :: get_seed => get_seed_lcg64
  procedure :: set_seed => set_seed_lcg64
  procedure :: reset => reset_lcg64
  procedure :: skip => skip_lcg64
  procedure :: rand_raw_real_1d_sp
  procedure :: rand_raw_real_1d_dp
end type Lcg64BasicRng

interface Lcg64BasicRng
  module procedure :: finit_lcg64_brng
end interface Lcg64BasicRng

contains

!******************************************************************************** 
!
! Initialization of the Lcg64Rng object
!
!******************************************************************************** 
subroutine init_lcg64_brng(self, seed, comm, mpi_distinct, omp_distinct)
  type(Lcg64BasicRng), intent(out)             :: self
  integer(i8), optional, intent(in)            :: seed
  type(mpi_communicator), optional, intent(in) :: comm
  logical, optional, intent(in)                :: mpi_distinct
  logical, optional, intent(in)                :: omp_distinct
  !local
  integer(i8)                                  :: seed_

  call self%base_init("LCG64", comm, mpi_distinct, omp_distinct)

  self%mult = def_lcg64_mult
  self%inc = def_lcg64_inc

  seed_ = def_brng_seed
  if (present(seed)) seed_ = seed

  call  self%set_seed(seed_)
end subroutine init_lcg64_brng


!******************************************************************************** 
!
! Lcg64Rng constructor
!
!******************************************************************************** 
function finit_lcg64_brng(seed, comm, mpi_distinct, omp_distinct) result(self)
  integer(i8), optional, intent(in)            :: seed
  type(mpi_communicator), optional, intent(in) :: comm
  logical, optional, intent(in)                :: mpi_distinct
  logical, optional, intent(in)                :: omp_distinct
  type(Lcg64BasicRng)                          :: self

  call init_lcg64_brng(self, seed, comm, mpi_distinct, omp_distinct)
end function finit_lcg64_brng


!******************************************************************************** 
!
! Update Lcg64Rng state:
!
!    x_{n+1} = a x_{n} + c mod m
!
! Relies on integer overflow, so x_n \in [-2^{63}, ..., 0, 1, ...m 2^{63}-1]
!
!******************************************************************************** 
subroutine next_state_lcg64(self)
  class(Lcg64BasicRng), intent(inout) :: self

  !mod arithmetic done automatically with integer overflow mechanism
  !gnu compiler may fail with integer overflow: try to use -fwrapv compiler option
  self%state = self%mult * self%state + self%inc
end subroutine next_state_lcg64


!******************************************************************************** 
!
! Clone an Lcg64Rng object
!
! Needed for safe drawing of random numbers within parallel OpenMP regions
!
!******************************************************************************** 
subroutine clone_lcg64(self, clone)
  class(Lcg64BasicRng), intent(in)          :: self
  class(BasicRng), allocatable, intent(out) :: clone
  !local
  type(Lcg64BasicRng), allocatable          :: clone_

  allocate(clone_, source=self)
  call move_alloc(clone_, clone)
end subroutine clone_lcg64


!******************************************************************************** 
!
! Retrieve a current seed
!
!******************************************************************************** 
function get_seed_lcg64(self) result(seed)
  class(Lcg64BasicRng), intent(in) :: self
  integer(i8)                      :: seed

  seed = self%seed
end function get_seed_lcg64


!******************************************************************************** 
!
! Setup seed and initialize the rng state
!
!******************************************************************************** 
subroutine set_seed_lcg64(self, seed) 
  class(Lcg64BasicRng), intent(inout) :: self
  integer(i8), intent(in)             :: seed

  self%seed = seed
  self%state = self%seed
end subroutine set_seed_lcg64


!******************************************************************************** 
!
! State is returned to the saved seed
!
!******************************************************************************** 
subroutine reset_lcg64(self)
  class(Lcg64BasicRng), intent(inout) :: self

  self%state = self%seed
end subroutine reset_lcg64


!******************************************************************************** 
!
! Skip n states without iterating state by state using 
!
!    x_{n} = a^{n} x_0 + c \frac{a^{n}-1}{a-1}  mod m
!
! Scales as O(log2(n)) instead of O(n)
!
!******************************************************************************** 
subroutine skip_lcg64(self, nskips)
  class(Lcg64BasicRng), intent(inout) :: self
  integer, intent(in)                 :: nskips
  !local
  integer(i8)                         :: acc_mult, acc_inc
  integer(i8)                         :: cur_mult, cur_inc
  integer(i8)                         :: n

  acc_mult = 1_i8
  acc_inc = 0_i8

  cur_mult = self%mult
  cur_inc = self%inc

  n = int(nskips, kind=i8)

  do while (n > 0)
    if (iand(n, 1_i8) == 1_i8) then
      acc_inc = acc_inc * cur_mult + cur_inc
      acc_mult = acc_mult * cur_mult
    end if

    cur_inc = cur_inc * (cur_mult + 1_i8)
    cur_mult = cur_mult * cur_mult

    n = ishft(n, -1)
  end do
  
  self%state = acc_mult * self%state + acc_inc
end subroutine skip_lcg64


!******************************************************************************** 
!
! Draw uniform 1d random array in signle precision
!
!******************************************************************************** 
subroutine rand_raw_real_1d_sp(self, sample)
  class(Lcg64BasicRng), intent(inout) :: self
  real(sp), intent(out)               :: sample(:)
  !local
  integer                             :: i

  do i = 1, size(sample)
    call self%next_state()
    sample(i) = sp64_resolution * real(self%state, sp)
  end do

  !Ensure that random numbers are in the range (0,1].
  !Due to the integer overflow arithmetic in Fortran,
  !state can have values from -2^{63} ... 2^{63} - 1 
  sample = 0.5_sp - sample
end subroutine rand_raw_real_1d_sp


!******************************************************************************** 
!
! Draw uniform 1d random array in double precision
!
!******************************************************************************** 
subroutine rand_raw_real_1d_dp(self, sample)
  class(Lcg64BasicRng), intent(inout) :: self
  real(dp), intent(out)               :: sample(:)
  !local  
  integer                             :: i

  do i = 1, size(sample)
    call self%next_state()
    sample(i) = dp64_resolution * real(self%state, dp)
  end do

  !Ensure that random numbers are in the range (0,1].
  !Due to the integer overflow arithmetic in Fortran,
  !state can have values from -2^{63} ... 2^{63} - 1 
  sample = 0.5_dp - sample
end subroutine rand_raw_real_1d_dp

end module lcg64_basic_rng