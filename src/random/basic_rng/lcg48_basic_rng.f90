! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module lcg48_basic_rng

use constants
use mpi, only: mpi_communicator
use rng_utils, only: sp48_resolution, dp48_resolution, def_brng_seed
use basic_rng, only: BasicRng

implicit none

private
public :: Lcg48BasicRng, init_lcg48_brng

integer(i8), parameter :: def_lcg48_mult = 44485709377909_i8
integer(i8), parameter :: def_lcg48_inc = 96309754297_i8

integer(i8), parameter :: mask24 = maskr(24, i8)
integer(i8), parameter :: mask48 = maskr(48, i8)

type, extends(BasicRng) :: Lcg48BasicRng
  integer(i8) :: state
  integer(i8) :: seed
  integer(i8) :: mult
  integer(i8) :: inc
contains
  procedure, private :: next_state => next_state_lcg48

  !deferred routines
  procedure :: clone => clone_lcg48
  procedure :: get_seed => get_seed_lcg48
  procedure :: set_seed => set_seed_lcg48
  procedure :: reset => reset_lcg48
  procedure :: skip => skip_lcg48
  procedure :: rand_raw_real_1d_sp
  procedure :: rand_raw_real_1d_dp
end type Lcg48BasicRng

interface Lcg48BasicRng
  module procedure :: finit_lcg48_brng
end interface Lcg48BasicRng

contains

!******************************************************************************** 
!
! Initialization of the Lcg48Rng object
!
!******************************************************************************** 
subroutine init_lcg48_brng(self, seed, comm, mpi_distinct, omp_distinct)
  type(Lcg48BasicRng), intent(out)             :: self
  integer(i8), optional, intent(in)            :: seed
  type(mpi_communicator), optional, intent(in) :: comm
  logical, optional, intent(in)                :: mpi_distinct
  logical, optional, intent(in)                :: omp_distinct
  !local
  integer(i8)                                  :: seed_

  call self%base_init("LCG48", comm, mpi_distinct, omp_distinct)

  self%mult = def_lcg48_mult
  self%inc = def_lcg48_inc

  seed_ = def_brng_seed
  if (present(seed)) seed_ = seed

  call  self%set_seed(seed_)
end subroutine init_lcg48_brng


!******************************************************************************** 
!
! Lcg48Rng constructor
!
!******************************************************************************** 
function finit_lcg48_brng(seed, comm, mpi_distinct, omp_distinct) result(self)
  integer(i8), optional, intent(in)            :: seed
  type(mpi_communicator), optional, intent(in) :: comm
  logical, optional, intent(in)                :: mpi_distinct
  logical, optional, intent(in)                :: omp_distinct
  type(Lcg48BasicRng)                          :: self

  call init_lcg48_brng(self, seed, comm, mpi_distinct, omp_distinct)
end function finit_lcg48_brng


!******************************************************************************** 
!
! Update Lcg48Rng state:
!
!    x_{n+1} = a x_{n} + c mod m
!
!******************************************************************************** 
subroutine next_state_lcg48(self)
  class(Lcg48BasicRng), intent(inout) :: self

  !emulates 48-bit integer multiplication
  call emulate_mult_48(self%mult, self%inc, self%state)
end subroutine next_state_lcg48


!******************************************************************************** 
!
! Clone an Lcg48Rng object
!
! Needed for safe drawing of random numbers within parallel OpenMP regions
!
!******************************************************************************** 
subroutine clone_lcg48(self, clone)
  class(Lcg48BasicRng), intent(in)          :: self
  class(BasicRng), allocatable, intent(out) :: clone
  !local
  type(Lcg48BasicRng), allocatable          :: clone_

  allocate(clone_, source=self)
  call move_alloc(clone_, clone)
end subroutine clone_lcg48


!******************************************************************************** 
!
! Retrieve a current seed
!
!******************************************************************************** 
function get_seed_lcg48(self) result(seed)
  class(Lcg48BasicRng), intent(in) :: self
  integer(i8)                      :: seed

  seed = self%seed
end function get_seed_lcg48


!******************************************************************************** 
!
! Setup seed and initialize the rng state
!
!******************************************************************************** 
subroutine set_seed_lcg48(self, seed) 
  class(Lcg48BasicRng), intent(inout) :: self
  integer(i8), intent(in)             :: seed

  call mvbits(seed, 0, 48, self%seed, 0)
  self%state = self%seed
end subroutine set_seed_lcg48


!******************************************************************************** 
!
! State is returned to the saved seed
!
!******************************************************************************** 
subroutine reset_lcg48(self)
  class(Lcg48BasicRng), intent(inout) :: self

  self%state = self%seed
end subroutine reset_lcg48


!******************************************************************************** 
!
! Skip n states without iterating state by state using 
!
!    x_{n} = a^{n} x_0 + c \frac{a^{n}-1}{a-1}  mod m
!
! Scales as O(log2(n)) instead of O(n)
!
!******************************************************************************** 
subroutine skip_lcg48(self, nskips)
  class(Lcg48BasicRng), intent(inout) :: self
  integer, intent(in)                 :: nskips
  !local
  integer(i8)                         :: acc_mult, acc_inc
  integer(i8)                         :: cur_mult, cur_inc
  integer(i8)                         :: temp_mult, temp_inc
  integer(i8)                         :: n

  acc_mult = 1_i8
  acc_inc = 0_i8

  cur_mult = self%mult
  cur_inc = self%inc

  n = int(nskips, kind=i8)

  do while (n > 0)
    if (iand(n, 1_i8) == 1_i8) then
      call emulate_mult_48(cur_mult, cur_inc, acc_inc)
      call emulate_mult_48(cur_mult, 0_i8, acc_mult)
    end if

    temp_mult = cur_mult
    temp_inc = cur_inc

    call emulate_mult_48(cur_mult, temp_inc, cur_inc)
    call emulate_mult_48(temp_mult, 0_i8, cur_mult)

    n = ishft(n, -1)
  end do
  
  call emulate_mult_48(acc_mult, acc_inc, self%state)
end subroutine skip_lcg48


!******************************************************************************** 
!
! Draw uniform 1d random array in signle precision
!
!******************************************************************************** 
subroutine rand_raw_real_1d_sp(self, sample)
  class(Lcg48BasicRng), intent(inout) :: self
  real(sp), intent(out)               :: sample(:)
  !local
  integer                             :: i

  do i = 1, size(sample)
    call self%next_state()
    sample(i) = sp48_resolution * self%state
  end do

  sample = 1.0_sp - sample
end subroutine rand_raw_real_1d_sp


!******************************************************************************** 
!
! Draw uniform 1d random array in double precision
!
!******************************************************************************** 
subroutine rand_raw_real_1d_dp(self, sample)
  class(Lcg48BasicRng), intent(inout) :: self
  real(dp), intent(out)               :: sample(:)
  !local  
  integer                             :: i

  do i = 1, size(sample)
    call self%next_state()
    sample(i) = dp48_resolution * self%state
  end do

  sample = 1.0_dp - sample
end subroutine rand_raw_real_1d_dp


!******************************************************************************** 
!
! Computes state update for LCG48:
!    new_state = mult * state + inc mod 2^{48}
!
! Emulates 48-bit integer multiplication to avoid overflow issues:
!    x ---> x_hi * 2^{24} + x_lo
!    y ---> y_hi * 2^{24} + y_lo
!
!    x*y + inc mod 2^{48} =  x_hi*y_hi*2^{48} + (x_hi*y_lo+x_lo*y_hi)*2^{24} +
!                            x_lo*y_lo + inc   mod 2^{48}
!                         =  (x_hi*y_lo+x_lo*y_hi)*2^{24} + x_lo*y_lo + inc
!
!******************************************************************************** 
subroutine emulate_mult_48(mult, inc, state)
  integer(i8), intent(in)    :: mult
  integer(i8), intent(in)    :: inc
  integer(i8), intent(inout) :: state
  !local
  integer(i8)                :: mult_hi, mult_lo
  integer(i8)                :: state_hi, state_lo

  mult_hi = iand(ishft(mult, -24), mask24)
  mult_lo = iand(mult, mask24)

  state_hi = iand(ishft(state, -24), mask24)
  state_lo = iand(state, mask24)

  state = iand(ishft(iand(mult_hi*state_lo+mult_lo*state_hi, mask24), 24) + mult_lo*state_lo + inc, mask48)
end subroutine emulate_mult_48

end module lcg48_basic_rng