! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module basic_rng

!$ use omp_lib

use constants
use mpi, only: mpi_communicator, comm_world

implicit none

private
public :: BasicRng

type, abstract :: BasicRng
  character(len=:), allocatable :: name
  type(mpi_communicator)        :: comm
  logical                       :: mpi_distinct
  logical                       :: omp_distinct
contains
  !deferred routines
  procedure(Iclone), deferred       :: clone
  procedure(Iget_seed), deferred    :: get_seed
  procedure(Iset_seed), deferred    :: set_seed
  procedure(Ireset), deferred       :: reset
  procedure(Iskip), deferred        :: skip
  procedure(Irand_raw_sp), deferred :: rand_raw_real_1d_sp
  procedure(Irand_raw_dp), deferred :: rand_raw_real_1d_dp
  
  procedure, public                 :: base_init 
  procedure, public                 :: prepare_stream
  procedure, public                 :: align_stream
  procedure, public                 :: report => report_brng

  procedure, private                :: get_nskips
  procedure, private                :: get_nskips_total

  generic, public                   :: rand_raw => rand_raw_real_sp, rand_raw_real_dp, &
                                                   rand_raw_real_1d_sp, rand_raw_real_1d_dp, &
                                                   rand_raw_real_2d_sp, rand_raw_real_2d_dp
  procedure, private                :: rand_raw_real_sp, rand_raw_real_dp
  procedure, private                :: rand_raw_real_2d_sp, rand_raw_real_2d_dp


  generic, public                   :: rand => rand_real_sp, rand_real_dp, &
                                               rand_real_1d_sp, rand_real_1d_dp, &
                                               rand_real_2d_sp, rand_real_2d_dp, &
                                               rand_real_3d_sp, rand_real_3d_dp
  procedure, private                :: rand_real_sp, rand_real_dp
  procedure, private                :: rand_real_1d_sp, rand_real_1d_dp
  procedure, private                :: rand_real_2d_sp, rand_real_2d_dp
  procedure, private                :: rand_real_3d_sp, rand_real_3d_dp
end type BasicRng

abstract interface
  subroutine Iclone(self, clone)
    import BasicRng
    class(BasicRng), intent(in)               :: self
    class(BasicRng), allocatable, intent(out) :: clone
  end subroutine Iclone

  function Iget_seed(self) result(seed)
    import BasicRng, i8
    class(BasicRng), intent(in) :: self
    integer(i8)                 :: seed
  end function Iget_seed
  
  subroutine Iset_seed(self, seed)
    import BasicRng, i8
    class(BasicRng), intent(inout) :: self
    integer(i8), intent(in)        :: seed
  end subroutine Iset_seed

  subroutine Ireset(self)
    import BasicRng
    class(BasicRng), intent(inout) :: self
  end subroutine Ireset

  subroutine Iskip(self, nskips) 
    import BasicRng
    class(BasicRng), intent(inout) :: self
    integer, intent(in)            :: nskips
  end subroutine Iskip

  subroutine Irand_raw_sp(self, sample)
    import BasicRng, sp
    class(BasicRng), intent(inout) :: self
    real(sp), intent(out)          :: sample(:)
  end subroutine Irand_raw_sp

  subroutine Irand_raw_dp(self, sample)
    import BasicRng, dp
    class(BasicRng), intent(inout) :: self
    real(dp), intent(out)          :: sample(:)
  end subroutine Irand_raw_dp
end interface

contains


!******************************************************************************** 
!
! Basic initialization of the components common to all concrete BasicRng types
!
!******************************************************************************** 
subroutine base_init(self, name, comm, mpi_distinct, omp_distinct)
  class(BasicRng), intent(inout)               :: self
  character(len=*), intent(in)                 :: name
  type(mpi_communicator), optional, intent(in) :: comm
  logical, optional, intent(in)                :: mpi_distinct
  logical, optional, intent(in)                :: omp_distinct

  self%name = name

  self%comm = comm_world
  if (present(comm)) self%comm = comm

  self%mpi_distinct = .true.
  if (present(mpi_distinct)) self%mpi_distinct = mpi_distinct

  self%omp_distinct = .false.
  if (present(omp_distinct)) self%omp_distinct = omp_distinct
end subroutine base_init


!******************************************************************************** 
!
! Write BasicRng  object to the file
!
!******************************************************************************** 
subroutine  report_brng(self, funit)
  class(BasicRng), intent(in) :: self
  integer, intent(in)         :: funit
  !local
  integer(i8)                 :: seed

  seed = self%get_seed()

  write(funit,*)   "Basic Random Number Generator:"
  write(funit,*)   "------------------------------"
  write(funit,100) "  brng       = ", self%name
  write(funit,101) "  brng_seed  = ", seed
  write(funit,*)

  100 format (1x,a,a)
  101 format (1x,a,i0) 
end subroutine report_brng


!******************************************************************************** 
!
! Preapare an rng stream so that each OpenMP thread / MPI proc. will draw
! unique random numbers depennding on mpi_distinct and omp_distinct variables
!
!******************************************************************************** 
subroutine prepare_stream(self, sample_size, stream)
  class(BasicRng), intent(in)               :: self
  integer, intent(in)                       :: sample_size
  class(BasicRng), allocatable, intent(out) :: stream
  !local
  integer                                   :: nskips

  call self%clone(stream)
  nskips = stream%get_nskips()
  call stream%skip(nskips*sample_size)
end subroutine prepare_stream


!******************************************************************************** 
!
! Align main stream after parallel drawing of random numbers
!
!******************************************************************************** 
subroutine align_stream(self, sample_size)
  class(BasicRng), intent(inout) :: self
  integer                        :: sample_size
  !local
  integer                        :: nskips_total

  nskips_total = self%get_nskips_total()

  !$omp barrier
  !$omp single
  call self%skip(nskips_total*sample_size)
  !$omp end single
end subroutine align_stream


!******************************************************************************** 
!
! Calculate amount of states given OMP thread / MPI proc. has to skip
! The number of states depends on mpi_distinct and omp_distinct variables
!
!******************************************************************************** 
function get_nskips(self) result(nskips)
  class(BasicRng), intent(in) :: self
  integer                     :: nskips
  !local
  integer                     :: thread_id

  thread_id = 0
  !$ thread_id = omp_get_thread_num()

  if (self%mpi_distinct) then
    if (self%omp_distinct) then
      nskips = self%comm%mpirank * self%comm%ompsize + thread_id
    else
      nskips = self%comm%mpirank
    end if
  else 
    if (self%omp_distinct) then
      nskips = thread_id
    else
      nskips = 0
    end if
  end if
end function get_nskips


!******************************************************************************** 
!
! Calculate amount of states main rng has to skip after drawing random numbers
!
! Note: random numbers are drawn by thread-local rng clones, and main 
! rng skips all states drawn by clones
!
!******************************************************************************** 
function get_nskips_total(self) result(nskips)
  class(BasicRng), intent(in) :: self
  integer                     :: nskips

  if (self%mpi_distinct) then
    if (self%omp_distinct) then
      nskips = self%comm%mpisize * self%comm%ompsize
    else
      nskips = self%comm%mpisize
    end if
  else 
    if (self%omp_distinct) then
      nskips = self%comm%ompsize
    else
      nskips = 1
    end if
  end if
end function get_nskips_total


!******************************************************************************** 
!
! Generate uniformly distributed random numbers in range (0,1]
!
! Does not include MPI / OMP logic (see also BasicRng%rand() routine)
!
! rand_raw_real_1d_sp, rand_raw_real_1d_dp are deferred routines 
!
!******************************************************************************** 
subroutine rand_raw_real_sp(self, sample)
  class(BasicRng), intent(inout) :: self
  real(sp), intent(out)          :: sample
  !local
  real(sp)                       :: sample_(1)

  call self%rand_raw(sample_)
  sample = sample_(1)
end subroutine rand_raw_real_sp

subroutine rand_raw_real_dp(self, sample)
  class(BasicRng), intent(inout) :: self
  real(dp), intent(out)          :: sample
  !local
  real(dp)                       :: sample_(1)

  call self%rand_raw(sample_)
  sample = sample_(1)
end subroutine rand_raw_real_dp

subroutine rand_raw_real_2d_sp(self, sample)
  class(BasicRng), intent(inout)            :: self
  real(sp), contiguous, target, intent(out) :: sample(:,:)
  !local
  real(sp), pointer                         :: sample_(:)

  sample_(1:size(sample)) => sample
  call self%rand_raw(sample_)
end subroutine rand_raw_real_2d_sp

subroutine rand_raw_real_2d_dp(self, sample)
  class(BasicRng), intent(inout)            :: self
  real(dp), contiguous, target, intent(out) :: sample(:,:)
  !local
  real(dp), pointer                         :: sample_(:)

  sample_(1:size(sample)) => sample
  call self%rand_raw(sample_)
end subroutine rand_raw_real_2d_dp


!******************************************************************************** 
!
! Main routine of the BasicRng class - crates uniformly distributed random
! numbers in range (0,1].
!
! Includes MPI / OMP logic
!
!******************************************************************************** 
subroutine rand_real_sp(self, sample)
  class(BasicRng), intent(inout) :: self
  real(sp), intent(out)          :: sample
  !local
  real(sp)                       :: sample_(1)

  call self%rand(sample_)
  sample = sample_(1)
end subroutine rand_real_sp

subroutine rand_real_dp(self, sample)
  class(BasicRng), intent(inout) :: self
  real(dp), intent(out)          :: sample
  !local
  real(dp)                       :: sample_(1)

  call self%rand(sample_)
  sample = sample_(1)
end subroutine rand_real_dp

subroutine rand_real_1d_sp(self, sample)
  class(BasicRng), intent(inout) :: self
  real(sp), intent(out)          :: sample(:)
  !local
  class(BasicRng), allocatable   :: brng_thread

  !create RNG stream for each OpenMP thread
  call self%prepare_stream(size(sample), brng_thread)
  call brng_thread%rand_raw(sample)
  call self%align_stream(size(sample))
end subroutine rand_real_1d_sp

subroutine rand_real_1d_dp(self, sample)
  class(BasicRng), intent(inout) :: self
  real(dp), intent(out)          :: sample(:)
  !local
  class(BasicRng), allocatable   :: brng_thread

  !create RNG stream for each OpenMP thread
  call self%prepare_stream(size(sample), brng_thread)
  call brng_thread%rand_raw(sample)
  call self%align_stream(size(sample))
end subroutine rand_real_1d_dp

subroutine rand_real_2d_sp(self, sample)
  class(BasicRng), intent(inout)            :: self
  real(sp), contiguous, target, intent(out) :: sample(:,:)
  !local
  real(sp), pointer                         :: sample_(:)

  sample_(1:size(sample)) => sample
  call self%rand(sample_)
end subroutine rand_real_2d_sp

subroutine rand_real_2d_dp(self, sample)
  class(BasicRng), intent(inout)            :: self
  real(dp), contiguous, target, intent(out) :: sample(:,:)
  !local
  real(dp), pointer                         :: sample_(:)

  sample_(1:size(sample)) => sample
  call self%rand(sample_)
end subroutine rand_real_2d_dp

subroutine rand_real_3d_sp(self, sample)
  class(BasicRng), intent(inout)            :: self
  real(sp), contiguous, target, intent(out) :: sample(:,:,:)
  !local
  real(sp), pointer                         :: sample_(:)

  sample_(1:size(sample)) => sample
  call self%rand(sample_)
end subroutine rand_real_3d_sp

subroutine rand_real_3d_dp(self, sample)
  class(BasicRng), intent(inout)            :: self
  real(dp), contiguous, target, intent(out) :: sample(:,:,:)
  !local
  real(dp), pointer                         :: sample_(:)

  sample_(1:size(sample)) => sample
  call self%rand(sample_)
end subroutine rand_real_3d_dp

end module basic_rng