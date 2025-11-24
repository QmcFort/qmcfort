! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module brng_factories

use constants
use string, only: lowercase
use mpi, only: mpi_communicator, comm_world
use rng_utils
use basic_rng, only: BasicRng
use intrinsic_basic_rng, only: IntrinsicBasicRng, init_intrinsic_brng
use lcg32_basic_rng, only: Lcg32BasicRng, init_lcg32_brng
use lcg48_basic_rng, only: Lcg48BasicRng, init_lcg48_brng
use lcg64_basic_rng, only: Lcg64BasicRng, init_lcg64_brng
#ifdef MKL
use vsl_basic_rng, only: VslBasicRng, init_vsl_brng
#endif

implicit none

public

!used to provide one brng for entire program
class(BasicRng), allocatable :: brng_glb

enum, bind(C)
  enumerator :: INTRINSIC_BRNG=1, LCG32_BRNG, LCG48_BRNG, LCG64_BRNG, VSL_BRNG
end enum

contains

!******************************************************************************** 
!
! Return brng type based on the given brng_name string
!
!******************************************************************************** 
function get_brng_type(brng_name) result(brng_type)
  character(len=*), intent(in) :: brng_name
  integer                      :: brng_type

  select case (lowercase(brng_name))
    case ("intrinsic") 
      brng_type = INTRINSIC_BRNG
    case ("lcg32") 
      brng_type = LCG32_BRNG
    case ("lcg48")
      brng_type = LCG48_BRNG
    case ("lcg64")
      brng_type = LCG64_BRNG
    case ("vsl_brng_mcg31", "vsl_brng_r250", "vsl_brng_mt19937")
      brng_type = VSL_BRNG
    case default
!debug: rng: warning that the brng is not recognized and the default one will be used
      brng_type = LCG48_BRNG
  end select
end function get_brng_type


!******************************************************************************** 
!
! Create a main brng of the qmcfort program
! Check the input flags to setup brng properly
!
!******************************************************************************** 
subroutine create_qmcfort_brng(brng)
  use qmcfort_in, only: add_input
  use qmcfort_io, only: io
  class(BasicRng), allocatable, intent(out) :: brng
  !local
  integer(i8)                               :: seed
  character(len=charlen)                    :: brng_name, seed_char
  logical                                   :: mpi_distinct, omp_distinct, found

  brng_name = "lcg48"
  call add_input("brng", brng_name)

  call add_input("brng_seed", seed_char, found=found) 
  if (found) then
    if (lowercase(seed_char) == "random") then 
      call get_random_seed(seed)
    else 
      call add_input("brng_seed", seed)
    end if
  else
    seed = def_brng_seed
  end if

  mpi_distinct = .true.
  call add_input("mpi_distinct", mpi_distinct)

  omp_distinct = .true.
  call add_input("omp_distinct", omp_distinct)

  call create_brng(brng_name, brng, seed, comm_world, mpi_distinct, omp_distinct)

  if (comm_world%mpirank == 0) call brng%report(io%qmcfort_out%funit)
end subroutine create_qmcfort_brng


!******************************************************************************** 
!
! Factory routine to construct a polymorphic BasicRng object
! based on the given brng_name string.
!
!******************************************************************************** 
subroutine create_brng(brng_name, brng, seed, comm, mpi_distinct, omp_distinct)
  character(len=*), intent(in)                 :: brng_name
  class(BasicRng), allocatable, intent(out)    :: brng
  integer(i8), optional, intent(in)            :: seed
  type(mpi_communicator), optional, intent(in) :: comm
  logical, optional, intent(in)                :: mpi_distinct
  logical, optional, intent(in)                :: omp_distinct
  !local
  integer                                      :: brng_type

  brng_type = get_brng_type(brng_name)

  select case (brng_type)
    case (INTRINSIC_BRNG)
      call create_intrinsic_brng(brng, seed, comm, mpi_distinct, omp_distinct)
    case (LCG32_BRNG)
      call create_lcg32_brng(brng, seed, comm, mpi_distinct, omp_distinct)
    case (LCG48_BRNG) 
      call create_lcg48_brng(brng, seed, comm, mpi_distinct, omp_distinct)
    case (LCG64_BRNG)
      call create_lcg64_brng(brng, seed, comm, mpi_distinct, omp_distinct)
#ifdef MKL
    case (VSL_BRNG)
      call create_vsl_brng(brng_name, brng, seed, comm, mpi_distinct, omp_distinct)
#endif
  end select
end subroutine create_brng


!******************************************************************************** 
!
! Allocate and initialize a IntrinsicRng instance, 
! and use it to create polymorphic BasicRng object.
!
!******************************************************************************** 
subroutine create_intrinsic_brng(brng, seed, comm, mpi_distinct, omp_distinct)
  class(BasicRng), allocatable, intent(out)    :: brng
  integer(i8), optional, intent(in)            :: seed
  type(mpi_communicator), optional, intent(in) :: comm
  logical, optional, intent(in)                :: mpi_distinct
  logical, optional, intent(in)                :: omp_distinct
  !local
  type(IntrinsicBasicRng), allocatable         :: intrinsic_brng

  allocate(IntrinsicBasicRng:: intrinsic_brng)
  call init_intrinsic_brng(intrinsic_brng, seed, comm, mpi_distinct, omp_distinct)
  call move_alloc(intrinsic_brng, brng)
end subroutine create_intrinsic_brng


!******************************************************************************** 
!
! Allocate and initialize a Lcg32Rng instance, 
! and use it to create polymorphic BasicRng object.
!
!******************************************************************************** 
subroutine create_lcg32_brng(brng, seed, comm, mpi_distinct, omp_distinct)
  class(BasicRng), allocatable, intent(out)    :: brng
  integer(i8), optional, intent(in)            :: seed
  type(mpi_communicator), optional, intent(in) :: comm
  logical, optional, intent(in)                :: mpi_distinct
  logical, optional, intent(in)                :: omp_distinct
  !local
  type(Lcg32BasicRng), allocatable             :: lcg32_brng

  allocate(Lcg32BasicRng:: lcg32_brng)
  call init_lcg32_brng(lcg32_brng, seed, comm, mpi_distinct, omp_distinct)
  call move_alloc(lcg32_brng, brng)
end subroutine create_lcg32_brng


!******************************************************************************** 
!
! Allocate and initialize a Lcg48Rng instance, 
! and use it to create polymorphic BasicRng object.
!
!******************************************************************************** 
subroutine create_lcg48_brng(brng, seed, comm, mpi_distinct, omp_distinct)
  class(BasicRng), allocatable, intent(out)    :: brng
  integer(i8), optional, intent(in)            :: seed
  type(mpi_communicator), optional, intent(in) :: comm
  logical, optional, intent(in)                :: mpi_distinct
  logical, optional, intent(in)                :: omp_distinct
  !local
  type(Lcg48BasicRng), allocatable             :: lcg48_brng

  allocate(Lcg48BasicRng:: lcg48_brng)
  call init_lcg48_brng(lcg48_brng, seed, comm, mpi_distinct, omp_distinct)
  call move_alloc(lcg48_brng, brng)
end subroutine create_lcg48_brng


!******************************************************************************** 
!
! Allocate and initialize a Lcg64Rng instance, 
! and use it to create polymorphic BasicRng object.
!
!******************************************************************************** 
subroutine create_lcg64_brng(brng, seed, comm, mpi_distinct, omp_distinct)
  class(BasicRng), allocatable, intent(out)    :: brng
  integer(i8), optional, intent(in)            :: seed
  type(mpi_communicator), optional, intent(in) :: comm
  logical, optional, intent(in)                :: mpi_distinct
  logical, optional, intent(in)                :: omp_distinct
  !local
  type(Lcg64BasicRng), allocatable             :: lcg64_brng

  allocate(Lcg64BasicRng:: lcg64_brng)
  call init_lcg64_brng(lcg64_brng, seed, comm, mpi_distinct, omp_distinct)
  call move_alloc(lcg64_brng, brng)
end subroutine create_lcg64_brng


!******************************************************************************** 
!
! Allocate and initialize a VslRng instance, 
! and use it to create polymorphic BasicRng object.
!
!******************************************************************************** 
#ifdef MKL
subroutine create_vsl_brng(name, brng, seed, comm, mpi_distinct, omp_distinct)
  character(len=*), intent(in)                 :: name
  class(BasicRng), allocatable, intent(out)    :: brng
  integer(i8), optional, intent(in)            :: seed
  type(mpi_communicator), optional, intent(in) :: comm
  logical, optional, intent(in)                :: mpi_distinct
  logical, optional, intent(in)                :: omp_distinct
  !local
  type(VslBasicRng), allocatable               :: vsl_brng

  allocate(VslBasicRng:: vsl_brng)
  call init_vsl_brng(name, vsl_brng, seed, comm, mpi_distinct, omp_distinct)
  call move_alloc(vsl_brng, brng)
end subroutine create_vsl_brng
#endif

end module brng_factories