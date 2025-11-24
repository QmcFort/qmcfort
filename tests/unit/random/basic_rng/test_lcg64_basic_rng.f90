! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module test_lcg64_basic_rng

use constants

use fruit, only: assert_equals, run_test_case, test_set_summary
use lcg64_basic_rng, only: Lcg64BasicRng, init_lcg64_brng
use mpi, only: comm_world

implicit none

private
public :: test_lcg64_basic_rng_set

contains

subroutine test_lcg64_basic_rng_set
  call run_test_case(test_lcg64_seed_sequence, "check whether lcg64 sequence is correct for a random seed")
  call run_test_case(test_lcg64_seed_skip_sequence, "check whether lcg64 skip works correctly")

  call test_set_summary("src/random/basic_rng/lcg64_basic_rng.f90")
end subroutine test_lcg64_basic_rng_set

subroutine test_lcg64_seed_sequence
  integer, parameter     :: len = 50
  real(sp), parameter    :: tol = 0.0001_sp
  integer(i8), parameter :: seed = 67729737732161_i8
  real(sp)               :: expected(len), result(len)
  type(Lcg64BasicRng)    :: brng

  !obtained from the ref. code that completely follows implementation here
  expected = [0.861333_sp, 0.930726_sp, 0.879357_sp, 0.325847_sp, 0.116008_sp, &
              0.822375_sp, 0.867145_sp, 0.519971_sp, 0.211038_sp, 0.516246_sp, &
              0.350910_sp, 0.312113_sp, 0.734684_sp, 0.533851_sp, 0.075433_sp, &
              0.631136_sp, 0.276698_sp, 0.422606_sp, 0.225962_sp, 0.053225_sp, &
              0.308299_sp, 0.066579_sp, 0.105867_sp, 0.760064_sp, 0.539040_sp, &
              0.178256_sp, 0.768535_sp, 0.323788_sp, 0.295879_sp, 0.953541_sp, &
              0.268560_sp, 0.543166_sp, 0.179785_sp, 0.992371_sp, 0.389952_sp, &
              0.798773_sp, 0.187829_sp, 0.464048_sp, 0.692571_sp, 0.570227_sp, &
              0.802387_sp, 0.787631_sp, 0.440805_sp, 0.032591_sp, 0.115752_sp, &
              0.593088_sp, 0.038361_sp, 0.532033_sp, 0.581456_sp, 0.936259_sp]

  if (.not. comm_world%initialized) call comm_world%init()
  brng = Lcg64BasicRng(seed=seed, comm=comm_world, mpi_distinct=.true., omp_distinct=.false.)
  call brng%rand(result)

  call assert_equals(expected, result, len, tol)
end subroutine test_lcg64_seed_sequence

subroutine test_lcg64_seed_skip_sequence
  integer, parameter     :: len=50, skip=200
  real(dp), parameter    :: tol = 0.0001_dp
  real(dp)               :: expected(len), result(len), temp(skip)
  type(Lcg64BasicRng)    :: brng

  if (.not. comm_world%initialized) call comm_world%init()
  brng = Lcg64BasicRng(comm=comm_world, mpi_distinct=.true., omp_distinct=.false.)

  call brng%rand(temp)
  call brng%rand(expected)

  call brng%reset()

  call brng%skip(skip)
  call brng%rand(result)

  call assert_equals(expected, result, len, tol)
end subroutine test_lcg64_seed_skip_sequence


end module test_lcg64_basic_rng