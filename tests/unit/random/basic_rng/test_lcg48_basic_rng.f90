! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module test_lcg48_basic_rng

use constants

use fruit, only: assert_equals, run_test_case, test_set_summary
use lcg48_basic_rng, only: Lcg48BasicRng, init_lcg48_brng
use mpi, only: comm_world

implicit none

private
public :: test_lcg48_basic_rng_set

contains

subroutine test_lcg48_basic_rng_set
  call run_test_case(test_lcg48_seed_sequence, "check whether lcg48 sequence is correct for a random seed")
  call run_test_case(test_lcg48_seed_skip_sequence, "check whether lcg48 skip works correctly")

  call test_set_summary("src/random/basic_rng/lcg48_basic_rng.f90")
end subroutine test_lcg48_basic_rng_set

subroutine test_lcg48_seed_sequence
  integer, parameter     :: len = 50
  integer(i8), parameter :: seed = 1_i8
  real(dp), parameter    :: tol = 0.0001_dp
  real(dp)               :: expected(len), result(len)
  type(Lcg48BasicRng)    :: brng

  !obtained from the ref. code that completely follows implementation here
  expected = [0.841613_dp, 0.022735_dp, 0.586973_dp, 0.918193_dp, 0.795010_dp, &
              0.972252_dp, 0.931519_dp, 0.730262_dp, 0.223280_dp, 0.702691_dp, &
              0.047234_dp, 0.482635_dp, 0.015613_dp, 0.622457_dp, 0.824324_dp, &
              0.562704_dp, 0.256889_dp, 0.160260_dp, 0.092974_dp, 0.614637_dp, &
              0.483149_dp, 0.093872_dp, 0.521472_dp, 0.994047_dp, 0.285505_dp, &
              0.211581_dp, 0.251202_dp, 0.893536_dp, 0.626515_dp, 0.435187_dp, &
              0.640469_dp, 0.072723_dp, 0.047121_dp, 0.071322_dp, 0.339910_dp, &
              0.255871_dp, 0.828645_dp, 0.123651_dp, 0.009056_dp, 0.130757_dp, &
              0.214095_dp, 0.825505_dp, 0.363805_dp, 0.985250_dp, 0.207900_dp, &
              0.427318_dp, 0.991630_dp, 0.672344_dp, 0.375364_dp, 0.685396_dp]

  if (.not. comm_world%initialized) call comm_world%init()
  brng = Lcg48BasicRng(comm=comm_world, seed=seed, mpi_distinct=.true., omp_distinct=.false.)
  call brng%rand(result)

  call assert_equals(expected, result, len, tol)
end subroutine test_lcg48_seed_sequence

subroutine test_lcg48_seed_skip_sequence
  integer, parameter     :: len=50, skip=200
  integer(i8), parameter :: seed = 7327371_i8
  real(dp), parameter    :: tol = 0.0001_dp
  real(dp)               :: expected(len), result(len), temp(skip)
  type(Lcg48BasicRng)    :: brng

  if (.not. comm_world%initialized) call comm_world%init()
  brng = Lcg48BasicRng(seed=seed, comm=comm_world, mpi_distinct=.true., omp_distinct=.false.)

  call brng%rand(temp)
  call brng%rand(expected)

  call brng%reset()

  call brng%skip(skip)
  call brng%rand(result)

  call assert_equals(expected, result, len, tol)
end subroutine test_lcg48_seed_skip_sequence

end module test_lcg48_basic_rng