! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module test_lcg32_basic_rng

use constants

use fruit, only: assert_equals, run_test_case, test_set_summary
use lcg32_basic_rng, only: Lcg32BasicRng, init_lcg32_brng
use mpi, only: comm_world

implicit none

private
public :: test_lcg32_basic_rng_set

contains

subroutine test_lcg32_basic_rng_set
  call run_test_case(test_lcg32_seed_sequence, "check whether lcg32 sequence is correct for a random seed")
  call run_test_case(test_lcg32_seed_skip_sequence, "check whether lcg32 sequence is correct after skipping")

  call test_set_summary("src/random/basic_rng/lcg32_basic_rng.f90")
end subroutine test_lcg32_basic_rng_set

subroutine test_lcg32_seed_sequence
  integer, parameter     :: len = 50
  real(sp), parameter    :: tol = 0.0001_sp
  integer(i8), parameter :: seed = 67773761_i8
  real(sp)               :: expected(len), result(len)
  type(Lcg32BasicRng)    :: brng

  !obtained from the ref. code that completely follows implementation here
  expected = [0.018963_sp, 0.996737_sp, 0.787193_sp, 0.150028_sp, 0.510439_sp, &
              0.413460_sp, 0.793932_sp, 0.384485_sp, 0.377907_sp, 0.511545_sp, &
              0.652976_sp, 0.050314_sp, 0.713354_sp, 0.564712_sp, 0.411202_sp, &
              0.820192_sp, 0.306006_sp, 0.820178_sp, 0.669508_sp, 0.819510_sp, &
              0.871761_sp, 0.723942_sp, 0.099839_sp, 0.164650_sp, 0.510154_sp, &
              0.033757_sp, 0.150693_sp, 0.193629_sp, 0.373247_sp, 0.741680_sp, &
              0.550481_sp, 0.436010_sp, 0.780799_sp, 0.978061_sp, 0.468110_sp, &
              0.618949_sp, 0.266105_sp, 0.318859_sp, 0.440162_sp, 0.256241_sp, &
              0.286870_sp, 0.168143_sp, 0.424064_sp, 0.617252_sp, 0.585504_sp, &
              0.137675_sp, 0.586596_sp, 0.122616_sp, 0.433196_sp, 0.595431_sp]

  if (.not. comm_world%initialized) call comm_world%init()
  brng = Lcg32BasicRng(seed=seed, comm=comm_world, mpi_distinct=.true., omp_distinct=.false.)
  call brng%rand(result)

  call assert_equals(expected, result, len, tol)
end subroutine test_lcg32_seed_sequence

subroutine test_lcg32_seed_skip_sequence
  integer, parameter     :: len=50, skip=200
  real(dp), parameter    :: tol = 0.0001_dp
  integer(i8), parameter :: seed = 212821_i8
  real(dp)               :: expected(len), result(len)
  type(Lcg32BasicRng)    :: brng

  !obtained from the ref. code that completely follows implementation here
  expected = [0.280969_dp, 0.786253_dp, 0.247276_dp, 0.572577_dp, 0.953909_dp, &
              0.493972_dp, 0.606844_dp, 0.207877_dp, 0.996355_dp, 0.507294_dp, &
              0.426425_dp, 0.761266_dp, 0.259397_dp, 0.224566_dp, 0.771554_dp, &
              0.726248_dp, 0.024096_dp, 0.040647_dp, 0.173711_dp, 0.593050_dp, &
              0.140088_dp, 0.368083_dp, 0.376844_dp, 0.784354_dp, 0.719859_dp, &
              0.186611_dp, 0.475490_dp, 0.299030_dp, 0.123063_dp, 0.792847_dp, &
              0.803876_dp, 0.860762_dp, 0.215247_dp, 0.672058_dp, 0.185367_dp, &
              0.618423_dp, 0.135369_dp, 0.663870_dp, 0.401647_dp, 0.695088_dp, &
              0.056479_dp, 0.843797_dp, 0.990994_dp, 0.132549_dp, 0.673938_dp, &
              0.792470_dp, 0.735036_dp, 0.028226_dp, 0.405554_dp, 0.174429_dp]

  if (.not. comm_world%initialized) call comm_world%init()
  brng = Lcg32BasicRng(seed=seed, comm=comm_world, mpi_distinct=.true., omp_distinct=.false.)
  call brng%skip(skip)
  call brng%rand(result)

  call assert_equals(expected, result, len, tol)
end subroutine test_lcg32_seed_skip_sequence

end module test_lcg32_basic_rng