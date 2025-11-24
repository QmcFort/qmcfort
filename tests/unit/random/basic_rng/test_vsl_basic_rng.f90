! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module test_vsl_basic_rng

use constants, only: dp
use fruit, only: assert_equals, run_test_case, test_set_summary
use mpi, only: comm_world

implicit none

private
public :: test_vsl_basic_rng_set

contains

subroutine test_vsl_basic_rng_set
#ifdef MKL
  call run_test_case(test_vsl_mt19937_skip_sequence, "check whether skip works wirh vsl_brng_mt19937")
  call run_test_case(test_vsl_r250_skip_sequence, "check whether skip works wirh vsl_brng_r250")

  call test_set_summary("src/random/basic_rng/vsl_basic_rng.f90")
#endif
end subroutine test_vsl_basic_rng_set

#ifdef MKL
subroutine test_vsl_mt19937_skip_sequence
  use vsl_basic_rng, only: VslBasicRng
  integer, parameter  :: len=50, skip=2000
  real(dp), parameter :: tol = 0.0001_dp
  real(dp)            :: expected(len), result(len), temp(skip)
  type(VslBasicRng)   :: brng

  if (.not. comm_world%initialized) call comm_world%init()
  brng = VslBasicRng(name="vsl_brng_mt19937", comm=comm_world, mpi_distinct=.true., omp_distinct=.false.)

  call brng%rand(temp)
  call brng%rand(expected)

  call brng%reset()

  call brng%skip(skip)
  call brng%rand(result)

  call assert_equals(expected, result, len, tol)
end subroutine test_vsl_mt19937_skip_sequence

subroutine test_vsl_r250_skip_sequence
  use vsl_basic_rng, only: VslBasicRng
  integer, parameter  :: len=50, skip=2000
  real(dp), parameter :: tol = 0.0001_dp
  real(dp)            :: expected(len), result(len), temp(skip)
  type(VslBasicRng)   :: brng

  if (.not. comm_world%initialized) call comm_world%init()
  brng = VslBasicRng(name="vsl_brng_r250", comm=comm_world, mpi_distinct=.true., omp_distinct=.false.)

  call brng%rand(temp)
  call brng%rand(expected)

  call brng%reset()

  call brng%skip(skip)
  call brng%rand(result)

  call assert_equals(expected, result, len, tol)
end subroutine test_vsl_r250_skip_sequence
#endif

end module test_vsl_basic_rng