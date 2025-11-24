! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module test_hamilton_operations

use hamilton_operations

use constants, only: wp
use fruit, only: assert_equals, run_test_case, test_set_summary

implicit none

private
public :: test_hamilton_operations_set

contains

subroutine test_hamilton_operations_set
  call run_test_case(test_sum_gres_potential, "Test summation of the g-resolved potential")

  call test_set_summary("src/hamilton/hamilton_operations.f90")
end subroutine test_hamilton_operations_set

subroutine test_sum_gres_potential
  integer, parameter    :: n=10, ng=20, ispin=2
  integer               :: g, spin
  real(wp), parameter   :: tol=1.0e-04_wp
  real(wp), allocatable :: pot_gres(:,:,:,:), pot_exp(:,:,:), pot_res(:,:,:)

  allocate(pot_gres(n,n,ng,ispin))
  call random_number(pot_gres)

  allocate(pot_exp(n,n,ispin))
  pot_exp = 0.0_wp

  do g = 1, ng
    pot_exp = pot_exp + pot_gres(:,:,g,:)
  end do

  call sum_gres_potential(pot_gres, pot_res)

  do spin = 1, ispin
    call assert_equals(pot_exp(:,:,spin), pot_res(:,:,spin), n, n, tol)
  end do
end subroutine test_sum_gres_potential

end module test_hamilton_operations