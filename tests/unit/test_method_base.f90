! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module test_method_base

use constants, only: charlen
use fruit, only: assert_equals, run_test_case, test_set_summary
use method_base, only: QmcFortMethod

implicit none

private
public :: test_method_base_set

contains

subroutine test_method_base_set
  call run_test_case(test_method_get_activity, "check activity of the defined method")
  call run_test_case(test_method_get_basis, "check basis of the defined method")

  call test_set_summary("src/method_base.f90")
end subroutine test_method_base_set

subroutine test_method_get_activity
  type(QmcFortMethod) :: method
  logical             :: expected
  logical             :: result_
  
  method = QmcFortMethod(method="afqmc", def_active=.false., basis="mean_field", integral_mode="cholesky")
  
  expected = .false.
  result_ = method%is_active()
  call assert_equals(expected, result_)
end subroutine test_method_get_activity

subroutine test_method_get_basis
  type(QmcFortMethod)           :: method
  character(len=:), allocatable :: expected
  character(len=:), allocatable :: result_
  
  method = QmcFortMethod(method="hfock", def_active=.true., basis="gauss", integral_mode="eri")

  expected = "gauss"
  result_ = method%get_basis()
  call assert_equals(expected, result_)
end subroutine test_method_get_basis

end module test_method_base
