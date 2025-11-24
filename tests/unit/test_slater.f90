! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module test_slater

use slater

use constants, only: wp
use fruit, only: assert_equals, run_test_case, test_set_summary

implicit none

private
public :: test_slater_set

contains

subroutine test_slater_set
  call run_test_case(test_str_len, "Test number of bits used to store Slater string")
  call run_test_case(test_str_exc1, "Single excitation of the Slater string")
  call run_test_case(test_str_content, "Content of the Slater string")
  call run_test_case(test_str_diff, "Difference between two strings")
  call run_test_case(test_str_phase, "Test phase between two Slater strings")
  call run_test_case(test_det_nel, "Test number of electrons in Slater determinant")
  call run_test_case(test_det_nmax, "Find highest occupied orbital in Slater determinant")
  call run_test_case(test_det_exc_ab, "Alpha-beta excitation of the Slater determinant")
  call run_test_case(test_det_exc_bb, "Beta-beta excitation of the Slater determinant")
  call run_test_case(test_det_content, "Content of the Slater determinant")
  call run_test_case(test_det_diff, "Difference between two Slater determinants")
  call run_test_case(test_det_phase_b, "Phase of the single beta excitation")
  call run_test_case(test_det_phase_ab, "Phase of the double alpha-beta excitation")
  call run_test_case(test_read_det_from_char, "Read Slater determinant from character string")

  call test_set_summary("src/slater.f90")
end subroutine test_slater_set

subroutine test_str_len
  type(slater_str) :: str
  integer          :: expected, result_

  str = slater_str(100_strkind)

  expected = 64
  result_ = str%len()

  call assert_equals(expected, result_)
end subroutine test_str_len

subroutine test_str_exc1
  type(slater_str) :: str, expected, result_
  integer          :: i, a
  
  i = 4
  a = 6
  
  str = slater_str(15_strkind)
  expected = slater_str(39_strkind)
  result_ = str%excite(i, a)
  
  call assert_equals(expected%str, result_%str)
end subroutine test_str_exc1

subroutine test_str_content
  type(slater_str)     :: str
  integer, allocatable :: expected(:), result_(:)
  
  expected = [1, 3, 5, 7, 9]
  str = slater_str(341_strkind)
  result_ = str%content()
  
  call assert_equals(expected, result_, size(expected))
end subroutine test_str_content

subroutine test_str_diff
  type(slater_str)     :: str_1, str_2
  integer, allocatable :: diff(:,:), diff_(:,:)
  
  str_1 = slater_str(55_strkind)
  str_2 = slater_str(199_strkind)
  
  diff = reshape([5, 7, 6, 8], shape=[2,2])
  diff_ = str_1%diff(str_2)
  
  call assert_equals(diff, diff_, size(diff,1), size(diff,2))
end subroutine test_str_diff

subroutine test_str_phase
  type(slater_str) :: str_1, str_2
  integer          :: expected, result_
  
  str_1 = slater_str(1465_strkind)
  str_2 = slater_str(3992_strkind)

  expected = -1
  result_ = str_1%phase(str_2)

  call assert_equals(expected, result_)
end subroutine test_str_phase

subroutine test_det_nel
  type(slater_det) :: det
  integer          :: expected(2), result_(2)

  det = slater_det(25_strkind, 99_strkind)
  expected = [3, 4]
  result_ = det%nel()

  call assert_equals(expected, result_, 2)
end subroutine test_det_nel

subroutine test_det_nmax
  type(slater_det) :: det
  integer          :: expected, result_

  det = slater_det(705_strkind, 171_strkind)
  expected = 10
  result_ = det%nmax()

  call assert_equals(expected, result_)
end subroutine test_det_nmax

subroutine test_det_exc_ab
  type(slater_det) :: expected, result_
  integer          :: i, a, s1, j, b, s2
  
  i = 3; a = 4; s1 = 1
  j = 3; b = 5; s2 = 2
  
  expected = slater_det(11_strkind, 19_strkind)
  result_ = slater_det(7_strkind, 7_strkind)
  result_ = result_%excite(i, a, s1, j, b, s2)
  
  call assert_equals(expected%str, result_%str, 2)
end subroutine test_det_exc_ab

subroutine test_det_exc_bb
  type(slater_det) :: expected, result_
  integer          :: i, a, s1, j, b, s2
  
  i = 2; a = 5; s1 = 2
  j = 3; b = 6; s2 = 2
  
  expected = slater_det(7_strkind, 49_strkind)
  result_ = slater_det(7_strkind, 7_strkind)
  result_ = result_%excite(i, a, s1, j, b, s2)
  
  call assert_equals(expected%str, result_%str, 2)
end subroutine test_det_exc_bb

subroutine test_det_content
  type(slater_det)     :: det
  integer, allocatable :: expected(:,:), result_(:,:)
  integer, parameter   :: n=5
  
  expected = reshape([1, 2, 3, 4, 6, 1, 1, 1, 2, 2], shape=[n, 2])
  det = slater_det(7_strkind, 40_strkind)
  result_ = det%content()
  
  call assert_equals(expected, result_, size(expected,1), size(expected,2))
end subroutine test_det_content

subroutine test_det_diff
  type(slater_det)     :: det_1, det_2
  integer, allocatable :: diff(:,:), diff_(:,:)
  
  det_1 = slater_det(13_strkind, 11_strkind)
  det_2 = slater_det(41_strkind, 19_strkind)
  
  diff = reshape([3, 6, 1, 4, 5, 2], shape=[3,2])
  diff_ = det_1%diff(det_2)
  
  call assert_equals(diff, diff_, size(diff,1), size(diff,2))
end subroutine test_det_diff

subroutine test_det_phase_b
  type(slater_det) :: det_1, det_2
  integer          :: expected, result_

  det_1 = slater_det(31_strkind, 597_strkind)
  det_2 = slater_det(31_strkind, 724_strkind)

  expected = -1
  result_ = det_1%phase(det_2)

  call assert_equals(expected, result_)
end subroutine test_det_phase_b

subroutine test_det_phase_ab
  type(slater_det) :: det_1, det_2
  integer          :: expected, result_

  det_1 = slater_det(103_strkind, 31_strkind)
  det_2 = slater_det(199_strkind, 155_strkind)

  expected = -1
  result_ = det_1%phase(det_2)

  call assert_equals(expected, result_)
end subroutine test_det_phase_ab

subroutine test_read_det_from_char
  type(slater_det)       :: det, det_
  complex(wp)            :: coeff, coeff_
  character(len=charlen) :: str
  real(wp), parameter    :: tol = 1.0e-06_wp

  det = slater_det(27_strkind, 39_strkind)
  coeff = (0.102080_wp, 0.001000_wp)

  str = " [0 1 3 4]  [ 0 1 2 5 ]  0.102080 0.001000"
  coeff_ = det_%read_from_char(str)

  call assert_equals(coeff, coeff_, tol)
  call assert_equals(det%str,  det_%str, 2)
end subroutine test_read_det_from_char

end module test_slater
