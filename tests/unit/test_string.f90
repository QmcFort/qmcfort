! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module test_string

use string

use constants, only: charlen
use fruit, only: assert_equals, run_test_case, test_set_summary

implicit none

private
public :: test_string_set

contains

subroutine test_string_set
  call run_test_case(test_lowercase_small_letters, "lowercase applied to the string with small letters only")
  call run_test_case(test_lowercase_capital_letters, "lowercase applied to the string with capital letters")
  call run_test_case(test_lowercase_mixed_letters, "lowercase applied to the string with small and capital letters")
  call run_test_case(test_lowercase_numbers, "lowercase applied to the string with numbers")
  call run_test_case(test_lowercase_general, "lowercase applied to the string with all possible characters")
  call run_test_case(test_str_integer, "Test str - convert integer into character")
  call run_test_case(test_split_tag_equal, "split tag - string that contains = character")
  call run_test_case(test_split_tag_space, "split tag - string that contains space character")
  call run_test_case(test_split_tag_longspace, "split tag - string that contains more than one space character")
  call run_test_case(test_isin, "check is the string part of the larger string")
  call run_test_case(test_isin_part, "check is the string part of the larger string")
  call run_test_case(test_extract_word_simple, "check first word from one word character")
  call run_test_case(test_extract_word, "check first word from multi word character")
  call run_test_case(test_rest_str_1, "extract rest of the string if matching is found 1")
  call run_test_case(test_rest_str_1, "extract rest of the string if matching is found 2")
  call run_test_case(test_read_int_arr_no_del, "read int array from string del=space")
  call run_test_case(test_read_int_arr_del, "read int array from string del=,")

  call test_set_summary("src/string.f90")
end subroutine test_string_set
  
subroutine test_lowercase_small_letters
  integer, parameter    :: strlen=50
  character(len=strlen) :: expected, result_
  !
  expected = "these are small letters"
  result_ = lowercase(expected)
  !
  call assert_equals(expected, result_)
end subroutine test_lowercase_small_letters

subroutine test_lowercase_capital_letters
  integer, parameter    :: strlen=50
  character(len=strlen) :: expected, result_
  !
  expected = "these are capital letters"
  result_  = "THESE ARE CAPITAL LETTERS"
  result_ = lowercase(result_)
  !
  call assert_equals(expected, result_)
end subroutine test_lowercase_capital_letters

subroutine test_lowercase_mixed_letters
  integer, parameter    :: strlen=50
  character(len=strlen) :: expected, result_
  !
  expected = "these are mixed letters"
  result_  = "These Are Mixed Letters"
  result_ = lowercase(result_)
  !
  call assert_equals(expected, result_)
end subroutine test_lowercase_mixed_letters

subroutine test_lowercase_numbers
  integer, parameter    :: strlen=50
  character(len=strlen) :: expected, result_
  !
  expected = "this is a string with numbers 12345"
  result_  = "This is a String with Numbers 12345"
  result_ = lowercase(result_)
  !
  call assert_equals(expected, result_)
end subroutine test_lowercase_numbers

subroutine test_lowercase_general
  integer, parameter    :: strlen = 50
  character(len=strlen) :: expected, result_
  !
  expected = "some.strange:string123"
  result_  = "Some.Strange:STRING123"
  result_ = lowercase(result_)
  !
  call assert_equals(expected, result_)
end subroutine test_lowercase_general

subroutine test_str_integer
  character(len=:), allocatable :: expected, result_
  integer                       :: ival
  !
  ival = 12345
  result_ = str(ival)
  allocate(character(len=len(result_)) :: expected)
  expected = "12345"
  !
  call assert_equals(expected, result_)
end subroutine test_str_integer

subroutine test_split_tag_equal
  integer, parameter    :: strlen = 50
  character(len=strlen) :: expected_1, expected_2, result_1, result_2
  !
  result_1   = "THIS IS KEY TO SPLIT = 100"
  expected_1 = "THIS IS KEY TO SPLIT"
  expected_2 = "100"
  call split_tag(result_1, result_2)
  !
  call assert_equals(expected_1, result_1)
  call assert_equals(expected_2, result_2)
end subroutine test_split_tag_equal

subroutine test_split_tag_space
  integer, parameter    :: strlen = 50
  character(len=strlen) :: expected_1, expected_2, result_1, result_2
  !
  result_1   = "MAXLEN 100"
  expected_1 = "MAXLEN"
  expected_2 = "100"
  call split_tag(result_1, result_2)
  !
  call assert_equals(expected_1, result_1)
  call assert_equals(expected_2, result_2)
end subroutine test_split_tag_space

subroutine test_split_tag_longspace
  integer, parameter    :: strlen = 50
  character(len=strlen) :: expected_1, expected_2, result_1, result_2
  !
  result_1   = "MAXLEN      100"
  expected_1 = "MAXLEN"
  expected_2 = "100"
  call split_tag(result_1, result_2)
  !
  call assert_equals(expected_1, result_1)
  call assert_equals(expected_2, result_2)
end subroutine test_split_tag_longspace

subroutine test_isin
  integer, parameter    :: strlen = 50
  character(len=strlen) :: a, b
  logical               :: expected, result_
  !
  a = "mean_field"
  b = "gauss_orth mean_field"
  !
  expected = .true.
  result_ = isin(a, b)
  !
  call assert_equals(expected, result_)
end subroutine test_isin

subroutine test_isin_part
  integer, parameter    :: strlen = 50
  character(len=strlen) :: a, b
  logical               :: expected, result_
  !
  a = "GAUSS"
  b = "gauss_orth mean_field"
  !
  expected = .false.
  result_ = isin(a, b)
  !
  call assert_equals(expected, result_)
end subroutine test_isin_part

subroutine test_extract_word_simple
  character(len=:), allocatable :: expected, result_
  !
  expected = "gauss_orth"
  result_ = extract_word(expected)
  expected = "gauss_orth"
  !
  call assert_equals(expected, result_)
end subroutine test_extract_word_simple  

subroutine test_extract_word
  integer, parameter    :: strlen = 50
  character(len=strlen) :: expected, result_
  !
  expected = "gauss_orth mean_field"
  result_ = extract_word(expected)
  expected = "gauss_orth"
  !
  call assert_equals(expected, result_)
end subroutine test_extract_word

subroutine test_rest_str_1
  integer, parameter    :: strlen = 50
  character(len=strlen) :: expected, result_
  expected = "'first': 1, 'second': 2, 'third': 3"
  result_ = rest_str("'third':", expected, sep=",")
  expected = "3"
  !
  call assert_equals(expected, result_)
end subroutine test_rest_str_1

subroutine test_rest_str_2
  integer, parameter    :: strlen = 50
  character(len=strlen) :: expected, result_
  expected = "'first': 1, 'second': 2, 'shape': (5 10 15"
  result_ = rest_str("'shape': (", expected, sep=",")
  expected = "5 10 15"
  !
  call assert_equals(expected, result_)
end subroutine test_rest_str_2

subroutine test_read_int_arr_no_del
  integer, parameter     :: n = 8
  integer                :: expected(n)
  integer, allocatable   :: result_(:)
  character(len=charlen) :: str

  expected = [1, 2, 3, 4, 5, 6, 7, 8]

  str = "   1 2   3 4 5   6  7  8   "
  result_ = int_arr_from_string(str)

  call assert_equals(expected, result_, n)
end subroutine test_read_int_arr_no_del

subroutine test_read_int_arr_del
  integer, parameter     :: n = 5
  integer                :: expected(n)
  integer, allocatable   :: result_(:)
  character(len=charlen) :: str

  expected = [12, 21, 30, -124, 100]

  str = " 12,21,  30,  -124,   100 "
  result_ = int_arr_from_string(str, del=",")

  call assert_equals(expected, result_, n)
end subroutine test_read_int_arr_del

end module test_string