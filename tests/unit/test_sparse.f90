! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module test_sparse

use sparse

use constants, only: sp, dp, wp, spprec
use fruit, only: assert_equals, run_test_case, test_set_summary

implicit none

private
public :: test_sparse_set

contains

subroutine test_sparse_set
  call run_test_case(test_init_sparse_matrix_trivial, "Initialization of the sparse matrix - manual comparison")
  call run_test_case(test_init_sparse_matrix_complex_trivial, "Initialization of the complex sparse matrix")
  call run_test_case(test_init_sparse_matrix_zero, "Initialization of zero sparse matrix")
  call run_test_case(test_init_sparse_matrix_sprec, "Sprse matrix initialization and get full matrix for single precision")
  call run_test_case(test_assign_sparse_matrix, "Check assigment of the sparse matrices")
  call run_test_case(test_assign_sparse_matrix_complex, "Check assignment of the complex sparse matrices")
  call run_test_case(test_add_sparse_trivial, "Add two simple sparse matrices")
  call run_test_case(test_add_sparse_complex_trivial, "Add two simple complex sparse matrices")
  call run_test_case(test_add_sparse_rc_trivial, "Add two real and complex sparse matrices")
  call run_test_case(test_add_assign_sparse, "Test a = a + b with sparse matrices")
  call run_test_case(test_add_sparse_random, "Add two random sparse matrices")
  call run_test_case(test_add_sparse_to_full, "Add sparse matrix to the full matrix")
  call run_test_case(test_add_sparse_to_full_c, "Add sparse matrix to the complex full matrix")
  call run_test_case(test_add_sparse_to_full_coeff, "Add sparse matrix to the full matrix with prefactor")
  call run_test_case(test_add_sparse_to_full_1d, "Add sparse matrix to the array of matrices")
  call run_test_case(test_sparse_left_scalar_sp, "Test multiplication of the sparse matrix by scalar from the left")
  call run_test_case(test_sparse_left_right_array, "Check consistency x * a = a * x for spasre matrices with array x")
  call run_test_case(test_sparse_left_right_cr_array, &
                     "Check consistency x * a = a * x for complex/real spasre matrices with array x")
  call run_test_case(test_sparse_add_mult, "Test c = a + x * b with sparse matrices")

  call test_set_summary("src/sparse.f90")
end subroutine test_sparse_set

subroutine test_init_sparse_matrix_trivial
  integer, parameter  :: i=5, j=3, k=3
  type(sparse_matrix) :: a_sp
  real(spprec)        :: a(i, j)
  real(spprec)        :: expected(k), expected_sparsity
  
  a = 0.0_dp
  a(1,1) = 1.0_dp
  a(1,3) = 1.0_dp
  a(5,3) = -1.0_dp
  expected(1) = 1.0_dp
  expected(2) = 1.0_dp
  expected(3) = -1.0_dp
  expected_sparsity = real(k, spprec) / real(i*j, spprec)
  
  a_sp = sparse_matrix(a)
  call assert_equals(expected, a_sp%val(1,:), k)
  call assert_equals(expected_sparsity, a_sp%sparsity)
end subroutine test_init_sparse_matrix_trivial

subroutine test_init_sparse_matrix_complex_trivial
  integer, parameter    :: i=5, j=3, k=3
  type(sparse_matrix_c) :: a_sp
  complex(spprec)       :: a(i, j)
  complex(spprec)       :: expected(k)
  real(spprec)          :: expected_sparsity
  
  a = (0.0_spprec, 0.0_spprec)
  a(1,1) = (1.0_spprec, 1.0_spprec)
  a(1,3) = 1.0_spprec
  a(5,3) = (-1.0_spprec, -1.0_spprec)
  expected(1) = (1.0_spprec, 1.0_spprec)
  expected(2) = (1.0_spprec, 0.0_spprec)
  expected(3) = (-1.0_spprec, -1.0_spprec)
  expected_sparsity = real(k, spprec) / real(i*j, spprec)
  
  a_sp = sparse_matrix_c(a)
  call assert_equals(expected, a_sp%val(1,:), k)
  call assert_equals(expected_sparsity, a_sp%sparsity)
end subroutine test_init_sparse_matrix_complex_trivial

subroutine test_init_sparse_matrix_zero
  integer, parameter        :: n=10
  real(spprec)              :: expected(n,n)
  real(spprec), allocatable :: result_(:,:)
  type(sparse_matrix)       :: a_
  
  expected = 0.0_spprec
  a_ = sparse_matrix(expected)
  call get_full_matrix(a_, result_)
  
  call assert_equals(expected, result_, n, n)
end subroutine test_init_sparse_matrix_zero

subroutine test_init_sparse_matrix_sprec
  integer, parameter        :: i=10, j=10
  real(spprec), parameter   :: tol=1.0e-04_spprec           
  type(sparse_matrix)       :: a_sp
  real(spprec)              :: a(i, j)
  real(spprec)              :: expected(i,j)
  real(spprec), allocatable :: result_(:,:)
  
  call random_number(a)
  where (abs(a) > tol)
    expected = a
  else where
    expected = 0.0_spprec
  end where
  
  a_sp = sparse_matrix(a, tol, 1)
  call get_full_matrix(a_sp, result_)
  call assert_equals(expected, result_, i, j)
end subroutine test_init_sparse_matrix_sprec

subroutine test_assign_sparse_matrix
  integer, parameter      :: n = 10
  real(spprec), parameter :: tol=1.0e-04_spprec
  type(sparse_matrix)     :: a, b
  real(spprec)            :: matrix(n, n)
  
  call random_number(matrix)
  a = sparse_matrix(matrix, tol)
  b = a 
  call assert_equals(a%size_, b%size_)
  call assert_equals(a%indx, b%indx, a%size_)
  call assert_equals(a%val, b%val, a%width, a%size_)
end subroutine test_assign_sparse_matrix

subroutine test_assign_sparse_matrix_complex
  integer, parameter      :: n = 10
  real(spprec), parameter :: tol=1.0e-04_spprec
  type(sparse_matrix_c)   :: a, b
  complex(spprec)         :: matrix(n,n)
  real(spprec)            :: matrix_r(n,n), matrix_i(n,n)
  !
  call random_number(matrix_r)
  call random_number(matrix_i)
  matrix = cmplx(matrix_r, matrix_i)
  !
  a = sparse_matrix_c(matrix, tol)
  b = a 
  call assert_equals(a%size_, b%size_)
  call assert_equals(a%indx, b%indx, a%size_)
  call assert_equals(a%val, b%val, a%width, a%size_)
end subroutine test_assign_sparse_matrix_complex

subroutine test_add_sparse_trivial
  integer, parameter      :: n=3, m=2
  real(spprec), parameter :: tol = 1.0e-09_dp
  real(spprec)            :: a(n,m), b(n,m)
  real(spprec)            :: expected(n,m)
  real(dp), allocatable   :: result_(:,:)
  type(sparse_matrix)     :: a_, b_, c_
  !
  a = reshape([1.0_spprec, 0.0_spprec, 0.0_spprec, 0.0_spprec, 0.0_spprec,  0.0_spprec], shape=[n, m]) 
  b = reshape([1.0_spprec, 0.0_spprec, 0.0_spprec, 3.0_spprec, 0.0_spprec, -1.0_spprec], shape=[n, m]) 
  expected = a + b
  a_ = sparse_matrix(a, tol)
  b_ = sparse_matrix(b, tol)
  !
  c_ = a_ + b_  
  call get_full_matrix(c_, result_)  
  !
  call assert_equals(expected, result_, n, m)
end subroutine test_add_sparse_trivial

subroutine test_add_sparse_complex_trivial
  integer, parameter           :: n=2, m=2
  real(spprec), parameter      :: tol = 1.0e-09_spprec
  complex(spprec)              :: a(n,m), b(n,m)
  complex(spprec)              :: expected(n,m)
  complex(spprec), allocatable :: result_(:,:)
  type(sparse_matrix_c)        :: a_, b_, c_
  !
  a = reshape([(1.0_spprec,-1.0_spprec), (0.0_spprec,0.0_spprec), (0.0_spprec, 0.0_spprec), & 
              (-1.0_spprec, 0.0_spprec)], shape=[n, m]) 
  b = reshape([(1.0_spprec, 0.0_spprec), (0.0_spprec, 3.0_spprec),(0.0_spprec, 0.0_spprec), & 
              (1.0_spprec, 0.0_spprec)], shape=[n, m]) 
  expected = a + b
  a_ = sparse_matrix_c(a, tol)
  b_ = sparse_matrix_c(b, tol)
  !
  c_ = a_ + b_  
  call get_full_matrix(c_, result_)  
  !
  call assert_equals(expected, result_, n, m)
end subroutine test_add_sparse_complex_trivial

subroutine test_add_sparse_rc_trivial
  integer, parameter           :: n=2, m=2
  real(spprec), parameter      :: tol = 1.0e-09_spprec
  real(spprec)                 :: a(n,m)
  complex(spprec)              :: b(n,m)
  complex(spprec)              :: expected(n,m)
  complex(spprec), allocatable :: result_(:,:)
  type(sparse_matrix)          :: a_
  type(sparse_matrix_c)        :: b_, c_
  !
  a = reshape([1.0_spprec, 0.0_spprec, 0.0_spprec, -1.0_spprec], shape=[n, m]) 
  b = reshape([(1.0_spprec, 0.0_spprec), (0.0_spprec, 3.0_spprec),(0.0_spprec, 0.0_spprec),  &
              (1.0_spprec, 0.0_spprec)], shape=[n, m]) 
  expected = a + b
  a_ = sparse_matrix(a, tol)
  b_ = sparse_matrix_c(b, tol)
  !
  c_ = a_ + b_  
  call get_full_matrix(c_, result_)  
  !
  call assert_equals(expected, result_, n, m)
end subroutine test_add_sparse_rc_trivial

subroutine test_add_assign_sparse
  integer, parameter        :: n=3, m=2
  real(spprec), parameter   :: tol = 1.0e-09_spprec
  real(spprec)              :: a(n,m), b(n,m)
  real(spprec)              :: expected(n,m)
  real(spprec), allocatable :: result_(:,:)
  type(sparse_matrix)       :: a_, b_
  !
  a = reshape([1.0_spprec, 0.0_spprec, 0.0_spprec, 0.0_spprec, 0.0_spprec,  3.0_spprec], shape=[n, m]) 
  b = reshape([1.0_spprec, 0.0_spprec, 0.0_spprec, 3.0_spprec, 0.0_spprec, -1.0_spprec], shape=[n, m]) 
  expected = a + b
  a_ = sparse_matrix(a, tol)
  b_ = sparse_matrix(b, tol)
  !
  a_ = a_ + b_  
  call get_full_matrix(a_, result_)  
  !
  call assert_equals(expected, result_, n, m)
end subroutine test_add_assign_sparse

subroutine test_add_sparse_random
  integer, parameter        :: n = 20
  real(spprec), parameter   :: tol=1.0e-03_spprec
  real(spprec)              :: a(n,n), b(n,n) 
  real(spprec)              :: expected(n,n)
  real(spprec), allocatable :: result_(:,:)
  type(sparse_matrix)       :: a_, b_, c_
  !
  call random_number(a)
  call random_number(b)
  where (abs(a) < tol)
    a = 0.0_spprec
  end where
  where (abs(b) < tol)
    b = 0.0_spprec
  end where
  expected = a + b
  !
  a_ = sparse_matrix(a, tol)
  b_ = sparse_matrix(b, tol)
  c_ = a_ + b_
  call get_full_matrix(c_, result_)
  !
  call assert_equals(expected, result_, n, n)
end subroutine test_add_sparse_random

subroutine test_add_sparse_to_full
  integer, parameter  :: n=2, m=3
  real(spprec)        :: a(n,m), expected(n,m), result_(n,m)
  type(sparse_matrix) :: a_
  !
  a = reshape([1.0_spprec, 0.0_spprec, -2.0_spprec, 0.0_spprec, 0.0_spprec, -1.0_spprec], shape=[n, m])
  expected  = reshape([4.0_spprec, 3.0_spprec,  2.0_spprec, 1.0_spprec, 0.0_spprec, -1.0_spprec], shape=[n, m])
  result_   = reshape([4.0_spprec, 3.0_spprec,  2.0_spprec, 1.0_spprec, 0.0_spprec, -1.0_spprec], shape=[n, m])
  !
  expected = expected + a
  a_ = sparse_matrix(a)
  call add_to(result_, a_)
  !
  call assert_equals(expected, result_, n, m)
end subroutine test_add_sparse_to_full

subroutine test_add_sparse_to_full_c
  integer, parameter  :: n=2, m=2
  real(spprec)        :: a(n,m)
  complex(spprec)     :: expected(n,m), result_(n,m)
  complex(spprec)     :: coeff
  type(sparse_matrix) :: a_
  !
  coeff = (1.0_spprec, 1.0_spprec)
  a = reshape([1.0_spprec, 0.0_spprec, -2.0_spprec, 0.0_spprec], shape=[n, m])
  expected  = reshape([4.0_spprec, 3.0_spprec,  2.0_spprec, 1.0_spprec], shape=[n, m])
  result_   = reshape([4.0_spprec, 3.0_spprec,  2.0_spprec, 1.0_spprec], shape=[n, m])
  !
  expected = expected + a * coeff
  a_ = sparse_matrix(a)
  call add_to(result_, a_, coeff)
  !
  call assert_equals(expected, result_, n, m)
end subroutine test_add_sparse_to_full_c

subroutine test_add_sparse_to_full_coeff
  integer, parameter  :: n=2, m=3
  real(spprec)        :: a(n,m), expected(n,m), result_(n,m)
  real(spprec)        :: coeff
  type(sparse_matrix) :: a_
  !
  coeff = 5.0_spprec
  a = reshape([1.0_spprec, 0.0_spprec, -2.0_spprec, 0.0_spprec, 0.0_spprec, -1.0_spprec], shape=[n, m])
  expected  = reshape([4.0_spprec, 3.0_spprec,  2.0_spprec, 1.0_spprec, 0.0_spprec, -1.0_spprec], shape=[n, m])
  result_   = reshape([4.0_spprec, 3.0_spprec,  2.0_spprec, 1.0_spprec, 0.0_spprec, -1.0_spprec], shape=[n, m])
  !
  expected = expected + coeff*a
  a_ = sparse_matrix(a)
  call add_to(result_, a_, coeff)
  !
  call assert_equals(expected, result_, n, m)
end subroutine test_add_sparse_to_full_coeff

subroutine test_add_sparse_to_full_1d
  integer, parameter  :: n=10, m=5
  integer             :: i, indx
  real(spprec)        :: a(m,m), expected(n,m,m), result_(n,m,m)
  real(spprec)        :: coeff(n)
  type(sparse_matrix) :: a_
  !
  call random_number(coeff)
  call random_number(a)
  call random_number(expected)
  result_ = expected
  !
  do i = 1, n
    expected(i,:,:) = expected(i,:,:) + coeff(i) * a 
  end do
  !
  a_ = sparse_matrix(a)
  call add_to(result_, a_, coeff)
  !
  indx = max(n-3, 1)
  call assert_equals(expected(indx,:,:), result_(indx,:,:), m, m)
end subroutine test_add_sparse_to_full_1d

subroutine test_sparse_left_scalar_sp
  integer, parameter        :: n=5
  real(spprec), parameter   :: tol=1.0e-04_spprec, coeff=3.0_spprec
  real(spprec)              :: a(n,n)
  real(spprec)              :: expected(n,n)
  real(spprec), allocatable :: result_(:,:)
  type(sparse_matrix)       :: a_
  !
  call random_number(a)
  where (abs(a) >= tol)
    expected = coeff * a
  else where
    expected = 0.0_spprec
  end where
  !
  a_ = sparse_matrix(a, tol, 1)
  a_ = a_ * coeff
  call get_full_matrix(a_, result_)
  !
  call assert_equals(expected, result_, n, n) 
end subroutine test_sparse_left_scalar_sp

subroutine test_sparse_left_right_array
  integer, parameter      :: n=20, m=5
  real(spprec), parameter :: tol=1.0e-03_spprec
  type(sparse_matrix)     :: a, al, ar
  real(spprec)            :: coeff(m), mat(n,n)
  !
  call random_number(coeff)
  call random_number(mat)
  a = sparse_matrix(mat, tol)
  !
  al = coeff * a
  ar = a * coeff
  !
  call assert_equals(al%val, ar%val, m, a%size_) 
end subroutine test_sparse_left_right_array

subroutine test_sparse_left_right_cr_array
  integer, parameter      :: n=20, m=5
  real(spprec), parameter :: tol=1.0e-03_dp
  type(sparse_matrix)     :: a
  type(sparse_matrix_c)   :: al, ar
  real(spprec)            :: coeff_r(m), coeff_i(m), mat(n,n)
  complex(spprec)         :: coeff(m)
  !
  call random_number(coeff_r)
  call random_number(coeff_i)
  coeff = cmplx(coeff_r, coeff_i)
  call random_number(mat)
  a = sparse_matrix(mat, tol)
  !
  al = coeff * a
  ar = a * coeff
  !
  call assert_equals(al%val, ar%val, m, a%size_) 
end subroutine test_sparse_left_right_cr_array

subroutine test_sparse_add_mult
  integer, parameter        :: n=10
  real(spprec), parameter   :: tol=1.0e-04_spprec
  real(spprec)              :: a(n,n), b(n,n), coeff 
  real(spprec)              :: expected(n,n)
  real(spprec), allocatable :: result_(:,:)
  type(sparse_matrix)       :: a_, b_, c_
  !
  coeff = 3.0_spprec
  call random_number(a)
  call random_number(b)
  where (abs(a) < tol)
    a = 0.0_spprec
  end where
  where (abs(b) < tol)
    b = 0.0_spprec
  end where
  a_ = sparse_matrix(a, tol)
  b_ = sparse_matrix(b, tol)
  !
  expected = a + coeff * b
  c_ = a_ + coeff * b_
  call get_full_matrix(c_, result_)
  !
  call assert_equals(expected, result_, n, n)
end subroutine test_sparse_add_mult

end module test_sparse
