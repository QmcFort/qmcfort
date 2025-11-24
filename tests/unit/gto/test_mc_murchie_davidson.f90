! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module test_mc_murchie_davidson

use mc_murchie_davidson

use constants, only: wp
use fruit, only: assert_equals, run_test_case, test_set_summary

implicit none

private
public :: test_mc_murchie_davidson_set

contains

subroutine test_mc_murchie_davidson_set
  call run_test_case(test_overlap_exp_coeff, "Compare new routine to calculate hermite expansion coefficients")
  call run_test_case(test_overlap_exp_coeff_diff, "Compare overlap expansion coefficients for same parameters but diff. l")
  call run_test_case(test_kinetic_exp_coeff, "Compare old and new routines to calculate kinetic expansion coefficients")
  call run_test_case(test_overlap_exp_coeff_table, "Test overlap_exp_coeff for whole shell with la = lb = 3")
  call run_test_case(test_hermite_integrals_shell, "Test values obtained via hermite_integral_shell for l=5")
  call run_test_case(test_hermite_integrals_shell_full, "Test hermite_integral_shell for whole shell with l=8")

  call test_set_summary("src/gto/mc_murchie_davidson.f90")
end subroutine test_mc_murchie_davidson_set

subroutine test_overlap_exp_coeff
  integer, parameter    :: n = 3
  integer               :: la, lb
  real(wp)              :: a, b
  real(wp)              :: Ax, Bx
  real(wp), allocatable :: E_ijt(:,:,:)
  real(wp)              :: expected(n), result_(n)
  integer               :: i1, j1, t1
  integer               :: i2, j2, t2
  integer               :: i3, j3, t3
  real(wp), parameter   :: tol = 1.0e-06_wp

  la = 3; lb = 2
  a = 1.00_wp; b = 0.75_wp
  Ax = 0.0_wp; Bx = 1.0_wp
  
  i1=3; j1=2; t1=0
  i2=1; j2=1; t2=1
  i3=2; j3=1; t3=1

  expected(1) = overlap_exp(i1, j1, t1, a, b, Ax, Bx)
  expected(2) = overlap_exp(i2, j2, t2, a, b, Ax, Bx)
  expected(3) = overlap_exp(i3, j3, t3, a, b, Ax, Bx)

  call overlap_exp_coeff(la, lb, a, b, Ax, Bx, E_ijt)
  result_(1) = E_ijt(i1, j1, t1)
  result_(2) = E_ijt(i2, j2, t2)
  result_(3) = E_ijt(i3, j3, t3)

  call assert_equals(expected, result_, n, tol)
end subroutine test_overlap_exp_coeff

subroutine test_overlap_exp_coeff_diff
  integer, parameter    :: n = 3
  integer               :: la, lb
  real(wp)              :: a, b
  real(wp)              :: Ax, Bx
  real(wp), allocatable :: E_ijt(:,:,:), E_ijt_new(:,:,:)
  real(wp)              :: expected(n), result_(n)
  integer               :: i1, j1, t1
  integer               :: i2, j2, t2
  integer               :: i3, j3, t3
  real(wp), parameter   :: tol = 1.0e-06_wp

  la = 4; lb = 2
  a = 1.00_wp; b = 1.0_wp
  Ax = 0.0_wp; Bx = 0.0_wp
  
  i1=1; j1=0; t1=1
  i2=2; j2=0; t2=0
  i3=2; j3=1; t3=1

  call overlap_exp_coeff(la, lb, a, b, Ax, Bx, E_ijt)
  expected(1) = E_ijt(i1, j1, t1)
  expected(2) = E_ijt(i2, j2, t2)
  expected(3) = E_ijt(i3, j3, t3)

  call overlap_exp_coeff(la+2, lb, a, b, Ax, Bx, E_ijt_new)
  result_(1) = E_ijt_new(i1, j1, t1)
  result_(2) = E_ijt_new(i2, j2, t2)
  result_(3) = E_ijt_new(i3, j3, t3)

  call assert_equals(expected, result_, n, tol)
end subroutine test_overlap_exp_coeff_diff

subroutine test_overlap_exp_coeff_table
  integer               :: la, lb
  real(wp)              :: a, b
  real(wp)              :: Ra(3), Rb(3)
  real(wp), allocatable :: E_ijt(:,:,:), E_klu(:,:,:), E_mnv(:,:,:)
  real(wp), allocatable :: E_ijt_(:,:,:), E_klu_(:,:,:), E_mnv_(:,:,:)
  real(wp), allocatable :: temp_ijt(:), temp_klu(:), temp_mnv(:)
  real(wp), allocatable :: temp_ijt_(:), temp_klu_(:), temp_mnv_(:)
  real(wp), parameter   :: tol = 1.0e-06_wp
  integer               :: i, j, t

  la = 3; lb = 3
  a = 1.428_wp; b = 1.428_wp
  Ra = [0.0_wp, 0.0_wp, 2.19232_wp]; Rb = [0.0_wp, 0.0_wp, 2.19232_wp]
  
  call overlap_exp_coeff(la, lb, a, b, Ra(1), Rb(1), E_ijt)
  call overlap_exp_coeff(la, lb, a, b, Ra(2), Rb(2), E_klu)
  call overlap_exp_coeff(la, lb, a, b, Ra(3), Rb(3), E_mnv)

  allocate(E_ijt_(0:la, 0:lb, 0:la+lb))
  allocate(E_klu_(0:la, 0:lb, 0:la+lb))
  allocate(E_mnv_(0:la, 0:lb, 0:la+lb))

  do t = 0, la+lb
    do j = 0, lb 
      do i = 0, la
        E_ijt_(i,j,t) = overlap_exp(i, j, t, a, b, Ra(1), Rb(1))
        E_klu_(i,j,t) = overlap_exp(i, j, t, a, b, Ra(2), Rb(2))
        E_mnv_(i,j,t) = overlap_exp(i, j, t, a, b, Ra(3), Rb(3))
      end do
    end do
  end do

  allocate(temp_ijt(size(E_ijt)), temp_ijt_(size(E_ijt_)))
  allocate(temp_klu(size(E_klu)), temp_klu_(size(E_klu_)))
  allocate(temp_mnv(size(E_mnv)), temp_mnv_(size(E_mnv_)))
  temp_ijt = reshape(E_ijt, shape=[size(E_ijt)])
  temp_ijt_ = reshape(E_ijt_, shape=[size(E_ijt_)])
  temp_klu = reshape(E_klu, shape=[size(E_klu)])
  temp_klu_ = reshape(E_klu_, shape=[size(E_klu_)])
  temp_mnv = reshape(E_mnv, shape=[size(E_mnv)])
  temp_mnv_ = reshape(E_mnv_, shape=[size(E_mnv_)])

  call assert_equals(temp_ijt, temp_ijt_, size(temp_ijt), tol)
  call assert_equals(temp_klu, temp_klu_, size(temp_klu), tol)
  call assert_equals(temp_mnv, temp_mnv_, size(temp_mnv), tol)
end subroutine test_overlap_exp_coeff_table

subroutine test_kinetic_exp_coeff
  integer               :: la, lb, i, j
  real(wp)              :: a, b
  real(wp)              :: Ax, Bx
  real(wp), allocatable :: T_ij(:,:), T_ij_(:,:)
  real(wp), allocatable :: E_ijt(:,:,:)
  real(wp), parameter   :: tol = 1.0e-06_wp

  la = 3
  lb = 4
  a = 0.800_wp
  b = 1.125_wp
  Ax = 0.0_wp
  Bx = 1.0_wp

  call overlap_exp_coeff(la, lb+2, a, b, Ax, Bx, E_ijt)
  call kin_exp_coeff(la, lb, a, b, E_ijt, T_ij)

  allocate(T_ij_(0:la,0:lb))
  T_ij_ = 0.0_wp

  do j = 0, lb 
    do i = 0, la
      T_ij_(i,j) = kinetic_exp(i, j, a, b, Ax, Bx)
    end do
  end do

  call assert_equals(T_ij,  T_ij_, la+1, lb+1, tol)
end subroutine test_kinetic_exp_coeff

subroutine test_hermite_integrals_shell
  integer, parameter    :: n = 5
  integer               :: l
  real(wp)              :: p
  real(wp)              :: R(3)
  real(wp), allocatable :: R_tuv(:,:,:)
  real(wp)              :: expected(n), result_(n)
  integer               :: t1, u1, v1
  integer               :: t2, u2, v2
  integer               :: t3, u3, v3
  integer               :: t4, u4, v4
  integer               :: t5, u5, v5
  real(wp), parameter   :: tol = 1.0e-06_wp

  l = 5
  p = 1.00_wp
  R = [0.5_wp, 0.5_wp, 1.0_wp]
  
  t1=5; u1=0; v1=0
  t2=1; u2=0; v2=1
  t3=3; u3=0; v3=2
  t4=3; u4=1; v4=1
  t5=1; u5=2; v5=2

  expected(1) = hermite_integral_rec(t1, u1, v1, 0, p, R) 
  expected(2) = hermite_integral_rec(t2, u2, v2, 0, p, R) 
  expected(3) = hermite_integral_rec(t3, u3, v3, 0, p, R) 
  expected(4) = hermite_integral_rec(t4, u4, v4, 0, p, R) 
  expected(5) = hermite_integral_rec(t5, u5, v5, 0, p, R) 

  call hermite_integral_shell(l, p, R, R_tuv)
  result_(1) = R_tuv(t1, u1, v1)
  result_(2) = R_tuv(t2, u2, v2)
  result_(3) = R_tuv(t3, u3, v3)
  result_(4) = R_tuv(t4, u4, v4)
  result_(5) = R_tuv(t5, u5, v5)

  call assert_equals(expected, result_, n, tol)
end subroutine test_hermite_integrals_shell

subroutine test_hermite_integrals_shell_full
  integer               :: l, t, u, v
  real(wp)              :: p
  real(wp)              :: R(3)
  real(wp), allocatable :: R_tuv(:,:,:), R_tuv_(:,:,:)
  real(wp), allocatable :: expected(:), result_(:)
  real(wp), parameter   :: tol = 1.0e-06_wp

  l = 8
  p = 14250.00_wp
  R = [0.0_wp, 0.0_wp, 0.0_wp]

  call hermite_integral_shell(l, p, R, R_tuv)
  
  allocate(R_tuv_(0:l, 0:l, 0:l))
  R_tuv_ = 0.0_wp
  do v = 0, l
    do u = 0, l
      do t = 0, l
        if (t+u+v <= l) then
          R_tuv_(t,u,v) = hermite_integral_rec(t, u, v, 0, p, R)
        end if
      end do
    end do
  end do

  allocate(expected(size(R_tuv)))
  expected = reshape(R_tuv, shape=[size(R_tuv)])

  allocate(result_(size(R_tuv_)))
  result_ = reshape(R_tuv_, shape=[size(R_tuv_)])

  call assert_equals(expected, result_, size(expected), tol)
end subroutine test_hermite_integrals_shell_full

end module test_mc_murchie_davidson