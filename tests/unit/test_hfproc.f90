! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module test_hfproc

use hfproc

use constants, only: wp
use fruit, only: assert_equals, run_test_case, test_set_summary

implicit none

private
public :: test_hfproc_set

contains

subroutine test_hfproc_set
  call run_test_case(test_get_degeneracies_p, "Find degeneracies in p states")
  call run_test_case(test_get_degeneracies_d, "Find degeneracies in d states")
  call run_test_case(test_get_degeneracies_d_pert, "Find degeneracies in slightly perturbed d states")
  call run_test_case(test_get_occ_groups_p, "Get bound for occupation numbers in 1s2s2p1 configuration")
  call run_test_case(test_get_occ_groups_pd, "Get bound for occupation numbers in 1s2s2p13d1 configuration")
  call run_test_case(test_get_occ_groups_empty, "Get bound for occupatien numbers with leading empty orbitals")
  call run_test_case(test_set_occupations_rhf, "Occupations of simple rhf model without degeneracies")
  call run_test_case(test_set_occupations_rhf_deg, "Occupations of simple rhf model with degeneracies")
  call run_test_case(test_set_occupations_spin_p1, "Occupations for 2p^1 configuration")
  call run_test_case(test_set_occupations_spin_p3, "Occupations for 2p^3 configuration")
  call run_test_case(test_set_occupations_spin_p5, "Occupations for 2p^5 configuration")
  call run_test_case(test_set_occupations_spin_p5_int, "Occupations for 2p^5 configuration with integer occ")
  call run_test_case(test_rdm_trafo, "Test spin-->trafo and trafo-->spin rdm transformation")
  call run_test_case(test_calculate_rdm_ispin1, "Calculate RDM for simple spin unpolarized system")
  call run_test_case(test_calculate_rdm_ispin2, "Calculate RDM for simple spin polarized system")

  call test_set_summary("src/hfproc.f90")
end subroutine test_hfproc_set

subroutine test_get_degeneracies_p
  integer, parameter :: n=6
  real(wp)           :: eigenval(n)
  integer            :: low, expected(2), result_(2)
  
  eigenval(1) = -1.0000_wp
  eigenval(2) = -0.8000_wp
  eigenval(3:5) = -0.5000_wp
  eigenval(6) = 0.2000_wp
  low = 3
  expected = [6, 3]
  
  call get_degeneracies(eigenval, low, result_(1), result_(2))
  
  call assert_equals(expected, result_, 2)
end subroutine test_get_degeneracies_p

subroutine test_get_degeneracies_d
  integer, parameter :: n=10
  real(wp)           :: eigenval(n)
  integer            :: low, expected(2), result_(2)
  
  eigenval(1) = -1.0000_wp
  eigenval(2) = -0.8000_wp
  eigenval(3:5) = -0.5000_wp
  eigenval(6:10) = 0.2000_wp
  low = 6
  expected = [11, 5]
  
  call get_degeneracies(eigenval, low, result_(1), result_(2))
  
  call assert_equals(expected, result_, 2)
end subroutine test_get_degeneracies_d

subroutine test_get_degeneracies_d_pert
  integer, parameter :: n=10
  real(wp)           :: eigenval(n)
  integer            :: low, expected(2), result_(2)
  
  eigenval(1) = -1.0000_wp
  eigenval(2) = -0.8000_wp
  eigenval(3:5) = -0.5000_wp
  eigenval(6) = 0.20000001_wp
  eigenval(7) = 0.20000002_wp
  eigenval(8) = 0.20000003_wp
  eigenval(9) = 0.20000004_wp
  eigenval(10) = 0.2000005_wp
  low = 6
  expected = [11, 5]
  
  call get_degeneracies(eigenval, low, result_(1), result_(2))
  
  call assert_equals(expected, result_, 2)
end subroutine test_get_degeneracies_d_pert

subroutine test_get_occ_groups_p
  integer, parameter   :: n=6
  real(wp), parameter  :: frac=0.33333_wp
  real(wp)             :: occ(n)
  integer, allocatable :: low(:), up(:)
  integer              :: low_(3), up_(3), ng_, ng
  
  occ = [1.0_wp, 1.0_wp, frac, frac, frac, 0.0_wp]
  ng_ = 3
  low_ = [1, 3, 6]
  up_ = [2, 5, 6]
  
  call get_occ_groups(occ, low, up, ng)
  
  call assert_equals(ng_, ng)
  call assert_equals(low_, low, ng)
  call assert_equals(up_, up, ng)
end subroutine test_get_occ_groups_p

subroutine test_get_occ_groups_pd
  integer, parameter   :: n=12
  real(wp), parameter  :: frac1=0.33333_wp, frac2=0.2_wp
  real(wp)             :: occ(n)
  integer, allocatable :: low(:), up(:)
  integer              :: low_(4), up_(4), ng_, ng
  
  occ = [1.0_wp, 1.0_wp, frac1, frac1, frac1, frac2, frac2, frac2, frac2, frac2, 0.0_wp, 0.0_wp]
  ng_ = 4
  low_ = [1, 3, 6, 11]
  up_ = [2, 5, 10, 12]
  
  call get_occ_groups(occ, low, up, ng)
  
  call assert_equals(ng_, ng)
  call assert_equals(low_, low, ng)
  call assert_equals(up_, up, ng)
end subroutine test_get_occ_groups_pd

subroutine test_get_occ_groups_empty
  integer, parameter   :: n=6
  real(wp)             :: occ(n)
  integer, allocatable :: low(:), up(:)
  integer              :: low_(2), up_(2), ng_, ng
  
  occ = [0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 1.0_wp, 1.0_wp]
  ng_ = 2
  low_ = [1, 4]
  up_ = [3, 6]
  
  call get_occ_groups(occ, low, up, ng)
  
  call assert_equals(ng_, ng)
  call assert_equals(low_, low, ng)
  call assert_equals(up_, up, ng)
end subroutine test_get_occ_groups_empty

subroutine test_set_occupations_rhf
  use hamilton_vars, only: hamil_descriptor
  integer, parameter     :: ispin=1, n=6
  real(wp), parameter    :: tol=1.0e-10_wp
  real(wp), parameter    :: frac=0.333333333333_wp
  type(hamil_descriptor) :: ham_des
  logical                :: lfrac
  real(wp), allocatable  :: eigenval(:,:)
  real(wp)               :: expected(n, 2)        
  real(wp), allocatable  :: result_(:,:)
  
  ham_des%n = n
  ham_des%nel = [3, 3]
  lfrac = .false.
  allocate(eigenval(n, ispin))
  eigenval(:,1) = [-1.0_wp, -0.8_wp, -0.5_wp, -0.3_wp, -0.1_wp, 0.0_wp]
  
  expected = 0.0_wp
  expected(1:3,1) = [1.0_wp, 1.0_wp, 1.0_wp]
  expected(1:3,2) = [1.0_wp, 1.0_wp, 1.0_wp]
  
  call set_occupations_eig(ham_des, result_, eigenval, lfrac)
  
  call assert_equals(expected, result_, n, ispin, tol)
end subroutine test_set_occupations_rhf

subroutine test_set_occupations_rhf_deg
  use hamilton_vars, only: hamil_descriptor
  integer, parameter     :: ispin=1, n=6
  real(wp), parameter    :: tol=1.0e-10_wp
  real(wp), parameter    :: frac=0.333333333333_wp
  type(hamil_descriptor) :: ham_des
  logical                :: lfrac
  real(wp), allocatable  :: eigenval(:,:)
  real(wp)               :: expected(n, 2)        
  real(wp), allocatable  :: result_(:,:)
  
  ham_des%n = n
  ham_des%nel = [3, 3]
  lfrac = .false.
  allocate(eigenval(n, ispin))
  eigenval(:,1) = [-1.0_wp, -0.8_wp, -0.5_wp, -0.5_wp, -0.5_wp, 0.0_wp]
  
  expected = 0.0_wp
  expected(1:4,1) = [1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp]
  expected(1:2,2) = [1.0_wp, 1.0_wp]
  
  call set_occupations_eig(ham_des, result_, eigenval, lfrac)
  
  call assert_equals(expected, result_, n, ispin, tol)
end subroutine test_set_occupations_rhf_deg

subroutine test_set_occupations_spin_p1
  use hamilton_vars, only: hamil_descriptor
  integer, parameter     :: ispin=2, n=8
  real(wp), parameter    :: tol=1.0e-10_wp
  real(wp), parameter    :: frac=0.333333333333_wp
  type(hamil_descriptor) :: ham_des
  logical                :: lfrac
  real(wp), allocatable  :: eigenval(:,:)
  real(wp)               :: expected(n, 2)        
  real(wp), allocatable  :: result_(:,:)
  
  ham_des%nel = [3, 2]
  lfrac = .true.
  allocate(eigenval(n, ispin))
  eigenval(:,1) = [-1.0_wp, -0.8_wp, -0.5_wp, -0.5_wp, -0.5_wp, -0.4_wp, -0.4_wp, -0.2_wp]
  eigenval(:,2) = [-1.0_wp, -0.8_wp, -0.5_wp, -0.5_wp, -0.5_wp, -0.4_wp, -0.4_wp, -0.2_wp]
  
  expected = 0.0_wp
  expected(1:5,1) = [1.0_wp, 1.0_wp, frac, frac, frac]
  expected(1:2,2) = [1.0_wp, 1.0_wp]
  
  call set_occupations_eig(ham_des, result_, eigenval, lfrac)
  
  call assert_equals(expected, result_, n, ispin, tol)
end subroutine test_set_occupations_spin_p1

subroutine test_set_occupations_spin_p3
  use hamilton_vars, only: hamil_descriptor
  integer, parameter     :: ispin=2, n=8
  real(wp), parameter    :: tol=1.0e-10_wp
  type(hamil_descriptor) :: ham_des
  logical                :: lfrac
  real(wp), allocatable  :: eigenval(:,:)
  real(wp)               :: expected(n, ispin)        
  real(wp), allocatable  :: result_(:,:)
  !
  ham_des%nel = [5, 2]
  lfrac = .true.
  allocate(eigenval(n, ispin))
  eigenval(:,1) = [-1.0_wp, -0.8_wp, -0.5_wp, -0.5_wp, -0.5_wp, -0.4_wp, -0.4_wp, -0.2_wp]
  eigenval(:,2) = [-1.0_wp, -0.8_wp, -0.5_wp, -0.5_wp, -0.5_wp, -0.4_wp, -0.4_wp, -0.2_wp]
  
  expected = 0.0_wp
  expected(1:5,1) = [1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp]
  expected(1:2,2) = [1.0_wp, 1.0_wp]
  
  call set_occupations_eig(ham_des, result_, eigenval, lfrac)
  
  call assert_equals(expected, result_, n, ispin, tol)
end subroutine test_set_occupations_spin_p3

subroutine test_set_occupations_spin_p5
  use hamilton_vars, only: hamil_descriptor
  integer, parameter     :: ispin=2, n=8
  real(wp), parameter    :: tol=1.0e-10_wp
  real(wp), parameter    :: frac=0.666666666666_wp
  type(hamil_descriptor) :: ham_des
  logical                :: lfrac
  real(wp), allocatable  :: eigenval(:,:)
  real(wp)               :: expected(n, ispin)        
  real(wp), allocatable  :: result_(:,:)
  
  ham_des%nel = [5, 4]
  lfrac = .true.
  allocate(eigenval(n, ispin))
  eigenval(:,1) = [-1.0_wp, -0.8_wp, -0.5_wp, -0.5_wp, -0.5_wp, -0.4_wp, -0.4_wp, -0.2_wp]
  eigenval(:,2) = [-1.0_wp, -0.8_wp, -0.5_wp, -0.5_wp, -0.5_wp, -0.4_wp, -0.4_wp, -0.2_wp]
  
  expected = 0.0_wp
  expected(1:5,1) = [1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp]
  expected(1:5,2) = [1.0_wp, 1.0_wp, frac, frac, frac]
  
  call set_occupations_eig(ham_des, result_, eigenval, lfrac)
  
  call assert_equals(expected, result_, n, ispin, tol)
end subroutine test_set_occupations_spin_p5


subroutine test_set_occupations_spin_p5_int
  use hamilton_vars, only: hamil_descriptor
  integer, parameter     :: ispin=2, n=8
  real(wp), parameter    :: tol=1.0e-10_wp
  real(wp), parameter    :: frac=0.666666666666_wp
  type(hamil_descriptor) :: ham_des
  logical                :: lfrac
  real(wp), allocatable  :: eigenval(:,:)
  real(wp)               :: expected(n, ispin)        
  real(wp), allocatable  :: result_(:,:)
  
  ham_des%nel = [5, 4]
  lfrac = .false.
  allocate(eigenval(n, ispin))
  eigenval(:,1) = [-1.0_wp, -0.8_wp, -0.5_wp, -0.5_wp, -0.5_wp, -0.4_wp, -0.4_wp, -0.2_wp]
  eigenval(:,2) = [-1.0_wp, -0.8_wp, -0.5_wp, -0.5_wp, -0.5_wp, -0.4_wp, -0.4_wp, -0.2_wp]
  
  expected = 0.0_wp
  expected(1:5,1) = [1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp]
  expected(1:4,2) = [1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp]
  
  call set_occupations_eig(ham_des, result_, eigenval, lfrac)
  
  call assert_equals(expected, result_, n, ispin, tol)
end subroutine test_set_occupations_spin_p5_int

subroutine test_rdm_trafo
  real(wp), allocatable :: rdm(:,:,:), rdm_(:,:,:)
  integer, parameter    :: n = 10
  real(wp), parameter   :: tol = 1.0e-10_wp
  
  allocate(rdm(n,n,4), rdm_(n,n,4))
  call random_number(rdm)
  rdm_ = rdm
  
  call spin_to_shell_rdm(rdm)
  call shell_to_spin_rdm(rdm)
  
  call assert_equals(rdm_(:,:,4), rdm(:,:,4), n, n, tol)
end subroutine test_rdm_trafo

subroutine test_calculate_rdm_ispin1
  integer, parameter     :: n=5, ispin=1
  real(wp), parameter    :: tol=1.0e-10_wp
  real(wp)               :: phi(n,n,ispin), occ(n,2), rdm_exp(n,n)
  real(wp), allocatable  :: rdm(:,:,:)
  
  occ = 0.0_wp
  occ(:,1) = [1.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
  occ(:,2) = [1.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
  
  phi = 0.0_wp
  phi(1,1,1) = 1.0_wp
  phi(2,2,1) = 1.0_wp
  
  rdm_exp = 0.0_wp
  rdm_exp(1,1) = 1.0_wp
  rdm_exp(2,2) = 1.0_wp
  
  call calculate_rdm(phi, occ, rdm)
  
  call assert_equals(rdm_exp, rdm(:,:,1), n, n, tol)
end subroutine test_calculate_rdm_ispin1

subroutine test_calculate_rdm_ispin2
  integer, parameter     :: n=5, ispin=2
  real(wp), parameter    :: tol=1.0e-10_wp, frac=0.3333333333_wp
  integer                :: i, spin
  real(wp)               :: phi(n,n,ispin), occ(n,2), rdm_exp(n,n)
  real(wp), allocatable  :: rdm(:,:,:)
  
  occ = 0.0_wp
  occ(:,1) = [1.0_wp, frac, frac, frac, 0.0_wp]
  occ(:,2) = [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
  
  phi = 0.0_wp
  do spin = 1, ispin
    do i = 1, n
      phi(i,i,spin) = 1.0_wp
    end do
  end do
  
  rdm_exp = 0.0_wp
  rdm_exp(1,1) = 2.0_wp
  rdm_exp(2,2) = frac
  rdm_exp(3,3) = frac
  rdm_exp(4,4) = frac
  
  call calculate_rdm(phi, occ, rdm)
  
  call assert_equals(rdm_exp(:,:), rdm(:,:,4), n, n, tol)
end subroutine test_calculate_rdm_ispin2

end module test_hfproc
