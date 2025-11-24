! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module test_slater_condon_rules

use lapack, only: gemm
use slater
use slater_condon_rules

use constants, only: wp
use fruit, only: assert_equals, run_test_case, test_set_summary

implicit none

private
public :: test_slater_condon_rules_set

contains

subroutine test_slater_condon_rules_set
  call run_test_case(test_hartree_consistency_ispin1, &
                     "Test Hartree vs g-resolved Hartree energy for relevant excitations with ispin=1")
  call run_test_case(test_hartree_consistency_ispin2, &
                     "Test Hartree vs g-resolved Hartree energy for relevant excitations with ispin=2")
  call run_test_case(test_exchange_consistency_ispin1, &
                     "Test exchange vs g-resolved exchange energy for relevant excitations with ispin=1")
  call run_test_case(test_exchange_consistency_ispin2, &
                     "Test exchange vs g-resolved exchange energy for relevant excitations with ispin=2")

  call test_set_summary("src/slater_condon_rules.f90")
end subroutine test_slater_condon_rules_set

subroutine test_hartree_consistency_ispin1
  use iso_c_binding, only: c_f_pointer, c_loc
  integer, parameter            :: n=10, ng=100, ispin1=1, ispin2=1, m=5
  integer                       :: g, spin, nn
  real(wp), parameter           :: tol = 1.0e-04_wp
  real(wp)                      :: eh(m), eh_(m), ex
  real(wp), allocatable         :: eh_gres(:), ex_gres(:)
  real(wp), allocatable, target :: h2(:,:,:,:,:), h2_gres(:,:,:,:)
  real(wp), pointer             :: h2_ptr(:,:,:)=>null(), h2_gres_ptr(:,:,:)=>null()
  type(slater_det)              :: det1, det2

  nn = n * n

  allocate(h2(n,n,n,n,ispin2))
  allocate(h2_gres(n,n,ng,ispin1))
  call c_f_pointer(c_loc(h2_gres), h2_gres_ptr, shape=[nn,ng,ispin1])
  call c_f_pointer(c_loc(h2), h2_ptr, shape=[nn, nn, ispin2])

  call random_number(h2_gres)
  !mimic g-resolved Hamiltonian 
  do spin = 1, ispin1
    do g = 1, ng
      h2_gres(:,:,g,spin) = 0.5_wp * (h2_gres(:,:,g,spin) + transpose(h2_gres(:,:,g,spin)))
    end do
  end do

  h2_ptr(:,:,1) = matmul(h2_gres_ptr(:,:,1), transpose(h2_gres_ptr(:,:,1)))

  !identical determinants
  det1 = slater_det(195_strkind, 60_strkind) !(1100001100, 0011110000)
  det2 = slater_det(195_strkind, 60_strkind) !(1100001100, 0011110000)
  call slater_condon_two(det1, det2, h2, eh(1), ex)
  call slater_condon_two_gres(det1, det2, h2_gres, eh_gres, ex_gres)
  eh_(1) = sum(eh_gres)

  !single excitation
  det1 = slater_det(340_strkind, 27_strkind)  !(0010101010, 1101100000)
  det2 = slater_det(340_strkind, 30_strkind)  !(0010101010, 0111100000)
  call slater_condon_two(det1, det2, h2, eh(2), ex)
  call slater_condon_two_gres(det1, det2, h2_gres, eh_gres, ex_gres)
  eh_(2) = sum(eh_gres)

  !double excitation
  det1 = slater_det(60_strkind, 15_strkind)  !(0011110000, 1111000000)
  det2 = slater_det(816_strkind, 15_strkind) !(0000110011, 1111000000)
  call slater_condon_two(det1, det2, h2, eh(3), ex)
  call slater_condon_two_gres(det1, det2, h2_gres, eh_gres, ex_gres)
  eh_(3) = sum(eh_gres)

  !double ab excitation
  det1 = slater_det(60_strkind, 340_strkind) !(0011110000, 0010101010)
  det2 = slater_det(92_strkind, 212_strkind) !(0011101000, 0010101100)
  call slater_condon_two(det1, det2, h2, eh(4), ex)
  call slater_condon_two_gres(det1, det2, h2_gres, eh_gres, ex_gres)
  eh_(4) = sum(eh_gres)

  !triple excitation
  det1 = slater_det(60_strkind, 340_strkind)  !(0011110000, 0010101010)
  det2 = slater_det(141_strkind, 212_strkind) !(1011000100, 0010101100)
  call slater_condon_two(det1, det2, h2, eh(5), ex)
  call slater_condon_two_gres(det1, det2, h2_gres, eh_gres, ex_gres)
  eh_(5) = sum(eh_gres)

  call assert_equals(eh, eh_, m, tol)
end subroutine test_hartree_consistency_ispin1


subroutine test_hartree_consistency_ispin2
  use iso_c_binding, only: c_f_pointer, c_loc
  integer, parameter            :: n=10, ng=10, ispin1=2, ispin2=3, m=5
  integer                       :: g, spin, nn
  real(wp), parameter           :: tol = 1.0e-04_wp
  real(wp)                      :: eh(m), eh_(m), ex
  real(wp), allocatable         :: eh_gres(:), ex_gres(:)
  real(wp), allocatable, target :: h2(:,:,:,:,:), h2_gres(:,:,:,:)
  real(wp), pointer             :: h2_ptr(:,:,:)=>null(), h2_gres_ptr(:,:,:)=>null()
  type(slater_det)              :: det1, det2

  nn = n * n

  allocate(h2(n,n,n,n,ispin2))
  allocate(h2_gres(n,n,ng,ispin1))
  call c_f_pointer(c_loc(h2_gres), h2_gres_ptr, shape=[nn,ng,ispin1])
  call c_f_pointer(c_loc(h2), h2_ptr, shape=[nn, nn, ispin2])

  call random_number(h2_gres)
  !mimic g-resolved Hamiltonian 
  do spin = 1, ispin1
    do g = 1, ng
      h2_gres(:,:,g,spin) = 0.5_wp * (h2_gres(:,:,g,spin) + transpose(h2_gres(:,:,g,spin)))
    end do
  end do

  h2_ptr(:,:,1) = matmul(h2_gres_ptr(:,:,1), transpose(h2_gres_ptr(:,:,1)))
  h2_ptr(:,:,2) = matmul(h2_gres_ptr(:,:,2), transpose(h2_gres_ptr(:,:,2)))
  h2_ptr(:,:,3) = matmul(h2_gres_ptr(:,:,1), transpose(h2_gres_ptr(:,:,2)))

  !identical determinants
  det1 = slater_det(771_strkind, 120_strkind) !(1100000011, 0001111000)
  det2 = slater_det(771_strkind, 120_strkind) !(1100000011, 0001111000)
  call slater_condon_two(det1, det2, h2, eh(1), ex)
  call slater_condon_two_gres(det1, det2, h2_gres, eh_gres, ex_gres)
  eh_(1) = sum(eh_gres)

  !single excitation
  det1 = slater_det(674_strkind, 108_strkind) !(0100010101, 0011011000)
  det2 = slater_det(674_strkind, 92_strkind)  !(0100010101, 0011101000)
  call slater_condon_two(det1, det2, h2, eh(2), ex)
  call slater_condon_two_gres(det1, det2, h2_gres, eh_gres, ex_gres)
  eh_(2) = sum(eh_gres)

  !double excitation
  det1 = slater_det(15_strkind, 277_strkind) !(1111000000, 1010100010)
  det2 = slater_det(15_strkind, 593_strkind) !(1111000000, 1000101001)
  call slater_condon_two(det1, det2, h2, eh(3), ex)
  call slater_condon_two_gres(det1, det2, h2_gres, eh_gres, ex_gres)
  eh_(3) = sum(eh_gres)

  !double ab excitation
  det1 = slater_det(408_strkind, 960_strkind) !(0001100110, 0000001111)
  det2 = slater_det(424_strkind, 706_strkind) !(0001010110, 0100001101)
  call slater_condon_two(det1, det2, h2, eh(4), ex)
  call slater_condon_two_gres(det1, det2, h2_gres, eh_gres, ex_gres)
  eh_(4) = sum(eh_gres)

  !triple excitation
  det1 = slater_det(58_strkind, 204_strkind)  !(0101110000, 0011001100)
  det2 = slater_det(170_strkind, 780_strkind) !(0101010100, 0011000011)
  call slater_condon_two(det1, det2, h2, eh(5), ex)
  call slater_condon_two_gres(det1, det2, h2_gres, eh_gres, ex_gres)
  eh_(5) = sum(eh_gres)

  call assert_equals(eh, eh_, m, tol)
end subroutine test_hartree_consistency_ispin2


subroutine test_exchange_consistency_ispin1
  integer, parameter            :: n=50, ng=50, ispin1=1, ispin2=1, m=4
  integer                       :: g, spin, nn
  real(wp), parameter           :: tol = 1.0e-04_wp
  real(wp)                      :: eh, ex(m), ex_(m)
  real(wp), allocatable         :: eh_gres(:), ex_gres(:)
  real(wp), allocatable, target :: h2(:,:,:,:,:), h2_gres(:,:,:,:)
  real(wp), pointer, contiguous :: h2_ptr(:,:,:)=>null(), h2_gres_ptr(:,:,:)=>null()
  type(slater_det)              :: det1, det2

  nn = n * n

  allocate(h2(n,n,n,n,ispin2))
  allocate(h2_gres(n,n,ng,ispin1))

  h2_ptr(1:nn,1:nn,1:ispin2) => h2
  h2_gres_ptr(1:nn,1:ng,1:ispin1) => h2_gres

  call random_number(h2_gres)
  !mimic g-resolved Hamiltonian 
  do spin = 1, ispin1
    do g = 1, ng
      h2_gres(:,:,g,spin) = 0.5_wp * (h2_gres(:,:,g,spin) + transpose(h2_gres(:,:,g,spin)))
    end do
  end do

  call gemm("n", "t", nn, nn, ng, 1.0_wp, h2_gres_ptr(:,:,1), nn, h2_gres_ptr(:,:,1), nn, 0.0_wp, h2_ptr(:,:,1), nn)

  !identical determinants
  det1 = slater_det(15_strkind, 240_strkind) !(11110000, 00001111)
  det2 = slater_det(15_strkind, 240_strkind) !(11110000, 00001111)
  call slater_condon_two(det1, det2, h2, eh, ex(1))
  call slater_condon_two_gres(det1, det2, h2_gres, eh_gres, ex_gres)
  ex_(1) = sum(ex_gres)

  !single excitation
  det1 = slater_det(15_strkind, 15_strkind)  !(111100000, 111100000)
  det2 = slater_det(15_strkind, 263_strkind) !(111100000, 111000001)
  call slater_condon_two(det1, det2, h2, eh, ex(2))
  call slater_condon_two_gres(det1, det2, h2_gres, eh_gres, ex_gres)
  ex_(2) = sum(ex_gres)

  !double excitation
  det1 = slater_det(60_strkind, 15_strkind)  !(0011110000, 1111000000)
  det2 = slater_det(816_strkind, 15_strkind) !(0000110011, 1111000000)
  call slater_condon_two(det1, det2, h2, eh, ex(3))
  call slater_condon_two_gres(det1, det2, h2_gres, eh_gres, ex_gres)
  ex_(3) = sum(ex_gres)

  !double ab excitation
  det1 = slater_det(60_strkind, 340_strkind) !(001111000, 001010101)
  det2 = slater_det(92_strkind, 212_strkind) !(001110100, 001010110)
  call slater_condon_two(det1, det2, h2, eh, ex(4))
  call slater_condon_two_gres(det1, det2, h2_gres, eh_gres, ex_gres)
  ex_(4) = sum(ex_gres)

  call assert_equals(ex, ex_, m, tol)
end subroutine test_exchange_consistency_ispin1


subroutine test_exchange_consistency_ispin2
  integer, parameter            :: n=50, ng=50, ispin1=2, ispin2=3, m=4
  integer                       :: g, spin, nn
  real(wp), parameter           :: tol = 1.0e-04_wp
  real(wp)                      :: eh, ex(m), ex_(m)
  real(wp), allocatable         :: eh_gres(:), ex_gres(:)
  real(wp), allocatable, target :: h2(:,:,:,:,:), h2_gres(:,:,:,:)
  real(wp), pointer, contiguous :: h2_ptr(:,:,:)=>null(), h2_gres_ptr(:,:,:)=>null()
  type(slater_det)              :: det1, det2

  nn = n * n

  allocate(h2(n,n,n,n,ispin2))
  allocate(h2_gres(n,n,ng,ispin1))

  h2_ptr(1:nn,1:nn,1:ispin2) => h2
  h2_gres_ptr(1:nn,1:ng,1:ispin1) => h2_gres

  call random_number(h2_gres)
  !mimic g-resolved Hamiltonian 
  do spin = 1, ispin1
    do g = 1, ng
      h2_gres(:,:,g,spin) = 0.5_wp * (h2_gres(:,:,g,spin) + transpose(h2_gres(:,:,g,spin)))
    end do
  end do

  call gemm("n", "t", nn, nn, ng, 1.0_wp, h2_gres_ptr(:,:,1), nn, h2_gres_ptr(:,:,1), nn, 0.0_wp, h2_ptr(:,:,1), nn)
  call gemm("n", "t", nn, nn, ng, 1.0_wp, h2_gres_ptr(:,:,2), nn, h2_gres_ptr(:,:,2), nn, 0.0_wp, h2_ptr(:,:,2), nn)
  call gemm("n", "t", nn, nn, ng, 1.0_wp, h2_gres_ptr(:,:,1), nn, h2_gres_ptr(:,:,2), nn, 0.0_wp, h2_ptr(:,:,3), nn)

  !identical determinants
  det1 = slater_det(15_strkind, 240_strkind) !(11110000, 00001111)
  det2 = slater_det(15_strkind, 240_strkind) !(11110000, 00001111)
  call slater_condon_two(det1, det2, h2, eh, ex(1))
  call slater_condon_two_gres(det1, det2, h2_gres, eh_gres, ex_gres)
  ex_(1) = sum(ex_gres)

  !single excitation
  det1 = slater_det(15_strkind, 15_strkind)  !(111100000, 111100000)
  det2 = slater_det(15_strkind, 263_strkind) !(111100000, 111000001)
  call slater_condon_two(det1, det2, h2, eh, ex(2))
  call slater_condon_two_gres(det1, det2, h2_gres, eh_gres, ex_gres)
  ex_(2) = sum(ex_gres)

  !double excitation
  det1 = slater_det(60_strkind, 15_strkind)  !(0011110000, 1111000000)
  det2 = slater_det(816_strkind, 15_strkind) !(0000110011, 1111000000)
  call slater_condon_two(det1, det2, h2, eh, ex(3))
  call slater_condon_two_gres(det1, det2, h2_gres, eh_gres, ex_gres)
  ex_(3) = sum(ex_gres)

  !double ab excitation
  det1 = slater_det(60_strkind, 340_strkind) !(001111000, 001010101)
  det2 = slater_det(92_strkind, 212_strkind) !(001110100, 001010110)
  call slater_condon_two(det1, det2, h2, eh, ex(4))
  call slater_condon_two_gres(det1, det2, h2_gres, eh_gres, ex_gres)
  ex_(4) = sum(ex_gres)

  call assert_equals(ex, ex_, m, tol)
end subroutine test_exchange_consistency_ispin2

end module test_slater_condon_rules