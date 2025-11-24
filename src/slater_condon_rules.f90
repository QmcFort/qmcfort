! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module slater_condon_rules
!******************************************************************************** 
!
! Routines to calculate various matrix element between orthogonal Slater determinants
!
!******************************************************************************** 

#include "preproc.inc"
use constants
use slater

use mpi, only: comm_world

implicit none

private
public :: slater_condon, slater_condon_one, slater_condon_one_array, slater_condon_two, &
          slater_condon_two_gres

contains


!******************************************************************************** 
!
! Calculate Hamiltonian matrix elements between orthogonal Slater determinants
!
! Calculate difference between determinants and call specific routine
!
! Matrix element is zero if triple or higher excitation
!
!******************************************************************************** 
subroutine slater_condon(det1, det2, h1, h2, e1, eh, ex)
  type(slater_det), intent(in) :: det1
  type(slater_det), intent(in) :: det2
  real(wp), intent(in)         :: h1(:,:,:)
  real(wp), intent(in)         :: h2(:,:,:,:,:)
  real(wp), intent(out)        :: e1
  real(wp), intent(out)        :: eh
  real(wp), intent(out)        :: ex
  !local variables        
  integer                      :: ndiff

  e1 = 0.0_wp
  eh = 0.0_wp
  ex = 0.0_wp

  ndiff = det1%exc_deg(det2)

  if (ndiff == 0) then
    call slater_condon_same(det1, h1, h2, e1, eh, ex)
  else if (ndiff == 1) then
    call slater_condon_single(det1, det2, h1, h2, e1, eh, ex)
  else if (ndiff == 2) then
    call slater_condon_double(det1, det2, h2, eh, ex)
  end if
end subroutine slater_condon


!******************************************************************************** 
!
! Matrix element between two identical Slater determinants
!
!    <Phi|H|Phi> = \sum_{i} h_{ii} + 0.5 \sum_{ij} <ij||ij> 
!
!******************************************************************************** 
subroutine slater_condon_same(det, h1, h2, e1, eh, ex)
  type(slater_det), intent(in) :: det
  real(wp), intent(in)         :: h1(:,:,:)
  real(wp), intent(in)         :: h2(:,:,:,:,:)
  real(wp), intent(inout)      :: e1
  real(wp), intent(inout)      :: eh
  real(wp), intent(inout)      :: ex
  !local variables
  integer                      :: n, i, j, nia, nib, nja, njb
  integer                      :: spin_a, spin_b, spin_ab
  real(wp), parameter          :: f = 0.5_wp

  n = det%nmax()

  if (size(h1,3) == 1) then
    spin_a = 1
    spin_b = 1
    spin_ab = 1    
  else 
    spin_a = 1
    spin_b = 2
    spin_ab = 3
  end if

  do j = 1, n
    nja = det%occ(j, 1)
    njb = det%occ(j, 2)

    e1 = e1 + nja*h1(j,j,spin_a) + njb*h1(j,j,spin_b)
    eh = eh + f*nja*h2(j,j,j,j,spin_a) + f*njb*h2(j,j,j,j,spin_b) + 2.0_wp*f*nja*njb*h2(j,j,j,j,spin_ab)
    ex = ex - f*nja*h2(j,j,j,j,spin_a) - f*njb*h2(j,j,j,j,spin_b)

    do i = j+1, n 
      nia = det%occ(i, 1)
      nib = det%occ(i, 2)

      eh = eh + nia*nja*h2(i,i,j,j,spin_a) + nib*njb*h2(i,i,j,j,spin_b) + &
                nia*njb*h2(i,i,j,j,spin_ab) + nib*nja*h2(j,j,i,i,spin_ab)
      ex = ex - nia*nja*h2(i,j,j,i,spin_a) - nib*njb*h2(i,j,j,i,spin_b)
    end do
  end do
end subroutine slater_condon_same


!******************************************************************************** 
!
! Matrix element between Slater determinant and its single excitation
!
!    <Phi|H|Phi_i^a> = h_{ia} + \sum_{j} <ij||aj> 
!
!******************************************************************************** 
subroutine slater_condon_single(det1, det2, h1, h2, e1, eh, ex)
  type(slater_det), intent(in) :: det1
  type(slater_det), intent(in) :: det2
  real(wp), intent(in)         :: h1(:,:,:)
  real(wp), intent(in)         :: h2(:,:,:,:,:)
  real(wp), intent(inout)      :: e1
  real(wp), intent(inout)      :: eh
  real(wp), intent(inout)      :: ex
  !local variables
  integer                      :: n, i, a, p, spin, spin_a, spin_ab
  integer                      :: na, nb, sgn
  integer, allocatable         :: diff(:,:)
  
  n = det1%nmax()

  diff = det1%diff(det2)

  i = diff(1,1)
  a = diff(2,1)
  spin = diff(3,1)

  if (size(h1,3) == 1) then
    spin_a = 1
    spin_ab = 1
  else
    spin_a = spin
    spin_ab = 3
  end if

  sgn = det1%phase(det2)

  e1 = 1 + sgn * h1(i,a,spin_a)

  do p = 1, n
    na = det1%occ(p, spin)
    nb = det1%occ(p, 3-spin)  

    eh = eh + sgn * na * h2(i,a,p,p,spin_a)
    if (spin == 1) then
      eh = eh + sgn * nb * h2(i,a,p,p,spin_ab)
    else 
      eh = eh + sgn * nb * h2(p,p,i,a,spin_ab)
    end if

    ex = ex - sgn * na*h2(i,p,p,a,spin_a)
  end do
end subroutine slater_condon_single


!******************************************************************************** 
!
! Matrix element between Slater determinant and its double excitation
!
!    <Phi|H|Phi_{ij}^{ab}> = <ij||ab> 
!
!******************************************************************************** 
subroutine slater_condon_double(det1, det2, h2, eh, ex)
  type(slater_det), intent(in) :: det1
  type(slater_det), intent(in) :: det2
  real(wp), intent(in)         :: h2(:,:,:,:,:)
  real(wp), intent(inout)      :: eh
  real(wp), intent(inout)      :: ex
  !local variables
  integer                      :: i, j, a, b, spin, spin1, spin2, sgn
  integer, allocatable         :: diff(:,:)

  diff = det1%diff(det2)

  i = diff(1,1)
  a = diff(2,1)
  spin1 = diff(3,1)

  j = diff(1,2)
  b = diff(2,2)
  spin2 = diff(3,2)

  if (size(h2,5) == 1) then
    spin = 1
  else
    if (spin1 == spin2) then
      spin = spin1
    else
      spin = 3
    end if
  end if

  sgn = det1%phase(det2)

  if (spin2==1 .and. spin1==2) then
    eh = eh + sgn*h2(j,b,i,a,spin)
  else 
    eh = eh + sgn*h2(i,a,j,b,spin)
  end if

  if (spin1 == spin2) ex = ex - sgn*h2(i,b,j,a,spin)
end subroutine slater_condon_double


!******************************************************************************** 
!
! Calculate one body matrix element between orthogonal Slater determinants
!
! Calculate difference between determinants and call specific routine
!
! Matrix element is zero if double or higher excitation
!
!******************************************************************************** 
subroutine slater_condon_one(det1, det2, h1, e)
  type(slater_det), intent(in) :: det1
  type(slater_det), intent(in) :: det2
  real(wp), intent(in)         :: h1(:,:,:)
  real(wp), intent(out)        :: e
  !local variables        
  integer                      :: ndiff

  e = 0.0_wp

  ndiff = det1%exc_deg(det2)

  if (ndiff == 0) then
    call slater_condon_one_same(det1, h1, e)
  else if (ndiff == 1) then
    call slater_condon_one_single(det1, det2, h1, e)
  end if
end subroutine slater_condon_one


!******************************************************************************** 
!
! One-body matrix element between two identical Slater determinants
!
!    <Phi|H|Phi> = \sum_{sig,i} h_{ii}^{sig}
!
!******************************************************************************** 
subroutine slater_condon_one_same(det, h1, e)
  type(slater_det), intent(in) :: det
  real(wp), intent(in)         :: h1(:,:,:)
  real(wp), intent(inout)      :: e
  !local variables
  integer                      :: n, i, nia, nib
  integer                      :: spin_a, spin_b

  n = det%nmax()

  if (size(h1,3) == 1) then
    spin_a = 1
    spin_b = 1
  else 
    spin_a = 1
    spin_b = 2
  end if

  do i = 1, n
    nia = det%occ(i, 1)
    nib = det%occ(i, 2)
    e = e + nia*h1(i,i,spin_a) + nib*h1(i,i,spin_b)
  end do
end subroutine slater_condon_one_same


!******************************************************************************** 
!
! One-body matrix element between Slater determinant and its single excitation
!
!    <Phi|H|Phi_i^a> = h_{ia}^{sig}
!
!******************************************************************************** 
subroutine slater_condon_one_single(det1, det2, h1, e)
  type(slater_det), intent(in) :: det1
  type(slater_det), intent(in) :: det2
  real(wp), intent(in)         :: h1(:,:,:)
  real(wp), intent(inout)      :: e
  !local variables
  integer                      :: i, a, spin
  integer, allocatable         :: diff(:,:)
  
  diff = det1%diff(det2)
  i = diff(1,1)
  a = diff(2,1)
  spin = diff(3,1)

  if (size(h1,3) == 1) spin = 1

  e = e + h1(i,a,spin) * det1%phase(det2)
end subroutine slater_condon_one_single


!******************************************************************************** 
!
! Calculate one-body matrix elements between orthogonal Slater determinants
! for an array of one-body operators
!
! First, calculate difference between determinants and call specific routine
!
! Matrix element is zero if double or higher excitation
!
!******************************************************************************** 
subroutine slater_condon_one_array(det1, det2, h1_arr, e)
  type(slater_det), intent(in)       :: det1
  type(slater_det), intent(in)       :: det2
  real(wp), intent(in)               :: h1_arr(:,:,:,:)
  real(wp), allocatable, intent(out) :: e(:)
  !local variables        
  integer                            :: ng, ndiff

  ng = size(h1_arr, 3)
  
  allocate(e(ng))
  e = 0.0_wp

  ndiff = det1%exc_deg(det2)

  if (ndiff == 0) then
    call slater_condon_one_array_same(det1, h1_arr, e)
  else if (ndiff == 1) then
    call slater_condon_one_array_single(det1, det2, h1_arr, e)
  end if
end subroutine slater_condon_one_array


!******************************************************************************** 
!
! One-body matrix element between two identical Slater determinants
! for an array of one-body operators
!
!    <Phi|H_g|Phi> = \sum_{sig,i} h_{g,ii}^{sig}
!
!******************************************************************************** 
subroutine slater_condon_one_array_same(det, h1_arr, e)
  type(slater_det), intent(in) :: det
  real(wp), intent(in)         :: h1_arr(:,:,:,:)
  real(wp), intent(inout)      :: e(:)
  !local variables
  integer                      :: n, i, nia, nib
  integer                      :: spin_a, spin_b

  n = det%nmax()

  if (size(h1_arr,4) == 1) then
    spin_a = 1
    spin_b = 1
  else 
    spin_a = 1
    spin_b = 2
  end if

  do i = 1, n
    nia = det%occ(i, 1)
    nib = det%occ(i, 2)
    e = e + nia*h1_arr(i,i,:,spin_a) + nib*h1_arr(i,i,:,spin_b)
  end do
end subroutine slater_condon_one_array_same


!******************************************************************************** 
!
! One-body matrix element between Slater determinant and its single excitation
! for an array of one-body operators
!
!    <Phi|H_g|Phi_i^a> = h_{g,ia}^{sig}
!
!******************************************************************************** 
subroutine slater_condon_one_array_single(det1, det2, h1_arr, e)
  type(slater_det), intent(in) :: det1
  type(slater_det), intent(in) :: det2
  real(wp), intent(in)         :: h1_arr(:,:,:,:)
  real(wp), intent(inout)      :: e(:)
  !local variables
  integer                      :: i, a, spin
  integer, allocatable         :: diff(:,:)
  
  diff = det1%diff(det2)
  i = diff(1,1)
  a = diff(2,1)
  spin = diff(3,1)

  if (size(h1_arr,4) == 1) spin = 1

  e = e + h1_arr(i,a,:,spin) * det1%phase(det2)
end subroutine slater_condon_one_array_single


!******************************************************************************** 
!
! Calculate Hartree and exchange energy between orthogonal Slater determinants
!
! Calculate difference between determinants and call specific routine
!
! Matrix element is zero if triple or higher excitation
!
!******************************************************************************** 
subroutine slater_condon_two(det1, det2, h2, eh, ex)
  type(slater_det), intent(in) :: det1
  type(slater_det), intent(in) :: det2
  real(wp), intent(in)         :: h2(:,:,:,:,:)
  real(wp), intent(out)        :: eh
  real(wp), intent(out)        :: ex
  !local variables        
  integer                      :: ndiff

  eh = 0.0_wp
  ex = 0.0_wp

  ndiff = det1%exc_deg(det2)

  if (ndiff == 0) then
    call slater_condon_two_same(det1, h2, eh, ex)
  else if (ndiff == 1) then
    call slater_condon_two_single(det1, det2, h2, eh, ex)
  else if (ndiff == 2) then
    call slater_condon_double(det1, det2, h2, eh, ex)
  end if
end subroutine slater_condon_two


!******************************************************************************** 
!
! Calculate Hartree and exchange energy between two identical Slater determinants
!
!    <Phi|H_2|Phi> = 0.5 \sum_{ij} <ij||ij> 
!
!******************************************************************************** 
subroutine slater_condon_two_same(det, h2, eh, ex)
  type(slater_det), intent(in) :: det
  real(wp), intent(in)         :: h2(:,:,:,:,:)
  real(wp), intent(inout)      :: eh
  real(wp), intent(inout)      :: ex
  !local variables
  integer                      :: n, i, j, nia, nib, nja, njb
  integer                      :: spin_a, spin_b, spin_ab
  real(wp), parameter          :: f = 0.5_wp

  n = det%nmax()

  if (size(h2,5) == 1) then
    spin_a = 1
    spin_b = 1
    spin_ab = 1    
  else 
    spin_a = 1
    spin_b = 2
    spin_ab = 3
  end if

  do j = 1, n
    nja = det%occ(j, 1)
    njb = det%occ(j, 2)

    eh = eh + f*nja*h2(j,j,j,j,spin_a) + f*njb*h2(j,j,j,j,spin_b) + 2.0_wp*f*nja*njb*h2(j,j,j,j,spin_ab)
    ex = ex - f*nja*h2(j,j,j,j,spin_a) - f*njb*h2(j,j,j,j,spin_b)

    do i = j+1, n 
      nia = det%occ(i, 1)
      nib = det%occ(i, 2)

      eh = eh + nia*nja*h2(i,i,j,j,spin_a) + nib*njb*h2(i,i,j,j,spin_b) + &
                nia*njb*h2(i,i,j,j,spin_ab) + nib*nja*h2(j,j,i,i,spin_ab)
      ex = ex - nia*nja*h2(i,j,j,i,spin_a) - nib*njb*h2(i,j,j,i,spin_b)
    end do
  end do
end subroutine slater_condon_two_same


!******************************************************************************** 
!
! Calculate Hartree and exchange energy between Slater determinant and 
! its single excitation
!
!    <Phi|H_2|Phi_i^a> = \sum_{p} <ip||ap> 
!
!******************************************************************************** 
subroutine slater_condon_two_single(det1, det2, h2, eh, ex)
  type(slater_det), intent(in) :: det1
  type(slater_det), intent(in) :: det2
  real(wp), intent(in)         :: h2(:,:,:,:,:)
  real(wp), intent(inout)      :: eh
  real(wp), intent(inout)      :: ex
  !local variables
  integer                      :: n, i, a, p, spin, spin_a, spin_b, spin_ab
  integer                      :: na, nb, sgn
  integer, allocatable         :: diff(:,:)
  
  n = det1%nmax()

  diff = det1%diff(det2)

  i = diff(1,1)
  a = diff(2,1)
  spin = diff(3,1)

  if (size(h2,5) == 1) then
    spin_a = 1
    spin_b = 1
    spin_ab = 1
  else
    spin_a = spin
    spin_b = 3 - spin
    spin_ab = 3
  end if

  sgn = det1%phase(det2)

  do p = 1, n
    na = det1%occ(p, spin)
    nb = det1%occ(p, 3-spin)  

    eh = eh + sgn * na * h2(i,a,p,p,spin_a)
    ex = ex - sgn * na * h2(i,p,p,a,spin_a)

    if (spin == 1) then
      eh = eh + sgn * nb * h2(i,a,p,p,spin_ab)
    else 
      eh = eh + sgn * nb * h2(p,p,i,a,spin_ab)
    end if
  end do
end subroutine slater_condon_two_single


!******************************************************************************** 
!
! Calculate g-resolved exchange energy between orthogonal Slater determinants
!
! Calculate difference between determinants and call specific routine
!
! Matrix element is zero if triple or higher excitation
!
!******************************************************************************** 
subroutine slater_condon_two_gres(det1, det2, h2_gres, eh, ex)
  type(slater_det), intent(in)       :: det1
  type(slater_det), intent(in)       :: det2
  real(wp), intent(in)               :: h2_gres(:,:,:,:)
  real(wp), allocatable, intent(out) :: eh(:)
  real(wp), allocatable, intent(out) :: ex(:)
  !local variables        
  integer                            :: ng, ndiff
  
  ng = size(h2_gres, 3)

  allocate(eh(ng), ex(ng))
  eh = 0.0_wp
  ex = 0.0_wp

  ndiff = det1%exc_deg(det2)

  if (ndiff == 0) then
    call slater_condon_two_gres_same(det1, h2_gres, eh, ex)
  else if (ndiff == 1) then
    call slater_condon_two_gres_single(det1, det2, h2_gres, eh, ex)
  else if (ndiff == 2) then
    call slater_condon_two_gres_double(det1, det2, h2_gres, eh, ex)
  end if
end subroutine slater_condon_two_gres


!******************************************************************************** 
!
! Calculate g-resolved Hartree and exchange energies 
! between two identical Slater determinants
!
!    <Phi|H_H,g|Phi> =  0.5 \sum_{sig,tau,ij} L_{g,ii}^{sig} L_{g,jj}^{tau}
!    <Phi|H_x,g|Phi> = -0.5 \sum_{sig,ij} L_{g,ij}^{sig} L_{g,ji}^{sig}
!
!******************************************************************************** 
subroutine slater_condon_two_gres_same(det, h2_gres, eh, ex)
  type(slater_det), intent(in) :: det
  real(wp), intent(in)         :: h2_gres(:,:,:,:)
  real(wp), intent(inout)      :: eh(:)
  real(wp), intent(inout)      :: ex(:)
  !local variables
  integer                      :: n, i, j, nia, nib, nja, njb
  integer                      :: spin_a, spin_b
  real(wp), parameter          :: f = 0.5_wp

  n = det%nmax()

  if (size(h2_gres,4) == 1) then
    spin_a = 1
    spin_b = 1
  else 
    spin_a = 1
    spin_b = 2
  end if

  do j = 1, n
    nja = det%occ(j, 1)
    njb = det%occ(j, 2)

    eh = eh + f*nja*h2_gres(j,j,:,spin_a)*h2_gres(j,j,:,spin_a) &    
            + f*njb*h2_gres(j,j,:,spin_b)*h2_gres(j,j,:,spin_b) &
            + 2.0_wp*f*nja*njb*h2_gres(j,j,:,spin_a)*h2_gres(j,j,:,spin_b)

    ex = ex - f*nja*h2_gres(j,j,:,spin_a)*h2_gres(j,j,:,spin_a) &
            - f*njb*h2_gres(j,j,:,spin_b)*h2_gres(j,j,:,spin_b)

    do i = j+1, n 
      nia = det%occ(i, 1)
      nib = det%occ(i, 2)

      eh = eh + nia*nja*h2_gres(i,i,:,spin_a)*h2_gres(j,j,:,spin_a) &
              + nib*njb*h2_gres(i,i,:,spin_b)*h2_gres(j,j,:,spin_b) &
              + nia*njb*h2_gres(i,i,:,spin_a)*h2_gres(j,j,:,spin_b) &
              + nib*nja*h2_gres(j,j,:,spin_a)*h2_gres(i,i,:,spin_b)

      ex = ex - nia*nja*h2_gres(i,j,:,spin_a)*h2_gres(j,i,:,spin_a) &
              - nib*njb*h2_gres(i,j,:,spin_b)*h2_gres(j,i,:,spin_b)
    end do
  end do
end subroutine slater_condon_two_gres_same


!******************************************************************************** 
!
! Calculate g-resolved Hartree and exchange energies 
! between Slater determinant and its single excitation
!
!    <Phi|H_H,g|Phi_i^a> =   \sum_{p,tau} L_{g,ia}^{sig} L_{g,pp}^{tau} 
!    <Phi|H_x,g|Phi_i^a> = - \sum_{p} L_{g,ip}^{sig} L_{g,pa}^{sig} 
!
!******************************************************************************** 
subroutine slater_condon_two_gres_single(det1, det2, h2_gres, eh, ex)
  type(slater_det), intent(in) :: det1
  type(slater_det), intent(in) :: det2
  real(wp), intent(in)         :: h2_gres(:,:,:,:)
  real(wp), intent(inout)      :: eh(:)
  real(wp), intent(inout)      :: ex(:)
  !local variables
  integer                      :: n, i, a, p, spin, spin_a, spin_b
  integer                      :: na, nb, sgn
  integer, allocatable         :: diff(:,:)
  
  n = det1%nmax()

  diff = det1%diff(det2)

  i = diff(1,1)
  a = diff(2,1)
  spin = diff(3,1)

  if (size(h2_gres,4) == 1) then
    spin_a = 1
    spin_b = 1
  else
    spin_a = spin
    spin_b = 3 - spin
  end if

  sgn = det1%phase(det2)

  do p = 1, n
    na = det1%occ(p, spin)
    nb = det1%occ(p, 3-spin)  

    eh = eh + sgn * na * h2_gres(i,a,:,spin_a) * h2_gres(p,p,:,spin_a) &
            + sgn * nb * h2_gres(i,a,:,spin_a) * h2_gres(p,p,:,spin_b)

    ex = ex - sgn * na * h2_gres(i,p,:,spin_a) * h2_gres(p,a,:,spin_a)
  end do
end subroutine slater_condon_two_gres_single


!******************************************************************************** 
!
! Calculate g-resolved exchange energy between Slater determinant 
! and its double excitation
!
!    <Phi|H_H,g|Phi_{ij}^{ab}> =   L_{g,ia}^{sig} L_{g,jb}^{tau}
!    <Phi|H_x,g|Phi_{ij}^{ab}> = - L_{g,ib}^{sig} L_{g,ja}^{sig} \delta_{sig,tau}
!
!******************************************************************************** 
subroutine slater_condon_two_gres_double(det1, det2, h2_gres, eh, ex)
  type(slater_det), intent(in) :: det1
  type(slater_det), intent(in) :: det2
  real(wp), intent(in)         :: h2_gres(:,:,:,:)
  real(wp), intent(inout)      :: eh(:)
  real(wp), intent(inout)      :: ex(:)
  !local variables
  integer                      :: i, j, a, b, spin1, spin2, spin_a, spin_b, sgn
  integer, allocatable         :: diff(:,:)

  diff = det1%diff(det2)

  i = diff(1,1)
  a = diff(2,1)
  spin1 = diff(3,1)

  j = diff(1,2)
  b = diff(2,2)
  spin2 = diff(3,2)

  if (size(h2_gres,4) == 1) then
    spin_a = 1
    spin_b = 1
  else 
    spin_a = spin1
    spin_b = spin2
  end if

  sgn = det1%phase(det2)

  eh = eh + sgn * h2_gres(i,a,:,spin_a) * h2_gres(j,b,:,spin_b)
  if (spin1 == spin2) ex = ex - sgn * h2_gres(i,b,:,spin_a) * h2_gres(j,a,:,spin_b)
end subroutine slater_condon_two_gres_double

end module slater_condon_rules