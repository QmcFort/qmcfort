! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module cas_hamilton

!********************************************************************************
!
! CasHamilton type is meant to be a comapact representation of the CAS
! Hamilonian for the CAS wave function in AFQMC
!
! Some requirements:
!    basis = mean_field (canonical basis assumed similar to FS approximation)
!    assumes nfrozen < 64 since Slater-Condon rules are used
!
!********************************************************************************

#include "../preproc.inc"
use constants
use slater
use slater_condon_rules
use hamilton_type
use hamilton_operations

implicit none

type, extends(Hamilton) :: CasHamilton
  real(wp)              :: e1_frozen
  real(wp), allocatable :: lgmean_frozen(:)
  real(wp), allocatable :: eh_frozen(:)
  real(wp), allocatable :: ex_frozen(:)
  real(wp), allocatable :: eself_frozen(:)
contains
  procedure, private :: project_cas_hamilton
end type CasHamilton

interface CasHamilton
  module procedure init_cas_hamilton
end interface CasHamilton

contains

!********************************************************************************
!
! CasHamilton constructor
! 
!********************************************************************************
function init_cas_hamilton(ham, nfrozen, nactive) result(self)
  type(Hamilton), intent(in) :: ham
  integer, intent(in)        :: nfrozen
  integer, intent(in)        :: nactive
  type(CasHamilton)          :: self
  !local variables
  integer                    :: i, n, ispin1
  integer(strkind)           :: slater_string
  real(wp), allocatable      :: hself_gres_frozen(:,:,:,:)
  real(wp), pointer          :: h2_gres_ob(:,:,:,:)
  type(slater_det)           :: frozen_det

  n = ham%des%n
  ispin1 = ham%des%ispin1

  !debug: note that the current method relies on determinants that can store up to 64 frozen orbitals
  slater_string = 0_strkind
  do i = 1, nfrozen
    slater_string = ibset(slater_string, i-1)
  end do
  frozen_det = slater_det(slater_string, slater_string)

  call ham%get_h2_gres_ob(h2_gres_ob)

  call slater_condon_one(frozen_det, frozen_det, ham%h1, self%e1_frozen)
  call slater_condon_one_array(frozen_det, frozen_det, h2_gres_ob, self%lgmean_frozen)
  call slater_condon_two_gres(frozen_det, frozen_det, h2_gres_ob, self%eh_frozen, self%ex_frozen)

  !calculate self energy Hamiltonian in frozen space only: max(1,nfrozen) needed to avoid strange behavior for nfrozen=0
  call calculate_gres_exchange_potential(h2_gres_ob, [1,1], [n, n], ispin1, hself_gres_frozen, 1, max(1,nfrozen))
  call slater_condon_one_array(frozen_det, frozen_det, hself_gres_frozen, self%eself_frozen)

  call self%project_cas_hamilton(ham, nfrozen, nactive)
end function init_cas_hamilton



!********************************************************************************
!
! Frozen-core like projection of the original Hamiltonian
!
!********************************************************************************
subroutine project_cas_hamilton(self, ham, nfrozen, nactive)
  class(CasHamilton), intent(inout) :: self
  type(Hamilton), intent(in)        :: ham
  integer, intent(in)               :: nfrozen
  integer, intent(in)               :: nactive
  !local variables
  integer                           :: n, ispin1
  integer                           :: first, last
  real(wp), pointer                 :: h2_gres_ob(:,:,:,:), hcas_h2_gres_ob(:,:,:,:)

  n = ham%des%n
  ispin1 = ham%des%ispin1

  first = nfrozen + 1
  last = nfrozen + nactive

  call ham%get_h2_gres_ob(h2_gres_ob)

  call calculate_gres_hartree_potential(h2_gres_ob, [1,1], [nfrozen,nfrozen], ispin1, ham%des%rspin, self%hpot_gres, first, last)
  call calculate_gres_exchange_potential(h2_gres_ob, [1,1], [nfrozen,nfrozen], ispin1, self%xpot_gres, first, last)
  call calculate_gres_exchange_potential(h2_gres_ob, [1,1], [n,n], ispin1, self%hself_gres, first, last)
  self%xpot_gres = -self%xpot_gres
  self%hself_gres = -self%hself_gres

  call sum_gres_potential(self%hpot_gres, self%hpot)
  call sum_gres_potential(self%xpot_gres, self%xpot)
  call sum_gres_potential(self%hself_gres, self%hself)

  self%des = ham%des
  self%des%n = nactive
  self%des%nbtot = nactive

  self%h_io = ham%h_io
  self%enuc = ham%enuc
  self%h0 = ham%h0

  if (allocated(ham%overlap)) then
    allocate(self%overlap(nactive,nactive,self%des%ispin1))
    self%overlap = ham%overlap(first:last,first:last,:)
  end if

  if (allocated(ham%ext_pot)) then
    allocate(self%ext_pot(nactive,nactive,self%des%ispin1))
    self%ext_pot = ham%ext_pot(first:last,first:last,:)
  end if

  if (allocated(ham%kinetic)) then
    allocate(self%kinetic(nactive,nactive,self%des%ispin1))
    self%kinetic = ham%kinetic(first:last,first:last,:)
  end if

  if (allocated(ham%hself)) then
    allocate(self%hself(nactive,nactive,self%des%ispin1))
    self%hself = ham%hself(first:last,first:last,:)
  end if

  if (allocated(ham%h1)) then
    allocate(self%h1(nactive,nactive,self%des%ispin1))
    self%h1 = ham%h1(first:last,first:last,:)
  end if

  if (allocated(ham%h1_orig)) then
    allocate(self%h1_orig(nactive,nactive,self%des%ispin1))
    self%h1_orig = ham%h1_orig(first:last,first:last,:)
  end if

  if (allocated(ham%h2)) then
    allocate(self%h2(nactive,nactive,nactive,nactive,self%des%ispin2))
    call self%associate_h2_ptr()
    self%h2 = ham%h2(first:last,first:last,first:last,first:last,:)
  end if

  if (allocated(ham%h2_gres)) then
    allocate(self%h2_gres(nactive*nactive, self%des%ng, self%des%ispin1))
    call self%get_h2_gres_ob(hcas_h2_gres_ob)
    call ham%get_h2_gres_ob(h2_gres_ob)
    hcas_h2_gres_ob = h2_gres_ob(first:last,first:last,:,:)
  end if

  if (allocated(ham%eigenval)) then
    allocate(self%eigenval(nactive,self%des%ispin))
    self%eigenval = ham%eigenval(first:last,:)
  end if

  if (allocated(ham%phi)) then
    allocate(self%phi(nactive,nactive,self%des%ispin))
    self%phi = ham%phi(first:last,first:last,:)
  end if

  if (allocated(ham%occ)) then
    allocate(self%occ(nactive,2))
    self%occ = ham%occ(first:last,:)
  end if

  if (allocated(ham%p1)) then
    allocate(self%p1(nactive,nactive,self%des%ispin1))
    self%p1 = ham%p1(first:last,first:last,:)
  end if
end subroutine project_cas_hamilton

end module cas_hamilton