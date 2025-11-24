! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module hamilton_vars

use constants
use spin_utils
use qmcfort_pos
use mpi, only: comm_world

implicit none

private
public hamil_descriptor

type hamil_descriptor
  integer                :: n
  integer                :: nbtot
  integer                :: nfrozen
  integer                :: ispin
  integer                :: ispin1
  integer                :: ispin2
  integer                :: ispin_fock
  integer                :: rspin
  !
  integer                :: nelect
  integer                :: nocc
  integer                :: nocc_
  integer                :: nel(2)
  integer                :: nell(2)
  integer                :: nel_space(2)
  integer                :: nell_space(2)
  integer                :: spin
  !
  real(wp)               :: h0
  !
  integer                :: ng
  integer                :: ng_
  logical                :: compress  = .false.
  real(wp)               :: chol_tol  = 1.0e-06_wp
  real(wp)               :: chol_res
  real(wp)               :: tol_spars = 1.0e-08_wp
  real(wp)               :: sparsity
  character(len=charlen) :: compress_mode  = "cholesky"
  !
  character(len=charlen) :: integral_mode
  character(len=charlen) :: basis
  character(len=charlen) :: basis_set
contains
  procedure              :: init          => init_ham_des
  generic, public        :: assignment(=) => assign_hamil_descriptor
  procedure, public      :: set_electrons_struc, set_electrons_nel, set_electrons_spin
  procedure, private     :: assign_hamil_descriptor
  procedure, public      :: get_spin, get_spin1, get_spin2
end type hamil_descriptor

contains

!******************************************************************************** 
!
!       Initialize default values in hamiltonian descriptor:
!
!******************************************************************************** 
subroutine init_ham_des(ham_des)
  class(hamil_descriptor), intent(inout) :: ham_des
  !
  ham_des%nbtot = unset
  ham_des%nfrozen = 0
  ham_des%integral_mode = "eri"
  ham_des%ispin = unset
  ham_des%ispin1 = 1
  ham_des%ispin2 = 1
  ham_des%ispin_fock = 1
  ham_des%nel = unset
  ham_des%spin = unset
  ham_des%nelect = unset

  ham_des%basis_set = ""
end subroutine init_ham_des


!******************************************************************************** 
!
!       Assign Hamiltonian descriptor
!
!******************************************************************************** 
subroutine assign_hamil_descriptor(to, from)
  class(hamil_descriptor), intent(inout) :: to
  class(hamil_descriptor), intent(in)    :: from
  !
  to%n = from%n
  to%nbtot = from%nbtot
  to%nfrozen = from%nfrozen
  to%nelect = from%nelect
  to%nocc = from%nocc
  to%nocc_ = from%nocc_
  to%nel = from%nel
  to%nell = from%nell
  to%nel_space = from%nel_space
  to%spin = from%spin
  to%rspin = from%rspin
  to%ispin = from%ispin
  to%ispin1 = from%ispin1 
  to%ispin2 = from%ispin2
  to%ispin_fock = from%ispin_fock
  to%ng = from%ng
  to%compress = from%compress
  to%integral_mode = from%integral_mode
  to%basis = from%basis
  to%basis_set = from%basis_set
end subroutine assign_hamil_descriptor


!******************************************************************************** 
!
!       Small routine to determine number of electrons
!
!******************************************************************************** 
subroutine set_electrons_struc(ham_des, struc)
  class(hamil_descriptor), intent(inout) :: ham_des
  type(Structure), intent(in)            :: struc
  !
  if (ham_des%nelect == unset) ham_des%nelect = nint(sum(struc%Zeff))
  if (ham_des%spin == unset) ham_des%spin = mod(ham_des%nelect, 2)
  if (ham_des%ispin==unset) then 
    if (mod(ham_des%nelect,2)==0) then
      ham_des%ispin = 1
    else
      ham_des%ispin = 2
    end if
  end if
  !
  if (ham_des%nelect <= 0) then
    if (comm_world%mpirank == 0) write(*,*) "Number of electrons can not be negative or zero - program stops "
    call exit
  end if
  !
  call ham_des%set_electrons_spin(ham_des%nelect, ham_des%spin)
end subroutine set_electrons_struc

subroutine set_electrons_spin(ham_des, nelect, spin)
  class(hamil_descriptor), intent(inout) :: ham_des
  integer, intent(in)                    :: nelect, spin
  
  ham_des%nelect = nelect
  ham_des%spin = spin
  ham_des%nel(1) = ceiling((nelect + spin) / 2.0_wp)
  ham_des%nel(2) = nelect - ham_des%nel(1)
  ham_des%nell(1) = 0
  ham_des%nell(2) = ham_des%nel(1)

  if (ham_des%ispin == 1) then
    ham_des%nocc = ham_des%nel(1)
  else
    ham_des%nocc = sum(ham_des%nel)
  end if
  
  ham_des%nocc_ = ham_des%nocc

  if (ham_des%ispin==1 .and. ham_des%spin==0) then
      ham_des%ispin_fock = 1
      ham_des%rspin = 2
  else
      ham_des%ispin_fock = 2
      ham_des%rspin = 1
  end if
end subroutine set_electrons_spin

subroutine set_electrons_nel(ham_des, nel_a, nel_b)
  class(hamil_descriptor), intent(inout) :: ham_des
  integer, intent(in)                    :: nel_a, nel_b
  
  ham_des%nel(1) = nel_a
  ham_des%nel(2) = nel_b
  ham_des%nell(1) = 0
  ham_des%nell(2) = ham_des%nel(1)
  ham_des%nelect = nel_a + nel_b
  ham_des%spin = nel_a - nel_b

  if (ham_des%ispin == 1) then
    ham_des%nocc = ham_des%nel(1)
  else
    ham_des%nocc = sum(ham_des%nel)
  end if
  
  ham_des%nocc_ = ham_des%nocc

  if (ham_des%ispin==1 .and. ham_des%spin==0) then
      ham_des%ispin_fock = 1
      ham_des%rspin = 2
  else
      ham_des%ispin_fock = 2
      ham_des%rspin = 1
  end if
end subroutine set_electrons_nel


!******************************************************************************** 
!
! Determine spin channel for wave functions
!
!******************************************************************************** 
integer function get_spin(ham_des, spin)
  class(hamil_descriptor), intent(in) :: ham_des
  integer, intent(in)                 :: spin
  
  get_spin = get_spin_channel(spin, ham_des%ispin)
end function get_spin


!******************************************************************************** 
!
! Determine spin channel for one-body operators
!
!******************************************************************************** 
integer function get_spin1(ham_des, spin)
  class(hamil_descriptor), intent(in) :: ham_des
  integer, intent(in)                 :: spin

  get_spin1 = get_spin1_channel(spin, ham_des%ispin1)
end function get_spin1


!******************************************************************************** 
!
! Determine spin channel for two-body operators
!
!******************************************************************************** 
integer function get_spin2(ham_des, spin1, spin2)
  class(hamil_descriptor), intent(in) :: ham_des
  integer, intent(in)                 :: spin1
  integer, intent(in)                 :: spin2
  
  get_spin2 = get_spin2_channel(spin1, spin2, ham_des%ispin2)
end function get_spin2

end module hamilton_vars