! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module frozen_core_approx

use mpi
use qmcfort_pos
use energy_types
use hamilton_layout
use hamilton_type
use method_base
use hfproc

private
public frozen_core

contains

!******************************************************************************** 
!
! Frozen-core approximation
!
!    Calculate core-valence interaction from frozen-core electrons and project
!    to the Hamiltonian in active space
!
! 
!******************************************************************************** 
subroutine frozen_core(struc, ham, add_pot, e_frozen)
  type(Structure), intent(in)         :: struc
  type(Hamilton), intent(inout)       :: ham
  logical, intent(in)                 :: add_pot
  type(Energy), optional, intent(out) :: e_frozen
  !local variables
  integer                             :: nfrozen
  real(wp), allocatable               :: occ(:,:)
  type(Energy)                        :: e_frozen_
  type(QmcFortMethod)                 :: method
  
  method = QmcFortMethod(method="frozen_core", def_active=.false., basis="mo", integral_mode="any")
  if (.not. method%is_active()) return
  call ham%meet_method_req(method)

  call detect_core_states(struc, ham%des, nfrozen)
  if (nfrozen == 0) return

  if (comm_world%mpirank == 0) call report_frozen_core_setup(nfrozen, add_pot)

  call set_occupations_simple(ham%des%n, [nfrozen, nfrozen], occ)

  e_frozen_ = Energy(enuc=0.0_wp, e0=ham%h0)
  call calculate_gse_full(ham, e_frozen_, occ)

  call core_valence_pot(ham, occ)
  call project_fs_ham(ham)

  if (add_pot) then      
    ham%h0 = ham%h0 + e_frozen_%electron_energy()
    ham%h1 = ham%h1 + ham%hpot - ham%xpot
    ham%h1_orig = ham%h1
  end if

  if (comm_world%mpirank == 0) call report_frozen_core(e_frozen_)

  if (present(e_frozen)) e_frozen = e_frozen_
end subroutine frozen_core


!******************************************************************************** 
!
! Determine the number of frozen core states to perform fc projection
! 
!    1. look at the number specified in the input file
!
!    2. calculate number of frozen states using cores of each atom
!       (defined in get_atomic_core_states)
!
!    3. if the nuclear charge does not agree with the number of electrons,
!       set ngrozen to 0
!
!******************************************************************************** 
subroutine detect_core_states(struc, hdes, nfrozen)
  type(Structure), intent(in)           :: struc
  type(hamil_descriptor), intent(inout) :: hdes
  integer, intent(out)                  :: nfrozen

  if (hdes%nfrozen > 0) then
    !nfrozen is explicitly requested 
    nfrozen = hdes%nfrozen
  else
    if (int(sum(struc%Zeff)) == sum(hdes%nel)) then
      nfrozen = struc%get_core_states()
    else
      nfrozen = 0
    end if

    !debug: this is not ideal on the long run
    hdes%nfrozen = nfrozen
  end if
end subroutine detect_core_states


!******************************************************************************** 
!
! Calculate core-valence interaction in full space (F+A) using HF routines
!
!******************************************************************************** 
subroutine core_valence_pot(ham, occ)
  type(Hamilton), intent(inout) :: ham
  real(wp), intent(in)          :: occ(:,:)
  !local variables
  real(wp), allocatable         :: pot(:,:,:), rdm_frozen(:,:,:)
  
  call calculate_rdm(ham%phi, occ, rdm_frozen)

  allocate(ham%hpot(ham%des%n, ham%des%n, ham%des%ispin1))
  call ham%hartree_potential_rdm(rdm_frozen, pot)
  ham%hpot = pot(:,:,1:ham%des%ispin1)

  allocate(ham%xpot(ham%des%n, ham%des%n, ham%des%ispin1))
  call ham%exchange_potential_rdm(rdm_frozen, pot)
  ham%xpot = pot(:,:,1:ham%des%ispin1)
end subroutine core_valence_pot


!******************************************************************************** 
!
! Project Hamiltonian onto active orbitals
!
! Trivial routines assuming the mean-field basis
!
!******************************************************************************** 
subroutine project_fs_ham(ham)
  use cholesky, only: compress_h2_gres
  type(Hamilton), intent(inout)    :: ham
  !local variables
  real(wp), allocatable            :: temp2(:,:), temp3(:,:,:)
  real(wp), allocatable            :: temp4(:,:,:,:), temp4_(:,:,:,:), temp5(:,:,:,:,:)
  real(wp), pointer                :: h2_gres_ob(:,:,:,:)
  integer                          :: nf, nold
  
  nold = ham%des%n
  nf = ham%des%nfrozen
  ham%des%n = ham%des%n - ham%des%nfrozen
  
  if (allocated(ham%eigenval)) then
    allocate(temp2(ham%des%n, ham%des%ispin))
    temp2 = ham%eigenval(nf+1:,:)
    call move_alloc(temp2, ham%eigenval)
  end if
  
  if (allocated(ham%occ)) then
    allocate(temp2(ham%des%n, ham%des%ispin))
    temp2 = ham%occ(nf+1:,:)
    call move_alloc(temp2, ham%occ)
  end if
  
  if (allocated(ham%h1)) then
    allocate(temp3(ham%des%n, ham%des%n, ham%des%ispin1))
    temp3 = ham%h1(nf+1:,nf+1:,:)
    call move_alloc(temp3, ham%h1)
  end if
  
  if (allocated(ham%h1_orig)) then
    allocate(temp3(ham%des%n, ham%des%n, ham%des%ispin1))
    temp3 = ham%h1_orig(nf+1:,nf+1:,:)
    call move_alloc(temp3, ham%h1_orig)
  end if
  
  if (allocated(ham%hpot)) then
    allocate(temp3(ham%des%n, ham%des%n, ham%des%ispin1))
    temp3 = ham%hpot(nf+1:,nf+1:,:)
    call move_alloc(temp3, ham%hpot)
  end if

  if (allocated(ham%xpot)) then
    allocate(temp3(ham%des%n, ham%des%n, ham%des%ispin1))
    temp3 = ham%xpot(nf+1:,nf+1:,:)
    call move_alloc(temp3, ham%xpot)
  end if

  if (allocated(ham%overlap)) then
    allocate(temp3(ham%des%n, ham%des%n, ham%des%ispin1))
    temp3 = ham%overlap(nf+1:,nf+1:,:)
    call move_alloc(temp3, ham%overlap)
  end if
  
  if (allocated(ham%phi)) then
    allocate(temp3(ham%des%n, ham%des%n, ham%des%ispin1))
    temp3 = ham%phi(nf+1:,nf+1:,:)
    call move_alloc(temp3, ham%phi)
  end if
  
  if (allocated(ham%h2)) then
    allocate(temp5(ham%des%n, ham%des%n, ham%des%n, ham%des%n, ham%des%ispin2))
    temp5 = ham%h2(nf+1:,nf+1:,nf+1:,nf+1:,:)
    call move_alloc(temp5, ham%h2)
    call ham%associate_h2_ptr()
  end if
  
  if (allocated(ham%h2_gres)) then
    allocate(temp3(ham%des%n*ham%des%n, ham%des%ng, ham%des%ispin1))
    call ham%get_h2_gres_ob(h2_gres_ob)
    temp3 = reshape(h2_gres_ob(nf+1:,nf+1:,:,:), shape=[ham%des%n*ham%des%n, ham%des%ng, ham%des%ispin1])
    call move_alloc(temp3, ham%h2_gres)
  end if
  
  ham%des%nbtot = ham%des%n
  ham%des%nel = ham%des%nel - ham%des%nfrozen
  ham%des%nfrozen = 0
  call ham%des%set_electrons_nel(ham%des%nel(1), ham%des%nel(2))
end subroutine project_fs_ham


!******************************************************************************** 
!
! Report the setup of the frozen-core approximation
!
!******************************************************************************** 
subroutine report_frozen_core_setup(nfrozen, add_pot)
  integer, intent(in) :: nfrozen
  logical, intent(in) :: add_pot

  call report_fc_setup_to_file(io%screen%funit)
  call report_fc_setup_to_file(io%qmcfort_log%funit)
  call report_fc_setup_to_file(io%qmcfort_out%funit)
contains
  subroutine report_fc_setup_to_file(funit)
    integer, intent(in) :: funit

    write(funit,*) "Frozen-core Approximation"
    write(funit,*) "-------------------------"
    write(funit,100) "number of frozen states", nfrozen
    write(funit,101) "add core-valence potential", add_pot
  
    100 format (1x,t5,a,t50,"= ",i6)
    101 format (1x,t5,a,t50,"= ",l)
  end subroutine report_fc_setup_to_file
end subroutine report_frozen_core_setup


!******************************************************************************** 
!
! Report the result of the frozen-core approximation
!
!******************************************************************************** 
subroutine report_frozen_core(frozen_energy)
  type(Energy), intent(in) :: frozen_energy

  call report_fc_to_file(io%screen)
  call report_fc_to_file(io%qmcfort_log)
  call report_fc_to_file(io%qmcfort_out)
contains
  subroutine report_fc_to_file(fh)
    type(FileHandle), intent(in) :: fh

    call frozen_energy%report(fh, method="frozen_core")
    write(fh%funit,*) starlong
  end subroutine report_fc_to_file
end subroutine report_frozen_core

end module frozen_core_approx
