! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module hfproc

#include "preproc.inc"

use constants
use hamilton_operations
use lapack
use method_base
use mpi
use profiling
use qmcfort_io
use qmcfort_pos
use standalone
use brng_factories, only: brng_glb
use diis_mixer, only: DiisConfig, DiisMixer
use energy_types, only: Energy
use file_handle, only: FileHandle
use hamilton_type, only: Hamilton
use hamilton_vars, only: hamil_descriptor
use orth, only: canonical_orthogonalization, symmetric_orthogonalization

implicit none

public

type scf_loop
  logical  :: iterate
  logical  :: converged
  integer  :: n
  integer  :: step
  integer  :: maxstep          = 50
  real(wp) :: en_tol           = 1.0E-10_wp
  real(wp) :: rdm_tol          = 1.0E-08_wp
  real(wp) :: energy
  real(wp) :: energy_old
  real(wp) :: ediff 
  real(wp) :: rdm_norm
  real(wp) :: rdm_norm_old
  real(wp) :: rdmdiff
end type scf_loop
  
type hf_descriptor
  type(scf_loop)         :: hf_loop
  character(len=charlen) :: mode           = "rhf"
  character(len=charlen) :: init_guess     = "none"
  integer                :: nshells
  integer                :: shell_low(3)
  integer                :: shell_len(3)
  real(wp)               :: shell_occ(3)
  real(wp)               :: uhf_mix        = 0.1_wp
  real(wp)               :: rdm_mix        = 0.8_wp
  real(wp)               :: s_square       = 0.0_wp
  real(wp)               :: s_mult         = 1.0_wp
  real(wp)               :: alpha_x        = 1.0_wp
  logical                :: gauss_proj     = .false.
  logical                :: lfrac          = .false.
  logical                :: spin_relax     = .false.
contains
  procedure              :: init   => init_hf_descriptor
  procedure              :: reader => hf_reader
end type hf_descriptor

type hf_hamil
  type(hf_descriptor)                      :: hf_des
  type(Energy)                             :: hf_en
  type(DiisMixer)                          :: diis
  type(Hamilton), pointer                  :: ham
  real(wp), allocatable                    :: occ_dens(:,:)
  real(wp), allocatable                    :: proj(:,:,:)
  procedure(init_guess_interface), pointer :: init_guess => null()
contains
  procedure                                :: scf_proc => hf_scf_proc
  procedure                                :: set => set_hf_hamil
  procedure                                :: print_header => print_hf_header
  procedure                                :: print_footer => print_hf_footer
  procedure                                :: print_info => print_hf_info
  procedure                                :: mulliken_analysis
  procedure                                :: associate_init_guess_ptr
  procedure                                :: get_hf_mode
end type hf_hamil

abstract interface
  subroutine init_guess_interface(hf_ham)
    import hf_hamil
    class(hf_hamil), intent(inout) :: hf_ham
  end subroutine init_guess_interface
end interface

contains

!********************************************************************************
!
!       Hartree Fock calculation - driver routine for hartree-fock method
!
!********************************************************************************
subroutine do_hartree_fock(ham, struc)
  type(Hamilton), intent(inout) :: ham
  type(structure), intent(in)   :: struc
  !local
  type(QmcFortMethod)           :: method
  type(hf_hamil)                :: hf_ham
  character(len=*), parameter   :: proc_name = "do_hartree_fock"
  
  method = QmcFortMethod(method="hfock", def_active=.false., basis="any", integral_mode="any")
  if (.not. method%is_active()) return
  call ham%meet_method_req(method)

  if (use_profiling) call start_profiling(proc_name)

  call hf_ham%set(struc, ham)
  call hf_ham%print_header
  call hf_ham%init_guess
  call hf_ham%scf_proc
  call hf_ham%print_footer

  if (use_profiling) call end_profiling(proc_name)
end subroutine do_hartree_fock


!******************************************************************************** 
!
!       Initialize hf descriptor from main
!
!******************************************************************************** 
subroutine init_hf_descriptor(hf_des, ham_des)
  class(hf_descriptor), intent(inout) :: hf_des
  type(hamil_descriptor), intent(in)  :: ham_des
  !
  call hf_des%reader
  !
  if (hf_des%hf_loop%maxstep <= 0) then
    hf_des%hf_loop%iterate = .false.
  else 
    hf_des%hf_loop%iterate = .true.
  end if
  hf_des%hf_loop%n = ham_des%n 
  hf_des%hf_loop%step = 0
  hf_des%hf_loop%energy = 1000.0_wp
  hf_des%hf_loop%energy_old = 0.0_wp
  !
  hf_des%nshells = 3
  hf_des%shell_occ = [2.0_wp, 1.0_wp, 0.0_wp]
  hf_des%shell_len(1) = min(ham_des%nel(1), ham_des%nel(2))
  hf_des%shell_len(2) = abs(ham_des%nel(1) - ham_des%nel(2))
  hf_des%shell_len(3) = ham_des%n - hf_des%shell_len(1) - hf_des%shell_len(2)
  hf_des%shell_low(1) = 0
  hf_des%shell_low(2) = hf_des%shell_len(1)
  hf_des%shell_low(3) = hf_des%shell_len(1) + hf_des%shell_len(2)
end subroutine init_hf_descriptor


!******************************************************************************** 
!
!       Read hf_des values from the qmcfort_in file
!
!******************************************************************************** 
subroutine hf_reader(hf_des)
  use qmcfort_in, only: add_input
  class(hf_descriptor), intent(inout) :: hf_des
   
  call add_input("hf_maxstep",      hf_des%hf_loop%maxstep)
  call add_input("hf_maxiter",      hf_des%hf_loop%maxstep)
  call add_input("uhf_mix",         hf_des%uhf_mix)
  call add_input("rdm_mix",         hf_des%rdm_mix)
  call add_input("rdm_eps",         hf_des%hf_loop%rdm_tol)
  call add_input("rdm_tol",         hf_des%hf_loop%rdm_tol)
  call add_input("hf_eps",          hf_des%hf_loop%en_tol)
  call add_input("hf_tol",          hf_des%hf_loop%en_tol)
  call add_input("proj",            hf_des%gauss_proj)
  call add_input("gauss_proj",      hf_des%gauss_proj)
  call add_input("projections",     hf_des%gauss_proj)
  call add_input("lfrac",           hf_des%lfrac)
  call add_input("spin_relax",      hf_des%spin_relax)
  call add_input("init_guess",      hf_des%init_guess)
  call add_input("scf_init_guess",  hf_des%init_guess)
  call add_input("alpha_x",         hf_des%alpha_x)
end subroutine hf_reader


!******************************************************************************** 
!
! Initializes hf_hamil and set pointers
!
!******************************************************************************** 
subroutine set_hf_hamil(hf_ham, struc, ham)
  class(hf_hamil), intent(inout)        :: hf_ham
  type(structure), intent(in)           :: struc
  type(Hamilton), target, intent(inout) :: ham
  !local
  integer                               :: param_size
  type(DiisConfig)                      :: diis_cfg

  
  if (.not. allocated(ham%eigenval)) allocate(ham%eigenval(ham%des%n, ham%des%ispin))
  hf_ham%ham => ham
  
  call hf_ham%get_hf_mode()
  call hf_ham%hf_des%init(ham%des)
  call hf_ham%associate_init_guess_ptr()

  !setup DiisConfig variable
  diis_cfg = DiisConfig()
  call diis_cfg%read_tags()

  !setup DiisMixer variable
  param_size = ham%des%ispin_fock * ham%des%n**2
  hf_ham%diis = DiisMixer(param_size, diis_cfg)
  
  hf_ham%hf_en%enuc = struc%calc_enuc()
end subroutine set_hf_hamil


!******************************************************************************** 
!
!       Associate scf_loop pointer
!
!******************************************************************************** 
subroutine get_hf_mode(hf_ham)
  class(hf_hamil), intent(inout) :: hf_ham
  !
  if (hf_ham%ham%des%ispin == 1) then
    if (hf_ham%ham%des%spin == 0) then
      hf_ham%hf_des%mode = "rhf"
    else
      hf_ham%hf_des%mode = "rohf"
    end if
  else
    hf_ham%hf_des%mode = "uhf"
  end if
end subroutine get_hf_mode


!******************************************************************************** 
!
!       Associate init_scf pointer
!
!******************************************************************************** 
subroutine associate_init_guess_ptr(hf_ham)
  class(hf_hamil), intent(inout) :: hf_ham
  !
  select case (hf_ham%hf_des%init_guess)
    case ("no", "none")
      hf_ham%init_guess => scf_init_guess_none
    case ("h1", "h1core", "core")
      hf_ham%init_guess => scf_init_guess_h1core 
    case ("diagonal", "diag")
      hf_ham%init_guess => scf_init_guess_diag
    case ("random")
      hf_ham%init_guess => scf_init_guess_random
    case default
      hf_ham%init_guess => scf_init_guess_none
  end select
end subroutine associate_init_guess_ptr


!********************************************************************************
!
!       Self-Consistent Hartree-Fock loop (rhf or uhf)
!
!********************************************************************************
subroutine hf_scf_proc(hf_ham)
  class(hf_hamil), intent(inout) :: hf_ham
  !local variables
  real(wp), allocatable          :: fock(:,:,:), fock_eff(:,:,:)
  real(wp), allocatable          :: rdm(:,:,:), rdm_old(:,:,:)
  real(wp), allocatable          :: coulomb(:,:,:), exchange(:,:,:)
  real(wp), allocatable          :: x(:,:,:,:)
  
  call canonical_orthogonalization(hf_ham%ham%overlap, x)
  call aufbau_principle(hf_ham%ham%des, hf_ham%ham%occ, hf_ham%ham%eigenval, hf_ham%hf_des%lfrac, hf_ham%hf_des%spin_relax)
  call initialize_rdm(hf_ham%ham%phi, hf_ham%ham%occ, rdm, hf_ham%hf_des%uhf_mix)
  call hf_ham%ham%hartree_potential_rdm(rdm, coulomb)
  call hf_ham%ham%exchange_potential_rdm(rdm, exchange)
  call calculate_gse(hf_ham%ham%des, hf_ham%ham%h1, coulomb, exchange, rdm, hf_ham%hf_en)
  call init_scf_convergence(hf_ham%hf_des%hf_loop, hf_ham%hf_en%total_energy())
  call hf_ham%print_info()
  
  do while (hf_ham%hf_des%hf_loop%iterate)
    call calculate_fock_matrix(hf_ham%ham%des, hf_ham%hf_des, hf_ham%ham%h1, coulomb, exchange, fock)
    
    !use diis to get the best 
    if (hf_ham%diis%cfg%enabled) call diis_interface(hf_ham, fock, rdm, x)

    call calculate_fock_eff(hf_ham%ham%des, hf_ham%ham%phi, hf_ham%ham%occ, rdm, fock, fock_eff)
    call transform_fock_general(hf_ham%ham%des, x, fock_eff)
    call diagonalize_fock_matrix(fock_eff, hf_ham%ham%eigenval)
    call transform_orbitals(x, fock_eff, hf_ham%ham%phi)    
    call normalize_orbitals(hf_ham%ham%des, hf_ham%ham%overlap, hf_ham%ham%phi)
    call aufbau_principle(hf_ham%ham%des, hf_ham%ham%occ, hf_ham%ham%eigenval, hf_ham%hf_des%lfrac, hf_ham%hf_des%spin_relax)
    call move_alloc(rdm, rdm_old)
    call calculate_rdm(hf_ham%ham%phi, hf_ham%ham%occ, rdm)
    call rdm_mixer(hf_ham%hf_des, hf_ham%diis%cfg%enabled, rdm, rdm_old)
    call hf_ham%ham%hartree_potential_rdm(rdm, coulomb)
    call hf_ham%ham%exchange_potential_rdm(rdm, exchange)
    call calculate_gse(hf_ham%ham%des, hf_ham%ham%h1, coulomb, exchange, rdm, hf_ham%hf_en)
    call check_scf_convergence(hf_ham%hf_des%hf_loop, hf_ham%hf_en%total_energy(), rdm(:,:,4))
    call hf_ham%print_info
  end do
  
  call hf_ham%mulliken_analysis(rdm)
  hf_ham%hf_des%s_square = estimate_s_square(hf_ham%ham%phi, hf_ham%ham%overlap, hf_ham%ham%des%nel)
  hf_ham%hf_des%s_mult = sqrt(1.0_wp + 4.0_wp*hf_ham%hf_des%s_square)
end subroutine hf_scf_proc


!******************************************************************************** 
!
!       Initializes scf loop after first energy evaluation
!
!******************************************************************************** 
subroutine init_scf_convergence(hf_loop, energy)
  type(scf_loop), intent(inout) :: hf_loop
  real(wp), intent(in)          :: energy
  !
  hf_loop%energy_old = hf_loop%energy
  hf_loop%energy = energy
  hf_loop%ediff = hf_loop%energy - hf_loop%energy_old  
  hf_loop%rdm_norm = 1.0_wp
  hf_loop%rdm_norm_old = 0.0_wp
end subroutine init_scf_convergence


!******************************************************************************** 
!
!       Checks convergence of the SCCF procedure 
!
!******************************************************************************** 
subroutine check_scf_convergence(hf_loop, energy, rdm)
  type(scf_loop), intent(inout) :: hf_loop
  real(wp), intent(in)          :: energy
  real(wp), intent(in)          :: rdm(:,:)
  !local variables
  logical                       :: try1, try2, try1_g, try2_g
  character(len=*), parameter   :: proc_name = "check_scf_convergence"
  
  if (use_profiling) call start_profiling(proc_name)

  hf_loop%step = hf_loop%step + 1
  hf_loop%energy_old = hf_loop%energy
  hf_loop%energy = energy
  hf_loop%ediff = abs(hf_loop%energy - hf_loop%energy_old)
  try1 = hf_loop%ediff < hf_loop%en_tol
  
  hf_loop%rdm_norm_old = hf_loop%rdm_norm
  hf_loop%rdm_norm = norm2(rdm)
  hf_loop%rdmdiff = abs(hf_loop%rdm_norm - hf_loop%rdm_norm_old)
  try2 = hf_loop%rdmdiff < hf_loop%rdm_tol

  call comm_world%reduce(try1, op=MPI_Land)
  call comm_world%reduce(try2, op=MPI_Land)

  hf_loop%converged = try1 .and. try2
  hf_loop%iterate = .not. hf_loop%converged
  if (hf_loop%step == hf_loop%maxstep) hf_loop%iterate = .false.

  if (use_profiling) call end_profiling(proc_name)
end subroutine check_scf_convergence


!******************************************************************************** 
!
!       Pass routine - doesn't change orbitals
!
!               if orbitals  are zero 
!                       equivalent to init_guess_h1core
!
!               if the file orbitals is present
!                       orbitals read from file
!
!******************************************************************************** 
subroutine scf_init_guess_none(hf_ham)
  class(hf_hamil), intent(inout) :: hf_ham
end subroutine scf_init_guess_none


!******************************************************************************** 
!
!       Orbitals set to zero - core hamiltonian used as the initial guess
!
!******************************************************************************** 
subroutine scf_init_guess_h1core(hf_ham)
  class(hf_hamil), intent(inout) :: hf_ham
  !
  hf_ham%ham%phi = 0.0_wp
end subroutine scf_init_guess_h1core


!******************************************************************************** 
!
!       Orbitals set to gaussian basis functions
!
!******************************************************************************** 
subroutine scf_init_guess_diag(hf_ham)
  class(hf_hamil), intent(inout) :: hf_ham
  !local variables
  integer                        :: spin
  !
  do spin = 1, size(hf_ham%ham%phi, 3)
    call eye(hf_ham%ham%phi(:,:,spin))
  end do
end subroutine scf_init_guess_diag


!******************************************************************************** 
!
!       Random orbitals used as the initial guess
!
!******************************************************************************** 
subroutine scf_init_guess_random(hf_ham)
  class(hf_hamil), intent(inout) :: hf_ham
  !local
  integer                        :: n, m, ispin, spin
  real(wp), allocatable          :: rand_phi(:,:,:), rmat(:,:)

  n = size(hf_ham%ham%phi, 1)
  m = size(hf_ham%ham%phi, 2)
  ispin = size(hf_ham%ham%phi, 3)

  do spin = 1, ispin
    call eye(hf_ham%ham%phi(:,:,spin))
  end do

  allocate(rand_phi, mold=hf_ham%ham%phi)
  call brng_glb%rand(rand_phi)
  rand_phi = hf_ham%ham%phi + 0.1_wp * (rand_phi-0.5_wp)

  do spin = 1, ispin
    if (allocated(rmat)) deallocate(rmat)
    allocate(rmat(m,m))
    call getqr(rand_phi(:,:,spin), rmat)
  end do

  hf_ham%ham%phi = rand_phi
  !call transform_orbitals(x, rand_phi, hf_ham%ham%phi)
end subroutine scf_init_guess_random


!******************************************************************************** 
!
!       Initializes reduced density matrix in rhf/uhf model
!
!******************************************************************************** 
subroutine initialize_rdm(phi, occ, rdm, uhf_mix)
  real(wp), intent(inout)             :: phi(:,:,:)
  real(wp), intent(in)                :: occ(:,:)
  real(wp), allocatable,  intent(out) :: rdm(:,:,:)
  real(wp), optional, intent(in)      :: uhf_mix
  !
  call calculate_rdm(phi, occ, rdm)
  !
  if (size(phi,3)==2 .and. present(uhf_mix)) then
    rdm(:,:,1) = rdm(:,:,1) + uhf_mix * rdm(:,:,2)
    rdm(:,:,2) = rdm(:,:,2) - uhf_mix * rdm(:,:,2)
    rdm(:,:,3) = rdm(:,:,1) - rdm(:,:,2)
  end if
end subroutine initialize_rdm


!********************************************************************************
!
!       Calculate Hartree-Fock reduced density matrix
!
!                      rdm = phi * phi^{T}
!
!                       rdm(:,:,1)      - spin up density matrix
!                       rdm(:,:,2)      - spin down density matrix
!                       rdm(:,:,3)      - spin density matix     rdm1-rdm2
!                       rdm(:,:,4)      - total density matrix   rdm1+rdm2
!
!********************************************************************************
subroutine calculate_rdm(phi, occ, rdm)
  real(wp), intent(in)                :: phi(:,:,:)
  real(wp), intent(in)                :: occ(:,:)
  real(wp), allocatable,  intent(out) :: rdm(:,:,:)
  !local variables 
  real(wp), parameter                 :: tol = 1.0e-06_wp
  real(wp)                            :: fact
  integer                             :: n, spin, spin_, nocc, g, ng
  integer, allocatable                :: low(:), up(:)
  character(len=*), parameter         :: proc_name = "calculate_rdm"
  
  if (use_profiling) call start_profiling(proc_name)

  n = size(phi, 1)

  if (.not. allocated(rdm)) allocate(rdm(n,n,4))
  rdm = 0.0_wp

  do spin = 1, 2
    spin_ = min(spin, size(phi,3))
    call get_occ_groups(occ(:,spin), low, up, ng)
    
    do g = 1, ng
      nocc = up(g) - low(g) + 1
      if (nocc == 0) cycle
      fact = occ(low(g), spin)
      if (fact < tol) cycle
      call gemm("n", "t", n, n, nocc, fact, phi(:,low(g):up(g),spin_), n, phi(:,low(g):up(g),spin_), n, 1.0_wp, rdm(:,:,spin), n)
    end do
  end do
  
  rdm(:,:,3) = rdm(:,:,1) - rdm(:,:,2)
  rdm(:,:,4) = rdm(:,:,1) + rdm(:,:,2)

  if (use_profiling) call end_profiling(proc_name)
end subroutine calculate_rdm


!******************************************************************************** 
!
!       Transofrm from spin to shell rdm representation
!
!               rdm_shell(1) = rdm_spin(2)
!               rdm_shell(2) = rdm_spin(1) - rdm_spin(2)
!               rdm_shell(3) = rdm_spin(3)
!               rdm_shell(4) = rdm_spin(4)
!
!******************************************************************************** 
subroutine spin_to_shell_rdm(rdm)
  real(wp), allocatable, intent(inout) :: rdm(:,:,:)
  !local variables
  real(wp), allocatable                :: rdm_(:,:,:)
  !
  allocate(rdm_, source=rdm)
  rdm_(:,:,1) = rdm(:,:,2)
  rdm_(:,:,2) = rdm(:,:,1) - rdm(:,:,2)
  call move_alloc(rdm_, rdm)
end subroutine spin_to_shell_rdm


!******************************************************************************** 
!
!       Transofrm from shell to spin rdm representation
!
!               rdm_spin(1) = rdm_shell(1) + rdm_shell(2)
!               rdm_spin(2) = rdm_shell(1)
!               rdm_spin(3) = rdm_shell(3)
!               rdm_spin(4) = rdm_shell(4)
!
!******************************************************************************** 
subroutine shell_to_spin_rdm(rdm)
  real(wp), allocatable, intent(inout) :: rdm(:,:,:)
  !local variables
  real(wp), allocatable                :: rdm_(:,:,:)
  !
  allocate(rdm_, source=rdm)
  rdm_(:,:,1) = rdm(:,:,1) + rdm(:,:,2)
  rdm_(:,:,2) = rdm(:,:,1)
  call move_alloc(rdm_, rdm)
end subroutine shell_to_spin_rdm


!******************************************************************************** 
!
!       Calculate rdm of the unoccupied (empty) orbitals 
!
!               rdm_empty = C_e C_e^{\dagger}
!
!       C_e - orbitals of the empty shell
!
!******************************************************************************** 
subroutine calculate_rdm_empty(phi, occ, rdm)
  real(wp), intent(in)               :: phi(:,:,:) 
  real(wp), allocatable, intent(in)  :: occ(:,:)
  real(wp), allocatable, intent(out) :: rdm(:,:)
  !local variables
  real(wp), allocatable              :: temp(:,:,:), occ_(:,:)
  integer                            :: n
  !
  n = size(phi, 1)
  if (.not. allocated(rdm)) allocate(rdm(n,n))
  allocate(occ_, source=occ)
  occ_ = 1.0_wp - occ_
  call calculate_rdm(phi, occ_, temp)
  rdm = temp(:,:,1)
end subroutine calculate_rdm_empty


!******************************************************************************** 
!
! DIIS interface subroutine - responsibilities:
! 
! 1) Pparameter vector x, and error vector e are prepared for DIIS
!    x, and e are defined as follows:
!        x = F
!        e = FRS - (FRS)^{\dager} = FRS - SRF  (in AO basis)
!
!        F - Fock matrix in the AO basis
!        R - Density matrix in the AO basis
!        S - AO overlap matrix
!
!    For better stability the error vector is transformed to the MO basis 
!        e_AO --> e_MO = X^{\dagger} e_AO X
!
!        X - orthogonlization matrix 
!
! 2) DIIS space is updated using x and e, and finally
!
! 3) DIIS equations are solved to obtain new optimal F = \sum_i c_i F_i
!
!******************************************************************************** 
subroutine diis_interface(hf_ham, fock, rdm, x)
  type(hf_hamil), intent(inout)               :: hf_ham
  real(wp), contiguous, target, intent(inout) :: fock(:,:,:)
  real(wp), allocatable, intent(inout)        :: rdm(:,:,:)
  real(wp), intent(in)                        :: x(:,:,:,:)
  !local
  integer                                     :: n, spin, spin1, ispin
  real(wp)                                    :: rspin(2)
  real(wp), allocatable, target               :: e(:,:,:)
  real(wp), allocatable                       :: temp(:,:), temp2(:,:)
  real(wp), pointer                           :: x_flat(:), e_flat(:)
  
  n = size(fock, 1)
  ispin = size(fock, 3)

  allocate(e(n,n,ispin), temp(n,n), temp2(n,n))
  e = 0.0_wp
  
  !For ROHF, go to the shell density matrices
  if (hf_ham%ham%des%ispin==1 .and. hf_ham%ham%des%spin/=0) then
    call spin_to_shell_rdm(rdm)
    rspin = [2.0_wp, 1.0_wp]
  else
    rspin = [1.0_wp, 1.0_wp]
  end if

  !calculate error vectors e = X^{\dagger} (FRS - SRF) X
  do spin = 1, ispin
    spin1 = hf_ham%ham%des%get_spin1(spin)
    call gemm("n", "n", n, n, n, rspin(spin), rdm(:,:,spin), n, hf_ham%ham%overlap(:,:,spin1), n, 0.0_wp, temp, n)
    call gemm("n", "n", n, n, n, 1.0_wp, fock(:,:,spin), n, temp, n, 0.0_wp, temp2, n)
    call gemm("n", "n", n, n, n, 1.0_wp, temp2, n, x(:,:,1,spin1), n, 0.0_wp, temp, n)
    call gemm("t", "n", n, n, n, 1.0_wp, x(:,:,1,spin1), n, temp, n, 0.0_wp, temp2, n)
    e(:,:,spin) = e(:,:,spin) + temp2 - transpose(temp2)
  end do

  !For ROHF, transform back to the spin density matrices
  if (hf_ham%ham%des%ispin==1 .and. hf_ham%ham%des%spin/=0) call shell_to_spin_rdm(rdm)
  
  !make 1d view of the parameter and error vectors 
  x_flat(1:size(fock)) => fock
  e_flat(1:size(e)) => e

  call hf_ham%diis%update(x_flat, e_flat)
  call hf_ham%diis%solve(x_flat)
end subroutine diis_interface


!******************************************************************************** 
!
!       Density mixing in SCF procedure
!
!******************************************************************************** 
subroutine rdm_mixer(hf_des, use_diis, rdm, rdm_old)
  type(hf_descriptor), intent(in) :: hf_des
  logical, intent(in)             :: use_diis
  real(wp), intent(inout)         :: rdm(:,:,:)
  real(wp), intent(in)            :: rdm_old(:,:,:)
  !local
  real(wp)                        :: mix

!debug: should have better solution on the long run
  if (use_diis) then 
    mix = 1.0_wp
  else
    mix = hf_des%rdm_mix
  end if

  rdm = mix*rdm + (1.0_wp-mix)*rdm_old
end subroutine rdm_mixer


!********************************************************************************
!
!       Normalize orbitals according to <phi|S|phi> = 1
!
!********************************************************************************
subroutine normalize_orbitals(ham_des, overlap, phi)
  type(hamil_descriptor), intent(in) :: ham_des
  real(wp), intent(in)               :: overlap(:,:,:)
  real(wp), intent(inout)            :: phi(:,:,:)
  !local variables
  integer                            :: p, n, spin, spin2
  real(wp)                           :: temp(ham_des%n)
  !
  n = ham_des%n
  !
  do spin = 1, ham_des%ispin
    spin2 = spin
    if (ham_des%ispin1 == 1) spin2 = 1
    do p = 1, n
      call gemv("N", n, n, 1.0_wp, overlap(:,:,spin2), n, phi(:,p,spin), 1, 0.0_wp, temp, 1)
      phi(:,p,spin) = phi(:,p,spin) / sqrt(dot_product(phi(:,p,spin), temp))
    end do
  end do
end subroutine normalize_orbitals


!******************************************************************************** 
!
!       Calculate Fock matrix
!
!       rhf/uhf model 
!
!               F1 = h + J(rdm1+rdm2) - K(rdm1)
!               F2 = h + J(rdm1+rdm2) - K(rdm2)
!
!               if spin unpolarized calculation F1 = F2
!
!       rohf model
!
!               F1 = h + J(rdm1+rdm2) - 0.5 * K
!               F2 = h + J(rdm1+rdm2) - K(rdm1)
!
!******************************************************************************** 
subroutine calculate_fock_matrix(ham_des, hf_des, ham1, coulomb, exchange, fock)
  type(hamil_descriptor), intent(in) :: ham_des
  type(hf_descriptor), intent(in)    :: hf_des
  real(wp), intent(in)               :: ham1(:,:,:)
  real(wp), intent(in)               :: coulomb(:,:,:)
  real(wp), intent(in)               :: exchange(:,:,:)
  real(wp), allocatable, intent(out) :: fock(:,:,:)
  !local
  integer                            :: spin, spin2
  !
  if (.not. allocated(fock)) allocate(fock(ham_des%n, ham_des%n, ham_des%ispin_fock))
  fock = 0.0_wp
  if (ham_des%ispin==1 .and. ham_des%spin/=0) then
    fock(:,:,1) = ham1(:,:,1) + coulomb(:,:,1) - 0.5_wp*exchange(:,:,1) - 0.5_wp*exchange(:,:,2)
    fock(:,:,2) = ham1(:,:,1) + coulomb(:,:,2) - exchange(:,:,1)
  else
    do spin = 1, size(fock, 3)
      spin2 = ham_des%get_spin1(spin)
      fock(:,:,spin) = ham1(:,:,spin2) + hf_des%alpha_x * (coulomb(:,:,spin) -  exchange(:,:,spin))
    end do
  end if
end subroutine calculate_fock_matrix


!********************************************************************************
!
!       Calculates single effective Fock Matrix in rohf model
!
!               F =  (R1+R2)(2F1-F2)(R1+R2) + (R1+R3)F1(R+R3) +
!                                 (R2+R3)F2(R2+R3)
!
!********************************************************************************
subroutine calculate_fock_eff(ham_des, phi, occ, rdm, fock, fock_eff)
  type(hamil_descriptor), intent(in)   :: ham_des
  real(wp), intent(in)                 :: phi(:,:,:)
  real(wp), allocatable, intent(in)    :: occ(:,:)
  real(wp), allocatable, intent(inout) :: rdm(:,:,:) 
  real(wp), allocatable, intent(in)    :: fock(:,:,:)
  real(wp), allocatable, intent(out)   :: fock_eff(:,:,:)
  !local variables
  real(wp), parameter                  :: tol = 1.0e-10_wp
  real(wp), allocatable                :: aux(:,:,:), prod(:,:), dfock(:,:,:), rdm_empty(:,:)
  integer                              :: i, n
  !
  if (allocated(fock_eff)) deallocate(fock_eff)
  !
  if (ham_des%ispin==1 .and. ham_des%spin/=0) then
    n = size(rdm, 1)
    allocate(fock_eff(n,n,1), aux(n,n,3), prod(n,n), dfock(n,n,3))
    call calculate_rdm_empty(phi, occ, rdm_empty)
    call spin_to_shell_rdm(rdm)
    fock_eff = 0.0_wp
    aux(:,:,1) = rdm(:,:,1) + rdm_empty
    aux(:,:,2) = rdm(:,:,2) + rdm_empty
    aux(:,:,3) = rdm(:,:,1) + rdm(:,:,2)
    call shell_to_spin_rdm(rdm)
    dfock(:,:,1) = fock(:,:,1)
    dfock(:,:,2) = fock(:,:,2)
    dfock(:,:,3) = 2.0_wp*fock(:,:,1) - fock(:,:,2)
    !
    do i = 1, 3
      call gemm("n","n",n,n,n,1.0_wp,dfock(:,:,i),n,aux(:,:,i),n,0.0_wp,prod,n)
      call gemm("n","n",n,n,n,1.0_wp,aux(:,:,i),n,prod,n,1.0_wp,fock_eff(:,:,1),n)
    end do
    !
    if (all(abs(fock_eff(:,:,1)) < tol)) fock_eff(:,:,1) = fock(:,:,1)
  else
    allocate(fock_eff, source=fock)
  end if
end subroutine calculate_fock_eff


!******************************************************************************** 
!
!       Diagonalize Fock matrix
!
!******************************************************************************** 
subroutine diagonalize_fock_matrix(fock, eigenval)
  real(wp), intent(inout)     :: fock(:,:,:)
  real(wp), intent(out)       :: eigenval(:,:)
  !local variables
  integer                     :: spin
  character(len=*), parameter :: proc_name = "diagonalize_fock_matrix"
  
  if (use_profiling) call start_profiling(proc_name)

  do spin = 1, size(fock, 3)
    call syev(fock(:,:,spin), eigenval(:,spin), "V")
  end do

  if (use_profiling) call end_profiling(proc_name)
end subroutine diagonalize_fock_matrix


!******************************************************************************** 
!
!       For RHF/UHF does the covariant transformation of the fock matrix with X
!
!       For ROHF method does contravariant trafo with X^{-1}
!
!******************************************************************************** 
subroutine transform_fock_general(ham_des, x, fock)
  type(hamil_descriptor), intent(in) :: ham_des
  real(wp), intent(in)               :: x(:,:,:,:)
  real(wp), intent(inout)            :: fock(:,:,:)
  !
  if (ham_des%ispin==1 .and. ham_des%spin/=0) then
    call transform_rdm(x, fock)
  else
    call transform_fock_matrix(x, fock)
  end if
end subroutine transform_fock_general


!******************************************************************************** 
!
!       Canonical transfromation of the Fock matrix:
!
!               F --> X^{+} F X
!
!******************************************************************************** 
subroutine transform_fock_matrix(x, fock)
  real(wp), intent(in)        :: x(:,:,:,:)
  real(wp), intent(inout)     :: fock(:,:,:)
  !local variables
  integer                     :: n, spin, spin2
  real(wp), allocatable       :: temp(:,:)
  character(len=*), parameter :: proc_name = "transform_fock_matrix"
  
  if (use_profiling) call start_profiling(proc_name)

  n = size(fock, 1)
  allocate(temp(n,n))
  
  do spin = 1, size(fock, 3)
    spin2 = spin
    if (size(x, 4) == 1) spin2 = 1
    call gemm("n", "n", n, n, n, 1.0_wp, fock(:,:,spin), n, x(:,:,1,spin2), n, 0.0_wp, temp, n)
    call gemm("t", "n", n, n, n, 1.0_wp, x(:,:,1,spin2), n, temp, n, 0.0_wp, fock(:,:,spin), n)
  end do

  if (use_profiling) call end_profiling(proc_name)
end subroutine transform_fock_matrix


!******************************************************************************** 
!
!       Canonical transfromation of the reduced density matrix
!
!             rdm --> X^{-1} rdm X^{-1 +}
!
!******************************************************************************** 
subroutine transform_rdm(x, rdm)
  real(wp), intent(in)        :: x(:,:,:,:)
  real(wp), intent(inout)     :: rdm(:,:,:)
  !local variables
  integer                     :: n, spin, spin2
  real(wp), allocatable       :: temp(:,:)
  character(len=*), parameter :: proc_name = "transform_rdm"
  
  if (use_profiling) call start_profiling(proc_name)

  n = size(rdm, 1)
  allocate(temp(n,n))
  
  do spin = 1, size(rdm, 3)
    spin2 = spin
    if (size(x, 4) == 1) spin2 = 1
    call gemm("n", "t", n, n, n, 1.0_wp, rdm(:,:,spin), n, x(:,:,2,spin2), n, 0.0_wp, temp, n)
    call gemm("n", "n", n, n, n, 1.0_wp, x(:,:,2,spin2), n, temp, n, 0.0_wp, rdm(:,:,spin), n)
  end do

  if (use_profiling) call end_profiling(proc_name)
end subroutine transform_rdm


!******************************************************************************** 
!
!       Canonical transfroamtion of the orbitals
!
!               phi --> X phi
!
!       Note initial orbitals are usually stored in the fock matrix
!
!******************************************************************************** 
subroutine transform_orbitals(x, fock, phi)
  real(wp), intent(in)        :: x(:,:,:,:)
  real(wp), intent(in)        :: fock(:,:,:)
  real(wp), intent(out)       :: phi(:,:,:)
  !local variables
  integer                     :: n, spin, spin2
  character(len=*), parameter :: proc_name = "transform_orbitals"
  
  if (use_profiling) call start_profiling(proc_name)

  n = size(fock, 1)

  do spin = 1, size(phi, 3)
    spin2 = spin
    if (size(x, 4) == 1) spin2 = 1
    call gemm("n", "n", n, n, n, 1.0_wp, x(:,:,1,spin2), n, fock(:,:,spin), n, 0.0_wp, phi(:,:,spin), n)
  end do

  if (use_profiling) call end_profiling(proc_name)
end subroutine transform_orbitals


!********************************************************************************
!
!       Calculate rhf/uhf ground state energy
!
!               E = Tr(h (rdm1+rdm2)) 
!                 + 0.5 * Tr((rdm1+rdm2) J(rdm1+rdm2))
!                 - 0.5 * Tr(rdm1 K(rdm1)) - 0.5 * Tr(rdm2 K(rdm2))
!
!********************************************************************************
subroutine calculate_gse(ham_des, ham1, coulomb, exchange, rdm, hf_en)
  type(hamil_descriptor), intent(in) :: ham_des
  real(wp), intent(in)               :: ham1(:,:,:)
  real(wp), intent(in)               :: coulomb(:,:,:)
  real(wp), intent(in)               :: exchange(:,:,:)
  real(wp), intent(in)               :: rdm(:,:,:)
  type(Energy), intent(inout)        :: hf_en
  !local variables
  integer                            :: spin, spin2
  character(len=*), parameter        :: proc_name = "calculate_gse"
  
  if (use_profiling) call start_profiling(proc_name)
  
  hf_en%e0 = ham_des%h0
  hf_en%e1 = 0.0_wp
  hf_en%eh = 0.0_wp
  hf_en%ex = 0.0_wp
  
  do spin = 1, ham_des%ispin_fock
    spin2 = ham_des%get_spin1(spin)
    hf_en%e1 = hf_en%e1 + ham_des%rspin * trace2(rdm(:,:,spin), ham1(:,:,spin2))
    hf_en%eh = hf_en%eh + ham_des%rspin * 0.5_wp * trace2(rdm(:,:,spin), coulomb(:,:,spin))
    hf_en%ex = hf_en%ex - ham_des%rspin * 0.5_wp * trace2(rdm(:,:,spin), exchange(:,:,spin))
  end do

  if (use_profiling) call end_profiling(proc_name)
end subroutine calculate_gse


!******************************************************************************** 
!
!       Calculate Energy of a determinant given by occupation numbers occ
!
!******************************************************************************** 
subroutine calculate_gse_full(ham, e, occ)
  type(Hamilton), intent(inout)               :: ham
  type(Energy), intent(inout)                 :: e
  real(wp), allocatable, optional, intent(in) :: occ(:,:)
  !local variables
  real(wp), allocatable                       :: rdm(:,:,:)
  real(wp), allocatable                       :: coulomb(:,:,:), exchange(:,:,:), occ_(:,:)
  
  if (present(occ)) then
    allocate(occ_, source=occ)
  else
    call aufbau_principle(ham%des, occ_, ham%eigenval)
  end if
  
  call calculate_rdm(ham%phi, occ_, rdm)
  call ham%hartree_potential_rdm(rdm, coulomb)
  call ham%exchange_potential_rdm(rdm, exchange)
  call calculate_gse(ham%des, ham%h1, coulomb, exchange, rdm, e)
end subroutine calculate_gse_full


!********************************************************************************
!
!       Mulliken Analysis of the Hartree Fock calculation
!
!******************************************************************************** 
subroutine mulliken_analysis(hf_ham, rdm)
  class(hf_hamil), intent(inout) :: hf_ham
  real(wp), intent(in)           :: rdm(:,:,:)
  
  call get_mulliken_occupations(hf_ham%ham%des, hf_ham%ham%overlap, rdm, hf_ham%occ_dens)
  if (hf_ham%hf_des%gauss_proj) call get_projections(hf_ham%ham%des, hf_ham%ham%phi, hf_ham%proj)
end subroutine mulliken_analysis


!******************************************************************************** 
!
!       Aufbau principle for orbitals
!               
!       Determines the spin state of the system, i.e. occupancies  according to
!               
!               eigenvalues if they are supplied
!
!               number of electrons assuming lowest states are populated
!
!       if lfrac = .true. fractional occ will be allowed (eigenval required)       
!
!       if relax = .true. spin state can change and ham_des is modified
!
!******************************************************************************** 
subroutine aufbau_principle(ham_des, occ, eigenval, lfrac, relax)
  type(hamil_descriptor), intent(inout)       :: ham_des
  real(wp), allocatable, intent(out)          :: occ(:,:)
  real(wp), allocatable, optional, intent(in) :: eigenval(:,:)
  logical, optional, intent(in)               :: lfrac
  logical, optional, intent(in)               :: relax
  !local variables
  real(wp), parameter                         :: tol = 1.0e-08_wp
  logical                                     :: lfrac_, relax_
  
  if (present(lfrac)) then
    lfrac_ = lfrac
  else
    lfrac_ = .false.
  end if
  
  if ( present(relax)) then
    relax_ = relax
  else
    relax_ = .false.
  end if
  
  if (present(eigenval) .and. relax_) then
    if (allocated(eigenval) .and. all(abs(eigenval) > tol)) then
      call set_occupations_eig(ham_des, occ, eigenval, lfrac_)
    else
      call set_occupations_simple(ham_des%n, ham_des%nel, occ)
    end if
  else
    call set_occupations_simple(ham_des%n, ham_des%nel, occ)
  end if
end subroutine aufbau_principle


!******************************************************************************** 
!
!       Simple routine to set occupancies according to the number of electrons 
!                                   per spin channel
!
!******************************************************************************** 
subroutine set_occupations_simple(n, nel, occ)
  integer, intent(in)                :: n
  integer, intent(in)                :: nel(2)
  real(wp), allocatable, intent(out) :: occ(:,:)
  !local variables
  integer                            :: i, spin
  
  if (.not. allocated(occ)) allocate(occ(n, 2))
  occ = 0.0_wp
  
  do spin = 1, 2
    do i = 1, nel(spin)
      occ(i, spin) = 1.0_wp
    end do
  end do
end subroutine set_occupations_simple


!******************************************************************************** 
!
!       Set occupation according to the supplied eigenvalues
!
!******************************************************************************** 
subroutine  set_occupations_eig(ham_des, occ, eigenval, lfrac)
  type(hamil_descriptor), intent(inout) :: ham_des
  real(wp), allocatable, intent(out)    :: occ(:,:)
  real(wp),  intent(in)                 :: eigenval(:,:)
  logical, intent(in)                   :: lfrac
  !local variables
  real(wp), parameter                   :: tol = 1.0e-06_wp
  integer                               :: low_a, up_a, low_b, up_b, counter, deg, spin
  integer                               :: spin_a, spin_b
  
  if (.not. allocated(occ)) allocate(occ(size(eigenval,1),2))
  occ = 0.0_wp
  
  counter = sum(ham_des%nel)
  spin_a = 1
  spin_b = size(eigenval, 2)
  low_a = 1
  low_b = 1
  
  do while (counter > 0)
    if (eigenval(low_a,spin_a)-eigenval(low_b,spin_b) <= tol) then
      call get_degeneracies(eigenval(:,spin_a), low_a, up_a, deg)
      if (deg <= counter) then
        occ(low_a:up_a-1,1) = 1.0_wp
        low_a = up_a
        counter = counter - deg
      else
        if (lfrac) then
          occ(low_a:up_a-1,1) = real(counter, wp) / real(deg, wp)
        else
          occ(low_a:low_a+counter-1,1) = 1.0_wp
        end if
        low_a = up_a
        counter = 0
      end if
    else
      call get_degeneracies(eigenval(:,spin_b), low_b, up_b, deg)
      if (deg <= counter) then
        occ(low_b:up_b-1,2) = 1.0_wp
        low_b = up_b
        counter = counter - deg
      else
        if (lfrac) then
          occ(low_b:up_b-1,2) = real(counter, wp) / real(deg, wp)
        else
          occ(low_b:low_b+counter-1,2) = 1.0_wp
        end if
        low_b = up_b
        counter = 0
      end if
    end if
  end do
   
  do spin = 1, 2
    ham_des%nel(spin) = nint(sum(occ(:,spin)))
  end do
  call ham_des%set_electrons_nel(ham_des%nel(1), ham_des%nel(2))
end subroutine set_occupations_eig


!******************************************************************************** 
!
!       Small helper routine to find the degenerate states
!
!******************************************************************************** 
subroutine get_degeneracies(eigenval, low, up, deg)
  real(wp), intent(in) :: eigenval(:)
  integer, intent(in)  :: low
  integer, intent(out) :: up
  integer, intent(out) :: deg
  !local
  real(wp), parameter  :: tol=1.0e-06_wp
  
  deg = 0
  up = low
  
  do
    if (up > size(eigenval)) exit
    if (abs(eigenval(low)-eigenval(up)) > tol) exit

    deg = deg + 1
    up = up + 1
  end do
end subroutine get_degeneracies


!******************************************************************************** 
!
!       Get groups of orbitals with same occupancies
!
!******************************************************************************** 
subroutine get_occ_groups(occ, low, up, ng)
  real(wp), intent(in)              :: occ(:)
  integer, allocatable, intent(out) :: low(:)
  integer, allocatable, intent(out) :: up(:)
  integer, intent(out)              :: ng
  !local variables
  real(wp), parameter               :: tol=1.0e-06_wp
  integer                           :: low_(size(occ)), up_(size(occ))
  integer                           :: p
  real(wp)                          :: old
  !
  ng = 1
  old = occ(1)
  low_(1) = 1
  p = 2
  !
  do while (p <= size(occ))
    if (abs(old - occ(p)) > tol) then
      up_(ng) = p - 1
      ng = ng + 1
      low_(ng) = p
      old = occ(p)
    end if
    p = p + 1
  end do 
  up_(ng) = p - 1 
  !
  if (allocated(low)) deallocate(low)
  if (allocated(up)) deallocate(up)
  allocate(low(ng), up(ng))
  low = low_(1:ng)
  up = up_(1:ng)
end subroutine get_occ_groups


!******************************************************************************** 
!
!       Estimates <S^2> for Hartrree-Fock wave functions
!
!           Remember RHF / ROHF wavefunctions are eigenfunctions of S^2,
!           while UHF wave function is only an approximation
!
!******************************************************************************** 
function estimate_s_square(phi, overlap, nel) result(s_square)
  real(wp), intent(in)  :: phi(:,:,:)
  real(wp), intent(in)  :: overlap(:,:,:)
  integer, intent(in)   :: nel(:)
  real(wp)              :: s_square
  !local variables
  real(wp), allocatable :: overlap_(:,:)
  integer               :: n, ispin, ispin1
  real(wp)              :: s

  n = size(phi, 1)
  ispin = size(phi, 3)
  ispin1 = size(overlap, 3)
  s = (nel(1)-nel(2))/2.0_wp
  s_square = s * (s+1.0_wp)

  if (ispin == 2) then
    allocate(overlap_(nel(1), nel(2)))
    overlap_ = matmul(transpose(phi(:,1:nel(1),1)), matmul(overlap(:,:,1), phi(:,1:nel(2),2)))
    s_square = s_square + nel(2)
    s_square = s_square - sum(overlap_**2)
  end if
end function estimate_s_square


!******************************************************************************** 
!
!       Calculates Mulliken occupation numbers
!
!******************************************************************************** 
subroutine  get_mulliken_occupations(ham_des, overlap, rdm, occ_dens)
  type(hamil_descriptor), intent(in) :: ham_des
  real(wp), intent(in)               :: overlap(:,:,:)
  real(wp), intent(in)               :: rdm(:,:,:)
  real(wp), allocatable, intent(out) :: occ_dens(:,:)
  !local variables
  integer                            :: i, n, spin, spin2, ispin
  !
  n = ham_des%n
  ispin = ham_des%ispin
  !
  if (.not. allocated(occ_dens)) allocate(occ_dens(n, ispin))
  do spin = 1, ispin
    spin2 = spin
    if (ham_des%ispin1 == 1) spin2 = 1
    do i = 1, n
      occ_dens(i,spin) = dot_product(overlap(:,i,spin2), rdm(:,i,spin)) * ham_des%rspin
    end do
  end do
end subroutine get_mulliken_occupations


!******************************************************************************** 
!
!       Get orbital character of the canonical orbitals
!
!******************************************************************************** 
subroutine get_projections(ham_des, phi, proj)
  type(hamil_descriptor), intent(in) :: ham_des
  real(wp), intent(in)               :: phi(:,:,:)
  real(wp), allocatable, intent(out) :: proj(:,:,:)
  !local variables
  integer                            :: i, n, spin, ispin
  !
  n = ham_des%n
  ispin = ham_des%ispin
  !
  if (.not. allocated(proj)) allocate(proj, mold=phi)
  !
  do spin = 1, ispin
    do i = 1, n
      proj(:,i,spin) = phi(:,i,spin) / norm2(phi(:,i,spin))
    end do
  end do
  proj = proj * proj * 100.0_wp
  !
  call print_projections(ham_des, proj)
end subroutine get_projections


!******************************************************************************** 
!
!       Print projections
!
!******************************************************************************** 
subroutine print_projections(ham_des, proj)
  type(hamil_descriptor), intent(in) :: ham_des
  real(wp), intent(in)               :: proj(:,:,:)
  !local
  integer                            :: block, block_size, nblocks, lb, ub
  integer                            :: spin
  type(FileHandle)                   :: fh
  character(len=*), parameter        :: fname = "qmcfort_proj"

  block_size = 20
  nblocks = ceiling(real(ham_des%n, wp) / real(block_size, wp))
  fh = FileHandle(fname)
  
  if (comm_world%mpirank == 0) then
    call fh%open(status="replace", action="write")
    write(fh%funit,*)
    write(fh%funit,*) "Projections of the orbitals on the Gaussian basis functions"
    do spin = 1, ham_des%ispin
      write(fh%funit,*) "Spin componenet ", spin
      do block = 1, nblocks
        lb = (block-1) * block_size + 1
        ub = min(block*block_size, ham_des%n)
        call print_projection_part(lb, proj(:,lb:ub,spin))
      end do
      write(fh%funit,*)
    end do
    call fh%close()
  end if
contains
  subroutine print_projection_part(lb, proj)
    integer, intent(in)  :: lb
    real(wp), intent(in) :: proj(:,:)
    !local variables
    integer              :: i, size_
    integer, allocatable :: index_(:)
    !
    size_ = size(proj, 2)
    allocate(index_(size_))
    do i = 1, size_
      index_(i) = lb + i - 1
    end do
    !
    write(fh%funit,100) "fun/orb", index_
    do i = 1, size(proj, 1)
      write(fh%funit,101) i, proj(i,:)
    end do
    write(fh%funit,*) 
    100 format(1x,a,t10,25i6)
    101 format(1x,i6,t10,25f6.2)   
  end subroutine print_projection_part
end subroutine print_projections


!******************************************************************************** 
!
!       Prints input of hf procedure
!
!******************************************************************************** 
subroutine print_hf_header(hf_ham)
  class(hf_hamil), intent(in) :: hf_ham
  
  if (comm_world%mpirank == 0) then
    write(io%screen%funit,*)       "SCF CALCULATION"
    write(io%screen%funit,*)       starshort
    write(io%screen%funit,*)       "Entering scf loop:"
    write(io%screen%funit,100)     " hf step  ", "     energy    ", "  energy diff  ", "    rdm norm   ", &
                                   "ldiis", "    diis_norm  ", " diis_size "
    
    write(io%qmcfort_log%funit,*)   "SCF CALCULATION"
    write(io%qmcfort_log%funit,*)   starshort
    write(io%qmcfort_log%funit,*)   "entering scf loop:"
    write(io%qmcfort_log%funit,100) " hf step  ", "     energy    ", "  energy diff  ", "   rdm norm   ", &
                                   "ldiis", "    diis_norm  ", " diis_size "
    
    write(io%qmcfort_out%funit,*)   "SCF CALCULATION"
    write(io%qmcfort_out%funit,*)   starshort
    write(io%qmcfort_out%funit,*)   "HF_descriptor:"
    write(io%qmcfort_out%funit,101) "maximal number of steps      maxstep       =", hf_ham%hf_des%hf_loop%maxstep
    write(io%qmcfort_out%funit,105) "convergence criterion energy en_tol        =", hf_ham%hf_des%hf_loop%en_tol
    write(io%qmcfort_out%funit,105) "convergence criterion rdm    rdm_tol       =", hf_ham%hf_des%hf_loop%rdm_tol
    write(io%qmcfort_out%funit,102) "mixing parameter for rdm     rdm_mix       =", hf_ham%hf_des%rdm_mix
    write(io%qmcfort_out%funit,102) "mixing parameter for uhf     uhf_mix       =", hf_ham%hf_des%uhf_mix
    write(io%qmcfort_out%funit,103) "scf initial guess            init_guess    =", hf_ham%hf_des%init_guess
    write(io%qmcfort_out%funit,*)   "entering scf loop:"
    write(io%qmcfort_out%funit,*)   starshort
  end if
  
  100 format (1x,t10,a,t20,a,t35,a,t50,a,t65,a,t70,a,t85,a) 
  101 format(1x,t5,a,t50,i3)
  102 format(1x,t5,a,t50,f10.6)
  103 format(1x,t5,a,t50,a)
  104 format(1x,t5,a,t50,l5)
  105 format(1x,t5,a,t50,es14.6)
end subroutine print_hf_header


!******************************************************************************** 
!
!       Prints output after hf procedure
!
!******************************************************************************** 
subroutine print_hf_footer(hf_ham)
  class(hf_hamil), intent(in) :: hf_ham
  
  if (comm_world%mpirank == 0) then
    write(io%screen%funit,*) starshort
    write(io%qmcfort_log%funit,*) starshort
   
    call print_hf_summary(hf_ham, io%screen)
    call print_hf_summary(hf_ham, io%qmcfort_log)
    call print_hf_summary(hf_ham, io%qmcfort_out)
  end if
  
  call hf_ham%hf_en%report(io%screen, "hf")
  call hf_ham%hf_en%report(io%qmcfort_log, "hf")
  call hf_ham%hf_en%report(io%qmcfort_out, "hf")
  call print_bands(hf_ham%ham%des, hf_ham%ham%eigenval, hf_ham%ham%occ)
  
  if (comm_world%mpirank == 0) then
    write(io%screen%funit,*)     starlong
    write(io%qmcfort_log%funit,*) starlong
    write(io%qmcfort_out%funit,*) starlong
  end if
end subroutine print_hf_footer


!******************************************************************************** 
!
!       Print info during the scf loop
!
!******************************************************************************** 
subroutine print_hf_info(hf_ham)
  class(hf_hamil), intent(in) :: hf_ham
  !local
  integer                     :: nbands, diis_size
  logical                     :: diis_ready
  real(wp)                    :: diis_err
  
  diis_ready = hf_ham%diis%is_ready()
  diis_size = hf_ham%diis%get_diis_size()
  diis_err = hf_ham%diis%rms_error()

  if (comm_world%mpirank == 0) then
    write(io%screen%funit,100)     "hf:", hf_ham%hf_des%hf_loop%step, hf_ham%hf_des%hf_loop%energy, hf_ham%hf_des%hf_loop%ediff, &
                 hf_ham%hf_des%hf_loop%rdmdiff, diis_ready, diis_err, diis_size
    write(io%qmcfort_log%funit,100) "hf:", hf_ham%hf_des%hf_loop%step, hf_ham%hf_des%hf_loop%energy, hf_ham%hf_des%hf_loop%ediff, &
                 hf_ham%hf_des%hf_loop%rdmdiff, diis_ready, diis_err, diis_size
    write(io%qmcfort_out%funit,101) "step number                  ", hf_ham%hf_des%hf_loop%step
    write(io%qmcfort_out%funit,102) "rdm norm                     ", hf_ham%hf_des%hf_loop%rdmdiff
    write(io%qmcfort_out%funit,103) "energy and energy difference ", hf_ham%hf_des%hf_loop%energy, hf_ham%hf_des%hf_loop%ediff
    write(io%qmcfort_out%funit,*)
    call hf_ham%hf_en%report(io%qmcfort_out)
    nbands = max(hf_ham%ham%des%nel(1), hf_ham%ham%des%nel(2))
    nbands = min(nbands+3, hf_ham%ham%des%n)
    call print_bands(hf_ham%ham%des, hf_ham%ham%eigenval, hf_ham%ham%occ, nbands)
    write(io%qmcfort_out%funit,*)   starshort
  end if
  !
  100 format (1x,a,t10,i6,t20,f14.8,t35,es14.6,t50,es14.6,t65,l5,t70,es14.6,t85,i6)
  101 format (1x,a,i3)
  102 format (1x,a,es14.6)
  103 format (1x,a,f14.8,es14.6)
end subroutine print_hf_info


!******************************************************************************** 
!
!       Print summary line for hf
!
!******************************************************************************** 
subroutine print_hf_summary(hf_ham, fh)
  type(hf_hamil), intent(in)    :: hf_ham
  type(FileHandle), intent(in)  :: fh
  
  write(fh%funit,101) "Eigenvalue of <S^{2}> operator  = ", hf_ham%hf_des%s_square
  write(fh%funit,101) "Spin multiplicity 2S + 1 = ", hf_ham%hf_des%s_mult
  write(fh%funit,100) "hf_summary: ", hf_ham%hf_des%hf_loop%step, hf_ham%hf_des%hf_loop%converged, &
                                 hf_ham%hf_en%total_energy(), hf_ham%hf_des%hf_loop%ediff, hf_ham%hf_des%hf_loop%rdmdiff
  write(fh%funit,*)
  100 format (1x,a,t15,i6,t30,l5,t45,f14.8,t60,es14.6,t75,es14.6)
  101 format (1x,a,f8.4)
end subroutine print_hf_summary


!******************************************************************************** 
!
!       Print energy bands info
!
!******************************************************************************** 
subroutine print_bands(ham_des, eigenval, occ, nbands)
  type(hamil_descriptor), intent(in) :: ham_des
  real(wp), allocatable,  intent(in) :: eigenval(:,:)
  real(wp), intent(in)               :: occ(:,:)
  integer, optional, intent(in)      :: nbands
  !local variables                   
  integer                            :: p, spin, nbands_
  real(wp)                           :: occup
  !
  if (comm_world%mpirank==0 .and. allocated(eigenval)) then
    do spin = 1, ham_des%ispin
      write(io%qmcfort_out%funit,*)
      write(io%qmcfort_out%funit,100) "eigenvalues spin component: ", spin
      write(io%qmcfort_out%funit,101) "No.", "eigenval", "occ"
      if (present(nbands)) then
        nbands_ = nbands
      else
        nbands_ = ham_des%n
      end if
      !
      do p = 1, nbands_
        if (ham_des%ispin==1) then
          occup = occ(p,1) + occ(p,2)
        else
          occup = occ(p,spin)
        end if
        write(io%qmcfort_out%funit, 102) p, eigenval(p,spin), occup
      end do
    end do
    100 format (1x,a,i3)
    101 format (1x,t10,a,t15,a,t30,a)
    102 format (1x,t10,i3,t15,f12.6,t30,f8.6)
  end if
end subroutine print_bands

end module hfproc