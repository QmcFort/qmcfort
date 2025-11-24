! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module hamilton_type

#include "../preproc.inc"
use constants
use file_handle, only: FileHandle
use array_file_io, only: ArrayFileIO
use array_io_factories, only: array_io_factory
use mpi, only: comm_world
use profiling
use method_base
use qmcfort_io
use qmcfort_pos
use qmcfort_in
use energy_types
use hamilton_vars
use hamilton_layout
use hamilton_operations
use hamilton_io
use gauss_base
use gauss_integrals
use cholesky

implicit none

private
public Hamilton, read_files, adjust_wavespin

type Hamilton
  type(hamil_descriptor)             :: des
  type(hamil_io)                     :: h_io
  type(gauss_basis_set)              :: basis_set
  
  real(wp)                           :: enuc = 0.0_wp
  real(wp)                           :: h0 = 0.0_wp

  real(wp), allocatable              :: overlap(:,:,:)
  real(wp), allocatable              :: ext_pot(:,:,:)
  real(wp), allocatable              :: kinetic(:,:,:)
  real(wp), allocatable              :: h1(:,:,:)
  real(wp), allocatable              :: h1_orig(:,:,:)
  real(wp), allocatable              :: p1(:,:,:)

  real(wp), allocatable              :: h2(:,:,:,:,:)
  real(wp), allocatable              :: h2_gres(:,:,:)
  real(wp), pointer                  :: h2_ptr(:,:,:) => null()

  real(wp), allocatable              :: hself(:,:,:)
  real(wp), allocatable              :: hself_gres(:,:,:,:)

  real(wp), allocatable              :: hpot(:,:,:)
  real(wp), allocatable              :: hpot_gres(:,:,:,:)

  real(wp), allocatable              :: xpot(:,:,:)
  real(wp), allocatable              :: xpot_gres(:,:,:,:)

  real(wp), allocatable              :: phi(:,:,:)
  real(wp), allocatable              :: eigenval(:,:)
  real(wp), allocatable              :: occ(:,:)

  procedure(Ipotential), pointer     :: hartree_potential => null()
  procedure(Ipotential), pointer     :: exchange_potential => null()
  procedure(Ipotential), pointer     :: self_potential => null()

  procedure(Ipotential_rdm), pointer :: hartree_potential_rdm => null()
  procedure(Ipotential_rdm), pointer :: exchange_potential_rdm => null()
contains
  procedure          :: set => set_hamiltonian
  procedure          :: reader => hamilton_reader
  procedure          :: write_ => write_ham_files
  procedure          :: sizeof => sizeof_hamilton
  procedure          :: reallocate => reallocate_spin_arrays

  procedure          :: meet_method_req
  procedure          :: check_basis, check_integral_mode

  procedure          :: associate_gres_pointers
  procedure          :: associate_h2_ptr
  procedure          :: get_h2_gres_ob

  procedure          :: decompose => decompose_hamiltonian
  procedure          :: assemble => assemble_hamiltonian

  procedure          :: transform_basis
  procedure          :: init_orbitals, init_overlap
  procedure          :: add_self_potential

  generic            :: subtract_mean_field => subtract_mean_field_r, subtract_mean_field_c
  procedure, private :: subtract_mean_field_r, subtract_mean_field_c

  procedure          :: get_cas_ham

  procedure, private :: hartree_potential_eri, hartree_potential_gres
  procedure, private :: exchange_potential_eri, exchange_potential_gres
  procedure, private :: self_potential_eri, self_potential_gres

  procedure, private :: hartree_potential_rdm_eri, hartree_potential_rdm_gres
  procedure, private :: exchange_potential_rdm_eri, exchange_potential_rdm_gres
end type Hamilton

abstract interface
  subroutine Ipotential(self, pot, first_state, last_state)
    import Hamilton, wp
    implicit none
    class(Hamilton), intent(inout)       :: self
    real(wp), allocatable, intent(inout) :: pot(:,:,:)
    integer, optional, intent(in)        :: first_state
    integer, optional, intent(in)        :: last_state
  end subroutine Ipotential

  subroutine Ipotential_rdm(self, rdm, pot, first_state, last_state)
    import Hamilton, wp
    implicit none
    class(Hamilton), intent(inout)       :: self
    real(wp), intent(in)                 :: rdm(:,:,:)
    real(wp), allocatable, intent(inout) :: pot(:,:,:)
    integer, optional, intent(in)        :: first_state
    integer, optional, intent(in)        :: last_state
  end subroutine Ipotential_rdm
end interface

contains

!******************************************************************************** 
!
! Set Hamiltonian matrix elements
!
!******************************************************************************** 
subroutine set_hamiltonian(ham, struc)
  class(Hamilton), intent(inout)  :: ham
  type(structure), intent(in)     :: struc
  !local variables
  integer                         :: rank
  logical                         :: lverb
  type(FileHandle)                :: fh
  class(ArrayFileIO), allocatable :: array_io

  fh = FileHandle("basis_set")
  
  call ham%des%init
  ham%enuc = struc%calc_enuc()
  call ham%reader
  call ham%h_io%init_files

  if (ham%h_io%if_files_exist()) then
    call ham%h_io%read_ham_info(ham%des)
    call array_io_factory(ham%h_io%file_format, array_io)
    call ham%des%set_electrons_nel(ham%des%nel(1), ham%des%nel(2))
    if (ham%h_io%h0%exists()) call array_io%load(ham%h_io%h0%fname, ham%h0)
    if (ham%h_io%h1%exists()) call array_io%load(ham%h_io%h1%fname, ham%h1)
    if (ham%h_io%overlap%exists()) call array_io%load(ham%h_io%overlap%fname, ham%overlap)
    call ham%init_overlap(array_io)
    call ham%h_io%read_h2(ham%des, ham%h2, ham%h2_gres, array_io)
    allocate(ham%h1_orig, source=ham%h1)
    if (ham%des%compress) call compress_h2(ham%des, ham%h2, ham%h2_gres)
  else if (fh%exists() .or. len_trim(ham%des%basis_set)>1) then
    if ((.not. fh%exists()) .and. len_trim(ham%des%basis_set)>1) then
      call generate_basis_set_file(ham%des%basis_set, struc%get_species_str())
    end if
    call ham%basis_set%read_mpi(ham%des, struc)
    ham%des%basis = "ao"
    call ham%des%set_electrons_struc(struc)
    ham%des%ispin1 = 1
    ham%des%ispin2 = 1
    call set_gauss_integrals(ham, struc)
  else
    if (comm_world%mpirank == 0) write(*,*) "Neither Hamiltonian files nor BASIS_SET file can not be found"
    if (comm_world%mpirank == 0) write(*,*) "There is no enough information to start the program - Program stops now"
    call exit
  end if
  
  call ham%associate_h2_ptr()
  call ham%associate_gres_pointers()

  call array_io_factory(ham%h_io%file_format, array_io)
  if (ham%h_io%eigenval%exists()) call array_io%load(ham%h_io%eigenval%fname, ham%eigenval)
  call ham%init_orbitals(array_io)
  ham%des%h0 = ham%h0
end subroutine set_hamiltonian


!******************************************************************************** 
!
! Reader for hamil_descriptor and hamil_io
!
!******************************************************************************** 
subroutine hamilton_reader(self)
  use qmcfort_in, only: add_input
  class(Hamilton), intent(inout) :: self
  
  call add_input("basis_set",        self%des%basis_set)
  call add_input("ispin",            self%des%ispin)
  call add_input("spin",             self%des%spin)
  call add_input("nfrozen",          self%des%nfrozen)
  call add_input("frozen",           self%des%nfrozen)
  call add_input("nelect",           self%des%nelect)
  call add_input("nelect_up",        self%des%nel(1))
  call add_input("nelect_alpha",     self%des%nel(1))
  call add_input("nelect_down",      self%des%nel(2))
  call add_input("nelect_beta",      self%des%nel(2))
  call add_input("compress",         self%des%compress)
  call add_input("compress_h2",      self%des%compress)
  call add_input("compress_mode",    self%des%compress_mode)
  call add_input("compress_h2_mode", self%des%compress_mode)
  call add_input("integral_mode",    self%des%integral_mode)
  call add_input("cholesky_tol",     self%des%chol_tol)
  call add_input("chol_tol",         self%des%chol_tol)
  call add_input("compress_tol",     self%des%chol_tol)
  call add_input("chol_eps",         self%des%chol_tol)
  call add_input("write_files",      self%h_io%iwrite)
  call add_input("iwrite",           self%h_io%iwrite)
  call add_input("file_format",      self%h_io%file_format)
  call add_input("fmt",              self%h_io%file_format)
end subroutine hamilton_reader


!******************************************************************************** 
!
! Storage size of the Hamilton type
!
!******************************************************************************** 
function sizeof_hamilton(self) result(mem)
  class(Hamilton), intent(in) :: self
  real(wp)                    :: mem

  mem = sizeof(self) + sizeof(self%overlap) + sizeof(self%ext_pot) + sizeof(self%kinetic) + sizeof(self%hself) &
      + sizeof(self%h1) + sizeof(self%h1_orig) + sizeof(self%h2) + sizeof(self%h2_gres) + sizeof(self%hpot) &
      + sizeof(self%xpot) + sizeof(self%eigenval) &
      + sizeof(self%phi) + sizeof(self%occ) + sizeof(self%p1)
end function sizeof_hamilton


!******************************************************************************** 
!
! Reallocate Hamiltoanian arrays
!        
!    Copies allocated hamiltonian arrays and changes spin dimensions
!    If basis=mean_field and ispin=2 spin dependent matrix elements are needed.  
!
!            ispin1=1 --> ispin1=2
!            ispin2=1 --> ispin2=3
!
!******************************************************************************** 
subroutine reallocate_spin_arrays(self)
  class(Hamilton), intent(inout) :: self
  !local variables
  type(Hamilton)                 :: temp
  integer                        :: spin
  
  if (self%des%ispin1==2 .and. self%des%ispin2==3) return

  if (trim(self%des%basis)=="mo" .and. self%des%ispin==2) then
    self%des%ispin1 = 2
    self%des%ispin2 = 3
    
    temp%des = self%des

    if (allocated(self%overlap)) then
      allocate(temp%overlap(temp%des%n, temp%des%n, temp%des%ispin1))
      do spin = 1, temp%des%ispin1
        temp%overlap(:,:,spin) = self%overlap(:,:,1)
      end do
      call move_alloc(temp%overlap, self%overlap)
    end if

    if (allocated(self%kinetic)) then
      allocate(temp%kinetic(temp%des%n, temp%des%n, temp%des%ispin1))
      do spin = 1, temp%des%ispin1
        temp%kinetic(:,:,spin) = self%kinetic(:,:,1)
      end do
      call move_alloc(temp%kinetic, self%kinetic)
    end if

    if (allocated(self%ext_pot)) then
      allocate(temp%ext_pot(temp%des%n, temp%des%n, temp%des%ispin1))
      do spin = 1, temp%des%ispin1
        temp%ext_pot(:,:,spin) = self%ext_pot(:,:,1)
      end do
      call move_alloc(temp%ext_pot, self%ext_pot)
    end if

    if (allocated(self%h1)) then
      allocate(temp%h1(temp%des%n, temp%des%n, temp%des%ispin1))
      do spin = 1, temp%des%ispin1
        temp%h1(:,:,spin) = self%h1(:,:,1)
      end do
      call move_alloc(temp%h1, self%h1)
    end if

    if (allocated(self%h1_orig)) then
      allocate(temp%h1_orig(temp%des%n, temp%des%n, temp%des%ispin1))
      do spin = 1, temp%des%ispin1
        temp%h1_orig(:,:,spin) = self%h1_orig(:,:,1)
      end do
      call move_alloc(temp%h1_orig, self%h1_orig)
    end if

    if (allocated(self%h2)) then
      allocate(temp%h2(temp%des%n, temp%des%n, temp%des%n, temp%des%n, temp%des%ispin2))
      do spin = 1, temp%des%ispin2
        temp%h2(:,:,:,:,spin) = self%h2(:,:,:,:,1)
      end do
      call move_alloc(temp%h2, self%h2)
      call self%associate_h2_ptr()
    end if

    if (allocated(self%h2_gres)) then
      allocate(temp%h2_gres(temp%des%n*temp%des%n, temp%des%ng, temp%des%ispin1))
      do spin = 1, temp%des%ispin1
        temp%h2_gres(:,:,spin) = self%h2_gres(:,:,1)
      end do
      call move_alloc(temp%h2_gres, self%h2_gres)
    end if
  end if
end subroutine reallocate_spin_arrays


!********************************************************************************
!
! Write Hamiltonian files according to h_io%iwrite
!
!    iwrite = 0 :
!        do not write files
!
!    iwrite = 1 :
!        orbitals
!        eigenval
!
!    iwrite = 2 :
!        orbitals 
!        eigenval 
!        overlap 
!        ham1
!        ham2
!
!********************************************************************************
subroutine write_ham_files(self)
  class(Hamilton), intent(inout)  :: self
  !local varaibles
  logical                         :: lverb
  character(len=*), parameter     :: proc_name = "write_ham_files"
  class(ArrayFileIO), allocatable :: array_io

  if (use_profiling) call start_profiling(proc_name)

  if (comm_world%mpirank == 0) then
    if (self%h_io%iwrite == 1) then
      call write_ham_io_header()
      call array_io_factory(self%h_io%file_format, array_io)
      call self%h_io%write_ham_info(self%des)
      call array_io%save(self%h_io%orbitals%fname, self%phi)
      if (allocated(self%eigenval)) call array_io%save(self%h_io%eigenval%fname, self%eigenval)
      call write_ham_io_footer
    else if (self%h_io%iwrite == 2) then
      call write_ham_io_header()
      call array_io_factory(self%h_io%file_format, array_io)
      call self%h_io%write_ham_info(self%des)
      call array_io%save(self%h_io%orbitals%fname, self%phi)
      if (allocated(self%eigenval)) call array_io%save(self%h_io%eigenval%fname, self%eigenval)
      call array_io%save(self%h_io%overlap%fname, self%overlap)
      call array_io%save(self%h_io%h0%fname, self%h0)
      call array_io%save(self%h_io%h1%fname, self%h1)
      call self%h_io%write_h2(self%des, self%h2, self%h2_gres, array_io)
      call write_ham_io_footer
    end if
  end if

  if (use_profiling) call end_profiling(proc_name)
contains
  subroutine write_ham_io_header
    if (comm_world%mpirank == 0) then
      write(*,*) "WRITING HAMILTONIAN FILES"
      write(*,*) starshort
    end if
  end subroutine write_ham_io_header
  
  subroutine write_ham_io_footer
    if (comm_world%mpirank == 0) then
      write(*,*) starlong
    end if
  end subroutine write_ham_io_footer
end subroutine write_ham_files


!******************************************************************************** 
!
! Associate h2_ptr 
!
! Original shapes:
!    h2(N,N,N,N,ispin2)
!    h2_ptr(N^2,N^2,ispin2)
!
!******************************************************************************** 
subroutine associate_h2_ptr(self)
  class(Hamilton), target, intent(inout) :: self
  !local variables
  integer                                :: nn, ispin

  if (allocated(self%h2)) then
    nn = size(self%h2,1) * size(self%h2,1)
    ispin = size(self%h2,5)
    self%h2_ptr(1:nn,1:nn,1:ispin) => self%h2
  else
    self%h2_ptr => null()
  end if
end subroutine associate_h2_ptr


!******************************************************************************** 
!
! Obtain g-resolved two-body Hamiltonian in "orbital basis" (N,N,Ng,ispin1)
!
! Note:
!    Original layout of the h2_gres is "pair-orbital basis" (N^2,Ng,ispin1)
!
!******************************************************************************** 
subroutine get_h2_gres_ob(self, h2_gres_ob)
  class(Hamilton), target, intent(in) :: self
  real(wp), pointer, intent(out)      :: h2_gres_ob(:,:,:,:)
  !local variables
  integer                             :: n, ng, ispin
  character(len=*), parameter         :: proc_name = "get_h2_gres_ob"

  if (use_profiling) call start_profiling(proc_name)

  if (allocated(self%h2_gres)) then
    n = nint(sqrt(real(size(self%h2_gres, 1))))
    ng = size(self%h2_gres,2)
    ispin = size(self%h2_gres,3)
    h2_gres_ob(1:n,1:n,1:ng,1:ispin) => self%h2_gres
  else
    h2_gres_ob => null()
  end if

  if (use_profiling) call end_profiling(proc_name)
end subroutine get_h2_gres_ob


!******************************************************************************** 
!
! Initialize orbitals
!
!    Read from file or
!    Set to unit matrix or
!    Set to zeros
!
!******************************************************************************** 
subroutine init_orbitals(self, array_io)
  use standalone, only: eye
  class(Hamilton), intent(inout) :: self
  class(ArrayFileIO), intent(in) :: array_io
  !local variables
  integer                        :: spin
  logical                        :: lverb
  real(wp), allocatable          :: eigenval(:,:)

  if (self%des%nbtot == unset) self%des%nbtot = self%des%n

  if (self%h_io%orbitals%exists()) then
    call array_io%load(self%h_io%orbitals%fname, self%phi)
    if (comm_world%mpirank == 0) write(*,*) "orbitals read from file"
    call adjust_wavespin(self%des, self%phi, self%eigenval)
  else
    if (.not. allocated(self%phi)) allocate(self%phi(self%des%n, self%des%nbtot, self%des%ispin))

    if (trim(self%des%basis) == "mo") then
      do spin = 1, size(self%phi, 3)
        call eye(self%phi(:,:,spin))
      end do
      if (comm_world%mpirank == 0) write(*,*) "orbitals set to unit matrix"
    else
      self%phi = 0.0_wp
      if (comm_world%mpirank == 0) write(*,*) "orbitals not found - set to zero"
    end if

    !adjust eigenvalues array with the orbitals
    allocate(eigenval(self%des%n, self%des%ispin))
    eigenval = 0.0_wp
    call move_alloc(eigenval, self%eigenval)
  end if
end subroutine init_orbitals


!******************************************************************************** 
!
! Initialize overlap matrix
!
!******************************************************************************** 
subroutine init_overlap(self, array_io)
  use standalone, only: eye
  class(Hamilton), intent(inout) :: self
  class(ArrayFileIO), intent(in) :: array_io
  !local variables
  integer                        :: spin
  logical                        :: lverb
  
  if (self%h_io%overlap%exists()) then
    call array_io%load(self%h_io%overlap%fname, self%overlap)
    if (comm_world%mpirank == 0) write(*,*) "overlap read from file"
  else
    if (.not. allocated(self%overlap)) allocate(self%overlap(self%des%n, self%des%n, self%des%ispin1))
    do spin = 1, size(self%overlap, 3)
      call eye(self%overlap(:,:,spin))
    end do
    if (comm_world%mpirank == 0) write(*,*) "overlap set to unit matrix"
  end if
end subroutine init_overlap


!******************************************************************************** 
!
!       If ham%des%ispin doesn't match with size(phi, 3), adjust spin variable
!
!******************************************************************************** 
subroutine adjust_wavespin(ham_des, phi, eigenval)
  type(hamil_descriptor), intent(in)             :: ham_des
  real(wp), allocatable, intent(inout)           :: phi(:,:,:)
  real(wp), allocatable, optional, intent(inout) :: eigenval(:,:)
  !local variables
  real(wp), allocatable                          :: phi_(:,:,:), eigenval_(:,:)
  
  if (ham_des%ispin == size(phi, 3)) return
  
  if (ham_des%ispin == 1) then
    allocate(phi_(size(phi,1), size(phi,2), ham_des%ispin))
    phi_(:,:,1) = phi(:,:,1)
    call move_alloc(phi_, phi)
    if (present(eigenval)) then
      if (allocated(eigenval) .and. size(eigenval,2)==2) then
        allocate(eigenval_(size(eigenval, 1), ham_des%ispin))
        eigenval_(:,1) = eigenval(:,1)
        call move_alloc(eigenval_, eigenval)
      end if
    end if
    if (comm_world%mpirank == 0) write(*,"(1x,a,i6)") "orbitals adjusted to ispin =", ham_des%ispin
  else
    allocate(phi_(size(phi,1), size(phi,2), ham_des%ispin))
    phi_(:,:,1) = phi(:,:,1)
    phi_(:,:,2) = phi(:,:,1)
    call move_alloc(phi_, phi)
    if (present(eigenval)) then
      if (allocated(eigenval) .and. size(eigenval,2)==1) then
        allocate(eigenval_(size(eigenval, 1), ham_des%ispin))
        eigenval_(:,1) = eigenval(:,1)
        eigenval_(:,2) = eigenval(:,1)
        call move_alloc(eigenval_, eigenval)
      end if
    end if
    if (comm_world%mpirank == 0) write(*,"(1x,a,i6)") "orbitals adjusted to ispin =", ham_des%ispin
  end if
end subroutine adjust_wavespin


!******************************************************************************** 
!
! Main routine for Hamiltonian decomposition
!
! Transform h2 to h2_gres
!
!******************************************************************************** 
subroutine decompose_hamiltonian(ham)
  class(Hamilton), intent(inout) :: ham
  
  if (.not. allocated(ham%h2) .or. allocated(ham%h2_gres)) return
  call compress_h2(ham%des, ham%h2, ham%h2_gres)
end subroutine decompose_hamiltonian


!******************************************************************************** 
!
! Transform h2_gres to h2
!
!******************************************************************************** 
subroutine assemble_hamiltonian(ham, dealloc)
  use lapack, only: gemm
  class(Hamilton), intent(inout)   :: ham
  logical, optional, intent(in) :: dealloc

  if (.not. allocated(ham%h2_gres) .or. allocated(ham%h2)) return
  call get_h2_full(ham%h2_gres, ham%h2, ham%des%n, ham%des%n, ham%des%ispin2, dealloc)
  ham%des%integral_mode = "eri"
end subroutine assemble_hamiltonian


!******************************************************************************** 
!
!       Main routine for the basis transformation
!
!******************************************************************************** 
subroutine transform_basis(ham, basis)
  use orth
  class(Hamilton), intent(inout)         :: ham
  character(len=*), optional, intent(in) :: basis
  !local
  type(basis_trafo_descriptor)           :: trafo_des
  real(wp), allocatable                  :: trafo(:,:,:,:)
 
  call trafo_des%init(ham%des, basis)
  if (.not. trafo_des%basis_changed) return

  call trafo_des%print_header()

  call calculate_trafo_matrix(trafo_des, ham%overlap, ham%phi, trafo)

  call ham%reallocate()
  call transform_ham1(trafo, ham%overlap)
  call transform_ham1(trafo, ham%h1)
  call transform_ham1(trafo, ham%h1_orig)
  call transform_ham2(trafo, ham%h2, ham%h2_gres)
  call transform_orbitals(trafo, ham%phi)

  call trafo_des%print_footer()
end subroutine transform_basis


!********************************************************************************
!
!       Calculates all molecular integrals over Gaussian functions stored in 
!       basis_set file and collect one- and two-electron integrals
!
!       output:         
!               (i,s|h|j,s)     = h1(i,j,s)            and
!               (is,js|kt,lt)   = h2(ijkl,st) == <is,kt|js,lt>  or
!               L_g(i,j)        = h2_gres(i,j,g)
!
!********************************************************************************
subroutine set_gauss_integrals(ham, struc)
  type(Hamilton), intent(inout) :: ham
  type(structure), intent(in)   :: struc
  !local variables
  character(len=*), parameter   :: proc_name = "set_gauss_integrals" 

  if (use_profiling) call start_profiling(proc_name)

  call print_gauss_integrals_header
  
  call collect_int_1e(ham%basis_set, struc, ham%overlap, ham%kinetic, ham%ext_pot, ham%des%ispin1)

  if (.not. allocated(ham%h1)) allocate(ham%h1, mold=ham%kinetic)
  ham%h1 = ham%kinetic + ham%ext_pot
  ham%h1_orig = ham%h1
  
  if (trim(ham%des%integral_mode) == "eri") then 
    call collect_int_2e(ham%basis_set, ham%h2, ham%des%ispin2)
  else if (trim(ham%des%integral_mode) == "cholesky") then
    call compress_h2_chol(ham%des, ham%basis_set, ham%h2_gres)
  end if
  
  call print_gauss_integrals_footer

  if (use_profiling) call end_profiling(proc_name)
end subroutine set_gauss_integrals


!******************************************************************************** 
!
! Check if the basis of the Hamiltonian matrix elements and integral_mode are 
! compatible with the QmcFort method requirements
!
! If not, appropriate transformations will be performed
!
!******************************************************************************** 
subroutine meet_method_req(self, method)
  class(Hamilton), intent(inout)  :: self
  type(QmcFortMethod), intent(in) :: method
  
  call self%check_basis(method)
  call self%check_integral_mode(method)
end subroutine meet_method_req


!******************************************************************************** 
!
! Check if the actual Hamiltonian basis meets the method requirements
!
! If not, Hamiltonian matrix elements will be transformed to the preferred basis
!
!******************************************************************************** 
subroutine check_basis(self, method)
  use string, only: isin, extract_word
  class(Hamilton), intent(inout)  :: self
  type(QmcFortMethod), intent(in) :: method
  !local
  character(len=:), allocatable   :: basis
  
  basis = method%get_basis()

  !if any basis allowed, return immediately
  if (basis == "any") return

  !if the actual basis is one of the allowed basis, return again
  if (isin(self%des%basis, basis)) return

  !otherwise, change the basis
  basis = extract_word(basis)
  call self%transform_basis(basis)
end subroutine check_basis


!******************************************************************************** 
!
! Check if the two-body Hamiltonian meets the method requirements
!
! If not, integrals will be transformed to the preferred form
!
!******************************************************************************** 
subroutine check_integral_mode(self, method)
  use string, only: isin, extract_word
  class(Hamilton), intent(inout)  :: self
  type(QmcFortMethod), intent(in) :: method
  !local
  character(len=:), allocatable   :: integral_mode
  
  integral_mode = method%get_integral_mode()

  !if any integral_mode allowed, return immediately
  if (trim(integral_mode) == "any") return

  !if the actual integral_mode is one of the allowed modes, return again
  if (isin(self%des%integral_mode, integral_mode)) return
  
  integral_mode = extract_word(integral_mode)
  
  select case (trim(integral_mode))
    case ("eri")
      call self%assemble()
    case ("cholesky")
      call self%decompose()
    case default
  end select
end subroutine check_integral_mode


!******************************************************************************** 
!
! Add self-energy Hamiltonian to the one-body Hamiltonian 
!
!******************************************************************************** 
subroutine add_self_potential(self)
  class(Hamilton), intent(inout) :: self

  if (.not. allocated(self%hself)) call self%self_potential(self%hself)
  self%h1 = self%h1 + self%hself
end subroutine add_self_potential


!******************************************************************************** 
!
! Header for gauss integrals calculation
!
!******************************************************************************** 
subroutine print_gauss_integrals_header()
  if (comm_world%mpirank == 0) then
    write(io%screen%funit,*)     "GAUSSIAN INTEGRALS CALCULATION"
    write(io%screen%funit,*)     starshort
    
    write(io%qmcfort_log%funit,*) "GAUSSIAN INTEGRALS CALCULATION"
    write(io%qmcfort_log%funit,*) starshort
    
    write(io%qmcfort_out%funit,*) "GAUSSIAN INTEGRALS CALCULATION"
    write(io%qmcfort_out%funit,*) starshort
  end if
end subroutine print_gauss_integrals_header


!******************************************************************************** 
!
! Footer for gauss integrals calculation
!
!******************************************************************************** 
subroutine print_gauss_integrals_footer()
  if (comm_world%mpirank == 0) then
    write(io%screen%funit,*)     starlong
    write(io%qmcfort_log%funit,*) starlong
    write(io%qmcfort_out%funit,*) starlong
  end if
end subroutine print_gauss_integrals_footer


!********************************************************************************
!       Read all files to initialize the program:
!               qmcfort_in
!               qmcfort_pos
!               BASIS_SET
!               ham1
!               ham2
!               orbitals
!               overlap
!               eigenval
!
!debug: shouldn't really be in hamilton.f
!********************************************************************************
subroutine read_files(struc, info, ham, sym)
  type(Structure), intent(in)    :: struc
  type(controll), intent(inout)  :: info
  type(Hamilton), intent(inout)  :: ham
  type(symmetry), intent(inout)  :: sym
  !local variables
  character(len=*), parameter    :: proc_name = "read_files"
  
  if (use_profiling) call start_profiling(proc_name)

  call print_read_files_header

  call read_input(info, sym)
  
  call ham%set(struc)
  
  call write_info(ham%des)
!debug:
call copy_to_info(info, ham%des)
  call print_read_files_footer

  if (use_profiling) call end_profiling(proc_name)
contains
  subroutine print_read_files_header
    if (comm_world%mpirank == 0) then
      write(*,*) "READING INPUT FILES"
      write(*,*) starshort
    end if
  end subroutine print_read_files_header

  subroutine print_read_files_footer
    if (comm_world%mpirank == 0) then 
      write(io%screen%funit,*) starlong
      write(io%qmcfort_out%funit,*) starlong
    end if
  end subroutine print_read_files_footer
  
  subroutine write_info(ham_des)
    type(hamil_descriptor), intent(inout) :: ham_des
    
    if (comm_world%mpirank == 0) then
      write(io%qmcfort_out%funit,*)   "Information about system:"
      write(io%qmcfort_out%funit,100) "Number of basis functions", "n          =", ham_des%n
      write(io%qmcfort_out%funit,100) "Number of electrons      ", "nelect     =", ham_des%nelect
      write(io%qmcfort_out%funit,100) "Number of up electrons   ", "nelect_a   =", ham_des%nel(1)
      write(io%qmcfort_out%funit,100) "Number of beta electrons ", "nelect_b   =", ham_des%nel(2)
      write(io%qmcfort_out%funit,100) "Number of occupied states", "nocc       =", ham_des%nocc
      write(io%qmcfort_out%funit,100) "Sz spin component        ", "spin       =", ham_des%spin
      write(io%qmcfort_out%funit,100) "Spin dimension           ", "ispin      =", ham_des%ispin
      write(io%qmcfort_out%funit,100) "Spin dimension 1         ", "ispin1     =", ham_des%ispin1
      write(io%qmcfort_out%funit,100) "Spin dimension 2         ", "ispin2     =", ham_des%ispin2
      write(io%qmcfort_out%funit,100) "Spin dimension FOCK      ", "ispin_fock =", ham_des%ispin_fock
      write(io%qmcfort_out%funit,100) "Spin multiplicity        ", "rspin      =", ham_des%rspin
    end if
    100 format(1x,t5,a,t35,a,t55,i6)
  end subroutine write_info
 
  subroutine copy_to_info(info, ham_des)
    type(controll), intent(inout)      :: info
    type(hamil_descriptor), intent(in) :: ham_des
    
    info%n = ham_des%n
    info%nelect = ham_des%nelect
    info%nocc = ham_des%nocc
    info%nel = ham_des%nel
    info%spin = ham_des%spin
    info%rspin = ham_des%rspin
    info%ispin = ham_des%ispin
    info%ispin1 = ham_des%ispin1
    info%ispin2 = ham_des%ispin2
    info%nelect_a = ham_des%nel(1)
    info%nelect_b = ham_des%nel(2)
  end subroutine copy_to_info
end subroutine read_files


!********************************************************************************
!
! This subroutine reads all flags from qmcfort_in file
!
!******************************************************************************** 
subroutine read_input(info, sym)
  type(controll), intent(inout)  :: info
  type(symmetry), intent(inout)  :: sym
  !local variables
  character(len=charlen)         :: char_mask
  
  call add_input("lci_vec", info%lci_vec)
  call add_input("ci_mode", info%ci_mode)
  call add_input("ci_diag_mode", info%ci_diag_mode)
  call add_input("ci_maxiter", info%ci_maxiter)
  call add_input("ci_eps", info%ci_eps)
  call add_input("sym", sym%sym)
  call add_input("eq_steps", info%eq_steps)
  call add_input("warm_up", info%warm_up)
  call add_input("num_mc_steps", info%num_mc_steps)
  call add_input("mc_tstep", info%mc_tstep)
  call add_input("print_step", info%print_step)
  call add_input("a_per", info%a_per)
  call add_input("damp", info%damp)
  call add_input("num_eigs", info%num_eigs)
  call add_input("write_ci", info%lwrite_ci)
  call add_input("randomize", info%lrandomize)
  call add_input("p_s", P_s)
  call add_input("w_start", info%w_start)
  call add_input("units", info%units)
  call add_input("lmask", info%lmask)
  call add_input("orb_mask", char_mask)
  call add_input("f_c", info%f_c)
  call add_input("fc", info%f_c)
  call add_input("spin_res", info%lspin_res)
  call add_input("ldmat", info%ldmat)
  call add_input("ldmet", info%ldmet)
  
  if (comm_world%mpirank == 0) write(*,*) "File qmcfort_in read successfully"
end subroutine read_input


!********************************************************************************
!
! Associate procedure pointers that depend on the integral_mode
!
!********************************************************************************
subroutine associate_gres_pointers(self)
  class(Hamilton), intent(inout) :: self

  select case (self%des%integral_mode)
    case ("eri")
      self%hartree_potential => hartree_potential_eri
      self%exchange_potential => exchange_potential_eri
      self%self_potential => self_potential_eri
      self%hartree_potential_rdm => hartree_potential_rdm_eri
      self%exchange_potential_rdm => exchange_potential_rdm_eri
    case ("cholesky")
      self%hartree_potential => hartree_potential_gres
      self%exchange_potential => exchange_potential_gres
      self%self_potential => self_potential_gres
      self%hartree_potential_rdm => hartree_potential_rdm_gres
      self%exchange_potential_rdm => exchange_potential_rdm_gres
    case default
      if (comm_world%mpirank == 0) write(*,*) "Unrecognized integral_mode value: program stops now!"
      stop
  end select
end subroutine associate_gres_pointers


!********************************************************************************
!
! Calculate Hartree potential from ERIs in the canonical basis
!
! Information about electrons provided in self%des
!
! Actual implementation in hamilton_operations.f
!
!********************************************************************************
subroutine hartree_potential_eri(self, pot, first_state, last_state)
  class(Hamilton), intent(inout)       :: self
  real(wp), allocatable, intent(inout) :: pot(:,:,:)
  integer, optional, intent(in)        :: first_state
  integer, optional, intent(in)        :: last_state

  call calculate_hartree_potential(self%h2, [1,1], self%des%nel, self%des%ispin_fock, self%des%rspin, pot, first_state, last_state)
end subroutine hartree_potential_eri


!********************************************************************************
!
! Calculate Hartree potential from Cholesky vectors in the canonical basis
!
! Information about electrons provided in self%des
!
! Actual implementation in hamilton_operations.f
!
!********************************************************************************
subroutine hartree_potential_gres(self, pot, first_state, last_state)
  class(Hamilton), intent(inout)       :: self
  real(wp), allocatable, intent(inout) :: pot(:,:,:)
  integer, optional, intent(in)        :: first_state
  integer, optional, intent(in)        :: last_state
  !local variables
  real(wp), pointer                    :: h2_gres_ob(:,:,:,:)

  call self%get_h2_gres_ob(h2_gres_ob)

  call calculate_hartree_potential_gres(h2_gres_ob, [1,1], self%des%nel, self%des%ispin_fock, self%des%rspin, pot, first_state, last_state)
end subroutine hartree_potential_gres


!********************************************************************************
!
! Calculate exchange potential from ERIs in the canonical basis
!
! Information about electrons provided in self%des
!
! Actual implementation in hamilton_operations.f
!
!********************************************************************************
subroutine exchange_potential_eri(self, pot, first_state, last_state)
  class(Hamilton), intent(inout)       :: self
  real(wp), allocatable, intent(inout) :: pot(:,:,:)
  integer, optional, intent(in)        :: first_state
  integer, optional, intent(in)        :: last_state

  call calculate_exchange_potential(self%h2, [1,1], self%des%nel, self%des%ispin_fock, pot, first_state, last_state)
end subroutine exchange_potential_eri


!********************************************************************************
!
! Calculate exchange potential from Cholesky vectors in the canonical basis
!
! Information about electrons provided in self%des
!
! Actual implementation in hamilton_operations.f
!
!********************************************************************************
subroutine exchange_potential_gres(self, pot, first_state, last_state)
  class(Hamilton), intent(inout)       :: self
  real(wp), allocatable, intent(inout) :: pot(:,:,:)
  integer, optional, intent(in)        :: first_state
  integer, optional, intent(in)        :: last_state
  !local variables
  real(wp), pointer                    :: h2_gres_ob(:,:,:,:)

  call self%get_h2_gres_ob(h2_gres_ob)

  call calculate_exchange_potential_gres(h2_gres_ob, [1,1], self%des%nel, self%des%ispin_fock, pot, first_state, last_state)
end subroutine exchange_potential_gres


!********************************************************************************
!
! Calculate self potential from ERIs
!
! Actual implementation in hamilton_operations.f
!
!********************************************************************************
subroutine self_potential_eri(self, pot, first_state, last_state)
  class(Hamilton), intent(inout)       :: self
  real(wp), allocatable, intent(inout) :: pot(:,:,:)
  integer, optional, intent(in)        :: first_state
  integer, optional, intent(in)        :: last_state
  !local variables
  integer                              :: n

  n = self%des%n

  call calculate_exchange_potential(self%h2, [1,1], [n,n], self%des%ispin1, pot, first_state, last_state)
  pot = -0.5_wp * pot
end subroutine self_potential_eri


!********************************************************************************
!
! Calculate self potential from Cholesky vectors
!
! Actual implementation in hamilton_operations.f
!
!********************************************************************************
subroutine self_potential_gres(self, pot, first_state, last_state)
  class(Hamilton), intent(inout)       :: self
  real(wp), allocatable, intent(inout) :: pot(:,:,:)
  integer, optional, intent(in)        :: first_state
  integer, optional, intent(in)        :: last_state
  !local variables
  integer                              :: n
  real(wp), pointer                    :: h2_gres_ob(:,:,:,:)


  n = self%des%n
  call self%get_h2_gres_ob(h2_gres_ob)

  call calculate_exchange_potential_gres(h2_gres_ob, [1,1], [n,n], self%des%ispin1, pot, first_state, last_state)
  pot = -0.5_wp * pot
end subroutine self_potential_gres


!********************************************************************************
!
! Calculate Hartree potential from ERIs in the canonical basis
!
! Information about electrons provided in self%des
!
! Actual implementation in hamilton_operations.f
!
!********************************************************************************
subroutine hartree_potential_rdm_eri(self, rdm, pot, first_state, last_state)
  class(Hamilton), intent(inout)       :: self
  real(wp), intent(in)                 :: rdm(:,:,:)
  real(wp), allocatable, intent(inout) :: pot(:,:,:)
  integer, optional, intent(in)        :: first_state
  integer, optional, intent(in)        :: last_state

  call calculate_hartree_potential_rdm(self%h2, rdm, self%des%ispin_fock, self%des%rspin, pot, first_state, last_state)
end subroutine hartree_potential_rdm_eri


!********************************************************************************
!
! Calculate Hartree potential from Cholesky vectors in the canonical basis
!
! Information about electrons provided in self%des
!
! Actual implementation in hamilton_operations.f
!
!********************************************************************************
subroutine hartree_potential_rdm_gres(self, rdm, pot, first_state, last_state)
  class(Hamilton), intent(inout)       :: self
  real(wp), intent(in)                 :: rdm(:,:,:)
  real(wp), allocatable, intent(inout) :: pot(:,:,:)
  integer, optional, intent(in)        :: first_state
  integer, optional, intent(in)        :: last_state
  !local variables
  real(wp), pointer                    :: h2_gres_ob(:,:,:,:)

  call self%get_h2_gres_ob(h2_gres_ob)

  call calculate_hartree_potential_rdm_gres(h2_gres_ob, rdm, self%des%ispin_fock, self%des%rspin, pot, first_state, last_state)
end subroutine hartree_potential_rdm_gres


!********************************************************************************
!
! Calculate exchange potential from ERIs in the canonical basis
!
! Information about electrons provided in self%des
!
! Actual implementation in hamilton_operations.f
!
!********************************************************************************
subroutine exchange_potential_rdm_eri(self, rdm, pot, first_state, last_state)
  class(Hamilton), intent(inout)       :: self
  real(wp), intent(in)                 :: rdm(:,:,:)
  real(wp), allocatable, intent(inout) :: pot(:,:,:)
  integer, optional, intent(in)        :: first_state
  integer, optional, intent(in)        :: last_state

  call calculate_exchange_potential_rdm(self%h2, rdm, self%des%ispin_fock, pot, first_state, last_state)
end subroutine exchange_potential_rdm_eri


!********************************************************************************
!
! Calculate exchange potential from Cholesky vectors in the canonical basis
!
! Information about electrons provided in self%des
!
! Actual implementation in hamilton_operations.f
!
!********************************************************************************
subroutine exchange_potential_rdm_gres(self, rdm, pot, first_state, last_state)
  class(Hamilton), intent(inout)       :: self
  real(wp), intent(in)                 :: rdm(:,:,:)
  real(wp), allocatable, intent(inout) :: pot(:,:,:)
  integer, optional, intent(in)        :: first_state
  integer, optional, intent(in)        :: last_state
  !local variables
  real(wp), pointer                    :: h2_gres_ob(:,:,:,:)

  call self%get_h2_gres_ob(h2_gres_ob)

  call calculate_exchange_potential_rdm_gres(h2_gres_ob, rdm, self%des%ispin_fock, pot, first_state, last_state)
end subroutine exchange_potential_rdm_gres


!******************************************************************************** 
!
! Mean-field subtraction:
!
!    Transformation of the two-body Hamiltonian
!        H2 = 0.5 \sum_g L_g^2  ==> 0.5 \sum_g (L_g - l_g)^2
!    by an arbitrary real-valued shift l_g
!
!    For a specific choice l_g = <Psi_HF|L_g|Psi_HF>, the mean-field
!    Hartree potential is shifted to the one-body Hamiltonian
!
! Original Hamiltonian:
!    H = H_0 + H_1 + 0.5 \sum_g L_g^2
!
! Rearrangerd Hamiltonian:
!    H = H'_0 + H'_1 + 0.5 \sum_g (L_g - l_g)^2
!
!    H'_0 = H_0 - E_H
!    H'_1 = H_1 + V_H
!    L_g  ==> L_g - l_g
!
!    E_H = 0.5 \sum_g l_g^2
!    V_H = \sum_g l_g L_g
!
! Note: this routine assumes g-resolved Hamiltonian,
!       and the full Hamiltonian H is not modified, just rearranged
!
!******************************************************************************** 
subroutine subtract_mean_field_r(self, lgmean, h0mean)
  class(Hamilton), intent(inout)  :: self
  real(wp), intent(in)            :: lgmean(:)
  real(wp), optional, intent(out) :: h0mean
  !local variables
  integer                         :: g, spin
  real(wp), pointer               :: h2_gres_ob(:,:,:,:)

  call self%get_h2_gres_ob(h2_gres_ob)

  do spin = 1, self%des%ispin1
    do g = 1, size(lgmean)
      self%h1(:,:,spin) = self%h1(:,:,spin) + lgmean(g) * h2_gres_ob(:,:,g,spin)
    end do
  end do

  if (present(h0mean)) h0mean = -0.5_wp * sum(lgmean**2)
end subroutine subtract_mean_field_r


!******************************************************************************** 
!
! Mean-field subtraction:
!
!    Transformation of the two-body Hamiltonian
!        H2 = 0.5 \sum_g L_g L_g^{+}  ==> 0.5 \sum_g (L_g - l_g) (L_g - l_g)^{+}
!    by an arbitrary complex-valued shift l_g
!
!    For a specific choice l_g = <Psi_HF|L_g|Psi_HF>, the mean-field
!    Hartree potential is shifted to the one-body Hamiltonian
!
! Original Hamiltonian:
!    H = H_0 + H_1 + 0.5 \sum_g L_g L_g^{+}
!
! Rearranged Hamiltonian:
!    H = H'_0 + H'_1 + 0.5 \sum_g (L_g - l_g)(L_g - l_g)^{+}
!
!    H'_0 = H_0 - E_H
!    H'_1 = H_1 + V_H
!    L_g  ==> L_g - l_g
!
!    E_H = 0.5 \sum_g |l_g|^2
!    V_H = \sum_g Rel_g L_g
!
! Note: this routine assumes g-resolved Hamiltonian,
!       and the full Hamiltonian H is not modified, just rearranged
!
!******************************************************************************** 
subroutine subtract_mean_field_c(self, lgmean, h0mean)
  class(Hamilton), intent(inout)  :: self
  complex(wp), intent(in)         :: lgmean(:)
  real(wp), optional, intent(out) :: h0mean
  !local variables
  integer                         :: g, spin
  real(wp), pointer               :: h2_gres_ob(:,:,:,:)

  call self%get_h2_gres_ob(h2_gres_ob)

  do spin = 1, self%des%ispin1
    do g = 1, size(lgmean)
      self%h1(:,:,spin) = self%h1(:,:,spin) + real(lgmean(g), wp) * h2_gres_ob(:,:,g,spin)
    end do
  end do

  if (present(h0mean)) h0mean = -0.5_wp * sum(abs(lgmean)**2)
end subroutine subtract_mean_field_c


!********************************************************************************
!
! Extract Hamiltonian in CAS space 
!
!    Input:
!        h - original Hamiltonian
!        n - dimension of the new Hamiltonian
!
!    Output:
!        hcas - Hamiltonian in CAS space
!
!******************************************************************************** 
subroutine get_cas_ham(h, hcas, n)
  class(Hamilton), intent(inout) :: h
  type(Hamilton), intent(inout)  :: hcas
  integer, intent(in)            :: n
  !local variables
  real(wp), pointer              :: h2_gres_ob(:,:,:,:), hcas_gres_ob(:,:,:,:)

  hcas%des = h%des
  hcas%des%n = n
  hcas%des%nbtot = n

  hcas%h_io = h%h_io
  hcas%enuc = h%enuc
  hcas%h0 = h%h0

  if (allocated(h%overlap)) then
    allocate(hcas%overlap(n,n,hcas%des%ispin1))
    hcas%overlap = h%overlap(1:n,1:n,:)
  end if

  if (allocated(h%ext_pot)) then
    allocate(hcas%ext_pot(n,n,hcas%des%ispin1))
    hcas%ext_pot = h%ext_pot(1:n,1:n,:)
  end if

  if (allocated(h%kinetic)) then
    allocate(hcas%kinetic(n,n,hcas%des%ispin1))
    hcas%kinetic = h%kinetic(1:n,1:n,:)
  end if

  if (allocated(h%hself)) then
    allocate(hcas%hself(n,n,hcas%des%ispin1))
    hcas%hself = h%hself(1:n,1:n,:)
  end if

  if (allocated(h%hpot)) then
    allocate(hcas%hpot(n,n,hcas%des%ispin1))
    hcas%hpot = h%hpot(1:n,1:n,:)
  end if

  if (allocated(h%xpot)) then
    allocate(hcas%xpot(n,n,hcas%des%ispin1))
    hcas%xpot = h%xpot(1:n,1:n,:)
  end if

  if (allocated(h%h1)) then
    allocate(hcas%h1(n,n,hcas%des%ispin1))
    hcas%h1 = h%h1(1:n,1:n,:)
  end if

  if (allocated(h%h1_orig)) then
    allocate(hcas%h1_orig(n,n,hcas%des%ispin1))
    hcas%h1_orig = h%h1_orig(1:n,1:n,:)
  end if

  if (allocated(h%h2)) then
    allocate(hcas%h2(n,n,n,n,hcas%des%ispin2))
    call hcas%associate_h2_ptr()
    hcas%h2 = h%h2(1:n,1:n,1:n,1:n,:)
  end if

  if (allocated(h%h2_gres)) then
    allocate(hcas%h2_gres(n*n, hcas%des%ng, hcas%des%ispin1))
  call h%get_h2_gres_ob(h2_gres_ob)
  call hcas%get_h2_gres_ob(hcas_gres_ob)
    hcas_gres_ob = h2_gres_ob(1:n,1:n,:,:)
  end if

  if (allocated(h%eigenval)) then
    allocate(hcas%eigenval(n,hcas%des%ispin))
    hcas%eigenval = h%eigenval(1:n,:)
  end if

  if (allocated(h%phi)) then
    allocate(hcas%phi(n,n,hcas%des%ispin))
    hcas%phi = h%phi(1:n,1:n,:)
  end if

  if (allocated(h%occ)) then
    allocate(hcas%occ(n,2))
    hcas%occ = h%occ(1:n,:)
  end if

  if (allocated(h%p1)) then
    allocate(hcas%p1(n,n,hcas%des%ispin1))
    hcas%p1 = h%p1(1:n,1:n,:)
  end if
end subroutine get_cas_ham

end module hamilton_type
