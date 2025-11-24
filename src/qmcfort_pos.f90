! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module qmcfort_pos

use constants
use qmcfort_io
use file_handle, only: FileHandle
use mpi, only: mpi_communicator, comm_world
use string, only: lowercase

implicit none

private
public :: Structure, get_atomic_number, get_atomic_mass
public :: read_qmcfort_structure, qmcfort_structure_factory

type Structure 
  character                         :: units                !Units of the coordinates
  character(len=atomc), allocatable :: atoms(:)             !List of symbols for each atom
  real(wp), allocatable             :: coords(:,:)          !Coordinates of each atom
  
  real(wp), allocatable             :: Z(:)                 !Atomic numbers
  real(wp), allocatable             :: Zeff(:)              !Effective atomic numbers
  real(wp), allocatable             :: mass(:)              !Atomic mass
contains
  procedure          :: get_natoms
  procedure          :: get_nspecies
  procedure          :: get_species
  procedure          :: get_species_str
  procedure          :: get_natoms_for_species
  procedure          :: get_coords_for_species
  procedure          :: get_core_states
  procedure          :: calc_enuc  => calculate_nuclear_energy
  procedure          :: mpi_bcast => mpi_bcast_structure
  procedure          :: report => report_structure
  procedure, private :: set_atomic_numbers
  procedure, private :: set_atomic_masses
end type structure

interface Structure
  module procedure finit_structure
end interface Structure

contains

!********************************************************************************
!
! Initialization of the Structure object
!
!********************************************************************************
subroutine init_structure(units, coords_str, self)
  character(len=*), intent(in) :: units
  character(len=*), intent(in) :: coords_str
  type(Structure), intent(out) :: self
  !local
  integer, parameter           :: natoms_max = 1000
  integer                      :: natoms, istart, iend, ierror
  character(len=atomc)         :: atom_symbol, atoms(natoms_max)
  character(len=charlen)       :: str_block
  real(wp)                     :: coords(3,natoms_max)
  real(wp)                     :: Zeff(natoms_max)

  natoms = 0
  istart = 1

  !read coord_str and extract atoms and coord
  do
    !early exit for empty structure
    if (len(trim(coords_str)) == 0) exit

    !identify part of the string belonging to a single atom
    iend = index(coords_str(istart:), ';')
    if (iend == 0) then
      str_block = adjustl(coords_str(istart:))
    else
      str_block = adjustl(coords_str(istart:istart+iend-2))
    end if

    !Read symbol and coords from the string
    natoms = natoms + 1
    read(str_block, *, iostat=ierror) atom_symbol, coords(1,natoms), coords(2,natoms), coords(3,natoms), Zeff(natoms)
    if (ierror == 0) then
      !to be consistent with cp correction setup for PyScf
      if (lowercase(atom_symbol(1:2)) == "x-") then
        atoms(natoms) = atom_symbol(3:)
      else
        atoms(natoms) = atom_symbol
      end if
    else
      read(str_block, *, iostat=ierror) atom_symbol, coords(1,natoms), coords(2,natoms), coords(3,natoms)
      if (ierror == 0) then
        !to be consistent with cp correction setup for PyScf
        if (lowercase(atom_symbol(1:2)) == "x-") then
          atoms(natoms) = atom_symbol(3:)
          Zeff(natoms) = 0.0_wp
        else
          atoms(natoms) = atom_symbol
          Zeff(natoms) = get_atomic_number(atoms(natoms))
        end if
      end if
    end if

    !Advance to the next atom string
    istart = istart + iend
    if (iend == 0) exit
  end do

  self%units = lowercase(units)

  allocate(self%atoms(natoms))
  if (natoms > 0) self%atoms = atoms(1:natoms)

  allocate(self%coords(3,natoms))
  if (natoms > 0) then
    self%coords = coords(:,1:natoms)
    if (self%units == "a")  self%coords = angstrom_to_bohr * self%coords
  end if

  allocate(self%Zeff(natoms))
  if (natoms > 0) self%Zeff = Zeff(1:natoms)

  call self%set_atomic_numbers()
  call self%set_atomic_masses()
end subroutine init_structure


!********************************************************************************
!
! Structure constructor
!
!********************************************************************************
function finit_structure(units, coords_str) result(self)
  character(len=*), intent(in) :: units
  character(len=*), intent(in) :: coords_str
  type(Structure)              :: self

  call init_structure(units, coords_str, self)
end function finit_structure


!********************************************************************************
!
! Get number of atoms in Structure 
!
!********************************************************************************
function get_natoms(self) result(natoms)
  class(Structure), intent(in) :: self
  integer                      :: natoms

  natoms = size(self%coords, 2)
end function get_natoms


!********************************************************************************
!
! Get number of species in Structure 
!
!********************************************************************************
function get_nspecies(self) result(nspecies)
  class(Structure), intent(in) :: self
  integer                      :: nspecies
  !local
  integer, parameter           :: nspecies_max = 20
  integer                      :: i, j
  integer                      :: Z(nspecies_max)
  logical                      :: unique

  nspecies = 0
   
  do i = 1, self%get_natoms()
    unique = .true.

    do j = 1, nspecies
       if (self%Z(i) == Z(j)) unique = .false.
    end do

    if (unique) then
      nspecies = nspecies + 1
      Z(nspecies) = self%Z(i)
    end if
  end do
end function get_nspecies


!********************************************************************************
!
! Get atomic species in Structure 
!
!********************************************************************************
subroutine get_species(self, species)
  class(Structure), intent(in)                   :: self
  character(len=atomc), allocatable, intent(out) :: species(:)
  !local
  integer, parameter                             :: nspecies_max = 20
  integer                                        :: i, j, nspecies
  integer                                        :: Z(nspecies_max)
  logical                                        :: unique
  character(len=atomc)                           :: species_(nspecies_max)

  nspecies = 0
   
  do i = 1, self%get_natoms()
    unique = .true.

    do j = 1, nspecies
       if (self%Z(i) == Z(j)) unique = .false.
    end do

    if (unique) then
      nspecies = nspecies + 1
      Z(nspecies) = self%Z(i)
      species_(nspecies) = self%atoms(i)
    end if
  end do

  allocate(species(nspecies))
  species = species_(1:nspecies)
end subroutine get_species


!********************************************************************************
!
! Get a string containing all atomic species in the system
!
!********************************************************************************
function get_species_str(self) result(species_str)
  class(Structure), intent(in)      :: self
  character(len=charlen)            :: species_str
  !local
  integer                           :: i
  character(len=atomc), allocatable :: species(:)
  
  call self%get_species(species)
  species_str = "" 

  do i = 1, size(species)
    species_str = trim(species_str) // " " // trim(species(i))
  end do
end function get_species_str


!********************************************************************************
!
! Get number of atoms for a given species
!
!********************************************************************************
function get_natoms_for_species(self, species) result(natoms)
  class(Structure), intent(in) :: self
  character(len=*), intent(in) :: species
  integer                      :: natoms
  !local
  integer                      :: i

  natoms = 0

  do i = 1, self%get_natoms()
    if (species == self%atoms(i)) natoms = natoms + 1
  end do
end function get_natoms_for_species


!********************************************************************************
!
! Get coords for a given species
!
!********************************************************************************
function get_coords_for_species(self, species) result(coords)
  class(Structure), intent(in) :: self
  character(len=*), intent(in) :: species
  real(wp), allocatable        :: coords(:,:)
  !local
  integer                      :: i, j, natoms

  natoms = self%get_natoms_for_species(species)
  allocate(coords(3, natoms))

  j = 0
  do i = 1, self%get_natoms()
    if (species == self%atoms(i)) then
      j = j +1
      coords(:,j) = self%coords(:,i)
    end if
  end do
end function get_coords_for_species


!********************************************************************************
!
! Calculate classical repulsion energy between nuclei in the system
!
!********************************************************************************
function calculate_nuclear_energy(self) result(enuc)
  class(Structure), intent(in) :: self
  real(wp)                     :: enuc
  !local
  integer                      :: i, j, natoms

  natoms = self%get_natoms()
  enuc = 0.0_wp

  do i = 1, natoms-1
    do j = i+1, natoms
      enuc = enuc + self%Zeff(i)*self%Zeff(j) / norm2(self%coords(:,j) - self%coords(:,i))
    end do
  end do
end function calculate_nuclear_energy


!********************************************************************************
!
! MPI Broadcast Structure object
!
!********************************************************************************
subroutine mpi_bcast_structure(self, comm, root)
  class(Structure), intent(inout)    :: self
  type(mpi_communicator), intent(in) :: comm
  integer, optional, intent(in)      :: root
  !local
  integer                            :: root_
  integer                            :: natoms

  root_ = 0
  if (present(root)) root_ = root

  if (comm%mpirank == root_) natoms = self%get_natoms()
  call comm%bcast(natoms, root_)

  if (comm%mpirank /= root_) then
    allocate(self%atoms(natoms))
    allocate(self%coords(3,natoms))
    allocate(self%Z(natoms))
    allocate(self%Zeff(natoms))
    allocate(self%mass(natoms))
  end if

  call comm%bcast(self%units, root_)
  call comm%bcast(self%atoms, root_)
  call comm%bcast(self%coords, root_)
  call comm%bcast(self%Z, root_)
  call comm%bcast(self%Zeff, root_)
  call comm%bcast(self%mass, root_)
end subroutine mpi_bcast_structure


!********************************************************************************
!
! Report Structure object
!
!********************************************************************************
subroutine report_structure(self, comm, fh)
  class(Structure), intent(in)       :: self
  type(mpi_communicator), intent(in) :: comm
  type(FileHandle), intent(in)       :: fh
  !local
  integer                            :: i
  character(len=charlen)             :: species_str
  
  if (comm%mpirank == 0) then
    species_str = self%get_species_str()

    write(fh%funit, 100) ""
    write(fh%funit, 100) "Molecular Geometry:"
    write(fh%funit, 100) "-------------------"
    write(fh%funit, 101) "  number of atoms", self%get_natoms()
    write(fh%funit, 101) "  number of atomic species", self%get_nspecies()
    write(fh%funit, 102) "  atomic species ", trim(species_str)
    write(fh%funit, 102) "  coordinate units", "A"
    write(fh%funit, 100) "  coordinates:"
    write(fh%funit, 103) "   atom   ", "     x    ", "     y    ", "     z    ", "   Zeff   "
    write(fh%funit, 103) "   ====   ", "     =    ", "     =    ", "     =    ", "   ====   "
    do i = 1, self%get_natoms()
      write(fh%funit, 104) self%atoms(i), self%coords(:,i) / angstrom_to_bohr, self%Zeff(i)
    end do
    write(fh%funit, 100) ""
  end if

  100 format (1x,a)
  101 format (1x,a,t50,"= ",i4)
  102 format (1x,a,t50,"= ",a)
  103 format (1x,t5,a,t15,a,t25,a,t35,a,t45,a)
  104 format (1x,t9,a,t15,f10.6,t25,f10.6,t35,f10.6,t45,f10.6)
end subroutine report_structure


!********************************************************************************
!
! Get the number of core states in the Molecular Structure
!
!********************************************************************************
function get_core_states(self) result(core_states)
  class(Structure), intent(in) :: self
  integer                      :: core_states
  !local
  integer                      :: i

  core_states = 0

  do i = 1, size(self%atoms)
    !don't count ghost states for the frozen-core
    if (int(self%Zeff(i)) == 0) cycle 

    core_states = core_states + get_atomic_core_states(self%atoms(i))
  end do
end function get_core_states


!********************************************************************************
!
! Set atomic numbers array in Structure object
!
!********************************************************************************
subroutine set_atomic_numbers(self)
  class(Structure), intent(inout) :: self
  !local
  integer                         :: i, natoms

  natoms = self%get_natoms()

  if (allocated(self%Z)) deallocate(self%Z)
  allocate(self%Z(natoms))

  do i = 1, natoms
    self%Z(i) = get_atomic_number(self%atoms(i))
  end do
end subroutine set_atomic_numbers


!********************************************************************************
!
! Set atomic masses array in Structure object
!
!********************************************************************************
subroutine set_atomic_masses(self)
  class(Structure), intent(inout) :: self
  !local
  integer                         :: i, natoms

  natoms = self%get_natoms()

  if (allocated(self%mass)) deallocate(self%mass)
  allocate(self%mass(natoms))

  do i = 1, natoms
    self%mass(i) = get_atomic_mass(self%atoms(i))
  end do
end subroutine set_atomic_masses


!********************************************************************************
!
! Read units and coords_str from qmcfort.xyz or qmcfort_pos files 
! and 
! create Structe object
!
!********************************************************************************
subroutine qmcfort_structure_factory(struc)
  type(Structure), intent(out) :: struc
  !local
  integer, parameter           :: root = 0
  character                    :: units
  character(len=longc)         :: coords_str
  type(FileHandle)             :: fh1, fh2

  if (comm_world%mpirank == root) then
    fh1 = FileHandle("qmcfort.xyz")
    fh2 = FileHandle("qmcfort_pos")

    if (fh1%exists()) then
      call read_qmcfort_xyz_file(units, coords_str)
    else if (fh2%exists()) then
      call read_qmcfort_pos_file(units, coords_str)
    else 
      units = "a"
      coords_str = ""
    end if

    struc = Structure(units, coords_str)
  end if

  call struc%mpi_bcast(comm_world, root)
end subroutine qmcfort_structure_factory


!********************************************************************************
!
! Read Structure from files and report it to the standard output
!
!********************************************************************************
subroutine read_qmcfort_structure(struc)
  type(Structure), intent(out) :: struc
  call qmcfort_structure_factory(struc)

  call struc%report(comm_world, io%screen)
  call struc%report(comm_world, io%qmcfort_log)
  call struc%report(comm_world, io%qmcfort_out)
end subroutine read_qmcfort_structure


!********************************************************************************
!
! qmcfort.xyz file reader
!
!    Creates two strings needed for the initialization of the Structure:
!        units,
!        coords_str
!        
!********************************************************************************
subroutine read_qmcfort_xyz_file(units, coords_str)
  character, intent(out)            :: units
  character(len=longc), intent(out) :: coords_str
  !local
  integer                           :: i, natoms, pos1, ios
  character(len=charlen)            :: line
  type(FileHandle)                  :: fh
  character(len=*), parameter       :: fname = "qmcfort.xyz"

  units = "a"
  fh = FileHandle(fname)

  !ensure the file is already present
  call fh%open(status="old", action="read")

  !read number of atoms
  read(fh%funit,"(a)") line
  read(line,*) natoms

  !read comment line
  read(fh%funit,"(a)") line
  pos1 = index(line, "units=")
  if (pos1 > 0) then
    pos1 = pos1 + len("units=")
    units = line(pos1:pos1+1)
  end if

  i = 0
  coords_str = ""

  !read atom by atom
  do
    read(fh%funit,"(a)", iostat=ios) line
    if (ios /= 0) exit
    if (len(trim(line)) <= 1) cycle

    i = i + 1
    if (i == natoms) then
      coords_str = trim(coords_str) // trim(line)
    else
      coords_str = trim(coords_str) // trim(line) // "; "
    end if
  end do
  
  call fh%close()
end subroutine read_qmcfort_xyz_file


!********************************************************************************
!
! qmcfort_pos file reader
!
!    Creates two strings needed for the initialization of the Structure:
!        units,
!        coords_str
!
!********************************************************************************
subroutine read_qmcfort_pos_file(units, coords_str)
  character, intent(out)            :: units
  character(len=longc), intent(out) :: coords_str
  !local
  integer                           :: i, j, natoms, nspecies, mult, ios
  character(len=atomc)              :: symbol
  character(len=charlen)            :: line, pos_str, z
  type(FileHandle)                  :: fh
  character(len=*), parameter       :: fname = "qmcfort_pos"
  
  units = "a"
  fh = FileHandle(fname)

  !ensure the file is already present
  call fh%open(status="old", action="read")

  !read and extract variables from the first line of the qmcfort_pos file
  read(fh%funit,"(a)") line
  read(line,*,iostat=ios) natoms, nspecies, units
  if (ios /= 0) read(line,*) natoms, nspecies

  coords_str = ""

  do i = 1, nspecies
    read(fh%funit,*) symbol, mult, z

    j = 0
    do 
      read(fh%funit,"(a)") pos_str
      if (len(trim(pos_str)) <= 1) cycle
      j = j + 1

      if (i==nspecies .and. j==mult) then
        coords_str = trim(coords_str) // symbol // " " // trim(pos_str) // " " // trim(z)
      else 
        coords_str = trim(coords_str) // symbol // " " // trim(pos_str) // " " // trim(z) // "; "
      end if

      if (j == mult) exit
    end do
  end do
  
  call fh%close()
end subroutine read_qmcfort_pos_file


!********************************************************************************
!
! Get atomic number for a given element
!
!********************************************************************************
function get_atomic_number(atom) result(z)
  character(len=*), intent(in) :: atom
  real(wp)                     :: z

  select case (lowercase(trim(atom)))
    case ("h")
      z = 1.0_wp
    case ("he")
      z = 2.0_wp
    case ("li")
      z = 3.0_wp
    case ("be")
      z = 4.0_wp
    case ("b")
      z = 5.0_wp
    case ("c")
      z = 6.0_wp
    case ("n")
      z = 7.0_wp
    case ("o")
      z = 8.0_wp
    case ("f")
      z = 9.0_wp
    case ("ne")
      z = 10.0_wp
    case ("na")
      z = 11.0_wp
    case ("mg")
      z = 12.0_wp
    case ("al")
      z = 13.0_wp
    case ("si")
      z = 14.0_wp
    case ("p")
      z = 15.0_wp
    case ("s")
      z = 16.0_wp
    case ("cl")
      z = 17.0_wp
    case ("ar")
      z = 18.0_wp
    case ("k")
      z = 19.0_wp
    case ("ca")
      z = 20.0_wp
    case ("sc")
      z = 21.0_wp
    case ("ti")
      z = 22.0_wp
    case ("v")
      z = 23.0_wp
    case ("cr")
      z = 24.0_wp
    case ("mn")
      z = 25.0_wp
    case ("fe")
      z = 26.0_wp
    case ("co")
      z = 27.0_wp
    case ("ni")
      z = 28.0_wp
    case ("cu")
      z = 29.0_wp
    case ("zn")
      z = 30.0_wp
    case ("ga")
      z = 31.0_wp
    case ("ge")
      z = 32.0_wp
    case ("as")
      z = 33.0_wp
    case ("se")
      z = 34.0_wp
    case ("br")
      z = 35.0_wp
    case ("kr")
      z = 36.0_wp
    case default
      write(io%screen%funit,*) "Atom is not recognized: only first four rows are introduced"
      z = 0.0_wp
  end select
end function get_atomic_number


!********************************************************************************
!
! Get number of core states for the atom
!
!********************************************************************************
function get_atomic_core_states(atom) result(zcore)
  character(len=*), intent(in) :: atom
  integer                      :: zcore

  select case (lowercase(trim(atom)))
    case ("h")
      zcore = 0
    case ("he")
      zcore = 0
    case ("li")
      zcore = 0
    case ("be")
      zcore = 0
    case ("b")
      zcore = 1
    case ("c")
      zcore = 1
    case ("n")
      zcore = 1
    case ("o")
      zcore = 1
    case ("f")
      zcore = 1
    case ("ne")
      zcore = 1
    case ("na")
      zcore = 1
    case ("mg")
      zcore = 1
    case ("al")
      zcore = 5
    case ("si")
      zcore = 5
    case ("p")
      zcore = 5
    case ("s")
      zcore = 5
    case ("cl")
      zcore = 5
    case ("ar")
      zcore = 5
    case ("k")
      zcore = 5
    case ("ca")
      zcore = 5
    case ("sc")
      zcore = 5
    case ("ti")
      zcore = 5
    case ("v")
      zcore = 5
    case ("cr")
      zcore = 5
    case ("mn")
      zcore = 5
    case ("fe")
      zcore = 5
    case ("co")
      zcore = 5
    case ("ni")
      zcore = 5
    case ("cu")
      zcore = 5
    case ("zn")
      zcore = 5
    case ("ga")
      zcore = 9
    case ("ge")
      zcore = 9
    case ("as")
      zcore = 9
    case ("se")
      zcore = 9
    case ("br")
      zcore = 9
    case ("kr")
      zcore = 9
    case default
      write(io%screen%funit,*) "Atom is not recognized: only first four rows are introduced"
      zcore = 0
  end select
end function get_atomic_core_states


!********************************************************************************
!
! Get atomic mass of the atom
!
!********************************************************************************
function get_atomic_mass(atom) result(u)
  character(len=*), intent(in) :: atom
  real(wp)                     :: u
 
  select case (lowercase(trim(atom)))
    case ("h")
      u = 1.0079_wp
    case ("he")
      u = 4.0026_wp
    case ("li")
      u = 6.941_wp
    case ("be")
      u = 9.01218_wp
    case ("b")
      u = 10.81_wp
    case ("c")
      u = 12.011_wp
    case ("n")
      u = 14.0067_wp
    case ("o")
      u = 15.9994_wp
    case ("f")
      u = 18.998403_wp
    case ("ne")
      u = 20.179_wp
    case ("na")
      u = 22.98977_wp
    case ("mg")
      u = 24.305_wp
    case ("al")
      u = 26.98154_wp
    case ("si")
      u = 28.0855_wp
    case ("p")
      u = 30.97376_wp
    case ("s")
      u = 32.064_wp
    case ("cl")
      u = 35.453_wp
    case ("ar")
      u = 39.948_wp
    case ("k")
      u = 39.0983_wp
    case ("ca")
      u = 40.08_wp
    case ("sc")
      u = 44.9559_wp
    case ("ti")
      u = 47.90_wp
    case ("v")
      u = 50.9415_wp
    case ("cr")
      u = 51.996_wp
    case ("mn")
      u = 54.9380_wp
    case ("fe")
      u = 55.847_wp
    case ("co")
      u = 58.9332_wp
    case ("ni")
      u = 58.6934_wp
    case ("cu")
      u = 63.546_wp
    case ("zn")
      u = 65.409_wp
    case ("ga")
      u = 69.72_wp
    case ("ge")
      u = 72.59_wp
    case ("as")
      u = 74.92159_wp
    case ("se")
      u = 78.96_wp
    case ("br")
      u = 79.904_wp
    case ("kr")
      u = 83.798_wp
    case default
      write(io%screen%funit,*) "Atom is not recognized: only first four rows are introduced"
      u = 0.0_wp
  end select
end function get_atomic_mass

end module qmcfort_pos