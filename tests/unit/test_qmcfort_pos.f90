! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module test_qmcfort_pos

use qmcfort_pos

use constants, only: wp, longc
use file_handle, only: FileHandle
use fruit, only: assert_equals, run_test_case, test_set_summary

implicit none

private
public :: test_qmcfort_pos_set

contains

subroutine test_qmcfort_pos_set
  call run_test_case(test_empty_structure, "Check natoms and nspecies in an empty structure")
  call run_test_case(test_structure_h2_natoms_nspecies, "Check natoms and nspecies for H2")
  call run_test_case(test_structure_h2_coords_species, "Check H coordinates in H2 molecule")
  call run_test_case(test_structure_h2o_zeff, "Check effective Z for H2O molecule with CP setting")
  call run_test_case(test_structure_unit_conv, "Check for the internal conversion of the coordinates units")
  call run_test_case(test_structure_z_ca2cl3, "Chek atomic numbers of Ca2Cl3")
  call run_test_case(test_structure_mass_caco2h3, "Check atomic masses of CaCO2H3")

  call test_set_summary("src/qmcfort_pos.f90")
end subroutine test_qmcfort_pos_set

subroutine test_empty_structure
  integer, allocatable :: expected(:)
  integer, allocatable :: result_(:)
  character            :: units
  character(len=longc) :: coords_str
  type(Structure)      :: struc
  
  units = "a"
  coords_str = ""

  struc = Structure(units, coords_str)

  expected = [0, 0]
  result_ = [struc%get_natoms(), struc%get_nspecies()]
  call assert_equals(expected, result_, size(expected))
end subroutine test_empty_structure

subroutine test_structure_h2_natoms_nspecies
  integer, allocatable :: expected(:)
  integer, allocatable :: result_(:)
  character            :: units
  character(len=longc) :: coords_str
  type(Structure)      :: struc
  
  units = "a"
  coords_str = " H 0.0 0.0 0.0;&
                 H 0.0 0.0 1.0"

  struc = Structure(units, coords_str)

  expected = [2, 1]
  result_ = [struc%get_natoms(), struc%get_nspecies()]
  call assert_equals(expected, result_, size(expected))
end subroutine test_structure_h2_natoms_nspecies

subroutine test_structure_h2_coords_species
  real(wp), allocatable :: expected(:,:)
  real(wp), allocatable :: result_(:,:)
  character             :: units
  character(len=longc)  :: coords_str
  type(Structure)       :: struc
  
  units = "a"
  coords_str = " H 0.0 0.0 0.0 ; & 
                 H 0.0 0.0 1.0  "
  
  struc = Structure(units, coords_str)
  
  expected = struc%coords
  result_ = struc%get_coords_for_species("H")

  call assert_equals(expected, result_, size(expected,1), size(expected,2))
end subroutine test_structure_h2_coords_species

subroutine test_structure_h2o_zeff
  integer, allocatable :: expected(:)
  integer, allocatable :: result_(:)
  character            :: units
  character(len=longc) :: coords_str
  type(Structure)      :: struc

  units = "a"
  coords_str = "H 0.0 0.0 0.0 0.0 ; &
                H 0.0 0.0 1.0 0.0 ; &
                O 0.0 0.0 2.0 "
  
  struc = Structure(units, coords_str)
  
  expected = [0.0_wp, 0.0_wp, 8.0_wp]
  result_ = struc%Zeff
  call assert_equals(expected, result_, size(expected))
end subroutine test_structure_h2o_zeff

subroutine test_structure_unit_conv
  real(wp), parameter   :: tol=1.0E-06_wp
  real(wp), allocatable :: expected(:,:)
  real(wp), allocatable :: result_(:,:)
  character             :: units
  character(len=longc)  :: coords_str
  type(Structure)       :: struc
  
  units = "a"
  coords_str =  "H 0.0 0.0 0.0; &
                 H 0.0 0.0 1.0; &
                 Be 0.0 0.0 2.0 "
  
  struc = Structure(units, coords_str)
  
  expected = reshape([0.0_wp, 0.0_wp, 3.779451977_wp], shape=[3,1])
  result_ = struc%get_coords_for_species("Be")
  call assert_equals(expected, result_, size(expected,1), size(expected,2), tol)
end subroutine test_structure_unit_conv

subroutine test_structure_z_ca2cl3
  integer, allocatable :: expected(:)
  integer, allocatable :: result_(:)
  type(Structure)      :: struc
  type(FileHandle)     :: fh
  
  fh = FileHandle("qmcfort_pos")
  call fh%open(status="replace", action="write")
  write(fh%funit,*) "5 2             "
  write(fh%funit,*) ""
  write(fh%funit,*) "Ca 2 20.0       "
  write(fh%funit,*) "  0.0 0.0 0.0   "
  write(fh%funit,*) "  0.0 0.0 1.0   "
  write(fh%funit,*) ""
  write(fh%funit,*) "Cl 3 17.0       "
  write(fh%funit,*) ""
  write(fh%funit,*) "  0.0 0.0 2.0 "
  write(fh%funit,*) ""
  write(fh%funit,*) "  0.0 0.0 3.0 "
  write(fh%funit,*) ""
  write(fh%funit,*) "  0.0 0.0 4.0 "
  call fh%close()
  
  call qmcfort_structure_factory(struc)
  
  expected = [20.0_wp, 20.0_wp, 17.0_wp, 17.0_wp, 17.0_wp]
  result_ = struc%Z
  call assert_equals(expected, result_, size(expected))

  call fh%delete()
end subroutine test_structure_z_ca2cl3

subroutine test_structure_mass_caco2h3
  integer, allocatable :: expected(:)
  integer, allocatable :: result_(:)
  type(structure)      :: struc
  type(FileHandle)     :: fh
  
  fh = FileHandle("qmcfort.xyz")
  call fh%open(status="replace", action="write")
  write(fh%funit,*) "7"
  write(fh%funit,*) "units=b"
  write(fh%funit,*) " "
  write(fh%funit,*) "Ca  0.0 0.0 0.0  "
  write(fh%funit,*) " "
  write(fh%funit,*) "C  0.0 0.0 1.0   "
  write(fh%funit,*) " "
  write(fh%funit,*) "O  0.0 0.0 2.0   "
  write(fh%funit,*) "O  0.0 0.0 3.0   "
  write(fh%funit,*) " "
  write(fh%funit,*) "H  0.0 0.0 4.0   "
  write(fh%funit,*) "H  0.0 0.0 5.0   "
  write(fh%funit,*) "H  0.0 0.0 6.0   "
  call fh%close()
  
  call qmcfort_structure_factory(struc)

  expected = [get_atomic_mass("Ca"), get_atomic_mass("C"), get_atomic_mass("O"), get_atomic_mass("O"), &
              get_atomic_mass("H"), get_atomic_mass("H"), get_atomic_mass("H")]
  result_ = struc%mass
  call assert_equals(expected, result_, size(expected))
  
  call fh%delete()
end subroutine test_structure_mass_caco2h3

end module test_qmcfort_pos
