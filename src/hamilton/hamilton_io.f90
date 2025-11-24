! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module hamilton_io

use constants
use file_handle, only: FileHandle
use mpi, only: comm_world
use array_file_io, only: ArrayFileIO
use hamilton_vars
use hamilton_layout

implicit none

private 
public hamil_io

type hamil_io
  type(FileHandle)           :: ham_info
  type(FileHandle)           :: h0
  type(FileHandle)           :: h1
  type(FileHandle)           :: h2
  type(FileHandle)           :: overlap
  type(FileHandle)           :: orbitals
  type(FileHandle)           :: eigenval
  type(FileHandle)           :: occup
  
  character(len=charlen)     :: file_format = "numpy"
  integer                    :: iwrite = 1
contains
  procedure                  :: init_files
  procedure                  :: if_files_exist
  procedure                  :: read_ham_info  => read_ham_info_mpi
  procedure                  :: write_ham_info
  procedure                  :: read_h2
  procedure                  :: write_h2
end type hamil_io

contains

!******************************************************************************** 
!
!       Initialize hamiltonian file handles
!
!******************************************************************************** 
subroutine init_files(h_io)
  class(hamil_io), intent(inout) :: h_io
  
  h_io%ham_info = FileHandle("ham_info")
  h_io%h0 = FileHandle("ham0")
  h_io%h1 = FileHandle("ham1")
  h_io%h2 = FileHandle("ham2")
  h_io%overlap = FileHandle("overlap")
  h_io%orbitals = FileHandle("orbitals")
  h_io%occup = FileHandle("occup")
  h_io%eigenval = FileHandle("eigenval")
end subroutine init_files


!******************************************************************************** 
!
!       Check whether hamiltonian files exist
!
!       If ham_info, ham1 and ham2 are present, hamiltonian files can be read
!
!       ham0 is not explicitly needed - it will be assumed to be zero
!
!******************************************************************************** 
logical function if_files_exist(h_io)
  class(hamil_io), intent(in) :: h_io
  
  if_files_exist = h_io%ham_info%exists() .and. h_io%h1%exists() .and. h_io%h2%exists()
end function if_files_exist


!******************************************************************************** 
!
!       Write ham_info file
!
!       Contains all info about Hamiltonian - has to be written with other files
!
!******************************************************************************** 
subroutine write_ham_info(h_io, ham_des)
  class(hamil_io), intent(inout)     :: h_io
  type(hamil_descriptor), intent(in) :: ham_des
  
  call h_io%ham_info%open(status="replace", action="write")
  write(h_io%ham_info%funit, "(1x,t1,a,a)")      "basis         = ", trim(ham_des%basis)
  write(h_io%ham_info%funit, "(1x,t1,a,a)")      "integral_mode = ", trim(ham_des%integral_mode)
  write(h_io%ham_info%funit, "(1x,t1,a,a)")      "file_format   = ", trim(h_io%file_format)
  write(h_io%ham_info%funit, "(1x,t1,a,i6)")     "n             = ", ham_des%n
  write(h_io%ham_info%funit, "(1x,t1,a,i6)")     "nbtot         = ", ham_des%nbtot
  write(h_io%ham_info%funit, "(1x,t1,a,i6)")     "ispin         = ", ham_des%ispin
  write(h_io%ham_info%funit, "(1x,t1,a,i6)")     "ispin1        = ", ham_des%ispin1
  write(h_io%ham_info%funit, "(1x,t1,a,i6)")     "ispin2        = ", ham_des%ispin2
  if (trim(ham_des%integral_mode) == "cholesky") then
    write(h_io%ham_info%funit, "(1x,t1,a,i6)")   "ng            = ", ham_des%ng
  end if 
  write(h_io%ham_info%funit, "(1x,t1,a,2i6)")    "nel           = ", ham_des%nel
  call h_io%ham_info%close()
end subroutine write_ham_info


!******************************************************************************** 
!
!       Read ham_info file
!     
!       First file read before all other hamiltonian files
!
!******************************************************************************** 
subroutine read_ham_info(h_io, ham_des)
  use string, only: split_tag
  class(hamil_io), intent(inout)        :: h_io
  type(hamil_descriptor), intent(inout) :: ham_des
  !local variables
  character(len=charlen)                :: key, value_
  integer                               :: ispin
  !
  if (.not. h_io%ham_info%exists()) return
  call h_io%ham_info%open(status="old", action="read")
  !
  read(h_io%ham_info%funit, "(a)") key
  call split_tag(key, value_)
  read(value_, *) ham_des%basis
  !
  read(h_io%ham_info%funit, "(a)") key
  call split_tag(key, value_)
  read(value_, *) ham_des%integral_mode
  !
  read(h_io%ham_info%funit, "(a)") key
  call split_tag(key, value_)
  read(value_, *) h_io%file_format
  !
  read(h_io%ham_info%funit, "(a)") key
  call split_tag(key, value_)
  read(value_, *) ham_des%n
  !
  read(h_io%ham_info%funit, "(a)") key
  call split_tag(key, value_)
  read(value_, *) ham_des%nbtot
  !
  read(h_io%ham_info%funit, "(a)") key
  call split_tag(key, value_)
  read(value_, *) ispin
  if (ham_des%ispin == unset) ham_des%ispin = ispin
  !
  read(h_io%ham_info%funit, "(a)") key
  call split_tag(key, value_)
  read(value_, *) ham_des%ispin1
  !
  read(h_io%ham_info%funit, "(a)") key
  call split_tag(key, value_)
  read(value_, *) ham_des%ispin2
  !
  if (trim(ham_des%integral_mode) == "cholesky") then
    read(h_io%ham_info%funit, "(a)") key
    call split_tag(key, value_)
    read(value_, *) ham_des%ng
  end if
  !
  read(h_io%ham_info%funit, "(a)") key
  call split_tag(key, value_)
  read(value_, *) ham_des%nel
  !
  call h_io%ham_info%close()
end subroutine read_ham_info


!******************************************************************************** 
!
!       MPI interface to read ham_info
!       
!       At the moment, sequential read by each proc is ensured
!
!******************************************************************************** 
subroutine read_ham_info_mpi(h_io, ham_des)
  class(hamil_io), intent(inout)        :: h_io
  type(hamil_descriptor), intent(inout) :: ham_des
  !local variables
  integer                               :: proc
  !
  do proc = 0, comm_world%mpisize-1
    if (proc == comm_world%mpirank) call read_ham_info(h_io, ham_des)
    call comm_world%barrier()
  end do
end subroutine read_ham_info_mpi


!******************************************************************************** 
!       
!       Read two-body Hamiltonian
!       
!       integral_mode read from ham_info or alternatively from the input file
!
!******************************************************************************** 
subroutine read_h2(h_io, ham_des, h2, h2_gres, array_io)
  class(hamil_io), intent(inout)       :: h_io
  type(hamil_descriptor), intent(in)   :: ham_des
  real(wp), allocatable, intent(inout) :: h2(:,:,:,:,:), h2_gres(:,:,:)
  class(ArrayFileIO), intent(in)       :: array_io
  !local variables
  real(wp), allocatable                :: h2_gres_(:,:,:,:)
  
  if (trim(ham_des%integral_mode) == "eri") then
    call array_io%load(h_io%h2%fname, h2)
  else if (trim(ham_des%integral_mode) == "cholesky") then
    call array_io%load(h_io%h2%fname, h2_gres_)
    call get_h2_gres_2d(h2_gres_, h2_gres, .true.)
  else 
    write(*,*) "If this message pops up - something very strange appeared - Program stops now"
    call exit
  end if
end subroutine read_h2


!******************************************************************************** 
!
!       Write two-body Hamiltonian
!
!******************************************************************************** 
subroutine write_h2(h_io, ham_des, h2, h2_gres, array_io)
  class(hamil_io), intent(inout)       :: h_io
  type(hamil_descriptor), intent(in)   :: ham_des
  real(wp), allocatable, intent(inout) :: h2(:,:,:,:,:), h2_gres(:,:,:)
  class(ArrayFileIO), intent(in)       :: array_io
  !local variables
  real(wp), allocatable                :: h2_gres_3d(:,:,:,:)
  
  if (trim(ham_des%integral_mode)=="eri" .or. allocated(h2)) then
    call array_io%save(h_io%h2%fname, h2)
  else if (trim(ham_des%integral_mode)=="cholesky" .and. allocated(h2_gres)) then
    call  get_h2_gres_3d(h2_gres, h2_gres_3d, ham_des%n, ham_des%n, .false.)
    call array_io%save(h_io%h2%fname, h2_gres_3d)
  else
    write(*,*) "If this message pops up - something very strange appeared - Program stops now"
    call exit
  end if
end subroutine write_h2

end module hamilton_io
