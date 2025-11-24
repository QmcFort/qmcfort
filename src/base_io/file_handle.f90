! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module file_handle

use iso_fortran_env, only: input_unit, output_unit, error_unit

use constants

implicit none

#include "../preproc.inc"

private 
public :: FileHandle

type FileHandle
  character(len=:), allocatable :: fname
  integer                       :: funit = unset
contains
  procedure :: open => open_file
  procedure :: close => close_file
  procedure :: delete => delete_file
  procedure :: is_open => is_file_open
  procedure :: exists => file_exists
end type FileHandle

interface  FileHandle
  procedure :: finit_file_handle
end interface FileHandle

contains

!******************************************************************************** 
!       
! Initialization of the FileHandle object
!
! Fortran IO unit is usually assigned automatically using open(neuwunit=...),
! However, a specific unit can be provided as an input argument for special
! cases, such as standard input, standard output and etc.
!
!******************************************************************************** 
subroutine init_file_handle(fname, self, funit)
  character(len=*), intent(in)  :: fname
  type(FileHandle), intent(out) :: self
  integer, optional, intent(in) :: funit
  
  self%fname = fname
  if (present(funit)) self%funit = funit
end subroutine init_file_handle


!******************************************************************************** 
!       
! FileHandle constructor
!
!******************************************************************************** 
function finit_file_handle(fname, funit) result(self)
  character(len=*), intent(in)  :: fname
  integer, optional, intent(in) :: funit
  type(FileHandle)              :: self

  call init_file_handle(fname, self, funit)
end function finit_file_handle


!******************************************************************************** 
!
! Open file using FileHandle
!
!******************************************************************************** 
subroutine open_file(self, access, action, status, form, position, recl, ierror, error_msg)
  class(FileHandle), intent(inout)                     :: self
  character(len=*), optional, intent(in)               :: access
  character(len=*), optional, intent(in)               :: action
  character(len=*), optional, intent(in)               :: status
  character(len=*), optional, intent(in)               :: form
  character(len=*), optional, intent(in)               :: position
  integer, optional, intent(in)                        :: recl
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: recl_
  integer                                              :: ierror_
  character(len=:), allocatable                        :: access_
  character(len=:), allocatable                        :: action_
  character(len=:), allocatable                        :: status_
  character(len=:), allocatable                        :: form_
  character(len=:), allocatable                        :: position_
  character(len=charlen)                               :: error_msg_

  if (self%funit==output_unit .or. self%funit==input_unit .or. self%funit==error_unit) then
    if (present(ierror)) ierror = ierror_
    return
  end if

!debug: nvidia
!call self%close()
  
  access_ = "sequential"
  action_ = "readwrite"
  status_ = "unknown"
  form_ = "formatted"
  position_ = "asis"
  recl_ = 0
  
  if (present(access)) access_ = access
  if (present(action)) action_ = action
  if (present(status)) status_ = status
  if (present(form)) form_ = form
  if (present(position)) position_ = position
  if (present(recl)) recl_ = recl
  
  if (trim(access_) == "direct") then
    open(newunit=self%funit, file=self%fname, access=access_, action=action_, form=form_, &
         status=status_, recl=recl_, iostat=ierror_, iomsg=error_msg_)
  else
    open(newunit=self%funit, file=self%fname, access=access_, action=action_, form=form_, &
         status=status_, position=position_, iostat=ierror_, iomsg=error_msg_)
  end if
  
  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine open_file


!******************************************************************************** 
!
! Close file 
!
!******************************************************************************** 
subroutine close_file(self, status, ierror, error_msg)
  class(FileHandle), intent(inout)                     :: self
  character(len=*), optional, intent(in)               :: status
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local         
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_

  ierror_ = 0
  
  if (.not. self%exists()) return
  
  if (present(status)) then
    if (self%is_open()) close(unit=self%funit, status=status, iostat=ierror_, iomsg=error_msg_)
  else 
    if (self%is_open()) close(unit=self%funit, iostat=ierror_, iomsg=error_msg_)
  end if 

  self%funit = unset

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine close_file


!******************************************************************************** 
!
! Delete file if exists
!
!******************************************************************************** 
subroutine delete_file(self)
  class(FileHandle), intent(inout) :: self
  !local
  character(len=:), allocatable    :: cmmd
  
  if (.not. self%exists()) return

  if (self%is_open()) then
    call self%close(status="delete")
  else 
    cmmd = "rm -f " // self%fname
    call execute_command_line(cmmd)
  end if 
end subroutine delete_file


!******************************************************************************** 
!
! Check whether the file is open
!
!******************************************************************************** 
logical function is_file_open(self)
  class(FileHandle), intent(in) :: self
  
  inquire(unit=self%funit, opened=is_file_open)
end function is_file_open


!******************************************************************************** 
!
! Check whether the file exists
!
!******************************************************************************** 
logical function file_exists(self)
  class(FileHandle), intent(in) :: self
  
  inquire(file=self%fname, exist=file_exists)
end function file_exists

end module file_handle