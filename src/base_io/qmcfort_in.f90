! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module qmcfort_in

use constants
use mpi, only: comm_world
use string
use qmcfort_io

implicit none

private
public add_input

interface add_input
  module procedure add_input_log, add_input_char
  module procedure add_input_int4, add_input_int4_1d_alloc
  module procedure add_input_int8, add_input_int8_1d_alloc
  module procedure add_input_real, add_input_real_1d_alloc
end interface add_input

contains

!******************************************************************************** 
!
!       Add input routines : reads variable with name var_name from the INPUT
!               fiel and stores it i value_
!
!******************************************************************************** 
subroutine add_input_log(var_name, value_, found)
  character(len=*), intent(in)   :: var_name
  logical, intent(inout)         :: value_
  logical, optional, intent(out) :: found
  !local variables
  character(len=charlen)         :: key, value_str
  integer                        :: ios, ios2, info
  logical                        :: found__

  if (comm_world%mpirank == 0) then
    if (.not. io%qmcfort_in%exists()) then
      found__ = .false.
      return
    end if

    call io%qmcfort_in%open(status="old", action="read")

    ios = 0
    found__ = .false.

    do while (ios == 0)
      read(io%qmcfort_in%funit,'(A)',iostat=ios) key
      
      if (ios == 0) then
        call split_tag(key, value_str, info)
        if (info == 0) cycle
        if (trim(lowercase(var_name)) == trim(lowercase(key))) then
          read(value_str,*,iostat=ios2) value_
          found__ = .true.
        end if
      end if
    end do

    call io%qmcfort_in%close()
  end if

  call comm_world%bcast(found__, 0)
  if (found__) call comm_world%bcast(value_, 0)
  if (present(found)) found = found__
end subroutine add_input_log

subroutine add_input_int4(var_name, value_, found)
  character(len=*), intent(in)   :: var_name
  integer(i4), intent(inout)     :: value_
  logical, optional, intent(out) :: found
  !local variables
  character(len=charlen)         :: key, value_str
  integer                        :: ios, ios2, info
  logical                        :: found__

  if (comm_world%mpirank == 0) then
    if (.not. io%qmcfort_in%exists()) then
      found__ = .false.
      return
    end if

    call io%qmcfort_in%open(status="old", action="read")

    ios = 0
    found__ = .false.

    do while (ios == 0)
      read(io%qmcfort_in%funit,'(A)',iostat=ios) key
      
      if (ios == 0) then
        call split_tag(key, value_str, info)
        if (info == 0) cycle
        if (trim(lowercase(var_name)) == trim(lowercase(key))) then
          read(value_str,*,iostat=ios2) value_
          found__ = .true.
        end if
      end if
    end do

    call io%qmcfort_in%close()
  end if

  call comm_world%bcast(found__, 0)
  if (found__) call comm_world%bcast(value_, 0)
  if (present(found)) found = found__
end subroutine add_input_int4

subroutine add_input_int4_1d(var_name, value_, found)
  character(len=*), intent(in)   :: var_name
  integer(i4), intent(inout)     :: value_(:)
  logical, optional, intent(out) :: found
  !local variables
  character(len=charlen)         :: key, value_str
  integer                        :: ios, ios2, info
  logical                        :: found__

  if (comm_world%mpirank == 0) then
    if (.not. io%qmcfort_in%exists()) then
      found__ = .false.
      return
    end if

    call io%qmcfort_in%open(status="old", action="read")

    ios = 0
    found__ = .false.

    do while (ios == 0)
      read(io%qmcfort_in%funit,'(A)',iostat=ios) key
      
      if (ios == 0) then
        call split_tag(key, value_str, info)
        if (info == 0) cycle
        if (trim(lowercase(var_name)) == trim(lowercase(key))) then
          read(value_str,*,iostat=ios2) value_
          found__ = .true.
        end if
      end if
    end do

    call io%qmcfort_in%close()
  end if

  call comm_world%bcast(found__, 0)
  if (found__) call comm_world%bcast(value_, 0)
  if (present(found)) found = found__
end subroutine add_input_int4_1d

subroutine add_input_int4_1d_alloc(var_name, value_, found)
  character(len=*), intent(in)            :: var_name
  integer(i4), allocatable, intent(inout) :: value_(:)
  logical, optional, intent(out)          :: found
  !local variables
  integer(i4), parameter                  :: const = 12345678
  integer                                 :: i, n
  integer(i4)                             :: temp(100)
  logical                                 :: found__

  temp = const
  call add_input_int4_1d(var_name, temp, found__)

  do i = 1, size(temp)
    if (temp(i) == const) then
      n = i - 1
      exit
    end if
  end do

  if (n > 0) then
    if (allocated(value_)) deallocate(value_)
    allocate(value_(n))
    value_(1:n) = temp(1:n)
  end if

  if (present(found)) found = found__
end subroutine add_input_int4_1d_alloc

subroutine add_input_int8(var_name, value_, found)
  character(len=*), intent(in)   :: var_name
  integer(i8), intent(inout)     :: value_
  logical, optional, intent(out) :: found
  !local variables
  character(len=charlen)         :: key, value_str
  integer                        :: ios, ios2, info
  logical                        :: found__

  if (comm_world%mpirank == 0) then
    if (.not. io%qmcfort_in%exists()) then
      found__ = .false.
      return
    end if

    call io%qmcfort_in%open(status="old", action="read")

    ios = 0
    found__ = .false.

    do while (ios == 0)
      read(io%qmcfort_in%funit,'(A)',iostat=ios) key
      
      if (ios == 0) then
        call split_tag(key, value_str, info)
        if (info == 0) cycle
        if (trim(lowercase(var_name)) == trim(lowercase(key))) then
          read(value_str,*,iostat=ios2) value_
          found__ = .true.
        end if
      end if
    end do

    call io%qmcfort_in%close()
  end if

  call comm_world%bcast(found__, 0)
  if (found__) call comm_world%bcast(value_, 0)
  if (present(found)) found = found__
end subroutine add_input_int8

subroutine add_input_int8_1d(var_name, value_, found)
  character(len=*), intent(in)   :: var_name
  integer(i8), intent(inout)     :: value_(:)
  logical, optional, intent(out) :: found
  !local variables
  character(len=charlen)         :: key, value_str
  integer                        :: ios, ios2, info
  logical                        :: found__

  if (comm_world%mpirank == 0) then
    if (.not. io%qmcfort_in%exists()) then
      found__ = .false.
      return
    end if

    call io%qmcfort_in%open(status="old", action="read")

    ios = 0
    found__ = .false.

    do while (ios == 0)
      read(io%qmcfort_in%funit,'(A)',iostat=ios) key
      
      if (ios == 0) then
        call split_tag(key, value_str, info)
        if (info == 0) cycle
        if (trim(lowercase(var_name)) == trim(lowercase(key))) then
          read(value_str,*,iostat=ios2) value_
          found__ = .true.
        end if
      end if
    end do

    call io%qmcfort_in%close()
  end if

  call comm_world%bcast(found__, 0)
  if (found__) call comm_world%bcast(value_, 0)
  if (present(found)) found = found__
end subroutine add_input_int8_1d

subroutine add_input_int8_1d_alloc(var_name, value_, found)
  character(len=*), intent(in)            :: var_name
  integer(i8), allocatable, intent(inout) :: value_(:)
  logical, optional, intent(out)          :: found
  !local variables
  integer(i8), parameter                  :: const = 12345678
  integer                                 :: i, n
  integer(i8)                             :: temp(100)
  logical                                 :: found__

  temp = const
  call add_input_int8_1d(var_name, temp, found__)

  do i = 1, size(temp)
    if (temp(i) == const) then
      n = i - 1
      exit
    end if
  end do

  if (n > 0) then
    if (allocated(value_)) deallocate(value_)
    allocate(value_(n))
    value_(1:n) = temp(1:n)
  end if

  if (present(found)) found = found__
end subroutine add_input_int8_1d_alloc

subroutine add_input_real(var_name, value_, found)
  character(len=*), intent(in)   :: var_name
  real(wp), intent(inout)        :: value_
  logical, optional, intent(out) :: found
  !local variables
  character(len=charlen)         :: key, value_str
  integer                        :: ios, ios2, info
  logical                        :: found__

  if (comm_world%mpirank == 0) then
    if (.not. io%qmcfort_in%exists()) then
      found__ = .false.
      return
    end if

    call io%qmcfort_in%open(status="old", action="read")

    ios = 0
    found__ = .false.

    do while (ios == 0)
      read(io%qmcfort_in%funit,'(A)',iostat=ios) key
      
      if (ios == 0) then
        call split_tag(key, value_str, info)
        if (info == 0) cycle
        if (trim(lowercase(var_name)) == trim(lowercase(key))) then
          read(value_str,*,iostat=ios2) value_
          found__ = .true.
        end if
      end if
    end do

    call io%qmcfort_in%close()
  end if

  call comm_world%bcast(found__, 0)
  if (found__) call comm_world%bcast(value_, 0)
  if (present(found)) found = found__
end subroutine add_input_real

subroutine add_input_real_1d(var_name, value_, found)
  character(len=*), intent(in)   :: var_name
  real(wp), intent(inout)        :: value_(:)
  logical, optional, intent(out) :: found
  !local variables
  character(len=charlen)         :: key, value_str
  integer                        :: ios, ios2, info
  logical                        :: found__

  if (comm_world%mpirank == 0) then
    if (.not. io%qmcfort_in%exists()) then
      found__ = .false.
      return
    end if

    call io%qmcfort_in%open(status="old", action="read")

    ios = 0
    found__ = .false.

    do while (ios == 0)
      read(io%qmcfort_in%funit,'(A)',iostat=ios) key
      
      if (ios == 0) then
        call split_tag(key, value_str, info)
        if (info == 0) cycle
        if (trim(lowercase(var_name)) == trim(lowercase(key))) then
          read(value_str,*,iostat=ios2) value_
          found__ = .true.
        end if
      end if
    end do

    call io%qmcfort_in%close()
  end if

  call comm_world%bcast(found__, 0)
  if (found__) call comm_world%bcast(value_, 0)
  if (present(found)) found = found__
end subroutine add_input_real_1d

subroutine add_input_real_1d_alloc(var_name, value_, found)
  character(len=*), intent(in)         :: var_name
  real(wp), allocatable, intent(inout) :: value_(:)
  logical, optional, intent(out)       :: found
  !local variables
  logical                              :: found__
  integer                              :: i, n
  real(wp), parameter                  :: const = -1234.5678_wp
  real(wp)                             :: temp(20)

  temp = const
  call add_input_real_1d(var_name, temp, found__)

  do i = 1, size(temp)
    if (temp(i) == const) then
      n = i - 1
      exit
    end if
  end do

  if (n > 0) then
    if (allocated(value_)) deallocate(value_)
    allocate(value_(n))
    value_(1:n) = temp(1:n)
  end if

  if (present(found)) found = found__
end subroutine add_input_real_1d_alloc

subroutine add_input_char(var_name, value_, found)
  character(len=*), intent(in)    :: var_name
  character(len=*), intent(inout) :: value_
  logical, optional, intent(out)  :: found
  !local variables
  character(len=charlen)          :: key, value_str
  integer                         :: ios, ios2, info
  logical                         :: found__

  if (comm_world%mpirank == 0) then
    if (.not. io%qmcfort_in%exists()) then
      found__ = .false.
      return
    end if

    call io%qmcfort_in%open(status="old", action="read")

    ios = 0
    found__ = .false.

    do while (ios == 0)
      read(io%qmcfort_in%funit,'(A)',iostat=ios) key
      
      if (ios == 0) then
        call split_tag(key, value_str, info)
        if (info == 0) cycle
        if (trim(lowercase(var_name)) == trim(lowercase(key))) then
          read(value_str,*,iostat=ios2) value_
!debug:
!          value_ = lowercase(value_)
          found__ = .true.
        end if
      end if
    end do

    call io%qmcfort_in%close()
  end if

  call comm_world%bcast(found__, 0)
  if (found__) call comm_world%bcast(value_, 0)
  if (present(found)) found = found__
end subroutine add_input_char

end module qmcfort_in 
