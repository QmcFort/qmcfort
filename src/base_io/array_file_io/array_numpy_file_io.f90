! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0


module array_numpy_file_io

use iso_c_binding, only: c_short

use constants
use string, only: str, rest_str, little_endian
use meta_data_utils, only: type_kind_label
use file_handle, only: FileHandle
use array_file_io, only: ArrayFileIO

implicit none

private
public :: ArrayNumpyFileIO

!Magic string at beginning of numpy file
integer, parameter :: magic_len = 8
character(len=magic_len), parameter :: magic = char(147) // "NUMPY" // char(1) // char(0)

type, extends(ArrayFileIO) :: ArrayNumpyFileIO
contains
  procedure, private :: numpy_byteorder
  procedure, private :: header_length
  procedure, private :: pad_header
  procedure, private :: get_header
  procedure, private :: write_header
  procedure, private :: read_header
  procedure, private :: read_shape_from_header

  !deferred procedures 
  procedure          :: write_int, write_int_1d, write_int_2d, write_int_3d, write_int_4d, write_int_5d
  procedure          :: write_real, write_real_1d, write_real_2d, write_real_3d, write_real_4d, write_real_5d
  procedure          :: write_cmplx, write_cmplx_1d, write_cmplx_2d, write_cmplx_3d, write_cmplx_4d, write_cmplx_5d

  procedure          :: read_int, read_int_1d, read_int_2d, read_int_3d, read_int_4d, read_int_5d
  procedure          :: read_real, read_real_1d, read_real_2d, read_real_3d, read_real_4d, read_real_5d
  procedure          :: read_cmplx, read_cmplx_1d, read_cmplx_2d, read_cmplx_3d, read_cmplx_4d, read_cmplx_5d

  procedure          :: save_int, save_int_1d, save_int_2d, save_int_3d, save_int_4d, save_int_5d
  procedure          :: save_real, save_real_1d, save_real_2d, save_real_3d, save_real_4d, save_real_5d
  procedure          :: save_cmplx, save_cmplx_1d, save_cmplx_2d, save_cmplx_3d, save_cmplx_4d, save_cmplx_5d

  procedure          :: load_int, load_int_1d, load_int_2d, load_int_3d, load_int_4d, load_int_5d
  procedure          :: load_real, load_real_1d, load_real_2d, load_real_3d, load_real_4d, load_real_5d
  procedure          :: load_cmplx, load_cmplx_1d, load_cmplx_2d, load_cmplx_3d, load_cmplx_4d, load_cmplx_5d
end type ArrayNumpyFileIO

contains

!******************************************************************************** 
!
! Returns "<" (little endian) or ">" (big endian)
!
!******************************************************************************** 
function numpy_byteorder(self) result(byteorder)
  class(ArrayNumpyFileIO), intent(in) :: self
  character                           :: byteorder

  if (little_endian()) then
    byteorder = "<"
  else
    byteorder = ">"
  end if
end function numpy_byteorder


!******************************************************************************** 
!
! Determine eligible length of header such that the header 
! and the extra bytes amount to a multiple of 64
!
!******************************************************************************** 
function header_length(self, header) result(length)
  class(ArrayNumpyFileIO), intent(in) :: self
  character(len=*), intent(in)        :: header
  integer(c_short)                    :: length
  !local
  integer(c_short), parameter         :: bytes_before_header = len(magic) + 2
  integer(c_short)                    :: multiple
  
  if (.not. little_endian()) stop "Big endian systems not supported"

  multiple = ((len_trim(header) + bytes_before_header) / 64 + 1)
  length = 64*multiple - bytes_before_header
end function header_length


!******************************************************************************** 
!
! Create a header for numpy file.
! The full header length with magic string should be multiple of 64.
!
!******************************************************************************** 
function get_header(self, shape_, dtype) result(header)
  class(ArrayNumpyFileIO), intent(in) :: self
  integer, intent(in)                 :: shape_(:)
  character(len=*), intent(in)        :: dtype
  character(len=:), allocatable       :: header
  !local
  integer                             :: i

  header = "{" &
        // "'descr': '" // self%numpy_byteorder() // dtype // "', " &
        // "'fortran_order': True, " &
        // "'shape': ("

  do i = 1, size(shape_)
    header = header // str(shape_(i)) // ","
  end do

  header = header // "), }"

  call self%pad_header(header)
end function get_header


!******************************************************************************** 
!
! Pad header with spaces to have a length that is multiple of 64 
!
!******************************************************************************** 
subroutine pad_header(self, header)
  class(ArrayNumpyFileIO), intent(in)          :: self
  character(len=:), allocatable, intent(inout) :: header
  !local
  integer                                      :: num_spaces 
  
  num_spaces = self%header_length(header) - len(header) - 1
  header = header // repeat(" ", num_spaces) // new_line("n")
end subroutine pad_header


!******************************************************************************** 
!
! Write header of the numpy file 
!
!******************************************************************************** 
subroutine write_header(self, funit, header, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)           :: self
  integer, intent(in)                           :: funit
  character(len=*), intent(in)                  :: header
  integer, optional, intent(out)                :: ierror
  character(len=charlen), optional, intent(out) :: error_msg
  !local
  integer                                       :: ierror_
  integer(c_short)                              :: len_header
  character(len=charlen)                        :: error_msg_
  
  len_header = len(header)

  write(funit, iostat=ierror_, iomsg=error_msg_) magic, len_header
  write(funit, iostat=ierror_, iomsg=error_msg_) header

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_header


!******************************************************************************** 
!
! Read header of the numpy file
!
!******************************************************************************** 
subroutine read_header(self, funit, header, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)           :: self
  integer, intent(in)                           :: funit
  character(len=:), allocatable, intent(out)    :: header
  integer, optional, intent(out)                :: ierror
  character(len=charlen), optional, intent(out) :: error_msg
  !local
  integer                                       :: ierror_
  integer(c_short)                              :: length
  character(len=magic_len)                      :: magic_str
  character(len=charlen)                        :: error_msg_
  
  read(funit, iostat=ierror_, iomsg=error_msg_) magic_str, length
  allocate(character(len=length) :: header)
  read(funit, iostat=ierror_, iomsg=error_msg_) header

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_header


!******************************************************************************** 
!
! Read header of the numpy file
!
!******************************************************************************** 
subroutine read_shape_from_header(self, header, shape_)
  class(ArrayNumpyFileIO), intent(in) :: self
  character(len=*), intent(in)        :: header
  integer, intent(out)                :: shape_(:)
  !local
  integer                             :: i, pos
  character(len=:), allocatable       :: temp

  temp = rest_str("'shape': (", header, ",")
  pos = scan(temp, ")")
  temp = temp(1:pos-1)
  
  do i = 1, size(shape_)
    read(temp, *) shape_(i)
    pos = scan(temp, ",")
    temp = temp(pos+1:)
  end do
end subroutine read_shape_from_header


!******************************************************************************** 
!
! Write an array to the numpy file
!
!******************************************************************************** 
subroutine write_int(self, funit, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self
  integer, intent(in)                                  :: funit
  integer, intent(in)                                  :: array
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: ierror_
  integer                                              :: temp(1)
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: header
  
  temp(1) = array
  header = self%get_header(shape(temp), type_kind_label(temp))

  call self%write_header(funit, header, ierror=ierror_, error_msg=error_msg_)
  write(funit, iostat=ierror_, iomsg=error_msg_) temp

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_int

subroutine write_int_1d(self, funit, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self
  integer, intent(in)                                  :: funit
  integer, intent(in)                                  :: array(:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: header
  
  header = self%get_header(shape(array), type_kind_label(array))

  call self%write_header(funit, header, ierror=ierror_, error_msg=error_msg_)
  write(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_int_1d

subroutine write_int_2d(self, funit, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self
  integer, intent(in)                                  :: funit
  integer, intent(in)                                  :: array(:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: header
  
  header = self%get_header(shape(array), type_kind_label(array))

  call self%write_header(funit, header, ierror=ierror_, error_msg=error_msg_)
  write(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_int_2d

subroutine write_int_3d(self, funit, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self
  integer, intent(in)                                  :: funit
  integer, intent(in)                                  :: array(:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: header
  
  header = self%get_header(shape(array), type_kind_label(array))

  call self%write_header(funit, header, ierror=ierror_, error_msg=error_msg_)
  write(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_int_3d

subroutine write_int_4d(self, funit, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self
  integer, intent(in)                                  :: funit
  integer, intent(in)                                  :: array(:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: header
  
  header = self%get_header(shape(array), type_kind_label(array))

  call self%write_header(funit, header, ierror=ierror_, error_msg=error_msg_)
  write(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_int_4d

subroutine write_int_5d(self, funit, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self
  integer, intent(in)                                  :: funit
  integer, intent(in)                                  :: array(:,:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: header
  
  header = self%get_header(shape(array), type_kind_label(array))

  call self%write_header(funit, header, ierror=ierror_, error_msg=error_msg_)
  write(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_int_5d

subroutine write_real(self, funit, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self
  integer, intent(in)                                  :: funit
  real(wp), intent(in)                                 :: array
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: ierror_
  real(wp)                                             :: temp(1)
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: header
  
  header = self%get_header(shape(temp), type_kind_label(temp))
  temp(1) = array

  call self%write_header(funit, header, ierror=ierror_, error_msg=error_msg_)
  write(funit, iostat=ierror_, iomsg=error_msg_) temp

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_real

subroutine write_real_1d(self, funit, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self
  integer, intent(in)                                  :: funit
  real(wp), intent(in)                                 :: array(:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: header
  
  header = self%get_header(shape(array), type_kind_label(array))

  call self%write_header(funit, header, ierror=ierror_, error_msg=error_msg_)
  write(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_real_1d

subroutine write_real_2d(self, funit, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self
  integer, intent(in)                                  :: funit
  real(wp), intent(in)                                 :: array(:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: header
  
  header = self%get_header(shape(array), type_kind_label(array))

  call self%write_header(funit, header, ierror=ierror_, error_msg=error_msg_)
  write(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_real_2d

subroutine write_real_3d(self, funit, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self
  integer, intent(in)                                  :: funit
  real(wp), intent(in)                                 :: array(:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: header
  
  header = self%get_header(shape(array), type_kind_label(array))

  call self%write_header(funit, header, ierror=ierror_, error_msg=error_msg_)
  write(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_real_3d

subroutine write_real_4d(self, funit, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self
  integer, intent(in)                                  :: funit
  real(wp), intent(in)                                 :: array(:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: header
  
  header = self%get_header(shape(array), type_kind_label(array))

  call self%write_header(funit, header, ierror=ierror_, error_msg=error_msg_)
  write(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_real_4d

subroutine write_real_5d(self, funit, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self
  integer, intent(in)                                  :: funit
  real(wp), intent(in)                                 :: array(:,:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: header
  
  header = self%get_header(shape(array), type_kind_label(array))

  call self%write_header(funit, header, ierror=ierror_, error_msg=error_msg_)
  write(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_real_5d

subroutine write_cmplx(self, funit, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self
  integer, intent(in)                                  :: funit
  complex(wp), intent(in)                              :: array
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: ierror_
  complex(wp)                                          :: temp(1)
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: header
  
  header = self%get_header(shape(temp), type_kind_label(temp))
  temp(1) = array

  call self%write_header(funit, header, ierror=ierror_, error_msg=error_msg_)
  write(funit, iostat=ierror_, iomsg=error_msg_) temp

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_cmplx

subroutine write_cmplx_1d(self, funit, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self
  integer, intent(in)                                  :: funit
  complex(wp), intent(in)                              :: array(:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: header
  
  header = self%get_header(shape(array), type_kind_label(array))

  call self%write_header(funit, header, ierror=ierror_, error_msg=error_msg_)
  write(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_cmplx_1d

subroutine write_cmplx_2d(self, funit, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self
  integer, intent(in)                                  :: funit
  complex(wp), intent(in)                              :: array(:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: header
  
  header = self%get_header(shape(array), type_kind_label(array))

  call self%write_header(funit, header, ierror=ierror_, error_msg=error_msg_)
  write(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_cmplx_2d

subroutine write_cmplx_3d(self, funit, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self
  integer, intent(in)                                  :: funit
  complex(wp), intent(in)                              :: array(:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: header
  
  header = self%get_header(shape(array), type_kind_label(array))

  call self%write_header(funit, header, ierror=ierror_, error_msg=error_msg_)
  write(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_cmplx_3d

subroutine write_cmplx_4d(self, funit, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self
  integer, intent(in)                                  :: funit
  complex(wp), intent(in)                              :: array(:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: header
  
  header = self%get_header(shape(array), type_kind_label(array))

  call self%write_header(funit, header, ierror=ierror_, error_msg=error_msg_)
  write(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_cmplx_4d

subroutine write_cmplx_5d(self, funit, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self
  integer, intent(in)                                  :: funit
  complex(wp), intent(in)                              :: array(:,:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: header
  
  header = self%get_header(shape(array), type_kind_label(array))

  call self%write_header(funit, header, ierror=ierror_, error_msg=error_msg_)
  write(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_cmplx_5d


!********************************************************************************
!
! Read an array from a numpy file
!
!******************************************************************************** 
subroutine read_int(self, funit, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  integer, intent(in)                                  :: funit
  integer, intent(out)                                 :: array
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 1
  integer                                              :: shape_(ndim), ierror_
  integer, allocatable                                 :: temp(:)
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: header

  call self%read_header(funit, header, ierror=ierror_, error_msg=error_msg_)
  call self%read_shape_from_header(header, shape_)
  allocate(temp(shape_(1)), stat=ierror_, errmsg=error_msg_)
  read(funit, iostat=ierror_, iomsg=error_msg_) temp
  array = temp(1)

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_int

subroutine read_int_1d(self, funit, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  integer, intent(in)                                  :: funit
  integer, allocatable, intent(out)                    :: array(:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 1
  integer                                              :: shape_(ndim), ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: header

  call self%read_header(funit, header, ierror=ierror_, error_msg=error_msg_)
  call self%read_shape_from_header(header, shape_)
  allocate(array(shape_(1)), stat=ierror_, errmsg=error_msg_)
  read(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_int_1d

subroutine read_int_2d(self, funit, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  integer, intent(in)                                  :: funit
  integer, allocatable, intent(out)                    :: array(:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 2
  integer                                              :: shape_(ndim), ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: header

  call self%read_header(funit, header, ierror=ierror_, error_msg=error_msg_)
  call self%read_shape_from_header(header, shape_)
  allocate(array(shape_(1), shape_(2)), stat=ierror_, errmsg=error_msg_)
  read(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_int_2d

subroutine read_int_3d(self, funit, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  integer, intent(in)                                  :: funit
  integer, allocatable, intent(out)                    :: array(:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 3
  integer                                              :: shape_(ndim), ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: header

  call self%read_header(funit, header, ierror=ierror_, error_msg=error_msg_)
  call self%read_shape_from_header(header, shape_)
  allocate(array(shape_(1), shape_(2), shape_(3)), stat=ierror_, errmsg=error_msg_)
  read(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_int_3d

subroutine read_int_4d(self, funit, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  integer, intent(in)                                  :: funit
  integer, allocatable, intent(out)                    :: array(:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 4
  integer                                              :: shape_(ndim), ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: header

  call self%read_header(funit, header, ierror=ierror_, error_msg=error_msg_)
  call self%read_shape_from_header(header, shape_)
  allocate(array(shape_(1), shape_(2), shape_(3), shape_(4)), stat=ierror_, errmsg=error_msg_)
  read(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_int_4d

subroutine read_int_5d(self, funit, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  integer, intent(in)                                  :: funit
  integer, allocatable, intent(out)                    :: array(:,:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 5
  integer                                              :: shape_(ndim), ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: header

  call self%read_header(funit, header, ierror=ierror_, error_msg=error_msg_)
  call self%read_shape_from_header(header, shape_)
  allocate(array(shape_(1), shape_(2), shape_(3), shape_(4), shape_(5)), stat=ierror_, errmsg=error_msg_)
  read(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_int_5d

subroutine read_real(self, funit, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  integer, intent(in)                                  :: funit
  real(wp), intent(out)                                :: array
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 1
  integer                                              :: shape_(ndim), ierror_
  real(wp), allocatable                                :: temp(:)
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: header

  call self%read_header(funit, header, ierror=ierror_, error_msg=error_msg_)
  call self%read_shape_from_header(header, shape_)
  allocate(temp(shape_(1)), stat=ierror_, errmsg=error_msg_)
  read(funit, iostat=ierror_, iomsg=error_msg_) temp
  array = temp(1)

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_real

subroutine read_real_1d(self, funit, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  integer, intent(in)                                  :: funit
  real(wp), allocatable, intent(out)                   :: array(:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 1
  integer                                              :: shape_(ndim), ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: header

  call self%read_header(funit, header, ierror=ierror_, error_msg=error_msg_)
  call self%read_shape_from_header(header, shape_)
  allocate(array(shape_(1)), stat=ierror_, errmsg=error_msg_)
  read(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_real_1d

subroutine read_real_2d(self, funit, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  integer, intent(in)                                  :: funit
  real(wp), allocatable, intent(out)                   :: array(:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 2
  integer                                              :: shape_(ndim), ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: header

  call self%read_header(funit, header, ierror=ierror_, error_msg=error_msg_)
  call self%read_shape_from_header(header, shape_)
  allocate(array(shape_(1), shape_(2)), stat=ierror_, errmsg=error_msg_)
  read(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_real_2d

subroutine read_real_3d(self, funit, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  integer, intent(in)                                  :: funit
  real(wp), allocatable, intent(out)                   :: array(:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 3
  integer                                              :: shape_(ndim), ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: header

  call self%read_header(funit, header, ierror=ierror_, error_msg=error_msg_)
  call self%read_shape_from_header(header, shape_)
  allocate(array(shape_(1), shape_(2), shape_(3)), stat=ierror_, errmsg=error_msg_)
  read(funit, iostat=ierror_, iomsg=error_msg_) array
  
  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_real_3d

subroutine read_real_4d(self, funit, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  integer, intent(in)                                  :: funit
  real(wp), allocatable, intent(out)                   :: array(:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 4
  integer                                              :: shape_(ndim), ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: header

  call self%read_header(funit, header, ierror=ierror_, error_msg=error_msg_)
  call self%read_shape_from_header(header, shape_)
  allocate(array(shape_(1), shape_(2), shape_(3), shape_(4)), stat=ierror_, errmsg=error_msg_)
  read(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_real_4d

subroutine read_real_5d(self, funit, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  integer, intent(in)                                  :: funit
  real(wp), allocatable, intent(out)                   :: array(:,:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 5
  integer                                              :: shape_(ndim), ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: header

  call self%read_header(funit, header, ierror=ierror_, error_msg=error_msg_)
  call self%read_shape_from_header(header, shape_)
  allocate(array(shape_(1), shape_(2), shape_(3), shape_(4), shape_(5)), stat=ierror_, errmsg=error_msg_)
  read(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_real_5d

subroutine read_cmplx(self, funit, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  integer, intent(in)                                  :: funit
  complex(wp), intent(out)                             :: array
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 1
  integer                                              :: shape_(ndim), ierror_
  complex(wp), allocatable                             :: temp(:)
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: header

  call self%read_header(funit, header, ierror=ierror_, error_msg=error_msg_)
  call self%read_shape_from_header(header, shape_)
  allocate(temp(shape_(1)), stat=ierror_, errmsg=error_msg_)
  read(funit, iostat=ierror_, iomsg=error_msg_) temp
  array = temp(1)

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_cmplx

subroutine read_cmplx_1d(self, funit, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  integer, intent(in)                                  :: funit
  complex(wp), allocatable, intent(out)                :: array(:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 1
  integer                                              :: shape_(ndim), ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: header

  call self%read_header(funit, header, ierror=ierror_, error_msg=error_msg_)
  call self%read_shape_from_header(header, shape_)
  allocate(array(shape_(1)), stat=ierror_, errmsg=error_msg_)
  read(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_cmplx_1d

subroutine read_cmplx_2d(self, funit, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  integer, intent(in)                                  :: funit
  complex(wp), allocatable, intent(out)                :: array(:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 2
  integer                                              :: shape_(ndim), ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: header

  call self%read_header(funit, header, ierror=ierror_, error_msg=error_msg_)
  call self%read_shape_from_header(header, shape_)
  allocate(array(shape_(1), shape_(2)), stat=ierror_, errmsg=error_msg_)
  read(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_cmplx_2d

subroutine read_cmplx_3d(self, funit, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  integer, intent(in)                                  :: funit
  complex(wp), allocatable, intent(out)                :: array(:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 3
  integer                                              :: shape_(ndim), ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: header

  call self%read_header(funit, header, ierror=ierror_, error_msg=error_msg_)
  call self%read_shape_from_header(header, shape_)
  allocate(array(shape_(1), shape_(2), shape_(3)), stat=ierror_, errmsg=error_msg_)
  read(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_cmplx_3d

subroutine read_cmplx_4d(self, funit, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  integer, intent(in)                                  :: funit
  complex(wp), allocatable, intent(out)                :: array(:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 4
  integer                                              :: shape_(ndim), ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: header

  call self%read_header(funit, header, ierror=ierror_, error_msg=error_msg_)
  call self%read_shape_from_header(header, shape_)
  allocate(array(shape_(1), shape_(2), shape_(3), shape_(4)), stat=ierror_, errmsg=error_msg_)
  read(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_cmplx_4d

subroutine read_cmplx_5d(self, funit, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  integer, intent(in)                                  :: funit
  complex(wp), allocatable, intent(out)                :: array(:,:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 5
  integer                                              :: shape_(ndim), ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: header

  call self%read_header(funit, header, ierror=ierror_, error_msg=error_msg_)
  call self%read_shape_from_header(header, shape_)
  allocate(array(shape_(1), shape_(2), shape_(3), shape_(4), shape_(5)), stat=ierror_, errmsg=error_msg_)
  read(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_cmplx_5d


!********************************************************************************
!
! Save an array to the numpy file
!
!*******************************************************************************
subroutine save_int(self, fname, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  character(len=*), intent(in)                         :: fname
  integer, intent(in)                                  :: array
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="write", status="replace", access="stream", form="unformatted")
  call self%write(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine save_int

subroutine save_int_1d(self, fname, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  character(len=*), intent(in)                         :: fname
  integer, intent(in)                                  :: array(:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="write", status="replace", access="stream", form="unformatted")
  call self%write(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine save_int_1d

subroutine save_int_2d(self, fname, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  character(len=*), intent(in)                         :: fname
  integer, intent(in)                                  :: array(:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="write", status="replace", access="stream", form="unformatted")
  call self%write(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine save_int_2d

subroutine save_int_3d(self, fname, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  character(len=*), intent(in)                         :: fname
  integer, intent(in)                                  :: array(:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="write", status="replace", access="stream", form="unformatted")
  call self%write(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine save_int_3d

subroutine save_int_4d(self, fname, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  character(len=*), intent(in)                         :: fname
  integer, intent(in)                                  :: array(:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="write", status="replace", access="stream", form="unformatted")
  call self%write(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine save_int_4d

subroutine save_int_5d(self, fname, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  character(len=*), intent(in)                         :: fname
  integer, intent(in)                                  :: array(:,:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="write", status="replace", access="stream", form="unformatted")
  call self%write(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine save_int_5d

subroutine save_real(self, fname, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  character(len=*), intent(in)                         :: fname
  real(wp), intent(in)                                 :: array
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="write", status="replace", access="stream", form="unformatted")
  call self%write(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine save_real

subroutine save_real_1d(self, fname, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  character(len=*), intent(in)                         :: fname
  real(wp), intent(in)                                 :: array(:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="write", status="replace", access="stream", form="unformatted")
  call self%write(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine save_real_1d

subroutine save_real_2d(self, fname, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  character(len=*), intent(in)                         :: fname
  real(wp), intent(in)                                 :: array(:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="write", status="replace", access="stream", form="unformatted")
  call self%write(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine save_real_2d

subroutine save_real_3d(self, fname, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  character(len=*), intent(in)                         :: fname
  real(wp), intent(in)                                 :: array(:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="write", status="replace", access="stream", form="unformatted")
  call self%write(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine save_real_3d

subroutine save_real_4d(self, fname, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  character(len=*), intent(in)                         :: fname
  real(wp), intent(in)                                 :: array(:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="write", status="replace", access="stream", form="unformatted")
  call self%write(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine save_real_4d

subroutine save_real_5d(self, fname, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  character(len=*), intent(in)                         :: fname
  real(wp), intent(in)                                 :: array(:,:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="write", status="replace", access="stream", form="unformatted")
  call self%write(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine save_real_5d

subroutine save_cmplx(self, fname, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  character(len=*), intent(in)                         :: fname
  complex(wp), intent(in)                              :: array
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="write", status="replace", access="stream", form="unformatted")
  call self%write(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine save_cmplx

subroutine save_cmplx_1d(self, fname, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  character(len=*), intent(in)                         :: fname
  complex(wp), intent(in)                              :: array(:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="write", status="replace", access="stream", form="unformatted")
  call self%write(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine save_cmplx_1d

subroutine save_cmplx_2d(self, fname, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  character(len=*), intent(in)                         :: fname
  complex(wp), intent(in)                              :: array(:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="write", status="replace", access="stream", form="unformatted")
  call self%write(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine save_cmplx_2d

subroutine save_cmplx_3d(self, fname, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  character(len=*), intent(in)                         :: fname
  complex(wp), intent(in)                              :: array(:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="write", status="replace", access="stream", form="unformatted")
  call self%write(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine save_cmplx_3d

subroutine save_cmplx_4d(self, fname, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  character(len=*), intent(in)                         :: fname
  complex(wp), intent(in)                              :: array(:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="write", status="replace", access="stream", form="unformatted")
  call self%write(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine save_cmplx_4d

subroutine save_cmplx_5d(self, fname, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  character(len=*), intent(in)                         :: fname
  complex(wp), intent(in)                              :: array(:,:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="write", status="replace", access="stream", form="unformatted")
  call self%write(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine save_cmplx_5d


!********************************************************************************
!
! Load an array from the file in txt format.
!
! Binary files are formatted as follows:
!    shape(array), array
!
!*******************************************************************************
subroutine load_int(self, fname, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  character(len=*), intent(in)                         :: fname
  integer, intent(out)                                 :: array
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="read", status="old", access="stream", form="unformatted")
  call self%read(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine load_int

subroutine load_int_1d(self, fname, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  character(len=*), intent(in)                         :: fname
  integer, allocatable, intent(out)                    :: array(:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="read", status="old", access="stream", form="unformatted")
  call self%read(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine load_int_1d

subroutine load_int_2d(self, fname, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  character(len=*), intent(in)                         :: fname
  integer, allocatable, intent(out)                    :: array(:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="read", status="old", access="stream", form="unformatted")
  call self%read(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine load_int_2d

subroutine load_int_3d(self, fname, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  character(len=*), intent(in)                         :: fname
  integer, allocatable, intent(out)                    :: array(:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="read", status="old", access="stream", form="unformatted")
  call self%read(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine load_int_3d

subroutine load_int_4d(self, fname, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  character(len=*), intent(in)                         :: fname
  integer, allocatable, intent(out)                    :: array(:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="read", status="old", access="stream", form="unformatted")
  call self%read(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine load_int_4d

subroutine load_int_5d(self, fname, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  character(len=*), intent(in)                         :: fname
  integer, allocatable, intent(out)                    :: array(:,:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="read", status="old", access="stream", form="unformatted")
  call self%read(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine load_int_5d

subroutine load_real(self, fname, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  character(len=*), intent(in)                         :: fname
  real(wp), intent(out)                                :: array
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="read", status="old", access="stream", form="unformatted")
  call self%read(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine load_real

subroutine load_real_1d(self, fname, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  character(len=*), intent(in)                         :: fname
  real(wp), allocatable, intent(out)                   :: array(:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="read", status="old", access="stream", form="unformatted")
  call self%read(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine load_real_1d

subroutine load_real_2d(self, fname, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  character(len=*), intent(in)                         :: fname
  real(wp), allocatable, intent(out)                   :: array(:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="read", status="old", access="stream", form="unformatted")
  call self%read(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine load_real_2d

subroutine load_real_3d(self, fname, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  character(len=*), intent(in)                         :: fname
  real(wp), allocatable, intent(out)                   :: array(:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="read", status="old", access="stream", form="unformatted")
  call self%read(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine load_real_3d

subroutine load_real_4d(self, fname, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  character(len=*), intent(in)                         :: fname
  real(wp), allocatable, intent(out)                   :: array(:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="read", status="old", access="stream", form="unformatted")
  call self%read(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine load_real_4d

subroutine load_real_5d(self, fname, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  character(len=*), intent(in)                         :: fname
  real(wp), allocatable, intent(out)                   :: array(:,:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="read", status="old", access="stream", form="unformatted")
  call self%read(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine load_real_5d

subroutine load_cmplx(self, fname, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  character(len=*), intent(in)                         :: fname
  complex(wp), intent(out)                             :: array
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="read", status="old", access="stream", form="unformatted")
  call self%read(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine load_cmplx

subroutine load_cmplx_1d(self, fname, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  character(len=*), intent(in)                         :: fname
  complex(wp), allocatable, intent(out)                :: array(:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="read", status="old", access="stream", form="unformatted")
  call self%read(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine load_cmplx_1d

subroutine load_cmplx_2d(self, fname, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  character(len=*), intent(in)                         :: fname
  complex(wp), allocatable, intent(out)                :: array(:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="read", status="old", access="stream", form="unformatted")
  call self%read(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine load_cmplx_2d

subroutine load_cmplx_3d(self, fname, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  character(len=*), intent(in)                         :: fname
  complex(wp), allocatable, intent(out)                :: array(:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="read", status="old", access="stream", form="unformatted")
  call self%read(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine load_cmplx_3d

subroutine load_cmplx_4d(self, fname, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  character(len=*), intent(in)                         :: fname
  complex(wp), allocatable, intent(out)                :: array(:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="read", status="old", access="stream", form="unformatted")
  call self%read(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine load_cmplx_4d

subroutine load_cmplx_5d(self, fname, array, ierror, error_msg)
  class(ArrayNumpyFileIO), intent(in)                  :: self 
  character(len=*), intent(in)                         :: fname
  complex(wp), allocatable, intent(out)                :: array(:,:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="read", status="old", access="stream", form="unformatted")
  call self%read(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine load_cmplx_5d

end module array_numpy_file_io