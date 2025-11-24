! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module array_bin_file_io

use constants
use file_handle, only: FileHandle
use array_file_io, only: ArrayFileIO

implicit none

private
public :: ArrayBinFileIO

type, extends(ArrayFileIO) :: ArrayBinFileIO
contains
  !deferred procedures 
  procedure :: write_int, write_int_1d, write_int_2d, write_int_3d, write_int_4d, write_int_5d
  procedure :: write_real, write_real_1d, write_real_2d, write_real_3d, write_real_4d, write_real_5d
  procedure :: write_cmplx, write_cmplx_1d, write_cmplx_2d, write_cmplx_3d, write_cmplx_4d, write_cmplx_5d

  procedure :: read_int, read_int_1d, read_int_2d, read_int_3d, read_int_4d, read_int_5d
  procedure :: read_real, read_real_1d, read_real_2d, read_real_3d, read_real_4d, read_real_5d
  procedure :: read_cmplx, read_cmplx_1d, read_cmplx_2d, read_cmplx_3d, read_cmplx_4d, read_cmplx_5d

  procedure :: save_int, save_int_1d, save_int_2d, save_int_3d, save_int_4d, save_int_5d
  procedure :: save_real, save_real_1d, save_real_2d, save_real_3d, save_real_4d, save_real_5d
  procedure :: save_cmplx, save_cmplx_1d, save_cmplx_2d, save_cmplx_3d, save_cmplx_4d, save_cmplx_5d

  procedure :: load_int, load_int_1d, load_int_2d, load_int_3d, load_int_4d, load_int_5d
  procedure :: load_real, load_real_1d, load_real_2d, load_real_3d, load_real_4d, load_real_5d
  procedure :: load_cmplx, load_cmplx_1d, load_cmplx_2d, load_cmplx_3d, load_cmplx_4d, load_cmplx_5d
end type ArrayBinFileIO

contains

!********************************************************************************
!
! Write an array in bin format.
!
! Binary files are formatted as follows:
!    shape(array), array
!
!*******************************************************************************
subroutine write_int(self, funit, array, ierror, error_msg)
  class(ArrayBinFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  integer, intent(in)                                  :: array
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: ierror_
  integer                                              :: temp(1)
  character(len=charlen)                               :: error_msg_

  temp(1) = array

  write(funit, iostat=ierror_, iomsg=error_msg_) shape(temp)
  write(funit, iostat=ierror_, iomsg=error_msg_) temp

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_int

subroutine write_int_1d(self, funit, array, ierror, error_msg)
  class(ArrayBinFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  integer, intent(in)                                  :: array(:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_

  write(funit, iostat=ierror_, iomsg=error_msg_) shape(array)
  write(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_int_1d

subroutine write_int_2d(self, funit, array, ierror, error_msg)
  class(ArrayBinFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  integer, intent(in)                                  :: array(:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_

  write(funit, iostat=ierror_, iomsg=error_msg_) shape(array)
  write(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_int_2d

subroutine write_int_3d(self, funit, array, ierror, error_msg)
  class(ArrayBinFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  integer, intent(in)                                  :: array(:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_

  write(funit, iostat=ierror_, iomsg=error_msg_) shape(array)
  write(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_int_3d

subroutine write_int_4d(self, funit, array, ierror, error_msg)
  class(ArrayBinFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  integer, intent(in)                                  :: array(:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_

  write(funit, iostat=ierror_, iomsg=error_msg_) shape(array)
  write(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_int_4d

subroutine write_int_5d(self, funit, array, ierror, error_msg)
  class(ArrayBinFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  integer, intent(in)                                  :: array(:,:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_

  write(funit, iostat=ierror_, iomsg=error_msg_) shape(array)

  write(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_int_5d

subroutine write_real(self, funit, array, ierror, error_msg)
  class(ArrayBinFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  real(wp), intent(in)                                 :: array
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: ierror_
  integer                                              :: temp(1)
  character(len=charlen)                               :: error_msg_

  temp(1) = array

  write(funit, iostat=ierror_, iomsg=error_msg_) shape(temp)
  write(funit, iostat=ierror_, iomsg=error_msg_) temp

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_real

subroutine write_real_1d(self, funit, array, ierror, error_msg)
  class(ArrayBinFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  real(wp), intent(in)                                 :: array(:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_

  write(funit, iostat=ierror_, iomsg=error_msg_) shape(array)
  write(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_real_1d

subroutine write_real_2d(self, funit, array, ierror, error_msg)
  class(ArrayBinFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  real(wp), intent(in)                                 :: array(:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_

  write(funit, iostat=ierror_, iomsg=error_msg_) shape(array)
  write(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_real_2d

subroutine write_real_3d(self, funit, array, ierror, error_msg)
  class(ArrayBinFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  real(wp), intent(in)                                 :: array(:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_

  write(funit, iostat=ierror_, iomsg=error_msg_) shape(array)
  write(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_real_3d

subroutine write_real_4d(self, funit, array, ierror, error_msg)
  class(ArrayBinFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  real(wp), intent(in)                                 :: array(:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_

  write(funit, iostat=ierror_, iomsg=error_msg_) shape(array)
  write(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_real_4d

subroutine write_real_5d(self, funit, array, ierror, error_msg)
  class(ArrayBinFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  real(wp), intent(in)                                 :: array(:,:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_

  write(funit, iostat=ierror_, iomsg=error_msg_) shape(array)
  write(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_real_5d

subroutine write_cmplx(self, funit, array, ierror, error_msg)
  class(ArrayBinFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  complex(wp), intent(in)                              :: array
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: ierror_
  integer                                              :: temp(1)
  character(len=charlen)                               :: error_msg_

  temp(1) = array

  write(funit, iostat=ierror_, iomsg=error_msg_) shape(temp)
  write(funit, iostat=ierror_, iomsg=error_msg_) temp

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_cmplx

subroutine write_cmplx_1d(self, funit, array, ierror, error_msg)
  class(ArrayBinFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  complex(wp), intent(in)                              :: array(:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_

  write(funit, iostat=ierror_, iomsg=error_msg_) shape(array)
  write(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_cmplx_1d

subroutine write_cmplx_2d(self, funit, array, ierror, error_msg)
  class(ArrayBinFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  complex(wp), intent(in)                              :: array(:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_

  write(funit, iostat=ierror_, iomsg=error_msg_) shape(array)
  write(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_cmplx_2d

subroutine write_cmplx_3d(self, funit, array, ierror, error_msg)
  class(ArrayBinFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  complex(wp), intent(in)                              :: array(:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_

  write(funit, iostat=ierror_, iomsg=error_msg_) shape(array)
  write(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_cmplx_3d

subroutine write_cmplx_4d(self, funit, array, ierror, error_msg)
  class(ArrayBinFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  complex(wp), intent(in)                              :: array(:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_

  write(funit, iostat=ierror_, iomsg=error_msg_) shape(array)
  write(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_cmplx_4d

subroutine write_cmplx_5d(self, funit, array, ierror, error_msg)
  class(ArrayBinFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  complex(wp), intent(in)                              :: array(:,:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_

  write(funit, iostat=ierror_, iomsg=error_msg_) shape(array)
  write(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_cmplx_5d


!********************************************************************************
!
! Read an array in bin format.
!
! Binary files are formatted as follows:
!    shape(array), array
!
!*******************************************************************************
subroutine read_int(self, funit, array, ierror, error_msg)
  class(ArrayBinFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  integer, intent(out)                                 :: array
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 1
  integer                                              :: shape_(ndim), ierror_
  integer, allocatable                                 :: temp(:)
  character(len=charlen)                               :: error_msg_

  read(funit, iostat=ierror_, iomsg=error_msg_) shape_
  allocate(temp(shape_(1)), stat=ierror_, errmsg=error_msg_)
  read(funit, iostat=ierror_, iomsg=error_msg_) temp

  array = temp(1)

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_int

subroutine read_int_1d(self, funit, array, ierror, error_msg)
  class(ArrayBinFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  integer, allocatable, intent(out)                    :: array(:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 1
  integer                                              :: shape_(ndim), ierror_
  character(len=charlen)                               :: error_msg_

  read(funit, iostat=ierror_, iomsg=error_msg_) shape_
  allocate(array(shape_(1)), stat=ierror_, errmsg=error_msg_)
  read(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_int_1d

subroutine read_int_2d(self, funit, array, ierror, error_msg)
  class(ArrayBinFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  integer, allocatable, intent(out)                    :: array(:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 2
  integer                                              :: shape_(ndim), ierror_
  character(len=charlen)                               :: error_msg_

  read(funit, iostat=ierror_, iomsg=error_msg_) shape_
  allocate(array(shape_(1), shape_(2)), stat=ierror_, errmsg=error_msg_)
  read(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_int_2d

subroutine read_int_3d(self, funit, array, ierror, error_msg)
  class(ArrayBinFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  integer, allocatable, intent(out)                    :: array(:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 3
  integer                                              :: shape_(ndim), ierror_
  character(len=charlen)                               :: error_msg_

  read(funit, iostat=ierror_, iomsg=error_msg_) shape_
  allocate(array(shape_(1), shape_(2), shape_(3)), stat=ierror_, errmsg=error_msg_)
  read(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_int_3d

subroutine read_int_4d(self, funit, array, ierror, error_msg)
  class(ArrayBinFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  integer, allocatable, intent(out)                    :: array(:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 4
  integer                                              :: shape_(ndim), ierror_
  character(len=charlen)                               :: error_msg_

  read(funit, iostat=ierror_, iomsg=error_msg_) shape_
  allocate(array(shape_(1), shape_(2), shape_(3), shape_(4)), stat=ierror_, errmsg=error_msg_)
  read(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_int_4d

subroutine read_int_5d(self, funit, array, ierror, error_msg)
  class(ArrayBinFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  integer, allocatable, intent(out)                    :: array(:,:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 5
  integer                                              :: shape_(ndim), ierror_
  character(len=charlen)                               :: error_msg_

  read(funit, iostat=ierror_, iomsg=error_msg_) shape_
  allocate(array(shape_(1), shape_(2), shape_(3), shape_(4), shape_(5)), stat=ierror_, errmsg=error_msg_)
  read(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_int_5d

subroutine read_real(self, funit, array, ierror, error_msg)
  class(ArrayBinFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  real(wp), intent(out)                                :: array
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 1
  integer                                              :: shape_(ndim), ierror_
  real(wp), allocatable                                :: temp(:)
  character(len=charlen)                               :: error_msg_

  read(funit, iostat=ierror_, iomsg=error_msg_) shape_
  allocate(temp(shape_(1)), stat=ierror_, errmsg=error_msg_)
  read(funit, iostat=ierror_, iomsg=error_msg_) temp

  array = temp(1)

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_real

subroutine read_real_1d(self, funit, array, ierror, error_msg)
  class(ArrayBinFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  real(wp), allocatable, intent(out)                   :: array(:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 1
  integer                                              :: shape_(ndim), ierror_
  character(len=charlen)                               :: error_msg_

  read(funit, iostat=ierror_, iomsg=error_msg_) shape_
  allocate(array(shape_(1)), stat=ierror_, errmsg=error_msg_)
  read(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_real_1d

subroutine read_real_2d(self, funit, array, ierror, error_msg)
  class(ArrayBinFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  real(wp), allocatable, intent(out)                   :: array(:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 2
  integer                                              :: shape_(ndim), ierror_
  character(len=charlen)                               :: error_msg_

  read(funit, iostat=ierror_, iomsg=error_msg_) shape_
  allocate(array(shape_(1), shape_(2)), stat=ierror_, errmsg=error_msg_)
  read(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_real_2d

subroutine read_real_3d(self, funit, array, ierror, error_msg)
  class(ArrayBinFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  real(wp), allocatable, intent(out)                   :: array(:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 3
  integer                                              :: shape_(ndim), ierror_
  character(len=charlen)                               :: error_msg_

  read(funit, iostat=ierror_, iomsg=error_msg_) shape_
  allocate(array(shape_(1), shape_(2), shape_(3)), stat=ierror_, errmsg=error_msg_)
  read(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_real_3d

subroutine read_real_4d(self, funit, array, ierror, error_msg)
  class(ArrayBinFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  real(wp), allocatable, intent(out)                   :: array(:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 4
  integer                                              :: shape_(ndim), ierror_
  character(len=charlen)                               :: error_msg_

  read(funit, iostat=ierror_, iomsg=error_msg_) shape_
  allocate(array(shape_(1), shape_(2), shape_(3), shape_(4)), stat=ierror_, errmsg=error_msg_)
  read(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_real_4d

subroutine read_real_5d(self, funit, array, ierror, error_msg)
  class(ArrayBinFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  real(wp), allocatable, intent(out)                   :: array(:,:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 5
  integer                                              :: shape_(ndim), ierror_
  character(len=charlen)                               :: error_msg_

  read(funit, iostat=ierror_, iomsg=error_msg_) shape_
  allocate(array(shape_(1), shape_(2), shape_(3), shape_(4), shape_(5)), stat=ierror_, errmsg=error_msg_)
  read(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_real_5d

subroutine read_cmplx(self, funit, array, ierror, error_msg)
  class(ArrayBinFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  complex(wp), intent(out)                             :: array
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 1
  integer                                              :: shape_(ndim), ierror_
  complex(wp), allocatable                             :: temp(:)
  character(len=charlen)                               :: error_msg_

  read(funit, iostat=ierror_, iomsg=error_msg_) shape_
  allocate(temp(shape_(1)), stat=ierror_, errmsg=error_msg_)
  read(funit, iostat=ierror_, iomsg=error_msg_) temp

  array = temp(1)

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_cmplx

subroutine read_cmplx_1d(self, funit, array, ierror, error_msg)
  class(ArrayBinFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  complex(wp), allocatable, intent(out)                :: array(:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 1
  integer                                              :: shape_(ndim), ierror_
  character(len=charlen)                               :: error_msg_

  read(funit, iostat=ierror_, iomsg=error_msg_) shape_
  allocate(array(shape_(1)), stat=ierror_, errmsg=error_msg_)
  read(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_cmplx_1d

subroutine read_cmplx_2d(self, funit, array, ierror, error_msg)
  class(ArrayBinFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  complex(wp), allocatable, intent(out)                :: array(:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 2
  integer                                              :: shape_(ndim), ierror_
  character(len=charlen)                               :: error_msg_

  read(funit, iostat=ierror_, iomsg=error_msg_) shape_
  allocate(array(shape_(1), shape_(2)), stat=ierror_, errmsg=error_msg_)
  read(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_cmplx_2d

subroutine read_cmplx_3d(self, funit, array, ierror, error_msg)
  class(ArrayBinFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  complex(wp), allocatable, intent(out)                :: array(:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 3
  integer                                              :: shape_(ndim), ierror_
  character(len=charlen)                               :: error_msg_

  read(funit, iostat=ierror_, iomsg=error_msg_) shape_
  allocate(array(shape_(1), shape_(2), shape_(3)), stat=ierror_, errmsg=error_msg_)
  read(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_cmplx_3d

subroutine read_cmplx_4d(self, funit, array, ierror, error_msg)
  class(ArrayBinFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  complex(wp), allocatable, intent(out)                :: array(:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 4
  integer                                              :: shape_(ndim), ierror_
  character(len=charlen)                               :: error_msg_

  read(funit, iostat=ierror_, iomsg=error_msg_) shape_
  allocate(array(shape_(1), shape_(2), shape_(3), shape_(4)), stat=ierror_, errmsg=error_msg_)
  read(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_cmplx_4d

subroutine read_cmplx_5d(self, funit, array, ierror, error_msg)
  class(ArrayBinFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  complex(wp), allocatable, intent(out)                :: array(:,:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 5
  integer                                              :: shape_(ndim), ierror_
  character(len=charlen)                               :: error_msg_

  read(funit, iostat=ierror_, iomsg=error_msg_) shape_
  allocate(array(shape_(1), shape_(2), shape_(3), shape_(4), shape_(5)), stat=ierror_, errmsg=error_msg_)
  read(funit, iostat=ierror_, iomsg=error_msg_) array

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_cmplx_5d


!********************************************************************************
!
! Save an array to the file in txt format.
!
! Binary files are formatted as follows:
!    shape(array), array
!
!*******************************************************************************
subroutine save_int(self, fname, array, ierror, error_msg)
  class(ArraybinFileIO), intent(in)                    :: self 
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
  class(ArraybinFileIO), intent(in)                    :: self 
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
  class(ArraybinFileIO), intent(in)                    :: self 
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
  class(ArraybinFileIO), intent(in)                    :: self 
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
  class(ArraybinFileIO), intent(in)                    :: self 
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
  class(ArraybinFileIO), intent(in)                    :: self 
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
  class(ArraybinFileIO), intent(in)                    :: self 
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
  class(ArraybinFileIO), intent(in)                    :: self 
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
  class(ArraybinFileIO), intent(in)                    :: self 
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
  class(ArraybinFileIO), intent(in)                    :: self 
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
  class(ArraybinFileIO), intent(in)                    :: self 
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
  class(ArraybinFileIO), intent(in)                    :: self 
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
  class(ArraybinFileIO), intent(in)                    :: self 
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
  class(ArraybinFileIO), intent(in)                    :: self 
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
  class(ArraybinFileIO), intent(in)                    :: self 
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
  class(ArraybinFileIO), intent(in)                    :: self 
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
  class(ArraybinFileIO), intent(in)                    :: self 
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
  class(ArraybinFileIO), intent(in)                    :: self 
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
  class(ArrayBinFileIO), intent(in)                    :: self 
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
  class(ArrayBinFileIO), intent(in)                    :: self 
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
  class(ArrayBinFileIO), intent(in)                    :: self 
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
  class(ArrayBinFileIO), intent(in)                    :: self 
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
  class(ArrayBinFileIO), intent(in)                    :: self 
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
  class(ArrayBinFileIO), intent(in)                    :: self 
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
  class(ArrayBinFileIO), intent(in)                    :: self 
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
  class(ArrayBinFileIO), intent(in)                    :: self 
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
  class(ArrayBinFileIO), intent(in)                    :: self 
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
  class(ArrayBinFileIO), intent(in)                    :: self 
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
  class(ArrayBinFileIO), intent(in)                    :: self 
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
  class(ArrayBinFileIO), intent(in)                    :: self 
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
  class(ArrayBinFileIO), intent(in)                    :: self 
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
  class(ArrayBinFileIO), intent(in)                    :: self 
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
  class(ArrayBinFileIO), intent(in)                    :: self 
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
  class(ArrayBinFileIO), intent(in)                    :: self 
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
  class(ArrayBinFileIO), intent(in)                    :: self 
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
  class(ArrayBinFileIO), intent(in)                    :: self 
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

end module array_bin_file_io