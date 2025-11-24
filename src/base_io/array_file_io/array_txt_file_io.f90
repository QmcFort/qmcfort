! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module array_txt_file_io

use constants
use string, only: str
use file_handle, only: FileHandle
use array_file_io, only: ArrayFileIO

implicit none

private
public :: ArrayTxtFileIO

character(len=*), parameter :: indx_format = "i6"
character(len=*), parameter :: int_format = "i10"
character(len=*), parameter :: real_format = "f18.10"
character(len=*), parameter :: cmplx_format = "2f18.10"

type, extends(ArrayFileIO) :: ArrayTxtFileIO
contains
  procedure, private :: get_shape_format
  procedure, private :: get_array_format

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
end type ArrayTxtFileIO

contains

!********************************************************************************
!
! Format descriptor for the shape of n-dim array
!
!********************************************************************************
function get_shape_format(self, n) result(shape_format)
  class(ArrayTxtFileIO), intent(in) :: self
  integer, intent(in)               :: n
  character(len=:), allocatable     :: shape_format

  shape_format = "(" // str(n) // indx_format // ")"
end function get_shape_format


!********************************************************************************
!
! Format descriptor for the n-dim array
!
!********************************************************************************
function get_array_format(self, n, val_format) result(array_format)
  class(ArrayTxtFileIO), intent(in) :: self
  integer, intent(in)               :: n
  character(len=*), intent(in)      :: val_format
  character(len=:), allocatable     :: array_format

  array_format = "(" // str(n) // indx_format // "," // val_format // ")"
end function get_array_format


!********************************************************************************
!
! Write an array in txt format.
!
! For example, for a 3d array, the format will look like:
!    ni nj nk
!    i j k array(i,j,k)      for all i, j, k
!
!*******************************************************************************
subroutine write_int(self, funit, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  integer, intent(in)                                  :: array
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: p
  integer                                              :: ierror_
  integer                                              :: temp(1)
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: shape_format, array_format

  shape_format = self%get_shape_format(rank(temp))
  array_format = self%get_array_format(rank(temp), int_format)

  temp(1) = array

  write(funit, shape_format, iostat=ierror_, iomsg=error_msg_) shape(temp)

  do p = 1, size(temp)
    write(funit, array_format, iostat=ierror_, iomsg=error_msg_) p, temp(p)
  end do

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_int

subroutine write_int_1d(self, funit, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  integer, intent(in)                                  :: array(:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: p
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: shape_format, array_format

  shape_format = self%get_shape_format(rank(array))
  array_format = self%get_array_format(rank(array), int_format)

  write(funit, shape_format, iostat=ierror_, iomsg=error_msg_) shape(array)

  do p = 1, size(array, 1)
    write(funit, array_format, iostat=ierror_, iomsg=error_msg_) p, array(p)
  end do

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_int_1d

subroutine write_int_2d(self, funit, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  integer, intent(in)                                  :: array(:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: p, q
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: shape_format, array_format

  shape_format = self%get_shape_format(rank(array))
  array_format = self%get_array_format(rank(array), int_format)

  write(funit, shape_format, iostat=ierror_, iomsg=error_msg_) shape(array)

  do q = 1, size(array, 2)
    do p = 1, size(array, 1)
      write(funit, array_format, iostat=ierror_, iomsg=error_msg_) p, q, array(p,q)
    end do
  end do

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_int_2d

subroutine write_int_3d(self, funit, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  integer, intent(in)                                  :: array(:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: p, q, r
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: shape_format, array_format

  shape_format = self%get_shape_format(rank(array))
  array_format = self%get_array_format(rank(array), int_format)

  write(funit, shape_format, iostat=ierror_, iomsg=error_msg_) shape(array)

  do r = 1, size(array, 3)
    do q = 1, size(array, 2)
      do p = 1, size(array, 1)
        write(funit, array_format, iostat=ierror_, iomsg=error_msg_) p, q, r, array(p,q,r)
      end do
    end do
  end do

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_int_3d

subroutine write_int_4d(self, funit, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  integer, intent(in)                                  :: array(:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: p, q, r, s
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: shape_format, array_format

  shape_format = self%get_shape_format(rank(array))
  array_format = self%get_array_format(rank(array), int_format)

  write(funit, shape_format, iostat=ierror_, iomsg=error_msg_) shape(array)

  do s = 1, size(array, 4)
    do r = 1, size(array, 3)
      do q = 1, size(array, 2)
        do p = 1, size(array, 1)
          write(funit, array_format, iostat=ierror_, iomsg=error_msg_) p, q, r, s, array(p,q,r,s)
        end do
      end do
    end do
  end do

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_int_4d

subroutine write_int_5d(self, funit, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  integer, intent(in)                                  :: array(:,:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: p, q, r, s, t
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: shape_format, array_format

  shape_format = self%get_shape_format(rank(array))
  array_format = self%get_array_format(rank(array), int_format)

  write(funit, shape_format, iostat=ierror_, iomsg=error_msg_) shape(array)

  do t = 1, size(array, 5)
    do s = 1, size(array, 4)
      do r = 1, size(array, 3)
        do q = 1, size(array, 2)
          do p = 1, size(array, 1)
            write(funit, array_format, iostat=ierror_, iomsg=error_msg_) p, q, r, s, t, array(p,q,r,s,t)
          end do
        end do
      end do
    end do
  end do

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_int_5d

subroutine write_real(self, funit, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  real(wp), intent(in)                                 :: array
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: p
  integer                                              :: ierror_
  real(wp)                                             :: temp(1)
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: shape_format, array_format

  shape_format = self%get_shape_format(rank(temp))
  array_format = self%get_array_format(rank(temp), real_format)

  temp(1) = array

  write(funit, shape_format, iostat=ierror_, iomsg=error_msg_) shape(temp)

  do p = 1, size(temp, 1)
    write(funit, array_format, iostat=ierror_, iomsg=error_msg_) p, temp(p)
  end do

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_real

subroutine write_real_1d(self, funit, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  real(wp), intent(in)                                 :: array(:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: p
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: shape_format, array_format

  shape_format = self%get_shape_format(rank(array))
  array_format = self%get_array_format(rank(array), real_format)

  write(funit, shape_format, iostat=ierror_, iomsg=error_msg_) shape(array)

  do p = 1, size(array, 1)
    write(funit, array_format, iostat=ierror_, iomsg=error_msg_) p, array(p)
  end do

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_real_1d

subroutine write_real_2d(self, funit, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  real(wp), intent(in)                                 :: array(:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: p, q
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: shape_format, array_format

  shape_format = self%get_shape_format(rank(array))
  array_format = self%get_array_format(rank(array), real_format)

  write(funit, shape_format, iostat=ierror_, iomsg=error_msg_) shape(array)

  do q = 1, size(array, 2)
    do p = 1, size(array, 1)
      write(funit, array_format, iostat=ierror_, iomsg=error_msg_) p, q, array(p,q)
    end do
  end do

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_real_2d

subroutine write_real_3d(self, funit, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  real(wp), intent(in)                                 :: array(:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: p, q, r
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: shape_format, array_format

  shape_format = self%get_shape_format(rank(array))
  array_format = self%get_array_format(rank(array), real_format)

  write(funit, shape_format, iostat=ierror_, iomsg=error_msg_) shape(array)

  do r = 1, size(array, 3)
    do q = 1, size(array, 2)
      do p = 1, size(array, 1)
        write(funit, array_format, iostat=ierror_, iomsg=error_msg_) p, q, r, array(p,q,r)
      end do
    end do
  end do

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_real_3d

subroutine write_real_4d(self, funit, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  real(wp), intent(in)                                 :: array(:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: p, q, r, s
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: shape_format, array_format

  shape_format = self%get_shape_format(rank(array))
  array_format = self%get_array_format(rank(array), real_format)

  write(funit, shape_format, iostat=ierror_, iomsg=error_msg_) shape(array)

  do s = 1, size(array, 4)
    do r = 1, size(array, 3)
      do q = 1, size(array, 2)
        do p = 1, size(array, 1)
          write(funit, array_format, iostat=ierror_, iomsg=error_msg_) p, q, r, s, array(p,q,r,s)
        end do
      end do
    end do
  end do

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_real_4d

subroutine write_real_5d(self, funit, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  real(wp), intent(in)                                 :: array(:,:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: p, q, r, s, t
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: shape_format, array_format

  shape_format = self%get_shape_format(rank(array))
  array_format = self%get_array_format(rank(array), real_format)

  write(funit, shape_format, iostat=ierror_, iomsg=error_msg_) shape(array)

  do t = 1, size(array, 5)
    do s = 1, size(array, 4)
      do r = 1, size(array, 3)
        do q = 1, size(array, 2)
          do p = 1, size(array, 1)
            write(funit, array_format, iostat=ierror_, iomsg=error_msg_) p, q, r, s, t, array(p,q,r,s,t)
          end do
        end do
      end do
    end do
  end do

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_real_5d

subroutine write_cmplx(self, funit, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  complex(wp), intent(in)                              :: array
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: p
  integer                                              :: ierror_
  real(wp)                                             :: temp(1)
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: shape_format, array_format

  shape_format = self%get_shape_format(rank(temp))
  array_format = self%get_array_format(rank(temp), cmplx_format)

  temp(1) = array

  write(funit, shape_format, iostat=ierror_, iomsg=error_msg_) shape(temp)

  do p = 1, size(temp, 1)
    write(funit, array_format, iostat=ierror_, iomsg=error_msg_) p, temp(p)
  end do

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_cmplx

subroutine write_cmplx_1d(self, funit, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  complex(wp), intent(in)                              :: array(:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: p
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: shape_format, array_format

  shape_format = self%get_shape_format(rank(array))
  array_format = self%get_array_format(rank(array), cmplx_format)

  write(funit, shape_format, iostat=ierror_, iomsg=error_msg_) shape(array)

  do p = 1, size(array, 1)
    write(funit, array_format, iostat=ierror_, iomsg=error_msg_) p, array(p)
  end do

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_cmplx_1d

subroutine write_cmplx_2d(self, funit, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  complex(wp), intent(in)                              :: array(:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: p, q
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: shape_format, array_format

  shape_format = self%get_shape_format(rank(array))
  array_format = self%get_array_format(rank(array), cmplx_format)

  write(funit, shape_format, iostat=ierror_, iomsg=error_msg_) shape(array)

  do q = 1, size(array, 2)
    do p = 1, size(array, 1)
      write(funit, array_format, iostat=ierror_, iomsg=error_msg_) p, q, array(p,q)
    end do
  end do

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_cmplx_2d

subroutine write_cmplx_3d(self, funit, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  complex(wp), intent(in)                              :: array(:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: p, q, r
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: shape_format, array_format

  shape_format = self%get_shape_format(rank(array))
  array_format = self%get_array_format(rank(array), cmplx_format)

  write(funit, shape_format, iostat=ierror_, iomsg=error_msg_) shape(array)

  do r = 1, size(array, 3)
    do q = 1, size(array, 2)
      do p = 1, size(array, 1)
        write(funit, array_format, iostat=ierror_, iomsg=error_msg_) p, q, r, array(p,q,r)
      end do
    end do
  end do

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_cmplx_3d

subroutine write_cmplx_4d(self, funit, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  complex(wp), intent(in)                              :: array(:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: p, q, r, s
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: shape_format, array_format

  shape_format = self%get_shape_format(rank(array))
  array_format = self%get_array_format(rank(array), cmplx_format)

  write(funit, shape_format, iostat=ierror_, iomsg=error_msg_) shape(array)

  do s = 1, size(array, 4)
    do r = 1, size(array, 3)
      do q = 1, size(array, 2)
        do p = 1, size(array, 1)
          write(funit, array_format, iostat=ierror_, iomsg=error_msg_) p, q, r, s, array(p,q,r,s)
        end do
      end do
    end do
  end do

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_cmplx_4d

subroutine write_cmplx_5d(self, funit, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  complex(wp), intent(in)                              :: array(:,:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer                                              :: p, q, r, s, t
  integer                                              :: ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: shape_format, array_format

  shape_format = self%get_shape_format(rank(array))
  array_format = self%get_array_format(rank(array), cmplx_format)

  write(funit, shape_format, iostat=ierror_, iomsg=error_msg_) shape(array)

  do t = 1, size(array, 5)
    do s = 1, size(array, 4)
      do r = 1, size(array, 3)
        do q = 1, size(array, 2)
          do p = 1, size(array, 1)
            write(funit, array_format, iostat=ierror_, iomsg=error_msg_) p, q, r, s, t, array(p,q,r,s,t)
          end do
        end do
      end do
    end do
  end do

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine write_cmplx_5d


!********************************************************************************
!
! Read array in txt format.
!
! For example, for a 3d array, the format will look like:
!    ni nj nk
!    i j k array(i,j,k)      for all i, j, k
!
!*******************************************************************************
subroutine read_int(self, funit, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  integer, intent(out)                                 :: array
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 1
  integer                                              :: p
  integer                                              :: shape_(1), indx(1), ierror_
  integer, allocatable                                 :: temp(:)
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: shape_format, array_format

  shape_format = self%get_shape_format(ndim)
  array_format = self%get_array_format(ndim, int_format)

  read(funit, shape_format, iostat=ierror_, iomsg=error_msg_) shape_
  allocate(temp(shape_(1)), stat=ierror_, errmsg=error_msg_)

  do p = 1, size(temp)
    read(funit, array_format, iostat=ierror_, iomsg=error_msg_) indx, temp(p)
  end do

  array = temp(1)

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_int

subroutine read_int_1d(self, funit, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  integer, allocatable, intent(out)                    :: array(:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 1
  integer                                              :: p
  integer                                              :: shape_(1), indx(1), ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: shape_format, array_format

  shape_format = self%get_shape_format(ndim)
  array_format = self%get_array_format(ndim, int_format)
  
  read(funit, shape_format, iostat=ierror_, iomsg=error_msg_) shape_
  allocate(array(shape_(1)), stat=ierror_, errmsg=error_msg_)

  do p = 1, size(array)
    read(funit, array_format, iostat=ierror_, iomsg=error_msg_) indx, array(p)
  end do

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_int_1d

subroutine read_int_2d(self, funit, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  integer, allocatable, intent(out)                    :: array(:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 2
  integer                                              :: p, q
  integer                                              :: shape_(ndim), indx(ndim), ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: shape_format, array_format

  shape_format = self%get_shape_format(ndim)
  array_format = self%get_array_format(ndim, int_format)
  
  read(funit, shape_format, iostat=ierror_, iomsg=error_msg_) shape_
  allocate(array(shape_(1), shape_(2)), stat=ierror_, errmsg=error_msg_)

  do q = 1, size(array, 2)
    do p = 1, size(array, 1)
      read(funit, array_format, iostat=ierror_, iomsg=error_msg_) indx, array(p,q)
    end do
  end do

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_int_2d

subroutine read_int_3d(self, funit, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  integer, allocatable, intent(out)                    :: array(:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 3
  integer                                              :: p, q, r
  integer                                              :: shape_(ndim), indx(ndim), ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: shape_format, array_format

  shape_format = self%get_shape_format(ndim)
  array_format = self%get_array_format(ndim, int_format)
  
  read(funit, shape_format, iostat=ierror_, iomsg=error_msg_) shape_
  allocate(array(shape_(1), shape_(2), shape_(3)), stat=ierror_, errmsg=error_msg_)

  do r = 1, size(array, 3)
    do q = 1, size(array, 2)
      do p = 1, size(array, 1)
        read(funit, array_format, iostat=ierror_, iomsg=error_msg_) indx, array(p,q,r)
      end do
    end do
  end do

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_int_3d

subroutine read_int_4d(self, funit, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  integer, allocatable, intent(out)                    :: array(:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 4
  integer                                              :: p, q, r, s
  integer                                              :: shape_(ndim), indx(ndim), ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: shape_format, array_format

  shape_format = self%get_shape_format(ndim)
  array_format = self%get_array_format(ndim, int_format)
  
  read(funit, shape_format, iostat=ierror_, iomsg=error_msg_) shape_
  allocate(array(shape_(1), shape_(2), shape_(3), shape_(4)), stat=ierror_, errmsg=error_msg_)

  do s = 1, size(array, 4)
    do r = 1, size(array, 3)
      do q = 1, size(array, 2)
        do p = 1, size(array, 1)
          read(funit, array_format, iostat=ierror_, iomsg=error_msg_) indx, array(p,q,r,s)
        end do
      end do
    end do
  end do

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_int_4d

subroutine read_int_5d(self, funit, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  integer, allocatable, intent(out)                    :: array(:,:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 5
  integer                                              :: p, q, r, s, t
  integer                                              :: shape_(ndim), indx(ndim), ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: shape_format, array_format

  shape_format = self%get_shape_format(ndim)
  array_format = self%get_array_format(ndim, int_format)
  
  read(funit, shape_format, iostat=ierror_, iomsg=error_msg_) shape_
  allocate(array(shape_(1), shape_(2), shape_(3), shape_(4), shape_(5)), stat=ierror_, errmsg=error_msg_)

  do t = 1, size(array, 5)
    do s = 1, size(array, 4)
      do r = 1, size(array, 3)
        do q = 1, size(array, 2)
          do p = 1, size(array, 1)
            read(funit, array_format, iostat=ierror_, iomsg=error_msg_) indx, array(p,q,r,s,t)
          end do
        end do
      end do
    end do
  end do

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_int_5d

subroutine read_real(self, funit, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  real(wp), intent(out)                                :: array
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 1
  integer                                              :: p
  integer                                              :: shape_(ndim), indx(ndim), ierror_
  real(wp), allocatable                                :: temp(:)
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: shape_format, array_format

  shape_format = self%get_shape_format(ndim)
  array_format = self%get_array_format(ndim, real_format)
  
  read(funit, shape_format, iostat=ierror_, iomsg=error_msg_) shape_
  allocate(temp(shape_(1)), stat=ierror_, errmsg=error_msg_)

  do p = 1, size(temp)
    read(funit, array_format, iostat=ierror_, iomsg=error_msg_) indx, temp(p)
  end do

  array = temp(1)

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_real

subroutine read_real_1d(self, funit, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  real(wp), allocatable, intent(out)                   :: array(:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 1
  integer                                              :: p
  integer                                              :: shape_(ndim), indx(ndim), ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: shape_format, array_format

  shape_format = self%get_shape_format(ndim)
  array_format = self%get_array_format(ndim, real_format)
  
  read(funit, shape_format, iostat=ierror_, iomsg=error_msg_) shape_
  allocate(array(shape_(1)), stat=ierror_, errmsg=error_msg_)

  do p = 1, size(array)
    read(funit, array_format, iostat=ierror_, iomsg=error_msg_) indx, array(p)
  end do

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_real_1d

subroutine read_real_2d(self, funit, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  real(wp), allocatable, intent(out)                   :: array(:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 2
  integer                                              :: p, q
  integer                                              :: shape_(ndim), indx(ndim), ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: shape_format, array_format

  shape_format = self%get_shape_format(ndim)
  array_format = self%get_array_format(ndim, real_format)
  
  read(funit, shape_format, iostat=ierror_, iomsg=error_msg_) shape_
  allocate(array(shape_(1), shape_(2)), stat=ierror_, errmsg=error_msg_)

  do q = 1, size(array, 2)
    do p = 1, size(array, 1)
      read(funit, array_format, iostat=ierror_, iomsg=error_msg_) indx, array(p,q)
    end do
  end do

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_real_2d

subroutine read_real_3d(self, funit, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  real(wp), allocatable, intent(out)                   :: array(:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 3
  integer                                              :: p, q, r
  integer                                              :: shape_(ndim), indx(ndim), ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: shape_format, array_format

  shape_format = self%get_shape_format(ndim)
  array_format = self%get_array_format(ndim, real_format)
  
  read(funit, shape_format, iostat=ierror_, iomsg=error_msg_) shape_
  allocate(array(shape_(1), shape_(2), shape_(3)), stat=ierror_, errmsg=error_msg_)

  do r = 1, size(array, 3)
    do q = 1, size(array, 2)
      do p = 1, size(array, 1)
        read(funit, array_format, iostat=ierror_, iomsg=error_msg_) indx, array(p,q,r)
      end do
    end do
  end do

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_real_3d

subroutine read_real_4d(self, funit, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  real(wp), allocatable, intent(out)                   :: array(:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 4
  integer                                              :: p, q, r, s
  integer                                              :: shape_(ndim), indx(ndim), ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: shape_format, array_format

  shape_format = self%get_shape_format(ndim)
  array_format = self%get_array_format(ndim, real_format)
  
  read(funit, shape_format, iostat=ierror_, iomsg=error_msg_) shape_
  allocate(array(shape_(1), shape_(2), shape_(3), shape_(4)), stat=ierror_, errmsg=error_msg_)

  do s = 1, size(array, 4)
    do r = 1, size(array, 3)
      do q = 1, size(array, 2)
        do p = 1, size(array, 1)
          read(funit, array_format, iostat=ierror_, iomsg=error_msg_) indx, array(p,q,r,s)
        end do
      end do
    end do
  end do

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_real_4d

subroutine read_real_5d(self, funit, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  real(wp), allocatable, intent(out)                   :: array(:,:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 5
  integer                                              :: p, q, r, s, t
  integer                                              :: shape_(ndim), indx(ndim), ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: shape_format, array_format

  shape_format = self%get_shape_format(ndim)
  array_format = self%get_array_format(ndim, real_format)
  
  read(funit, shape_format, iostat=ierror_, iomsg=error_msg_) shape_
  allocate(array(shape_(1), shape_(2), shape_(3), shape_(4), shape_(5)), stat=ierror_, errmsg=error_msg_)

  do t = 1, size(array, 5)
    do s = 1, size(array, 4)
      do r = 1, size(array, 3)
        do q = 1, size(array, 2)
          do p = 1, size(array, 1)
            read(funit, array_format, iostat=ierror_, iomsg=error_msg_) indx, array(p,q,r,s,t)
          end do
        end do
      end do
    end do
  end do

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_real_5d

subroutine read_cmplx(self, funit, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  complex(wp), intent(out)                             :: array
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 1
  integer                                              :: p
  integer                                              :: shape_(ndim), indx(ndim), ierror_
  complex(wp), allocatable                             :: temp(:)
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: shape_format, array_format

  shape_format = self%get_shape_format(ndim)
  array_format = self%get_array_format(ndim, cmplx_format)
  
  read(funit, shape_format, iostat=ierror_, iomsg=error_msg_) shape_
  allocate(temp(shape_(1)), stat=ierror_, errmsg=error_msg_)

  do p = 1, size(temp)
    read(funit, array_format, iostat=ierror_, iomsg=error_msg_) indx, temp(p)
  end do

  array = temp(1)

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_cmplx

subroutine read_cmplx_1d(self, funit, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  complex(wp), allocatable, intent(out)                :: array(:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 1
  integer                                              :: p
  integer                                              :: shape_(ndim), indx(ndim), ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: shape_format, array_format

  shape_format = self%get_shape_format(ndim)
  array_format = self%get_array_format(ndim, cmplx_format)
  
  read(funit, shape_format, iostat=ierror_, iomsg=error_msg_) shape_
  allocate(array(shape_(1)), stat=ierror_, errmsg=error_msg_)

  do p = 1, size(array)
    read(funit, array_format, iostat=ierror_, iomsg=error_msg_) indx, array(p)
  end do

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_cmplx_1d

subroutine read_cmplx_2d(self, funit, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  complex(wp), allocatable, intent(out)                :: array(:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 2
  integer                                              :: p, q
  integer                                              :: shape_(ndim), indx(ndim), ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: shape_format, array_format

  shape_format = self%get_shape_format(ndim)
  array_format = self%get_array_format(ndim, cmplx_format)
  
  read(funit, shape_format, iostat=ierror_, iomsg=error_msg_) shape_
  allocate(array(shape_(1), shape_(2)), stat=ierror_, errmsg=error_msg_)

  do q = 1, size(array, 2)
    do p = 1, size(array, 1)
      read(funit, array_format, iostat=ierror_, iomsg=error_msg_) indx, array(p,q)
    end do
  end do

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_cmplx_2d

subroutine read_cmplx_3d(self, funit, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  complex(wp), allocatable, intent(out)                :: array(:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 3
  integer                                              :: p, q, r
  integer                                              :: shape_(ndim), indx(ndim), ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: shape_format, array_format

  shape_format = self%get_shape_format(ndim)
  array_format = self%get_array_format(ndim, cmplx_format)
  
  read(funit, shape_format, iostat=ierror_, iomsg=error_msg_) shape_
  allocate(array(shape_(1), shape_(2), shape_(3)), stat=ierror_, errmsg=error_msg_)

  do r = 1, size(array, 3)
    do q = 1, size(array, 2)
      do p = 1, size(array, 1)
        read(funit, array_format, iostat=ierror_, iomsg=error_msg_) indx, array(p,q,r)
      end do
    end do
  end do

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_cmplx_3d

subroutine read_cmplx_4d(self, funit, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  complex(wp), allocatable, intent(out)                :: array(:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 4
  integer                                              :: p, q, r, s
  integer                                              :: shape_(ndim), indx(ndim), ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: shape_format, array_format

  shape_format = self%get_shape_format(ndim)
  array_format = self%get_array_format(ndim, cmplx_format)
  
  read(funit, shape_format, iostat=ierror_, iomsg=error_msg_) shape_
  allocate(array(shape_(1), shape_(2), shape_(3), shape_(4)), stat=ierror_, errmsg=error_msg_)

  do s = 1, size(array, 4)
    do r = 1, size(array, 3)
      do q = 1, size(array, 2)
        do p = 1, size(array, 1)
          read(funit, array_format, iostat=ierror_, iomsg=error_msg_) indx, array(p,q,r,s)
        end do
      end do
    end do
  end do

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_cmplx_4d

subroutine read_cmplx_5d(self, funit, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  integer, intent(in)                                  :: funit
  complex(wp), allocatable, intent(out)                :: array(:,:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  integer, parameter                                   :: ndim = 5
  integer                                              :: p, q, r, s, t
  integer                                              :: shape_(ndim), indx(ndim), ierror_
  character(len=charlen)                               :: error_msg_
  character(len=:), allocatable                        :: shape_format, array_format

  shape_format = self%get_shape_format(ndim)
  array_format = self%get_array_format(ndim, cmplx_format)
  
  read(funit, shape_format, iostat=ierror_, iomsg=error_msg_) shape_
  allocate(array(shape_(1), shape_(2), shape_(3), shape_(4), shape_(5)), stat=ierror_, errmsg=error_msg_)

  do t = 1, size(array, 5)
    do s = 1, size(array, 4)
      do r = 1, size(array, 3)
        do q = 1, size(array, 2)
          do p = 1, size(array, 1)
            read(funit, array_format, iostat=ierror_, iomsg=error_msg_) indx, array(p,q,r,s,t)
          end do
        end do
      end do
    end do
  end do

  if (present(ierror)) ierror = ierror_
  if (present(error_msg)) error_msg = error_msg_
end subroutine read_cmplx_5d


!********************************************************************************
!
! Save an array to the file in txt format.
!
! For example, for a 3d array, the format will look like:
!    ni nj nk
!    i j k array(i,j,k)      for all i, j, k
!
!*******************************************************************************
subroutine save_int(self, fname, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  character(len=*), intent(in)                         :: fname
  integer, intent(in)                                  :: array
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="write", status="replace", access="sequential", form="formatted")
  call self%write(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine save_int

subroutine save_int_1d(self, fname, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  character(len=*), intent(in)                         :: fname
  integer, intent(in)                                  :: array(:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="write", status="replace", access="sequential", form="formatted")
  call self%write(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine save_int_1d

subroutine save_int_2d(self, fname, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  character(len=*), intent(in)                         :: fname
  integer, intent(in)                                  :: array(:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="write", status="replace", access="sequential", form="formatted")
  call self%write(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine save_int_2d

subroutine save_int_3d(self, fname, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  character(len=*), intent(in)                         :: fname
  integer, intent(in)                                  :: array(:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="write", status="replace", access="sequential", form="formatted")
  call self%write(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine save_int_3d

subroutine save_int_4d(self, fname, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  character(len=*), intent(in)                         :: fname
  integer, intent(in)                                  :: array(:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="write", status="replace", access="sequential", form="formatted")
  call self%write(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine save_int_4d

subroutine save_int_5d(self, fname, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  character(len=*), intent(in)                         :: fname
  integer, intent(in)                                  :: array(:,:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="write", status="replace", access="sequential", form="formatted")
  call self%write(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine save_int_5d

subroutine save_real(self, fname, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  character(len=*), intent(in)                         :: fname
  real(wp), intent(in)                                 :: array
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="write", status="replace", access="sequential", form="formatted")
  call self%write(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine save_real

subroutine save_real_1d(self, fname, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  character(len=*), intent(in)                         :: fname
  real(wp), intent(in)                                 :: array(:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="write", status="replace", access="sequential", form="formatted")
  call self%write(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine save_real_1d

subroutine save_real_2d(self, fname, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  character(len=*), intent(in)                         :: fname
  real(wp), intent(in)                                 :: array(:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="write", status="replace", access="sequential", form="formatted")
  call self%write(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine save_real_2d

subroutine save_real_3d(self, fname, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  character(len=*), intent(in)                         :: fname
  real(wp), intent(in)                                 :: array(:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="write", status="replace", access="sequential", form="formatted")
  call self%write(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine save_real_3d

subroutine save_real_4d(self, fname, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  character(len=*), intent(in)                         :: fname
  real(wp), intent(in)                                 :: array(:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="write", status="replace", access="sequential", form="formatted")
  call self%write(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine save_real_4d

subroutine save_real_5d(self, fname, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  character(len=*), intent(in)                         :: fname
  real(wp), intent(in)                                 :: array(:,:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="write", status="replace", access="sequential", form="formatted")
  call self%write(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine save_real_5d

subroutine save_cmplx(self, fname, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  character(len=*), intent(in)                         :: fname
  complex(wp), intent(in)                              :: array
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="write", status="replace", access="sequential", form="formatted")
  call self%write(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine save_cmplx

subroutine save_cmplx_1d(self, fname, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  character(len=*), intent(in)                         :: fname
  complex(wp), intent(in)                              :: array(:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="write", status="replace", access="sequential", form="formatted")
  call self%write(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine save_cmplx_1d

subroutine save_cmplx_2d(self, fname, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  character(len=*), intent(in)                         :: fname
  complex(wp), intent(in)                              :: array(:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="write", status="replace", access="sequential", form="formatted")
  call self%write(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine save_cmplx_2d

subroutine save_cmplx_3d(self, fname, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  character(len=*), intent(in)                         :: fname
  complex(wp), intent(in)                              :: array(:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="write", status="replace", access="sequential", form="formatted")
  call self%write(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine save_cmplx_3d

subroutine save_cmplx_4d(self, fname, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  character(len=*), intent(in)                         :: fname
  complex(wp), intent(in)                              :: array(:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="write", status="replace", access="sequential", form="formatted")
  call self%write(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine save_cmplx_4d

subroutine save_cmplx_5d(self, fname, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  character(len=*), intent(in)                         :: fname
  complex(wp), intent(in)                              :: array(:,:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="write", status="replace", access="sequential", form="formatted")
  call self%write(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine save_cmplx_5d


!********************************************************************************
!
! Load an array from the file in txt format.
!
! For example, for a 3d array, the format will look like:
!    ni nj nk
!    i j k array(i,j,k)      for all i, j, k
!
!*******************************************************************************
subroutine load_int(self, fname, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  character(len=*), intent(in)                         :: fname
  integer, intent(out)                                 :: array
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="read", status="old", access="sequential", form="formatted")
  call self%read(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine load_int

subroutine load_int_1d(self, fname, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  character(len=*), intent(in)                         :: fname
  integer, allocatable, intent(out)                    :: array(:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="read", status="old", access="sequential", form="formatted")
  call self%read(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine load_int_1d

subroutine load_int_2d(self, fname, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  character(len=*), intent(in)                         :: fname
  integer, allocatable, intent(out)                    :: array(:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="read", status="old", access="sequential", form="formatted")
  call self%read(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine load_int_2d

subroutine load_int_3d(self, fname, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  character(len=*), intent(in)                         :: fname
  integer, allocatable, intent(out)                    :: array(:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="read", status="old", access="sequential", form="formatted")
  call self%read(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine load_int_3d

subroutine load_int_4d(self, fname, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  character(len=*), intent(in)                         :: fname
  integer, allocatable, intent(out)                    :: array(:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="read", status="old", access="sequential", form="formatted")
  call self%read(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine load_int_4d

subroutine load_int_5d(self, fname, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  character(len=*), intent(in)                         :: fname
  integer, allocatable, intent(out)                    :: array(:,:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="read", status="old", access="sequential", form="formatted")
  call self%read(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine load_int_5d

subroutine load_real(self, fname, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  character(len=*), intent(in)                         :: fname
  real(wp), intent(out)                                :: array
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="read", status="old", access="sequential", form="formatted")
  call self%read(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine load_real

subroutine load_real_1d(self, fname, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  character(len=*), intent(in)                         :: fname
  real(wp), allocatable, intent(out)                   :: array(:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="read", status="old", access="sequential", form="formatted")
  call self%read(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine load_real_1d

subroutine load_real_2d(self, fname, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  character(len=*), intent(in)                         :: fname
  real(wp), allocatable, intent(out)                   :: array(:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="read", status="old", access="sequential", form="formatted")
  call self%read(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine load_real_2d

subroutine load_real_3d(self, fname, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  character(len=*), intent(in)                         :: fname
  real(wp), allocatable, intent(out)                   :: array(:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="read", status="old", access="sequential", form="formatted")
  call self%read(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine load_real_3d

subroutine load_real_4d(self, fname, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  character(len=*), intent(in)                         :: fname
  real(wp), allocatable, intent(out)                   :: array(:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="read", status="old", access="sequential", form="formatted")
  call self%read(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine load_real_4d

subroutine load_real_5d(self, fname, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  character(len=*), intent(in)                         :: fname
  real(wp), allocatable, intent(out)                   :: array(:,:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="read", status="old", access="sequential", form="formatted")
  call self%read(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine load_real_5d

subroutine load_cmplx(self, fname, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  character(len=*), intent(in)                         :: fname
  complex(wp), intent(out)                             :: array
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="read", status="old", access="sequential", form="formatted")
  call self%read(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine load_cmplx

subroutine load_cmplx_1d(self, fname, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  character(len=*), intent(in)                         :: fname
  complex(wp), allocatable, intent(out)                :: array(:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="read", status="old", access="sequential", form="formatted")
  call self%read(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine load_cmplx_1d

subroutine load_cmplx_2d(self, fname, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  character(len=*), intent(in)                         :: fname
  complex(wp), allocatable, intent(out)                :: array(:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="read", status="old", access="sequential", form="formatted")
  call self%read(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine load_cmplx_2d

subroutine load_cmplx_3d(self, fname, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  character(len=*), intent(in)                         :: fname
  complex(wp), allocatable, intent(out)                :: array(:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="read", status="old", access="sequential", form="formatted")
  call self%read(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine load_cmplx_3d

subroutine load_cmplx_4d(self, fname, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  character(len=*), intent(in)                         :: fname
  complex(wp), allocatable, intent(out)                :: array(:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="read", status="old", access="sequential", form="formatted")
  call self%read(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine load_cmplx_4d

subroutine load_cmplx_5d(self, fname, array, ierror, error_msg)
  class(ArrayTxtFileIO), intent(in)                    :: self 
  character(len=*), intent(in)                         :: fname
  complex(wp), allocatable, intent(out)                :: array(:,:,:,:,:)
  integer, optional, intent(out)                       :: ierror
  character(len=:), allocatable, optional, intent(out) :: error_msg
  !local
  type(FileHandle)                                     :: fh

  fh = FileHandle(fname)
  call fh%open(action="read", status="old", access="sequential", form="formatted")
  call self%read(fh%funit, array, ierror=ierror, error_msg=error_msg)
  call fh%close()
end subroutine load_cmplx_5d

end module array_txt_file_io