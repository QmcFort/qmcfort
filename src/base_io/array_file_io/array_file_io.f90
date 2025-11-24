! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module array_file_io

use constants

implicit none

private
public :: ArrayFileIO

type, abstract :: ArrayFileIO
contains
  generic, public     :: write => write_int, write_int_1d, write_int_2d, write_int_3d, write_int_4d, write_int_5d, &
                                  write_real, write_real_1d, write_real_2d, write_real_3d, write_real_4d, write_real_5d, &
                                  write_cmplx, write_cmplx_1d, write_cmplx_2d, write_cmplx_3d, write_cmplx_4d, write_cmplx_5d

  generic, public     :: read => read_int, read_int_1d, read_int_2d, read_int_3d, read_int_4d, read_int_5d, &
                                 read_real, read_real_1d, read_real_2d, read_real_3d, read_real_4d, read_real_5d, &
                                 read_cmplx, read_cmplx_1d, read_cmplx_2d, read_cmplx_3d, read_cmplx_4d, read_cmplx_5d

  generic, public     :: save => save_int, save_int_1d, save_int_2d, save_int_3d, save_int_4d, save_int_5d, &
                                 save_real, save_real_1d, save_real_2d, save_real_3d, save_real_4d, save_real_5d, &
                                 save_cmplx, save_cmplx_1d, save_cmplx_2d, save_cmplx_3d, save_cmplx_4d, save_cmplx_5d

  generic, public     :: load => load_int, load_int_1d, load_int_2d, load_int_3d, load_int_4d, load_int_5d, &
                                 load_real, load_real_1d, load_real_2d, load_real_3d, load_real_4d, load_real_5d, &
                                 load_cmplx, load_cmplx_1d, load_cmplx_2d, load_cmplx_3d, load_cmplx_4d, load_cmplx_5d

  procedure(Iwrite_int), deferred      :: write_int
  procedure(Iwrite_int_1d), deferred   :: write_int_1d
  procedure(Iwrite_int_2d), deferred   :: write_int_2d
  procedure(Iwrite_int_3d), deferred   :: write_int_3d
  procedure(Iwrite_int_4d), deferred   :: write_int_4d
  procedure(Iwrite_int_5d), deferred   :: write_int_5d
  procedure(Iwrite_real), deferred     :: write_real
  procedure(Iwrite_real_1d), deferred  :: write_real_1d
  procedure(Iwrite_real_2d), deferred  :: write_real_2d
  procedure(Iwrite_real_3d), deferred  :: write_real_3d
  procedure(Iwrite_real_4d), deferred  :: write_real_4d
  procedure(Iwrite_real_5d), deferred  :: write_real_5d
  procedure(Iwrite_cmplx), deferred    :: write_cmplx
  procedure(Iwrite_cmplx_1d), deferred :: write_cmplx_1d
  procedure(Iwrite_cmplx_2d), deferred :: write_cmplx_2d
  procedure(Iwrite_cmplx_3d), deferred :: write_cmplx_3d
  procedure(Iwrite_cmplx_4d), deferred :: write_cmplx_4d
  procedure(Iwrite_cmplx_5d), deferred :: write_cmplx_5d

  procedure(Iread_int), deferred      :: read_int
  procedure(Iread_int_1d), deferred   :: read_int_1d
  procedure(Iread_int_2d), deferred   :: read_int_2d
  procedure(Iread_int_3d), deferred   :: read_int_3d
  procedure(Iread_int_4d), deferred   :: read_int_4d
  procedure(Iread_int_5d), deferred   :: read_int_5d
  procedure(Iread_real), deferred     :: read_real
  procedure(Iread_real_1d), deferred  :: read_real_1d
  procedure(Iread_real_2d), deferred  :: read_real_2d
  procedure(Iread_real_3d), deferred  :: read_real_3d
  procedure(Iread_real_4d), deferred  :: read_real_4d
  procedure(Iread_real_5d), deferred  :: read_real_5d
  procedure(Iread_cmplx), deferred    :: read_cmplx
  procedure(Iread_cmplx_1d), deferred :: read_cmplx_1d
  procedure(Iread_cmplx_2d), deferred :: read_cmplx_2d
  procedure(Iread_cmplx_3d), deferred :: read_cmplx_3d
  procedure(Iread_cmplx_4d), deferred :: read_cmplx_4d
  procedure(Iread_cmplx_5d), deferred :: read_cmplx_5d

  procedure(Isave_int), deferred      :: save_int
  procedure(Isave_int_1d), deferred   :: save_int_1d
  procedure(Isave_int_2d), deferred   :: save_int_2d
  procedure(Isave_int_3d), deferred   :: save_int_3d
  procedure(Isave_int_4d), deferred   :: save_int_4d
  procedure(Isave_int_5d), deferred   :: save_int_5d
  procedure(Isave_real), deferred     :: save_real
  procedure(Isave_real_1d), deferred  :: save_real_1d
  procedure(Isave_real_2d), deferred  :: save_real_2d
  procedure(Isave_real_3d), deferred  :: save_real_3d
  procedure(Isave_real_4d), deferred  :: save_real_4d
  procedure(Isave_real_5d), deferred  :: save_real_5d
  procedure(Isave_cmplx), deferred    :: save_cmplx
  procedure(Isave_cmplx_1d), deferred :: save_cmplx_1d
  procedure(Isave_cmplx_2d), deferred :: save_cmplx_2d
  procedure(Isave_cmplx_3d), deferred :: save_cmplx_3d
  procedure(Isave_cmplx_4d), deferred :: save_cmplx_4d
  procedure(Isave_cmplx_5d), deferred :: save_cmplx_5d

  procedure(Iload_int), deferred      :: load_int
  procedure(Iload_int_1d), deferred   :: load_int_1d
  procedure(Iload_int_2d), deferred   :: load_int_2d
  procedure(Iload_int_3d), deferred   :: load_int_3d
  procedure(Iload_int_4d), deferred   :: load_int_4d
  procedure(Iload_int_5d), deferred   :: load_int_5d
  procedure(Iload_real), deferred     :: load_real
  procedure(Iload_real_1d), deferred  :: load_real_1d
  procedure(Iload_real_2d), deferred  :: load_real_2d
  procedure(Iload_real_3d), deferred  :: load_real_3d
  procedure(Iload_real_4d), deferred  :: load_real_4d
  procedure(Iload_real_5d), deferred  :: load_real_5d
  procedure(Iload_cmplx), deferred    :: load_cmplx
  procedure(Iload_cmplx_1d), deferred :: load_cmplx_1d
  procedure(Iload_cmplx_2d), deferred :: load_cmplx_2d
  procedure(Iload_cmplx_3d), deferred :: load_cmplx_3d
  procedure(Iload_cmplx_4d), deferred :: load_cmplx_4d
  procedure(Iload_cmplx_5d), deferred :: load_cmplx_5d
end type ArrayFileIO

abstract interface
  subroutine Iwrite_int(self, funit, array, ierror, error_msg)
    import ArrayFileIO
    class(ArrayFileIO), intent(in)                       :: self
    integer, intent(in)                                  :: funit
    integer, intent(in)                                  :: array
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iwrite_int

  subroutine Iwrite_int_1d(self, funit, array, ierror, error_msg)
    import ArrayFileIO
    class(ArrayFileIO), intent(in)                       :: self
    integer, intent(in)                                  :: funit
    integer, intent(in)                                  :: array(:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iwrite_int_1d

  subroutine Iwrite_int_2d(self, funit, array, ierror, error_msg)
    import ArrayFileIO
    class(ArrayFileIO), intent(in)                       :: self
    integer, intent(in)                                  :: funit
    integer, intent(in)                                  :: array(:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iwrite_int_2d

  subroutine Iwrite_int_3d(self, funit, array, ierror, error_msg)
    import ArrayFileIO
    class(ArrayFileIO), intent(in)                       :: self
    integer, intent(in)                                  :: funit
    integer, intent(in)                                  :: array(:,:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iwrite_int_3d

  subroutine Iwrite_int_4d(self, funit, array, ierror, error_msg)
    import ArrayFileIO
    class(ArrayFileIO), intent(in)                       :: self
    integer, intent(in)                                  :: funit
    integer, intent(in)                                  :: array(:,:,:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iwrite_int_4d

  subroutine Iwrite_int_5d(self, funit, array, ierror, error_msg)
    import ArrayFileIO
    class(ArrayFileIO), intent(in)                       :: self
    integer, intent(in)                                  :: funit
    integer, intent(in)                                  :: array(:,:,:,:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iwrite_int_5d

  subroutine Iwrite_real(self, funit, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    integer, intent(in)                                  :: funit
    real(wp), intent(in)                                 :: array
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iwrite_real

  subroutine Iwrite_real_1d(self, funit, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    integer, intent(in)                                  :: funit
    real(wp), intent(in)                                 :: array(:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iwrite_real_1d

  subroutine Iwrite_real_2d(self, funit, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    integer, intent(in)                                  :: funit
    real(wp), intent(in)                                 :: array(:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iwrite_real_2d

  subroutine Iwrite_real_3d(self, funit, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    integer, intent(in)                                  :: funit
    real(wp), intent(in)                                 :: array(:,:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iwrite_real_3d

  subroutine Iwrite_real_4d(self, funit, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    integer, intent(in)                                  :: funit
    real(wp), intent(in)                                 :: array(:,:,:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iwrite_real_4d

  subroutine Iwrite_real_5d(self, funit, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    integer, intent(in)                                  :: funit
    real(wp), intent(in)                                 :: array(:,:,:,:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iwrite_real_5d

  subroutine Iwrite_cmplx(self, funit, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    integer, intent(in)                                  :: funit
    complex(wp), intent(in)                              :: array
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iwrite_cmplx

  subroutine Iwrite_cmplx_1d(self, funit, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    integer, intent(in)                                  :: funit
    complex(wp), intent(in)                              :: array(:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iwrite_cmplx_1d

  subroutine Iwrite_cmplx_2d(self, funit, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    integer, intent(in)                                  :: funit
    complex(wp), intent(in)                              :: array(:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iwrite_cmplx_2d

  subroutine Iwrite_cmplx_3d(self, funit, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    integer, intent(in)                                  :: funit
    complex(wp), intent(in)                              :: array(:,:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iwrite_cmplx_3d

  subroutine Iwrite_cmplx_4d(self, funit, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    integer, intent(in)                                  :: funit
    complex(wp), intent(in)                              :: array(:,:,:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iwrite_cmplx_4d

  subroutine Iwrite_cmplx_5d(self, funit, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    integer, intent(in)                                  :: funit
    complex(wp), intent(in)                              :: array(:,:,:,:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iwrite_cmplx_5d

  subroutine Iread_int(self, funit, array, ierror, error_msg)
    import ArrayFileIO
    class(ArrayFileIO), intent(in)                       :: self
    integer, intent(in)                                  :: funit
    integer, intent(out)                                 :: array
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iread_int

  subroutine Iread_int_1d(self, funit, array, ierror, error_msg)
    import ArrayFileIO
    class(ArrayFileIO), intent(in)                       :: self
    integer, intent(in)                                  :: funit
    integer, allocatable, intent(out)                    :: array(:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iread_int_1d

  subroutine Iread_int_2d(self, funit, array, ierror, error_msg)
    import ArrayFileIO
    class(ArrayFileIO), intent(in)                       :: self
    integer, intent(in)                                  :: funit
    integer, allocatable, intent(out)                    :: array(:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iread_int_2d

  subroutine Iread_int_3d(self, funit, array, ierror, error_msg)
    import ArrayFileIO
    class(ArrayFileIO), intent(in)                       :: self
    integer, intent(in)                                  :: funit
    integer, allocatable, intent(out)                    :: array(:,:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iread_int_3d

  subroutine Iread_int_4d(self, funit, array, ierror, error_msg)
    import ArrayFileIO
    class(ArrayFileIO), intent(in)                       :: self
    integer, intent(in)                                  :: funit
    integer, allocatable, intent(out)                    :: array(:,:,:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iread_int_4d

  subroutine Iread_int_5d(self, funit, array, ierror, error_msg)
    import ArrayFileIO
    class(ArrayFileIO), intent(in)                       :: self
    integer, intent(in)                                  :: funit
    integer, allocatable, intent(out)                    :: array(:,:,:,:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iread_int_5d

  subroutine Iread_real(self, funit, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    integer, intent(in)                                  :: funit
    real(wp), intent(out)                                :: array
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iread_real

  subroutine Iread_real_1d(self, funit, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                      :: self
    integer, intent(in)                                 :: funit
    real(wp), allocatable, intent(out)                  :: array(:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iread_real_1d

  subroutine Iread_real_2d(self, funit, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    integer, intent(in)                                  :: funit
    real(wp), allocatable, intent(out)                   :: array(:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iread_real_2d

  subroutine Iread_real_3d(self, funit, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    integer, intent(in)                                  :: funit
    real(wp), allocatable, intent(out)                   :: array(:,:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iread_real_3d

  subroutine Iread_real_4d(self, funit, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    integer, intent(in)                                  :: funit
    real(wp), allocatable, intent(out)                   :: array(:,:,:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iread_real_4d

  subroutine Iread_real_5d(self, funit, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    integer, intent(in)                                  :: funit
    real(wp), allocatable, intent(out)                   :: array(:,:,:,:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iread_real_5d

  subroutine Iread_cmplx(self, funit, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    integer, intent(in)                                  :: funit
    complex(wp), intent(out)                             :: array
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iread_cmplx

  subroutine Iread_cmplx_1d(self, funit, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    integer, intent(in)                                  :: funit
    complex(wp), allocatable, intent(out)                :: array(:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iread_cmplx_1d

  subroutine Iread_cmplx_2d(self, funit, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    integer, intent(in)                                  :: funit
    complex(wp), allocatable, intent(out)                :: array(:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iread_cmplx_2d

  subroutine Iread_cmplx_3d(self, funit, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    integer, intent(in)                                  :: funit
    complex(wp), allocatable, intent(out)                :: array(:,:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iread_cmplx_3d

  subroutine Iread_cmplx_4d(self, funit, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    integer, intent(in)                                  :: funit
    complex(wp), allocatable, intent(out)                :: array(:,:,:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iread_cmplx_4d

  subroutine Iread_cmplx_5d(self, funit, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    integer, intent(in)                                  :: funit
    complex(wp), allocatable, intent(out)                :: array(:,:,:,:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iread_cmplx_5d

  subroutine Isave_int(self, fname, array, ierror, error_msg)
    import ArrayFileIO
    class(ArrayFileIO), intent(in)                       :: self
    character(len=*), intent(in)                         :: fname
    integer, intent(in)                                  :: array
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Isave_int

  subroutine Isave_int_1d(self, fname, array, ierror, error_msg)
    import ArrayFileIO
    class(ArrayFileIO), intent(in)                       :: self
    character(len=*), intent(in)                         :: fname
    integer, intent(in)                                  :: array(:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Isave_int_1d

  subroutine Isave_int_2d(self, fname, array, ierror, error_msg)
    import ArrayFileIO
    class(ArrayFileIO), intent(in)                       :: self
    character(len=*), intent(in)                         :: fname
    integer, intent(in)                                  :: array(:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Isave_int_2d

  subroutine Isave_int_3d(self, fname, array, ierror, error_msg)
    import ArrayFileIO
    class(ArrayFileIO), intent(in)                       :: self
    character(len=*), intent(in)                         :: fname
    integer, intent(in)                                  :: array(:,:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Isave_int_3d

  subroutine Isave_int_4d(self, fname, array, ierror, error_msg)
    import ArrayFileIO
    class(ArrayFileIO), intent(in)                       :: self
    character(len=*), intent(in)                         :: fname
    integer, intent(in)                                  :: array(:,:,:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Isave_int_4d

  subroutine Isave_int_5d(self, fname, array, ierror, error_msg)
    import ArrayFileIO
    class(ArrayFileIO), intent(in)                       :: self
    character(len=*), intent(in)                         :: fname
    integer, intent(in)                                  :: array(:,:,:,:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Isave_int_5d

  subroutine Isave_real(self, fname, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    character(len=*), intent(in)                         :: fname
    real(wp), intent(in)                                 :: array
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Isave_real

  subroutine Isave_real_1d(self, fname, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    character(len=*), intent(in)                         :: fname
    real(wp), intent(in)                                 :: array(:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Isave_real_1d

  subroutine Isave_real_2d(self, fname, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    character(len=*), intent(in)                         :: fname
    real(wp), intent(in)                                 :: array(:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Isave_real_2d

  subroutine Isave_real_3d(self, fname, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    character(len=*), intent(in)                         :: fname
    real(wp), intent(in)                                 :: array(:,:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Isave_real_3d

  subroutine Isave_real_4d(self, fname, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    character(len=*), intent(in)                         :: fname
    real(wp), intent(in)                                 :: array(:,:,:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Isave_real_4d

  subroutine Isave_real_5d(self, fname, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    character(len=*), intent(in)                         :: fname
    real(wp), intent(in)                                 :: array(:,:,:,:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Isave_real_5d

  subroutine Isave_cmplx(self, fname, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    character(len=*), intent(in)                         :: fname
    complex(wp), intent(in)                              :: array
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Isave_cmplx

  subroutine Isave_cmplx_1d(self, fname, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    character(len=*), intent(in)                         :: fname
    complex(wp), intent(in)                              :: array(:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Isave_cmplx_1d

  subroutine Isave_cmplx_2d(self, fname, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    character(len=*), intent(in)                         :: fname
    complex(wp), intent(in)                              :: array(:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Isave_cmplx_2d

  subroutine Isave_cmplx_3d(self, fname, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    character(len=*), intent(in)                         :: fname
    complex(wp), intent(in)                              :: array(:,:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Isave_cmplx_3d

  subroutine Isave_cmplx_4d(self, fname, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    character(len=*), intent(in)                         :: fname
    complex(wp), intent(in)                              :: array(:,:,:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Isave_cmplx_4d

  subroutine Isave_cmplx_5d(self, fname, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    character(len=*), intent(in)                         :: fname
    complex(wp), intent(in)                              :: array(:,:,:,:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Isave_cmplx_5d

  subroutine Iload_int(self, fname, array, ierror, error_msg)
    import ArrayFileIO
    class(ArrayFileIO), intent(in)                       :: self
    character(len=*), intent(in)                         :: fname
    integer, intent(out)                                 :: array
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iload_int

  subroutine Iload_int_1d(self, fname, array, ierror, error_msg)
    import ArrayFileIO
    class(ArrayFileIO), intent(in)                       :: self
    character(len=*), intent(in)                         :: fname
    integer, allocatable, intent(out)                    :: array(:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iload_int_1d

  subroutine Iload_int_2d(self, fname, array, ierror, error_msg)
    import ArrayFileIO
    class(ArrayFileIO), intent(in)                       :: self
    character(len=*), intent(in)                         :: fname
    integer, allocatable, intent(out)                    :: array(:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iload_int_2d

  subroutine Iload_int_3d(self, fname, array, ierror, error_msg)
    import ArrayFileIO
    class(ArrayFileIO), intent(in)                       :: self
    character(len=*), intent(in)                         :: fname
    integer, allocatable, intent(out)                    :: array(:,:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iload_int_3d

  subroutine Iload_int_4d(self, fname, array, ierror, error_msg)
    import ArrayFileIO
    class(ArrayFileIO), intent(in)                       :: self
    character(len=*), intent(in)                         :: fname
    integer, allocatable, intent(out)                    :: array(:,:,:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iload_int_4d

  subroutine Iload_int_5d(self, fname, array, ierror, error_msg)
    import ArrayFileIO
    class(ArrayFileIO), intent(in)                       :: self
    character(len=*), intent(in)                         :: fname
    integer, allocatable, intent(out)                    :: array(:,:,:,:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iload_int_5d

  subroutine Iload_real(self, fname, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    character(len=*), intent(in)                         :: fname
    real(wp), intent(out)                                :: array
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iload_real

  subroutine Iload_real_1d(self, fname, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    character(len=*), intent(in)                         :: fname
    real(wp), allocatable, intent(out)                   :: array(:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iload_real_1d

  subroutine Iload_real_2d(self, fname, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    character(len=*), intent(in)                         :: fname
    real(wp), allocatable, intent(out)                   :: array(:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iload_real_2d

  subroutine Iload_real_3d(self, fname, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    character(len=*), intent(in)                         :: fname
    real(wp), allocatable, intent(out)                   :: array(:,:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iload_real_3d

  subroutine Iload_real_4d(self, fname, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    character(len=*), intent(in)                         :: fname
    real(wp), allocatable, intent(out)                   :: array(:,:,:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iload_real_4d

  subroutine Iload_real_5d(self, fname, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    character(len=*), intent(in)                         :: fname
    real(wp), allocatable, intent(out)                   :: array(:,:,:,:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iload_real_5d

  subroutine Iload_cmplx(self, fname, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    character(len=*), intent(in)                         :: fname
    complex(wp), intent(out)                             :: array
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iload_cmplx

  subroutine Iload_cmplx_1d(self, fname, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    character(len=*), intent(in)                         :: fname
    complex(wp), allocatable, intent(out)                :: array(:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iload_cmplx_1d

  subroutine Iload_cmplx_2d(self, fname, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    character(len=*), intent(in)                         :: fname
    complex(wp), allocatable, intent(out)                :: array(:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iload_cmplx_2d

  subroutine Iload_cmplx_3d(self, fname, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    character(len=*), intent(in)                         :: fname
    complex(wp), allocatable, intent(out)                :: array(:,:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iload_cmplx_3d

  subroutine Iload_cmplx_4d(self, fname, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    character(len=*), intent(in)                         :: fname
    complex(wp), allocatable, intent(out)                :: array(:,:,:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iload_cmplx_4d

  subroutine Iload_cmplx_5d(self, fname, array, ierror, error_msg)
    import ArrayFileIO, wp
    class(ArrayFileIO), intent(in)                       :: self
    character(len=*), intent(in)                         :: fname
    complex(wp), allocatable, intent(out)                :: array(:,:,:,:,:)
    integer, optional, intent(out)                       :: ierror
    character(len=:), allocatable, optional, intent(out) :: error_msg
  end subroutine Iload_cmplx_5d
end interface

end module array_file_io