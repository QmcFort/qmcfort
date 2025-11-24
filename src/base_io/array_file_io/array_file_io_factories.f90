! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module array_io_factories

use array_file_io, only: ArrayFileIO
use array_txt_file_io, only: ArrayTxtFileIO
use array_bin_file_io, only: ArrayBinFileIO
use array_numpy_file_io, only: ArrayNumpyFileIO

implicit none

enum, bind(c)
  enumerator :: NUMPY = 0, TXT, BIN
end enum

private
public :: array_io_factory

contains

function select_file_format(file_format_request) result(file_format)
  character(len=*), intent(in) :: file_format_request
  integer                      :: file_format

  select case (file_format_request)
    case ("text", "txt")
      file_format = TXT
    case ("binary", "bin")
      file_format = BIN
    case ("numpy", "np")
      file_format = NUMPY
    case default
      file_format = TXT
  end select
end function select_file_format

subroutine array_io_factory(file_format_request, array_io)
  character(len=*), intent(in)        :: file_format_request
  class(ArrayFileIO), allocatable     :: array_io
  !local
  integer                             :: file_format
  type(ArrayTxtFileIO), allocatable   :: array_io_txt
  type(ArrayBinFileIO), allocatable   :: array_io_bin
  type(ArrayNumpyFileIO), allocatable :: array_io_numpy

  file_format = select_file_format(file_format_request)

  select case (file_format)
    case (TXT)
      allocate(ArrayTxtFileIO :: array_io_txt)
      call move_alloc(array_io_txt, array_io)
    case (BIN)
      allocate(ArrayBinFileIO :: array_io_bin)
      call move_alloc(array_io_bin, array_io)
    case (NUMPY)
      allocate(ArrayNumpyFileIO :: array_io_numpy)
      call move_alloc(array_io_numpy, array_io)
  end select
end subroutine array_io_factory

end module array_io_factories