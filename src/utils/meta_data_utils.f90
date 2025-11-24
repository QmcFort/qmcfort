! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module meta_data_utils

use constants

implicit none

private
public :: type_kind_label

interface type_kind_label
  module procedure type_kind_label_int_4, type_kind_label_int_8
  module procedure type_kind_label_real_sp, type_kind_label_real_dp
  module procedure type_kind_label_cmplx_sp, type_kind_label_cmplx_dp
end interface type_kind_label

contains

!******************************************************************************** 
!
! Return a character label indicating the type and kind of the given scalar
! argument.
!
!******************************************************************************** 
function type_kind_label_int_4(arg) result(label)
  integer(i4), intent(in)       :: arg(..)
  character(len=:), allocatable :: label

  label = "i4"
end function type_kind_label_int_4

function type_kind_label_int_8(arg) result(label)
  integer(i8), intent(in)       :: arg(..)
  character(len=:), allocatable :: label

  label = "i8"
end function type_kind_label_int_8

function type_kind_label_real_sp(arg) result(label)
  real(sp), intent(in)          :: arg(..)
  character(len=:), allocatable :: label

  label = "f4"
end function type_kind_label_real_sp

function type_kind_label_real_dp(arg) result(label)
  real(dp), intent(in)          :: arg(..)
  character(len=:), allocatable :: label

  label = "f8"
end function type_kind_label_real_dp

function type_kind_label_cmplx_sp(arg) result(label)
  complex(sp), intent(in)       :: arg(..)
  character(len=:), allocatable :: label

  label = "c8"
end function type_kind_label_cmplx_sp

function type_kind_label_cmplx_dp(arg) result(label)
  complex(dp), intent(in)       :: arg(..)
  character(len=:), allocatable :: label

  label = "c16"
end function type_kind_label_cmplx_dp

end module meta_data_utils