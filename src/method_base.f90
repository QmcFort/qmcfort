! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module method_base

use constants
use qmcfort_in, only: add_input

implicit none

private
public QmcFortMethod

type QmcFortMethod
  logical                       :: active
  character(len=:), allocatable :: method
  character(len=:), allocatable :: basis
  character(len=:), allocatable :: integral_mode
contains
  generic, public    :: assignment(=) => assign_qmcfort_method
  procedure, private :: assign_qmcfort_method

  procedure, public  :: is_active => is_qmcfort_method_active
  procedure, public  :: get_basis => get_qmcfort_method_basis
  procedure, public  :: get_integral_mode => get_qmcfort_method_integral_mode
end type QmcFortMethod

interface QmcFortMethod
  module procedure :: finit_qmcfort_method
end interface QmcFortMethod

contains

!******************************************************************************** 
!
! Initialization of the QmcFortMethod objects
!
!******************************************************************************** 
subroutine init_qmcfort_method(self, method, def_active, basis, integral_mode)
  type(QmcFortMethod), intent(out) :: self
  character(len=*), intent(in)     :: method
  logical, intent(in)              :: def_active
  character(len=*), intent(in)     :: basis
  character(len=*), intent(in)     :: integral_mode
  
  self%method = method
  self%active = def_active
  self%basis = basis
  self%integral_mode = integral_mode

  !read the method's activity from the input file
  call add_input(self%method, self%active)
end subroutine init_qmcfort_method


!******************************************************************************** 
!
! QmcFortMethod constructor
!
!******************************************************************************** 
function finit_qmcfort_method(method, def_active, basis, integral_mode) result(self)
  character(len=*), intent(in) :: method
  logical, intent(in)          :: def_active
  character(len=*), intent(in) :: basis
  character(len=*), intent(in) :: integral_mode
  type(QmcFortMethod)          :: self

  call init_qmcfort_method(self, method, def_active, basis, integral_mode)
end function finit_qmcfort_method


!******************************************************************************** 
!
! Assignment operator for the QmcfortMethod type
!
!******************************************************************************** 
subroutine assign_qmcfort_method(to, from)
  class(QmcFortMethod), intent(inout) :: to
  type(QmcFortMethod), intent(in)     :: from
  
  to%method = from%method
  to%active = from%active
  to%basis = from%basis
  to%integral_mode = from%integral_mode
end subroutine assign_qmcfort_method


!******************************************************************************** 
!
! Is QmcFortMethod active
!
!******************************************************************************** 
function is_qmcfort_method_active(self) result(active)
  class(QmcFortMethod), intent(in) :: self
  logical                          :: active

  active = self%active
end function is_qmcfort_method_active


!******************************************************************************** 
!
! Get basis variable from the QmcFortMethod object   
!
!******************************************************************************** 
function get_qmcfort_method_basis(self) result(basis)
  class(QmcFortMethod), intent(in) :: self
  character(len=:), allocatable    :: basis

  basis = self%basis
end function get_qmcfort_method_basis


!******************************************************************************** 
!
! Get integral_mode variable from the QmcFortMethod object   
!
!******************************************************************************** 
function get_qmcfort_method_integral_mode(self) result(integral_mode)
  class(QmcFortMethod), intent(in) :: self
  character(len=:), allocatable    :: integral_mode

  integral_mode = self%integral_mode
end function get_qmcfort_method_integral_mode

end module method_base
