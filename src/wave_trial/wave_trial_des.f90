! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module wave_trial_des

use constants
use file_handle, only: FileHandle
use hamilton_vars, only: hamil_descriptor

implicit none

private
public :: WaveTrialDes

type WaveTrialDes
  character(len=:), allocatable :: exchange_mode
  integer                       :: h_block_size = unset
  integer                       :: x_block_size = unset
  logical                       :: compress_x
  real(wp)                      :: compress_x_tol

  integer                       :: nfrozen
  integer                       :: nactive
  integer                       :: nel_active(2)
  integer                       :: nell_active(2)
  real(wp)                      :: delta_x
contains
  generic            :: assignment(=) => assign_wave_trial_des
  procedure          :: report => report_wave_trial_des

  procedure, private :: assign_wave_trial_des
end type WaveTrialDes

interface WaveTrialDes
  module procedure finit_wave_trial_des
end interface WaveTrialDes

contains

!******************************************************************************** 
!
! Initialization of the WaveTrialDes type
!
!******************************************************************************** 
subroutine init_wave_trial_des(hdes, exchange_mode, h_block_size, x_block_size, compress_x, compress_x_tol, self)
  type(hamil_descriptor), intent(in) :: hdes
  character(len=*), intent(in)       :: exchange_mode
  integer, intent(in)                :: h_block_size
  integer, intent(in)                :: x_block_size
  logical, intent(in)                :: compress_x
  real(wp), intent(in)               :: compress_x_tol
  type(WaveTrialDes), intent(out)    :: self

  self%exchange_mode = exchange_mode
  self%h_block_size = h_block_size
  self%x_block_size = x_block_size
  self%compress_x = compress_x
  self%compress_x_tol = compress_x_tol

  self%delta_x = 0.0_wp

  self%nfrozen = 0
  self%nactive = hdes%n
  self%nel_active = hdes%nel
end subroutine init_wave_trial_des


!******************************************************************************** 
!
! Assignment operator for the WaveTrialDes
!
!******************************************************************************** 
subroutine assign_wave_trial_des(to, from)
  class(WaveTrialDes), intent(inout) :: to
  type(WaveTrialDes), intent(in)     :: from

  to%exchange_mode = from%exchange_mode
  to%h_block_size = from%h_block_size
  to%x_block_size = from%x_block_size
  to%compress_x = from%compress_x
  to%compress_x_tol = from%compress_x_tol

  to%nfrozen = from%nfrozen
  to%nactive = from%nactive
  to%nel_active = from%nel_active
  to%nell_active = from%nell_active
  to%delta_x = from%delta_x
end subroutine assign_wave_trial_des


!******************************************************************************** 
!
! WaveTrialDes constructor       
!
! Note:
!    The part that requires input info is initialized here,
!    remaining variables are initialized during the construction of the WaveTrial
!
!******************************************************************************** 
function finit_wave_trial_des(hdes, exchange_mode, h_block_size, x_block_size, compress_x, compress_x_tol) result(self)
  type(hamil_descriptor), intent(in) :: hdes
  character(len=*), intent(in)       :: exchange_mode
  integer, intent(in)                :: h_block_size
  integer, intent(in)                :: x_block_size
  logical, intent(in)                :: compress_x
  real(wp), intent(in)               :: compress_x_tol
  type(WaveTrialDes)                 :: self

  call init_wave_trial_des(hdes, exchange_mode, h_block_size, x_block_size, compress_x, compress_x_tol, self)
end function finit_wave_trial_des


!******************************************************************************** 
!
! Report WaveTrialDesc object
!
!******************************************************************************** 
subroutine report_wave_trial_des(self, fh)
  class(WaveTrialDes), intent(in) :: self
  type(FileHandle), intent(in)    :: fh

  write(fh%funit,*)
  write(fh%funit,*)
  write(fh%funit,*)     "Trial Wave function descriptor" 
  write(fh%funit,*)     "-------------------------------" 
  write(fh%funit,100) "number of frozen orbitals       ", self%nfrozen
  write(fh%funit,100) "number of active orbitals       ", self%nactive
  write(fh%funit,100) "number of active electrons      ", self%nel_active
  write(fh%funit,100) "block size for Hartree terms    ", self%h_block_size
  write(fh%funit,100) "block size for exchange terms   ", self%x_block_size
  write(fh%funit,102) "Exchange energy mode            ", self%exchange_mode 

  if (self%exchange_mode == "cholesky") then
    write(fh%funit,101)  "Compressed Cholesky vectors X   ", self%compress_x
    write(fh%funit,103)  "Tolerance for the compression   ", self%compress_x_tol
    write(fh%funit,103)  "Exchange energy correction      ", self%delta_x
  end if

  100 format(1x, t5, a, t50, "= ", 10i8)
  101 format(1x, t5, a, t50, "= ", 10l6)
  102 format(1x, t5, a, t50, "= ", a)
  103 format(1x, t5, a, t50, "= ", 10es12.4)
end subroutine report_wave_trial_des

end module wave_trial_des