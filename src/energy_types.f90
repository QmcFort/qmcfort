! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module energy_types

#include "preproc.inc"
use constants
use mpi
use file_handle, only: FileHandle
use statistics
use qmcfort_io

implicit none

public
!public Energy, CEnergy, EnergyStd, read_c_energies, write_c_energies

!length of energy array == [e%e1, e%eh, e%ex]
integer, parameter :: en_len = 3

type Energy
  real(wp) :: enuc
  real(wp) :: e0

  real(wp) :: e1
  real(wp) :: eh
  real(wp) :: ex
contains
  generic            :: assignment(=)   => assign_energy_to_energy, assign_c_energy_to_energy
  generic            :: operator(+)     => add_energy
  generic            :: operator(*)     => energy_left_scalar, energy_right_scalar
  generic            :: operator(/)     => energy_scalar_division
  generic            :: report          => report_energy, report_std_energy
  procedure          :: reset           => reset_energy
  procedure          :: to_array        => energy_to_array
  procedure          :: electron_energy => get_electron_energy 
  procedure          :: total_energy    => get_total_energy
  procedure          :: mpisum          => energy_mpisum
  procedure, private :: assign_energy_to_energy, assign_c_energy_to_energy
  procedure, private :: add_energy, energy_right_scalar, energy_scalar_division
  procedure, private :: report_energy, report_std_energy
  procedure, pass(self), private :: energy_left_scalar
end type Energy

type CEnergy
  real(wp)    :: enuc
  real(wp)    :: e0

  complex(wp) :: e1
  complex(wp) :: eh
  complex(wp) :: ex
contains
  generic            :: assignment(=)   => assign_c_energy_to_c_energy, assign_energy_to_c_energy
  generic            :: operator(+)     => add_c_energy
  generic            :: operator(*)     => c_energy_real_left_scalar, c_energy_real_right_scalar, c_energy_cmplx_left_scalar, c_energy_cmplx_right_scalar
  generic            :: operator(/)     => c_energy_real_scalar_division, c_energy_cmplx_scalar_division
  generic            :: report          => report_c_energy, report_c_std_energy
  procedure          :: reset           => reset_c_energy
  procedure          :: to_array        => c_energy_to_array
  procedure          :: electron_energy => get_electron_c_energy
  procedure          :: total_energy    => get_total_c_energy
  procedure          :: mpisum          => c_energy_mpisum
  procedure, private :: assign_c_energy_to_c_energy, assign_energy_to_c_energy
  procedure, private :: add_c_energy, c_energy_real_right_scalar, c_energy_cmplx_right_scalar
  procedure, private :: c_energy_real_scalar_division, c_energy_cmplx_scalar_division
  procedure, private :: report_c_energy, report_c_std_energy
  procedure, pass(self), private :: c_energy_real_left_scalar, c_energy_cmplx_left_scalar
end type CEnergy

type EnergyStd
  real(wp) :: e1
  real(wp) :: eh
  real(wp) :: ex
  real(wp) :: e

  real(wp) :: e1_corr_len
  real(wp) :: eh_corr_len
  real(wp) :: ex_corr_len
  real(wp) :: e_corr_len
end type EnergyStd

type CorrStddev
  real(wp) :: stddev
  real(wp) :: corr_len
end type CorrStddev

interface Energy
  module procedure init_energy
end interface Energy

interface CEnergy
  module procedure init_c_energy
end interface CEnergy

interface EnergyStd
  module procedure init_energy_std, init_energy_std_val
end interface EnergyStd

contains

!******************************************************************************** 
!
! Energy type constructor
!
!******************************************************************************** 
function init_energy(enuc, e0, energy_array) result(self)
  real(wp), optional, intent(in) :: enuc
  real(wp), optional, intent(in) :: e0
  real(wp), optional, intent(in) :: energy_array(:)
  type(Energy)                   :: self

  self%enuc = zeror
  if (present(enuc)) self%enuc = enuc

  self%e0 = zeror
  if (present(e0)) self%e0 = e0

  self%e1 = zeror
  self%eh = zeror
  self%ex = zeror
  if (present(energy_array)) then
    self%e1 = energy_array(1)
    self%eh = energy_array(2)
    self%ex = energy_array(3)
  end if
end function init_energy


!******************************************************************************** 
!
! Reset variable part of the Energy type
!
!******************************************************************************** 
subroutine reset_energy(self)
  class(Energy), intent(inout) :: self

  self%e1 = zeror
  self%eh = zeror
  self%ex = zeror
end subroutine reset_energy


!******************************************************************************** 
!
! Assignment operator:  Energy = Energy
!
!******************************************************************************** 
subroutine assign_energy_to_energy(to, from)
  class(Energy), intent(inout) :: to
  class(Energy), intent(in)    :: from
  
  to%enuc = from%enuc
  to%e0 = from%e0

  to%e1 = from%e1
  to%eh = from%eh
  to%ex = from%ex
end subroutine assign_energy_to_energy


!******************************************************************************** 
!
! Assignment operator: Energy = CEnergy
!
!******************************************************************************** 
subroutine assign_c_energy_to_energy(to, from)
  class(Energy), intent(inout) :: to
  class(CEnergy), intent(in)   :: from
  
  to%enuc = from%enuc
  to%e0 = from%e0

  to%e1 = real(from%e1, wp)
  to%eh = real(from%eh, wp)
  to%ex = real(from%ex, wp)
end subroutine assign_c_energy_to_energy


!******************************************************************************** 
!
! Add two Energy types (only e1 eh ex e etot are modified)
!
!******************************************************************************** 
function add_energy(e1, e2) result(e)
  class(Energy), intent(in) :: e1
  class(Energy), intent(in) :: e2
  type(Energy)              :: e
  
  e = e1

  e%e1 = e1%e1 + e2%e1 
  e%eh = e1%eh + e2%eh 
  e%ex = e1%ex + e2%ex 
end function add_energy


!******************************************************************************** 
!
! Left scalar multiplication of the Energy type
!
!******************************************************************************** 
function energy_left_scalar(x, self) result(en)
  class(energy), intent(in) :: self
  real(wp), intent(in)      :: x
  type(energy)              :: en
  
  en = self

  en%e1 = x * self%e1
  en%eh = x * self%eh
  en%ex = x * self%ex
end function energy_left_scalar


!******************************************************************************** 
!
! Right scalar multiplication of the Energy type
!
!******************************************************************************** 
function energy_right_scalar(self, x) result(en)
  class(energy), intent(in) :: self
  real(wp), intent(in)      :: x
  type(energy)              :: en
  
  en = self

  en%e1 = x * self%e1
  en%eh = x * self%eh
  en%ex = x * self%ex
end function energy_right_scalar


!******************************************************************************** 
!
! Division of the Energy type by a scalar value
!
!******************************************************************************** 
function energy_scalar_division(self, x) result(en)
  class(Energy), intent(in) :: self
  real(wp), intent(in)      :: x
  type(Energy)              :: en
  
  en = self

  en%e1 = self%e1 / x
  en%eh = self%eh / x
  en%ex = self%ex / x
end function energy_scalar_division


!******************************************************************************** 
!
! Extract e1 eh ex from the Energy type to the array
!
!******************************************************************************** 
function energy_to_array(self) result(energy_array)
  class(Energy), intent(in) :: self
  real(wp)                  :: energy_array(en_len)
  
  energy_array(1) = self%e1
  energy_array(2) = self%eh
  energy_array(3) = self%ex
end function energy_to_array


!******************************************************************************** 
!
! Compute electronic Energy: e0 + e1 + eh + ex
!
!******************************************************************************** 
elemental function get_electron_energy(self) result(el_en)
  class(Energy), intent(in) :: self
  real(wp)                  :: el_en

  el_en = self%e0 + self%e1 + self%eh + self%ex
end function get_electron_energy


!******************************************************************************** 
!
! Compute total Energy: enuc + e0 + e1 + eh + ex
!
!******************************************************************************** 
elemental function get_total_energy(self) result(tot_en)
  class(Energy), intent(in) :: self
  real(wp)                  :: tot_en

  tot_en = self%enuc + self%e0 + self%e1 + self%eh + self%ex
end function get_total_energy


!******************************************************************************** 
!
! MPI sum of the Energy type
!
!******************************************************************************** 
subroutine energy_mpisum(self, comm)
  class(Energy), intent(inout)       :: self
  type(mpi_communicator), intent(in) :: comm

  call comm%mpisum(self%e1)
  call comm%mpisum(self%eh)
  call comm%mpisum(self%ex)
end subroutine energy_mpisum


!******************************************************************************** 
!
! report the energy type
!
!******************************************************************************** 
subroutine report_energy(self, fh, method)
  class(Energy), intent(in)              :: self
  type(FileHandle), intent(in)           :: fh
  character(len=*), optional, intent(in) :: method
  !local
  character(len=:), allocatable          :: method_
  
  if (present(method)) then
    method_ = trim(method) // "_"
  else
    method_ = ""
  end if
  
  write(fh%funit,*)
  write(fh%funit,*)   method_ // "energy of the electronic-ionic system:"
  write(fh%funit,*)   "---------------------------------------------------------"
  write(fh%funit,100)
  write(fh%funit,100) "contribution", "     Re E      "
  write(fh%funit,100) "------------", "     ----      "
  write(fh%funit,101) method_ // "enuc =", self%enuc
  write(fh%funit,101) method_ // "e0   =", self%e0
  write(fh%funit,101) method_ // "e1   =", self%e1
  write(fh%funit,101) method_ // "eh   =", self%eh
  write(fh%funit,101) method_ // "ex   =", self%ex
  write(fh%funit,101) method_ // "e    =", self%electron_energy()
  write(fh%funit,101) method_ // "etot =", self%total_energy()

  100 format (1x,a,t30,a)
  101 format (1x,a,t30,f14.8)
end subroutine report_energy


!******************************************************************************** 
!
! Report the Energy type with standard deviation
!
!******************************************************************************** 
subroutine report_std_energy(self, dself, fh, method)
  class(Energy), intent(in)              :: self
  type(EnergyStd), intent(in)            :: dself
  type(FileHandle), intent(in)           :: fh
  character(len=*), optional, intent(in) :: method
  !local variables
  real(wp)                               :: el_en, tot_en
  character(len=:), allocatable          :: method_
  
  el_en = self%electron_energy()
  tot_en = self%total_energy()

  if (present(method)) then
    method_ = trim(method) // "_"
  else
    method_ = ""
  end if

  write(fh%funit,*)
  write(fh%funit,*)   method_ // "energy of the electronic-ionic system:"
  write(fh%funit,*)   "---------------------------------------------------------"
  write(fh%funit,100) "contribution", "     Re E      ", "    stddev    ", "corr. length", " stddev final "
  write(fh%funit,100) "------------", "     ----      ", "    ------    ", "------------", " ------------ "
  write(fh%funit,101) method_ // "enuc =", self%enuc
  write(fh%funit,101) method_ // "e0   =", self%e0
  write(fh%funit,101) method_ // "e1   =", self%e1,          dself%e1,   dself%e1_corr_len, dself%e1*sqrt(dself%e1_corr_len)
  write(fh%funit,101) method_ // "eh   =", self%eh,          dself%eh,   dself%eh_corr_len, dself%eh*sqrt(dself%eh_corr_len)
  write(fh%funit,101) method_ // "ex   =", self%ex,          dself%ex,   dself%ex_corr_len, dself%ex*sqrt(dself%ex_corr_len)
  write(fh%funit,101) method_ // "e    =", el_en,            dself%e,    dself%e_corr_len,  dself%e*sqrt(dself%e_corr_len)
  write(fh%funit,101) method_ // "etot =", tot_en,           dself%e,    dself%e_corr_len,  dself%e*sqrt(dself%e_corr_len)
  
  100 format (1x,a,t30,4(a,2x))
  101 format (1x,a,t30,f14.8,2x,es14.6,2x,2x,f8.2,2x,2x,es14.6)
end subroutine report_std_energy


!******************************************************************************** 
!
! CEnergy type constructor
!
!******************************************************************************** 
function init_c_energy(enuc, e0, energy_array) result(self)
  real(wp), optional, intent(in)    :: enuc
  real(wp), optional, intent(in)    :: e0
  complex(wp), optional, intent(in) :: energy_array(:)
  type(CEnergy)                     :: self

  self%enuc = zeror
  if (present(enuc)) self%enuc = enuc

  self%e0 = zeror
  if (present(e0)) self%e0 = e0

  self%e1 = zeroc
  self%eh = zeroc
  self%ex = zeroc
  if (present(energy_array)) then
    self%e1 = energy_array(1)
    self%eh = energy_array(2)
    self%ex = energy_array(3)
  end if
end function init_c_energy


!******************************************************************************** 
!
! Reset variable part of the CEnergy type
!
!******************************************************************************** 
subroutine reset_c_energy(self)
  class(CEnergy), intent(inout) :: self

  self%e1 = zeroc
  self%eh = zeroc
  self%ex = zeroc
end subroutine reset_c_energy


!******************************************************************************** 
!
! Assignment operator: CEnergy = CEnergy
!
!******************************************************************************** 
subroutine assign_c_energy_to_c_energy(to, from)
  class(CEnergy), intent(inout) :: to
  class(CEnergy), intent(in)    :: from
  
  to%enuc = from%enuc
  to%e0 = from%e0

  to%e1 = from%e1
  to%eh = from%eh
  to%ex = from%ex
end subroutine assign_c_energy_to_c_energy


!******************************************************************************** 
!
! Assignment operator: CEnergy = Energy
!
!******************************************************************************** 
subroutine assign_energy_to_c_energy(to, from)
  class(CEnergy), intent(inout) :: to
  class(Energy), intent(in)     :: from
  
  to%enuc = from%enuc
  to%e0 = from%e0

  to%e1 = from%e1
  to%eh = from%eh
  to%ex = from%ex
end subroutine assign_energy_to_c_energy


!******************************************************************************** 
!
! Add two CEnergy types (only e1 eh ex e etot are modified)
!
!******************************************************************************** 
function add_c_energy(en1, en2) result(en)
  class(CEnergy), intent(in) :: en1
  class(CEnergy), intent(in) :: en2
  type(CEnergy)              :: en
  
  en = en1

  en%e1 = en1%e1 + en2%e1 
  en%eh = en1%eh + en2%eh 
  en%ex = en1%ex + en2%ex 
end function add_c_energy


!******************************************************************************** 
!
! Real left scalar multiplication of the CEnergy type
!
!******************************************************************************** 
function c_energy_real_left_scalar(x, self) result(en)
  real(wp), intent(in)       :: x
  class(CEnergy), intent(in) :: self
  type(CEnergy)              :: en
  
  en = self

  en%e1 = x * self%e1
  en%eh = x * self%eh
  en%ex = x * self%ex
end function c_energy_real_left_scalar


!******************************************************************************** 
!
! Real right scalar multiplication of the CEnergy type
!
!******************************************************************************** 
function c_energy_real_right_scalar(self, x) result(en)
  class(CEnergy), intent(in) :: self
  real(wp), intent(in)       :: x
  type(CEnergy)              :: en
  
  en = self

  en%e1 = x * self%e1
  en%eh = x * self%eh
  en%ex = x * self%ex
end function c_energy_real_right_scalar


!******************************************************************************** 
!
! Complex left scalar multiplication of the CEnergy type
!
!******************************************************************************** 
function c_energy_cmplx_left_scalar(x , self) result(en)
  complex(wp), intent(in)    :: x 
  class(CEnergy), intent(in) :: self
  type(CEnergy)              :: en
  
  en = self

  en%e1 = x * self%e1
  en%eh = x * self%eh
  en%ex = x * self%ex
end function c_energy_cmplx_left_scalar


!******************************************************************************** 
!
! Complex right scalar multiplication of the CEnergy type
!
!******************************************************************************** 
function c_energy_cmplx_right_scalar(self, x) result(en)
  class(CEnergy), intent(in) :: self
  complex(wp), intent(in)    :: x
  type(CEnergy)              :: en
  
  en = self

  en%e1 = x * self%e1
  en%eh = x * self%eh
  en%ex = x * self%ex
end function c_energy_cmplx_right_scalar


!******************************************************************************** 
!
! Division of the CEnergy type by a real scalar value
!
!******************************************************************************** 
function c_energy_real_scalar_division(self, x) result(en)
  class(CEnergy), intent(in) :: self
  real(wp), intent(in)       :: x
  type(CEnergy)              :: en
  
  en = self

  en%e1 = self%e1 / x
  en%eh = self%eh / x
  en%ex = self%ex / x
end function c_energy_real_scalar_division


!******************************************************************************** 
!
! Division of the CEnergy type by a complex scalar value
!
!******************************************************************************** 
function c_energy_cmplx_scalar_division(self, x) result(en)
  class(CEnergy), intent(in) :: self
  complex(wp), intent(in)    :: x
  type(CEnergy)              :: en
  
  en = self

  en%e1 = self%e1 / x
  en%eh = self%eh / x
  en%ex = self%ex / x
end function c_energy_cmplx_scalar_division


!******************************************************************************** 
!
! Extract e1 eh ex from the CEnergy type to the array
!
!******************************************************************************** 
function c_energy_to_array(self) result(energy_array)
  class(CEnergy), intent(in) :: self
  complex(wp)                :: energy_array(en_len)
  
  energy_array(1) = self%e1
  energy_array(2) = self%eh
  energy_array(3) = self%ex
end function c_energy_to_array


!******************************************************************************** 
!
! Compute electronic Energy: e0 + e1 + eh + ex
!
!******************************************************************************** 
elemental function get_electron_c_energy(self) result(el_en)
  class(CEnergy), intent(in) :: self
  complex(wp)                :: el_en

  el_en = self%e0 + self%e1 + self%eh + self%ex
end function get_electron_c_energy


!******************************************************************************** 
!
! Compute total CEnergy: enuc + e0 + e1 + eh + ex
!
!******************************************************************************** 
elemental function get_total_c_energy(self) result(tot_en)
  class(CEnergy), intent(in) :: self
  complex(wp)                :: tot_en

  tot_en = self%enuc + self%e0 + self%e1 + self%eh + self%ex
end function get_total_c_energy


!******************************************************************************** 
!
! MPI sum of the CEnergy type
!
!******************************************************************************** 
subroutine c_energy_mpisum(self, comm)
  class(CEnergy), intent(inout)      :: self
  type(mpi_communicator), intent(in) :: comm

  call comm%mpisum(self%e1)
  call comm%mpisum(self%eh)
  call comm%mpisum(self%ex)
end subroutine c_energy_mpisum
 

!******************************************************************************** 
!
! Report the CEnergy type
!
!******************************************************************************** 
subroutine report_c_energy(self, fh, method)
  class(CEnergy), intent(in)             :: self
  type(FileHandle), intent(in)           :: fh
  character(len=*), optional, intent(in) :: method
  !local variables
  complex(wp)                            :: el_en, tot_en
  character(len=:), allocatable          :: method_
  
  el_en = self%electron_energy()
  tot_en = self%total_energy()

  if (present(method)) then
    method_ = trim(method) // "_"
  else
    method_ = ""
  end if
  
  write(fh%funit,*)
  write(fh%funit,*)   method_ // "energy of the electronic-ionic system:"
  write(fh%funit,*)   "---------------------------------------------------------"
  write(fh%funit,100) "contribution", "     Re E      ", "    Im E    "
  write(fh%funit,100) "------------", "     ----      ", "    ----    "
  write(fh%funit,101) method_ // "enuc =", self%enuc
  write(fh%funit,101) method_ // "e0   =", self%e0
  write(fh%funit,101) method_ // "e1   =", real(self%e1,wp),  aimag(self%e1)
  write(fh%funit,101) method_ // "eh   =", real(self%eh,wp),  aimag(self%eh)
  write(fh%funit,101) method_ // "ex   =", real(self%ex,wp),  aimag(self%ex)
  write(fh%funit,101) method_ // "e    =", real(el_en,wp),    aimag(el_en)
  write(fh%funit,101) method_ // "etot =", real(tot_en,wp),   aimag(tot_en)

  100 format (1x,a,t30,2(a,2x))
  101 format (1x,a,t30,f14.8,2x,f12.8)
end subroutine report_c_energy


!******************************************************************************** 
!
! Report the CEnergy type with standard deviation
!
!******************************************************************************** 
subroutine report_c_std_energy(self, dself, fh, method)
  class(CEnergy), intent(in)             :: self
  type(EnergyStd), intent(in)            :: dself
  type(FileHandle), intent(in)           :: fh
  character(len=*), optional, intent(in) :: method
  !local variables
  complex(wp)                            :: el_en, tot_en
  character(len=:), allocatable          :: method_
  
  el_en = self%electron_energy()
  tot_en = self%total_energy()

  if (present(method)) then
    method_ = trim(method) // "_"
  else
    method_ = ""
  end if

  write(fh%funit,*)
  write(fh%funit,*)   method_ // "energy of the electronic-ionic system:"
  write(fh%funit,*)   "---------------------------------------------------------"
  write(fh%funit,100) "contribution", "     Re E      ", "    Im E    ", "    stddev    ", "corr. length", " stddev final "
  write(fh%funit,100) "------------", "     ----      ", "    ----    ", "    ------    ", "------------", " ------------ "
  write(fh%funit,101) method_ // "enuc =", self%enuc
  write(fh%funit,101) method_ // "e0   =", self%e0
  write(fh%funit,101) method_ // "e1   =", real(self%e1, wp), aimag(self%e1), dself%e1/sqrt(dself%e1_corr_len), dself%e1_corr_len, dself%e1
  write(fh%funit,101) method_ // "eh   =", real(self%eh, wp), aimag(self%eh), dself%eh/sqrt(dself%eh_corr_len), dself%eh_corr_len, dself%eh
  write(fh%funit,101) method_ // "ex   =", real(self%ex, wp), aimag(self%ex), dself%ex/sqrt(dself%ex_corr_len), dself%ex_corr_len, dself%ex
  write(fh%funit,101) method_ // "e    =", real(el_en, wp),   aimag(el_en),   dself%e/sqrt(dself%e_corr_len), dself%e_corr_len,  dself%e
  write(fh%funit,101) method_ // "etot =", real(tot_en, wp),  aimag(tot_en),  dself%e/sqrt(dself%e_corr_len), dself%e_corr_len,  dself%e
  
  100 format (1x,a,t30,5(a,2x))
  101 format (1x,a,t30,f14.8,2x,f12.8,2x,es14.6,2x,2x,f8.2,2x,2x,es14.6)
end subroutine report_c_std_energy


!******************************************************************************** 
!
! EnergyStd constructor
!
!******************************************************************************** 
function init_energy_std(e1, eh, ex, e) result(self)
  real(wp), intent(in) :: e1
  real(wp), intent(in) :: eh
  real(wp), intent(in) :: ex
  real(wp), intent(in) :: e
  type(EnergyStd)      :: self

  self%e1 = e1
  self%eh = eh
  self%ex = ex
  self%e = e

  self%e1_corr_len = 1.0_wp
  self%eh_corr_len = 1.0_wp
  self%ex_corr_len = 1.0_wp
  self%e_corr_len = 1.0_wp
end function init_energy_std


!******************************************************************************** 
!
! EnergyStd constructor (single value option)
!
!******************************************************************************** 
function init_energy_std_val(e) result(self)
  real(wp), intent(in) :: e
  type(EnergyStd)      :: self

  self%e1 = e
  self%eh = e
  self%ex = e
  self%e = e

  self%e1_corr_len = 1.0_wp
  self%eh_corr_len = 1.0_wp
  self%ex_corr_len = 1.0_wp
  self%e_corr_len = 1.0_wp
end function init_energy_std_val


!******************************************************************************** 
!
! Read Energy array from the file (file unit given)
!
!******************************************************************************** 
subroutine read_c_energies(en, funit)
  type(CEnergy), intent(inout) :: en(:)
  integer, intent(in)          :: funit
  !local variables            
  integer                      :: i
  
  do i = 1, size(en)
    read(funit) en(i)%enuc
  end do

  do i = 1, size(en)
    read(funit) en(i)%e0
  end do

  do i = 1, size(en)
    read(funit) en(i)%e1
  end do

  do i = 1, size(en)
    read(funit) en(i)%eh
  end do

  do i = 1, size(en)
    read(funit) en(i)%ex
  end do
end subroutine read_c_energies


!******************************************************************************** 
!
! Write CEnergy array to the file
!
!******************************************************************************** 
subroutine write_c_energies(en, funit)
  type(CEnergy), intent(in) :: en(:)
  integer, intent(in)       :: funit
  !local variables            
  integer                   :: i
  
  do i = 1, size(en)
    write(funit) en(i)%enuc
  end do

  do i = 1, size(en)
    write(funit) en(i)%e0
  end do

  do i = 1, size(en)
    write(funit) en(i)%e1
  end do

  do i = 1, size(en)
    write(funit) en(i)%eh
  end do

  do i = 1, size(en)
    write(funit) en(i)%ex
  end do
end subroutine write_c_energies

end module energy_types