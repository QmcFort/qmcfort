! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module wave_trial_factories

#include "../preproc.inc"
use constants
use mpi, only: comm_world
use file_handle, only: FileHandle
use array_numpy_file_io, only: ArrayNumpyFileIO
use qmcfort_io, only: io
use hamilton_type, only: Hamilton, adjust_wavespin
use slater, only: slater_det
use wave_trial_des, only: WaveTrialDes
use wave_trial, only: WaveTrial
use wave_trial_sd, only: WaveTrialSD, init_wave_trial_sd
use wave_trial_cas, only: WaveTrialCas, init_wave_trial_cas
use wave_trial_noci, only: WaveTrialNOCI, init_wave_trial_noci

implicit none

private 
public :: wave_trial_factory_afqmc, wave_trial_factory_noci

enum, bind(c)
  enumerator :: SD = 0, CAS, NOCI
end enum

contains

!******************************************************************************** 
!
! Select type of the WaveTrial object by inspecting the existence of
! required files 
!
!******************************************************************************** 
function select_wave_trial_type(wave_trial_request) result(wave_trial_type)
  character(len=*), intent(in) :: wave_trial_request
  integer                      :: wave_trial_type
  !local variables
  logical                      :: lexist

  select case (wave_trial_request)
    case ("sd")
      wave_trial_type = SD
      return
    case ("cas")
      wave_trial_type = CAS
      return
    case ("noci")
      wave_trial_type = NOCI
      return
  end select

  inquire(file="orbitals_trial_noci", exist=lexist)
  if (lexist) then
    wave_trial_type = NOCI
    return
  end if

  inquire(file="cas_space", exist=lexist)
  if (lexist) then
    wave_trial_type = CAS
    return
  end if

  wave_trial_type = SD
end function select_wave_trial_type


!******************************************************************************** 
!
! WaveTrial Factory to use in AFQMC module
!
!******************************************************************************** 
subroutine wave_trial_factory_afqmc(wave_trial_request, wtdes, ham, phi_t)
  character(len=*), intent(in)               :: wave_trial_request
  type(WaveTrialDes), intent(inout)          :: wtdes
  type(Hamilton), intent(inout)              :: ham
  class(WaveTrial), allocatable, intent(out) :: phi_t
  !logical variables 
  integer                                    :: wave_trial_type

  wave_trial_type = select_wave_trial_type(wave_trial_request)

  select case (wave_trial_type)
    case (SD)
      call poly_wave_trial_factory_sd(wtdes, ham, phi_t)
    case (CAS)
      call poly_wave_trial_factory_cas(wtdes, ham, phi_t)
    case (NOCI)
      call poly_wave_trial_factory_noci(wtdes, ham, phi_t)
  end select
end subroutine wave_trial_factory_afqmc


!******************************************************************************** 
!
! WaveTrial Factory for the WaveTrialSD wave functions
!
!******************************************************************************** 
subroutine poly_wave_trial_factory_sd(wtdes, ham, phi_t)
  type(WaveTrialDes), intent(inout)          :: wtdes
  type(Hamilton), intent(inout)              :: ham
  class(WaveTrial), allocatable, intent(out) :: phi_t
  !local variables
  integer                                    :: nocc
  real(wp), allocatable                      :: coeff(:,:,:)
  type(WaveTrialSD), allocatable             :: phi_t_sd

  !debug: nel_space and nell_space should already be by construction equal to nel and nell
  ham%des%nel_space = ham%des%nel
  ham%des%nell_space = ham%des%nell
  nocc = maxval(ham%des%nel_space)

  call read_trial_orbitals_real(ham, coeff)

  allocate(WaveTrialSD :: phi_t_sd)
  call init_wave_trial_sd(wtdes, ham, coeff(:,1:nocc,:), phi_t_sd)
  call move_alloc(from=phi_t_sd, to=phi_t)
end subroutine poly_wave_trial_factory_sd


!******************************************************************************** 
!
! WaveTrial Factory for the WaveTrialCas wave functions
!
!******************************************************************************** 
subroutine poly_wave_trial_factory_cas(wtdes, ham, phi_t)
  type(WaveTrialDes), intent(inout)          :: wtdes
  type(Hamilton), intent(inout)              :: ham
  class(WaveTrial), allocatable, intent(out) :: phi_t
  !local variables
  integer                                    :: nocc, nfrozen(2)
  logical                                    :: lexist
  real(wp), allocatable                      :: coeff(:,:,:), ci_coeff(:)
  type(slater_det), allocatable              :: ci_det(:)
  type(WaveTrialCas), allocatable            :: phi_t_cas

  !debug: nel_space and nell_space should already be by construction equal to nel and nell
  ham%des%nel_space = ham%des%nel
  ham%des%nell_space = ham%des%nell
  nocc = maxval(ham%des%nel_space)

  call read_trial_orbitals_real(ham, coeff)

  inquire(file="cas_space", exist=lexist)
  if (lexist) then
    call read_cas_space(wtdes, ci_coeff, ci_det)

    !debug: this should also be reconsidered
    wtdes%nell_active(1) = 0
    wtdes%nell_active(2) = wtdes%nel_active(1)

    !debug: this should also be reconsidered
    nfrozen(1) = ham%des%nel(1) - wtdes%nel_active(1)
    nfrozen(2) = ham%des%nel(2) - wtdes%nel_active(2)
    if (nfrozen(1) == nfrozen(2)) then
      wtdes%nfrozen = nfrozen(1)
    else
      if (comm_world%mpirank == 0) write(io%screen%funit,*)      "number of frozen electrons inconsistent over two spin channels"
      if (comm_world%mpirank == 0) write(io%qmcfort_log%funit,*) "number of frozen electrons inconsistent over two spin channels"
      call exit
    end if

    !debug: this should also be reconsidered
    ham%des%nel_space = [wtdes%nfrozen+wtdes%nactive, wtdes%nfrozen+wtdes%nactive]
    ham%des%nell_space(1) = 0
    ham%des%nell_space(2) = ham%des%nel_space(1)

    nocc = maxval(ham%des%nel_space)

    allocate(WaveTrialCas :: phi_t_cas)
    call init_wave_trial_cas(wtdes, ham, coeff(:,1:nocc,:), ci_det, ci_coeff, phi_t_cas)
    call move_alloc(from=phi_t_cas, to=phi_t)
  else
    if (comm_world%mpirank == 0) write(io%screen%funit, *)      "WaveTrialCas can not be initialized without cas_space file - program stops now"
    if (comm_world%mpirank == 0) write(io%qmcfort_log%funit, *) "WaveTrialCas can not be initialized without cas_space file - program stops now"
    call comm_world%barrier()
    call exit()
  end if
end subroutine poly_wave_trial_factory_cas


!******************************************************************************** 
!
! WaveTrial Factory for the WaveTrialNOCI wave functions
!
!******************************************************************************** 
subroutine poly_wave_trial_factory_noci(wtdes, ham, phi_t)
  type(WaveTrialDes), intent(inout)          :: wtdes
  type(Hamilton), intent(inout)              :: ham
  class(WaveTrial), allocatable, intent(out) :: phi_t
  !local variables
  logical                                    :: lexist1,lexist2
  complex(wp), allocatable                   :: coeff(:,:,:,:), ci_coeff(:)
  type(WaveTrialNOCI), allocatable           :: phi_t_noci

  !debug: nel_space and nell_space should already be by construction equal to nel and nell
  ham%des%nel_space = ham%des%nel
  ham%des%nell_space = ham%des%nell

  inquire(file="orbitals_trial_noci", exist=lexist1)
  inquire(file="ci_coeff_noci", exist=lexist2)

  if (lexist1 .and. lexist2) then
    call read_trial_orbitals_noci(coeff, ci_coeff)
  else 
    if (comm_world%mpirank == 0) write(io%screen%funit, *)      "WaveTrialNoci can not be initialized without orbitals_trial_noci, ci_coeff_noci files - program stops now"
    if (comm_world%mpirank == 0) write(io%qmcfort_log%funit, *) "WaveTrialNoci can not be initialized without orbitals_trial_noci, ci_coeff_noci files - program stops now"
    call comm_world%barrier()
    call exit()
  end if

  allocate(WaveTrialNOCI :: phi_t_noci)
  call init_wave_trial_noci(wtdes, ham, coeff, ci_coeff, phi_t_noci)
  call move_alloc(from=phi_t_noci, to=phi_t)
end subroutine poly_wave_trial_factory_noci


!******************************************************************************** 
!
! Factory for the WaveTrialNOCI for NOCI selection
!
!******************************************************************************** 
subroutine wave_trial_factory_noci(wtdes, ham, phi_t_noci)
  type(WaveTrialDes), intent(inout)     :: wtdes
  type(Hamilton), target, intent(inout) :: ham
  type(WaveTrialNOCI)                   :: phi_t_noci
  !local
  integer                               :: nocc
  logical                               :: lexist
  real(wp), allocatable                 :: coeff_real(:,:,:)
  complex(wp), allocatable              :: coeff(:,:,:,:), ci_coeff(:)

  !debug: nel_space and nell_space should already be by construction equal to nel and nell
  ham%des%nel_space = ham%des%nel
  ham%des%nell_space = ham%des%nell
  nocc = maxval(ham%des%nel_space)

  inquire(file="orbitals_trial_noci", exist=lexist)

  if (lexist) then
    call read_trial_orbitals_noci(coeff, ci_coeff)
  else
    call read_trial_orbitals_real(ham, coeff_real)
    allocate(coeff(size(coeff_real, 1), nocc, size(coeff_real, 3), 1))
    coeff(:,:,:,1) = coeff_real(:,1:nocc,:)

    allocate(ci_coeff(1))
    ci_coeff = onec
  end if

  call init_wave_trial_noci(wtdes, ham, coeff, ci_coeff, phi_t_noci)
end subroutine wave_trial_factory_noci


!******************************************************************************** 
!
! Read real trial orbitals from the file orbitals_trial
! or use orbitals stored in Hamilton if the file is not present
!
!******************************************************************************** 
subroutine read_trial_orbitals_real(ham, coeff)
  type(Hamilton), intent(in)         :: ham
  real(wp), allocatable, intent(out) :: coeff(:,:,:)
  !local variables
  logical                            :: lexist
  type(FileHandle)                   :: fh
  type(ArrayNumpyFileIO)             :: numpy_io
  character(len=*), parameter        :: fname = "orbitals_trial"

  fh = FileHandle(fname)

  if (fh%exists()) then
    call numpy_io%load(fname, coeff)
    call adjust_wavespin(ham%des, coeff) 

    if (comm_world%mpirank == 0) write(io%screen%funit, *)      "trial orbitals read from the file"
    if (comm_world%mpirank == 0) write(io%qmcfort_log%funit, *) "trial orbitals read from the file"
  else
    allocate(coeff, source=ham%phi)
    if (comm_world%mpirank == 0) write(io%screen%funit, *)      "trial orbitals same as the initial orbitals"
    if (comm_world%mpirank == 0) write(io%qmcfort_log%funit, *) "trial orbitals same as the initial orbitals"
  end if
end subroutine read_trial_orbitals_real


!******************************************************************************** 
!
! Read real trial orbitals from the file orbitals_trial
! or use orbitals stored in Hamilton if the file is not present
!
!******************************************************************************** 
subroutine read_trial_orbitals_cmplx(ham, coeff)
  type(Hamilton), intent(in)            :: ham
  complex(wp), allocatable, intent(out) :: coeff(:,:,:)
  !local variables
  logical                               :: lexist
  type(FileHandle)                      :: fh
  type(ArrayNumpyFileIO)                :: numpy_io
  character(len=*), parameter           :: fname = "orbitals_trial"

  fh = FileHandle(fname)

  if (fh%exists()) then
    call numpy_io%load(fname, coeff)

    if (comm_world%mpirank == 0) write(io%screen%funit, *)      "trial orbitals read from the file"
    if (comm_world%mpirank == 0) write(io%qmcfort_log%funit, *) "trial orbitals read from the file"
  else
    allocate(coeff(size(ham%phi,1), size(ham%phi,2), size(ham%phi,3)))
    coeff = ham%phi
    if (comm_world%mpirank == 0) write(io%screen%funit, *)      "trial orbitals same as the initial orbitals"
    if (comm_world%mpirank == 0) write(io%qmcfort_log%funit, *) "trial orbitals same as the initial orbitals"
  end if
end subroutine read_trial_orbitals_cmplx


!******************************************************************************** 
!
! Read real trial orbitals from the file orbitals_trial_noci
!
!******************************************************************************** 
subroutine read_trial_orbitals_noci(coeff, ci_coeff)
  complex(wp), allocatable, intent(out) :: coeff(:,:,:,:)
  complex(wp), allocatable, intent(out) :: ci_coeff(:)
  !local
  type(FileHandle)                      :: fh_coeff, fh_ci_coeff
  type(ArrayNumpyFileIO)                :: numpy_io
  character(len=*), parameter           :: fname_coeff = "orbitals_trial_noci"
  character(len=*), parameter           :: fname_ci_coeff = "ci_coeff_noci"

  fh_coeff = FileHandle(fname_coeff)
  fh_ci_coeff = FileHandle(fname_ci_coeff)

  if (.not. (fh_coeff%exists() .and. fh_ci_coeff%exists())) then
    if (comm_world%mpirank == 0) write(io%screen%funit, *)      "Either orbitals_trial_noci or ci_coeff_noci are not present - program stops now"
    if (comm_world%mpirank == 0) write(io%qmcfort_log%funit, *) "Either orbitals_trial_noci or ci_coeff_noci are not present - program stops now"
    call comm_world%barrier()
    call exit()
  end if

  call numpy_io%load(fname_coeff, coeff)
  call numpy_io%load(fname_ci_coeff, ci_coeff)
end subroutine read_trial_orbitals_noci


!******************************************************************************** 
!
! Read CAS space from the file cas_space
!
! Required for the WaveTrialCas object
!
!******************************************************************************** 
subroutine read_cas_space(wtdes, ci_coeff, ci_det)
  type(WaveTrialDes), intent(inout)          :: wtdes
  real(wp), allocatable, intent(out)         :: ci_coeff(:)
  type(slater_det), allocatable, intent(out) :: ci_det(:)
  !local variables
  integer                                    :: d, ndet
  character(len=charlen)                     :: fline
  type(FileHandle)                           :: fh

  fh = FileHandle("cas_space")
  call fh%open(status="old", action="read")

  read(fh%funit,*) fline
  read(fh%funit,*) ndet, wtdes%nactive, wtdes%nel_active
  read(fh%funit,*) fline
  allocate(ci_coeff(ndet), ci_det(ndet))
  do d = 1, ndet
    read(fh%funit,"(a)") fline
    ci_coeff(d) = ci_det(d)%read_from_char(fline)
  end do

  call fh%close()

  if (comm_world%mpirank == 0) write(io%screen%funit,*)      "cas space read from the file"
  if (comm_world%mpirank == 0) write(io%qmcfort_log%funit,*) "cas space read from the file"
end subroutine read_cas_space

end module wave_trial_factories