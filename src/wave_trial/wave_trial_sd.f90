! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module wave_trial_sd

#include "../preproc.inc"
use constants
use file_handle, only: FileHandle
use mpi, only: mpi_communicator, comm_world
use profiling
use qmcfort_io
use energy_types
use hamilton_vars, only: hamil_descriptor
use hamilton_layout
use hamilton_operations, only: sum_gres_potential
use hamilton_type
use nosd
use wave_trial_des, only: WaveTrialDes
use wave_trial, only: WaveTrial

implicit none

private
public :: WaveTrialSD, init_wave_trial_sd

type, extends(WaveTrial) :: WaveTrialSD
  real(wp), allocatable :: coeff(:,:,:)
  real(wp), allocatable :: coeff_occ(:,:)
  real(wp), allocatable :: h1phi(:,:)
  real(wp), allocatable :: h1phi_s(:,:)
  real(wp), allocatable :: h2phi_h(:,:)
  real(wp), allocatable :: h2phi_s(:,:)
  real(wp), allocatable :: h2phi_x_chol(:,:)
  real(wp), allocatable :: h2phi_x_int(:,:,:)
  real(wp), allocatable :: h2phi_full(:,:,:,:)
  real(wp), allocatable :: hselfphi(:,:,:)
  real(wp), allocatable :: hselfphi_gres(:,:,:,:)
contains
  procedure :: setup_coeff => setup_coeff_sd
  procedure :: setup_h1 => setup_h1_sd
  procedure :: setup_h2 => setup_h2_sd
  procedure :: setup_hartree => setup_hartree_sd
  procedure :: setup_exchange => setup_exchange_sd
  procedure :: setup_self => setup_self_sd
  procedure :: exchange_calculator => exchange_calculator_sd
  procedure :: exchange_calculator_gres => exchange_calculator_gres_sd

  !deferred routines
  procedure :: assign_wave_trial => assign_wave_trial_sd
  procedure :: setup => setup_wave_trial_sd
  procedure :: report => report_wave_trial_sd
  procedure :: sizeof => sizeof_wave_trial_sd

  procedure :: trial_energy_full => trial_energy_sd
  procedure :: trial_h2_gres_mean => trial_h2_gres_mean_sd
  procedure :: initial_energy_full => trial_energy_sd
  procedure :: initial_h2_gres_mean => trial_h2_gres_mean_sd

  procedure :: get_ovlp_batch => get_ovlp_batch_sd
  procedure :: h1_mean_batch => h1_mean_batch_sd
  procedure :: h2_gres_mean_batch => h2_gres_mean_batch_sd
  procedure :: exchange_energy_batch => exchange_energy_batch_sd
  procedure :: exchange_energy_batch_gres => exchange_energy_batch_gres_sd
  procedure :: self_energy_batch => self_energy_batch_sd
  procedure :: self_energy_batch_gres => self_energy_batch_gres_sd
  procedure :: local_energy_batch => local_energy_batch_sd
  procedure :: local_energy_batch_gres => local_energy_batch_gres_sd
end type WaveTrialSD

interface WaveTrialSD
  module procedure finit_wave_trial_sd
end interface WaveTrialSD

contains

!******************************************************************************** 
!
! Initialization of the WaveTrialSD object
!
!******************************************************************************** 
subroutine init_wave_trial_sd(wtdes, ham, coeff, self)
  type(WaveTrialDes), intent(in)     :: wtdes
  type(Hamilton), target, intent(in) :: ham
  real(wp), intent(in)               :: coeff(:,:,:)
  type(WaveTrialSD), intent(out)     :: self

  self%des = wtdes
  allocate(self%coeff, source=coeff)

  call self%base_init(ham)

  call self%setup_coeff(self%ham%des)

!debug: wave_trial : this feauter temporarily removed
  !call self%estimate_delta_x()
end subroutine init_wave_trial_sd


!******************************************************************************** 
!
! WaveTrialSD constructor
!
!******************************************************************************** 
function finit_wave_trial_sd(wtdes, ham, coeff) result(self)
  type(WaveTrialDes), intent(in)     :: wtdes
  type(Hamilton), target, intent(in) :: ham
  real(wp), intent(in)               :: coeff(:,:,:)
  type(WaveTrialSD)                  :: self

  call init_wave_trial_sd(wtdes, ham, coeff, self)
end function finit_wave_trial_sd


!******************************************************************************** 
!
! Assignment operator for the WaveTrialSD
!
!******************************************************************************** 
subroutine assign_wave_trial_sd(to, from)
  class(WaveTrialSD), intent(inout) :: to
  class(WaveTrial), intent(in)      :: from

  select type (from)
    type is (WaveTrialSD)
      call to%base_assign(from)

      if (allocated(from%coeff)) then
        if (allocated(to%coeff)) deallocate(to%coeff)
        allocate(to%coeff, source=from%coeff)
      end if

      if (allocated(from%coeff_occ)) then
        if (allocated(to%coeff_occ)) deallocate(to%coeff_occ)
        allocate(to%coeff_occ, source=from%coeff_occ)
      end if

      if (allocated(from%h1phi)) then
        if (allocated(to%h1phi)) deallocate(to%h1phi)
        allocate(to%h1phi, source=from%h1phi)
      end if

      if (allocated(from%h1phi_s)) then
        if (allocated(to%h1phi_s)) deallocate(to%h1phi_s)
        allocate(to%h1phi_s, source=from%h1phi_s)
      end if

      if (allocated(from%h2phi_h)) then
        if (allocated(to%h2phi_h)) deallocate(to%h2phi_h)
        allocate(to%h2phi_h, source=from%h2phi_h)
      end if

      if (allocated(from%h2phi_s)) then
        if (allocated(to%h2phi_s)) deallocate(to%h2phi_s)
        allocate(to%h2phi_s, source=from%h2phi_s)
      end if

      if (allocated(from%h2phi_x_chol)) then
        if (allocated(to%h2phi_x_chol)) deallocate(to%h2phi_x_chol)
        allocate(to%h2phi_x_chol, source=from%h2phi_x_chol)
      end if

      if (allocated(from%h2phi_x_int)) then
        if (allocated(to%h2phi_x_int)) deallocate(to%h2phi_x_int)
        allocate(to%h2phi_x_int, source=from%h2phi_x_int)
      end if

      if (allocated(from%h2phi_full)) then
        if (allocated(to%h2phi_full)) deallocate(to%h2phi_full)
        allocate(to%h2phi_full, source=from%h2phi_full)
      end if

      if (allocated(from%hselfphi)) then
        if (allocated(to%hselfphi)) deallocate(to%hselfphi)
        allocate(to%hselfphi, source=from%hselfphi)
      end if

      if (allocated(from%hselfphi_gres)) then
        if (allocated(to%hselfphi_gres)) deallocate(to%hselfphi_gres)
        allocate(to%hselfphi_gres, source=from%hselfphi_gres)
      end if
    class default
      if (comm_world%mpirank == 0) write(*,*) "both dynamic types have to be WaveTrialSD in the assignment operator"
      stop
  end select
end subroutine assign_wave_trial_sd


!******************************************************************************** 
!
! Setup WaveTrialSD object for the energy evaluation
!
!******************************************************************************** 
subroutine setup_wave_trial_sd(self, lself)
  class(WaveTrialSD), intent(inout) :: self
  logical, optional, intent(in)     :: lself
  !local
  logical                           :: lself_

  lself_ = .true.
  if (present(lself)) lself_ = lself

  call self%setup_h1(self%ham%h1_orig)
  call self%setup_h2(self%ham%h2_gres, lself=lself_)
  call self%setup_hartree()
  call self%setup_exchange()
  if (lself_) call self%setup_self()
end subroutine setup_wave_trial_sd


!******************************************************************************** 
!
! Setup h1 related variables: h1phi
!
!******************************************************************************** 
subroutine setup_h1_sd(self, h1)
  class(WaveTrialSD), intent(inout) :: self
  real(wp), intent(in)              :: h1(:,:,:)
  !local
  real(wp), allocatable             :: h1phi(:,:,:)

  call get_h1phi(self%coeff, h1, h1phi)
  call get_h1phi_occ(self%ham%des, h1phi, self%h1phi)
end subroutine setup_h1_sd


!******************************************************************************** 
!
! Setup h2 related variables: h2phi_full, hselfphi, hselfphi_gres
!
!******************************************************************************** 
subroutine setup_h2_sd(self, h2_gres, lself)
  class(WaveTrialSD), intent(inout) :: self
  real(wp), intent(in)              :: h2_gres(:,:,:)
  logical, intent(in)               :: lself
  
  call get_h2phi(self%coeff, h2_gres, self%h2phi_full)

  if (lself) then
    call get_hselfphi(self%coeff, h2_gres, self%hselfphi_gres)
    call sum_gres_potential(self%hselfphi_gres, self%hselfphi)
  end if
end subroutine setup_h2_sd


!******************************************************************************** 
!
! Setup hartree related variables: h2phi_h
!
!******************************************************************************** 
subroutine setup_hartree_sd(self)
  class(WaveTrialSD), intent(inout) :: self

  call get_h2phi_h_occ(self%ham%des, self%h2phi_full, self%h2phi_h)
end subroutine setup_hartree_sd


!******************************************************************************** 
!
! Setup exchage energy related variables:
!
!    h2phi_x_chol
!    h2phi_x_int (only if exchange_mode = "eri")
!
!******************************************************************************** 
subroutine setup_exchange_sd(self)
  class(WaveTrialSD), intent(inout) :: self

  select case (trim(self%des%exchange_mode))
    case ("cholesky")
      call get_h2phi_x_occ_chol(self%ham%des, self%h2phi_full, self%h2phi_x_chol, self%des%compress_x, self%des%compress_x_tol)
    case ("eri")
      call get_h2phi_x_occ_chol(self%ham%des, self%h2phi_full, self%h2phi_x_chol, self%des%compress_x, self%des%compress_x_tol)
      call get_h2phi_x_occ_int(self%ham%des, self%h2phi_full, self%h2phi_x_int)
  end select
end subroutine setup_exchange_sd


!******************************************************************************** 
!
! Setup self energy related variables: h2phi_s
!
!******************************************************************************** 
subroutine setup_self_sd(self)
  class(WaveTrialSD), intent(inout) :: self

  call get_h2phi_h_occ(self%ham%des, self%hselfphi_gres, self%h2phi_s)
  call get_h1phi_occ(self%ham%des, self%hselfphi, self%h1phi_s)
end subroutine setup_self_sd


!******************************************************************************** 
!
! Report WaveTrialSD object
!
!******************************************************************************** 
subroutine report_wave_trial_sd(self, comm, fh)
  class(WaveTrialSD), intent(inout)  :: self
  type(mpi_communicator), intent(in) :: comm
  type(FileHandle), intent(in)       :: fh
  !local
  real(dp)                           :: mem
  type(CEnergy)                      :: e_init, e_trial

  call self%initial_energy(e_init)
  call self%trial_energy(e_trial)

  mem = self%sizeof()

  if (comm%mpirank == 0) then
    write(fh%funit,*)
    write(fh%funit,*)   "Trial wave function descriptor" 
    write(fh%funit,*)   "-------------------------------" 
    write(fh%funit,100) "trial wave function type        ", "WaveTrialSD"
    write(fh%funit,100) "size of the object              ", write_bytes(mem)
    write(fh%funit,100) "Exchange energy mode            ", self%des%exchange_mode 
    write(fh%funit,101) "block size for Hartree terms    ", self%des%h_block_size
    write(fh%funit,101) "block size for exchange terms   ", self%des%x_block_size
  
!debug: wave_trial : return when estimate_delta_x reimplemented
    !if (self%exchange_mode == "cholesky") then
    !  write(fh%funit,102)  "Compressed Cholesky vectors X   ", self%compress_x
    !  write(fh%funit,103)  "Tolerance for the compression   ", self%compress_x_tol
    !  write(fh%funit,103)  "Exchange energy correction      ", self%delta_x
    !end if

    call e_trial%report(fh, method="trial")
    call e_init%report(fh, method="init")
  end if

  100 format(1x, t5, a, t50, "= ", a)
  101 format(1x, t5, a, t50, "= ", 10i8)
  102 format(1x, t5, a, t50, "= ", 10l6)
  103 format(1x, t5, a, t50, "= ", 10es12.4)
end subroutine report_wave_trial_sd


!******************************************************************************** 
!
! Estimate memory consumption by WaveTrialSD object
!
!******************************************************************************** 
function sizeof_wave_trial_sd(self) result(mem)
  class(WaveTrialSD), intent(in) :: self
  real(dp)                       :: mem

  mem = sizeof(self) + sizeof(self%coeff) + sizeof(self%coeff_occ) + sizeof(self%h1phi) + sizeof(self%h1phi_s) + &
        sizeof(self%h2phi_h) + sizeof(self%h2phi_x_chol) + sizeof(self%h2phi_x_int) + sizeof(self%h2phi_s) +  &
        sizeof(self%h2phi_full) + sizeof(self%hselfphi) + sizeof(self%hselfphi_gres)
end function sizeof_wave_trial_sd


!******************************************************************************** 
!
! Obtain occupied orbitals self%coeff_occ from all orbitals stored in self%coeff
!
!******************************************************************************** 
subroutine setup_coeff_sd(self, hdes)
  class(WaveTrialSD), intent(inout)  :: self
  type(hamil_descriptor), intent(in) :: hdes
  !local
  integer                            :: spin, sp, nel, w1, w2

  if (allocated(self%coeff_occ)) deallocate(self%coeff_occ)
  allocate(self%coeff_occ(hdes%n, hdes%nocc))

  do spin = 1, hdes%ispin
    sp = hdes%get_spin(spin) 
    nel = hdes%nel(spin)
    w1 = hdes%nell(spin) + 1
    w2 = hdes%nell(spin) + nel
    self%coeff_occ(:,w1:w2) = self%coeff(:,1:nel,sp)
  end do
end subroutine setup_coeff_sd


!******************************************************************************** 
!
! Estimate trial energy using local_energy with trial orbitals
!
!******************************************************************************** 
subroutine trial_energy_sd(self, e1, lgmean, eh, ex, eself)
  class(WaveTrialSD), intent(inout)     :: self
  complex(wp), intent(out)              :: e1
  complex(wp), allocatable, intent(out) :: lgmean(:)
  complex(wp), allocatable, intent(out) :: eh(:)
  complex(wp), allocatable, intent(out) :: ex(:)
  complex(wp), allocatable, intent(out) :: eself(:)

  call self%local_energy_gres(self%coeff_occ*onec, e1, lgmean, eh, ex, eself)
  ex(1) = ex(1) + self%des%delta_x
end subroutine trial_energy_sd


!******************************************************************************** 
!
! Calculate trial exectation vlaue of the Cholesky vectors
!
!    <L_g> = <Psi_T|L_g|Psi_T> / <Psi_T|Psi_T>
!
!******************************************************************************** 
subroutine trial_h2_gres_mean_sd(self, lgmean)
  class(WaveTrialSD), intent(inout)     :: self
  complex(wp), allocatable, intent(out) :: lgmean(:)
  !local
  complex(wp), allocatable              :: eh(:)
  
  call self%h2_gres_mean(self%coeff_occ*onec, lgmean, eh)
end subroutine trial_h2_gres_mean_sd


!******************************************************************************** 
!
! Calculate overlap between trial determinant and walker determinants
!
! Actual low-level implementation in nosd.f: get_ovlp_batch
!
!******************************************************************************** 
subroutine get_ovlp_batch_sd(self, coeff, ovlp)
  class(WaveTrialSD), intent(inout)     :: self
  complex(wp), intent(in)               :: coeff(:,:)
  complex(wp), allocatable, intent(out) :: ovlp(:)
  !local
  character(len=*), parameter           :: proc_name = "get_overlap_sd"
  
  if (profile_code) call start_profiling(proc_name)

  call get_ovlp_batch(self%ham%des, self%coeff_occ, coeff, ovlp)

  if (profile_code) call end_profiling(proc_name)
end subroutine get_ovlp_batch_sd


!******************************************************************************** 
!
! Calculate expectation value of the one-electron Hamiltonian for the batch of walkers
!
!    E1 = <Psi_T | H1 | Phi>
!
! Actual low-level implementation in nosd.f: h1_mean_no_batch
!
!******************************************************************************** 
subroutine h1_mean_batch_sd(self, coeff, e1)
  class(WaveTrialSD), intent(inout)     :: self
  complex(wp), intent(in)               :: coeff(:,:)
  complex(wp), allocatable, intent(out) :: e1(:)
  !local
  complex(wp), allocatable              :: coeff_bi(:,:)
  character(len=*), parameter           :: proc_name = "h1_mean_sd"

  if (profile_code) call start_profiling(proc_name)

  call biorth_orb_batch(self%ham%des, self%coeff_occ, coeff, coeff_bi)
  call h1_mean_no_batch(self%ham%des, self%h1phi, coeff_bi, e1)

  if (profile_code) call end_profiling(proc_name)
end subroutine h1_mean_batch_sd


!******************************************************************************** 
!
! Calculate expectation value of the Cholesky vectors for the batch of walkers
!
!    l_{g,w} = <Psi_T | L_g | Phi_w>
!
!    coeff --> coeff_bi = coeff * (coeff_T * coeff)^{-1}
!    l_{g,w} = h2phi_h{g,I} coeff_bi_{I,w} 
!   
! Actual low-level implementation in nosd.f: h2_gres_mean_no_batch
!
!******************************************************************************** 
subroutine h2_gres_mean_batch_sd(self, coeff, lgmean, eh)
  class(WaveTrialSD), intent(inout)     :: self
  complex(wp), intent(in)               :: coeff(:,:)
  complex(wp), allocatable, intent(out) :: lgmean(:,:)
  complex(wp), allocatable, intent(out) :: eh(:,:)
  !local
  complex(wp), allocatable              :: coeff_bi(:,:)
  character(len=*), parameter           :: proc_name = "h2_gres_mean_sd"
  
  if (profile_code) call start_profiling(proc_name)

  call biorth_orb_batch(self%ham%des, self%coeff_occ, coeff, coeff_bi)
  call h2_gres_mean_no_batch(self%ham%des, self%h2phi_h, coeff_bi, lgmean, self%des%h_block_size)

  allocate(eh, mold=lgmean)
  eh = 0.5_wp * lgmean**2

  if (profile_code) call end_profiling(proc_name)
end subroutine h2_gres_mean_batch_sd


!******************************************************************************** 
!
! Calculate exchange energy for the batch of walkers
!
!    Ex_w = <Psi_T | Hx | Phi_w>
! 
! exchange_calculator chooses the method to calculate exchange energy
!
! Actual low-level implementations in nosd.f
!
!******************************************************************************** 
subroutine exchange_energy_batch_sd(self, coeff, ex)
  class(WaveTrialSD), intent(inout)     :: self
  complex(wp), intent(in)               :: coeff(:,:)
  complex(wp), allocatable, intent(out) :: ex(:)
  !local
  complex(wp), allocatable              :: coeff_bi(:,:)
  character(len=*), parameter           :: proc_name = "exchange_sd"
  
  if (profile_code) call start_profiling(proc_name)

  call biorth_orb_batch(self%ham%des, self%coeff_occ, coeff, coeff_bi)
  call self%exchange_calculator(self%ham%des, coeff_bi, ex)
  ex = ex + self%des%delta_x

  if (profile_code) call end_profiling(proc_name)
end subroutine exchange_energy_batch_sd


!******************************************************************************** 
!
! Calculate g-resolved exchange energy for the batch of walkers
!
!    Ex_{g,w} = <Psi_T | Hx_g | Phi_w>
! 
! Actual low-level implementations in nosd.f : exchange_energy_no_chol_batch
!
!******************************************************************************** 
subroutine exchange_energy_batch_gres_sd(self, coeff, ex)
  class(WaveTrialSD), intent(inout)     :: self
  complex(wp), intent(in)               :: coeff(:,:)
  complex(wp), allocatable, intent(out) :: ex(:,:)
  !local
  complex(wp), allocatable              :: coeff_bi(:,:)
  character(len=*), parameter           :: proc_name = "exchange_sd"
  
  if (profile_code) call start_profiling(proc_name)

  call biorth_orb_batch(self%ham%des, self%coeff_occ, coeff, coeff_bi)
  call self%exchange_calculator_gres(self%ham%des, coeff_bi, ex)

  if (profile_code) call end_profiling(proc_name)
end subroutine exchange_energy_batch_gres_sd


!******************************************************************************** 
!
! Setup the method for the exchange energy evaluation
!
! Note: input orbitals are biorthogonal orbitals
!
!******************************************************************************** 
subroutine exchange_calculator_sd(self, hdes, coeff_bi, ex)
  class(WaveTrialSD), intent(inout)     :: self
  type(hamil_descriptor), intent(in)    :: hdes
  complex(wp), intent(in)               :: coeff_bi(:,:)
  complex(wp), allocatable, intent(out) :: ex(:)
  !local
  integer                               :: w
  complex(wp), allocatable              :: ex_(:,:)

  select case (self%des%exchange_mode)
    case ("cholesky")
      call exchange_energy_no_chol_batch(hdes, self%h2phi_x_chol, coeff_bi, ex_, self%des%x_block_size)
      allocate(ex(size(ex_, 2)))
      do w = 1, size(ex_, 2)
        ex(w) = sum(ex_(:,w))
      end do
    case ("eri")
      call exchange_energy_no_int_batch(hdes, self%h2phi_x_int , coeff_bi, ex, self%des%x_block_size)
  end select
end subroutine exchange_calculator_sd


!******************************************************************************** 
!
! Setup the method for the g-resolved exchange energy evaluation
!
! Note: input orbitals are biorthogonal orbitals
!
!******************************************************************************** 
subroutine exchange_calculator_gres_sd(self, hdes, coeff_bi, ex)
  class(WaveTrialSD), intent(inout)     :: self
  type(hamil_descriptor), intent(in)    :: hdes
  complex(wp), intent(in)               :: coeff_bi(:,:)
  complex(wp), allocatable, intent(out) :: ex(:,:)

  call exchange_energy_no_chol_batch(hdes, self%h2phi_x_chol, coeff_bi, ex, self%des%x_block_size)
end subroutine exchange_calculator_gres_sd


!******************************************************************************** 
!
! Calculate self energy for the batch of walkers
!
!    Eself_w = <Psi_T | Hself | Phi_w>
! 
! Actual low-level implementations in nosd.f: h2_gres_mean_no_batch
!
!******************************************************************************** 
subroutine self_energy_batch_sd(self, coeff, eself)
  class(WaveTrialSD), intent(inout)     :: self
  complex(wp), intent(in)               :: coeff(:,:)
  complex(wp), allocatable, intent(out) :: eself(:)
  !local
  complex(wp), allocatable              :: coeff_bi(:,:)
  character(len=*), parameter           :: proc_name = "self_energy_sd"
  
  if (profile_code) call start_profiling(proc_name)

  call biorth_orb_batch(self%ham%des, self%coeff_occ, coeff, coeff_bi)
  call h1_mean_no_batch(self%ham%des, self%h1phi_s, coeff_bi, eself)

  if (profile_code) call end_profiling(proc_name)
end subroutine self_energy_batch_sd


!******************************************************************************** 
!
! Calculate g-resolved self energy for the batch of walkers
!
!    Eself_{g,w} = <Psi_T | Hself_g | Phi_w>
! 
! Actual low-level implementations in nosd.f: h2_gres_mean_no_batch
!
!******************************************************************************** 
subroutine self_energy_batch_gres_sd(self, coeff, eself)
  class(WaveTrialSD), intent(inout)     :: self
  complex(wp), intent(in)               :: coeff(:,:)
  complex(wp), allocatable, intent(out) :: eself(:,:)
  !local
  complex(wp), allocatable              :: coeff_bi(:,:)
  character(len=*), parameter           :: proc_name = "self_energy_gres_sd"
  
  if (profile_code) call start_profiling(proc_name)

  call biorth_orb_batch(self%ham%des, self%coeff_occ, coeff, coeff_bi)
  call h2_gres_mean_no_batch(self%ham%des, self%h2phi_s, coeff_bi, eself, self%des%h_block_size)

  if (profile_code) call end_profiling(proc_name)
end subroutine self_energy_batch_gres_sd


!******************************************************************************** 
!
! Calculate local energy contributions for a batch of walkers
!
!    E1 = <Psi_T | H1 | Phi>
!    l_g = <Psi_T | L_g | Phi>
!    Ex = <Psi_T | Hx | Phi>
!    Eself = <Psi_T | Hself | Phi>
!
! Calculates one-body, Hartree, and exchange energy accodrding to routines
! in nosd.f
!
!******************************************************************************** 
subroutine local_energy_batch_sd(self, coeff, e1, lgmean, eh, ex, eself, lsamplex, lself, en)
  class(WaveTrialSD), intent(inout)                 :: self
  complex(wp), intent(in)                           :: coeff(:,:)
  complex(wp), allocatable, intent(out)             :: e1(:)
  complex(wp), allocatable, intent(out)             :: lgmean(:,:)
  complex(wp), allocatable, intent(out)             :: eh(:)
  complex(wp), allocatable, intent(out)             :: ex(:)
  complex(wp), allocatable, intent(out)             :: eself(:)
  logical, optional, intent(in)                     :: lsamplex
  logical, optional, intent(in)                     :: lself
  type(CEnergy), allocatable, optional, intent(out) :: en(:)
  !local
  integer                                           :: w, nw
  logical                                           :: lsamplex_, lself_
  complex(wp), allocatable                          :: coeff_bi(:,:)
  complex(wp), allocatable                          :: eself_(:,:)
  character(len=*), parameter                       :: proc_name = "local_energy_sd"
  
  if (profile_code) call start_profiling(proc_name)

  lsamplex_ = .true.
  if (present(lsamplex)) lsamplex_ = lsamplex

  lself_ = .true.
  if (present(lself)) lself_ = lself

  nw = size(coeff, 2) / self%ham%des%nocc

  call biorth_orb_batch(self%ham%des, self%coeff_occ, coeff, coeff_bi)
  call h1_mean_no_batch(self%ham%des, self%h1phi, coeff_bi, e1)
  call h2_gres_mean_no_batch(self%ham%des, self%h2phi_h, coeff_bi, lgmean)

  if (lself_) then
    call h2_gres_mean_no_batch(self%ham%des, self%h2phi_s, coeff_bi, eself_)
  else
    allocate(eself_, mold=lgmean)
    eself_ = zeroc
  end if

  if (lsamplex_) then
    call self%exchange_calculator(self%ham%des, coeff_bi, ex)
    ex = ex + self%des%delta_x
  else
    allocate(ex(nw))
    ex = zeroc
  end if

  allocate(eh(nw))
  allocate(eself(nw))
  do w = 1, nw
    eh(w) = 0.5_wp * sum(lgmean(:,w)**2)
    eself(w) = sum(eself_(:,w))
  end do

  if (present(en)) then
    allocate(en(nw))
    do w = 1, nw
      en(w) = self%pack_energy(e1(w), eh(w), ex(w))
    end do
  end if
  
  if (profile_code) call end_profiling(proc_name)
end subroutine local_energy_batch_sd


!******************************************************************************** 
!
! Calculate g-resolved local energy contributions for a single walker
!
!    E1_w = <Psi_T | H1 | Phi_w>
!    l_{g,w} = <Psi_T | L_g | Phi_w>
!    Ex_{g,w} = <Psi_T | Hx_g | Phi_w>
!    Eself_{g,w} = <Psi_T | L_g^2 | Phi_w>
!
! Calculates one-body, Hartree, and exchange energy accodrding to routines
! in nosd.f
!
!******************************************************************************** 
subroutine local_energy_batch_gres_sd(self, coeff, e1, lgmean, eh, ex, eself, lsamplex, lself, en)
  class(WaveTrialSD), intent(inout)                 :: self
  complex(wp), intent(in)                           :: coeff(:,:)
  complex(wp), allocatable, intent(out)             :: e1(:)
  complex(wp), allocatable, intent(out)             :: lgmean(:,:)
  complex(wp), allocatable, intent(out)             :: eh(:,:)
  complex(wp), allocatable, intent(out)             :: ex(:,:)
  complex(wp), allocatable, intent(out)             :: eself(:,:)
  logical, optional, intent(in)                     :: lsamplex
  logical, optional, intent(in)                     :: lself
  type(CEnergy), allocatable, optional, intent(out) :: en(:)
  !local
  integer                                           :: w, nw
  logical                                           :: lsamplex_, lself_
  complex(wp), allocatable                          :: coeff_bi(:,:)
  character(len=*), parameter                       :: proc_name = "local_energy_gres_sd"
  
  if (profile_code) call start_profiling(proc_name)

  lsamplex_ = .true.
  if (present(lsamplex)) lsamplex_ = lsamplex

  lself_ = .true.
  if (present(lself)) lself_ = lself

  nw = size(coeff, 2) / self%ham%des%nocc

  call biorth_orb_batch(self%ham%des, self%coeff_occ, coeff, coeff_bi)
  call h1_mean_no_batch(self%ham%des, self%h1phi, coeff_bi, e1)
  call h2_gres_mean_no_batch(self%ham%des, self%h2phi_h, coeff_bi, lgmean)

  if (lself_) then
    call h2_gres_mean_no_batch(self%ham%des, self%h2phi_s, coeff_bi, eself)
  else
    allocate(eself, mold=lgmean)
    eself = zeroc
  end if

  if (lsamplex_) then
    call self%exchange_calculator_gres(self%ham%des, coeff_bi, ex)
  else
    allocate(ex, mold=lgmean)
    ex = zeroc
  end if

  allocate(eh, mold=lgmean)
  eh = 0.5_wp * lgmean**2

  if (present(en)) then
    allocate(en(nw))
    do w = 1, nw
      en(w) = self%pack_energy_gres(e1(w), eh(:,w), ex(:,w))
    end do
  end if
  
  if (profile_code) call end_profiling(proc_name)
end subroutine local_energy_batch_gres_sd

end module wave_trial_sd