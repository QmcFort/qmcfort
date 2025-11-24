! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module wave_trial_cas

#include "../preproc.inc"
use constants
use file_handle, only: FileHandle
use energy_types
use mpi, only: mpi_communicator, comm_world
use profiling
use qmcfort_io
use standalone, only: trace, trace2, get_nblocks
use lapack
use hamilton_vars, only: hamil_descriptor
use hamilton_layout
use hamilton_type
use cas_hamilton, only: CasHamilton
use slater
use slater_condon_rules
use nosd
use wave_trial_des, only: WaveTrialDes
use wave_trial, only: WaveTrial
use wave_trial_sd, only: WaveTrialSD

implicit none

private
public :: WaveTrialCAS, init_wave_trial_cas

type, extends(WaveTrialSD) :: WaveTrialCAS
  type(slater_det), allocatable       :: det(:)
  real(wp), allocatable               :: ci_coeff(:)
  type(CasHamilton)                   :: hcas

  real(wp), allocatable               :: h1phi_occ(:,:)
  real(wp), allocatable               :: h1phi_s_occ(:,:)
  real(wp), allocatable               :: h2phi_h_occ(:,:)
  real(wp), allocatable               :: h2phi_s_occ(:,:)
  real(wp), allocatable               :: h2phi_x_chol_occ(:,:)
  real(wp), allocatable               :: h2phi_x_int_occ(:,:,:)
  type(hamil_descriptor), allocatable :: hdes_cp
contains
  procedure          :: get_ispin => get_det_ispin
  procedure          :: setup_coeff => setup_coeff_cas
  procedure          :: setup_coeff_occ => setup_coeff_occ_cas
  procedure          :: setup_h1phi_occ => setup_h1phi_occ_cas
  procedure          :: setup_h1phi_s_occ => setup_h1phi_s_occ_cas
  procedure          :: setup_h2phi_h_occ => setup_h2phi_h_occ_cas
  procedure          :: setup_h2phi_s_occ => setup_h2phi_s_occ_cas
  procedure          :: setup_h2phi_x_occ => setup_h2phi_x_occ_cas
  procedure, private :: setup_h2phi_x_occ_chol_cas, setup_h2phi_x_occ_int_cas

  !overridden routines
  procedure          :: exchange_calculator => exchange_calculator_cas
  procedure          :: exchange_calculator_gres => exchange_calculator_gres_cas

  !deferred routines
  procedure          :: assign_wave_trial => assign_wave_trial_cas
  procedure          :: report => report_wave_trial_cas
  procedure          :: sizeof => sizeof_wave_trial_cas

  procedure          :: trial_energy_full => trial_energy_cas
  procedure          :: trial_h2_gres_mean => trial_h2_gres_mean_cas
  procedure          :: initial_energy_full => initial_energy_cas
  procedure          :: initial_h2_gres_mean => initial_h2_gres_mean_cas

  procedure          :: get_ovlp_batch => get_ovlp_batch_cas
  procedure          :: block_h1_mean => h1_mean_batch_cas
  procedure          :: h2_gres_mean_batch => h2_gres_mean_batch_cas
  procedure          :: exchange_energy_batch => exchange_energy_batch_cas
  procedure          :: exchange_energy_batch_gres => exchange_energy_batch_gres_cas
  procedure          :: self_energy_batch => self_energy_batch_cas
  procedure          :: self_energy_batch_gres => self_energy_batch_gres_cas
  procedure          :: local_energy_batch => local_energy_batch_cas
  procedure          :: local_energy_batch_gres => local_energy_batch_gres_cas
end type WaveTrialCAS

interface WaveTrialCAS
  module procedure finit_wave_trial_cas
end interface WaveTrialCAS

contains

!******************************************************************************** 
!
! Initialization of the WaveTrialCAS object
!
!******************************************************************************** 
subroutine init_wave_trial_cas(wtdes, ham, coeff, ci_det, ci_coeff, self)
  type(WaveTrialDes), intent(in)            :: wtdes
  type(Hamilton), target, intent(in)        :: ham
  real(wp), intent(in)                      :: coeff(:,:,:)
  type(slater_det), allocatable, intent(in) :: ci_det(:)
  real(wp), allocatable, intent(in)         :: ci_coeff(:)
  type(WaveTrialCAS), intent(out)           :: self
  
  self%des = wtdes
  allocate(self%coeff, source=coeff)
  allocate(self%ci_coeff, source=ci_coeff)
  allocate(self%det, source=ci_det)

  call self%base_init(ham)
  allocate(self%hdes_cp, source=ham%des)

  call self%setup_coeff(self%ham%des)

!debug: wave_trial : this feature temporarily removed
  !call self%estimate_delta_x()

  self%hcas = CasHamilton(ham, self%des%nfrozen, self%des%nactive)
end subroutine init_wave_trial_cas


!******************************************************************************** 
!
! WaveTrialCAS constructor
!
!******************************************************************************** 
function finit_wave_trial_cas(wtdes, ham, coeff, ci_det, ci_coeff) result(self)
  type(WaveTrialDes), intent(in)            :: wtdes
  type(Hamilton), target, intent(in)        :: ham
  real(wp), intent(in)                      :: coeff(:,:,:)
  type(slater_det), allocatable, intent(in) :: ci_det(:)
  real(wp), allocatable, intent(in)         :: ci_coeff(:)
  type(WaveTrialCAS)                        :: self

  call init_wave_trial_cas(wtdes, ham, coeff, ci_det, ci_coeff, self)
end function finit_wave_trial_cas


!******************************************************************************** 
!
! Assignment operator for the WaveTrialCAS
!
!******************************************************************************** 
subroutine assign_wave_trial_cas(to, from)
  class(WaveTrialCAS), intent(inout) :: to
  class(WaveTrial), intent(in)       :: from

  select type (from)
    type is (WaveTrialCAS)
      call to%base_assign(from)
      to%WaveTrialSD = from%WaveTrialSD
      to%hcas = from%hcas

      if (allocated(from%det)) then
        if (allocated(to%det)) deallocate(to%det)
        allocate(to%det, source=from%det)
      end if

      if (allocated(from%ci_coeff)) then
        if (allocated(to%ci_coeff)) deallocate(to%ci_coeff)
        allocate(to%ci_coeff, source=from%ci_coeff)
      end if
    class default
      if (comm_world%mpirank == 0) write(*,*) "both dynamic types have to be WaveTrialSD in the assignment operator"
      stop
  end select
end subroutine assign_wave_trial_cas


!******************************************************************************** 
!
! Report WaveTrialCAS object
!
!******************************************************************************** 
subroutine report_wave_trial_cas(self, comm, fh)
  class(WaveTrialCAS), intent(inout) :: self
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
    write(fh%funit,100) "trial wave function type        ", "WaveTrialCAS"
    write(fh%funit,100) "size of the object              ", write_bytes(mem)
    write(fh%funit,101) "number of determinants          ", size(self%ci_coeff)
    write(fh%funit,101) "number of frozen orbitals       ", self%des%nfrozen
    write(fh%funit,101) "number of active orbitals       ", self%des%nactive
    write(fh%funit,101) "number of active electrons      ", self%des%nel_active
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
end subroutine report_wave_trial_cas


!******************************************************************************** 
!
! Estimate memory consumption by WaveTrialCAS object
!
!******************************************************************************** 
function sizeof_wave_trial_cas(self) result(mem)
  class(WaveTrialCAS), intent(in) :: self
  real(dp)                        :: mem

  mem = self%WaveTrialSD%sizeof()
  mem = mem + self%hcas%sizeof()
  mem = mem + sizeof(self%h1phi_occ) + sizeof(self%h1phi_s_occ) +  sizeof(self%h2phi_h_occ) + sizeof(self%h2phi_s_occ) + &
        sizeof(self%h2phi_x_chol_occ) + sizeof(self%h2phi_x_int_occ) + sizeof(self%det) + sizeof(self%ci_coeff)
end function sizeof_wave_trial_cas


!******************************************************************************** 
!
! Obtain occupied orbitals from all orbitals stored in self%coeff
!
! Note:
!    The coeff_occ stores explicitly both spin channels even in the 
!    restricted spin regime. This is needed since some of the determinants
!    will have  different spin-up and spind-down strings.
!
!******************************************************************************** 
subroutine setup_coeff_cas(self, hdes)
  class(WaveTrialCAS), intent(inout) :: self
  type(hamil_descriptor), intent(in) :: hdes
  !local variables
  integer                            :: spin, sp, nocc, w1, w2

  if (allocated(self%coeff_occ)) deallocate(self%coeff_occ)
  allocate(self%coeff_occ(hdes%n, sum(hdes%nel)))

  do spin = 1, size(hdes%nel)
    sp = min(spin, size(self%coeff, 3)) 
    nocc = hdes%nel(spin)
    w1 = hdes%nell(spin) + 1
    w2 = hdes%nell(spin) + nocc
    self%coeff_occ(:,w1:w2) = self%coeff(:,1:nocc,sp)
  end do
end subroutine setup_coeff_cas


!******************************************************************************** 
!
! Set active orbitals according to det (from self%ceoff to self%coeff_occ)
!
!******************************************************************************** 
subroutine setup_coeff_occ_cas(self, hdes, det)
  class(WaveTrialCAS), intent(inout) :: self
  type(hamil_descriptor), intent(in) :: hdes
  type(slater_det), intent(in)       :: det
  !local variables
  integer                            :: i, ii, nocc, indx, indx_, spin, sp
  integer, allocatable               :: orb(:,:)

  if (.not. allocated(self%coeff_occ)) then
    allocate(self%coeff_occ(hdes%n, sum(hdes%nel)))

    do spin = 1, size(hdes%nel)
      sp = min(spin, size(self%coeff, 3)) 
      nocc = hdes%nel(spin)
      indx = hdes%nell(spin) + 1
      indx_ = hdes%nell(spin) + nocc
      self%coeff_occ(:,indx:indx_) = self%coeff(:,1:nocc,sp)
    end do
  end if

  orb = det%content()

  do i = 1, size(orb, 1)
    indx = self%des%nfrozen + orb(i,1)
    spin = orb(i,2)
    sp = min(spin, size(self%coeff, 3))
    ii = i - self%des%nell_active(spin)
    indx_ = self%des%nfrozen + hdes%nell(spin) + ii
    self%coeff_occ(:,indx_) = self%coeff(:,indx,sp)
  end do
end subroutine setup_coeff_occ_cas


!******************************************************************************** 
!
! Set active orbitals (stored in Slater det.) from h1phi to h1phi_occ
! To check data layout of h1phi look at nosd.f: get_h1phi_occ
!
!******************************************************************************** 
subroutine setup_h1phi_occ_cas(self, hdes, det)
  class(WaveTrialCAS), intent(inout) :: self
  type(hamil_descriptor), intent(in) :: hdes
  type(slater_det), intent(in)       :: det
  !local variables
  integer                            :: i, ii, indx1, indx2, indx1_, indx2_, spin
  integer, allocatable               :: orb(:,:)

  if (.not. allocated(self%h1phi_occ)) then
    allocate(self%h1phi_occ(hdes%n, sum(hdes%nel)))

    do spin = 1, size(hdes%nel)
      indx1 = hdes%nell(spin) + 1
      indx2 = hdes%nell(spin) + hdes%nel(spin)
      indx1_ = hdes%nell_space(spin) + 1
      indx2_ = hdes%nell_space(spin) + hdes%nel(spin)
      self%h1phi_occ(:,indx1:indx2) = self%h1phi(:,indx1_:indx2_)
    end do
  end if

  orb = det%content()

  do i = 1, size(orb, 1)
    spin = orb(i,2)
    ii = i - self%des%nell_active(spin)
    indx1 = hdes%nell(spin) + self%des%nfrozen + ii
    indx1_ = hdes%nell_space(spin) + self%des%nfrozen + orb(i,1)
    self%h1phi_occ(:,indx1) = self%h1phi(:,indx1_)
  end do
end subroutine setup_h1phi_occ_cas


!******************************************************************************** 
!
! Set active orbitals (stored in Slater det.) from h1phi_s to h1phi_s_occ
! To check data layout of h1phi look at nosd.f: get_h1phi_occ
!
!******************************************************************************** 
subroutine setup_h1phi_s_occ_cas(self, hdes, det)
  class(WaveTrialCAS), intent(inout) :: self
  type(hamil_descriptor), intent(in) :: hdes
  type(slater_det), intent(in)       :: det
  !local variables
  integer                            :: i, ii, indx1, indx2, indx1_, indx2_, spin
  integer, allocatable               :: orb(:,:)

  if (.not. allocated(self%h1phi_s_occ)) then
    allocate(self%h1phi_s_occ(hdes%n, sum(hdes%nel)))

    do spin = 1, size(hdes%nel)
      indx1 = hdes%nell(spin) + 1
      indx2 = hdes%nell(spin) + hdes%nel(spin)
      indx1_ = hdes%nell_space(spin) + 1
      indx2_ = hdes%nell_space(spin) + hdes%nel(spin)
      self%h1phi_s_occ(:,indx1:indx2) = self%h1phi_s(:,indx1_:indx2_)
    end do
  end if

  orb = det%content()

  do i = 1, size(orb, 1)
    spin = orb(i,2)
    ii = i - self%des%nell_active(spin)
    indx1 = hdes%nell(spin) + self%des%nfrozen + ii
    indx1_ = hdes%nell_space(spin) + self%des%nfrozen + orb(i,1)
    self%h1phi_s_occ(:,indx1) = self%h1phi_s(:,indx1_)
  end do
end subroutine setup_h1phi_s_occ_cas


!******************************************************************************** 
!
! Set active orbitals (stored in Slater det.) from h2phi_h to h2phi_h_occ
! To check data layout of h2phi_h look at nosd.f: get_h2phi_h_occ
!
!******************************************************************************** 
subroutine setup_h2phi_h_occ_cas(self, hdes, det)
  class(WaveTrialCAS), intent(inout) :: self
  type(hamil_descriptor), intent(in) :: hdes
  type(slater_det), intent(in)       :: det
  !local variables
  integer                            :: i, ii, n, ng, indx1, indx2, indx1_, indx2_, spin
  integer, allocatable               :: orb(:,:)

  n = hdes%n
  ng = size(self%h2phi_h, 1)

  if (.not. allocated(self%h2phi_h_occ)) then
    allocate(self%h2phi_h_occ(ng, n*sum(hdes%nel)))

    do spin = 1, size(hdes%nel)
      indx1 = n*hdes%nell(spin) + 1
      indx2 = n*hdes%nell(spin) + n*hdes%nel(spin)
      indx1_ = n*hdes%nell_space(spin) + 1
      indx2_ = n*hdes%nell_space(spin) + n*hdes%nel(spin)
      self%h2phi_h_occ(:,indx1:indx2) = self%h2phi_h(:,indx1_:indx2_)
    end do
  end if

  orb = det%content()

  do i = 1, size(orb, 1)
    spin = orb(i,2)
    ii = i - self%des%nell_active(spin)
    indx1 = (hdes%nell(spin) + self%des%nfrozen + ii - 1) * n + 1
    indx2 = (hdes%nell(spin) + self%des%nfrozen + ii) * n
    indx1_ = (hdes%nell_space(spin) + self%des%nfrozen + orb(i,1) - 1) * n + 1
    indx2_ = (hdes%nell_space(spin) + self%des%nfrozen + orb(i,1)) * n
    self%h2phi_h_occ(:,indx1:indx2) = self%h2phi_h(:,indx1_:indx2_)
  end do
end subroutine setup_h2phi_h_occ_cas


!******************************************************************************** 
!
! Set active orbitals (stored in Slater det.) from h2phi_h to h2phi_h_occ
! To check data layout of h2phi_h look at nosd.f: get_h2phi_h_occ
!
!******************************************************************************** 
subroutine setup_h2phi_s_occ_cas(self, hdes, det)
  class(WaveTrialCAS), intent(inout) :: self
  type(hamil_descriptor), intent(in) :: hdes
  type(slater_det), intent(in)       :: det
  !local variables
  integer                            :: i, ii, n, ng, indx1, indx2, indx1_, indx2_, spin
  integer, allocatable               :: orb(:,:)

  n = hdes%n
  ng = size(self%h2phi_s, 1)

  if (.not. allocated(self%h2phi_s_occ)) then
    allocate(self%h2phi_s_occ(ng, n*sum(hdes%nel)))

    do spin = 1, size(hdes%nel)
      indx1 = n*hdes%nell(spin) + 1
      indx2 = n*hdes%nell(spin) + n*hdes%nel(spin)
      indx1_ = n*hdes%nell_space(spin) + 1
      indx2_ = n*hdes%nell_space(spin) + n*hdes%nel(spin)
      self%h2phi_s_occ(:,indx1:indx2) = self%h2phi_s(:,indx1_:indx2_)
    end do
  end if

  orb = det%content()

  do i = 1, size(orb, 1)
    spin = orb(i,2)
    ii = i - self%des%nell_active(spin)
    indx1 = (hdes%nell(spin) + self%des%nfrozen + ii - 1) * n + 1
    indx2 = (hdes%nell(spin) + self%des%nfrozen + ii) * n
    indx1_ = (hdes%nell_space(spin) + self%des%nfrozen + orb(i,1) - 1) * n + 1
    indx2_ = (hdes%nell_space(spin) + self%des%nfrozen + orb(i,1)) * n
    self%h2phi_s_occ(:,indx1:indx2) = self%h2phi_s(:,indx1_:indx2_)
  end do
end subroutine setup_h2phi_s_occ_cas


!******************************************************************************** 
!
! Set active orbitals in det from h2phi_x_chol to h2phi_x_chol_occ
! To check data layout of h2phi_x_chol look at nosd.f: get_h2phi_x_chol_occ
!
!******************************************************************************** 
subroutine setup_h2phi_x_occ_cas(self, hdes, det)
  class(WaveTrialCAS), intent(inout) :: self
  type(hamil_descriptor), intent(in) :: hdes
  type(slater_det), intent(in)       :: det

  select case (self%des%exchange_mode)
    case ("cholesky")
      call self%setup_h2phi_x_occ_chol_cas(hdes, det)
    case ("eri")
      call self%setup_h2phi_x_occ_int_cas(hdes, det)
  end select
end subroutine setup_h2phi_x_occ_cas


!******************************************************************************** 
!
! Set active orbitals in det from h2phi_x_chol to h2phi_x_chol_occ
! To check data layout of h2phi_x_chol look at nosd.f: get_h2phi_x_chol_occ
!
!******************************************************************************** 
subroutine setup_h2phi_x_occ_chol_cas(self, hdes, det)
  class(WaveTrialCAS), intent(inout) :: self
  type(hamil_descriptor), intent(in) :: hdes
  type(slater_det), intent(in)       :: det
  !local variables
  integer                            :: i, ii, g, n, ng, indx1, indx2, indx1_, indx2_, spin
  integer, allocatable               :: orb(:,:)

  n = hdes%n
  ng = size(self%h2phi_x_chol, 2) / sum(hdes%nel_space)

  if (.not. allocated(self%h2phi_x_chol_occ)) then
    allocate(self%h2phi_x_chol_occ(n, ng*sum(hdes%nel)))

    do g = 1, ng
      do spin = 1, size(hdes%nel)
        indx1 = ng*hdes%nell(spin) + (g-1)*hdes%nel(spin) + 1
        indx2 = ng*hdes%nell(spin) + g*hdes%nel(spin)
        indx1_ = ng*hdes%nell_space(spin) + (g-1)*hdes%nel_space(spin) + 1
        indx2_ = ng*hdes%nell_space(spin) + (g-1)*hdes%nel_space(spin) + hdes%nel(spin)
        self%h2phi_x_chol_occ(:,indx1:indx2) = self%h2phi_x_chol(:,indx1_:indx2_)
      end do
    end do
  end if

  orb = det%content()

  do g = 1, ng
    do i = 1, size(orb, 1)
      spin = orb(i,2)
      ii = i - self%des%nell_active(spin)
      indx1 = ng*hdes%nell(spin) + (g-1)*hdes%nel(spin) + self%des%nfrozen + ii
      indx1_ = ng*hdes%nell_space(spin) + (g-1)*hdes%nel_space(spin) + self%des%nfrozen + orb(i,1)
      self%h2phi_x_chol_occ(:,indx1) = self%h2phi_x_chol(:,indx1_)
    end do
  end do
end subroutine setup_h2phi_x_occ_chol_cas


!******************************************************************************** 
!
! Set active orbitals in det from h2phi_x_int to h2phi_x_int_occ
! To check data layout of h2phi_x_int look at nosd.f: get_h2phi_x_int_occ
!
!******************************************************************************** 
subroutine setup_h2phi_x_occ_int_cas(self, hdes, det)
  class(WaveTrialCAS), intent(inout) :: self
  type(hamil_descriptor), intent(in) :: hdes
  type(slater_det), intent(in)       :: det
  !local variables
  integer                            :: i, ii, j, jj, n, noccmax
  integer                            :: indx1, indx2, indx3, indx4, indx1_, indx2_, indx3_, indx4_
  integer                            :: spin, spin1, spin2, sp, ispin
  integer, allocatable               :: orb(:,:)

  n = hdes%n
  ispin = size(self%h2phi_x_int, 3)
  noccmax = maxval(hdes%nel)

  if (.not. allocated(self%h2phi_x_int_occ)) then
    allocate(self%h2phi_x_int_occ(n*noccmax, n*noccmax, size(hdes%nel)))

    do spin = 1, size(hdes%nel)
      sp = min(spin, ispin)
      indx1 = 1
      indx2 = n*hdes%nel(spin)
      self%h2phi_x_int_occ(indx1:indx2,indx1:indx2,spin) = self%h2phi_x_int(indx1:indx2,indx1:indx2,sp)
    end do
  end if

  orb = det%content()

  do j = 1, size(orb,1)
    spin2 = orb(j,2)
    sp = min(spin2, ispin)
    jj = j - self%des%nell_active(spin2)
    indx3 = n*(self%des%nfrozen + jj - 1) + 1
    indx4 = n*(self%des%nfrozen + jj)
    indx3_ = n*(self%des%nfrozen + orb(j,1) - 1) + 1
    indx4_ = n*(self%des%nfrozen + orb(j,1))

    indx1 = 1
    indx2 = n * self%des%nfrozen
    self%h2phi_x_int_occ(indx3:indx4,indx1:indx2,spin2) = self%h2phi_x_int(indx3_:indx4_,indx1:indx2,sp)
    self%h2phi_x_int_occ(indx1:indx2,indx3:indx4,spin2) = self%h2phi_x_int(indx1:indx2,indx3_:indx4_,sp)

    do i = 1, size(orb,1)
      spin1 = orb(i,2)
      ii = i - self%des%nell_active(spin1) 
      if (spin1 /= spin2) cycle
      indx1 = n*(self%des%nfrozen + ii - 1) + 1
      indx2 = n*(self%des%nfrozen + ii)
      indx1_ = n*(self%des%nfrozen + orb(i,1) - 1) + 1
      indx2_ = n*(self%des%nfrozen + orb(i,1))
      self%h2phi_x_int_occ(indx1:indx2,indx3:indx4,spin1) = self%h2phi_x_int(indx1_:indx2_,indx3_:indx4_,sp)
    end do
  end do
end subroutine setup_h2phi_x_occ_int_cas


!******************************************************************************** 
!
! Estimate trial energy: 
!
!    E_T = <Psi_T|H|Psi_T> / <Psi_T|Psi_T>
!        = \sum_{ij} c_i* c_j <Phi_i|H|Phi_j>  / \sum_i |c_i|^2
!
!******************************************************************************** 
subroutine trial_energy_cas(self, e1, lgmean, eh, ex, eself)
  class(WaveTrialCAS), intent(inout)    :: self
  complex(wp), intent(out)              :: e1
  complex(wp), allocatable, intent(out) :: lgmean(:)
  complex(wp), allocatable, intent(out) :: eh(:)
  complex(wp), allocatable, intent(out) :: ex(:)
  complex(wp), allocatable, intent(out) :: eself(:)
  !local variables
  integer                               :: i, j, ng
  real(wp)                              :: den, weight, e1_
  real(wp), allocatable                 :: lgmean_(:), eh_(:), ex_(:), eself_(:)
  real(wp), allocatable                 :: eh_af(:), ex_af(:)
  real(wp), pointer                     :: h2_gres_ob(:,:,:,:)

  ng = size(self%hcas%lgmean_frozen)

  allocate(lgmean(ng))
  allocate(eh(ng))
  allocate(ex(ng))
  allocate(eself(ng))

  call self%hcas%get_h2_gres_ob(h2_gres_ob)

  den = zeror
  e1 = zeroc
  lgmean = zeroc
  eh = zeroc
  ex = zeroc
  eself = zeroc

  do i = 1, size(self%det)
    do j = 1, size(self%det)
      weight = self%ci_coeff(i) * self%ci_coeff(j)

      call slater_condon_one(self%det(i), self%det(j), self%hcas%h1, e1_)
      call slater_condon_one_array(self%det(i), self%det(j), h2_gres_ob, lgmean_)
      call slater_condon_two_gres(self%det(i), self%det(j), h2_gres_ob, eh_, ex_)
      call slater_condon_one_array(self%det(i), self%det(j), self%hcas%hself_gres, eself_)

      !active-frozen interaction energy
      call slater_condon_one_array(self%det(i), self%det(j), self%hcas%hpot_gres, eh_af)
      call slater_condon_one_array(self%det(i), self%det(j), self%hcas%xpot_gres, ex_af)

      den = den + weight * self%det(i)%ovlp(self%det(j))
      e1 = e1 + weight * e1_
      lgmean = lgmean + weight * lgmean_
      eh = eh + weight * (eh_ + eh_af)
      ex = ex + weight * (ex_ + ex_af)
      eself = eself + weight * eself_
    end do
  end do

  e1 = e1/den + self%hcas%e1_frozen
  lgmean = lgmean/den + self%hcas%lgmean_frozen
  eh = eh/den + self%hcas%eh_frozen
  ex = ex/den + self%hcas%ex_frozen
  eself = eself/den + self%hcas%eself_frozen
end subroutine trial_energy_cas


!******************************************************************************** 
!
! Calculate trial exectation vlaue of the Cholesky vectors
!
!    l_g = <Psi_T | L_g | Psi_T> / <Psi_T | Psi_T>
!        = \sum_ij c_i* c_j <Phi_i | L_g | Phi_j> / \sum_i |c_i|^2
!
! Assumes MO basis and uses Slater-Condon rules
!
!******************************************************************************** 
subroutine trial_h2_gres_mean_cas(self, lgmean)
  class(WaveTrialCAS), intent(inout)    :: self
  complex(wp), allocatable, intent(out) :: lgmean(:)
  !local variables
  integer                               :: i, j
  real(wp)                              :: den, weight
  real(wp), allocatable                 :: lgmean_(:)
  real(wp), pointer                     :: h2_gres_ob(:,:,:,:)

  allocate(lgmean(size(self%hcas%lgmean_frozen)))

  den = zeror
  lgmean = zeroc

  call self%hcas%get_h2_gres_ob(h2_gres_ob)

  do i = 1, size(self%det)
    do j = 1, size(self%det)
      weight = self%ci_coeff(i) * self%ci_coeff(j)
      call slater_condon_one_array(self%det(i), self%det(j), h2_gres_ob, lgmean_)

      den = den + weight * self%det(i)%ovlp(self%det(j))
      lgmean = lgmean + weight * lgmean_
    end do
  end do

  lgmean = lgmean/den + self%hcas%lgmean_frozen
end subroutine trial_h2_gres_mean_cas


!******************************************************************************** 
!
! Estimate initial energy: 
!
!    E_I = <Psi_T|H|Phi_1> / <Psi_T|Phi_1>
!        = \sum_{d} c_d* <Phi_d|H|Phi_1>  / c_1*
!
!******************************************************************************** 
subroutine initial_energy_cas(self, e1, lgmean, eh, ex, eself)
  class(WaveTrialCAS), intent(inout)    :: self
  complex(wp), intent(out)              :: e1
  complex(wp), allocatable, intent(out) :: lgmean(:)
  complex(wp), allocatable, intent(out) :: eh(:)
  complex(wp), allocatable, intent(out) :: ex(:)
  complex(wp), allocatable, intent(out) :: eself(:)
  !local variables
  integer                               :: d, ng
  real(wp)                              :: den, e1_
  real(wp), allocatable                 :: lgmean_(:), eh_(:), ex_(:), eself_(:)
  real(wp), allocatable                 :: eh_af(:), ex_af(:)
  real(wp), pointer                     :: h2_gres_ob(:,:,:,:)

  ng = size(self%hcas%lgmean_frozen)

  allocate(lgmean(ng))
  allocate(eh(ng))
  allocate(ex(ng))
  allocate(eself(ng))

  den = zeror
  e1 = zeroc
  lgmean = zeroc
  eh = zeroc
  ex = zeroc
  eself = zeroc

  call self%hcas%get_h2_gres_ob(h2_gres_ob)

  do d = 1, size(self%det)
    call slater_condon_one(self%det(d), self%det(1), self%hcas%h1, e1_)
    call slater_condon_one_array(self%det(d), self%det(1), h2_gres_ob, lgmean_)
    call slater_condon_two_gres(self%det(d), self%det(1), h2_gres_ob, eh_, ex_)
    call slater_condon_one_array(self%det(d), self%det(1), self%hcas%hself_gres, eself_)

    !frozen-active components
    call slater_condon_one_array(self%det(d), self%det(1), self%hcas%hpot_gres, eh_af)
    call slater_condon_one_array(self%det(d), self%det(1), self%hcas%xpot_gres, ex_af)

    den = den + self%ci_coeff(d) * self%det(d)%ovlp(self%det(1))
    e1 = e1 + self%ci_coeff(d) *  e1_
    lgmean = lgmean + self%ci_coeff(d) * lgmean_
    eself = eself + self%ci_coeff(d) * eself_
    eh = eh + self%ci_coeff(d) * (eh_ + eh_af)
    ex = ex + self%ci_coeff(d) * (ex_ + ex_af)
  end do

  e1 = e1/den + self%hcas%e1_frozen
  lgmean = lgmean/den + self%hcas%lgmean_frozen
  eh = eh/den + self%hcas%eh_frozen
  ex = ex/den + self%hcas%ex_frozen
  eself = eself/den + self%hcas%eself_frozen
end subroutine initial_energy_cas


!******************************************************************************** 
!
! Calculates initial exectation vlaue of the Cholesky vectors
!
!    <L_g> = <Psi_T|L_g|Phi_1> / <Psi_T|Phi_1>
!          = \sum_d c_d* <Phi_d|L_g|Phi_1> / c_1*
!
!******************************************************************************** 
subroutine initial_h2_gres_mean_cas(self, lgmean)
  class(WaveTrialCAS), intent(inout)    :: self
  complex(wp), allocatable, intent(out) :: lgmean(:)
  !local variables
  integer                               :: d
  real(wp)                              :: den
  real(wp), allocatable                 :: lgmean_(:)
  real(wp), pointer                     :: h2_gres_ob(:,:,:,:)

  allocate(lgmean(size(self%hcas%lgmean_frozen)))

  den = zeror
  lgmean = zeroc

  call self%hcas%get_h2_gres_ob(h2_gres_ob)

  do d = 1, size(self%det)
    call slater_condon_one_array(self%det(d), self%det(1), h2_gres_ob, lgmean_)

    den = den + self%ci_coeff(d) * self%det(d)%ovlp(self%det(1))
    lgmean = lgmean + self%ci_coeff(d) * lgmean_
  end do
  
  lgmean  = lgmean/den + self%hcas%lgmean_frozen
end subroutine initial_h2_gres_mean_cas


!******************************************************************************** 
!
!  Determines ispin and rspin values for given determinant
!
!      spin-polarized case                     : ispin=2, rspin=1
!      No spin-polarization and a_str /= b_str : ispin=2, rspin=1 
!      No spin polarization and a_str == b_str : ispin=1, rspin=2
!  
!******************************************************************************** 
subroutine get_det_ispin(self, det)
  class(WaveTrialCAS), intent(inout) :: self
  type(slater_det), intent(in)       :: det
  
  if (det%str(1)==det%str(2) .and. self%ham%des%ispin==1) then
    self%hdes_cp%ispin = 1
    self%hdes_cp%rspin = 2
    self%hdes_cp%nocc_ = self%ham%des%nel(1)
  else
    self%hdes_cp%ispin = 2
    self%hdes_cp%rspin = 1
    self%hdes_cp%nocc_ = sum(self%ham%des%nel)
  end if
end subroutine get_det_ispin


!******************************************************************************** 
!
! Calculate overlap between trial determinant and walker determinants
!
!    <Psi_T|Phi> = \sum_i c_i* <Psi_i|Phi>
!
! where <Psi_i|Phi> is the overlap between two non-orthogonal
!
! Actual low-level implementation in nosd.f: get_ovlp_batch
!
!******************************************************************************** 
subroutine get_ovlp_batch_cas(self, coeff, ovlp)
  class(WaveTrialCAS), intent(inout)    :: self
  complex(wp), intent(in)               :: coeff(:,:)
  complex(wp), allocatable, intent(out) :: ovlp(:)
  !local variables
  integer                               :: d, nw
  complex(wp), allocatable              :: ovlp_(:)
  character(len=*), parameter           :: proc_name = "get_ovlp_cas"

  if (profile_code) call start_profiling(proc_name)

  nw = size(coeff, 2) / self%ham%des%nocc

  allocate(ovlp(nw))
  ovlp = zeroc

  do d = 1, size(self%det)
    call self%get_ispin(self%det(d))
    call self%setup_coeff_occ(self%hdes_cp, self%det(d))

    call get_ovlp_batch(self%hdes_cp, self%coeff_occ, coeff, ovlp_)
    ovlp = ovlp + self%ci_coeff(d) * ovlp_
  end do 

  if (profile_code) call end_profiling(proc_name)
end subroutine get_ovlp_batch_cas


!******************************************************************************** 
!
! Calculate one-electron energy for a batch of walkers 
!
!    E1 = <Psi_T | H1 | Phi_w> / <Psi_T | Phi_w>
!       = \sum_d c_d* <Phi_d | H1 | Phi_w>  / \sum_d c_d* <Phi_d | Phi_w>
!
! Actual implementation of <Phi_d | H1 | Phi_w> / <Phi_d | Phi_w> in
! nosd.f: h1_mean_no_batch
!            
!******************************************************************************** 
subroutine h1_mean_batch_cas(self, coeff, e1)
  class(WaveTrialCAS), intent(inout)    :: self
  complex(wp), intent(in)               :: coeff(:,:)
  complex(wp), allocatable, intent(out) :: e1(:)
  !local variables
  integer                               :: nw, d
  complex(wp), allocatable              :: ovlp(:), e1_(:), den(:)
  complex(wp), allocatable              :: coeff_bi(:,:)
  character(len=*), parameter           :: proc_name = "h1_mean_cas"

  if (profile_code) call start_profiling(proc_name)

  nw = size(coeff, 2) / self%ham%des%nocc

  allocate(e1(nw))
  e1 = zeroc

  allocate(den(nw))
  den = zeroc

  do d = 1, size(self%det)
    call self%get_ispin(self%det(d))
    call self%setup_coeff_occ(self%hdes_cp, self%det(d))

    call get_ovlp_batch(self%hdes_cp, self%coeff_occ, coeff, ovlp)
    ovlp = self%ci_coeff(d) * ovlp
    den = den + ovlp

    call biorth_orb_batch(self%hdes_cp, self%coeff_occ, coeff, coeff_bi)

    call self%setup_h1phi_occ(self%hdes_cp, self%det(d))
    call h1_mean_no_batch(self%hdes_cp, self%h1phi_occ, coeff_bi, e1_)
    e1 = e1 + ovlp * e1_
  end do

  e1 = e1 / den

  if (profile_code) call end_profiling(proc_name)
end subroutine h1_mean_batch_cas


!******************************************************************************** 
!
! Calculate expectation value of the Cholesky vectors for a batch of walkers
!
!    l_{g,w} = <Psi_T |L_g |Phi_w> / <Psi_T | Phi_w>
!            = \sum_d c_d* <Phi_d | L_g | Phi_w> / c_d*
!
! Actual implementation of <Phi_d | L_g | Phi_w> / <Phi_d | Phi_w> in
! nosd.f: h2_gres_mean_no_batch
!
!******************************************************************************** 
subroutine h2_gres_mean_batch_cas(self, coeff, lgmean, eh)
  class(WaveTrialCAS), intent(inout)    :: self
  complex(wp), intent(in)               :: coeff(:,:)
  complex(wp), allocatable, intent(out) :: lgmean(:,:)
  complex(wp), allocatable, intent(out) :: eh(:,:)
  !local variables
  integer                               :: ng, nw, w, d
  complex(wp), allocatable              :: ovlp(:), lgmean_(:,:), den(:)
  complex(wp), allocatable              :: coeff_bi(:,:)
  character(len=*), parameter           :: proc_name = "h2_gres_mean_cas"

  if (profile_code) call start_profiling(proc_name)

  ng = size(self%hcas%lgmean_frozen)
  nw = size(coeff, 2) / self%ham%des%nocc

  allocate(lgmean(ng,nw))
  lgmean = zeroc

  allocate(eh(ng,nw))
  eh = zeroc

  allocate(den(nw))
  den = zeroc

  do d = 1, size(self%det)
    call self%get_ispin(self%det(d))
    call self%setup_coeff_occ(self%hdes_cp, self%det(d))

    call get_ovlp_batch(self%hdes_cp, self%coeff_occ, coeff, ovlp)
    ovlp = self%ci_coeff(d) * ovlp
    den = den + ovlp

    call biorth_orb_batch(self%hdes_cp, self%coeff_occ, coeff, coeff_bi)

    call self%setup_h2phi_h_occ(self%hdes_cp, self%det(d))
    call h2_gres_mean_no_batch(self%hdes_cp, self%h2phi_h_occ, coeff_bi, lgmean_, self%des%h_block_size)

    do w = 1, nw
      lgmean(:,w) = lgmean(:,w) + ovlp(w) * lgmean_(:,w)
      eh(:,w) = eh(:,w) + 0.5_wp * ovlp(w) * lgmean_(:,w)**2
    end do
  end do

  do w = 1, nw
    lgmean(:,w) = lgmean(:,w) / den(w)
    eh(:,w) = eh(:,w) / den(w)
  end do

  if (profile_code) call end_profiling(proc_name)
end subroutine h2_gres_mean_batch_cas


!******************************************************************************** 
!
! Calculate exchange energy for the batch of walkers
!
!    Ex_w = \sum_d c_d <Psi_d | Hx | Phi_w> / \sum_d c_d <Psi_d | Phi_w>
! 
! exchange_calculator chooses the method to calculate exchange energy
!
! Actual low-level implementations in nosd.f
!
!******************************************************************************** 
subroutine exchange_energy_batch_cas(self, coeff, ex)
  class(WaveTrialCAS), intent(inout)    :: self
  complex(wp), intent(in)               :: coeff(:,:)
  complex(wp), allocatable, intent(out) :: ex(:)
  !local variables
  integer                               :: nw, d
  complex(wp), allocatable              :: den(:), ovlp(:), ex_(:)
  complex(wp), allocatable              :: coeff_bi(:,:)
  character(len=*), parameter           :: proc_name = "exchange_cas"
  
  if (profile_code) call start_profiling(proc_name)

  nw = size(coeff,2) / self%ham%des%nocc

  allocate(ex(nw))
  ex = zeroc

  allocate(den(nw))
  den = zeroc

  do d = 1, size(self%det)
    call self%get_ispin(self%det(d))
    call self%setup_coeff_occ(self%hdes_cp, self%det(d))

    call get_ovlp_batch(self%hdes_cp, self%coeff_occ, coeff, ovlp)
    ovlp = self%ci_coeff(d) * ovlp
    den = den + ovlp

    call biorth_orb_batch(self%hdes_cp, self%coeff_occ, coeff, coeff_bi)

    call self%setup_h2phi_x_occ(self%hdes_cp, self%det(d))
    call self%exchange_calculator(self%hdes_cp, coeff_bi, ex_)

    ex = ex + ovlp * ex_
  end do

  ex = ex / den

!debug: delta_x still not functional with WaveTrialCAS
  !ex = ex + self%des%delta_x

  if (profile_code) call end_profiling(proc_name)
end subroutine exchange_energy_batch_cas


!******************************************************************************** 
!
! Calculate exchange energy for the batch of walkers
!
!    Ex_{g,w} = \sum_d c_d <Psi_d | Hx_g | Phi_w> / \sum_d c_d <Psi_d | Phi_w>
! 
! exchange_calculator chooses the method to calculate exchange energy
!
! Actual low-level implementations in nosd.f
!
!******************************************************************************** 
subroutine exchange_energy_batch_gres_cas(self, coeff, ex)
  class(WaveTrialCAS), intent(inout)    :: self
  complex(wp), intent(in)               :: coeff(:,:)
  complex(wp), allocatable, intent(out) :: ex(:,:)
  !local variables
  integer                               :: ng, nw, w, d
  complex(wp), allocatable              :: den(:), ovlp(:), ex_(:,:)
  complex(wp), allocatable              :: coeff_bi(:,:)
  character(len=*), parameter           :: proc_name = "exchange_cas"
  
  if (profile_code) call start_profiling(proc_name)

  ng = size(self%hcas%lgmean_frozen)
  nw = size(coeff,2) / self%ham%des%nocc

  allocate(ex(ng, nw))
  ex = zeroc

  allocate(den(nw))
  den = zeroc

  do d = 1, size(self%det)
    call self%get_ispin(self%det(d))
    call self%setup_coeff_occ(self%hdes_cp, self%det(d))

    call get_ovlp_batch(self%hdes_cp, self%coeff_occ, coeff, ovlp)
    ovlp = self%ci_coeff(d) * ovlp
    den = den + ovlp

    call biorth_orb_batch(self%hdes_cp, self%coeff_occ, coeff, coeff_bi)

    call self%setup_h2phi_x_occ(self%hdes_cp, self%det(d))
    call self%exchange_calculator_gres(self%hdes_cp, coeff_bi, ex_)

    do w = 1, nw
      ex(:,w) = ex(:,w) + ovlp(w) * ex_(:,w)
    end do
  end do

  do w = 1, nw
    ex(:,w) = ex(:,w) / den(w)
  end do

!debug: delta_x still not functional with WaveTrialCAS
  !ex = ex + self%des%delta_x

  if (profile_code) call end_profiling(proc_name)
end subroutine exchange_energy_batch_gres_cas


!******************************************************************************** 
!
! Setup the method for the exchange energy evaluation
!
! Note: input orbitals are biorthogonal orbitals
!
!******************************************************************************** 
subroutine exchange_calculator_cas(self, hdes, coeff_bi, ex)
  class(WaveTrialCAS), intent(inout)    :: self
  type(hamil_descriptor), intent(in)    :: hdes
  complex(wp), intent(in)               :: coeff_bi(:,:)
  complex(wp), allocatable, intent(out) :: ex(:)
  !local variables
  integer                               :: w
  complex(wp), allocatable              :: ex_(:,:)

  select case (self%des%exchange_mode)
    case ("cholesky")
      call exchange_energy_no_chol_batch(hdes, self%h2phi_x_chol_occ, coeff_bi, ex_, self%des%x_block_size)
      allocate(ex(size(ex_, 2)))
      do w = 1, size(ex_, 2)
        ex(w) = sum(ex_(:,w))
      end do
    case ("eri")
      call exchange_energy_no_int_batch(hdes, self%h2phi_x_int_occ , coeff_bi, ex, self%des%x_block_size)
  end select
end subroutine exchange_calculator_cas


!******************************************************************************** 
!
! Setup the method for the g-resolved exchange energy evaluation
!
! Note: input orbitals are biorthogonal orbitals
!
!******************************************************************************** 
subroutine exchange_calculator_gres_cas(self, hdes, coeff_bi, ex)
  class(WaveTrialCAS), intent(inout)    :: self
  type(hamil_descriptor), intent(in)    :: hdes
  complex(wp), intent(in)               :: coeff_bi(:,:)
  complex(wp), allocatable, intent(out) :: ex(:,:)

  call exchange_energy_no_chol_batch(hdes, self%h2phi_x_chol, coeff_bi, ex, self%des%x_block_size)
end subroutine exchange_calculator_gres_cas


!******************************************************************************** 
!
! Calculate self energy for a batch of walkers
!
!    Eself_{w} = \sum_g <Psi_T | L_g^2 | Phi_w> / <Psi_T | Phi_w>
!              = \sum_g \sum_d c_d* <Phi_d|L_g^2| Phi_w> / c_d*
!
! Actual implementation of <Phi_d | L_g^2 | Phi_w> / <Phi_d | Phi_w> in
! nosd.f: h2_gres_mean_no_batch
!
!******************************************************************************** 
subroutine self_energy_batch_cas(self, coeff, eself)
  class(WaveTrialCAS), intent(inout)    :: self
  complex(wp), intent(in)               :: coeff(:,:)
  complex(wp), allocatable, intent(out) :: eself(:)
  !local variables
  integer                               :: nw, d
  complex(wp), allocatable              :: ovlp(:), eself_(:), den(:)
  complex(wp), allocatable              :: coeff_bi(:,:)
  character(len=*), parameter           :: proc_name = "self_energy_cas"

  if (profile_code) call start_profiling(proc_name)

  nw = size(coeff, 2) / self%ham%des%nocc

  allocate(eself(nw))
  eself = zeroc

  allocate(den(nw))
  den = zeroc

  do d = 1, size(self%det)
    call self%get_ispin(self%det(d))
    call self%setup_coeff_occ(self%hdes_cp, self%det(d))

    call get_ovlp_batch(self%hdes_cp, self%coeff_occ, coeff, ovlp)
    ovlp = self%ci_coeff(d) * ovlp
    den = den + ovlp

    call biorth_orb_batch(self%hdes_cp, self%coeff_occ, coeff, coeff_bi)

    call self%setup_h1phi_s_occ(self%hdes_cp, self%det(d))
    call h1_mean_no_batch(self%hdes_cp, self%h1phi_s_occ, coeff_bi, eself_)

    eself = eself + ovlp * eself_
  end do

  eself = eself / den

  if (profile_code) call end_profiling(proc_name)
end subroutine self_energy_batch_cas


!******************************************************************************** 
!
! Calculate g-resolved self energy for a batch of walkers
!
!    Eself_{g,w} = <Psi_T | L_g^2 | Phi_w> / <Psi_T | Phi_w>
!                = \sum_d c_d* <Phi_d|L_g^2| Phi_w> / c_d*
!
! Actual implementation of <Phi_d | L_g^2 | Phi_w> / <Phi_d | Phi_w> in
! nosd.f: h2_gres_mean_no_batch
!
!******************************************************************************** 
subroutine self_energy_batch_gres_cas(self, coeff, eself)
  class(WaveTrialCAS), intent(inout)    :: self
  complex(wp), intent(in)               :: coeff(:,:)
  complex(wp), allocatable, intent(out) :: eself(:,:)
  !local variables
  integer                               :: ng, nw, w, d
  complex(wp), allocatable              :: ovlp(:), eself_(:,:), den(:)
  complex(wp), allocatable              :: coeff_bi(:,:)
  character(len=*), parameter           :: proc_name = "self_energy_gres_cas"

  if (profile_code) call start_profiling(proc_name)

  ng = size(self%hcas%lgmean_frozen)
  nw = size(coeff, 2) / self%ham%des%nocc

  allocate(eself(ng,nw))
  eself = zeroc

  allocate(den(nw))
  den = zeroc

  do d = 1, size(self%det)
    call self%get_ispin(self%det(d))
    call self%setup_coeff_occ(self%hdes_cp, self%det(d))

    call get_ovlp_batch(self%hdes_cp, self%coeff_occ, coeff, ovlp)
    ovlp = self%ci_coeff(d) * ovlp
    den = den + ovlp

    call biorth_orb_batch(self%hdes_cp, self%coeff_occ, coeff, coeff_bi)

    call self%setup_h2phi_s_occ(self%hdes_cp, self%det(d))
    call h2_gres_mean_no_batch(self%hdes_cp, self%h2phi_s_occ, coeff_bi, eself_, self%des%h_block_size)

    do w = 1, nw
      eself(:,w) = eself(:,w) + ovlp(w) * eself_(:,w)
    end do
  end do

  do w = 1, nw
    eself(:,w) = eself(:,w) / den(w)
  end do

  if (profile_code) call end_profiling(proc_name)
end subroutine self_energy_batch_gres_cas


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
subroutine local_energy_batch_cas(self, coeff, e1, lgmean, eh, ex, eself, lsamplex, lself, en)
  class(WaveTrialCAS), intent(inout)                :: self
  complex(wp), intent(in)                           :: coeff(:,:)
  complex(wp), allocatable, intent(out)             :: e1(:)
  complex(wp), allocatable, intent(out)             :: lgmean(:,:)
  complex(wp), allocatable, intent(out)             :: eh(:)
  complex(wp), allocatable, intent(out)             :: ex(:)
  complex(wp), allocatable, intent(out)             :: eself(:)
  logical, optional, intent(in)                     :: lsamplex
  logical, optional, intent(in)                     :: lself
  type(CEnergy), allocatable, optional, intent(out) :: en(:)
  !local variables
  integer                                           :: ng, nw, w, d
  logical                                           :: lsamplex_, lself_
  complex(wp), allocatable                          :: coeff_bi(:,:)
  complex(wp), allocatable                          :: ovlp(:), den(:), lgmean_(:,:), e1_(:), ex_(:), eself_(:)
  character(len=*), parameter                       :: proc_name = "local_energy_cas"
  
  if (profile_code) call start_profiling(proc_name)

  lsamplex_ = .true.
  if (present(lsamplex)) lsamplex_ = lsamplex

  lself_ = .true.
  if (present(lself)) lself_ = lself

  ng = size(self%hcas%lgmean_frozen)
  nw = size(coeff, 2) / self%ham%des%nocc

  allocate(e1(nw))
  allocate(lgmean(ng, nw))
  allocate(eh(nw))
  allocate(ex(nw))
  allocate(eself(nw))
  allocate(den(nw))

  e1 = zeroc
  lgmean = zeroc
  eh = zeroc
  ex = zeroc
  eself = zeroc
  den = zeroc

  do d = 1, size(self%det)
    call self%get_ispin(self%det(d))
    call self%setup_coeff_occ(self%hdes_cp, self%det(d))

    call get_ovlp_batch(self%hdes_cp, self%coeff_occ, coeff, ovlp)
    ovlp = self%ci_coeff(d) * ovlp
    den = den + ovlp

    call biorth_orb_batch(self%hdes_cp, self%coeff_occ, coeff, coeff_bi)

    call self%setup_h1phi_occ(self%hdes_cp, self%det(d))
    call h1_mean_no_batch(self%hdes_cp, self%h1phi_occ, coeff_bi, e1_)
    e1 = e1 + ovlp * e1_

    call self%setup_h2phi_h_occ(self%hdes_cp, self%det(d))
    call h2_gres_mean_no_batch(self%hdes_cp, self%h2phi_h_occ, coeff_bi, lgmean_, self%des%h_block_size)
    do w = 1, nw
      lgmean(:,w) = lgmean(:,w) + ovlp(w) * lgmean_(:,w)
      eh(w) = eh(w) + 0.5_wp * ovlp(w) * sum(lgmean_(:,w)**2)
    end do 

    if (lself_) then
      call self%setup_h1phi_s_occ(self%hdes_cp, self%det(d))
      call h1_mean_no_batch(self%hdes_cp, self%h1phi_s_occ, coeff_bi, eself_)
      eself = eself + ovlp * eself_
    end if

    if (lsamplex_) then 
      call self%setup_h2phi_x_occ(self%hdes_cp, self%det(d))
      call self%exchange_calculator(self%hdes_cp, coeff_bi, ex_)
      ex = ex + ovlp * ex_
    end if
  end do
  
  e1 = e1 / den
  eself = eself / den
  eh = eh / den 
  ex = ex / den
  
  do w = 1, nw
    lgmean(:,w) = lgmean(:,w) / den(w)
  end do

  if (present(en)) then
    allocate(en(nw))
    do w = 1, nw
      en(w) = self%pack_energy(e1(w), eh(w), ex(w))
    end do
  end if

  if (profile_code) call end_profiling(proc_name)
end subroutine local_energy_batch_cas


!******************************************************************************** 
!
! Calculate g-resolved local energy contributions for a batch of walkers
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
subroutine local_energy_batch_gres_cas(self, coeff, e1, lgmean, eh, ex, eself, lsamplex, lself, en)
  class(WaveTrialCAS), intent(inout)                :: self
  complex(wp), intent(in)                           :: coeff(:,:)
  complex(wp), allocatable, intent(out)             :: e1(:)
  complex(wp), allocatable, intent(out)             :: lgmean(:,:)
  complex(wp), allocatable, intent(out)             :: eh(:,:)
  complex(wp), allocatable, intent(out)             :: ex(:,:)
  complex(wp), allocatable, intent(out)             :: eself(:,:)
  logical, optional, intent(in)                     :: lsamplex
  logical, optional, intent(in)                     :: lself
  type(CEnergy), allocatable, optional, intent(out) :: en(:)
  !local variables
  integer                                           :: w, ng, nw, d
  logical                                           :: lsamplex_, lself_
  complex(wp), allocatable                          :: coeff_bi(:,:)
  complex(wp), allocatable                          :: ovlp(:), den(:), lgmean_(:,:), e1_(:), ex_(:,:), eself_(:,:)
  character(len=*), parameter                       :: proc_name = "local_energy_gres_cas"
  
  if (profile_code) call start_profiling(proc_name)

  lsamplex_ = .true.
  if (present(lsamplex)) lsamplex_ = lsamplex

  lself_ = .true.
  if (present(lself)) lself_ = lself

  ng = size(self%hcas%lgmean_frozen)
  nw = size(coeff, 2) / self%ham%des%nocc

  allocate(e1(nw))
  allocate(lgmean(ng, nw))
  allocate(eh(ng, nw))
  allocate(eself(ng, nw))
  allocate(ex(ng, nw))
  allocate(den(nw))

  e1 = zeroc
  lgmean = zeroc
  eself = zeroc
  ex = zeroc
  den = zeroc

  do d = 1, size(self%det)
    call self%get_ispin(self%det(d))
    call self%setup_coeff_occ(self%hdes_cp, self%det(d))

    call get_ovlp_batch(self%hdes_cp, self%coeff_occ, coeff, ovlp)
    ovlp = self%ci_coeff(d) * ovlp
    den = den + ovlp

    call biorth_orb_batch(self%hdes_cp, self%coeff_occ, coeff, coeff_bi)

    call self%setup_h1phi_occ(self%hdes_cp, self%det(d))
    call h1_mean_no_batch(self%hdes_cp, self%h1phi_occ, coeff_bi, e1_)

    call self%setup_h2phi_h_occ(self%hdes_cp, self%det(d))
    call h2_gres_mean_no_batch(self%hdes_cp, self%h2phi_h_occ, coeff_bi, lgmean_, self%des%h_block_size)

    if (lself_) then
      call self%setup_h2phi_s_occ(self%hdes_cp, self%det(d))
      call h2_gres_mean_no_batch(self%hdes_cp, self%h2phi_s_occ, coeff_bi, eself_, self%des%h_block_size)
    end if

    if (lsamplex_) then 
      call self%setup_h2phi_x_occ(self%hdes_cp, self%det(d))
      call self%exchange_calculator_gres(self%hdes_cp, coeff_bi, ex_)
    end if

    e1 = e1 + ovlp * e1_

    do w = 1, nw
      lgmean(:,w) = lgmean(:,w) + ovlp(w) * lgmean_(:,w)
      eh(:,w) = eh(:,w) + 0.5_wp * ovlp(w) * lgmean_(:,w)**2
      if (lsamplex_) ex(:,w) = ex(:,w) + ovlp(w) * ex_(:,w)
      if (lself_) eself(:,w) = eself(:,w) + ovlp(w) * eself_(:,w)
    end do 
  end do
  
  e1 = e1 / den
  
  do w = 1, nw
    lgmean(:,w) = lgmean(:,w) / den(w)
    eh(:,w) = eh(:,w) / den(w)
    ex(:,w) = ex(:,w) / den(w)
    eself(:,w) = eself(:,w) / den(w)
  end do

  if (present(en)) then
    allocate(en(nw))
    do w = 1, nw
      en(w) = self%pack_energy_gres(e1(w), eh(:,w), ex(:,w))
    end do
  end if

  if (profile_code) call end_profiling(proc_name)
end subroutine local_energy_batch_gres_cas

end module wave_trial_cas