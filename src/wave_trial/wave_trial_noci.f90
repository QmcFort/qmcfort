! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module wave_trial_noci

#include "../preproc.inc"
use constants
use file_handle, only: FileHandle
use mpi, only: mpi_communicator, comm_world
use array_numpy_file_io, only: ArrayNumpyFileIO
use profiling
use qmcfort_io
use energy_types
use hamilton_vars, only: hamil_descriptor
use hamilton_layout
use hamilton_operations, only: sum_gres_potential
use hamilton_type
use wave_trial_des, only: WaveTrialDes
use wave_trial, only: WaveTrial
use c_wave_trial_sd, only: CWaveTrialSD, init_c_wave_trial_sd

implicit none

private
public :: WaveTrialNOCI, init_wave_trial_noci

integer, parameter :: MAX_INPLACE_DET = 2000

type, extends(WaveTrial) :: WaveTrialNOCI
  complex(wp), allocatable        :: coeff(:,:,:,:)
  complex(wp), allocatable        :: coeff_occ(:,:)
  complex(wp), allocatable        :: ci_coeff(:)
  type(CWaveTrialSD), allocatable :: det(:)
contains
  procedure :: setup_coeff => setup_coeff_noci
  procedure :: is_inplace_energy_eval => is_inplace_energy_eval_noci
  procedure :: add_det => add_det_to_noci
  procedure :: setup_det_ptr => setup_det_ptr_noci
  procedure :: report_orbital_character => report_noci_orbital_character

  procedure :: write_to_file => write_to_file_noci

  !deferred routines
  procedure :: assign_wave_trial => assign_wave_trial_noci
  procedure :: setup => setup_wave_trial_noci
  procedure :: report => report_wave_trial_noci
  procedure :: sizeof => sizeof_wave_trial_noci

  procedure :: trial_energy_full => trial_energy_noci
  procedure :: trial_h2_gres_mean => trial_h2_gres_mean_noci
  procedure :: initial_energy_full => initial_energy_noci
  procedure :: initial_h2_gres_mean => initial_h2_gres_mean_noci

  procedure :: get_ovlp_batch => get_ovlp_batch_noci
  procedure :: h1_mean_batch => h1_mean_batch_noci
  procedure :: h2_gres_mean_batch => h2_gres_mean_batch_noci
  procedure :: exchange_energy_batch => exchange_energy_batch_noci
  procedure :: exchange_energy_batch_gres => exchange_energy_batch_gres_noci
  procedure :: self_energy_batch => self_energy_batch_noci
  procedure :: self_energy_batch_gres => self_energy_batch_gres_noci
  procedure :: local_energy_batch => local_energy_batch_noci
  procedure :: local_energy_batch_gres => local_energy_batch_gres_noci
end type WaveTrialNOCI

interface WaveTrialNOCI
  module procedure finit_wave_trial_noci
end interface WaveTrialNOCI

contains

!******************************************************************************** 
!
! Initialization of the WaveTrialNOCI object
!
!******************************************************************************** 
subroutine init_wave_trial_noci(wtdes, ham, coeff, ci_coeff, self)
  type(WaveTrialDes), intent(in)       :: wtdes
  type(Hamilton), target, intent(in)   :: ham
  complex(wp), allocatable, intent(in) :: coeff(:,:,:,:)
  complex(wp), allocatable, intent(in) :: ci_coeff(:)
  type(WaveTrialNOCI), intent(out)     :: self

  self%des = wtdes
  allocate(self%coeff, source=coeff)
  allocate(self%ci_coeff, source=ci_coeff)

  call self%base_init(ham)
  call self%setup_coeff(self%ham%des)
end subroutine init_wave_trial_noci


!******************************************************************************** 
!
! WaveTrialNOCI constructor
!
!******************************************************************************** 
function finit_wave_trial_noci(wtdes, ham, coeff, ci_coeff) result(self)
  type(WaveTrialDes), intent(in)          :: wtdes
  type(Hamilton), target, intent(in)      :: ham
  complex(wp), allocatable, intent(inout) :: coeff(:,:,:,:)
  complex(wp), allocatable, intent(inout) :: ci_coeff(:)
  type(WaveTrialNOCI)                     :: self

  call init_wave_trial_noci(wtdes, ham, coeff, ci_coeff, self)
end function finit_wave_trial_noci


!******************************************************************************** 
!
! Assignment operator for the WaveTrialNOCI
!
!******************************************************************************** 
subroutine assign_wave_trial_noci(to, from)
  class(WaveTrialNOCI), intent(inout) :: to
  class(WaveTrial), intent(in)        :: from
  !local
  integer                             :: d

  select type (from)
    type is (WaveTrialNOCI)
      call to%base_assign(from)

      if (allocated(from%coeff)) then
        if (allocated(to%coeff)) deallocate(to%coeff)
        allocate(to%coeff, source=from%coeff)
      end if

      if (allocated(from%coeff_occ)) then 
        if (allocated(to%coeff_occ)) deallocate(to%coeff_occ)
        allocate(to%coeff_occ, source=from%coeff_occ)
      end if

      if (allocated(from%ci_coeff)) then
        if (allocated(to%ci_coeff)) deallocate(to%ci_coeff)
        allocate(to%ci_coeff, source=from%ci_coeff)
      end if

      if (allocated(from%det)) then
        if (allocated(to%det)) deallocate(to%det)
        allocate(to%det(size(from%det)))
        do d = 1, size(to%det)
          to%det(d) = from%det(d)
        end do
      end if
    class default
      if (comm_world%mpirank == 0) write(*,*) "both dynamic types have to be WaveTrialSD in the assignment operator"
      stop
  end select
end subroutine assign_wave_trial_noci


!******************************************************************************** 
!
! Setup WaveTrialNOCI object for the energy evaluation
!
!******************************************************************************** 
subroutine setup_wave_trial_noci(self, lself)
  class(WaveTrialNOCI), intent(inout) :: self
  logical, optional, intent(in)       :: lself
  !local
  integer                             :: d, ndet

  ndet = size(self%ci_coeff)

  if (allocated(self%det)) deallocate(self%det)

  if (self%is_inplace_energy_eval()) then 
    allocate(self%det(ndet))
    do d = 1 , ndet
      self%det(d) = CWaveTrialSD(self%des, self%ham, self%coeff(:,:,:,d))
      call self%det(d)%setup(lself=lself)
    end do
  end if
end subroutine setup_wave_trial_noci


!******************************************************************************** 
!
! Report WaveTrialNOCI object
!
!******************************************************************************** 
subroutine report_wave_trial_noci(self, comm, fh)
  class(WaveTrialNOCI), intent(inout) :: self
  type(mpi_communicator), intent(in)  :: comm
  type(FileHandle), intent(in)        :: fh
  !local
  real(dp)                            :: mem
  type(CEnergy)                       :: e_init, e_trial

  call self%initial_energy(e_init)
  call self%trial_energy(e_trial)

  mem = self%sizeof()

  if (comm%mpirank == 0) then
    write(fh%funit,*)
    write(fh%funit,*)   "Trial wave function descriptor" 
    write(fh%funit,*)   "-------------------------------" 
    write(fh%funit,100) "trial wave function type        ", "WaveTrialNOCI"
    write(fh%funit,100) "size of the object              ", write_bytes(mem)
    write(fh%funit,101) "number of determinants          ", size(self%ci_coeff)
    write(fh%funit,100) "Exchange energy mode            ", self%des%exchange_mode 
    write(fh%funit,101) "block size for Hartree terms    ", self%des%h_block_size
    write(fh%funit,101) "block size for exchange terms   ", self%des%x_block_size
    write(fh%funit,102) "inplce energy evaluation        ", self%is_inplace_energy_eval()
  
!debug: wave_trial : return when estimate_delta_x reimplemented
    !if (self%exchange_mode == "cholesky") then
    !  write(fh%funit,102)  "Compressed Cholesky vectors X   ", self%compress_x
    !  write(fh%funit,103)  "Tolerance for the compression   ", self%compress_x_tol
    !  write(fh%funit,101)  "Number of Cholesky vectors H    ", self%ng
    !  write(fh%funit,101)  "Number of Cholesky vectors X    ", self%ng_comp
    !  write(fh%funit,103)  "Exchange energy correction      ", self%delta_x
    !end if

    call e_trial%report(fh, method="trial")
    call e_init%report(fh, method="init")
  end if

  100 format(1x, t5, a, t50, "= ", a)
  101 format(1x, t5, a, t50, "= ", 10i8)
  102 format(1x, t5, a, t50, "= ", 10l6)
  103 format(1x, t5, a, t50, "= ", 10es12.4)
end subroutine report_wave_trial_noci


!******************************************************************************** 
!
! Report orbital character of the determinants stored in WaveTrialNoci
!
!******************************************************************************** 
subroutine report_noci_orbital_character(self)
  use standalone, only: dump_array
  class(WaveTrialNOCI), intent(in) :: self
  !local 
  integer                          :: spin, ispin, d, nd, nel

  ispin = size(self%coeff, 3)
  nd = size(self%coeff, 4)

  if (comm_world%mpirank==0) then
    do d = 1, nd
      do spin = 1, ispin
        write(*,100) "determinant no. = ", d, " spin component = ", spin
        nel = self%ham%des%nel(spin)
        call dump_array(abs(self%coeff(:,1:nel,spin,d))**2, -1, -1)
      end do
    end do
  end if

  100 format(1x,a,i6,a,i6)
end subroutine report_noci_orbital_character


!******************************************************************************** 
!
! Write WaveTrialNOCI to the file
!
!******************************************************************************** 
subroutine write_to_file_noci(self)
  class(WaveTrialNOCI), intent(in) :: self
  !local
  type(ArrayNumpyFileIO)           :: numpy_io
  character(len=*), parameter      :: coeff_fname = "orbitals_trial_noci"
  character(len=*), parameter      :: ci_coeff_fname = "ci_coeff_noci"
  
  call numpy_io%save(coeff_fname, self%coeff)
  call numpy_io%save(ci_coeff_fname, self%ci_coeff)
end subroutine write_to_file_noci


!******************************************************************************** 
!
! Estimate memory consumption by WaveTrialNOCI object
!
!******************************************************************************** 
function sizeof_wave_trial_noci(self) result(mem)
  class(WaveTrialNOCI), intent(in) :: self
  real(dp)                         :: mem
  !local
  integer                          :: d

  mem = sizeof(self) + sizeof(self%coeff) + sizeof(self%coeff_occ) + sizeof(self%ci_coeff)

  do d = 1, size(self%det)
    mem = mem + self%det(d)%sizeof()
  end do
end function sizeof_wave_trial_noci


!******************************************************************************** 
!
! Obtain occupied orbitals self%coeff_occ from all orbitals stored in self%coeff
!
!******************************************************************************** 
subroutine setup_coeff_noci(self, hdes)
  class(WaveTrialNOCI), intent(inout) :: self
  type(hamil_descriptor), intent(in)  :: hdes
  !local
  integer                             :: spin, sp, ispin
  integer                             :: d, nd, nel, w1, w2

  ispin = size(self%coeff, 3)
  nd = size(self%coeff, 4)

  if (allocated(self%coeff_occ)) deallocate(self%coeff_occ)
  allocate(self%coeff_occ(hdes%n, hdes%nocc*nd))

  do spin = 1, size(self%coeff, 3)
    nel = hdes%nel(spin)
    do d = 1, size(self%coeff, 4)
      w1 = hdes%nell(spin)*nd + (d-1)*nel + 1
      w2 = hdes%nell(spin)*nd + d*nel
      self%coeff_occ(:,w1:w2) = self%coeff(:,1:nel,spin,d)
    end do
  end do
end subroutine setup_coeff_noci


!******************************************************************************** 
!
! Get inplace_energy_eval variable for a specified number of determinants n
!
!******************************************************************************** 
function is_inplace_energy_eval_noci(self, n) result(inplace_energy_eval)
  class(WaveTrialNOCI), intent(in) :: self
  integer, optional, intent(in)    :: n
  logical                          :: inplace_energy_eval
  !local
  integer                          :: n_

  n_ = size(self%ci_coeff)
  if (present(n)) n_ = n

!debug: for the moment use only false
  !inplace_energy_eval = n_ <= MAX_INPLACE_DET
  inplace_energy_eval = .true.
end function is_inplace_energy_eval_noci


!******************************************************************************** 
!
! Add determinant to the WaveTrialnoci obejct, and update CI coefficients
!
!******************************************************************************** 
subroutine add_det_to_noci(self, det, ci_coeff)
  class(WaveTrialNOCI), intent(inout)     :: self
  type(CWaveTrialSD), intent(in)          :: det
  complex(wp), allocatable, intent(inout) :: ci_coeff(:)
  !local
  integer                                 :: d, n
  logical                                 :: inplace_energy_eval
  complex(wp), allocatable                :: coeff(:,:,:,:)
  type(CWaveTrialSD), allocatable         :: det_new(:)

  n = size(ci_coeff)

  call move_alloc(ci_coeff, self%ci_coeff)

  allocate(coeff(size(self%coeff,1), size(self%coeff,2), size(self%coeff,3), n))
  coeff(:,:,:,1:n-1) = self%coeff
  coeff(:,:,:,n) = det%coeff
  call move_alloc(from=coeff, to=self%coeff)

  call self%setup_coeff(self%ham%des)

  inplace_energy_eval = self%is_inplace_energy_eval(n)

  if (inplace_energy_eval) then
    allocate(det_new(n))
    do d = 1, n-1
      det_new(d) = self%det(d)
    end do
    det_new(n) = det

    call move_alloc(from=det_new, to=self%det)
  else 
    if (allocated(self%det)) deallocate(self%det)
  end if
end subroutine add_det_to_noci


!******************************************************************************** 
!
! Setup CWaveTrialSD pointer for the energy evaluation depending on the
! inplace_energy_eval flag
!
!******************************************************************************** 
subroutine setup_det_ptr_noci(self, d, det, lself)
  class(WaveTrialNOCI), target, intent(inout) :: self
  integer, intent(in)                         :: d
  type(CWaveTrialSD), pointer, intent(out)    :: det
  logical, optional, intent(in)               :: lself
  !local
  logical                                     :: lself_
  type(CWaveTrialSD), target, save            :: det_help

  lself_ = .true.
  if (present(lself)) lself_ = lself

  if (self%is_inplace_energy_eval()) then
    det => self%det(d)
  else
    det_help = CWaveTrialSD(self%des, self%ham, self%coeff(:,:,:,d))
    call det_help%setup(lself=lself_)
    det => det_help
  end if
end subroutine setup_det_ptr_noci


!******************************************************************************** 
!
! Estimate trial energy of the WaveTrialNOCI object using 
! prepared CWaveTriaSD obejcts
!
!******************************************************************************** 
subroutine trial_energy_noci(self, e1, lgmean, eh, ex, eself)
  class(WaveTrialNOCI), intent(inout)   :: self
  complex(wp), intent(out)              :: e1
  complex(wp), allocatable, intent(out) :: lgmean(:)
  complex(wp), allocatable, intent(out) :: eh(:)
  complex(wp), allocatable, intent(out) :: ex(:)
  complex(wp), allocatable, intent(out) :: eself(:)
  !local varaibles
  integer                               :: i, j, ng
  complex(wp)                           :: den
  complex(wp), allocatable              :: ovlp(:), e1_(:)
  complex(wp), allocatable              :: lgmean_(:,:), eh_(:,:), ex_(:,:), eself_(:,:)
  type(CWaveTrialSD), pointer           :: det

  ng = self%ham%des%ng

  allocate(lgmean(ng))
  allocate(eh(ng))
  allocate(ex(ng))
  allocate(eself(ng))

  e1 = zeroc
  lgmean = zeroc
  eh = zeroc
  ex = zeroc
  eself = zeroc
  den = zeroc

  do i = 1, size(self%ci_coeff)
    call self%setup_det_ptr(i, det, lself=.true.)

    call det%get_ovlp(self%coeff_occ, ovlp)
    ovlp = conjg(self%ci_coeff(i)) * self%ci_coeff * ovlp
    den = den + sum(ovlp)

    call det%local_energy_gres(self%coeff_occ, e1_, lgmean_, eh_, ex_, eself_)

    e1 = e1 + sum(ovlp * e1_)
    do j = 1, size(self%ci_coeff)
      lgmean = lgmean + ovlp(j) * lgmean_(:,j)
      eh = eh + ovlp(j) * eh_(:,j)
      ex = ex + ovlp(j) * ex_(:,j)
      eself = eself + ovlp(j) * eself_(:,j)
    end do
  end do

  e1 = e1 / den
  lgmean = lgmean / den
  eh = eh / den
  ex = ex / den
  eself = eself / den
end subroutine trial_energy_noci


!******************************************************************************** 
!
! Trial exectation vlaue of the Cholesky vectors for WaveTrialNOCI
!
!    <L_g> = <Psi_T|L_g|Psi_T> / <Psi_T|Psi_T>
!          = \sum_{ij} <Psi_i|L_g|Psi_j> / \sum_{ij} <Psi_i|Psi_j>
!
!******************************************************************************** 
subroutine trial_h2_gres_mean_noci(self, lgmean)
  class(WaveTrialNOCI), intent(inout)   :: self
  complex(wp), allocatable, intent(out) :: lgmean(:)
  !local
  integer                               :: i, j, ng
  complex(wp)                           :: den
  complex(wp), allocatable              :: ovlp(:), lgmean_(:,:), eh_(:,:)
  type(CWaveTrialSD), pointer           :: det

  ng = self%ham%des%ng
  allocate(lgmean(ng))
  lgmean = zeroc
  den = zeroc
  
  do i = 1, size(self%ci_coeff)
    call self%setup_det_ptr(i, det, lself=.false.)

    call det%get_ovlp(self%coeff_occ, ovlp)
    ovlp = conjg(self%ci_coeff(i)) * self%ci_coeff * ovlp
    den = den + sum(ovlp)

    call det%h2_gres_mean(self%coeff_occ, lgmean_, eh_)

    do j = 1, size(self%ci_coeff)
      lgmean = lgmean + ovlp(j) * lgmean_(:,j)
    end do
  end do

  lgmean = lgmean / den
end subroutine trial_h2_gres_mean_noci


!******************************************************************************** 
!
! Initial energy of the WaveTrialNOCI object using WaveTrialSD objects and
! asumming first determinant is the initial determinant
!
! Note:
!    Don't forget that the initial determinant is the first determinant in the 
!    WaveTrialNOCI object
!
!******************************************************************************** 
subroutine initial_energy_noci(self, e1, lgmean, eh, ex, eself)
  class(WaveTrialNOCI), intent(inout)   :: self
  complex(wp), intent(out)              :: e1
  complex(wp), allocatable, intent(out) :: lgmean(:)
  complex(wp), allocatable, intent(out) :: eh(:)
  complex(wp), allocatable, intent(out) :: ex(:)
  complex(wp), allocatable, intent(out) :: eself(:)
  !local varaibles
  integer                               :: d, ng
  complex(wp)                           :: e1_, den, ovlp
  complex(wp), allocatable              :: lgmean_(:), eh_(:), ex_(:), eself_(:)
  type(CWaveTrialSD), pointer           :: det, det1

  ng = self%ham%des%ng 

  allocate(lgmean(ng))
  allocate(eh(ng))
  allocate(ex(ng))
  allocate(eself(ng))

  e1 = zeroc
  lgmean = zeroc
  eh = zeroc
  ex = zeroc
  eself = zeroc
  den = zeroc

  call self%setup_det_ptr(1, det1, lself=.false.)

  do d = 1, size(self%ci_coeff)
    call self%setup_det_ptr(d, det, lself=.true.)

    call det%get_ovlp(det1%coeff_occ, ovlp)
    ovlp = conjg(self%ci_coeff(d)) * ovlp
    den = den + ovlp

    call det%local_energy_gres(det1%coeff_occ, e1_, lgmean_, eh_, ex_, eself_)

    e1 = e1 + ovlp * e1_
    lgmean = lgmean + ovlp * lgmean_
    eh = eh + ovlp * eh_
    ex = ex + ovlp * ex_
    eself = eself + ovlp * eself_
  end do

  e1 = e1 / den
  lgmean = lgmean / den
  eh = eh / den
  ex = ex / den
  eself = eself / den
end subroutine initial_energy_noci


!******************************************************************************** 
!
! Initial exectation vlaue of the Cholesky vectors for WaveTrialNOCI
!
!    <L_g> = <Psi_T|L_g|Psi_I> / <Psi_T|Psi_I>
!          = \sum_i c_i* <Psi_i|L_g|Psi_1> / \sum_i c_i* <Psi_i|Psi_1>
!
!******************************************************************************** 
subroutine initial_h2_gres_mean_noci(self, lgmean)
  class(WaveTrialNOCI), intent(inout)   :: self
  complex(wp), allocatable, intent(out) :: lgmean(:)
  !local
  integer                               :: d, ng
  complex(wp)                           :: den, ovlp
  complex(wp), allocatable              :: lgmean_(:), eh_(:)
  type(CWaveTrialSD), pointer           :: det, det1

  ng = self%ham%des%ng
  allocate(lgmean(ng))
  lgmean = zeroc
  den = zeroc

  call self%setup_det_ptr(1, det1, lself=.false.)
  !debug: optimization : this loop should be removed by calculating everything 
  !with respect to the det1
  do d = 1, size(self%ci_coeff)
    call self%setup_det_ptr(d, det, lself=.false.)

    call det%get_ovlp(det1%coeff_occ, ovlp)
    ovlp = conjg(self%ci_coeff(d)) * ovlp
    den = den + ovlp

    call det%h2_gres_mean(det1%coeff_occ, lgmean_, eh_)
    lgmean = lgmean + ovlp * lgmean_
  end do

  lgmean = lgmean / den
end subroutine initial_h2_gres_mean_noci


!******************************************************************************** 
!
! Calculate overlap between WaveTrialNOCI and walker determinants
!
!    <Psi_T|Phi> = \sum_i c_i* <Psi_i|Phi>
!
! where <Psi_i|Phi> is the overlap between two non-orthogonal determinants
!
!******************************************************************************** 
subroutine get_ovlp_batch_noci(self, coeff, ovlp)
  class(WaveTrialNOCI), intent(inout)   :: self
  complex(wp), intent(in)               :: coeff(:,:)
  complex(wp), allocatable, intent(out) :: ovlp(:)
  !local
  integer                               :: d, nw
  complex(wp), allocatable              :: ovlp_(:)
  type(CWaveTrialSD)                    :: det
  character(len=*), parameter           :: proc_name = "get_ovlp_noci"

  if (profile_code) call start_profiling(proc_name)

  nw = size(coeff, 2) / self%ham%des%nocc

  allocate(ovlp(nw))
  ovlp = zeroc

  do d = 1, size(self%ci_coeff)
    det = CWaveTrialSD(self%des, self%ham, self%coeff(:,:,:,d))
    call det%get_ovlp(coeff, ovlp_)
    ovlp = ovlp + conjg(self%ci_coeff(d)) * ovlp_
  end do 

  if (profile_code) call end_profiling(proc_name)
end subroutine get_ovlp_batch_noci


!******************************************************************************** 
!
! Calculate one-electron energy for a batch of walkers 
!
!    E1 = <Psi_T | H1 | Phi_w> / <Psi_T | Phi_w>
!       = \sum_d c_d* <Phi_d | H1 | Phi_w>  / \sum_d c_d* <Phi_d | Phi_w>
!
!******************************************************************************** 
subroutine h1_mean_batch_noci(self, coeff, e1)
  class(WaveTrialNOCI), intent(inout)   :: self
  complex(wp), intent(in)               :: coeff(:,:)
  complex(wp), allocatable, intent(out) :: e1(:)
  !local
  integer                               :: nw, d
  complex(wp), allocatable              :: ovlp(:), e1_(:), den(:)
  type(CWaveTrialSD), pointer           :: det
  character(len=*), parameter           :: proc_name = "h1_mean_noci"

  if (profile_code) call start_profiling(proc_name)

  nw = size(coeff, 2) / self%ham%des%nocc

  allocate(e1(nw))
  allocate(den(nw))

  e1 = zeroc
  den = zeroc

  do d = 1, size(self%ci_coeff)
    call self%setup_det_ptr(d, det, lself=.false.)

    call det%get_ovlp(coeff, ovlp)
    ovlp = conjg(self%ci_coeff(d)) * ovlp
    den = den + ovlp

    call det%h1_mean(coeff, e1_)
    e1 = e1 + ovlp * e1_
  end do

  e1 = e1 / den

  if (profile_code) call end_profiling(proc_name)
end subroutine h1_mean_batch_noci


!******************************************************************************** 
!
! Calculate expectation value of the Cholesky vectors for a batch of walkers
!
!    l_{g,w} = <Psi_T |L_g |Phi_w> / <Psi_T | Phi_w>
!            = \sum_d c_d* <Phi_d | L_g | Phi_w> / \sum_d c_d*
!
! <Phi_d | Phi_w>    and
! <Phi_d | L_g | Phi_w> / <Phi_d | Phi_w>
! computed using CWaveTrialSD
!
!******************************************************************************** 
subroutine h2_gres_mean_batch_noci(self, coeff, lgmean, eh)
  class(WaveTrialNOCI), intent(inout)   :: self
  complex(wp), intent(in)               :: coeff(:,:)
  complex(wp), allocatable, intent(out) :: lgmean(:,:)
  complex(wp), allocatable, intent(out) :: eh(:,:)
  !local
  integer                               :: ng, nw, w, d
  complex(wp), allocatable              :: ovlp(:), den(:), lgmean_(:,:), eh_(:,:)
  type(CWaveTrialSD), pointer           :: det
  character(len=*), parameter           :: proc_name = "h2_gres_mean_noci"

  if (profile_code) call start_profiling(proc_name)

  ng = self%ham%des%ng
  nw = size(coeff, 2) / self%ham%des%nocc

  allocate(lgmean(ng,nw))
  allocate(eh(ng,nw))
  allocate(den(nw))

  lgmean = zeroc
  eh = zeroc
  den = zeroc

  do d = 1, size(self%ci_coeff)
    call self%setup_det_ptr(d, det, lself=.false.)

    call det%get_ovlp(coeff, ovlp)
    ovlp = conjg(self%ci_coeff(d)) * ovlp
    den = den + ovlp

    call det%h2_gres_mean(coeff, lgmean_, eh_)

    do w = 1, nw
      lgmean(:,w) = lgmean(:,w) + ovlp(w) * lgmean_(:,w)
      eh(:,w) = eh(:,w) + ovlp(w) * eh_(:,w)
    end do
  end do

  do w = 1, nw
    lgmean(:,w) = lgmean(:,w) / den(w)
    eh(:,w) = eh(:,w) / den(w)
  end do

  if (profile_code) call end_profiling(proc_name)
end subroutine h2_gres_mean_batch_noci


!******************************************************************************** 
!
! Calculate exchange energy for the batch of walkers
!
!    Ex_w = \sum_d c_d* <Psi_d | Hx | Phi_w> / \sum_d c_d* <Psi_d | Phi_w>
! 
! <Phi_d | Phi_w>    and
! <Phi_d | Hx | Phi_w> / <Phi_d | Phi_w>
! computed using CWaveTrialSD
!
!******************************************************************************** 
subroutine exchange_energy_batch_noci(self, coeff, ex)
  class(WaveTrialNOCI), intent(inout)   :: self
  complex(wp), intent(in)               :: coeff(:,:)
  complex(wp), allocatable, intent(out) :: ex(:)
  !local
  integer                               :: nw, d
  complex(wp), allocatable              :: den(:), ovlp(:), ex_(:)
  type(CWaveTrialSD), pointer           :: det
  character(len=*), parameter           :: proc_name = "exchange_noci"
  
  if (profile_code) call start_profiling(proc_name)

  nw = size(coeff,2) / self%ham%des%nocc

  allocate(ex(nw))
  allocate(den(nw))

  ex = zeroc
  den = zeroc

  do d = 1, size(self%ci_coeff)
    call self%setup_det_ptr(d, det, lself=.false.)

    call det%get_ovlp(coeff, ovlp)
    ovlp = conjg(self%ci_coeff(d)) * ovlp
    den = den + ovlp

    call det%exchange_energy(coeff, ex_)
    ex = ex + ovlp * ex_
  end do

  ex = ex / den

  if (profile_code) call end_profiling(proc_name)
end subroutine exchange_energy_batch_noci


!******************************************************************************** 
!
! Calculate exchange energy for the batch of walkers
!
!    Ex_{g,w} = \sum_d c_d <Psi_d | Hx_g | Phi_w> / \sum_d c_d <Psi_d | Phi_w>
! 
! <Phi_d | Phi_w>    and
! <Phi_d | Hx_g | Phi_w> / <Phi_d | Phi_w>
! computed using CWaveTrialSD
!
!******************************************************************************** 
subroutine exchange_energy_batch_gres_noci(self, coeff, ex)
  class(WaveTrialNOCI), intent(inout)   :: self
  complex(wp), intent(in)               :: coeff(:,:)
  complex(wp), allocatable, intent(out) :: ex(:,:)
  !local
  integer                               :: ng, nw, w, d
  complex(wp), allocatable              :: den(:), ovlp(:), ex_(:,:)
  type(CWaveTrialSD), pointer           :: det
  character(len=*), parameter           :: proc_name = "exchange_noci"
  
  if (profile_code) call start_profiling(proc_name)

  ng = self%ham%des%ng
  nw = size(coeff,2) / self%ham%des%nocc

  allocate(ex(ng, nw))
  allocate(den(nw))

  ex = zeroc
  den = zeroc

  do d = 1, size(self%ci_coeff)
    call self%setup_det_ptr(d, det, lself=.false.)

    call det%get_ovlp(coeff, ovlp)
    ovlp = conjg(self%ci_coeff(d)) * ovlp
    den = den + ovlp

    call det%exchange_energy_gres(coeff, ex_)

    do w = 1, nw
      ex(:,w) = ex(:,w) + ovlp(w) * ex_(:,w)
    end do
  end do

  do w = 1, nw
    ex(:,w) = ex(:,w) / den(w)
  end do

  if (profile_code) call end_profiling(proc_name)
end subroutine exchange_energy_batch_gres_noci


!******************************************************************************** 
!
! Calculate self energy for a batch of walkers
!
!    Eself_{w} = \sum_g <Psi_T | L_g^2 | Phi_w> / <Psi_T | Phi_w>
!              = \sum_g \sum_d c_d* <Phi_d|L_g^2| Phi_w> / c_d*
!
! <Phi_d | Phi_w>    and
! <Phi_d | L_g^2 | Phi_w> / <Phi_d | Phi_w>
! computed using CWaveTrialSD
!
!******************************************************************************** 
subroutine self_energy_batch_noci(self, coeff, eself)
  class(WaveTrialNOCI), intent(inout)   :: self
  complex(wp), intent(in)               :: coeff(:,:)
  complex(wp), allocatable, intent(out) :: eself(:)
  !local
  integer                               :: nw, d
  complex(wp), allocatable              :: ovlp(:), den(:), eself_(:)
  type(CWaveTrialSD), pointer           :: det
  character(len=*), parameter           :: proc_name = "self_energy_noci"

  if (profile_code) call start_profiling(proc_name)

  nw = size(coeff, 2) / self%ham%des%nocc

  allocate(eself(nw))
  allocate(den(nw))

  eself = zeroc
  den = zeroc

  do d = 1, size(self%ci_coeff)
    call self%setup_det_ptr(d, det, lself=.true.)

    call det%get_ovlp(coeff, ovlp)
    ovlp = conjg(self%ci_coeff(d)) * ovlp
    den = den + ovlp

    call det%self_energy(coeff, eself_)
    eself = eself + ovlp * eself_
  end do

  eself = eself / den

  if (profile_code) call end_profiling(proc_name)
end subroutine self_energy_batch_noci


!******************************************************************************** 
!
! Calculate g-resolved self energy for a batch of walkers
!
!    Eself_{g,w} = <Psi_T | L_g^2 | Phi_w> / <Psi_T | Phi_w>
!                = \sum_d c_d* <Phi_d|L_g^2| Phi_w> / c_d*
!
!******************************************************************************** 
subroutine self_energy_batch_gres_noci(self, coeff, eself)
  class(WaveTrialNOCI), intent(inout)   :: self
  complex(wp), intent(in)               :: coeff(:,:)
  complex(wp), allocatable, intent(out) :: eself(:,:)
  !local
  integer                               :: ng, nw, w, d
  complex(wp), allocatable              :: ovlp(:), den(:), eself_(:,:)
  type(CWaveTrialSD), pointer           :: det
  character(len=*), parameter           :: proc_name = "self_energy_gres_noci"

  if (profile_code) call start_profiling(proc_name)

  ng = self%ham%des%ng 
  nw = size(coeff, 2) / self%ham%des%nocc

  allocate(eself(ng,nw))
  allocate(den(nw))

  eself = zeroc
  den = zeroc

  do d = 1, size(self%ci_coeff)
    call self%setup_det_ptr(d, det, lself=.true.)

    call det%get_ovlp(coeff, ovlp)
    ovlp = conjg(self%ci_coeff(d)) * ovlp
    den = den + ovlp

    call det%self_energy_gres(coeff, eself_)

    do w = 1, nw
      eself(:,w) = eself(:,w) + ovlp(w) * eself_(:,w)
    end do
  end do

  do w = 1, nw
    eself(:,w) = eself(:,w) / den(w)
  end do

  if (profile_code) call end_profiling(proc_name)
end subroutine self_energy_batch_gres_noci


!******************************************************************************** 
!
! Calculate local energy contributions for a batch of walkers
!
!    E1 = <Psi_T | H1 | Phi>
!    l_g = <Psi_T | L_g | Phi>
!    Ex = <Psi_T | Hx | Phi>
!    Eself = <Psi_T | Hself | Phi>
!
!
!******************************************************************************** 
subroutine local_energy_batch_noci(self, coeff, e1, lgmean, eh, ex, eself, lsamplex, lself, en)
  class(WaveTrialNOCI), intent(inout)               :: self
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
  integer                                           :: ng, nw, w, d
  complex(wp), allocatable                          :: ovlp(:), den(:), e1_(:), lgmean_(:,:), eh_(:), ex_(:), eself_(:)
  type(CWaveTrialSD), pointer                       :: det
  character(len=*), parameter                       :: proc_name = "local_energy_noci"
  
  if (profile_code) call start_profiling(proc_name)

  ng = self%ham%des%ng 
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

  do d = 1, size(self%ci_coeff)
    call self%setup_det_ptr(d, det, lself=lself)

    call det%get_ovlp(coeff, ovlp)
    ovlp = conjg(self%ci_coeff(d)) * ovlp
    den = den + ovlp

    call det%local_energy(coeff, e1_, lgmean_, eh_, ex_, eself_, lsamplex, lself)

    e1 = e1 + ovlp * e1_
    eh = eh + ovlp * eh_
    ex = ex + ovlp * ex_
    eself = eself + ovlp * eself_

    do w = 1, nw
      lgmean(:,w) = lgmean(:,w) + ovlp(w) * lgmean_(:,w)
    end do 
  end do
  
  e1 = e1 / den
  eh = eh / den 
  ex = ex / den
  eself = eself / den
  
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
end subroutine local_energy_batch_noci


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
subroutine local_energy_batch_gres_noci(self, coeff, e1, lgmean, eh, ex, eself, lsamplex, lself, en)
  class(WaveTrialNOCI), intent(inout)               :: self
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
  integer                                           :: w, ng, nw, d
  complex(wp), allocatable                          :: ovlp(:), den(:), e1_(:), lgmean_(:,:), eh_(:,:), ex_(:,:), eself_(:,:)
  type(CWaveTrialSD), pointer                       :: det
  character(len=*), parameter                       :: proc_name = "local_energy_gres_noci"
  
  if (profile_code) call start_profiling(proc_name)

  ng = self%ham%des%ng 
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

  do d = 1, size(self%ci_coeff)
    call self%setup_det_ptr(d, det, lself=lself)

    call det%get_ovlp(coeff, ovlp)
    ovlp = conjg(self%ci_coeff(d)) * ovlp
    den = den + ovlp

    call det%local_energy_gres(coeff, e1_, lgmean_, eh_, ex_, eself_, lsamplex, lself)

    e1 = e1 + ovlp * e1_

    do w = 1, nw
      lgmean(:,w) = lgmean(:,w) + ovlp(w) * lgmean_(:,w)
      eh(:,w) = eh(:,w) + ovlp(w) * eh_(:,w)
      ex(:,w) = ex(:,w) + ovlp(w) * ex_(:,w)
      eself(:,w) = eself(:,w) + ovlp(w) * eself_(:,w)
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
end subroutine local_energy_batch_gres_noci

end module wave_trial_noci