! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module noci_selector

#include "../preproc.inc"
use constants
use profiling
use file_handle, only: FileHandle
use lapack
use energy_types, only: CEnergy
use c_wave_trial_sd, only: CWaveTrialSD, init_c_wave_trial_sd
use wave_trial_noci, only: WaveTrialNOCI
use noci_selector_des, only: NOCISelectorDes

implicit none

private
public :: NOCISelector, init_noci_selector

type, private :: NOCIMAtrixElements
  complex(wp)              :: diag_ovlp    !<mu|mu>
  complex(wp)              :: diag_ham     !<mu|H|mu>
  complex(wp), allocatable :: col_ovlp(:)  !<p|mu>
  complex(wp), allocatable :: col_ham(:)   !<p|H|mu>
end type NOCIMAtrixElements

type NOCISelector
  type(NOCISelectorDes)    :: des
  type(WaveTrialNOCI)      :: phi

  real(wp), allocatable    :: old_ci_energy(:)
  real(wp), allocatable    :: ci_energy(:)

  complex(wp), allocatable :: ovlp(:,:)
  complex(wp), allocatable :: inv_ovlp(:,:)
  complex(wp), allocatable :: ham(:,:)

  real(wp)                 :: ovlp_residual
  real(wp)                 :: energy_residual
contains
  generic            :: assignment(=) => assign
  procedure, private :: assign => assign_noci
  procedure          :: report => report_noci_selector
  procedure          :: report_update => report_noci_update
  procedure          :: init_ovlp_matrix => init_ovlp_matrix_noci
  procedure          :: init_ham_matrix => init_ham_matrix_noci
  procedure          :: noci_matrix_elements => get_noci_matrix_elements
  procedure          :: get_weights => get_noci_weights
  procedure          :: get_ci_weights => get_noci_ci_weights
  procedure          :: compress => compress_noci_selector
  procedure          :: update => update_noci_selector
  procedure          :: expand => expand_noci_selector
  procedure          :: overlap_test => overlap_test_noci
  procedure          :: energy_test => energy_test_noci
end type NOCISelector

interface NOCISelector
  module procedure finit_noci_selector
end interface NOCISelector

contains

!******************************************************************************** 
!
! Initialization of the NOCIMatrixElements object
!
! For a new determinant |mu>, and a current |NOCI> = \sum_p c_p |p>, calcualtes:
!    <mu|mu> = s,     <mu|H|mu> = e,
!    <p|mu> = s_p,    <p|H|mu> = e_p
!
!******************************************************************************** 
subroutine get_noci_matrix_elements(self, det, ci_el)
  class(NOCISelector), intent(in)         :: self
  type(CWaveTrialSD), intent(inout)       :: det
  type(NOCIMAtrixElements), intent(inout) :: ci_el
  !local
  type(CEnergy)                           :: trial_en
  type(CEnergy), allocatable              :: ens(:)

  call det%get_ovlp(det%coeff_occ, ci_el%diag_ovlp)

  call det%local_energy(det%coeff_occ, trial_en, lself=.false.)
  ci_el%diag_ham = trial_en%total_energy() * ci_el%diag_ovlp

  call det%get_ovlp(self%phi%coeff_occ, ci_el%col_ovlp)
  ci_el%col_ovlp = conjg(ci_el%col_ovlp)

  call det%local_energy(self%phi%coeff_occ, ens, lself=.false.)
  ci_el%col_ham = conjg(ens%total_energy()) * ci_el%col_ovlp
end subroutine get_noci_matrix_elements


!******************************************************************************** 
!
! Initialization of the NOCISelector object
!
!******************************************************************************** 
subroutine init_noci_selector(des, phi_t_noci, self)
  type(NOCISelectorDes), intent(in) :: des
  type(WaveTrialNOCI), intent(in)   :: phi_t_noci
  type(NOCISelector), intent(out)   :: self
  !local
  complex(wp), allocatable          :: ci_coeff(:)

  self%des = des
  self%phi = phi_t_noci

  call self%init_ovlp_matrix()
  call self%init_ham_matrix()
  call diagonalize_noci(self%ovlp, self%ham, self%phi%ci_coeff, self%ci_energy)
end subroutine init_noci_selector


!******************************************************************************** 
!
! NOCISelector constructor
!
!******************************************************************************** 
function finit_noci_selector(des, phi_t_noci) result(self)
  type(NOCISelectorDes), intent(in) :: des
  type(WaveTrialNOCI), intent(in)   :: phi_t_noci
  type(NOCISelector)                :: self

  call init_noci_selector(des, phi_t_noci, self)
end function finit_noci_selector


!******************************************************************************** 
!
! Assignment for the NOCISelectro
!
!******************************************************************************** 
subroutine assign_noci(to, from)
  class(NOCISelector), intent(inout) :: to
  type(NOCISelector), intent(in)     :: from

  to%des = from%des
  to%phi = from%phi

  to%ci_energy = from%ci_energy
  allocate(to%ovlp, source=from%ovlp)
  allocate(to%inv_ovlp, source=from%inv_ovlp)
  allocate(to%ham, source=from%ham)
end subroutine assign_noci


!******************************************************************************** 
!
! Report NOCISelector type
!
!******************************************************************************** 
subroutine report_noci_selector(self, fh)
  class(NOCISelector), intent(in) :: self
  type(FileHandle), intent(in)    :: fh
  !local
  integer                         :: w
  real(wp), allocatable           :: weights(:), ci_weights(:)
  complex(wp)                     :: det
  complex(wp), allocatable        :: ovlp(:,:)

  if (comm_world%mpirank == 0) then
    weights = self%get_weights()
    ci_weights =  self%get_ci_weights()

    write(fh%funit,*)
    write(fh%funit,*) "Final number of NOCI determinants ", size(weights)
    write(fh%funit,100) "    det no.    ", "     |c_i|^2    ", "     weight    ", "      H_ii     ",  "      E_i      " 
    write(fh%funit,100) "    =======    ", "     =======    ", "     ======    ", "      ====     ",  "      ===      "
    do w = 1, size(weights)
      write(fh%funit,101) " noci:f: ", w, ci_weights(w), weights(w), real(self%ham(w,w), wp), self%ci_energy(w)
    end do
    write(fh%funit,*)

    allocate(ovlp, source=self%ovlp)
    call determinant(ovlp, det)
    write(fh%funit,*) "determinant of the overlap matrix is = ", det
  end if

  100 format (1x,t10,a,t25,a,t40,a,t55,a,t70,a)
  101 format (1x,a,t10,i6,t30,f6.4,t45,f6.4,t57,f10.4,t72,f10.4)
end subroutine report_noci_selector


!******************************************************************************** 
!
! Report NOCISelector update
!
!******************************************************************************** 
subroutine report_noci_update(self)
  class(NOCISelector), intent(in) :: self
  !local
  integer                         :: n

  if (comm_world%mpirank == 0) then
    n = size(self%ci_energy)
    write(io%qmcfort_out%funit, 100) "noci:s:", n, self%old_ci_energy(1), self%ci_energy(1), self%energy_residual, self%ovlp_residual
    write(io%qmcfort_log%funit, 100) "noci:s:", n, self%old_ci_energy(1), self%ci_energy(1), self%energy_residual, self%ovlp_residual
    write(io%screen%funit, 100)      "noci:s:", n, self%old_ci_energy(1), self%ci_energy(1), self%energy_residual, self%ovlp_residual
  end if

  100 format (1x,t10,a,t20,i6,t35,f12.6,t54,f12.6,t72,es12.4,t87,es12.4)
end subroutine report_noci_update


!******************************************************************************** 
!
! Initialize overlap matrix for the determiants stored in self%phi
!
!******************************************************************************** 
subroutine init_ovlp_matrix_noci(self)
  class(NOCISelector), intent(inout) :: self
  !local
  integer                            :: d, n
  complex(wp), allocatable           :: col_ovlp(:)
  type(CWaveTrialSD)                 :: det

  n = size(self%phi%ci_coeff)

  allocate(self%ovlp(n,n))

  do d = 1, n
    call init_c_wave_trial_sd(self%phi%des, self%phi%ham, self%phi%coeff(:,:,:,d), det)
    call det%get_ovlp(self%phi%coeff_occ, col_ovlp)
    self%ovlp(:,d) = conjg(col_ovlp)
  end do

  allocate(self%inv_ovlp(n,n))
  self%inv_ovlp = self%ovlp
  call inverse(self%inv_ovlp)
end subroutine init_ovlp_matrix_noci


!******************************************************************************** 
!
! Initialize Hamilton matrix for the determiants stored in self%phi
!
!******************************************************************************** 
subroutine init_ham_matrix_noci(self)
  class(NOCISelector), intent(inout) :: self
  !local
  integer                            :: d, n
  type(CEnergy), allocatable         :: col_ham(:)
  type(CWaveTrialSD), pointer        :: det

  n = size(self%phi%ci_coeff)

  allocate(self%ham(n,n))

  do d = 1, n
    call self%phi%setup_det_ptr(d, det, lself=.false.)
    call det%local_energy(self%phi%coeff_occ, col_ham, lself=.false.)
    self%ham(:,d) = conjg(col_ham%total_energy())
  end do

  self%ham = self%ham * self%ovlp
end subroutine init_ham_matrix_noci


!******************************************************************************** 
!
! Main driver routine for NOCISelector:
!
!    perform overlap and energy test, and if both passed add this determinant
!
!******************************************************************************** 
subroutine update_noci_selector(self, coeff)
  class(NOCISelector), intent(inout) :: self
  complex(wp), intent(in)            :: coeff(:,:,:)
  !local
  logical                            :: lovlp, lenergy
  type(CWaveTrialSD)                 :: det
  type(NOCIMAtrixElements)           :: ci_el

  if (use_profiling) call start_profiling("init_det")
  call init_c_wave_trial_sd(self%phi%des, self%phi%ham, coeff, det)
  if (use_profiling) call end_profiling("init_det")

  if (use_profiling) call start_profiling("setup_det")
  call det%setup(lself=.false.)
  if (use_profiling) call end_profiling("setup_det")

  if (use_profiling) call start_profiling("noci_elements")
  call self%noci_matrix_elements(det, ci_el)
  if (use_profiling) call end_profiling("noci_elements")

  lovlp = self%overlap_test(ci_el)
  if (.not. lovlp) return

  lenergy = self%energy_test(ci_el)

  if (lovlp .and. lenergy) call self%expand(det, ci_el)
  if (lovlp .and. lenergy) call self%report_update()
end subroutine update_noci_selector


!******************************************************************************** 
!
! Initialize Hamilton matrix for the determiants stored in self%phi
!
!******************************************************************************** 
subroutine expand_noci_selector(self, det, ci_el)
  class(NOCISelector), intent(inout)   :: self
  type(CWaveTrialSD), intent(in)       :: det
  type(NOCIMAtrixElements), intent(in) :: ci_el
  !local
  integer                              :: n
  complex(wp), allocatable             :: ci_coeff(:)
  complex(wp), allocatable             :: ovlp(:,:), inv_ovlp(:,:), ham(:,:)
  character(len=*), parameter          :: proc_name="expand_noci"

  if (use_profiling) call start_profiling(proc_name)

  n = size(self%ovlp, 1)

  allocate(ovlp(n+1,n+1))
  allocate(inv_ovlp(n+1,n+1))
  allocate(ham(n+1,n+1))

  ovlp(1:n,1:n) = self%ovlp
  ovlp(1:n,n+1) = ci_el%col_ovlp
  ovlp(n+1,1:n) = conjg(ci_el%col_ovlp)
  ovlp(n+1,n+1) = ci_el%diag_ovlp

  inv_ovlp = ovlp
  call inverse(inv_ovlp)

  ham(1:n,1:n) = self%ham
  ham(1:n,n+1) = ci_el%col_ham
  ham(n+1,1:n) = conjg(ci_el%col_ham)
  ham(n+1,n+1) = ci_el%diag_ham

  call move_alloc(self%ci_energy, self%old_ci_energy)

  if (use_profiling) call start_profiling("diag_noci")
  call diagonalize_noci(ovlp, ham, ci_coeff, self%ci_energy)
  if (use_profiling) call end_profiling("diag_noci")

  call move_alloc(ovlp, self%ovlp)
  call move_alloc(inv_ovlp, self%inv_ovlp)
  call move_alloc(ham, self%ham)

  if (use_profiling) call start_profiling("add_det_noci")
  call self%phi%add_det(det, ci_coeff)
  if (use_profiling) call end_profiling("add_det_noci")

  if (use_profiling) call end_profiling(proc_name)
end subroutine expand_noci_selector


!******************************************************************************** 
!
! Solve generalized eigenvalue problem for the |NOCI> = \sum_{p} c_p |p>
!
!   HC = ESC
!
!******************************************************************************** 
subroutine diagonalize_noci(ovlp, ham, ci_coeff, ci_energy)
  complex(wp), intent(in)               :: ovlp(:,:)
  complex(wp), intent(in)               :: ham(:,:)
  complex(wp), allocatable, intent(out) :: ci_coeff(:)
  real(wp), allocatable, intent(out)    :: ci_energy(:)
  !local
  integer                               :: n
  complex(wp), allocatable              :: ci_coeff_all(:,:), ovlp_(:,:)

  n = size(ovlp, 1)

  allocate(ci_coeff(n))  
  allocate(ci_energy(n))
  allocate(ci_coeff_all, source=ham)
  allocate(ovlp_, source=ovlp)

  call gsyev(ci_coeff_all, ovlp_, ci_energy)
  ci_coeff = ci_coeff_all(:,1)
end subroutine diagonalize_noci


!******************************************************************************** 
!
! Perform overlap test on the candidate determinant |mu>
!
! Overlap test compares ||Q|mu>|| / || |mu> ||,
! where Q is an appropriate projection operator
!
! See following reference for more details:
!
!******************************************************************************** 
function overlap_test_noci(self, ci_el) result(lovlp)
  class(NOCISelector), intent(inout)   :: self
  type(NOCIMAtrixElements), intent(in) :: ci_el
  logical                              :: lovlp
  !local  
  character(len=*), parameter          :: proc_name = "overlap_test_noci"

  if (profile_code) call start_profiling(proc_name)  

  self%ovlp_residual = dot_product(ci_el%col_ovlp, matmul(self%inv_ovlp, ci_el%col_ovlp))
  self%ovlp_residual = sqrt(1.0_wp - real(self%ovlp_residual/ci_el%diag_ovlp, wp))
  lovlp = self%ovlp_residual > self%des%ovlp_thresh

  if (profile_code) call end_profiling(proc_name)  
end function overlap_test_noci


!******************************************************************************** 
!
! Perform energy test on the candidate determinant |mu>
!
! Energy test estimates the variaional energy of:
!
!    |Phi> = c_1 |NOCI> + c_2 |mu>
!
! See following reference for more details:
!
!******************************************************************************** 
function energy_test_noci(self, ci_el) result(lenergy)
  class(NOCISelector), intent(inout)   :: self
  type(NOCIMAtrixElements), intent(in) :: ci_el
  logical                              :: lenergy
  !local
  real(wp), allocatable                :: ci_energy(:)
  complex(wp)                          :: ovlp(2,2), ham(2,2)
  complex(wp), allocatable             :: ci_coeff(:)
  character(len=*), parameter          :: proc_name = "energy_test_noci"

  if (profile_code) call start_profiling(proc_name)

  ovlp(1,1) = dot_product(self%phi%ci_coeff, matmul(self%ovlp, self%phi%ci_coeff))
  ovlp(2,1) = zeroc
  ovlp(1,2) = zeroc
  ovlp(2,2) = ci_el%diag_ovlp - dot_product(ci_el%col_ovlp, matmul(self%inv_ovlp, ci_el%col_ovlp))

  ham(1,1) = self%ci_energy(1) * ovlp(1,1)
  ham(2,1) = conjg(dot_product(self%phi%ci_coeff, ci_el%col_ham-matmul(self%ham, matmul(self%inv_ovlp, ci_el%col_ovlp))))
  ham(1,2) = conjg(ham(2,1))
  ham(2,2) = ci_el%diag_ham - 2.0_wp*real(dot_product(ci_el%col_ovlp, matmul(self%inv_ovlp, ci_el%col_ham)), wp) + &
                 dot_product(ci_el%col_ovlp, matmul(self%inv_ovlp, matmul(self%ham, matmul(self%inv_ovlp, ci_el%col_ovlp))))

  if (use_profiling) call start_profiling("diag_noci")
  call diagonalize_noci(ovlp, ham, ci_coeff, ci_energy)
  if (use_profiling) call end_profiling("diag_noci")

  self%energy_residual = abs((ci_energy(1)-self%ci_energy(1)) / self%ci_energy(1))
  lenergy = self%energy_residual > self%des%energy_thresh

  if (profile_code) call end_profiling(proc_name)
end function energy_test_noci


!******************************************************************************** 
!
! Calculate the weights of the determinants in the NOCI wave function
!
!    a_i = (S^{1/2} c)_i
!    w_i = |a_i|^2
!
!    \sum_i w_i = 1 automatically fulfilled
!
! The alternative is to compute
!
!    a_i = (Sc)_i
!    w_i = |a_i|^2 / \sum_j |a_j|^2
! 
!    this version relies on a_i = <i|Psi>
!
!******************************************************************************** 
function get_noci_weights(self) result(weights)
  class(NOCISelector), intent(in) :: self
  real(wp), allocatable           :: weights(:)
  !local
  complex(wp), allocatable        :: ovlp_sqrt(:,:)

  allocate(weights(size(self%ovlp, 1)))
  weights = abs(matmul(self%ovlp, self%phi%ci_coeff))**2
  weights = weights / sum(weights)

  !allocate(weights(size(self%ovlp, 1)))
  !allocate(ovlp_sqrt, source=self%ovlp)
  !call matrix_power(ovlp_sqrt, 0.5_wp)
  !weights = abs(matmul(ovlp_sqrt, self%phi%ci_coeff))**2
end function get_noci_weights


!******************************************************************************** 
!
! Calculate the CI weights of the determinants in the NOCI wave function
!
!    w_i = |c_i|^2 / \sum_j |c_j|^2
!
!******************************************************************************** 
function get_noci_ci_weights(self) result(weights)
  class(NOCISelector), intent(in) :: self
  real(wp), allocatable           :: weights(:)

  allocate(weights(size(self%ovlp, 1)))
  weights = abs(self%phi%ci_coeff)**2
  weights = weights / sum(weights)
end function get_noci_ci_weights


!******************************************************************************** 
!
! Compress NOCI wave function by leaving only the determinants with 
!
!    weight > tol
!
!******************************************************************************** 
subroutine compress_noci_selector(self, tol)
  class(NOCISelector), intent(inout) :: self
  real(wp), intent(in)               :: tol
  !local
  integer                            :: i, ii, j, jj, nd, nd_new
  integer, allocatable               :: indices(:)
  real(wp), allocatable              :: weights(:)
  complex(wp), allocatable           :: ovlp(:,:), ham(:,:), coeff(:,:,:,:)

  if (tol < 0.0_wp) return

  weights = self%get_ci_weights()

  nd = size(weights)
  nd_new = count(weights > tol)

  allocate(indices(nd_new))
  indices = pack([(i, i=1,nd)], weights>tol)

  allocate(ovlp(nd_new, nd_new), ham(nd_new, nd_new))
  allocate(coeff(size(self%phi%coeff,1), size(self%phi%coeff,2), size(self%phi%coeff,3), nd_new))

  do j = 1, nd_new
    jj = indices(j)
    do i = 1, nd_new
      ii = indices(i)
      ovlp(i,j) = self%ovlp(ii,jj)
      ham(i,j) = self%ham(ii,jj)
    end do
    coeff(:,:,:,j) = self%phi%coeff(:,:,:,jj)
  end do

  call move_alloc(ovlp, self%ovlp)
  call move_alloc(ham, self%ham)
  call move_alloc(coeff, self%phi%coeff)

  call diagonalize_noci(self%ovlp, self%ham, self%phi%ci_coeff, self%ci_energy)
end subroutine compress_noci_selector

end module noci_selector