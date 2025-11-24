! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module nosd

!******************************************************************************** 
!
! Module implemetns routines needed to calculate various matrix element
! between non-orthogonal Slater determinants. 
!
! It is counterpart of the slater_condon_rules which implement matrix elements
! between orthogonal Slater determinants
!
!******************************************************************************** 

#include "preproc.inc"
use constants
use lapack 
use profiling
use mpi, only: mpi_communicator, comm_world
use qmcfort_io
use standalone, only: trace, trace2, get_nblocks, dump_array
use hamilton_vars
use hamilton_layout

implicit none

public

interface get_ovlp_matrix
  module procedure get_ovlp_matrix_rc, get_ovlp_matrix_c
end interface get_ovlp_matrix

interface get_ovlp_batch
  module procedure get_ovlp_batch_rc, get_ovlp_batch_c
end interface get_ovlp_batch

interface biorth_orb
  module procedure biorth_orb_rc, biorth_orb_c
end interface biorth_orb

interface biorth_orb_batch
  module procedure biorth_orb_batch_rc, biorth_orb_batch_c
end interface biorth_orb_batch

interface h1_mean_no_batch
  module procedure h1_mean_no_batch_rc, h1_mean_no_batch_c
end interface h1_mean_no_batch

interface h2_gres_mean_no_batch
  module procedure h2_gres_mean_no_batch_rc, h2_gres_mean_no_batch_c
end interface h2_gres_mean_no_batch

interface exchange_energy_no_chol_batch
  module procedure exchange_energy_no_chol_batch_r, exchange_energy_no_chol_batch_c
end interface exchange_energy_no_chol_batch

interface exchange_energy_no_int_batch
  module procedure exchange_energy_no_int_batch_r, exchange_energy_no_int_batch_c
end interface exchange_energy_no_int_batch

interface get_coeff_occ
  module procedure get_coeff_occ_r, get_coeff_occ_c
end interface get_coeff_occ

interface get_h1phi
  module procedure get_h1phi_r, get_h1phi_cr
end interface get_h1phi

interface get_h2phi
  module procedure get_h2phi_r, get_h2phi_cr
end interface get_h2phi

interface get_hselfphi
  module procedure get_hselfphi_r, get_hselfphi_cr
end interface get_hselfphi

interface get_h1phi_occ
  module procedure get_h1phi_occ_r, get_h1phi_occ_c
end interface get_h1phi_occ

interface get_h2phi_h_occ
  module procedure get_h2phi_h_occ_r, get_h2phi_h_occ_c
end interface get_h2phi_h_occ

interface get_h2phi_x_occ_chol
  module procedure get_h2phi_x_occ_chol_r, get_h2phi_x_occ_chol_c
end interface get_h2phi_x_occ_chol

interface get_h2phi_x_occ_int
  module procedure get_h2phi_x_occ_int_r, get_h2phi_x_occ_int_c
end interface get_h2phi_x_occ_int

contains 

!******************************************************************************** 
!
! Reorder orbitals from the batch of the determinants in the following manner
!
!    From:
!        -------------------------------------
!        |    |   |    |   |    |   |    |   |
!        |    |   |    |   |    |   |    |   |
!        |1up |1dn|2up |2dn|3up |3dn|4up |4dn|
!        |    |   |    |   |    |   |    |   |
!        |    |   |    |   |    |   |    |   |
!        -------------------------------------
!
!    To:
!        -------------------------------------
!        |    |    |    |    |   |   |   |   |
!        |    |    |    |    |   |   |   |   |
!        |1up |2up |3up |4up |1dn|2dn|3dn|4dn|
!        |    |    |    |    |   |   |   |   |
!        |    |    |    |    |   |   |   |   |
!        -------------------------------------
!
!******************************************************************************** 
subroutine reorder_coeff(hdes, coeff, coeff_)
  type(hamil_descriptor), intent(in)      :: hdes
  complex(wp), intent(in)                 :: coeff(:,:)
  complex(wp), allocatable, intent(inout) :: coeff_(:,:)
  !local variables
  integer                                 :: nw, nocc, nel
  integer                                 :: spin, w, i1, i2, j1, j2

  nocc = hdes%nocc
  nw = size(coeff, 2) / nocc

  if (hdes%ispin == 1) then
    allocate(coeff_, source=coeff)
  else 
    allocate(coeff_, mold=coeff)

    do w = 1, nw 
      do spin = 1, hdes%ispin
        nel = hdes%nel(spin)

        i1 = (w-1)*nocc + hdes%nell(spin) + 1
        i2 = (w-1)*nocc + hdes%nell(spin) + nel

        j1 = nw * hdes%nell(spin) + (w-1)*nel + 1 
        j2 = nw * hdes%nell(spin) + w*nel 

        coeff_(:,j1:j2) = coeff(:,i1:i2)
      end do
    end do
  end if 
end subroutine reorder_coeff


!******************************************************************************** 
!
! Allocate space and double spin if needed for biorthogonal orbitals
!
!******************************************************************************** 
subroutine get_coeff_spin(hdes, coeff, coeff_bi)
  type(hamil_descriptor), intent(in)    :: hdes
  complex(wp), intent(in)               :: coeff(:,:)
  complex(wp), allocatable, intent(out) :: coeff_bi(:,:)

  if (hdes%nocc == hdes%nocc_) then
    allocate(coeff_bi(size(coeff,1), size(coeff,2)))
    coeff_bi = coeff
  else
    allocate(coeff_bi(size(coeff,1), 2*size(coeff,2)))
    coeff_bi(:,1:size(coeff,2)) = coeff
    coeff_bi(:,size(coeff,2)+1:2*size(coeff,2)) = coeff
  end if
end subroutine get_coeff_spin


!******************************************************************************** 
!
! Calculates overlap matrix for nw set of orbitals
!
!    Input:
!        coeff_l - orbitals in the left Slater determinant
!        coeff_r  - nw set of oribtals in the right Slater determinant
!
!    Output:
!        ovlp - nw overlap matrices (ne_spin x ne_spin*nw)
!
!    ovlp = coeff_l^{+} coeff_r
!
!******************************************************************************** 
subroutine get_ovlp_matrix_rc(coeff_l, coeff_r, ovlp)
  real(wp), intent(in)                  :: coeff_l(:,:)
  complex(wp), intent(in)               :: coeff_r(:,:)
  complex(wp), allocatable, intent(out) :: ovlp(:,:)
  !local variables
  integer                               :: n, ne, ne_
  character(len=*), parameter           :: proc_name = "ovlp_matrix"
  
 if (profile_code) call start_profiling(proc_name)

  n = size(coeff_r, 1)
  ne = size(coeff_r, 2)
  ne_ = size(coeff_l, 2)
  
  allocate(ovlp(ne_,ne))
  
  call gemm("t", "n", ne_, ne, n, onec, coeff_l, n, coeff_r, n, zeroc, ovlp, ne_)

  if (profile_code) call end_profiling(proc_name, ne_, ne, n, "rc")
end subroutine get_ovlp_matrix_rc

subroutine get_ovlp_matrix_c(coeff_l, coeff_r, ovlp)
  complex(wp), intent(in)               :: coeff_l(:,:)
  complex(wp), intent(in)               :: coeff_r(:,:)
  complex(wp), allocatable, intent(out) :: ovlp(:,:)
  !local variables
  integer                               :: n, ne, ne_
  character(len=*), parameter           :: proc_name = "ovlp_matrix"
  
  if (profile_code) call start_profiling(proc_name)

  n = size(coeff_r, 1)
  ne = size(coeff_r, 2)
  ne_ = size(coeff_l, 2)
  
  allocate(ovlp(ne_,ne))
  
  call gemm("c", "n", ne_, ne, n, onec, coeff_l, n, coeff_r, n, zeroc, ovlp, ne_)

  if (profile_code) call end_profiling(proc_name, ne_, ne, n, "c")
end subroutine get_ovlp_matrix_c


!******************************************************************************** 
!
! Calculate overlap for nw set of orbitals stored in coeff
!
!    <Phi_l|Phi_r^w>  = det(coeff_l^{+} coeff_r^w)
!
!******************************************************************************** 
subroutine get_ovlp_batch_rc(hdes, coeff_l, coeff_r, ovlp)
  type(hamil_descriptor), intent(in)    :: hdes
  real(wp), intent(in)                  :: coeff_l(:,:)
  complex(wp), intent(in)               :: coeff_r(:,:)
  complex(wp), allocatable, intent(out) :: ovlp(:)
  !local variables
  integer                               :: nw, w, w1, w2, w1_, w2_, spin, nel
  complex(wp)                           :: ovlp_
  complex(wp), allocatable              :: coeff_bi(:,:), ovlp_matrix(:,:)
  character(len=*), parameter           :: proc_name = "get_overlap"

  if (profile_code) call start_profiling(proc_name)

  nw = size(coeff_r, 2) / hdes%nocc

  allocate(ovlp(nw))
  ovlp = onec

  call get_coeff_spin(hdes, coeff_r, coeff_bi)

  do spin = 1, hdes%ispin
    nel = hdes%nel(spin)
    if (nel == 0) cycle
    w1 = hdes%nell(spin) + 1
    w2 = hdes%nell(spin) + nel
    w1_ = hdes%nell(spin)*nw + 1
    w2_ = hdes%nell(spin)*nw + nel*nw
    call get_ovlp_matrix(coeff_l(:,w1:w2), coeff_bi(:,w1_:w2_), ovlp_matrix)
    do w = 1, nw
      w1 = (w-1)*nel + 1
      w2 = w*nel
      call determinant(ovlp_matrix(1:nel,w1:w2), ovlp_)
      ovlp(w) = ovlp(w) * ovlp_**(hdes%rspin)
    end do
  end do

  if (profile_code) call end_profiling(proc_name)
end subroutine get_ovlp_batch_rc

subroutine get_ovlp_batch_c(hdes, coeff_l, coeff_r, ovlp)
  type(hamil_descriptor), intent(in)    :: hdes
  complex(wp), intent(in)               :: coeff_l(:,:)
  complex(wp), intent(in)               :: coeff_r(:,:)
  complex(wp), allocatable, intent(out) :: ovlp(:)
  !local variables
  integer                               :: nw, w, w1, w2, w1_, w2_, spin, nel
  complex(wp)                           :: ovlp_
  complex(wp), allocatable              :: coeff_bi(:,:), ovlp_matrix(:,:)
  character(len=*), parameter           :: proc_name = "get_overlap"

  if (profile_code) call start_profiling(proc_name)

  nw = size(coeff_r, 2) / hdes%nocc

  allocate(ovlp(nw))
  ovlp = onec

  call get_coeff_spin(hdes, coeff_r, coeff_bi)

  do spin = 1, hdes%ispin
    nel = hdes%nel(spin)
    if (nel == 0) cycle
    w1 = hdes%nell(spin) + 1
    w2 = hdes%nell(spin) + nel
    w1_ = hdes%nell(spin)*nw + 1
    w2_ = hdes%nell(spin)*nw + nel*nw
    call get_ovlp_matrix(coeff_l(:,w1:w2), coeff_bi(:,w1_:w2_), ovlp_matrix)
    do w = 1, nw
      w1 = (w-1)*nel + 1
      w2 = w*nel
      call determinant(ovlp_matrix(1:nel,w1:w2), ovlp_)
      ovlp(w) = ovlp(w) * ovlp_**(hdes%rspin)
    end do
  end do

  if (profile_code) call end_profiling(proc_name)
end subroutine get_ovlp_batch_c


!******************************************************************************** 
!
! Biorthogonalize right set of orbitals:
!
!    Input:
!        hdes     - Hamiltonian descriptor
!        coeff_l  - orbitals in the left Slater determinant
!        coeff_r  - oribtals in the right Slater determinant
!
!    Output:
!        coeff_bi - biorthogonalized orbitals with resepect to coeff_l
!            
!        coeff_bi = coeff_r * S^{-1} with
!        S = coeff_l^{+} coeff_r
!    
!******************************************************************************** 
subroutine biorth_orb_rc(hdes, coeff_l, coeff_r, coeff_bi)
  type(hamil_descriptor), intent(in)    :: hdes
  real(wp), intent(in)                  :: coeff_l(:,:)
  complex(wp), intent(in)               :: coeff_r(:,:)
  complex(wp), allocatable, intent(out) :: coeff_bi(:,:)
  !local variables
  integer                               :: spin, w1, w2, n, nel
  complex(wp), allocatable              :: ovlp_matrix(:,:), coeff_r_sp(:,:)
  character(len=*), parameter           :: proc_name = "biorth_orb"
  
  if (profile_code) call start_profiling(proc_name)

  call get_coeff_spin(hdes, coeff_r, coeff_r_sp)
  allocate(coeff_bi, mold=coeff_r_sp)

  n = hdes%n
  do spin = 1, hdes%ispin
    nel = hdes%nel(spin)
    if (nel == 0) cycle
    w1 = hdes%nell(spin) + 1 
    w2 = hdes%nell(spin) + nel
    call get_ovlp_matrix(coeff_l(:,w1:w2), coeff_r_sp(:,w1:w2), ovlp_matrix)
    call inverse(ovlp_matrix)
    call gemm("n", "n", n, nel, nel, onec, coeff_r_sp(:,w1:w2), n, ovlp_matrix, nel, zeroc, coeff_bi(:,w1:w2), n)
  end do

  if (profile_code) call end_profiling(proc_name)
end subroutine biorth_orb_rc

subroutine biorth_orb_c(hdes, coeff_l, coeff_r, coeff_bi)
  type(hamil_descriptor), intent(in)    :: hdes
  complex(wp), intent(in)               :: coeff_l(:,:)
  complex(wp), intent(in)               :: coeff_r(:,:)
  complex(wp), allocatable, intent(out) :: coeff_bi(:,:)
  !local variables
  integer                               :: spin, w1, w2, n, nel
  complex(wp), allocatable              :: ovlp_matrix(:,:), coeff_r_sp(:,:)
  character(len=*), parameter           :: proc_name = "biorth_orb"
  
  if (profile_code) call start_profiling(proc_name)

  call get_coeff_spin(hdes, coeff_r, coeff_r_sp)
  allocate(coeff_bi, mold=coeff_r_sp)

  n = hdes%n
  do spin = 1, hdes%ispin
    nel = hdes%nel(spin)
    if (nel == 0) cycle
    w1 = hdes%nell(spin) + 1 
    w2 = hdes%nell(spin) + nel
    call get_ovlp_matrix(coeff_l(:,w1:w2), coeff_r_sp(:,w1:w2), ovlp_matrix)
    call inverse(ovlp_matrix)
    call gemm("n", "n", n, nel, nel, onec, coeff_r_sp(:,w1:w2), n, ovlp_matrix, nel, zeroc, coeff_bi(:,w1:w2), n)
  end do

  if (profile_code) call end_profiling(proc_name)
end subroutine biorth_orb_c


!******************************************************************************** 
!
! Biorthogonalize nw set of orbitals:
!
!    Input:
!        hdes     - Hamiltonian descriptor
!        coeff_l  - orbitals in the left Slater determinant
!        coeff_r  - Nw set of oribtals in the right Slater determinant
!        mode     - determines the data layout of coeff_r
!
!    Output:
!        coeff_bi - biorthogonalized orbitals with respect to coeff_l
!            
!        coeff_bi = coeff * S^{-1} with
!        S = coeff_l^{+} coeff_r
!
! Data layout:
!
!   default mode (mode not present) 
!       -------------------------------------------------------------
!       |       |       |     |       |       |       |     |       |
!       |       |       |     |       |       |       |     |       |
!       | phi_1 | phi_2 | ... | phi_n | phi_1 | phi_2 | ... | phi_n |
!       |  _up  |  _up  |     |  _up  |  _dn  |  _dn  |     |  _dn  |
!       |       |       |     |       |       |       |     |       |
!       -------------------------------------------------------------
!
!   if present mode
!       -------------------------------------------------------
!       |       |       |       |       |     |       |       |
!       |       |       |       |       |     |       |       |
!       | phi_1 | phi_1 | phi_2 | phi_2 | ... | phi_n | phi_n |
!       |  _up  |  _dn  |  _up  |  _dn  |     |  _up  |  _dn  |
!       |       |       |       |       |     |       |       |
!       -------------------------------------------------------
!    
!******************************************************************************** 
subroutine biorth_orb_batch_rc(hdes, coeff_l, coeff_r, coeff_bi, mode)
  type(hamil_descriptor), intent(in)    :: hdes
  real(wp), intent(in)                  :: coeff_l(:,:)
  complex(wp), intent(in)               :: coeff_r(:,:)
  complex(wp), allocatable, intent(out) :: coeff_bi(:,:)
  integer, optional, intent(in)         :: mode
  !local variables
  integer                               :: nw, spin, w, w1, w2, w1_, w2_, n, nel
  complex(wp), allocatable              :: ovlp_matrix(:,:), coeff_r_sp(:,:)
  character(len=*), parameter           :: proc_name = "biorth_orb"
  
  if (profile_code) call start_profiling(proc_name)

  nw = size(coeff_r, 2) / hdes%nocc
  n = hdes%n

  call get_coeff_spin(hdes, coeff_r, coeff_r_sp)
  allocate(coeff_bi, mold=coeff_r_sp)

  if (present(mode)) then
    w1_ = 0
    w2_ = 0
    do w = 1, nw
      do spin = 1, hdes%ispin
        nel = hdes%nel(spin)
        if (nel == 0) cycle
        w1 = hdes%nell(spin) + 1
        w2 = hdes%nell(spin) + nel
        w1_ = w2_ + 1
        w2_ = w2_ + nel
        call get_ovlp_matrix(coeff_l(:,w1:w2), coeff_r_sp(:,w1_:w2_), ovlp_matrix)
        call inverse(ovlp_matrix)
        call gemm("n", "n", n, nel, nel, onec, coeff_r_sp(:,w1_:w2_), n, ovlp_matrix, nel, zeroc, coeff_bi(:,w1_:w2_), n)
      end do
    end do
  else
    w1_ = 0 
    w2_ = 0
    do spin = 1, hdes%ispin
      nel = hdes%nel(spin)
      if (nel == 0) cycle
      w1 = hdes%nell(spin) + 1 
      w2 = hdes%nell(spin) + nel
      do w = 1, nw
        w1_ = w2_ + 1
        w2_ = w2_ + nel
        call get_ovlp_matrix(coeff_l(:,w1:w2), coeff_r_sp(:,w1_:w2_), ovlp_matrix)
        if (profile_code) call start_profiling("inverse_ovlp")
!debug: pseudoinv
        call inverse(ovlp_matrix)
        !call pseudoinverse(ovlp_matrix, 0.01_wp)
        if (profile_code) call end_profiling("inverse_ovlp")
        if (profile_code) call start_profiling("gemm_biorth")
        call gemm("n", "n", n, nel, nel, onec, coeff_r_sp(:,w1_:w2_), n, ovlp_matrix, nel, zeroc, coeff_bi(:,w1_:w2_), n)
        if (profile_code) call end_profiling("gemm_biorth", n, nel, nel, "c")
      end do
    end do
  end if

  if (profile_code) call end_profiling(proc_name)
end subroutine biorth_orb_batch_rc

subroutine biorth_orb_batch_c(hdes, coeff_l, coeff_r, coeff_bi, mode)
  type(hamil_descriptor), intent(in)    :: hdes
  complex(wp), intent(in)               :: coeff_l(:,:)
  complex(wp), intent(in)               :: coeff_r(:,:)
  complex(wp), allocatable, intent(out) :: coeff_bi(:,:)
  integer, optional, intent(in)         :: mode
  !local variables
  integer                               :: nw, spin, w, w1, w2, w1_, w2_, n, nel
  complex(wp), allocatable              :: ovlp_matrix(:,:), coeff_r_sp(:,:)
  character(len=*), parameter           :: proc_name = "biorth_orb"

  
  if (profile_code) call start_profiling(proc_name)

  nw = size(coeff_r, 2) / hdes%nocc
  n = hdes%n

  call get_coeff_spin(hdes, coeff_r, coeff_r_sp)
  allocate(coeff_bi, mold=coeff_r_sp)

  if (present(mode)) then
    w1_ = 0
    w2_ = 0
    do w = 1, nw
      do spin = 1, hdes%ispin
        nel = hdes%nel(spin)
        if (nel == 0) cycle
        w1 = hdes%nell(spin) + 1
        w2 = hdes%nell(spin) + nel
        w1_ = w2_ + 1
        w2_ = w2_ + nel
        call get_ovlp_matrix(coeff_l(:,w1:w2), coeff_r_sp(:,w1_:w2_), ovlp_matrix)
        call inverse(ovlp_matrix)
        call gemm("n", "n", n, nel, nel, onec, coeff_r_sp(:,w1_:w2_), n, ovlp_matrix, nel, zeroc, coeff_bi(:,w1_:w2_), n)
      end do
    end do
  else
    w1_ = 0 
    w2_ = 0
    do spin = 1, hdes%ispin
      nel = hdes%nel(spin)
      if (nel == 0) cycle
      w1 = hdes%nell(spin) + 1 
      w2 = hdes%nell(spin) + nel
      do w = 1, nw
        w1_ = w2_ + 1
        w2_ = w2_ + nel
        call get_ovlp_matrix(coeff_l(:,w1:w2), coeff_r_sp(:,w1_:w2_), ovlp_matrix)
        if (profile_code) call start_profiling("inverse_ovlp")
        call inverse(ovlp_matrix)
        if (profile_code) call end_profiling("inverse_ovlp")
        if (profile_code) call start_profiling("gemm_biorth")
        call gemm("n", "n", n, nel, nel, onec, coeff_r_sp(:,w1_:w2_), n, ovlp_matrix, nel, zeroc, coeff_bi(:,w1_:w2_), n)
        if (profile_code) call end_profiling("gemm_biorth", n, nel, nel, "c")
      end do
    end do
  end if

  if (profile_code) call end_profiling(proc_name)
end subroutine biorth_orb_batch_c


!******************************************************************************** 
!
! Calculation of the one-body energy for nw set of orbitals
!
!    mean_w = sum_{pi} h1phi_{pi} coeff_{pi}(w)
!
!    Input:
!        hdes  - Hamilton descriptor
!        h1phi - intermediate one-body array (see get_h1phi_occ)
!        coeff - nw sets of orbitals (biorthogonalized orbitals)
!
!    Output:
!        mean - array of expectation values
!
!******************************************************************************** 
subroutine h1_mean_no_batch_rc(hdes, h1phi, coeff_bi, mean)
  type(hamil_descriptor), intent(in)          :: hdes
  real(wp), target, contiguous, intent(in)    :: h1phi(:,:)
  complex(wp), target, contiguous, intent(in) :: coeff_bi(:,:)
  complex(wp), allocatable, intent(out)       :: mean(:)
  !local variables
  integer                                     :: n, nw, nel, nnel
  integer                                     :: w1, w2, i1, i2, spin
  complex(wp)                                 :: alpha
  complex(wp), allocatable                    :: mean_(:,:)
  real(wp), pointer                           :: h1phi_ptr(:,:)
  complex(wp), pointer                        :: coeff_ptr(:,:)
  character(len=*), parameter                 :: proc_name = "h1_mean_no"

  if (profile_code) call start_profiling(proc_name)

  n = size(h1phi, 1)
  nw = size(coeff_bi, 2) / hdes%nocc_

  allocate(mean_(1,nw))
  mean_ = zeroc

  do spin = 1, hdes%ispin
    nel = hdes%nel(spin)
    if (nel == 0) cycle
    nnel = n * nel
    i1 = hdes%nell(spin) + 1
    i2 = hdes%nell(spin) + nel
    w1 = hdes%nell(spin)*nw + 1
    w2 = (hdes%nell(spin)+nel) * nw
    h1phi_ptr(1:1,1:nnel) => h1phi(:,i1:i2)
    coeff_ptr(1:nnel,1:nw) => coeff_bi(:,w1:w2)

    alpha = onec * hdes%rspin
    call gemm("n", "n", 1, nw, nnel, alpha, h1phi_ptr, 1, coeff_ptr, nnel, onec, mean_, 1)
  end do

  allocate(mean(nw))
  mean = mean_(1,:)

  if (profile_code) call end_profiling(proc_name)
end subroutine h1_mean_no_batch_rc

subroutine h1_mean_no_batch_c(hdes, h1phi, coeff_bi, mean)
  type(hamil_descriptor), intent(in)          :: hdes
  complex(wp), target, contiguous, intent(in) :: h1phi(:,:)
  complex(wp), target, contiguous, intent(in) :: coeff_bi(:,:)
  complex(wp), allocatable, intent(out)       :: mean(:)
  !local variables
  integer                                     :: n, nw, nel, nnel
  integer                                     :: w1, w2, i1, i2, spin
  complex(wp)                                 :: alpha
  complex(wp), allocatable                    :: mean_(:,:)
  complex(wp), pointer                        :: h1phi_ptr(:,:), coeff_ptr(:,:)
  character(len=*), parameter                 :: proc_name = "h1_mean_no"

  if (profile_code) call start_profiling(proc_name)

  n = size(h1phi, 1)
  nw = size(coeff_bi, 2) / hdes%nocc_

  allocate(mean_(1,nw))
  mean_ = zeroc

  do spin = 1, hdes%ispin
    nel = hdes%nel(spin)
    if (nel == 0) cycle
    nnel = n * nel
    i1 = hdes%nell(spin) + 1
    i2 = hdes%nell(spin) + nel
    w1 = hdes%nell(spin)*nw + 1
    w2 = (hdes%nell(spin)+nel) * nw
    h1phi_ptr(1:1,1:nnel) => h1phi(:,i1:i2)
    coeff_ptr(1:nnel,1:nw) => coeff_bi(:,w1:w2)

    alpha = onec * hdes%rspin
    call gemm("n", "n", 1, nw, nnel, alpha, h1phi_ptr, 1, coeff_ptr, nnel, onec, mean_, 1)
  end do

  allocate(mean(nw))
  mean = mean_(1,:)

  if (profile_code) call end_profiling(proc_name)
end subroutine h1_mean_no_batch_c


!******************************************************************************** 
!
! Calculates exp. values of ng one-body operators for nw set of orbitals
!
!    force_{gw} = \sum_{I} h2phi_h{gI} coeff_{Iw}
! 
! with composite index I = {p,i}
!
!    Input:
!        hdes     - Hamilton descriptor
!        h2phi_h  - intermediate two-body array (see get_h2phi_occ_h)
!        coeff_bi - nw sets of orbitals (biorthogonalized orbitals)
!
!    Output:
!        lgmean - array of expectation values
!
!******************************************************************************** 
subroutine h2_gres_mean_no_batch_rc(hdes, h2phi_h, coeff_bi, lgmean, bsize)
  type(hamil_descriptor), intent(in)    :: hdes
  real(wp), intent(in)                  :: h2phi_h(:,:)
  complex(wp), intent(in)               :: coeff_bi(:,:)
  complex(wp), allocatable, intent(out) :: lgmean(:,:)
  integer, optional, intent(in)         :: bsize
  !local variables
  integer                               :: n, nel, ng, nw, nw_, bsize_, nblocks
  integer                               :: w1, w2, i1, i2, spin, hblock
  complex(wp), allocatable              :: coeff_(:,:)
  character(len=*), parameter           :: proc_name = "lgmean_no"

  if (profile_code) call start_profiling(proc_name)

  n = size(coeff_bi, 1)
  ng = size(h2phi_h, 1)
  nw = size(coeff_bi, 2) / hdes%nocc_

  if (present(bsize)) then
    bsize_ = bsize
  else
    bsize_ = nw
  end if 

  allocate(lgmean(ng, nw))
  lgmean = zeroc

  do spin = 1, hdes%ispin
    nel = hdes%nel(spin)
    if (nel == 0) cycle
    i1 = hdes%nell(spin)*n + 1
    i2 = (hdes%nell(spin) + hdes%nel(spin)) * n
    w1 = hdes%nell(spin)*nw + 1
    w2 = (hdes%nell(spin) + hdes%nel(spin)) * nw
!debug: opt - get rid of the allocation and reshaping
    if (allocated(coeff_)) deallocate(coeff_)
    allocate(coeff_(n*nel, nw))
    coeff_ = reshape(coeff_bi(:,w1:w2), shape=[n*nel, nw])

    if (is_blocking(nw, bsize_)) then
      nblocks = get_nblocks(nw, bsize_)
      !$omp parallel do private(w1, w2, nw_)
      do hblock = 1, nblocks
        w1 = (hblock-1)*bsize_ + 1
        w2 = min(hblock*bsize_, nw)
        nw_ = w2 - w1 + 1
        if (profile_code) call start_profiling("gemm_lgmean")
        call gemm("n", "n", ng, nw_, n*nel, onec*hdes%rspin, h2phi_h(:,i1:i2), ng,  & 
                  coeff_(:,w1:w2), n*nel, onec, lgmean(:,w1:w2), ng)
        if (profile_code) call end_profiling("gemm_lgmean", ng, n*nel, nw_, "rc")
      end do
      !$omp end parallel do
    else
      if (profile_code) call start_profiling("gemm_lgmean")
      call gemm("n", "n", ng, nw, n*nel, onec*hdes%rspin, h2phi_h(:,i1:i2), ng,  & 
                coeff_, n*nel, onec, lgmean, ng)
      if (profile_code) call end_profiling("gemm_lgmean", ng, n*nel, nw, "rc")
    end if
  end do

  if (profile_code) call end_profiling(proc_name)
end subroutine h2_gres_mean_no_batch_rc

subroutine h2_gres_mean_no_batch_c(hdes, h2phi_h, coeff_bi, lgmean, bsize)
  type(hamil_descriptor), intent(in)    :: hdes
  complex(wp), intent(in)               :: h2phi_h(:,:)
  complex(wp), intent(in)               :: coeff_bi(:,:)
  complex(wp), allocatable, intent(out) :: lgmean(:,:)
  integer, optional, intent(in)         :: bsize
  !local variables
  integer                               :: n, nel, ng, nw, nw_, bsize_, nblocks
  integer                               :: w1, w2, i1, i2, spin, hblock
  complex(wp), allocatable              :: coeff_(:,:)
  character(len=*), parameter           :: proc_name = "lgmean_no"

  if (profile_code) call start_profiling(proc_name)

  n = size(coeff_bi, 1)
  ng = size(h2phi_h, 1)
  nw = size(coeff_bi, 2) / hdes%nocc_

  if (present(bsize)) then
    bsize_ = bsize
  else
    bsize_ = nw
  end if 

  allocate(lgmean(ng, nw))
  lgmean = zeroc

  do spin = 1, hdes%ispin
    nel = hdes%nel(spin)
    if (nel == 0) cycle
    i1 = hdes%nell(spin)*n + 1
    i2 = (hdes%nell(spin) + hdes%nel(spin)) * n
    w1 = hdes%nell(spin)*nw + 1
    w2 = (hdes%nell(spin) + hdes%nel(spin)) * nw
!debug: opt - get rid of the allocation and reshaping
    if (allocated(coeff_)) deallocate(coeff_)
    allocate(coeff_(n*nel, nw))
    coeff_ = reshape(coeff_bi(:,w1:w2), shape=[n*nel, nw])

    if (is_blocking(nw, bsize_)) then
      nblocks = get_nblocks(nw, bsize_)
      !$omp parallel do private(w1, w2, nw_)
      do hblock = 1, nblocks
        w1 = (hblock-1)*bsize_ + 1
        w2 = min(hblock*bsize_, nw)
        nw_ = w2 - w1 + 1
        if (profile_code) call start_profiling("gemm_lgmean")
        call gemm("n", "n", ng, nw_, n*nel, onec*hdes%rspin, h2phi_h(:,i1:i2), ng,  & 
                  coeff_(:,w1:w2), n*nel, onec, lgmean(:,w1:w2), ng)
        if (profile_code) call end_profiling("gemm_lgmean", ng, n*nel, nw_, "c")
      end do
      !$omp end parallel do
    else
      if (profile_code) call start_profiling("gemm_lgmean")
      call gemm("n", "n", ng, nw, n*nel, onec*hdes%rspin, h2phi_h(:,i1:i2), ng,  & 
                coeff_, n*nel, onec, lgmean, ng)
      if (profile_code) call end_profiling("gemm_lgmean", ng, n*nel, nw, "c")
    end if
  end do

  if (profile_code) call end_profiling(proc_name)
end subroutine h2_gres_mean_no_batch_c


!******************************************************************************** 
!
! Calculates exchange energy for nw set of orbitals (using Cholesky vectors)
!
!    theta_{g,ij}(w) = \sum_{n} h2phi_x{gi,p} coeff_{p,j}
!    E_x(w) = -0.5 \sum_g \sum_{ij} theta_{g,ij}(w)
! 
!    Input:
!        hdes    - Hamilton descriptor
!        h2phi_x - intermediate two-body array (see get_h2phi_occ_x_chol)
!        coeff_bi- nw sets of orbitals (biorthogonalized orbitals)
!        bsize   - specifies how many walkers are treated simultaneously
!
!    Output:
!        ex - array of exchange energies
!
!******************************************************************************** 
subroutine exchange_energy_no_chol_batch_r(hdes, h2phi_x, coeff_bi, ex, bsize)
  type(hamil_descriptor), intent(in)    :: hdes
  real(wp), intent(in)                  :: h2phi_x(:,:)
  complex(wp), intent(in)               :: coeff_bi(:,:)
  complex(wp), allocatable, intent(out) :: ex(:,:)
  integer, optional, intent(in)         :: bsize
  !local variables
  integer                               :: ispin, n, nocc, nel, nel_, ng, nw, nw_, bsize_
  integer                               :: spin, g, g1, g2, w, w1, w2, i1, i2, j1, j2, k1, k2, xblock
  complex(wp), allocatable              :: help(:,:)
  character(len=*), parameter           :: proc_name = "exchange_energy_no_chol"

  if (profile_code) call start_profiling(proc_name)

  n = hdes%n
  ispin = hdes%ispin
  nocc = hdes%nocc_
  ng = size(h2phi_x, 2) / sum(hdes%nel)
  nw = size(coeff_bi, 2) / nocc

  if (present(bsize)) then
    bsize_ = bsize
  else
    bsize_ = nw
  end if 

  allocate(ex(ng,nw))
  ex = zeroc

  do spin = 1, ispin
    nel = hdes%nel(spin) 
    if (nel == 0) cycle
    i1 = ng * hdes%nell(spin) + 1
    i2 = ng * (hdes%nell(spin) + nel)

    if (is_blocking(nw, bsize_)) then
      if (.not. allocated(help)) allocate(help(ng*maxval(hdes%nel), bsize_*maxval(hdes%nel)))
      help = zeroc

      !$omp parallel do private(w, w1, w2, g, g1, g2, nw_, nel_, j1, j2, k1, k2) firstprivate(help)
      do xblock = 1, get_nblocks(nw, bsize_)
        w1 = (xblock-1) * bsize_
        w2 = min(xblock*bsize_, nw)
        nw_ = w2 - w1
        nel_ = nel * nw_

        j1 = hdes%nell(spin)*nw + w1*nel + 1
        j2 = hdes%nell(spin)*nw + w2*nel

        if (profile_code) call start_profiling("gemm_en")
        call gemm("t", "n", ng*nel, nel_, n, onec, h2phi_x(:,i1:i2), n, coeff_bi(:,j1:j2), n, &
                  zeroc, help(1:ng*nel,1:nel_), ng*nel)
        if (profile_code) call end_profiling("gemm_en", ng*nel, nel_, n, "rc")

        if (profile_code) call start_profiling("alpha_g_ij")
        do w = w1+1, w2
          k1 = (w-w1-1)*nel + 1
          k2 = (w-w1)*nel 
          do g = 1, ng
            g1 = nel * (g-1) + 1
            g2 = nel * g
            !be careful with e_x_g = \sum_ij \alpha_g_ij \alpha_g_ji
            ex(g,w) = ex(g,w) - 0.5_wp * hdes%rspin * trace2(help(g1:g2,k1:k2), help(g1:g2,k1:k2))
          end do
        end do
        if (profile_code) call end_profiling("alpha_g_ij", 8*nel**2*ng*nw_)
      end do
      !$omp end parallel do
    else 
      if (.not. allocated(help)) allocate(help(ng*maxval(hdes%nel), nw*maxval(hdes%nel)))
      help = zeroc

      nel_ = nel * nw
      j1 = hdes%nell(spin)*nw + 1
      j2 = hdes%nell(spin)*nw + nel_

      if (profile_code) call start_profiling("gemm_en")
      call gemm("t", "n", ng*nel, nel_, n, onec, h2phi_x(:,i1:i2), n, coeff_bi(:,j1:j2), n, &
                zeroc, help(1:ng*nel,1:nel_), ng*nel)
      if (profile_code) call end_profiling("gemm_en", ng*nel, nel_, n, "rc")

      if (profile_code) call start_profiling("alpha_g_ij")
      do w = 1, nw
        k1 = (w-1)*nel + 1
        k2 = w*nel 
        do g = 1, ng
          g1 = nel * (g-1) + 1
          g2 = nel * g
          !be careful with e_x_g = \sum_ij \alpha_g_ij \alpha_g_ji
          ex(g,w) = ex(g,w) - 0.5_wp * hdes%rspin * trace2(help(g1:g2,k1:k2), help(g1:g2,k1:k2))
        end do
      end do
      if (profile_code) call end_profiling("alpha_g_ij", 8*nel**2*ng*nw)
    end if 
  end do

  if (profile_code) call end_profiling(proc_name)
end subroutine exchange_energy_no_chol_batch_r

subroutine exchange_energy_no_chol_batch_c(hdes, h2phi_x, coeff_bi, ex, bsize)
  type(hamil_descriptor), intent(in)    :: hdes
  complex(wp), intent(in)               :: h2phi_x(:,:)
  complex(wp), intent(in)               :: coeff_bi(:,:)
  complex(wp), allocatable, intent(out) :: ex(:,:)
  integer, optional, intent(in)         :: bsize
  !local variables
  integer                               :: ispin, n, nocc, nel, nel_, ng, nw, nw_, bsize_
  integer                               :: spin, g, g1, g2, w, w1, w2, i1, i2, j1, j2, k1, k2, xblock
  complex(wp), allocatable              :: help(:,:)
  character(len=*), parameter           :: proc_name = "exchange_energy_no_chol"

  if (profile_code) call start_profiling(proc_name)

  n = hdes%n
  ispin = hdes%ispin
  nocc = hdes%nocc_
  ng = size(h2phi_x, 2) / sum(hdes%nel)
  nw = size(coeff_bi, 2) / nocc

  if (present(bsize)) then
    bsize_ = bsize
  else
    bsize_ = nw
  end if 

  allocate(ex(ng,nw))
  ex = zeroc

  do spin = 1, ispin
    nel = hdes%nel(spin) 
    if (nel == 0) cycle
    i1 = ng * hdes%nell(spin) + 1
    i2 = ng * (hdes%nell(spin) + nel)

    if (is_blocking(nw, bsize_)) then
      if (.not. allocated(help)) allocate(help(ng*maxval(hdes%nel), bsize_*maxval(hdes%nel)))
      help = zeroc

      !$omp parallel do private(w, w1, w2, g, g1, g2, nw_, nel_, j1, j2, k1, k2) firstprivate(help)
      do xblock = 1, get_nblocks(nw, bsize_)
        w1 = (xblock-1) * bsize_
        w2 = min(xblock*bsize_, nw)
        nw_ = w2 - w1
        nel_ = nel * nw_

        j1 = hdes%nell(spin)*nw + w1*nel + 1
        j2 = hdes%nell(spin)*nw + w2*nel

        if (profile_code) call start_profiling("gemm_en")
!debug: cmplx nosd - careful analysis whether we need "t", "n" or "c"
        call gemm("t", "n", ng*nel, nel_, n, onec, h2phi_x(:,i1:i2), n, coeff_bi(:,j1:j2), n, &
                  zeroc, help(1:ng*nel,1:nel_), ng*nel)
        if (profile_code) call end_profiling("gemm_en", ng*nel, nel_, n, "c")

        if (profile_code) call start_profiling("alpha_g_ij")
        do w = w1+1, w2
          k1 = (w-w1-1)*nel + 1
          k2 = (w-w1)*nel 
          do g = 1, ng
            g1 = nel * (g-1) + 1
            g2 = nel * g
            !be careful with e_x_g = \sum_ij \alpha_g_ij \alpha_g_ji
            ex(g,w) = ex(g,w) - 0.5_wp * hdes%rspin * trace2(help(g1:g2,k1:k2), help(g1:g2,k1:k2))
          end do
        end do
        if (profile_code) call end_profiling("alpha_g_ij", 8*nel**2*ng*nw_)
      end do
      !$omp end parallel do
    else 
      if (.not. allocated(help)) allocate(help(ng*maxval(hdes%nel), nw*maxval(hdes%nel)))
      help = zeroc

      nel_ = nel * nw
      j1 = hdes%nell(spin)*nw + 1
      j2 = hdes%nell(spin)*nw + nel_

      if (profile_code) call start_profiling("gemm_en")
!debug: cmplx nosd - careful analysis whether we need "t", "n" or "c"
      call gemm("t", "n", ng*nel, nel_, n, onec, h2phi_x(:,i1:i2), n, coeff_bi(:,j1:j2), n, &
                zeroc, help(1:ng*nel,1:nel_), ng*nel)
      if (profile_code) call end_profiling("gemm_en", ng*nel, nel_, n, "c")

      if (profile_code) call start_profiling("alpha_g_ij")
      do w = 1, nw
        k1 = (w-1)*nel + 1
        k2 = w*nel 
        do g = 1, ng
          g1 = nel * (g-1) + 1
          g2 = nel * g
          !be careful with e_x_g = \sum_ij \alpha_g_ij \alpha_g_ji
          ex(g,w) = ex(g,w) - 0.5_wp * hdes%rspin * trace2(help(g1:g2,k1:k2), help(g1:g2,k1:k2))
        end do
      end do
      if (profile_code) call end_profiling("alpha_g_ij", 8*nel**2*ng*nw)
    end if 
  end do

  if (profile_code) call end_profiling(proc_name)
end subroutine exchange_energy_no_chol_batch_c


!******************************************************************************** 
!
! Calculates exchange energy for nw set of orbitals (using half-rotated integrals)
!
!    E_x(w) = -0.5 \sum_{IJ} coeff_{Iw} h2phi_x_{IJ} coeff_{Jw}
! 
!    Input:
!        hdes    - Hamilton descriptor
!        h2phi_x - intermediate two-body array (see get_h2phi_occ_x_int)
!        coeff   - nw sets of orbitals (biorthogonalized orbitals)
!        bsize   - specifies how many walkers are treated simultaneously
!
!    Output:
!        ex - array of exchange energies
!
!******************************************************************************** 
subroutine exchange_energy_no_int_batch_r(hdes, h2phi_x, coeff_bi, ex, bsize)
  type(hamil_descriptor), intent(in)    :: hdes
  real(wp), intent(in)                  :: h2phi_x(:,:,:)
  complex(wp), intent(in)               :: coeff_bi(:,:)
  complex(wp), allocatable, intent(out) :: ex(:)
  integer, optional, intent(in)         :: bsize
  !local variables
  integer                               :: ispin, n, nocc, nel, nel_, ng, nw, nw_, bsize_
  integer                               :: spin, w, ww, w1, w2, j1, j2, xblock
  complex(wp), allocatable              :: coeff_(:,:), help(:,:)
  character(len=*), parameter           :: proc_name = "exchange_energy_no_eri"

  if (profile_code) call start_profiling(proc_name)

  n = hdes%n
  ispin = hdes%ispin
  nocc = hdes%nocc_
  nw = size(coeff_bi, 2) / nocc

  if (present(bsize)) then
    bsize_ = bsize
  else
    bsize_ = nw
  end if 

  allocate(ex(nw))
  ex = zeroc

  do spin = 1, ispin
    nel = hdes%nel(spin) 
    if (nel == 0) cycle

    j1 = hdes%nell(spin)*nw + 1
    j2 = hdes%nell(spin)*nw + hdes%nel(spin)*nw
!debug: opt - try to get rid of that reshaping
    if (allocated(coeff_)) deallocate(coeff_)
    allocate(coeff_(n*nel,nw))
    coeff_ = reshape(coeff_bi(:,j1:j2), shape=[n*nel,nw])

    if (is_blocking(nw, bsize_)) then
      if (.not. allocated(help)) allocate(help(n*maxval(hdes%nel), bsize_))
      help = zeroc

      !$omp parallel do private(w, ww, w1, w2, nw_, nel_, j1, j2) firstprivate(help)
      do xblock = 1, get_nblocks(nw, bsize_)
        w1 = (xblock-1) * bsize_
        w2 = min(xblock*bsize_, nw)
        nw_ = w2 - w1

        j1 = w1 + 1
        j2 = w2

        if (profile_code) call start_profiling("gemm_en")
        call gemm("n", "n", n*nel, nw_, n*nel, onec, h2phi_x(1:n*nel,1:n*nel,spin), n*nel, coeff_(:,j1:j2), n*nel, &
                  zeroc, help(1:n*nel,1:nw_), n*nel)
        if (profile_code) call end_profiling("gemm_en", n*nel, n*nel, nw_, "rc")

        if (profile_code) call start_profiling("theta_I")
        do w = w1+1, w2
          ww = w - w1
          ex(w) = ex(w) - 0.5_wp * hdes%rspin * sum(coeff_(1:n*nel,w)*help(1:n*nel,ww))
        end do
        if (profile_code) call end_profiling("theta_I")
      end do
      !$omp end parallel do
    else 
      if (.not. allocated(help)) allocate(help(n*maxval(hdes%nel), nw))
      help = zeroc

      if (profile_code) call start_profiling("gemm_en")
      call gemm("n", "n", n*nel, nw, n*nel, onec, h2phi_x(1:n*nel,1:n*nel,spin), n*nel, coeff_, n*nel, &
                zeroc, help(1:n*nel,1:nw), n*nel)
      if (profile_code) call end_profiling("gemm_en", n*nel, n*nel, nw, "rc")

      if (profile_code) call start_profiling("theta_I")
      do w = 1, nw
        ex(w) = ex(w) - 0.5_wp * hdes%rspin * sum(coeff_(1:n*nel,w)*help(1:n*nel,w))
      end do
      if (profile_code) call end_profiling("theta_I")
    end if 
  end do

  if (profile_code) call end_profiling(proc_name)
end subroutine exchange_energy_no_int_batch_r

subroutine exchange_energy_no_int_batch_c(hdes, h2phi_x, coeff_bi, ex, bsize)
  type(hamil_descriptor), intent(in)    :: hdes
  complex(wp), intent(in)               :: h2phi_x(:,:,:)
  complex(wp), intent(in)               :: coeff_bi(:,:)
  complex(wp), allocatable, intent(out) :: ex(:)
  integer, optional, intent(in)         :: bsize
  !local variables
  integer                               :: ispin, n, nocc, nel, nel_, ng, nw, nw_, bsize_
  integer                               :: spin, w, ww, w1, w2, j1, j2, xblock
  complex(wp), allocatable              :: coeff_(:,:), help(:,:)
  character(len=*), parameter           :: proc_name = "exchange_energy_no_eri"

  if (profile_code) call start_profiling(proc_name)

  n = hdes%n
  ispin = hdes%ispin
  nocc = hdes%nocc_
  nw = size(coeff_bi, 2) / nocc

  if (present(bsize)) then
    bsize_ = bsize
  else
    bsize_ = nw
  end if 

  allocate(ex(nw))
  ex = zeroc

  do spin = 1, ispin
    nel = hdes%nel(spin) 
    if (nel == 0) cycle

    j1 = hdes%nell(spin)*nw + 1
    j2 = hdes%nell(spin)*nw + hdes%nel(spin)*nw
!debug: opt - try to get rid of that reshaping
    if (allocated(coeff_)) deallocate(coeff_)
    allocate(coeff_(n*nel,nw))
    coeff_ = reshape(coeff_bi(:,j1:j2), shape=[n*nel,nw])

    if (is_blocking(nw, bsize_)) then
      if (.not. allocated(help)) allocate(help(n*maxval(hdes%nel), bsize_))
      help = zeroc

      !$omp parallel do private(w, ww, w1, w2, nw_, nel_, j1, j2) firstprivate(help)
      do xblock = 1, get_nblocks(nw, bsize_)
        w1 = (xblock-1) * bsize_
        w2 = min(xblock*bsize_, nw)
        nw_ = w2 - w1

        j1 = w1 + 1
        j2 = w2

        if (profile_code) call start_profiling("gemm_en")
        call gemm("n", "n", n*nel, nw_, n*nel, onec, h2phi_x(1:n*nel,1:n*nel,spin), n*nel, coeff_(:,j1:j2), n*nel, &
                  zeroc, help(1:n*nel,1:nw_), n*nel)
        if (profile_code) call end_profiling("gemm_en", n*nel, n*nel, nw_, "c")

        if (profile_code) call start_profiling("theta_I")
        do w = w1+1, w2
          ww = w - w1
          ex(w) = ex(w) - 0.5_wp * hdes%rspin * sum(coeff_(1:n*nel,w)*help(1:n*nel,ww))
        end do
        if (profile_code) call end_profiling("theta_I")
      end do
      !$omp end parallel do
    else 
      if (.not. allocated(help)) allocate(help(n*maxval(hdes%nel), nw))
      help = zeroc

      if (profile_code) call start_profiling("gemm_en")
      call gemm("n", "n", n*nel, nw, n*nel, onec, h2phi_x(1:n*nel,1:n*nel,spin), n*nel, coeff_, n*nel, &
                zeroc, help(1:n*nel,1:nw), n*nel)
      if (profile_code) call end_profiling("gemm_en", n*nel, n*nel, nw, "c")

      if (profile_code) call start_profiling("theta_I")
      do w = 1, nw
        ex(w) = ex(w) - 0.5_wp * hdes%rspin * sum(coeff_(1:n*nel,w)*help(1:n*nel,w))
      end do
      if (profile_code) call end_profiling("theta_I")
    end if 
  end do

  if (profile_code) call end_profiling(proc_name)
end subroutine exchange_energy_no_int_batch_c


!******************************************************************************** 
!
! Extracts occupied manifold in Slater determinant from all orbitals
! 
!******************************************************************************** 
subroutine get_coeff_occ_r(hdes, coeff, coeff_occ)
  type(hamil_descriptor), intent(in) :: hdes
  real(wp), intent(in)               :: coeff(:,:,:)
  real(wp), allocatable, intent(out) :: coeff_occ(:,:)
  !local variables
  integer                            :: spin, sp, nel, w1, w2
  
  allocate(coeff_occ(hdes%n, hdes%nocc))

  do spin = 1, hdes%ispin
    sp = hdes%get_spin(spin) 
    nel = hdes%nel(spin)
    w1 = hdes%nell(spin) + 1
    w2 = hdes%nell(spin) + nel
    coeff_occ(:,w1:w2) = coeff(:,1:nel,sp)
  end do
end subroutine get_coeff_occ_r

subroutine get_coeff_occ_c(hdes, coeff, coeff_occ)
  type(hamil_descriptor), intent(in)    :: hdes
  complex(wp), intent(in)               :: coeff(:,:,:)
  complex(wp), allocatable, intent(out) :: coeff_occ(:,:)
  !local variables
  integer                               :: spin, sp, nel, w1, w2
  
  allocate(coeff_occ(hdes%n, hdes%nocc))

  do spin = 1, hdes%ispin
    sp = hdes%get_spin(spin) 
    nel = hdes%nel(spin)
    w1 = hdes%nell(spin) + 1
    w2 = hdes%nell(spin) + nel
    coeff_occ(:,w1:w2) = coeff(:,1:nel,sp)
  end do
end subroutine get_coeff_occ_c


!******************************************************************************** 
!
! Applies orbitals from left SD on the one-body hamiltonian
!
!    Input:
!        phi - set of orbitals
!        h   - one-body Hamiltonian
!
!    Output:
!        hphi = h phi            for real h, phi
!        hphi = h phi^{*}        for real h, cmplx phi
!        hphi = h^{*} phi^{*}    for cmplx h, phi
!
!******************************************************************************** 
subroutine get_h1phi_r(phi, h, hphi)
  real(wp), intent(in)               :: phi(:,:,:)
  real(wp), intent(in)               :: h(:,:,:)
  real(wp), allocatable, intent(out) :: hphi(:,:,:)
  !local variables
  integer                            :: n, m, ispin, spin, spin2
  
  n = size(phi, 1)
  m = size(phi, 2)
  ispin = size(phi, 3)

  allocate(hphi(n,m,ispin))
  
  do spin = 1, ispin
    spin2 = spin
    if (size(h, 3) == 1) spin2 = 1
    call gemm("n", "n", n, m, n, oner, h(:,:,spin2), n, phi(:,:,spin), n, zeror, hphi(:,:,spin), n)
  end do
end subroutine get_h1phi_r

subroutine get_h1phi_cr(phi, h, hphi)
  complex(wp), intent(in)               :: phi(:,:,:)
  real(wp), intent(in)                  :: h(:,:,:)
  complex(wp), allocatable, intent(out) :: hphi(:,:,:)
  !local variables
  integer                               :: n, m, ispin, spin, spin1
  complex(wp), allocatable              :: phi_c(:,:,:)
  
  n = size(phi, 1)
  m = size(phi, 2)
  ispin = size(phi, 3)

  allocate(hphi(n,m,ispin))
  allocate(phi_c, mold=phi)
  phi_c = conjg(phi)
  
  do spin = 1, ispin
    spin1 = spin
    if (size(h, 3) == 1) spin1 = 1
    call gemm("n", "n", n, m, n, onec, h(:,:,spin1), n, phi_c(:,:,spin), n, zeroc, hphi(:,:,spin), n)
  end do
end subroutine get_h1phi_cr


!******************************************************************************** 
!
! Apply orbitals from left SD on the Cholesky vectors
!
!    Input:
!        phi     - set of orbitals
!        h2_gres - Cholesky vectors (nn,ng,ispin1)
!
!    Output:
!        h2phi_gres(:,:,g,:) = h2_gres(:,:,g,:) phi           for real h2_gres, phi
!        h2phi_gres(:,:,g,:) = h2_gres(:,:,g,:) phi^{*}       for real h2_gres, cmplx phi
!        h2phi_gres(:,:,g,:) = h2_gres(:,:,g,:)^{*} phi^{*}   for cmplx h2_gres, phi
!
!******************************************************************************** 
subroutine get_h2phi_r(phi, h2_gres, h2phi_gres)
  real(wp), intent(in)               :: phi(:,:,:)
  real(wp), target, intent(in)       :: h2_gres(:,:,:)
  real(wp), allocatable, intent(out) :: h2phi_gres(:,:,:,:)
  !local variables
  integer                            :: g, ng, n, m, ispin, ispin1, spin, spin1
  real(wp), pointer                  :: h2_gres_g(:,:)
  character(len=*), parameter        :: proc_name = "get_h2phi"
  character(len=*), parameter        :: gemm_name = "gemm_h2phi"

  if (use_profiling) call start_profiling(proc_name)
  
  n = size(phi, 1)
  m = size(phi, 2)
  ispin = size(phi, 3)
  ng = size(h2_gres, 2)
  ispin1 = size(h2_gres, 3)

  allocate(h2phi_gres(n,m,ng,ispin))
  
  do spin = 1, ispin
    spin1 = spin
    if (ispin1 == 1) spin1 = 1
    do g = 1, ng
      h2_gres_g(1:n,1:n) => h2_gres(:,g,spin1)
      if (profile_code) call start_profiling(gemm_name)
      call gemm("n", "n", n, m, n, oner, h2_gres_g, n, phi(:,:,spin), n, zeror, h2phi_gres(:,:,g,spin), n)
      if (profile_code) call end_profiling(gemm_name, n, m, n, "r")
    end do
  end do

  if (use_profiling) call end_profiling(proc_name)
end subroutine get_h2phi_r

subroutine get_h2phi_cr(phi, h2_gres, h2phi_gres)
  complex(wp), intent(in)               :: phi(:,:,:)
  real(wp), target, intent(in)          :: h2_gres(:,:,:)
  complex(wp), allocatable, intent(out) :: h2phi_gres(:,:,:,:)
  !local variables
  integer                               :: g, ng, n, m, ispin, ispin1, spin, spin1
  complex(wp), allocatable              :: phi_c(:,:,:)
  real(wp), pointer                     :: h2_gres_g(:,:)
  character(len=*), parameter           :: proc_name = "get_h2phi"
  character(len=*), parameter           :: gemm_name = "gemm_h2phi"

  if (use_profiling) call start_profiling(proc_name)
  
  n = size(phi, 1)
  m = size(phi, 2)
  ispin = size(phi, 3)
  ng = size(h2_gres, 2)
  ispin1 = size(h2_gres, 3)

  allocate(h2phi_gres(n,m,ng,ispin))
  allocate(phi_c, mold=phi)
  phi_c = conjg(phi)
  
  do spin = 1, ispin
    spin1 = spin
    if (ispin1 == 1) spin1 = 1
    do g = 1, ng
      h2_gres_g(1:n,1:n) => h2_gres(:,g,spin1)
      if (profile_code) call start_profiling(gemm_name)
      call gemm("n", "n", n, m, n, onec, h2_gres_g, n, phi_c(:,:,spin), n, zeroc, h2phi_gres(:,:,g,spin), n)
      if (profile_code) call end_profiling(gemm_name, n, m, n, "rc")
    end do
  end do

  if (use_profiling) call end_profiling(proc_name)
end subroutine get_h2phi_cr


!******************************************************************************** 
!
! Apply orbitals from left SD on the g-resolve self energy 
!
!    Input:
!        phi     - set of orbitals
!        h2_gres - Cholesky vectors (nn,ng,ispin1)
!
!    Output:
!        hself_gres(:,g,:) = phi^{+} h2_gres^2(:,g,:)
!
!******************************************************************************** 
subroutine get_hselfphi_r(phi, h2_gres, hselfphi_gres)
  real(wp), intent(in)               :: phi(:,:,:)
  real(wp), target, intent(in)       :: h2_gres(:,:,:)
  real(wp), allocatable, intent(out) :: hselfphi_gres(:,:,:,:)
  !local variables
  integer                            :: g, ng, n, m, spin, spin1, ispin, ispin1
  real(wp), allocatable              :: hself_g(:,:)
  real(wp), pointer                  :: h2_gres_g(:,:)
  character(len=*), parameter        :: proc_name = "get_hselfphi"
  character(len=*), parameter        :: gemm_name = "gemm_hselfphi", gemm2_name = "gemm_hself"
  
  if (use_profiling) call start_profiling(proc_name)

  n = size(phi, 1)
  m = size(phi, 2)
  ispin = size(phi, 3)
  ng = size(h2_gres, 2)
  ispin1 = size(h2_gres, 3)

  allocate(hselfphi_gres(n, m, ng, ispin))
  allocate(hself_g(n,n))
  
  do spin = 1, ispin
    spin1 = spin
    if (ispin1 == 1) spin1 = 1
    do g = 1, ng
      h2_gres_g(1:n,1:n) => h2_gres(:,g,spin1)
      if (profile_code) call start_profiling(gemm2_name)
      call gemm("n", "n", n, n, n, -0.5_wp, h2_gres_g, n, h2_gres_g, n, zeror, hself_g, n)
      if (profile_code) call end_profiling(gemm2_name, n, n, n, "r")
      if (profile_code) call start_profiling(gemm_name)
      call gemm("n", "n", n, m, n, oner, hself_g, n, phi(:,:,spin), n, zeror, hselfphi_gres(:,:,g,spin), n)
      if (profile_code) call end_profiling(gemm_name, n, m, n, "r")
    end do
  end do

  if (use_profiling) call end_profiling(proc_name)
end subroutine get_hselfphi_r

subroutine get_hselfphi_cr(phi, h2_gres, hselfphi_gres)
  complex(wp), intent(in)               :: phi(:,:,:)
  real(wp), target, intent(in)          :: h2_gres(:,:,:)
  complex(wp), allocatable, intent(out) :: hselfphi_gres(:,:,:,:)
  !local variables
  integer                               :: g, ng, n, m, spin, spin1, ispin, ispin1
  real(wp), allocatable                 :: hself_g(:,:)
  real(wp), pointer                     :: h2_gres_g(:,:)
  complex(wp), allocatable              :: phi_c(:,:,:)
  character(len=*), parameter           :: proc_name = "get_hselfphi"
  character(len=*), parameter           :: gemm_name = "gemm_hselfphi", gemm2_name = "gemm_hself"
  
  if (use_profiling) call start_profiling(proc_name)
  
  n = size(phi, 1)
  m = size(phi, 2)
  ispin = size(phi, 3)
  ng = size(h2_gres, 2)
  ispin1 = size(h2_gres, 3)

  allocate(hselfphi_gres(n, m, ng, ispin))
  allocate(hself_g(n,n))
  allocate(phi_c, mold=phi)
  phi_c = conjg(phi)
  
  do spin = 1, ispin
    spin1 = spin
    if (ispin1 == 1) spin1 = 1
    do g = 1, ng
      h2_gres_g(1:n,1:n) => h2_gres(:,g,spin1)
      if (profile_code) call start_profiling(gemm2_name)
      call gemm("n", "n", n, n, n, -0.5_wp, h2_gres_g, n, h2_gres_g, n, zeror, hself_g, n)
      if (profile_code) call end_profiling(gemm2_name, n, n, n, "r")
      if (profile_code) call start_profiling(gemm_name)
      call gemm("n", "n", n, m, n, onec, hself_g, n, phi_c(:,:,spin), n, zeroc, hselfphi_gres(:,:,g,spin), n)
      if (profile_code) call end_profiling(gemm_name, n, m, n, "rc")
    end do
  end do

  if (use_profiling) call end_profiling(proc_name)
end subroutine get_hselfphi_cr


!******************************************************************************** 
!
! Extract h1phi_occ from  h1phi
!
! Both spin channels stored explicitly due to the simplicity
!
!    Input:
!        hdes  - Hamiltonian descriptor that defines occupied space
!        h1phi(n,m,ispin)
!
!    Output:
!        h1phi_occ(n,ne)
!
!    -    -- -- -- -- -- -- --
!    |   |  |  |  |  |  |  |  |
!    |   |  |  |  |  |  |  |  |
!    |   |  |  |  |  |  |  |  |
!    n   |  |  |  |  |  |  |  |
!    |   |  |  |  |  |  |  |  |
!    |   |  |  |  |  |  |  |  |
!    -    -- -- -- -- -- -- --
!        <---ne_up---><-ne_dn->
!
! Effectively, extracts the occupied manifold from h1phi
!
!******************************************************************************** 
subroutine get_h1phi_occ_r(hdes, h1phi, h1phi_occ)
  type(hamil_descriptor), intent(in) :: hdes
  real(wp), intent(in)               :: h1phi(:,:,:)
  real(wp), allocatable, intent(out) :: h1phi_occ(:,:)
  !local variables
  integer                            :: n, nocc, nel
  integer                            :: spin, sp, w1, w2
  character(len=*), parameter        :: proc_name = "get_h1phi_occ"

  if (use_profiling) call start_profiling(proc_name)
  
  n = size(h1phi, 1)
  nocc = sum(hdes%nel_space)
  
  allocate(h1phi_occ(n, nocc))

  do spin = 1, size(hdes%nel_space)
    sp = hdes%get_spin(spin) 
    nel = hdes%nel_space(spin)
    if (nel == 0) cycle
    w1 = hdes%nell_space(spin) + 1
    w2 = hdes%nell_space(spin) + nel
    h1phi_occ(:,w1:w2) = h1phi(:,1:nel,sp)
  end do

  if (use_profiling) call end_profiling(proc_name)
end subroutine get_h1phi_occ_r

subroutine get_h1phi_occ_c(hdes, h1phi, h1phi_occ)
  type(hamil_descriptor), intent(in)    :: hdes
  complex(wp), intent(in)               :: h1phi(:,:,:)
  complex(wp), allocatable, intent(out) :: h1phi_occ(:,:)
  !local variables
  integer                               :: n, nocc, nel
  integer                               :: spin, sp, w1, w2
  character(len=*), parameter           :: proc_name = "get_h1phi_occ"

  if (use_profiling) call start_profiling(proc_name)
  
  n = size(h1phi, 1)
  nocc = sum(hdes%nel_space)
  
  allocate(h1phi_occ(n, nocc))

  do spin = 1, size(hdes%nel_space)
    sp = hdes%get_spin(spin) 
    nel = hdes%nel_space(spin)
    if (nel == 0) cycle
    w1 = hdes%nell_space(spin) + 1
    w2 = hdes%nell_space(spin) + nel
    h1phi_occ(:,w1:w2) = h1phi(:,1:nel,sp)
  end do

  if (use_profiling) call end_profiling(proc_name)
end subroutine get_h1phi_occ_c


!******************************************************************************** 
!
! Extract h2phi_h from h2phi (for Hartree energy evaluation)                    
!
! Both spin channels stored explicitly due to the simplicity
!
!               ---------
!               | |     |
!         h2phi(n,m,ng,ispin)  ---> h2phi_h(ng,n*nocc)
!
!        |s=1 ne=1|s=1 ne=2|s=1 ne=3|s=2 ne=1|s=2 ne=2
!    -    -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!    |   |  !  !  |  !  !  |  !  !  |  !  !  |  !  !  |
!    |   |  !  !  |  !  !  |  !  !  |  !  !  |  !  !  |
!    |   |  !  !  |  !  !  |  !  !  |  !  !  |  !  !  |
!    ng  |  !  !  |  !  !  |  !  !  |  !  !  |  !  !  |
!    |   |  !  !  |  !  !  |  !  !  |  !  !  |  !  !  |
!    |   |  !  !  |  !  !  |  !  !  |  !  !  |  !  !  |
!    -    -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!        |---n----|
!         
!******************************************************************************** 
subroutine get_h2phi_h_occ_r(hdes, h2phi, h2phi_h)
  type(hamil_descriptor), intent(in) :: hdes
  real(wp), intent(in)               :: h2phi(:,:,:,:)
  real(wp), allocatable, intent(out) :: h2phi_h(:,:)
  !local variables
  integer                            :: n, ng, nel, nocc, ispin
  integer                            :: spin, sp, off, i, p, ii 
  character(len=*), parameter        :: proc_name = "get_h2phi_h_occ"

  if (use_profiling) call start_profiling(proc_name)

  n = size(h2phi, 1)
  ng = size(h2phi, 3)
  nocc = sum(hdes%nel_space)

  allocate(h2phi_h(ng, n*nocc))
  h2phi_h = 0.0_wp

  do spin = 1, size(hdes%nel_space)
    sp = hdes%get_spin(spin)
    nel = hdes%nel_space(spin)
    off = hdes%nell_space(spin)
    do i = 1, nel
      do p = 1, n
        ii = (off+i-1)*n + p
        h2phi_h(:,ii) = h2phi(p,i,:,sp)
      end do
    end do
  end do

  if (use_profiling) call end_profiling(proc_name)
end subroutine get_h2phi_h_occ_r

subroutine get_h2phi_h_occ_c(hdes, h2phi, h2phi_h)
  type(hamil_descriptor), intent(in)    :: hdes
  complex(wp), intent(in)               :: h2phi(:,:,:,:)
  complex(wp), allocatable, intent(out) :: h2phi_h(:,:)
  !local variables
  integer                               :: n, ng, nel, nocc, ispin
  integer                               :: spin, sp, off, i, p, ii 
  character(len=*), parameter           :: proc_name = "get_h2phi_h_occ"

  if (use_profiling) call start_profiling(proc_name)

  n = size(h2phi, 1)
  ng = size(h2phi, 3)
  nocc = sum(hdes%nel_space)

  allocate(h2phi_h(ng, n*nocc))
  h2phi_h = 0.0_wp

  do spin = 1, size(hdes%nel_space)
    sp = hdes%get_spin(spin)
    nel = hdes%nel_space(spin)
    off = hdes%nell_space(spin)
    do i = 1, nel
      do p = 1, n
        ii = (off+i-1)*n + p
        h2phi_h(:,ii) = h2phi(p,i,:,sp)
      end do
    end do
  end do

  if (use_profiling) call end_profiling(proc_name)
end subroutine get_h2phi_h_occ_c


!******************************************************************************** 
!
! Extract h2phi_x from h2phi (for exchange energy evaluation)                    
!
! Both spin channels stored explicitly due to the simplicity
!
!               ---------
!               |   |   |
!         h2phi(n,m,ng,ispin)  ---> h2phi_x(n,ng*nocc)
!
!        |s=1  g=1|s=1  g=2|s=1 g=ng|s2 g2|s2 g2|s2 ng|
!    -    -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!    |   |  !  !  |  !  !  |  !  !  |  !  |  !  |  !  |
!    |   |  !  !  |  !  !  |  !  !  |  !  |  !  |  !  |
!    |   |  !  !  |  !  !  |  !  !  |  !  |  !  |  !  |
!    n   |  !  !  |  !  !  |  !  !  |  !  |  !  |  !  |
!    |   |  !  !  |  !  !  |  !  !  |  !  |  !  |  !  |
!    |   |  !  !  |  !  !  |  !  !  |  !  |  !  |  !  |
!    -    -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!        |-ne_up--|                 |ne_dn|
!         
!******************************************************************************** 
subroutine get_h2phi_x_occ_chol_r(hdes, h2phi, h2phi_x, compress, tol)
  use cholesky, only: iterative_cholesky_compress
  type(hamil_descriptor), intent(in) :: hdes
  real(wp), intent(in)               :: h2phi(:,:,:,:)
  real(wp), allocatable, intent(out) :: h2phi_x(:,:)
  logical, intent(in)                :: compress
  real(wp), intent(in)               :: tol
  !local variables
  integer                            :: n, ng, ispin, nocc, noccmax, nel
  integer                            :: spin, sp, w1, w2, g
  real(wp), allocatable              :: h2phi_occ(:,:,:,:), h2phi_occ_(:,:,:)
  character(len=*), parameter        :: proc_name = "get_h2phi_x_occ_chol"

  if (use_profiling) call start_profiling(proc_name)
  
  n = size(h2phi, 1)
  ispin = size(h2phi, 4)
  nocc = sum(hdes%nel_space)
  noccmax = maxval(hdes%nel_space)
  
  allocate(h2phi_occ(n, noccmax, size(h2phi,3), ispin))
  h2phi_occ = h2phi(:,1:noccmax,:,:)

!debug:
!I am not sure that this line is a good idea
!  if (ispin == 2) h2phi_occ(:,hdes%nel_space(2)+1:,:,2) = zeror
  
  if (compress) then
    call get_h2_gres_2d(h2phi_occ, h2phi_occ_, .true.)
    call iterative_cholesky_compress(h2phi_occ_, tol, .false.)
    call get_h2_gres_3d(h2phi_occ_, h2phi_occ, n, noccmax, .true.)
  end if
  
  ng = size(h2phi_occ, 3)
  
  allocate(h2phi_x(n,ng*nocc))

  do spin = 1, size(hdes%nel_space)
    sp = hdes%get_spin(spin)
    nel = hdes%nel_space(spin)
    if (nel == 0) cycle
    do g = 1, ng
      w1 = hdes%nell_space(spin)*ng + (g-1)*nel + 1
      w2 = hdes%nell_space(spin)*ng + g*nel
      h2phi_x(:,w1:w2) = h2phi_occ(:,1:nel,g,sp)
    end do
  end do

!debug: optimization : simpler but did not show speedup at the moment
  !!do spin = 1, size(hdes%nel_space)
  !!  sp = hdes%get_spin(spin)
  !!  nel = hdes%nel_space(spin)
  !!  if (nel == 0) cycle
  !!  w1 = hdes%nell_space(spin)*ng + 1
  !!  w2 = (hdes%nell_space(spin)+nel) * ng
  !!  h2phi_x(:,w1:w2) = reshape(h2phi_occ(:,1:nel,:,sp), shape=[n,ng*nel])
  !!end do

  if (use_profiling) call end_profiling(proc_name)
end subroutine get_h2phi_x_occ_chol_r

subroutine get_h2phi_x_occ_chol_c(hdes, h2phi, h2phi_x, compress, tol)
  use cholesky, only: iterative_cholesky_compress
  type(hamil_descriptor), intent(in)    :: hdes
  complex(wp), intent(in)               :: h2phi(:,:,:,:)
  complex(wp), allocatable, intent(out) :: h2phi_x(:,:)
  logical, intent(in)                   :: compress
  real(wp), intent(in)                  :: tol
  !local variables
  integer                               :: n, ng, ispin, nocc, noccmax, nel
  integer                               :: spin, sp, w1, w2, g
  complex(wp), allocatable              :: h2phi_occ(:,:,:,:), h2phi_occ_(:,:,:)
  character(len=*), parameter           :: proc_name = "get_h2phi_x_occ_chol"

  if (use_profiling) call start_profiling(proc_name)
  
  n = size(h2phi, 1)
  ispin = size(h2phi, 4)
  nocc = sum(hdes%nel_space)
  noccmax = maxval(hdes%nel_space)
  
!debug: optimization : as long as the compression is note used,
!h2phi_occ is not necessary
  !!allocate(h2phi_occ(n, noccmax, size(h2phi,3), ispin))
  !!h2phi_occ = h2phi(:,1:noccmax,:,:)

!debug:
!I am not sure that this line is a good idea
!  if (ispin == 2) h2phi_occ(:,hdes%nel_space(2)+1:,:,2) = zeroc
  
!debug: cmplx nosd - this does not work for complex-valued Cholesky vectors at the moment
  !if (compress) then
  !  call get_h2_gres_2d(h2phi_occ, h2phi_occ_, .true.)
  !  call iterative_cholesky_compress(h2phi_occ_, tol, .false.)
  !  call get_h2_gres_3d(h2phi_occ_, h2phi_occ, n, noccmax, .true.)
  !end if
  
!debug: optimization
  ng = size(h2phi, 3)
  !!ng = size(h2phi_occ, 3)

  allocate(h2phi_x(n,ng*nocc))

  do spin = 1, size(hdes%nel_space)
    sp = hdes%get_spin(spin)
    nel = hdes%nel_space(spin)
    if (nel == 0) cycle
    do g = 1, ng
      w1 = hdes%nell_space(spin)*ng + (g-1)*nel + 1
      w2 = hdes%nell_space(spin)*ng + g*nel
      h2phi_x(:,w1:w2) = h2phi(:,1:nel,g,sp)
      !h2phi_x(:,w1:w2) = h2phi_occ(:,1:nel,g,sp)
    end do
  end do

!debug: optimization : simpler but did not show speedup at the moment
  !!do spin = 1, size(hdes%nel_space)
  !!  sp = hdes%get_spin(spin)
  !!  nel = hdes%nel_space(spin)
  !!  if (nel == 0) cycle
  !!  w1 = hdes%nell_space(spin)*ng + 1
  !!  w2 = (hdes%nell_space(spin)+nel) * ng
  !!  h2phi_x(:,w1:w2) = reshape(h2phi_occ(:,1:nel,:,sp), shape=[n,ng*nel])
  !!end do

  if (use_profiling) call end_profiling(proc_name)
end subroutine get_h2phi_x_occ_chol_c


!******************************************************************************** 
!
! Precalculate array for the exchange energy evaluation using ERIs
!
!    h2phi_x{I[pj]J[qi]}  = [pi|qj] = V_{II[pi]JJ[qj]}
!
! h2phi_x tensor has the following shape: h2phi_x(n x nel_space x ispin)
! Exchange energy is evaluated as:
!
!    E_x = \sum_{IJ} \theta_I h2phi_x{IJ} \theta_J
!
! Note: After Coulomb integrals V_{II[pi]JJ[qj]} are constructed, indices have to
!       be reordered (II, JJ) --> (I, J)
!
!******************************************************************************** 
subroutine get_h2phi_x_occ_int_r(hdes, h2phi, h2phi_x)
  type(hamil_descriptor), intent(in) :: hdes
  real(wp), intent(in)               :: h2phi(:,:,:,:)
  real(wp), allocatable, intent(out) :: h2phi_x(:,:,:)
  !local variables
  real(wp), allocatable              :: temp(:,:,:,:), temp_(:,:,:)
  integer                            :: n, ng, ispin, nel, noccmax
  integer                            :: i, j, ii, jj, iii, jjj, p, q, spin
  character(len=*), parameter        :: proc_name = "get_h2phi_x_occ_int"

  if (use_profiling) call start_profiling(proc_name)

  n = size(h2phi, 1)
  ng = size(h2phi, 3)
  ispin = size(h2phi, 4)
  noccmax = maxval(hdes%nel_space)

  allocate(temp(n, noccmax, ng, ispin))
  temp = h2phi(:,1:noccmax,:,:)
  if (ispin == 2) temp(:,hdes%nel_space(2)+1:,:,2) = zeror
  call get_h2_gres_2d(temp, temp_, .true.)
  call get_h2_full(temp_, ispin)

  allocate(h2phi_x(n*noccmax,n*noccmax,ispin))
  h2phi_x = zeror

  do spin = 1, ispin
    nel = hdes%nel_space(spin)
    do jj = 1, n*nel
      p = 1 + mod(jj-1, n)
      i = 1 + (jj-1) / n
      do ii = 1, n*nel
        q = 1 + mod(ii-1, n)
        j = 1 + (ii-1) / n

        iii = (j-1)*n + p
        jjj = (i-1)*n + q

        h2phi_x(ii,jj,spin) = temp_(iii,jjj,spin)
      end do
    end do
  end do

  if (use_profiling) call end_profiling(proc_name)
end subroutine get_h2phi_x_occ_int_r

subroutine get_h2phi_x_occ_int_c(hdes, h2phi, h2phi_x)
  type(hamil_descriptor), intent(in)    :: hdes
  complex(wp), intent(in)               :: h2phi(:,:,:,:)
  complex(wp), allocatable, intent(out) :: h2phi_x(:,:,:)
  !local variables
  complex(wp), allocatable              :: temp(:,:,:,:), temp_(:,:,:)
  integer                               :: n, ng, ispin, nel, noccmax
  integer                               :: i, j, ii, jj, iii, jjj, p, q, spin
  character(len=*), parameter           :: proc_name = "get_h2phi_x_occ_int"

  if (use_profiling) call start_profiling(proc_name)

  n = size(h2phi, 1)
  ng = size(h2phi, 3)
  ispin = size(h2phi, 4)
  noccmax = maxval(hdes%nel_space)

  allocate(temp(n, noccmax, ng, ispin))
  temp = zeroc
  temp = h2phi(:,1:noccmax,:,:)
  if (ispin == 2) temp(:,hdes%nel_space(2)+1:,:,2) = zeroc
  call get_h2_gres_2d(temp, temp_, .true.)
  call get_h2_full(temp_, ispin)

  allocate(h2phi_x(n*noccmax,n*noccmax,ispin))
  h2phi_x = zeroc

  do spin = 1, ispin
    nel = hdes%nel_space(spin)
    do jj = 1, n*nel
      p = 1 + mod(jj-1, n)
      i = 1 + (jj-1) / n
      do ii = 1, n*nel
        q = 1 + mod(ii-1, n)
        j = 1 + (ii-1) / n

        iii = (j-1)*n + p
        jjj = (i-1)*n + q

        h2phi_x(ii,jj,spin) = temp_(iii,jjj,spin)
      end do
    end do
  end do

  if (use_profiling) call end_profiling(proc_name)
end subroutine get_h2phi_x_occ_int_c


!******************************************************************************** 
!
! Determine whether omp is activated
!
!******************************************************************************** 
logical function is_omp(size_, factor)
  integer, intent(in) :: size_
  integer, intent(in) :: factor

  is_omp =  comm_world%ompsize*factor <= size_
end function is_omp


!******************************************************************************** 
!
! Determine whether blocking is activated
!
!******************************************************************************** 
logical function is_blocking(size_, block_size)
  integer, intent(in)         :: size_
  integer, intent(in)         :: block_size

!debug:
  !is_blocking =  (comm_world%ompsize-1)*block_size <= size_
  is_blocking =  comm_world%ompsize*block_size <= size_
end function is_blocking

end module nosd