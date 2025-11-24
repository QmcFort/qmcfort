! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module cholesky

#include "preproc.inc"
use constants
use mpi, only: comm_world
use profiling
use lapack
use hamilton_vars, only: hamil_descriptor
use hamilton_layout
use qmcfort_io
use gauss_base
use gauss_integrals

implicit none

interface compress_h2_chol
  module procedure compress_h2_chol_3d, compress_h2_chol_5d, compress_h2_chol_gauss
end interface compress_h2_chol

interface compress_h2_eig
  module procedure compress_h2_eig_3d, compress_h2_eig_5d
end interface compress_h2_eig

interface compress_h2_gres_chol 
  module procedure compress_h2_gres_chol_3d, compress_h2_gres_chol_4d
end interface compress_h2_gres_chol

interface compress_h2_gres_eig
  module procedure compress_h2_gres_eig_3d, compress_h2_gres_eig_4d
end interface compress_h2_gres_eig

interface h2_eig
  module procedure h2_eig_2d, h2_eig_3d
end interface h2_eig

interface h2_gres_eig
  module procedure h2_gres_eig_2d, h2_gres_eig_3d
end interface h2_gres_eig

interface iterative_cholesky
  module procedure iterative_cholesky_2d, iterative_cholesky_3d, iterative_cholesky_gauss
end interface iterative_cholesky

interface iterative_cholesky_compress
  module procedure iterative_cholesky_compress_2d, iterative_cholesky_compress_3d
end interface iterative_cholesky_compress

public

type cholesky_buffer
  integer               :: buffer_size
  integer               :: buffer_indx
  integer, allocatable  :: pivot_indx(:)
  real(wp), allocatable :: buffer(:,:)
  real(wp), allocatable :: chol_new(:,:)
contains
  procedure             :: init => init_chol_buffer
  procedure             :: is_full => is_chol_buffer_full
  procedure             :: append => append_chol_buffer
  procedure             :: add_to_h2 => add_chol_buffer_to_h2
end type cholesky_buffer

contains

!******************************************************************************** 
!
! Initialize buffer for sotrage of Cholesky vectors that are still not 
!    added to the main Cholesky vectors
!
!******************************************************************************** 
subroutine init_chol_buffer(self, size_, buffer_size)
  class(cholesky_buffer), intent(inout) :: self
  integer, intent(in)                   :: size_
  integer, intent(in)                   :: buffer_size

  self%buffer_size = buffer_size
  self%buffer_indx = 0
  allocate(self%buffer(size_, self%buffer_size))
  allocate(self%pivot_indx(self%buffer_size))
end subroutine init_chol_buffer


!******************************************************************************** 
!
! Chekc whether the Cholesky buffer is full
!
!******************************************************************************** 
logical function is_chol_buffer_full(self)
  class(cholesky_buffer), intent(in) :: self

  is_chol_buffer_full = self%buffer_indx == self%buffer_size
end function is_chol_buffer_full


!******************************************************************************** 
!
! Append Cholesky buffer:
!    Store new Coulomb row
!    Store pivot indx of that Coulomb row
!
!******************************************************************************** 
subroutine append_chol_buffer(self, coulomb_row, indx)
  class(cholesky_buffer), intent(inout) :: self
  real(wp), intent(in)                  :: coulomb_row(:)
  integer, intent(in)                   :: indx

  self%buffer_indx = self%buffer_indx + 1
  self%buffer(:,self%buffer_indx) = coulomb_row
  self%pivot_indx(self%buffer_indx) = indx
end subroutine append_chol_buffer


!******************************************************************************** 
!
! Empty Cholesky buffers: 
!
!    Reorthogonalize full buffer with respect to existing Cholesky vectors 
!
!******************************************************************************** 
subroutine add_chol_buffer_to_h2(self, h2, pivots, ng, chol_tol)
  class(cholesky_buffer), intent(inout) :: self
  real(wp), allocatable, intent(inout)  :: h2(:,:)
  real(wp), intent(inout)               :: pivots(:)
  integer, intent(inout)                :: ng
  real(wp), intent(in)                  :: chol_tol
  !local variables
  integer                              :: nbuff, size_, temp_ng, ng_, i, j, indx
  real(wp)                             :: res_val
  real(wp), allocatable                :: h2_v(:,:), h2_(:,:), h2_v_(:)
  logical, allocatable                 :: accept(:)
  character(len=*), parameter          :: proc_name = "add_to_h2"

  if (use_profiling) call start_profiling(proc_name)

#ifdef MKL
  call mkl_set_num_threads(48)
#endif

  nbuff = self%buffer_indx
  size_ = size(self%buffer, 1)
  temp_ng = size(h2, 2)

  allocate(h2_v(ng, nbuff))
  allocate(h2_(size_, nbuff))
  allocate(h2_v_(nbuff))
  allocate(accept(nbuff))

  do i = 1, self%buffer_indx
    indx = self%pivot_indx(i)
    h2_v(:,i) =  h2(indx,1:ng)
  end do

  if (profile_code) call start_profiling("chol_gemm")
  if (ng > 0) call gemm("n", "n", size_, nbuff, ng, -1.0_wp, h2(:,1:ng), size_, h2_v, ng, 1.0_wp, self%buffer, size_)
  if (profile_code) call end_profiling("chol_gemm")


  if (profile_code) call start_profiling("chol_rest")
  ng_ = 0
  do i = 1, nbuff
    indx = self%pivot_indx(i)

    if (ng_ > 0) then
      do j = 1, ng_
        h2_v_(j) = h2_(indx,j)
      end do
      call gemv("n", size_, ng_, -1.0_wp, h2_(:,1:ng_), size_, h2_v_(1:ng_), 1, 1.0_wp, self%buffer(:,i), 1)
    end if 

    res_val = self%buffer(indx,i)
    if (res_val < chol_tol) cycle       

    ng_ = ng_ + 1
    h2_(:,ng_) = self%buffer(:,i) / sqrt(res_val)
    pivots = pivots - h2_(:,ng_)*h2_(:,ng_)
    pivots(indx) = 0.0_wp
  end do
  if (profile_code) call end_profiling("chol_rest")

  if (ng+ng_ > temp_ng) then 
    temp_ng = temp_ng + max(100, ng_)
    call resize_h2_gres(temp_ng, h2)
  end if

  do i = 1, ng_
    h2(:,ng+i) = h2_(:,i)
  end do
  ng = ng + ng_

  self%buffer_indx = 0

#ifdef MKL
  call mkl_set_num_threads(1)
#endif

  if (use_profiling) call end_profiling(proc_name)
end subroutine add_chol_buffer_to_h2


!******************************************************************************** 
!
!       Driver for the compression of the Hamiltonian
!
!       If full Coulomb matrix elements are passed as input
!           Usual Cholesky decomposition is  done
! 
!       If auxiliary fields are passes as input
!           Compression of auxiliary fields is done
!           if Ng > N:
!               first calculate integrals then decompose
!           else:
!               direct compression of auxiliary fields
!
!******************************************************************************** 
subroutine compress_h2(hdes, h2, h2_gres)
  type(hamil_descriptor), intent(inout) :: hdes
  real(wp), allocatable, intent(inout)  :: h2(:,:,:,:,:)
  real(wp), allocatable, intent(inout)  :: h2_gres(:,:,:)
  
  if (allocated(h2)) then
    select case(trim(hdes%compress_mode))
      case ("cholesky", "chol")
        call compress_h2_chol(hdes, h2, h2_gres)
      case ("eig", "svd")
        call compress_h2_eig(hdes, h2, h2_gres)
      case default 
        call compress_h2_chol(hdes, h2, h2_gres)
    end select
  else if (allocated(h2_gres)) then
    call compress_h2_gres(hdes, h2_gres)
  else 
    write(*,*) "Problem  - h2 / h2_gres is not allocated"
  end if
end subroutine compress_h2


!******************************************************************************** 
!
!       Driver for the compression of the auxiliary fields
!
!******************************************************************************** 
subroutine compress_h2_gres(hdes, h2_gres)
  type(hamil_descriptor), intent(inout) :: hdes
  real(wp), allocatable, intent(inout)  :: h2_gres(:,:,:)
  !local variables
  real(wp), allocatable                 :: h2_spin(:,:)

  select case (trim(hdes%compress_mode))
    case ("cholesky", "chol")
      call compress_h2_gres_chol(hdes, h2_gres)
    case ("svd", "eig")
      call compress_h2_gres_eig(hdes, h2_gres)
  end select
end subroutine compress_h2_gres


!********************************************************************************
!
!       Modified Cholesky decompoisition for input h2(n*m,n*m,ispin2)
!
!********************************************************************************
subroutine compress_h2_chol_3d(hdes, h2, h2_gres)
  type(hamil_descriptor), intent(inout) :: hdes
  real(wp), allocatable, intent(inout)  :: h2(:,:,:)
  real(wp), allocatable, intent(inout)  :: h2_gres(:,:,:)
  !local variables
  real(wp), allocatable                 :: h2_spin(:,:)
  
  call move_alloc(h2, h2_gres)
  call iterative_cholesky(h2_gres, hdes%chol_tol, .true.)
  hdes%ng = size(h2_gres, 2)
  hdes%integral_mode = "cholesky"
end subroutine compress_h2_chol_3d

!******************************************************************************** 
!
!       Modified Cholesky decompoisition for input h2(n,m,n,m,ispin2)
!
!********************************************************************************
subroutine compress_h2_chol_5d(hdes, h2, h2_gres)
  type(hamil_descriptor), intent(inout) :: hdes
  real(wp), allocatable, intent(inout)  :: h2(:,:,:,:,:)
  real(wp), allocatable, intent(inout)  :: h2_gres(:,:,:)
  !local variables
  real(wp), allocatable                 :: h2_spin(:,:)
  
  call get_h2_fullspin(h2, h2_spin, dealloc=.true.)
  call iterative_cholesky(h2_spin, hdes%chol_tol, lverb=.true.)
  call get_h2_gres_from_fullspin(h2_gres, h2_spin, hdes%ispin1, dealloc=.true.)
  hdes%ng = size(h2_gres, 2)
  hdes%integral_mode = "cholesky"
end subroutine compress_h2_chol_5d

!******************************************************************************** 
!
!       Modified Cholesky decompoisition directly for gaussian basis 
!
!********************************************************************************
subroutine compress_h2_chol_gauss(hdes, bs, h2_gres)
  type(hamil_descriptor), intent(inout) :: hdes
  type(gauss_basis_set), intent(in)     :: bs
  real(wp), allocatable, intent(inout)  :: h2_gres(:,:,:)
  !local variables
  real(wp), allocatable                 :: h2_spin(:,:)
  
  call iterative_cholesky_gauss_batch(bs, h2_spin, hdes%chol_tol)
  call get_h2_gres_from_fullspin(h2_gres, h2_spin, hdes%ispin1, dealloc=.true.)
  hdes%ng = size(h2_gres, 2)
  hdes%integral_mode = "cholesky"
end subroutine compress_h2_chol_gauss


!******************************************************************************** 
!
!       Spectral decomposition of ERIs for input h2(n*m,n*m,ispin2)
!
!********************************************************************************
subroutine compress_h2_eig_3d(hdes, h2, h2_gres)
  type(hamil_descriptor), intent(inout) :: hdes
  real(wp), allocatable, intent(inout)  :: h2(:,:,:)
  real(wp), allocatable, intent(inout)  :: h2_gres(:,:,:)

  call move_alloc(h2, h2_gres)
  call h2_eig(h2_gres, hdes%chol_tol)
  hdes%ng = size(h2_gres, 2)
  hdes%integral_mode = "cholesky"
end subroutine compress_h2_eig_3d

!******************************************************************************** 
!
!       Modified Cholesky decompoisition for input h2(n,m,n,m,ispin2)
!
!********************************************************************************
subroutine compress_h2_eig_5d(hdes, h2, h2_gres)
  type(hamil_descriptor), intent(inout) :: hdes
  real(wp), allocatable, intent(inout)  :: h2(:,:,:,:,:)
  real(wp), allocatable, intent(inout)  :: h2_gres(:,:,:)
  !local variables
  real(wp), allocatable                 :: h2_spin(:,:)
  
  call get_h2_fullspin(h2, h2_spin, dealloc=.true.)
  call h2_eig(h2_spin, hdes%chol_tol)
  call get_h2_gres_from_fullspin(h2_gres, h2_spin, hdes%ispin1, dealloc=.true.)
  hdes%ng = size(h2_gres, 2)
  hdes%integral_mode = "cholesky"
end subroutine compress_h2_eig_5d


!******************************************************************************** 
!
!       Further compress h2_gres using mCD for input h2_gres(n*m,ng,ispin1)
!
!******************************************************************************** 
subroutine compress_h2_gres_chol_3d(hdes, h2_gres)
  type(hamil_descriptor), intent(inout)  :: hdes
  real(wp), allocatable, intent(inout)   :: h2_gres(:,:,:)
  !local variables
  real(wp), allocatable                  :: h2_spin(:,:)
  
  call iterative_cholesky_compress(h2_gres, hdes%chol_tol, lverb=.true.)
  hdes%ng = size(h2_gres, 2)
  hdes%integral_mode = "cholesky"
end subroutine compress_h2_gres_chol_3d

!******************************************************************************** 
!
!       Further compress h2_gres using mCD for input h2_gres(n,m,ng,ispin1)
!
!******************************************************************************** 
subroutine compress_h2_gres_chol_4d(hdes, h2_gres)
  type(hamil_descriptor), intent(inout)  :: hdes
  real(wp), allocatable, intent(inout)   :: h2_gres(:,:,:,:)
  !local variables
  real(wp), allocatable                  :: h2_spin(:,:)
  integer                                :: n, ispin

  n = size(h2_gres, 1)
  ispin = size(h2_gres, 4)
  
  call get_h2_gres_fullspin(h2_gres, h2_spin, dealloc=.true.)
  call iterative_cholesky_compress(h2_spin, hdes%chol_tol, lverb=.true.)
  call get_h2_gres_from_fullspin(h2_gres, h2_spin, n, ispin, dealloc=.true.)
  hdes%ng = size(h2_gres, 3)
  hdes%integral_mode = "cholesky"
end subroutine compress_h2_gres_chol_4d


!******************************************************************************** 
!
!       Further compress h2_gres using SVD/EIG for input h2_gres(n*m,ng,ispin1)
!
!******************************************************************************** 
subroutine compress_h2_gres_eig_3d(hdes, h2_gres)
  type(hamil_descriptor), intent(inout)  :: hdes
  real(wp), allocatable, intent(inout)   :: h2_gres(:,:,:)
  !local variables
  real(wp), allocatable                  :: h2_spin(:,:)
  
  call h2_gres_eig(h2_gres, hdes%chol_tol)
  hdes%ng = size(h2_gres, 2)
  hdes%integral_mode = "cholesky"
end subroutine compress_h2_gres_eig_3d

!******************************************************************************** 
!
!       Further compress h2_gres using SVD/EIG for input h2_gres(n,m,ng,ispin1)
!
!******************************************************************************** 
subroutine compress_h2_gres_eig_4d(hdes, h2_gres)
  type(hamil_descriptor), intent(inout)  :: hdes
  real(wp), allocatable, intent(inout)   :: h2_gres(:,:,:,:)
  !local variables
  real(wp), allocatable                  :: h2_spin(:,:)
  integer                                :: n, ispin

  n = size(h2_gres, 1)
  ispin = size(h2_gres, 4)
  
  call get_h2_gres_fullspin(h2_gres, h2_spin, dealloc=.true.)
  call h2_gres_eig(h2_spin, hdes%chol_tol)
  call get_h2_gres_from_fullspin(h2_gres, h2_spin, n, ispin, dealloc=.true.)
  hdes%ng = size(h2_gres, 3)
  hdes%integral_mode = "cholesky"
end subroutine compress_h2_gres_eig_4d


!******************************************************************************** 
!
!       Spectral decomposition of the ERIs
!
!           h2(N^2 x N^2)  ---> h2_gres(N^2 x N_g)
!
!******************************************************************************** 
subroutine h2_eig_2d(h2, tol)
  real(wp), allocatable, intent(inout) :: h2(:,:)
  real(wp), intent(in)                 :: tol
  !local variables
  real(wp), allocatable                :: sigma(:)
  integer                              :: size_, ng, g, i, j
  
  size_ = size(h2, 1)
  allocate(sigma(size_))
  call syev(h2, sigma)

  ng = 0
  do g = 1, size_
    if (sigma(g) >= tol) then
      ng = ng + 1
      h2(:,ng) = h2(:,g) * sqrt(sigma(g))
    end if
  end do

  call resize_h2_gres(ng, h2)
end subroutine h2_eig_2d


!******************************************************************************** 
!
!       Spectral decomposition of the ERIs
!
!           h2(N^2 x N^2, ispin2)  ---> h2_gres(N^2 x N_g, ispin1)
!
!******************************************************************************** 
subroutine h2_eig_3d(h2, tol)
  real(wp), allocatable, intent(inout) :: h2(:,:,:)
  real(wp), intent(in)                 :: tol
  !local variables
  integer                              :: ispin            
  real(wp), allocatable                :: h2_spin(:,:)
  
  ispin = size(h2, 3)
  call get_h2_fullspin(h2, h2_spin, dealloc=.true.)
  call h2_eig(h2_spin, tol)
  call get_h2_gres_from_fullspin(h2, h2_spin, ispin, dealloc=.true.)
end subroutine h2_eig_3d


!******************************************************************************** 
!
!       SVD decomposition of h2_gres
!
!           h2_gres(N^2 x Ng)  ---> h2_gres(N^2 x N_gamma)   with N_gamma <= Ng
!
!******************************************************************************** 
subroutine h2_gres_eig_2d(h2_gres, tol)
  real(wp), allocatable, intent(inout) :: h2_gres(:,:)
  real(wp), intent(in)                 :: tol
  !local variables
  real(wp), allocatable                :: sigma(:), temp(:,:), h2_gres_temp(:,:)
  integer                              :: size_, ng, ng_, g
  
  size_ = size(h2_gres, 1)
  ng = size(h2_gres, 2)

  if (size_ <= ng) then
    call get_h2_full(h2_gres)
    call h2_eig(h2_gres, tol)
  else 
    allocate(temp(ng,ng))
    call gemm("t", "n", ng, ng, size_, 1.0_wp, h2_gres, size_, h2_gres, size_, 0.0_wp, temp, ng)
    allocate(sigma(ng))
    call syev(temp, sigma)

    ng_ = 0
    do g = 1, ng
      if (sigma(g) >= tol) then
        ng_= ng_ + 1
        temp(:,ng_) = temp(:,g)
      end if 
    end do

    allocate(h2_gres_temp(size_, ng_))
    call gemm("n", "n", size_, ng_, ng, 1.0_wp, h2_gres, size_, temp(:,1:ng_), ng, 0.0_wp, h2_gres_temp, size_)
    call move_alloc(h2_gres_temp, h2_gres)
  end if
end subroutine h2_gres_eig_2d

!******************************************************************************** 
!
!       SVD decomposition of h2_gres
!
!           h2_gres(N^2 x Ng, ispin1)  ---> h2_gres(N^2 x N_gamma, ispin1)
!
!******************************************************************************** 
subroutine h2_gres_eig_3d(h2_gres, tol)
  real(wp), allocatable, intent(inout) :: h2_gres(:,:,:)
  real(wp), intent(in)                 :: tol
  !local variables
  integer                              :: ispin            
  real(wp), allocatable                :: h2_spin(:,:)
  
  ispin = size(h2_gres, 3)
  call get_h2_gres_fullspin(h2_gres, h2_spin, dealloc=.true.)
  call h2_gres_eig(h2_spin, tol)
  call get_h2_gres_from_fullspin(h2_gres, h2_spin, ispin, dealloc=.true.)
end subroutine h2_gres_eig_3d


!******************************************************************************** 
!
!       Iterative Cholesky decomposition
!
!             H(NxN)  --->  H(NxNg) where Ng <= N
!
!******************************************************************************** 
subroutine iterative_cholesky_2d(h2, chol_tol, lverb)
  real(wp), allocatable, intent(inout) :: h2(:,:)
  real(wp), intent(in)                 :: chol_tol
  logical, optional, intent(in)        :: lverb
  !local variables
  integer                              :: frac, n, ng, ng_est, ind, i, j, size_, temp_size
  logical                              :: lverb_
  real(wp)                             :: chol_res
  real(wp), allocatable                :: residual_matrix(:), temp(:,:)
  character(len=*), parameter          :: proc_name = "iterative_cholesky"

  if (use_profiling) call start_profiling(proc_name)

  lverb_ = .false.
  if (present(lverb)) lverb_ = lverb

  size_ = size(h2, 1)
  frac = int(2 - log10(chol_tol))
  n = nint(sqrt(real(size_)))
  ng_est = frac * n
  temp_size = ng_est

  allocate(residual_matrix(size_), temp(size_, temp_size))
  temp = 0.0_wp

  do i = 1, size_
    residual_matrix(i) = h2(i,i)
  end do

  if (lverb_) call print_chol_header(size_, 0, ng_est)
    
  ng = 0
  do i = 1, size_
    call find_max(residual_matrix, ind, chol_res)
    if (chol_res <= chol_tol) exit
    call schwarz_screen(residual_matrix, chol_res, chol_tol)
    ng = ng + 1
    if (ng > temp_size) then
      temp_size = temp_size + n
      call resize_h2_gres(temp_size, temp)
    end if
    temp(:,ng) = h2(:,ind)
    do j = 1, ng-1
      temp(:,ng) = temp(:,ng) - temp(:,j) * temp(ind,j)
    end do
    temp(:,ng) = temp(:,ng) / sqrt(chol_res)
    call update_residual_matrix(residual_matrix, temp, ng)
    if (lverb_) call print_chol_info(ng, ng_est, chol_res, chol_tol)
  end do
    
  if (.not. chol_res<=chol_tol) then 
    if (comm_world%mpirank == 0) write(*,*) "Iterative Cholesky decomposition did not converge"
  end if

  if (lverb_) call print_chol_footer(ng, chol_res)
  
  call resize_h2_gres(ng, temp)
  call move_alloc(temp, h2)

  if (use_profiling) call end_profiling(proc_name)
end subroutine iterative_cholesky_2d

subroutine iterative_cholesky_3d(h2, chol_tol, lverb)
  real(wp), allocatable, intent(inout) :: h2(:,:,:)
  real(wp), intent(in)                 :: chol_tol
  logical, optional, intent(in)        :: lverb
  !local variables
  integer                              :: ispin
  real(wp), allocatable                :: h2_spin(:,:)

  ispin = size(h2, 3)
  call get_h2_fullspin(h2, h2_spin, dealloc=.true.)
  call iterative_cholesky(h2_spin, chol_tol, lverb)
  call get_h2_gres_from_fullspin(h2, h2_spin, ispin, dealloc=.true.)
end subroutine iterative_cholesky_3d


!******************************************************************************** 
!
!       Compression of the partially done Cholesky decomposition
!
!             H(NxNg)  --->   H(NxN_gamma) where N_gamma <= N_g
!
!******************************************************************************** 
subroutine iterative_cholesky_compress_2d(h2_gres, chol_tol, lverb)
  real(wp), allocatable, intent(inout) :: h2_gres(:,:)
  real(wp), intent(in)                 :: chol_tol
  logical, optional, intent(in)        :: lverb
  !local variables
  integer                              :: frac, n, ng, ng_est, size_, temp_size
  integer                              :: index_, i, j, g, spin
  logical                              :: lverb_
  real(wp)                             :: chol_res
  real(wp), allocatable                :: pivots(:), coulomb_row(:), temp(:,:)
  character(len=*), parameter          :: proc_name = "iterative_cholesky_compress"
  
  if (use_profiling) call start_profiling(proc_name)

  size_ = size(h2_gres, 1)
  ng = size(h2_gres, 2)

  if (size_ <= ng) then
    call get_h2_full(h2_gres)
    call iterative_cholesky(h2_gres, chol_tol, lverb)
  else
    lverb_ = .false.
    if (present(lverb)) lverb_ = lverb

    frac = int(2 - log10(chol_tol))
    n = nint(sqrt(real(size_)))
    ng_est = frac * n
    temp_size = ng_est

    allocate(temp(size_, temp_size))
    temp = 0.0_wp
  
    if (lverb_) call print_chol_header(size_, size(h2_gres, 2), ng_est)
    call calculate_coulomb_integral_pivot(h2_gres, pivots)
       
    ng = 0
    do i = 1, size_
      call find_max(pivots, index_, chol_res)
      if (chol_res <= chol_tol) exit
      call schwarz_screen(pivots, chol_res, chol_tol)
      call calculate_coulomb_integral_row(h2_gres, coulomb_row, index_)
      do j = 1, size_
        do g = 1, ng
          coulomb_row(j) = coulomb_row(j) - temp(j,g) * temp(index_,g) 
        end do
      end do

      coulomb_row = coulomb_row / sqrt(chol_res)
      ng = ng + 1
      if (ng > temp_size) then
        temp_size = temp_size + n
        call resize_h2_gres(temp_size, temp)
      end if

      temp(:,ng) = coulomb_row
      pivots = pivots - coulomb_row**2
      if (lverb_) call print_chol_info(ng, ng_est, chol_res, chol_tol)
    end do

    if (lverb_) call print_chol_footer(ng, chol_res)
  
    call resize_h2_gres(ng, temp)
    call move_alloc(temp, h2_gres)
  end if

  call end_profiling(proc_name)
end subroutine iterative_cholesky_compress_2d

subroutine iterative_cholesky_compress_3d(h2_gres, chol_tol, lverb)
  real(wp), allocatable, intent(inout) :: h2_gres(:,:,:)
  real(wp), intent(in)                 :: chol_tol
  logical, optional, intent(in)        :: lverb
  !local variables
  integer                              :: ispin
  real(wp), allocatable                :: h2_spin(:,:)

  ispin = size(h2_gres, 3)
  call get_h2_gres_fullspin(h2_gres, h2_spin, dealloc=.true.)
  call iterative_cholesky_compress(h2_spin, chol_tol, lverb)
  call get_h2_gres_from_fullspin(h2_gres, h2_spin, ispin, dealloc=.true.)
end subroutine iterative_cholesky_compress_3d


!******************************************************************************** 
!
! Direct construction of the Cholesky vectors from the Gaussian basis set
!
!******************************************************************************** 
subroutine iterative_cholesky_gauss(bs, h2_gres, chol_tol)
  type(gauss_basis_set), intent(in)     :: bs
  real(wp), allocatable, intent(inout)  :: h2_gres(:,:)
  real(wp), intent(in)                  :: chol_tol
  !local variables
  logical, parameter                    :: reuse_shell = .true.
  integer                               :: frac, nsh, nint, nshells, nscr
  integer                               :: indx, indx_, indx_sh, ng, ng_est, n, size_, temp_size
  integer                               :: a, b, sh_start, sh_end
  real(wp)                              :: chol_res, shell_res
  logical, allocatable                  :: lscreen(:)
  real(wp), allocatable                 :: V_ab(:,:), pivots(:), coulomb_row(:)
  character(len=*), parameter           :: proc_name = "cholesky_gauss"
  
  if (use_profiling) call start_profiling(proc_name)

  n = bs%get_norbitals()
  nshells = bs%get_nshells()**2
  size_ = n * n
  frac = int(2 - log10(chol_tol))
  ng_est = frac * n
  temp_size = ng_est

  if (allocated(h2_gres)) deallocate(h2_gres)
  allocate(h2_gres(size_, temp_size))
  allocate(coulomb_row(size_))

  call print_chol_header_gauss(size_, 0, ng_est, chol_tol)
 
  call collect_int_2e_diag(bs, pivots)

  ng = 0
  nint = 0
  nsh = 0 
  nscr = 0 
  chol_res = 1.0_wp
  do while (chol_res > chol_tol)
    call find_max(pivots, indx_, chol_res)
    call schwarz_screen_gauss(pivots, indx_, chol_res, chol_tol, lscreen, nscr)

    call bs%get_shells_from_indx(indx_, a, b)
    call collect_int_2e_2sh(bs, a, b, V_ab)
    sh_start = bs%get_start_indx(a,b) + 1
    sh_end = bs%get_start_indx(a,b) + bs%shell(a)%get_shell_size()*bs%shell(b)%get_shell_size()
    nint = nint + bs%shell(a)%get_shell_size()*bs%shell(b)%get_shell_size()
    nsh = nsh + 1

    shell_res = chol_res
    do while (shell_res > max(chol_tol, chol_res/1000.0_wp))
      call find_max(pivots(sh_start:sh_end), indx, shell_res)
      indx_ = indx + sh_start - 1
      coulomb_row = V_ab(:,indx)
      call update_coulomb_row(coulomb_row, h2_gres(:,1:ng), h2_gres(indx_,1:ng), shell_res)
      lscreen(indx_) = .true.
      ng = ng + 1

      if (ng > temp_size) then 
        temp_size = temp_size + n
        call resize_h2_gres(temp_size, h2_gres)
      end if

      h2_gres(:,ng) = coulomb_row
      pivots = pivots - coulomb_row**2
      call find_max(pivots(sh_start:sh_end), indx, shell_res)

      call print_chol_info_gauss(ng, ng_est, chol_res, real(nsh,wp)/nshells, real(ng,wp)/size_, &
                                 real(nint,wp)/size_, real(nscr,wp)/size_)
      if (.not. reuse_shell) shell_res = 0.0_wp
    end do
  end do
  
  call resize_h2_gres(ng, h2_gres)
  call shell_reorder(bs, h2_gres)
  call print_chol_footer(ng, chol_res)

  if (use_profiling) call end_profiling(proc_name)
end subroutine iterative_cholesky_gauss


!******************************************************************************** 
!
! Direct construction of the Cholesky vectors from the Gaussian basis set
!
!******************************************************************************** 
subroutine iterative_cholesky_gauss_batch(bs, h2_gres, chol_tol)
  type(gauss_basis_set), intent(in)     :: bs
  real(wp), allocatable, intent(inout)  :: h2_gres(:,:)
  real(wp), intent(in)                  :: chol_tol
  !local variables
  logical, parameter                    :: reuse_shell = .true.
  integer                               :: batch_size, frac, nsh, nint, nshells, nscr
  integer                               :: indx, indx_, indx_sh, ng, ng_est, n, size_, temp_size
  integer                               :: a, b, sh_start, sh_end
  real(wp)                              :: chol_res, shell_res
  logical, allocatable                  :: lscreen(:)
  real(wp), allocatable                 :: V_ab(:,:), pivots(:)
  type(cholesky_buffer)                 :: chol_buffer
  character(len=*), parameter           :: proc_name = "cholesky_gauss"
  
  if (use_profiling) call start_profiling(proc_name)

  n = bs%get_norbitals()
  nshells = bs%get_nshells()**2
  size_ = n * n
  frac = int(2 - log10(chol_tol))
  ng_est = frac * n
  temp_size = ng_est
  batch_size = 50
  batch_size = min(n, batch_size)

  if (allocated(h2_gres)) deallocate(h2_gres)
  allocate(h2_gres(size_, temp_size))
  call chol_buffer%init(size_, batch_size)

  call print_chol_header_gauss(size_, 0, ng_est, chol_tol)
 
  call collect_int_2e_diag(bs, pivots)

  ng = 0
  nint = 0
  nsh = 0 
  nscr = 0 
  chol_res = 1.0_wp
  do while (chol_res > chol_tol)
    call find_max(pivots, indx_, chol_res)
    call schwarz_screen_gauss(pivots, indx_, chol_res, chol_tol, lscreen, nscr)

    call bs%get_shells_from_indx(indx_, a, b)
    call collect_int_2e_2sh(bs, a, b, V_ab)
    sh_start = bs%get_start_indx(a,b) + 1
    sh_end = bs%get_start_indx(a,b) + bs%shell(a)%get_shell_size()*bs%shell(b)%get_shell_size()
    nint = nint + bs%shell(a)%get_shell_size()*bs%shell(b)%get_shell_size()
    nsh = nsh + 1

    shell_res = chol_res
    do while (shell_res > max(chol_tol, chol_res/1000.0_wp))
      call find_max(pivots(sh_start:sh_end), indx, shell_res)
      indx_ = indx + sh_start - 1
      call chol_buffer%append(V_ab(:,indx), indx_)
      pivots(indx_) = 0.0_wp
      lscreen(indx_) = .true.
      if (chol_buffer%is_full()) call chol_buffer%add_to_h2(h2_gres, pivots, ng, chol_tol)

      call find_max(pivots(sh_start:sh_end), indx, shell_res)

      call print_chol_info_gauss(ng, ng_est, chol_res, real(nsh,wp)/nshells, real(ng,wp)/size_, &
                                   real(nint,wp)/size_, real(nscr,wp)/size_)

      if (.not. reuse_shell) shell_res = 0.0_wp
    end do
  end do
  
  call resize_h2_gres(ng, h2_gres)
  call shell_reorder(bs, h2_gres)
  call print_chol_footer(ng, chol_res)

  if (use_profiling) call end_profiling(proc_name)
end subroutine iterative_cholesky_gauss_batch


!******************************************************************************** 
!
! Remove contribution of the already calculated Cholesky vectors 
!
!******************************************************************************** 
subroutine update_coulomb_row(coulomb_row, h2, h2_v, chol_res)
  !use mkl_service
  real(wp), intent(inout)     :: coulomb_row(:)
  real(wp), intent(in)        :: h2(:,:)
  real(wp), intent(in)        :: h2_v(:)
  real(wp), intent(in)        :: chol_res
  !local variables
  integer                     :: i, j
  integer                     :: ng, size_, chunk_size, chunk_size_
  character(len=*), parameter :: proc_name = "update_coulomb_row"

  if (profile_code) call start_profiling(proc_name)

  size_ = size(h2, 1)
  ng = size(h2, 2)

  !$omp parallel private(i,j,chunk_size,chunk_size_)
  chunk_size = ceiling(real(size_,wp)/omp_get_num_threads())
  i = omp_get_thread_num()*chunk_size + 1
  j = min((1+omp_get_thread_num())*chunk_size, size_)
  chunk_size_ = j - i + 1
  coulomb_row(i:j) = (coulomb_row(i:j) - matmul(h2(i:j,:), h2_v)) / sqrt(chol_res)
  !$omp end parallel

!debug: for threaded blas probably the best version
  !!call mkl_set_num_threads(48)
  !!call gemv("n", size_, ng, -1.0_wp, h2, size_, h2_v, 1, 1.0_wp, coulomb_row, 1)
  !!call mkl_set_num_threads(1)
  !!coulomb_row = coulomb_row / sqrt(chol_res)
  if (profile_code) call end_profiling(proc_name) 
end subroutine update_coulomb_row


!******************************************************************************** 
!
!       Find maximal element in residual matrix
!
!debug:
! Be careful: choice of the algorithm can change the order of g vectors
! which can affect afqmc results within statistics
!
!******************************************************************************** 
subroutine find_max(array, index_, element)
  real(wp), intent(in)  :: array(:)
  integer, intent(out)  :: index_
  real(wp), intent(out) :: element
  !local variables
  integer               :: i 
  real(wp), parameter   :: tol = 1.0e-12_wp
  
  index_ = 1 
  element = array(1)
  do i = 2, size(array)
    if (array(i)-element > tol) then
      element = array(i)
      index_ = i    
    end if
  end do
end subroutine find_max


!******************************************************************************** 
!
!       Cauchy-Schwarz like prescreening: all elements on diagaonal with 
!                       sqrt(|matrix(i,i)|*maxEL)<tresh
!       are set equal to 0
!
!******************************************************************************** 
subroutine schwarz_screen(array, maxel, tol)
  real(wp), intent(inout) :: array(:)
  real(wp), intent(in)    :: maxel
  real(wp), intent(in)    :: tol
  !local variables
  real(wp)                :: tol_
  
  if (tol < 1.0e-08_wp) then
    tol_ = tol
  else
    tol_ = 1.0e-08_wp
  end if

  where (sqrt(abs(array*maxel)) <= tol_)
    array = 0.0_wp
  end where
end subroutine schwarz_screen

subroutine schwarz_screen_gauss(array, indx, maxel, tol, lscreen, nscr)
  real(wp), intent(inout)             :: array(:)
  integer, intent(in)                 :: indx
  real(wp), intent(in)                :: maxel
  real(wp), intent(in)                :: tol
  logical, allocatable, intent(inout) :: lscreen(:)
  integer, intent(inout)              :: nscr 
  !local variables
  integer                             :: i 
  real(wp)                            :: tol_
  
  if (tol < 1.0e-09_wp) then
    tol_ = tol
  else
    tol_ = 1.0e-09_wp
  end if

  if (.not. allocated(lscreen)) then
    allocate(lscreen(size(array)))
    lscreen = .false.
  end if

  lscreen(indx) = .true.

  do i = 1, size(array)
    if (sqrt(abs(array(i)*maxel))<=tol_ .and. .not. lscreen(i)) then
      array(i) = 0.0_wp
      lscreen(i) = .true.
      nscr = nscr + 1
    end if
  end do
end subroutine schwarz_screen_gauss


!******************************************************************************** 
!
!       Update residual matrix according to
!                       M(i,j,N+1) = M(i,j,N) - L(i,N+1)L(j,N+1)
!       where 
!               i,j     - matrix indices
!               N       - iteration index
!               M       - current residual matrix
!               L       - Cholesky vector
!
!************************************************** 
subroutine update_residual_matrix(residual_matrix, matrix2, ind)
  real(wp), intent(inout) :: residual_matrix(:)
  real(wp), intent(in)    :: matrix2(:,:)
  integer,intent(in)      :: ind
  
  residual_matrix = residual_matrix - matrix2(:,ind) * matrix2(:,ind)
end subroutine update_residual_matrix


!******************************************************************************** 
!
!       Print header for the compression of the Cholesky vectors
!
!******************************************************************************** 
subroutine print_chol_header(size_, ng, ng_est)
  integer, intent(in)           :: size_
  integer, intent(in)           :: ng
  integer, intent(in)           :: ng_est
  
  if (comm_world%mpirank == 0) then
    write(io%screen%funit,*)
    write(io%screen%funit,*)       "Cholesky compresion"
    write(io%screen%funit,*)
    write(io%screen%funit,*)       "current number of the Cholesky vectors             = ", ng
    write(io%screen%funit,*)       "maximal number of the Cholesky vectors             = ", size_
    write(io%screen%funit,*)       "estimated number of the Cholesky vectors           = ", ng_est
    write(io%screen%funit,*)
    write(io%screen%funit,100)     "  # of Chol. vectors  ", "      residual value     ", "     tolerance      "
    write(io%screen%funit,100)     "  ------------------  ", "      --------------     ", "     ---------      "
    
    write(io%qmcfort_log%funit,*)
    write(io%qmcfort_log%funit,*)   "Cholesky compression"
    write(io%qmcfort_log%funit,*)
    write(io%qmcfort_log%funit,*)   "current number of the Cholesky vectors             = ", ng
    write(io%qmcfort_log%funit,*)   "maximal number of the Cholesky vectors             = ", size_
    write(io%qmcfort_log%funit,*)   "estimated number of the Cholesky vectors           = ", ng_est
    write(io%qmcfort_log%funit,*)
    write(io%qmcfort_log%funit,100) "  # of Chol. vectors  ", "      residual value     ", "     tolerance      "
    write(io%qmcfort_log%funit,100) "  ------------------  ", "      --------------     ", "     ---------      "
    
    write(io%qmcfort_out%funit,*) 
    write(io%qmcfort_out%funit,*)   "Cholesky compression"
    write(io%qmcfort_out%funit,*) 
    write(io%qmcfort_out%funit,*)   "current number of the Cholesky vectors             = ", ng
    write(io%qmcfort_out%funit,*)   "maximal number of the Cholesky vectors             = ", size_ 
    write(io%qmcfort_out%funit,*)   "upper limit for the number of the Cholesky vectors = ", ng_est
    write(io%qmcfort_out%funit,*)
    write(io%qmcfort_out%funit,100) "  # of Chol. vectors  ", "    residual value   ", "   tolerance      "
    write(io%qmcfort_out%funit,100) "  ------------------  ", "    --------------   ", "   ---------      "
  end if
 
  100 format(1x,t5,a,t26,a,t49,a)
end subroutine print_chol_header


!******************************************************************************** 
!
!       Print header for the compression of the Cholesky vectors
!
!******************************************************************************** 
subroutine print_chol_header_gauss(size_, ng, ng_est, chol_tol)
  integer, intent(in)           :: size_
  integer, intent(in)           :: ng
  integer, intent(in)           :: ng_est
  real(wp), intent(in)          :: chol_tol
  
  if (comm_world%mpirank == 0) then
    write(io%screen%funit,*)
    write(io%screen%funit,103)     "Cholesky compresion"
    write(io%screen%funit,*)
    write(io%screen%funit,101)     "current number of the Cholesky vectors             = ", ng
    write(io%screen%funit,101)     "maximal number of the Cholesky vectors             = ", size_
    write(io%screen%funit,101)     "estimated number of the Cholesky vectors           = ", ng_est
    write(io%screen%funit,102)     "Cholesky threshold                                 = ", chol_tol 
    write(io%screen%funit,*)
    write(io%screen%funit,100)     "  # Chol. vectors ", "  residual value   ", " # sh pairs ", " # acc. int ", & 
                                   " # calc. int ", " # scr. int "
    write(io%screen%funit,100)     "  --------------- ", "  --------------   ", " ---------- ", " ---------- ", &
                                   " ----------- ", " ---------- "
    
    write(io%qmcfort_log%funit,*)
    write(io%qmcfort_log%funit,103) "Cholesky compresion"
    write(io%qmcfort_log%funit,*)
    write(io%qmcfort_log%funit,101) "current number of the Cholesky vectors             = ", ng
    write(io%qmcfort_log%funit,101) "maximal number of the Cholesky vectors             = ", size_
    write(io%qmcfort_log%funit,101) "estimated number of the Cholesky vectors           = ", ng_est
    write(io%qmcfort_log%funit,102) "Cholesky threshold                                 = ", chol_tol 
    write(io%qmcfort_log%funit,*)
    write(io%qmcfort_log%funit,100) "  # Chol. vectors ", "  residual value   ", " # sh pairs ", " # acc. int ", & 
                                    " # calc. int ", " # scr. int "
    write(io%qmcfort_log%funit,100) "  --------------- ", "  --------------   ", " ---------- ", " ---------- ", &
                                    " ----------- ", " ---------- "
    
    write(io%qmcfort_out%funit,*)
    write(io%qmcfort_out%funit,103) "Cholesky compresion"
    write(io%qmcfort_out%funit,*)
    write(io%qmcfort_out%funit,101) "current number of the Cholesky vectors             = ", ng
    write(io%qmcfort_out%funit,101) "maximal number of the Cholesky vectors             = ", size_
    write(io%qmcfort_out%funit,101) "estimated number of the Cholesky vectors           = ", ng_est
    write(io%qmcfort_out%funit,102) "Cholesky threshold                                 = ", chol_tol 
    write(io%qmcfort_out%funit,*)
    write(io%qmcfort_out%funit,100) "  # Chol. vectors ", "  residual value   ", " # sh pairs ", " # acc. int ", & 
                                    " # calc. int ", " # scr. int "
    write(io%qmcfort_out%funit,100) "  --------------- ", "  --------------   ", " ---------- ", " ---------- ", &
                                    " ----------- ", " ---------- "
  end if
 
  100 format(1x,t4,a,t22,a,t42,a,t57,a,t69,a,t82,a)
  101 format(1x,a,i8)
  102 format(1x,a,es14.6)
  103 format(1x,a)
end subroutine print_chol_header_gauss


!******************************************************************************** 
!
! Print info for the cholesky_decomposition_direct
!
!******************************************************************************** 
subroutine print_chol_info(ng, ng_est, chol_res, chol_tol)
  use standalone, only: modul
  integer, intent(in)  :: ng
  integer, intent(in)  :: ng_est
  real(wp), intent(in) :: chol_res
  real(wp), intent(in) :: chol_tol
  !local variables
  integer, parameter   :: nprints = 20
  integer              :: print_screen
  
  print_screen = get_print_size(ng_est, nprints)

  if (comm_world%mpirank == 0) then
    if (ng == 0) return
    if (.not. modul(ng, print_screen)) return
        
    write(io%screen%funit,100)      ng, chol_res, chol_tol
    write(io%qmcfort_log%funit,100) ng, chol_res, chol_tol
    write(io%qmcfort_out%funit,100) ng, chol_res, chol_tol
  end if
  
  100 format(1x,t10,i6,t28,es14.6,t50,es14.6)
end subroutine print_chol_info


!******************************************************************************** 
!
! Print info for the cholesky_decomposition_direct
!
!******************************************************************************** 
subroutine print_chol_info_gauss(ng, ng_est, chol_res, nsh, nint, nint_, nscr)
  use standalone, only: modul
  integer, intent(in)  :: ng
  integer, intent(in)  :: ng_est
  real(wp), intent(in) :: chol_res
  real(wp), intent(in) :: nsh
  real(wp), intent(in) :: nint
  real(wp), intent(in) :: nint_
  real(wp), intent(in) :: nscr
  !local variables
  integer, parameter   :: nprints = 20
  integer              :: print_screen, ng_
  integer, save        :: ratio = 0
  
  print_screen = get_print_size(ng_est, nprints)

  if (comm_world%mpirank == 0) then
    if (ng == 0) return
    if (ratio /= ng/print_screen) then
      ratio = ng / print_screen
    else 
      return
    end if
        
    ng_ = ratio * print_screen

    write(io%screen%funit,100)      ng_, chol_res, write_perc(nsh), write_perc(nint), write_perc(nint_), write_perc(nscr)
    write(io%qmcfort_log%funit,100) ng_, chol_res, write_perc(nsh), write_perc(nint), write_perc(nint_), write_perc(nscr)
    write(io%qmcfort_out%funit,100) ng_, chol_res, write_perc(nsh), write_perc(nint), write_perc(nint_), write_perc(nscr)
  end if
  
  100 format(1x,t6,i8,t23,es14.6,t44,a,t59,a,t71,a,t84,a)
end subroutine print_chol_info_gauss


!******************************************************************************** 
!
!       Print info for the cholesky_decomposition_direct
!
!******************************************************************************** 
subroutine print_chol_footer(ng, chol_res)
  integer, intent(in)  :: ng
  real(wp), intent(in) :: chol_res
  
  if (comm_world%mpirank == 0) then
    write(io%screen%funit,*) 
    write(io%screen%funit,100)      "number of cholesky vecotrs found   = ", ng
    write(io%screen%funit,101)      "residual value                     = ", chol_res
  
    write(io%qmcfort_log%funit,*) 
    write(io%qmcfort_log%funit,100) "number of cholesky vecotrs found   = ", ng
    write(io%qmcfort_log%funit,101) "residual value                     = ", chol_res
    
    write(io%qmcfort_out%funit,*) 
    write(io%qmcfort_out%funit,100) "number of cholesky vecotrs found   = ", ng
    write(io%qmcfort_out%funit,101) "residual value                     = ", chol_res
  end if
  
  100 format(1x,a,i6)
  101 format(1x,a,es14.6)
end subroutine print_chol_footer


!******************************************************************************** 
!
! Calculate Coulomb integral diagonal elements (pivots) from the Cholesky vectors
!
!******************************************************************************** 
subroutine calculate_coulomb_integral_pivot(h2_gres, pivots, pivots_sq)
  real(wp), intent(in)                           :: h2_gres(:,:)
  real(wp), allocatable, intent(inout)           :: pivots(:)
  real(wp), allocatable, optional, intent(inout) :: pivots_sq(:)
  !local variables
  integer                                        :: size_, a
  character(len=*), parameter                    :: proc_name = "coulomb_integral_pivot"
  
  if (profile_code) call start_profiling(proc_name)
  
  size_ = size(h2_gres, 1)
  if (.not. allocated(pivots)) allocate(pivots(size_))
  pivots = 0.0_wp
  
  do a = 1, size_
    pivots(a) = dot_product(h2_gres(a,:),h2_gres(a,:))
  end do
   
  if (present(pivots_sq)) then
    if (.not. allocated(pivots_sq)) allocate(pivots_sq, mold=pivots)
    pivots_sq = sqrt(pivots)
  end if

  if (profile_code) call end_profiling(proc_name)
end subroutine calculate_coulomb_integral_pivot


!******************************************************************************** 
!
! Calculate Coulomb integral row from the Cholesky vectors
!
!******************************************************************************** 
subroutine calculate_coulomb_integral_row(h2_gres, coulomb_row, b)
  real(wp), intent(in)                 :: h2_gres(:,:)
  real(wp), allocatable, intent(inout) :: coulomb_row(:)
  integer , intent(in)                 :: b
  !local variables
  integer                              :: size_, a
  character(len=*), parameter          :: proc_name = "calculate_coulomb_integral_row"
  
  if (profile_code) call start_profiling(proc_name)

  size_ = size(h2_gres, 1)
  if (.not. allocated(coulomb_row)) allocate(coulomb_row(size_))
  coulomb_row = 0.0_wp
  
  !$omp parallel do 
  do a = 1, size_
    coulomb_row(a) = dot_product(h2_gres(a,:), h2_gres(b,:))
  end do
  !$omp end parallel do

  if (profile_code) call end_profiling(proc_name)
end subroutine calculate_coulomb_integral_row

end module cholesky