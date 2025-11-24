! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module hamilton_layout


#include "../preproc.inc"
use constants
use mpi, only: comm_world
use profiling
use hamilton_vars
use lapack

implicit none

interface get_h2_full 
  module procedure :: get_h2_full_2d, get_h2_full_3d, get_h2_full_3d_c, get_h2_full_5d
end interface get_h2_full

interface get_h2_fullspin
  module procedure :: get_h2_fullspin_3d, get_h2_fullspin_5d
end interface

interface get_h2_gres_2d
  module procedure get_h2_gres_2d_r, get_h2_gres_2d_c
end interface get_h2_gres_2d

interface get_h2_gres_fullspin
  module procedure :: get_h2_gres_fullspin_3d, get_h2_gres_fullspin_4d
end interface

interface get_h2_gres_from_fullspin
  module procedure :: get_h2_gres_from_fullspin_3d, get_h2_gres_from_fullspin_4d
end interface

interface resize_h2_gres
  module procedure :: resize_h2_gres_2d, resize_h2_gres_3d
end interface 

contains

!******************************************************************************** 
!
!       Assemble full h2(n*m,n*m) from h2_gres(n*m,ng)
!
!******************************************************************************** 
subroutine get_h2_full_2d(h2_gres)
  real(wp), allocatable, intent(inout) :: h2_gres(:,:)
  !local variables
  integer                              :: size_, ng
  real(wp), allocatable                :: temp(:,:)

  size_ = size(h2_gres, 1)
  ng = size(h2_gres, 2)

  allocate(temp(size_, size_))

  call gemm("n", "t", size_, size_, ng, 1.0_wp, h2_gres, size_, & 
              h2_gres, size_, 0.0_wp, temp, size_)

  call move_alloc(temp, h2_gres)
end subroutine get_h2_full_2d

!******************************************************************************** 
!
!       Assemble full h2(n*m,n*m,ispin2) from h2_gres(n*m,ng,ispin1)
!
!******************************************************************************** 
subroutine get_h2_full_3d(h2_gres, ispin2)
  real(wp), allocatable, intent(inout) :: h2_gres(:,:,:)
  integer, intent(in)                  :: ispin2
  !local variables
  integer                              :: size_, ng
  integer                              :: spin, spin1, spin2
  real(wp), allocatable                :: temp(:,:,:)
  character(len=*), parameter          :: proc_name = "get_h2_full"
  character(len=*), parameter          :: gemm_name = "gemm_h2_full"

  if (use_profiling) call start_profiling(proc_name)

  size_ = size(h2_gres, 1)
  ng = size(h2_gres, 2)

  allocate(temp(size_, size_, ispin2))

  do spin = 1, ispin2
    if (spin == 3) then
      spin1 = 1
      spin2 = 2
    else
      spin1 = spin 
      spin2 = spin
    end if
    if (profile_code) call start_profiling(gemm_name)
    call gemm("n", "t", size_, size_, ng, oner, h2_gres(:,:,spin1), size_, & 
              h2_gres(:,:,spin2), size_, zeror, temp(:,:,spin), size_)
    if (profile_code) call end_profiling(gemm_name, size_, size_, ng, "r")
  end do

  call move_alloc(temp, h2_gres)

  if (use_profiling) call end_profiling(proc_name)
end subroutine get_h2_full_3d

subroutine get_h2_full_3d_c(h2_gres, ispin2)
  complex(wp), allocatable, intent(inout) :: h2_gres(:,:,:)
  integer, intent(in)                     :: ispin2
  !local variables
  integer                                 :: size_, ng
  integer                                 :: spin, spin1, spin2
  complex(wp), allocatable                :: temp(:,:,:)
  character(len=*), parameter             :: proc_name = "get_h2_full"
  character(len=*), parameter             :: gemm_name = "gemm_h2_full"

  if (use_profiling) call start_profiling(proc_name)

  size_ = size(h2_gres, 1)
  ng = size(h2_gres, 2)

  allocate(temp(size_, size_, ispin2))

  do spin = 1, ispin2
    if (spin == 3) then
      spin1 = 1
      spin2 = 2
    else
      spin1 = spin 
      spin2 = spin
    end if
    if (profile_code) call start_profiling(gemm_name)
    call gemm("n", "t", size_, size_, ng, onec, h2_gres(:,:,spin1), size_, & 
              h2_gres(:,:,spin2), size_, zeroc, temp(:,:,spin), size_)
    if (profile_code) call end_profiling(gemm_name, size_, size_, ng, "c")
  end do

  call move_alloc(temp, h2_gres)

  if (use_profiling) call end_profiling(proc_name)
end subroutine get_h2_full_3d_c


!******************************************************************************** 
!
!       Assemble full h2(n,m,n,m,ispin2) from h2_gres(n*m,ng,ispin1)
!
!******************************************************************************** 
subroutine get_h2_full_5d(h2_gres, h2, n, m, ispin2, dealloc)
  real(wp), allocatable, intent(inout) :: h2_gres(:,:,:)
  real(wp), allocatable, intent(inout) :: h2(:,:,:,:,:)
  integer, intent(in)                  :: n
  integer, intent(in)                  :: m
  integer, intent(in)                  :: ispin2
  logical, optional, intent(in)        :: dealloc
  !local variables
  logical                              :: dealloc_
  
  dealloc_ = .true.
  if (present(dealloc)) dealloc_ = dealloc

  call get_h2_full(h2_gres, ispin2)
  
  if (allocated(h2)) deallocate(h2)
  allocate(h2(n,m,n,m,ispin2))

  h2 = reshape(h2_gres, shape=[n,m,n,m,ispin2])
  if (dealloc_) deallocate(h2_gres)
end subroutine get_h2_full_5d


!******************************************************************************** 
!
!       Increase number of cholesky vectors by inc
! 
!             h2_gres(n,temp_size,ispin1)  --> h2_gres(n,temp_size+inc,ispin1)
!
!******************************************************************************** 
subroutine increase_h2_gres_size(inc, h2_gres)
  integer, intent(in)                  :: inc
  real(wp), allocatable, intent(inout) :: h2_gres(:,:,:)
  !local varaibles
  integer                              :: size_
  real(wp), allocatable                :: temp(:,:,:)

  size_ = size(h2_gres,2)
  allocate(temp(size(h2_gres,1), size_+inc, size(h2_gres,3)))
  temp = 0.0_wp
  temp(:,1:size_,:) = h2_gres
  call move_alloc(temp, h2_gres)
end subroutine increase_h2_gres_size


!******************************************************************************** 
!
!       Resize Cholesky vectors to desired length
! 
!             h2_gres(size_,temp_size,ispin)  --> h2_gres(size_,ng,ispin)
!
!******************************************************************************** 
subroutine resize_h2_gres_2d(ng, h2_gres)
  integer, intent(in)                  :: ng
  real(wp), allocatable, intent(inout) :: h2_gres(:,:)
  !local varaibles
  integer                              :: size_
  real(wp), allocatable                :: temp(:,:)

  size_ = size(h2_gres,2)
  allocate(temp(size(h2_gres,1), ng))

  temp = 0.0_wp
  if (ng <= size_) then
    temp = h2_gres(:,1:ng)
  else
    temp(:,1:size_) = h2_gres
  end if
  call move_alloc(temp, h2_gres)
end subroutine resize_h2_gres_2d

subroutine resize_h2_gres_3d(ng, h2_gres)
  integer, intent(in)                  :: ng
  real(wp), allocatable, intent(inout) :: h2_gres(:,:,:)
  !local varaibles
  integer                              :: size_
  real(wp), allocatable                :: temp(:,:,:)

  size_ = size(h2_gres,2)
  allocate(temp(size(h2_gres,1), ng, size(h2_gres,3)))

  temp = 0.0_wp
  if (ng <= size_) then
    temp = h2_gres(:,1:ng,:)
  else
    temp(:,1:size_,:) = h2_gres
  end if
  call move_alloc(temp, h2_gres)
end subroutine resize_h2_gres_3d


!******************************************************************************** 
!
!       Transforms h2_gres from orbital basis to pair orbital basis
!
!               h(n,m,ng,ispin)  ---> h(nm,ng,ispin)
!
!******************************************************************************** 
subroutine get_h2_gres_2d_r(h2_gres_3d, h2_gres_2d, dealloc)
  real(wp), allocatable, intent(inout) :: h2_gres_3d(:,:,:,:)
  real(wp), allocatable, intent(inout) :: h2_gres_2d(:,:,:) 
  logical, optional, intent(in)        :: dealloc
  !local variables
  logical                              :: dealloc_
  integer                              :: n, m, ng, ispin
  integer                              :: j, g, lo, up, spin
  character(len=*), parameter          :: proc_name = "get_h2_gres_2d"

  if (use_profiling) call start_profiling(proc_name)
  
  dealloc_ = .true.
  if (present(dealloc)) dealloc_ = dealloc

  n = size(h2_gres_3d, 1)
  m = size(h2_gres_3d, 2)
  ng = size(h2_gres_3d, 3)
  ispin = size(h2_gres_3d, 4)

  if (allocated(h2_gres_2d)) deallocate(h2_gres_2d)
  allocate(h2_gres_2d(n*m,ng,ispin))
  
  do spin = 1, ispin
    do g = 1, ng
      do j = 1, m
        lo = n * (j-1) + 1
        up = n * j
        h2_gres_2d(lo:up,g,spin) = h2_gres_3d(:,j,g,spin)
      end do
    end do
  end do
  
  if (dealloc_) deallocate(h2_gres_3d)

  if (use_profiling) call end_profiling(proc_name)
end subroutine get_h2_gres_2d_r

subroutine get_h2_gres_2d_c(h2_gres_3d, h2_gres_2d, dealloc)
  complex(wp), allocatable, intent(inout) :: h2_gres_3d(:,:,:,:)
  complex(wp), allocatable, intent(inout) :: h2_gres_2d(:,:,:) 
  logical, optional, intent(in)           :: dealloc
  !local variables
  logical                                 :: dealloc_
  integer                                 :: n, m, ng, ispin
  integer                                 :: j, g, lo, up, spin
  character(len=*), parameter             :: proc_name = "get_h2_gres_2d"

  if (use_profiling) call start_profiling(proc_name)
  
  dealloc_ = .true.
  if (present(dealloc)) dealloc_ = dealloc

  n = size(h2_gres_3d, 1)
  m = size(h2_gres_3d, 2)
  ng = size(h2_gres_3d, 3)
  ispin = size(h2_gres_3d, 4)

  if (allocated(h2_gres_2d)) deallocate(h2_gres_2d)
  allocate(h2_gres_2d(n*m,ng,ispin))
  do spin = 1, ispin
    do g = 1, ng
      do j = 1, m
        lo = n * (j-1) + 1
        up = n * j
        h2_gres_2d(lo:up,g,spin) = h2_gres_3d(:,j,g,spin)
      end do
    end do
  end do
  
  if (dealloc_) deallocate(h2_gres_3d)

  if (use_profiling) call end_profiling(proc_name)
end subroutine get_h2_gres_2d_c


!******************************************************************************** 
!
!       Transforms h2_gres_2d from pair orbital basis to orbital basis
!
!               h(nm,ng,ispin1)  ---> h(n,m,ng,ispin1)
!
!******************************************************************************** 
subroutine get_h2_gres_3d(h2_gres_2d, h2_gres_3d, n, m, dealloc)
  real(wp), allocatable, intent(inout) :: h2_gres_2d(:,:,:) 
  real(wp), allocatable, intent(inout) :: h2_gres_3d(:,:,:,:)
  integer, intent(in)                  :: n, m
  logical, optional, intent(in)        :: dealloc
  !local variables
  integer                              :: ng, ispin
  integer                              :: j, g, lo, up, spin
  
  ng = size(h2_gres_2d, 2)
  ispin = size(h2_gres_2d, 3)

  if (allocated(h2_gres_3d)) deallocate(h2_gres_3d)
  allocate(h2_gres_3d(n,m,ng,ispin))
  
  do spin = 1, ispin
    do g = 1, ng
      do j = 1, m
        lo = n * (j-1) + 1
        up = n * j
        h2_gres_3d(:,j,g,spin) = h2_gres_2d(lo:up,g,spin)
      end do
    end do
  end do
  
  if (present(dealloc)) then
    if (dealloc) deallocate(h2_gres_2d)
  end if
end subroutine get_h2_gres_3d


!******************************************************************************** 
!
!       Transform h2 from format h2(N^2,N^2,ispin3) to h2(2N^2, 2N^2)
!
!******************************************************************************** 
subroutine get_h2_fullspin_3d(h2, h2_spin, dealloc)
  real(wp), allocatable, intent(inout) :: h2(:,:,:)
  real(wp), allocatable, intent(inout) :: h2_spin(:,:)
  logical, optional, intent(in)        :: dealloc
  !local variables
  integer                              :: n, m, ispin, ispin_

  n = size(h2, 1)
  m = size(h2, 2)
  ispin = size(h2, 3)
  ispin_ = min(ispin, 2)

  if (allocated(h2_spin)) deallocate(h2_spin)
  allocate(h2_spin(ispin_*n, ispin_*m))

  h2_spin(1:n,1:m) = h2(:,:,1)

  if (ispin_ == 2) then
    h2_spin(n+1:,m+1:) = h2(:,:,2)
    h2_spin(n+1:,1:m) = h2(:,:,3)
    h2_spin(1:n,m+1:) = h2(:,:,3)
  end if

  if (present(dealloc)) then
    if (dealloc) deallocate(h2)
  end if
end subroutine get_h2_fullspin_3d

!******************************************************************************** 
!
!       Transform h2 from format h2(N,N,M,M,ispin3) to h2(2N^2, 2N^2)
!
!******************************************************************************** 
subroutine get_h2_fullspin_5d(h2, h2_spin, dealloc)
  real(wp), allocatable, intent(inout) :: h2(:,:,:,:,:)
  real(wp), allocatable, intent(inout) :: h2_spin(:,:)
  logical, optional, intent(in)        :: dealloc
  !local variables
  integer                              :: n, m, nm, ispin, ispin_, i , j, p, q, r, s

  n = size(h2, 1)
  m = size(h2, 2)
  nm = n * m
  ispin = size(h2, 5)
  ispin_ = min(ispin, 2)

  if (allocated(h2_spin)) deallocate(h2_spin)
  allocate(h2_spin(ispin_*nm, ispin_*nm))

  h2_spin(1:nm,1:nm) = reshape(h2(:,:,:,:,1), shape=[nm, nm])

  if (ispin_ == 2) then
    h2_spin(nm+1:,nm+1:) = reshape(h2(:,:,:,:,2), shape=[nm, nm])
    h2_spin(nm+1:,1:nm) = reshape(h2(:,:,:,:,3), shape=[nm, nm])
    h2_spin(1:nm,nm+1:) = reshape(h2(:,:,:,:,3), shape=[nm, nm])
    h2_spin(nm+1:,1:nm) = transpose(h2_spin(nm+1:,1:nm))
  end if

  if (present(dealloc)) then
    if (dealloc) deallocate(h2)
  end if
end subroutine get_h2_fullspin_5d


!******************************************************************************** 
!
!       Transform h2_gres from format h2(N^2,Ng,ispin1) to h2(2N^2, Ng)
!
!******************************************************************************** 
subroutine get_h2_gres_fullspin_3d(h2, h2_spin, dealloc)
  real(wp), allocatable, intent(inout) :: h2(:,:,:)
  real(wp), allocatable, intent(inout) :: h2_spin(:,:)
  logical, optional, intent(in)        :: dealloc
  !local variables
  integer                              :: n, ng, ispin

  n = size(h2, 1)
  ng = size(h2, 2)
  ispin = size(h2, 3)

  if (allocated(h2_spin)) deallocate(h2_spin)
  allocate(h2_spin(ispin*n, ng))

  h2_spin(1:n,:) = h2(:,:,1)
  if (ispin == 2) h2_spin(n+1:,:) = h2(:,:,2)

  if (present(dealloc)) then
    if (dealloc) deallocate(h2)
  end if
end subroutine get_h2_gres_fullspin_3d

!******************************************************************************** 
!
!       Transform h2_gres from format h2(N,M,Ng,ispin1) to h2(2NM, Ng)
!
!******************************************************************************** 
subroutine get_h2_gres_fullspin_4d(h2, h2_spin, dealloc)
  real(wp), allocatable, intent(inout) :: h2(:,:,:,:)
  real(wp), allocatable, intent(inout) :: h2_spin(:,:)
  logical, optional, intent(in)        :: dealloc
  !local variables
  integer                              :: n, m, nm, ng, ispin

  n = size(h2, 1)
  m = size(h2, 2)
  ng = size(h2, 3)
  ispin = size(h2, 4)
  nm = n * m

  if (allocated(h2_spin)) deallocate(h2_spin)
  allocate(h2_spin(ispin*nm, ng))

  h2_spin(1:nm,:) = reshape(h2(:,:,:,1), shape=[nm, ng]) 
  if (ispin == 2) h2_spin(nm+1:,:) = reshape(h2(:,:,:,2), shape=[nm,ng])

  if (present(dealloc)) then
    if (dealloc) deallocate(h2)
  end if
end subroutine get_h2_gres_fullspin_4d


!******************************************************************************** 
!
!       Obtain h2_gres from h2_gres_fulspin:
!
!               h2_gres_fullspin(ispin*NM,Ng) --> h2_gres(NM,Ng,ispin)
!
!******************************************************************************** 
subroutine get_h2_gres_from_fullspin_3d(h2, h2_spin, ispin, dealloc)
  real(wp), allocatable, intent(inout) :: h2(:,:,:)
  real(wp), allocatable, intent(inout) :: h2_spin(:,:)
  integer, intent(in)                  :: ispin
  logical, optional, intent(in)        :: dealloc
  !local variables
  integer                              :: n, ng

  n = size(h2_spin, 1) / ispin
  ng = size(h2_spin, 2)

  if (allocated(h2)) deallocate(h2)
  allocate(h2(n, ng, ispin))

  h2(:,:,1) = h2_spin(1:n,:)
  if (ispin == 2) h2(:,:,2) = h2_spin(n+1:,:)

  if (present(dealloc)) then
    if (dealloc) deallocate(h2_spin)
  end if
end subroutine get_h2_gres_from_fullspin_3d

!******************************************************************************** 
!
!       Obtain h2_gres from h2_gres_fulspin:
!
!               h2_gres_fullspin(ispin*NM,Ng) --> h2_gres(N,M,Ng,ispin)
!
!******************************************************************************** 
subroutine get_h2_gres_from_fullspin_4d(h2, h2_spin, n, ispin, dealloc)
  real(wp), allocatable, intent(inout) :: h2(:,:,:,:)
  real(wp), allocatable, intent(inout) :: h2_spin(:,:)
  integer, intent(in)                  :: n
  integer, intent(in)                  :: ispin
  logical, optional, intent(in)        :: dealloc
  !local variables
  integer                              :: nm, m, ng

  nm = size(h2_spin, 1)
  ng = size(h2_spin, 2)
  m = nm / n / ispin

  if (allocated(h2)) deallocate(h2)
  allocate(h2(n, m, ng, ispin))

  h2(:,:,:,1) = reshape(h2_spin(1:nm,:), shape=[n, m, ng])
  if (ispin == 2) h2(:,:,:,2) = reshape(h2_spin(nm+1:,:), shape=[n, m, ng])

  if (present(dealloc)) then
    if (dealloc) deallocate(h2_spin)
  end if
end subroutine get_h2_gres_from_fullspin_4d


!******************************************************************************** 
!
!       Transpose h2_gres_2d
!
!               h2_gres_2d(n*n, ng, ispin) --> h2_gres_2d(ng, n*n, ispin)
!       or
!               h2_gres_2d(ng, n*n, ispin) --> h2_gres_2d(n*n, ng, ispin)
!
!******************************************************************************** 
subroutine transpose_h2_gres_2d(h2_gres)
  real(wp), allocatable, intent(inout) :: h2_gres(:,:,:)
  !local variables
  integer                              :: n, m, ispin, spin
  real(wp), allocatable                :: temp(:,:,:)

  n = size(h2_gres, 1)
  m = size(h2_gres, 2)
  ispin = size(h2_gres, 3)

  allocate(temp(m,n,ispin))

  do spin  = 1, ispin
    temp(:,:,spin) = transpose(h2_gres(:,:,spin))
  end do

  call move_alloc(temp, h2_gres)
end subroutine transpose_h2_gres_2d

end module hamilton_layout