! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module statistics

#include "preproc.inc"
use constants
use mpi
use qmcfort_io

implicit none 

interface mean
  module procedure mean, mean_c, mean_w, mean_w_rc, mean_w_c  
  module procedure mean_2d, mean_2d_c, mean_w_2d, mean_w_2d_rc, mean_w_2d_c
  module procedure mean_2d_ax, mean_2d_ax_c, mean_w_2d_ax, mean_w_2d_ax_rc, mean_w_2d_ax_c
end interface mean

interface std
  module procedure std, std_c, std_w, std_w_rc, std_w_c
  module procedure std_2d, std_2d_c, std_w_2d,  std_w_2d_rc, std_w_2d_c
  module procedure std_2d_ax, std_2d_ax_c, std_w_2d_ax, std_w_2d_ax_rc, std_w_2d_ax_c
end interface std

interface std_mean
  module procedure std_mean, std_mean_c, std_mean_w, std_mean_w_c, std_mean_w_rc
  module procedure std_mean_2d, std_mean_2d_c, std_mean_w_2d, std_mean_w_2d_rc, std_mean_w_2d_c
  module procedure std_mean_2d_ax, std_mean_2d_ax_c, std_mean_w_2d_ax, std_mean_w_2d_ax_rc, std_mean_w_2d_ax_c
end interface std_mean

interface cov
  module procedure cov, cov_c, cov_w, cov_w_rc, cov_w_c 
  module procedure cov_2d, cov_2d_c, cov_w_2d, cov_w_2d_rc, cov_w_2d_c
  module procedure cov_2d_ax, cov_2d_ax_c, cov_w_2d_ax, cov_w_2d_ax_rc, cov_w_2d_ax_c
end interface cov

interface autocov
  module procedure autocov_r, autocov_c
end interface autocov

interface autocorr_len
  module procedure autocorr_len_r, autocorr_len_c
end interface autocorr_len

interface blocking
  module procedure blocking, blocking_array, blocking_w, blocking_w_array
  module procedure blocking_c, blocking_array_c, blocking_w_c, blocking_w_rc, blocking_w_array_c, blocking_w_array_rc
end interface blocking 

interface blocking_output
  module procedure blocking_output_real, blocking_output_cmplx
end interface blocking_output

interface histogram
  module procedure histogram_1d, histogram_2d, histogram_w_2d
end interface histogram

contains

!******************************************************************************** 
!
!       Calculate mean value of the array
!
!               \bar x = \sum_i x_i / N
!
!******************************************************************************** 
function mean(x, comm) result(mean_)
  real(wp), intent(in)                         :: x(:)
  type(mpi_communicator), optional, intent(in) :: comm
  real(wp)                                     :: mean_
  !local variables                          
  integer                                      :: size_
  
  mean_ = sum(x) 
  size_ = size(x)

  if (present(comm)) call comm%mpisum(mean_)
  if (present(comm)) call comm%mpisum(size_)

  mean_ = mean_ / size_
end function mean

function mean_c(z, comm) result(mean_)
  complex(wp), intent(in)                      :: z(:)
  type(mpi_communicator), optional, intent(in) :: comm
  complex(wp)                                  :: mean_
  !local variables                          
  integer                                      :: size_
  
  mean_ = sum(z)
  size_ = size(z)

  if (present(comm)) call comm%mpisum(mean_)
  if (present(comm)) call comm%mpisum(size_)

  mean_ = mean_ / size_
end function mean_c

function mean_2d(x, comm) result(mean_) 
  real(wp), intent(in)                         :: x(:,:)
  type(mpi_communicator), optional, intent(in) :: comm
  real(wp)                                     :: mean_
  !local variables:
  integer                                      :: size_
  
  mean_ = sum(x) 
  size_ = size(x)

  if (present(comm)) call comm%mpisum(mean_)
  if (present(comm)) call comm%mpisum(size_)

  mean_ = mean_ / size_
end function mean_2d

function mean_2d_ax(x, ax, comm) result(mean_) 
  real(wp), intent(in)                         :: x(:,:)
  integer, intent(in)                          :: ax
  type(mpi_communicator), optional, intent(in) :: comm
  real(wp), allocatable                        :: mean_(:)
  !local variables:
  integer                                      :: i
  
  if (ax == 1) then
    if (.not. allocated(mean_)) allocate(mean_(size(x, 2)))
    do i = 1, size(mean_)
      mean_(i) = mean(x(:,i), comm)
    end do
  else if (ax == 2) then
    if (.not. allocated(mean_)) allocate(mean_(size(x, 1)))
    do i = 1, size(mean_)
      mean_(i) = mean(x(i,:), comm)
    end do
  else
    write(*,*) "wrong ax value in mean_arrray - program stops now"
    stop
  end if
end function mean_2d_ax

function mean_2d_c(x, comm) result(mean_) 
  complex(wp), intent(in)                      :: x(:,:)
  type(mpi_communicator), optional, intent(in) :: comm
  complex(wp)                                  :: mean_
  !local variables:
  integer                                      :: size_
  
  mean_ = sum(x) 
  size_ = size(x)

  if (present(comm)) call comm%mpisum(mean_)
  if (present(comm)) call comm%mpisum(size_)

  mean_ = mean_ / size_
end function mean_2d_c

function mean_2d_ax_c(z, ax, comm) result(mean_)
  complex(wp), intent(in)                      :: z(:,:)
  integer, intent(in)                          :: ax
  type(mpi_communicator), optional, intent(in) :: comm
  complex(wp), allocatable                     :: mean_(:)
  !local variables:
  integer                                      :: i
  
  if (ax == 1) then
    if (.not. allocated(mean_)) allocate(mean_(size(z, 2)))
    do i = 1, size(mean_)
      mean_(i) = mean(z(:, i), comm)
    end do
  else if (ax == 2) then
    if (.not. allocated(mean_)) allocate(mean_(size(z, 1)))
    do i = 1, size(mean_)
      mean_(i) = mean(z(i,:), comm)
    end do
  else
    write(*,*) "wrong ax value in mean_arrray_c - program stops now"
    stop
  end if
end function mean_2d_ax_c


!******************************************************************************** 
!
!       Calculate mean value for weighted array
!
!               \bar x = \sum_i w_i x_i / v1
!
!               v1 = \sum_i w_i
!
!******************************************************************************** 
function mean_w(w, x, comm) result(mean_)
  real(wp), intent(in)                         :: w(:)
  real(wp), intent(in)                         :: x(:)
  type(mpi_communicator), optional, intent(in) :: comm
  real(wp)                                     :: mean_
  !local variables
  real(wp)                                     :: size_
  
  mean_ = sum(x * w)
  size_ = sum(w)

  if (present(comm)) call comm%mpisum(mean_)
  if (present(comm)) call comm%mpisum(size_)

  mean_ = mean_ / size_
end function mean_w

function mean_w_c(w, z, comm) result(mean_)
  complex(wp), intent(in)                      :: w(:)
  complex(wp), intent(in)                      :: z(:)
  type(mpi_communicator), optional, intent(in) :: comm
  complex(wp)                                  :: mean_
  !local variables
  complex(wp)                                  :: size_
  
  mean_ = sum(z * w)
  size_ = sum(w)

  if (present(comm)) call comm%mpisum(mean_)
  if (present(comm)) call comm%mpisum(size_)
  
  mean_ = mean_ / size_
end function mean_w_c

function mean_w_rc(w, z, comm) result(mean_)
  real(wp), intent(in)                         :: w(:)
  complex(wp), intent(in)                      :: z(:)
  type(mpi_communicator), optional, intent(in) :: comm
  complex(wp)                                  :: mean_
  !local variables
  real(wp)                                     :: size_
  
  mean_ = sum(z * w)
  size_ = sum(w)

  if (present(comm)) call comm%mpisum(mean_)
  if (present(comm)) call comm%mpisum(size_)
  
  mean_ = mean_ / size_
end function mean_w_rc

function mean_w_2d(w, x, comm) result(mean_)
  real(wp), intent(in)                         :: w(:,:)
  real(wp), intent(in)                         :: x(:,:)
  type(mpi_communicator), optional, intent(in) :: comm
  real(wp)                                     :: mean_
  !local variables
  real(wp)                                     :: size_
  
  mean_ = sum(x * w)
  size_ = sum(w)

  if (present(comm)) call comm%mpisum(mean_)
  if (present(comm)) call comm%mpisum(size_)

  mean_ = mean_ / size_
end function mean_w_2d

function mean_w_2d_ax(w, x, ax, comm) result(mean_)
  real(wp), intent(in)                         :: w(:)
  real(wp), intent(in)                         :: x(:,:)
  integer, intent(in)                          :: ax
  type(mpi_communicator), optional, intent(in) :: comm
  real(wp), allocatable                        :: mean_(:)
  !local variables
  integer                                      :: i
  
  if (ax == 1) then
    if (.not. allocated(mean_)) allocate(mean_(size(x, 2)))
    do i = 1, size(mean_)
      mean_(i) = mean(w, x(:,i), comm)
    end do
  else if (ax == 2) then
    if (.not. allocated(mean_)) allocate(mean_(size(x, 1)))
    do i = 1, size(mean_)
      mean_(i) = mean(w, x(i,:), comm)
    end do
  else
    write(*,*) "wrong ax value in mean_w_arrray - program stops now"
    stop
  end if
end function mean_w_2d_ax

function mean_w_2d_c(w, x, comm) result(mean_)
  complex(wp), intent(in)                      :: w(:,:)
  complex(wp), intent(in)                      :: x(:,:)
  type(mpi_communicator), optional, intent(in) :: comm
  complex(wp)                                  :: mean_
  !local variables
  real(wp)                                     :: size_
  
  mean_ = sum(x * w)
  size_ = sum(w)

  if (present(comm)) call comm%mpisum(mean_)
  if (present(comm)) call comm%mpisum(size_)

  mean_ = mean_ / size_
end function mean_w_2d_c

function mean_w_2d_ax_c(w, z, ax, comm) result(mean_)
  complex(wp), intent(in)                      :: w(:)
  complex(wp), intent(in)                      :: z(:,:)
  integer, intent(in)                          :: ax
  type(mpi_communicator), optional, intent(in) :: comm
  complex(wp), allocatable                     :: mean_(:)
  !local variables
  integer                                      :: i
  
  if (ax == 1) then
    if (.not. allocated(mean_)) allocate(mean_(size(z, 2)))
    do i = 1, size(mean_)
      mean_(i) = mean(w, z(:,i), comm)
    end do
  else if (ax == 2) then
    if (.not. allocated(mean_)) allocate(mean_(size(z, 1)))
    do i = 1, size(mean_)
      mean_(i) = mean(w, z(i,:), comm)
    end do
  else
    write(*,*) "wrong ax value in mean_w_arrray_c - program stops now"
    stop
  end if
end function mean_w_2d_ax_c

function mean_w_2d_rc(w, x, comm) result(mean_)
  real(wp), intent(in)                         :: w(:,:)
  complex(wp), intent(in)                      :: x(:,:)
  type(mpi_communicator), optional, intent(in) :: comm
  complex(wp)                                  :: mean_
  !local variables
  real(wp)                                     :: size_
  
  mean_ = sum(x * w)
  size_ = sum(w)

  if (present(comm)) call comm%mpisum(mean_)
  if (present(comm)) call comm%mpisum(size_)

  mean_ = mean_ / size_
end function mean_w_2d_rc

function mean_w_2d_ax_rc(w, z, ax, comm) result(mean_)
  real(wp), intent(in)                         :: w(:)
  complex(wp), intent(in)                      :: z(:,:)
  integer, intent(in)                          :: ax
  type(mpi_communicator), optional, intent(in) :: comm
  complex(wp), allocatable                     :: mean_(:)
  !local variables
  integer                                      :: i
  
  if (ax == 1) then
    if (.not. allocated(mean_)) allocate(mean_(size(z, 2)))
    do i = 1, size(mean_)
      mean_(i) = mean(w, z(:,i), comm)
    end do
  else if (ax == 2) then
    if (.not. allocated(mean_)) allocate(mean_(size(z, 1)))
    do i = 1, size(mean_)
      mean_(i) = mean(w, z(i,:), comm)
    end do
  else
    write(*,*) "wrong ax value in mean_w_arrray_rc - program stops now"
    stop
  end if
end function mean_w_2d_ax_rc


!******************************************************************************** 
!
!       Calculates standard deviation of the array
!
!               std = sqrt( \sum_i (x_i - \bar x)^2 / (N-1) )
!
!               std = sqrt(\sum_i w_i (x_i - \bar x)^2) / sqrt(V1- V2/V1)
!
!       Factors for std and std_mean
!
!                       unweighted              weighted
!
!       std             1 / sqrt(N-1)           1 / sqrt(V1 - V2/V1)
!
!       std_mean        1 / sqrt(N)             sqrt(V2/V1^2)
!
!******************************************************************************** 
function std(x, comm) result(std_)
  real(wp), intent(in)                         :: x(:)
  type(mpi_communicator), optional, intent(in) :: comm
  real(wp)                                     :: std_
  !local variables
  integer                                      :: size_

  if (size(x) == 1) then
    std_ = 0.0_wp
    return
  end if
  
  std_ = sum((x-mean(x,comm))**2)
  size_ = size(x)

  if (present(comm)) call comm%mpisum(std_)
  if (present(comm)) call comm%mpisum(size_)

  std_ = sqrt(std_ / real(size_-1, wp))
end function std

function std_c(z, comm) result(std_)
  complex(wp), intent(in)                      :: z(:)
  type(mpi_communicator), optional, intent(in) :: comm
  real(wp)                                     :: std_
  !local varaibles       
  integer                                      :: size_

  if (size(z) == 1) then
    std_ = 0.0_wp
    return
  end if
  
  std_ = sum(abs(z-mean(z,comm))**2)
  size_ = size(z)

  if (present(comm)) call comm%mpisum(std_)
  if (present(comm)) call comm%mpisum(size_)

  std_ = sqrt(std_ / real(size_-1, wp))
end function std_c

function std_2d(x, comm) result(std_)
  real(wp), intent(in)                         :: x(:,:)
  type(mpi_communicator), optional, intent(in) :: comm
  !local variables
  real(wp)                                     :: std_
  integer                                      :: size_

  if (size(x) == 1) then
    std_ = 0.0_wp
    return
  end if
  
  std_ = sum((x-mean(x,comm))**2)
  size_ = size(x)

  if (present(comm)) call comm%mpisum(std_)
  if (present(comm)) call comm%mpisum(size_)

  std_ = sqrt(std_ / real(size_-1, wp))
end function std_2d

function std_2d_ax(x, ax, comm) result(std_)
  real(wp), intent(in)                         :: x(:,:)
  integer, intent(in)                          :: ax
  type(mpi_communicator), optional, intent(in) :: comm
  real(wp), allocatable                        :: std_(:)
  !local variables
  integer                                      :: i
  
  if (ax == 1) then
    if (.not. allocated(std_)) allocate(std_(size(x, 2)))
    do i = 1, size(std_)
      std_(i) = std(x(:,i), comm)
    end do
  else if (ax == 2) then
    if (.not. allocated(std_)) allocate(std_(size(x, 1)))
    do i = 1, size(std_)
      std_(i) = std(x(i,:), comm)
    end do
  else
    write(*,*) "wrong ax value in std_arrray - program stops now"
    stop
  end if
end function std_2d_ax

function std_2d_c(z, comm) result(std_)
  complex(wp), intent(in)                      :: z(:,:)
  type(mpi_communicator), optional, intent(in) :: comm
  real(wp)                                     :: std_
  !local variables
  integer                                      :: size_

  if (size(z) == 1) then
    std_ = 0.0_wp
    return
  end if
  
  std_ = sum((z-mean(z,comm))**2)
  size_ = size(z)

  if (present(comm)) call comm%mpisum(std_)
  if (present(comm)) call comm%mpisum(size_)

  std_ = sqrt(std_ / real(size_-1, wp))
end function std_2d_c

function std_2d_ax_c(z, ax, comm) result(std_)
  complex(wp), intent(in)                      :: z(:,:)
  integer, intent(in)                          :: ax
  type(mpi_communicator), optional, intent(in) :: comm
  real(wp), allocatable                        :: std_(:)
  !local variables
  integer                                      :: i
  
  if (ax == 1) then
    if (.not. allocated(std_)) allocate(std_(size(z, 2)))
    do i = 1, size(std_)
      std_(i) = std(z(:,i), comm)
    end do
  else if (ax == 2) then
    if (.not. allocated(std_)) allocate(std_(size(z, 1)))
    do i = 1, size(std_)
      std_(i) = std(z(i,:), comm)
    end do
  else
    write(*,*) "wrong ax value in std_array_c - program stops now"
    stop
  end if
end function std_2d_ax_c

function std_w(w, x, comm) result(std_)
  real(wp), intent(in)                         :: w(:)
  real(wp), intent(in)                         :: x(:)
  type(mpi_communicator), optional, intent(in) :: comm
  real(wp)                                     :: std_
  !loclal variables
  real(wp)                                     :: size_, size_2

  if (size(x) == 1) then
    std_ = 0.0_wp
    return
  end if
  
  std_ = sum(w * (x-mean(w,x,comm))**2)
  size_ = sum(w)
  size_2 = sum(w**2)

  if (present(comm)) call comm%mpisum(std_)
  if (present(comm)) call comm%mpisum(size_)
  if (present(comm)) call comm%mpisum(size_2)

  std_ = sqrt(std_ / size_) / sqrt(1.0_wp - size_2/size_**2)
end function std_w

function std_w_c(w, z, comm) result(std_)
  complex(wp), intent(in)                      :: w(:)
  complex(wp), intent(in)                      :: z(:)
  type(mpi_communicator), optional, intent(in) :: comm
  real(wp)                                     :: std_
  !loclal variables
  complex(wp)                                  :: std_c, size_, size_2

  if (size(z) == 1) then
    std_ = 0.0_wp
    return
  end if
  
  std_c = sum(w * abs(z-mean(w,z,comm))**2)
  size_ = sum(w)
  size_2 = sum(w**2)

  if (present(comm)) call comm%mpisum(std_c)
  if (present(comm)) call comm%mpisum(size_)
  if (present(comm)) call comm%mpisum(size_2)

  std_ = real(sqrt(std_c / size_)/sqrt(1.0_wp - size_2/size_**2), wp)
end function std_w_c

function std_w_rc(w, z, comm) result(std_)
  real(wp), intent(in)                         :: w(:)
  complex(wp), intent(in)                      :: z(:)
  type(mpi_communicator), optional, intent(in) :: comm
  real(wp)                                     :: std_
  !loclal variables
  real(wp)                                     :: size_, size_2

  if (size(z) == 1) then
    std_ = 0.0_wp
    return
  end if
  
  std_ = sum(w * abs(z-mean(w,z,comm))**2)
  size_ = sum(w)
  size_2 = sum(w**2)

  if (present(comm)) call comm%mpisum(std_)
  if (present(comm)) call comm%mpisum(size_)
  if (present(comm)) call comm%mpisum(size_2)

  std_ = sqrt(std_ / size_) / sqrt(1.0_wp - size_2/size_**2)
end function std_w_rc

function std_w_2d(w, x, comm) result(std_)
  real(wp), intent(in)                         :: w(:,:)
  real(wp), intent(in)                         :: x(:,:)
  type(mpi_communicator), optional, intent(in) :: comm
  real(wp)                                     :: std_
  !loclal variables
  real(wp)                                     :: size_, size_2

  if (size(x) == 1) then
    std_ = 0.0_wp
    return
  end if
  
  std_ = sum(w * (x-mean(w,x,comm))**2)
  size_ = sum(w)
  size_2 = sum(w**2)

  if (present(comm)) call comm%mpisum(std_)
  if (present(comm)) call comm%mpisum(size_)
  if (present(comm)) call comm%mpisum(size_2)

  std_ = sqrt(std_ / size_) / sqrt(1.0_wp - size_2/size_**2)
end function std_w_2d

function std_w_2d_ax(w, x, ax, comm) result(std_)
  real(wp), intent(in)                         :: w(:)
  real(wp), intent(in)                         :: x(:,:)
  integer, intent(in)                          :: ax
  type(mpi_communicator), optional, intent(in) :: comm
  real(wp), allocatable                        :: std_(:)
  !local variables
  integer                                      :: i, ax_
  
  if (ax == 1) then
    if (.not. allocated(std_)) allocate(std_(size(x, 2)))
    do i = 1, size(std_)
      std_(i) = std(w, x(:, i), comm)
    end do
  else if (ax == 2) then
    if (.not. allocated(std_)) allocate(std_(size(x, 1)))
    do i = 1, size(std_)
      std_(i) = std(w, x(i, :), comm)
    end do
  else
    write(*,*) "wrong ax value in std_w_array - program stops now"
    stop
  end if
end function std_w_2d_ax

function std_w_2d_c(w, z, comm) result(std_)
  complex(wp), intent(in)                      :: w(:,:)
  complex(wp), intent(in)                      :: z(:,:)
  type(mpi_communicator), optional, intent(in) :: comm
  real(wp)                                     :: std_
  !loclal variables
  complex(wp)                                  :: std_c, size_, size_2

  if (size(z) == 1) then
    std_ = 0.0_wp
    return
  end if
  
  std_c = sum(w * abs(z-mean(w,z,comm))**2)
  size_ = sum(w)
  size_2 = sum(w**2)

  if (present(comm)) call comm%mpisum(std_c)
  if (present(comm)) call comm%mpisum(size_)
  if (present(comm)) call comm%mpisum(size_2)

  std_ = real(sqrt(std_c / size_)/sqrt(1.0_wp - size_2/size_**2), wp)
end function std_w_2d_c

function std_w_2d_ax_c(w, z, ax, comm) result(std_)
  complex(wp), intent(in)                      :: w(:)
  complex(wp), intent(in)                      :: z(:,:)
  integer, intent(in)                          :: ax
  type(mpi_communicator), optional, intent(in) :: comm
  real(wp), allocatable                        :: std_(:)
  !local variables
  integer                                      :: i, ax_
  
  if (ax == 1) then
    if (.not. allocated(std_)) allocate(std_(size(z, 2)))
    do i = 1, size(std_)
      std_(i) = std(w, z(:, i), comm)
    end do
  else if (ax == 2) then
    if (.not. allocated(std_)) allocate(std_(size(z, 1)))
    do i = 1, size(std_)
      std_(i) = std(w, z(i, :), comm)
    end do
  else
    write(*,*) "wrong ax value in std_w_array_c - program stops now"
    stop
  end if
end function std_w_2d_ax_c

function std_w_2d_rc(w, z, comm) result(std_)
  real(wp), intent(in)                         :: w(:,:)
  complex(wp), intent(in)                      :: z(:,:)
  type(mpi_communicator), optional, intent(in) :: comm
  real(wp)                                     :: std_
  !loclal variables
  real(wp)                                     :: size_, size_2

  if (size(z) == 1) then
    std_ = 0.0_wp
    return
  end if
  
  std_ = sum(w * abs(z-mean(w,z,comm))**2)
  size_ = sum(w)
  size_2 = sum(w**2)

  if (present(comm)) call comm%mpisum(std_)
  if (present(comm)) call comm%mpisum(size_)
  if (present(comm)) call comm%mpisum(size_2)

  std_ = sqrt(std_ / size_) / sqrt(1.0_wp - size_2/size_**2)
end function std_w_2d_rc

function std_w_2d_ax_rc(w, z, ax, comm) result(std_)
  real(wp), intent(in)                         :: w(:)
  complex(wp), intent(in)                      :: z(:,:)
  integer, intent(in)                          :: ax
  type(mpi_communicator), optional, intent(in) :: comm
  real(wp), allocatable                        :: std_(:)
  !local variables
  integer                                      :: i
  
  if (ax == 1) then
    if (.not. allocated(std_)) allocate(std_(size(z, 2)))
    do i = 1, size(std_)
      std_(i) = std(w, z(:, i), comm)
    end do
  else if (ax == 2) then
    if (.not. allocated(std_)) allocate(std_(size(z, 1)))
    do i = 1, size(std_)
      std_(i) = std(w, z(i, :), comm)
    end do
  else
    write(*,*) "wrong ax value in std_w_array_rc - program stops now"
    stop
  end if
end function std_w_2d_ax_rc


!******************************************************************************** 
!
!       Standard deviation of the mean of the array x
!
!               std_mean = std / sqrt(N)
!
!               std_mean = std * sqrt(V2 / V1^2)
!
!******************************************************************************** 
function std_mean(x, comm) result(std_)
  real(wp), intent(in)                         :: x(:)
  type(mpi_communicator), optional, intent(in) :: comm
  real(wp)                                     :: std_
  !local variables
  integer                                      :: size_
  
  std_ = std(x, comm) 
  size_ = size(x)

  if (present(comm)) call comm%mpisum(size_)
  
  std_ = std_ / sqrt(real(size_, wp))
end function std_mean

function std_mean_c(z, comm) result(std_)
  complex(wp), intent(in)                      :: z(:)
  type(mpi_communicator), optional, intent(in) :: comm
  real(wp)                                     :: std_
  !local variables
  integer                                      :: size_
  
  std_ = std(z, comm) 
  size_ = size(z)

  if (present(comm)) call comm%mpisum(size_)
  
  std_ = std_ / sqrt(real(size_, wp))
end function std_mean_c

function std_mean_2d(x, comm) result(std_)
  real(wp), intent(in)                         :: x(:,:)
  type(mpi_communicator), optional, intent(in) :: comm
  real(wp)                                     :: std_
  !local variables
  integer                                      :: size_
  
  std_ = std(x, comm) 
  size_ = size(x)

  if (present(comm)) call comm%mpisum(size_)
  
  std_ = std_ / sqrt(real(size_, wp))
end function std_mean_2d

function std_mean_2d_ax(x, ax, comm) result(std_)
  real(wp), intent(in)                         :: x(:,:)
  integer, intent(in)                          :: ax
  type(mpi_communicator), optional, intent(in) :: comm
  real(wp), allocatable                        :: std_(:)
  !local variables
  integer                                      :: i
  
  if (ax == 1) then
    if (.not. allocated(std_)) allocate(std_(size(x, 2)))
    do i = 1, size(std_)
      std_(i) = std_mean(x(:,i), comm)
    end do
  else if (ax == 2) then
    if (.not. allocated(std_)) allocate(std_(size(x, 1)))
    do i = 1, size(std_)
      std_(i) = std_mean(x(i,:), comm)
    end do
  else
    write(*,*) "wrong ax value in std_mean_array - program stops now"
    stop
  end if
end function std_mean_2d_ax

function std_mean_2d_c(z, comm) result(std_)
  complex(wp), intent(in)                      :: z(:,:)
  type(mpi_communicator), optional, intent(in) :: comm
  real(wp)                                     :: std_
  !local variables
  integer                                      :: size_
  
  std_ = std(z, comm) 
  size_ = size(z)

  if (present(comm)) call comm%mpisum(size_)
  
  std_ = std_ / sqrt(real(size_, wp))
end function std_mean_2d_c

function std_mean_2d_ax_c(z, ax, comm) result(std_)
  complex(wp), intent(in)                      :: z(:,:)
  integer, intent(in)                          :: ax
  type(mpi_communicator), optional, intent(in) :: comm
  real(wp), allocatable                        :: std_(:)
  !local variables
  integer                                      :: i
  
  if (ax == 1) then
    if (.not. allocated(std_)) allocate(std_(size(z, 2)))
    do i = 1, size(std_)
      std_(i) = std_mean(z(:,i), comm)
    end do
  else if (ax == 2) then
    if (.not. allocated(std_)) allocate(std_(size(z, 1)))
    do i = 1, size(std_)
      std_(i) = std_mean(z(i,:), comm)
    end do
  else
    write(*,*) "wrong ax value in std_mean_array_c - program stops now"
    stop
  end if
end function std_mean_2d_ax_c
  
function std_mean_w(w, x, comm) result(std_)
  real(wp), intent(in)                         :: w(:)
  real(wp), intent(in)                         :: x(:)
  type(mpi_communicator), optional, intent(in) :: comm
  real(wp)                                     :: std_
  !local variables
  real(wp)                                     :: size_, size_2
 
  std_ = std(w, x, comm) 
  size_ = sum(w)
  size_2 = sum(w**2)

  if (present(comm)) call comm%mpisum(size_)
  if (present(comm)) call comm%mpisum(size_2)
  
  std_ = std_ * sqrt(size_2/size_**2)
end function std_mean_w

function std_mean_w_c(w, z, comm) result(std_)
  complex(wp), intent(in)                      :: w(:)
  complex(wp), intent(in)                      :: z(:)
  type(mpi_communicator), optional, intent(in) :: comm
  real(wp)                                     :: std_
  !local variables
  complex(wp)                                  :: size_, size_2
  
  std_ = std(w, z, comm) 
  size_ = sum(w)
  size_2 = sum(w**2)

  if (present(comm)) call comm%mpisum(size_)
  if (present(comm)) call comm%mpisum(size_2)
  
  std_ = std_ * real(sqrt(size_2/size_**2), wp)
end function std_mean_w_c

function std_mean_w_rc(w, z, comm) result(std_)
  real(wp), intent(in)                         :: w(:)
  complex(wp), intent(in)                      :: z(:)
  type(mpi_communicator), optional, intent(in) :: comm
  real(wp)                                     :: std_
  !local variables
  real(wp)                                     :: size_, size_2
  
  std_ = std(w, z, comm) 
  size_ = sum(w)
  size_2 = sum(w**2)

  if (present(comm)) call comm%mpisum(size_)
  if (present(comm)) call comm%mpisum(size_2)
  
  std_ = std_ * sqrt(size_2/size_**2)
end function std_mean_w_rc

function std_mean_w_2d(w, x, comm) result(std_)
  real(wp), intent(in)                         :: w(:,:)
  real(wp), intent(in)                         :: x(:,:)
  type(mpi_communicator), optional, intent(in) :: comm
  real(wp)                                     :: std_
  !local variables
  real(wp)                                     :: size_, size_2
 
  std_ = std(w, x, comm) 
  size_ = sum(w)
  size_2 = sum(w**2)

  if (present(comm)) call comm%mpisum(size_)
  if (present(comm)) call comm%mpisum(size_2)
  
  std_ = std_ * sqrt(size_2/size_**2)
end function std_mean_w_2d

function std_mean_w_2d_ax(w, x, ax, comm) result(std_)
  real(wp), intent(in)                         :: w(:)
  real(wp), intent(in)                         :: x(:,:)
  integer, intent(in)                          :: ax
  type(mpi_communicator), optional, intent(in) :: comm
  real(wp), allocatable                        :: std_(:)
  !local variables
  integer                                      :: i
  
  if (ax == 1) then
    if (.not. allocated(std_)) allocate(std_(size(x, 2)))
    do i = 1, size(std_)
      std_(i) = std_mean(w, x(:,i), comm)
    end do
  else if (ax == 2) then
    if (.not. allocated(std_)) allocate(std_(size(x, 1)))
    do i = 1, size(std_)
      std_(i) = std_mean(w, x(i,:), comm)
    end do
  else
    write(*,*) "wrong ax value in std_mean_w_array - program stops now"
    stop
  end if
end function std_mean_w_2d_ax

function std_mean_w_2d_c(w, z, comm) result(std_)
  complex(wp), intent(in)                      :: w(:,:)
  complex(wp), intent(in)                      :: z(:,:)
  type(mpi_communicator), optional, intent(in) :: comm
  real(wp)                                     :: std_
  !local variables
  complex(wp)                                  :: size_, size_2
  
  std_ = std(w, z, comm) 
  size_ = sum(w)
  size_2 = sum(w**2)

  if (present(comm)) call comm%mpisum(size_)
  if (present(comm)) call comm%mpisum(size_2)
  
  std_ = std_ * real(sqrt(size_2/size_**2), wp)
end function std_mean_w_2d_c

function std_mean_w_2d_ax_c(w, z, ax, comm) result(std_)
  complex(wp), intent(in)                      :: w(:)
  complex(wp), intent(in)                      :: z(:,:)
  integer, intent(in)                          :: ax
  type(mpi_communicator), optional, intent(in) :: comm
  real(wp), allocatable                        :: std_(:)
  !local variables
  integer                                      :: i
  
  if (ax == 1) then
    if (.not. allocated(std_)) allocate(std_(size(z, 2)))
    do i = 1, size(std_)
      std_(i) = std_mean(w, z(:,i), comm)
    end do
  else if (ax == 2) then
    if (.not. allocated(std_)) allocate(std_(size(z, 1)))
    do i = 1, size(std_)
      std_(i) = std_mean(w, z(i,:), comm)
    end do
  else
    write(*,*) "wrong ax value in std_mean_w_array_c - program stops now"
    stop
  end if
end function std_mean_w_2d_ax_c

function std_mean_w_2d_rc(w, z, comm) result(std_)
  real(wp), intent(in)                         :: w(:,:)
  complex(wp), intent(in)                      :: z(:,:)
  type(mpi_communicator), optional, intent(in) :: comm
  real(wp)                                     :: std_
  !local variables
  real(wp)                                     :: size_, size_2
  
  std_ = std(w, z, comm) 
  size_ = sum(w)
  size_2 = sum(w**2)

  if (present(comm)) call comm%mpisum(size_)
  if (present(comm)) call comm%mpisum(size_2)
  
  std_ = std_ * sqrt(size_2/size_**2)
end function std_mean_w_2d_rc

function std_mean_w_2d_ax_rc(w, z, ax, comm) result(std_)
  real(wp), intent(in)                         :: w(:)
  complex(wp), intent(in)                      :: z(:,:)
  integer, intent(in)                          :: ax
  type(mpi_communicator), optional, intent(in) :: comm
  real(wp), allocatable                        :: std_(:)
  !local variables
  integer                                      :: i
  
  if (ax == 1) then
    if (.not. allocated(std_)) allocate(std_(size(z, 2)))
    do i = 1, size(std_)
      std_(i) = std_mean(w, z(:,i), comm)
    end do
  else if (ax == 2) then
    if (.not. allocated(std_)) allocate(std_(size(z, 1)))
    do i = 1, size(std_)
      std_(i) = std_mean(w, z(i,:), comm)
    end do
  else
    write(*,*) "wrong ax value in std_mean_w_array_rc - program stops now"
    stop
  end if
end function std_mean_w_2d_ax_rc


!******************************************************************************** 
!
!       Covariance of the samples x and y (optionally weighted by w)
!
!               Cov(x, y) = (\sum_i (x_i - mx) (y_u - my)) / (N-1) 
!
!               Cov(w, x, y) = (\sum_i w_i (x_i-mx) (y_i-my)) / (V_1 - V_2/V_1)
!
!******************************************************************************** 
function cov(x, y, comm) result(cov_)
  real(wp), intent(in)                         :: x(:)
  real(wp), intent(in)                         :: y(:)
  type(mpi_communicator), optional, intent(in) :: comm
  real(wp)                                     :: cov_
  !local variables
  integer                                      :: size_
  
  if (size(x) == 1) then
    cov_ = 0.0_wp
    return
  end if

  cov_ = sum((x-mean(x,comm)) * (y-mean(y,comm)))
  size_ = size(x)

  if (present(comm)) call comm%mpisum(cov_)
  if (present(comm)) call comm%mpisum(size_)

  cov_ = sqrt(cov_ / real(size_-1, wp))
end function cov

function cov_c(x, y, comm) result(cov_)
  complex(wp), intent(in)                      :: x(:)
  complex(wp), intent(in)                      :: y(:)
  type(mpi_communicator), optional, intent(in) :: comm
  complex(wp)                                  :: cov_
  !local variables
  integer                                      :: size_
  
  if (size(x) == 1) then
    cov_ = 0.0_wp
    return
  end if
  
  cov_ = sum(conjg(x-mean(x,comm)) * (y-mean(y,comm)))
  size_ = size(x)

  if (present(comm)) call comm%mpisum(cov_)
  if (present(comm)) call comm%mpisum(size_)
  
  cov_ = sqrt(cov_ / real(size_-1, wp))
end function cov_c

function cov_2d(x, y, comm) result(cov_)
  real(wp), intent(in)                         :: x(:,:)
  real(wp), intent(in)                         :: y(:,:)
  type(mpi_communicator), optional, intent(in) :: comm
  real(wp)                                     :: cov_
  !local variables
  integer                                      :: size_
  
  if (size(x) == 1) then
    cov_ = 0.0_wp
    return
  end if

  cov_ = sum((x-mean(x,comm)) * (y-mean(y,comm)))
  size_ = size(x)

  if (present(comm)) call comm%mpisum(cov_)
  if (present(comm)) call comm%mpisum(size_)

  cov_ = sqrt(cov_ / real(size_-1, wp))
end function cov_2d

function cov_2d_ax(x, ax, comm) result(cov_)
  real(wp), intent(in)                         :: x(:,:)
  integer, intent(in)                :: ax
  type(mpi_communicator), optional, intent(in) :: comm
  real(wp), allocatable                        :: cov_(:,:)
  !local variables
  integer                                      :: i, j
  
  if (ax == 1) then
    if (.not. allocated(cov_)) allocate(cov_(size(x,2), size(x,2)))
    do i = 1, size(cov_, 1)
      do j = 1, i
        cov_(j,i) = cov(x(:,j), x(:,i), comm)
        cov_(i,j) = cov_(j,i)
      end do
    end do
  else if (ax == 2) then
    if (.not. allocated(cov_)) allocate(cov_(size(x,1), size(x,1)))
    do i = 1, size(cov_, 1)
      do j = 1, i
        cov_(j,i) = cov(x(j,:), x(i,:), comm)
        cov_(j,i) = cov_(i,j)
      end do
    end do
  else
    write(*,*) "wrong ax value in cov_arrray - program stops now"
    stop
  end if
end function cov_2d_ax

function cov_2d_c(x, y, comm) result(cov_)
  complex(wp), intent(in)                      :: x(:,:)
  complex(wp), intent(in)                      :: y(:,:)
  type(mpi_communicator), optional, intent(in) :: comm
  complex(wp)                                  :: cov_
  !local variables
  integer                                      :: size_
  
  if (size(x) == 1) then
    cov_ = 0.0_wp
    return
  end if
  
  cov_ = sum(conjg(x-mean(x,comm)) * (y-mean(y,comm)))
  size_ = size(x)

  if (present(comm)) call comm%mpisum(cov_)
  if (present(comm)) call comm%mpisum(size_)
  
  cov_ = sqrt(cov_ / real(size_-1, wp))
end function cov_2d_c

function cov_2d_ax_c(x, ax, comm) result(cov_)
  complex(wp), intent(in)                      :: x(:,:)
  integer, intent(in)                          :: ax
  type(mpi_communicator), optional, intent(in) :: comm
  complex(wp), allocatable                     :: cov_(:,:)
  !local variables
  integer                                      :: i, j
  
  if (ax == 1) then
    if (.not. allocated(cov_)) allocate(cov_(size(x,2), size(x,2)))
    do i = 1, size(cov_)
      do j = 1, i
        cov_(j,i) = cov(x(:,j), x(:,i), comm)
        cov_(i,j) = conjg(cov_(j,i))
      end do
    end do
  else if (ax == 2) then
    if (.not. allocated(cov_)) allocate(cov_(size(x,1), size(x,1)))
    do i = 1, size(cov_)
      do j = 1, i
        cov_(j,i) = cov(x(j,:), x(i,:), comm)
        cov_(i,j) = conjg(cov_(j,i))
      end do
    end do
  else
    write(*,*) "wrong ax value in cov_array_c - program stops now"
    stop
  end if
end function cov_2d_ax_c

function cov_w(w, x, y, comm) result(cov_)
  real(wp), intent(in)                         :: w(:)
  real(wp), intent(in)                         :: x(:)
  real(wp), intent(in)                         :: y(:)
  type(mpi_communicator), optional, intent(in) :: comm
  real(wp)                                     :: cov_
  !loclal variables
  real(wp)                                     :: size_, size_2
  
  if (size(x) ==  1) then
    cov_ = 0.0_wp
    return
  end if 

  cov_ = sum(w * (x-mean(w,x,comm)) * (y-mean(w,y,comm)))
  size_ = sum(w)
  size_2 = sum(w**2)

  if (present(comm)) call comm%mpisum(cov_)
  if (present(comm)) call comm%mpisum(size_)
  if (present(comm)) call comm%mpisum(size_2)

  cov_ = cov_ / (size_ - size_2/size_)
end function cov_w

function cov_w_c(w, x, y, comm) result(cov_)
  complex(wp), intent(in)                      :: w(:)
  complex(wp), intent(in)                      :: x(:)
  complex(wp), intent(in)                      :: y(:)
  type(mpi_communicator), optional, intent(in) :: comm
  complex(wp)                                  :: cov_
  !loclal variables
  complex(wp)                                  :: size_, size_2

  if (size(x) ==  1) then
    cov_ = 0.0_wp
    return
  end if 
  
  cov_ = sum(w * conjg(x-mean(w,x,comm)) * (y-mean(w,y,comm)))
  size_ = sum(w)
  size_2 = sum(w**2)

  if (present(comm)) call comm%mpisum(cov_)
  if (present(comm)) call comm%mpisum(size_)
  if (present(comm)) call comm%mpisum(size_2)

  cov_ = cov_ / (size_ - size_2/size_)
end function cov_w_c

function cov_w_rc(w, x, y, comm) result(cov_)
  real(wp), intent(in)                         :: w(:)
  complex(wp), intent(in)                      :: x(:)
  complex(wp), intent(in)                      :: y(:)
  type(mpi_communicator), optional, intent(in) :: comm
  complex(wp)                                  :: cov_
  !loclal variables
  real(wp)                                     :: size_, size_2

  if (size(x) ==  1) then
    cov_ = 0.0_wp
    return
  end if 
  
  cov_ = sum(w * conjg(x-mean(w,x,comm)) * (y-mean(w,y,comm)))
  size_ = sum(w)
  size_2 = sum(w**2)

  if (present(comm)) call comm%mpisum(cov_)
  if (present(comm)) call comm%mpisum(size_)
  if (present(comm)) call comm%mpisum(size_2)

  cov_ = cov_ / (size_ - size_2/size_)
end function cov_w_rc

function cov_w_2d(w, x, y, comm) result(cov_)
  real(wp), intent(in)                         :: w(:,:)
  real(wp), intent(in)                         :: x(:,:)
  real(wp), intent(in)                         :: y(:,:)
  type(mpi_communicator), optional, intent(in) :: comm
  real(wp)                                     :: cov_
  !loclal variables
  real(wp)                                     :: size_, size_2
  
  if (size(x) ==  1) then
    cov_ = 0.0_wp
    return
  end if 

  cov_ = sum(w * (x-mean(w,x,comm)) * (y-mean(w,y,comm)))
  size_ = sum(w)
  size_2 = sum(w**2)

  if (present(comm)) call comm%mpisum(cov_)
  if (present(comm)) call comm%mpisum(size_)
  if (present(comm)) call comm%mpisum(size_2)

  cov_ = cov_ / (size_ - size_2/size_)
end function cov_w_2d

function cov_w_2d_ax(w, x, ax, comm) result(cov_)
  real(wp), intent(in)                         :: w(:)
  real(wp), intent(in)                         :: x(:,:)
  integer, intent(in)                          :: ax
  type(mpi_communicator), optional, intent(in) :: comm
  real(wp), allocatable                        :: cov_(:,:)
  !local variables
  integer                                      :: i, j
  
  if (ax == 1) then
    if (.not. allocated(cov_)) allocate(cov_(size(x,2), size(x,2)))
    do i = 1, size(cov_, 1)
      do j = 1, i
        cov_(j,i) = cov(w, x(:,j), x(:,i), comm)
        cov_(i,j) = cov_(j,i)
      end do
    end do
  else if (ax == 2) then
    if (.not. allocated(cov_)) allocate(cov_(size(x,1), size(x,1)))
    do i = 1, size(cov_, 1)
      do j = 1, i
        cov_(j,i) = cov(w, x(j,:), x(i,:), comm)
        cov_(i,j) = cov_(j,i)
      end do
    end do
  else
    write(*,*) "wrong ax value in cov_w_array - program stops now"
    stop
  end if
end function cov_w_2d_ax

function cov_w_2d_c(w, x, y, comm) result(cov_)
  complex(wp), intent(in)                      :: w(:,:)
  complex(wp), intent(in)                      :: x(:,:)
  complex(wp), intent(in)                      :: y(:,:)
  type(mpi_communicator), optional, intent(in) :: comm
  complex(wp)                                  :: cov_
  !loclal variables
  complex(wp)                                  :: size_, size_2

  if (size(x) ==  1) then
    cov_ = 0.0_wp
    return
  end if 
  
  cov_ = sum(w * conjg(x-mean(w,x,comm)) * (y-mean(w,y,comm)))
  size_ = sum(w)
  size_2 = sum(w**2)

  if (present(comm)) call comm%mpisum(cov_)
  if (present(comm)) call comm%mpisum(size_)
  if (present(comm)) call comm%mpisum(size_2)

  cov_ = cov_ / (size_ - size_2/size_)
end function cov_w_2d_c

function cov_w_2d_ax_c(w, x, ax, comm) result(cov_)
  complex(wp), intent(in)                      :: w(:)
  complex(wp), intent(in)                      :: x(:,:)
  integer, intent(in)                          :: ax
  type(mpi_communicator), optional, intent(in) :: comm
  complex(wp), allocatable                     :: cov_(:,:)
  !local variables
  integer                                      :: i, j
  
  if (ax == 1) then
    if (.not. allocated(cov_)) allocate(cov_(size(x,2), size(x,2)))
    do i = 1, size(cov_, 1)
      do j = 1, size(cov_, 1)
        cov_(j,i) = cov(w, x(:,j), x(:,i), comm)
      end do
    end do
  else if (ax == 2) then
    if (.not. allocated(cov_)) allocate(cov_(size(x,1), size(x,1)))
    do i = 1, size(cov_, 1)
      do j = 1, size(cov_, 1)
        cov_(j,i) = cov(w, x(j,:), x(i,:), comm)
      end do
    end do
  else
    write(*,*) "wrong ax value in cov_w_array_c - program stops now"
    stop
  end if
end function cov_w_2d_ax_c

function cov_w_2d_rc(w, x, y, comm) result(cov_)
  real(wp), intent(in)                         :: w(:,:)
  complex(wp), intent(in)                      :: x(:,:)
  complex(wp), intent(in)                      :: y(:,:)
  type(mpi_communicator), optional, intent(in) :: comm
  complex(wp)                                  :: cov_
  !loclal variables
  real(wp)                                     :: size_, size_2

  if (size(x) ==  1) then
    cov_ = 0.0_wp
    return
  end if 
  
  cov_ = sum(w * conjg(x-mean(w,x,comm)) * (y-mean(w,y,comm)))
  size_ = sum(w)
  size_2 = sum(w**2)

  if (present(comm)) call comm%mpisum(cov_)
  if (present(comm)) call comm%mpisum(size_)
  if (present(comm)) call comm%mpisum(size_2)

  cov_ = cov_ / (size_ - size_2/size_)
end function cov_w_2d_rc

function cov_w_2d_ax_rc(w, x, ax, comm) result(cov_)
  real(wp), intent(in)                         :: w(:)
  complex(wp), intent(in)                      :: x(:,:)
  integer, intent(in)                          :: ax
  type(mpi_communicator), optional, intent(in) :: comm
  complex(wp), allocatable                     :: cov_(:,:)
  !local variables
  integer                                      :: i, j
  
  if (ax == 1) then
    if (.not. allocated(cov_)) allocate(cov_(size(x,2), size(x,2)))
    do i = 1, size(cov_, 1)
      do j = 1, i
        cov_(j,i) = cov(w, x(:,j), x(:,i), comm)
        cov_(i,j) = conjg(cov_(j,i))
      end do
    end do
  else if (ax == 2) then
    if (.not. allocated(cov_)) allocate(cov_(size(x,1), size(x,1)))
    do i = 1, size(cov_, 1)
      do j = 1, i
        cov_(j,i) = cov(w, x(j,:), x(i,:), comm)
        cov_(i,j) = conjg(cov_(j,i))
      end do
    end do
  else
    write(*,*) "wrong ax value in cov_w_array_rc - program stops now"
    stop
  end if
end function cov_w_2d_ax_rc


!******************************************************************************** 
!
! Auto-covariance of the array for a given lag time
!
!    \gamma_a = 1/(N-a) \sum_k (O_k - <O>) (O_{k+a} - <O>)
!
!******************************************************************************** 
function autocov_r(x, lag, weights) result(autocov_)
  real(wp), intent(in)           :: x(:)
  integer, intent(in)            :: lag
  real(wp), optional, intent(in) :: weights(:)
  real(wp)                       :: autocov_
  !local variables
  integer                        :: n
  real(wp)                       :: mu
  real(wp), allocatable          :: weights_(:)

  n = size(x)

  if (present(weights)) then
    allocate(weights_(size(weights))) 
    weights_ = weights
  else
    allocate(weights_(n))
    weights_ = oner
  end if

  mu =  mean(weights, x)

  autocov_ = sum(weights_(1:n-lag) * (x(1:n-lag)-mu) * (x(lag+1:n)-mu)) / sum(weights_(1:n-lag))
end function autocov_r

function autocov_c(x, lag, weights) result(autocov_)
  complex(wp), intent(in)           :: x(:)
  integer, intent(in)               :: lag
  complex(wp), optional, intent(in) :: weights(:)
  real(wp)                          :: autocov_
  !local variables
  integer                           :: n
  complex(wp)                       :: mu
  real(wp), allocatable             :: weights_(:)

  n = size(x)

  if (present(weights)) then
    allocate(weights_(size(weights))) 
    weights_ = real(weights,wp)
  else
    allocate(weights_(n))
    weights_ = oner
  end if

  mu =  mean(weights, x)

  autocov_ = sum(weights_(1:n-lag) * real(x(1:n-lag)-mu, wp) * real(x(lag+1:n)-mu, wp)) / sum(weights_(1:n-lag))
end function autocov_c


!******************************************************************************** 
!
! Estimate autocorrelation length  a given lag time
!
! Condition to stop summation:
!
!    lag >= \sum_{a=-lag, lag} \rho_a = 1 + 2 \sum_{a=1,lag} \rho_a
!
!    \rho_a = \Gamma_a / \Gamma_0
!
!******************************************************************************** 
function autocorr_len_r(x, weights) result(corr_len)
  real(wp), intent(in)           :: x(:)
  real(wp), optional, intent(in) :: weights(:)
  real(wp)                       :: corr_len
  !local variables
  integer                        :: n, lag, maxlag
  real(wp)                       :: sigma
  real(wp), allocatable          :: weights_(:)
  
  n = size(x)
  maxlag = size(x) / 2

  if (present(weights)) then
    allocate(weights_(size(weights))) 
    weights_ = weights
  else
    allocate(weights_(n))
    weights_ = oner
  end if

  sigma = autocov(x, 0, weights_)

  corr_len = 1.0_wp

  do lag = 1, maxlag
    corr_len = corr_len + 2.0_wp * autocov(x, lag, weights_) / sigma
    !if (lag >= 2.0_wp*corr_len-1.0_wp) then
    if (lag >= corr_len) then
      corr_len = max(1.0_wp, corr_len)
      exit
    end if
  end do
end function autocorr_len_r

function autocorr_len_c(x, weights) result(corr_len)
  complex(wp), intent(in)           :: x(:)
  complex(wp), optional, intent(in) :: weights(:)
  real(wp)                          :: corr_len
  !local variables
  integer                           :: n, lag, maxlag
  complex(wp)                       :: sigma
  complex(wp), allocatable          :: weights_(:)
  
  n = size(x)
  maxlag = size(x) / 2

  if (present(weights)) then
    allocate(weights_(size(weights))) 
    weights_ = weights
  else
    allocate(weights_(n))
    weights_ = onec
  end if

  sigma = autocov(x, 0, weights_)

  corr_len = 1.0_wp

  do lag = 1, maxlag
    corr_len = corr_len + 2.0_wp * autocov(x, lag, weights_) / sigma
    !if (lag >= 2.0_wp*corr_len-1.0_wp) then 
    if (lag >= corr_len) then
      corr_len = max(1.0_wp, corr_len)
      exit
    end if
  end do
end function autocorr_len_c


!******************************************************************************** 
!
! Determine the block sizes for the blocking analysis
!
! Use Flyvbjerg-Petersen analysis:
!   block sizes are chosen as 2^i, for i=0,1,...,nmax-1
!
! Criterion for nmax:
!    2^{nmax - 1} < n/4
!
!******************************************************************************** 
subroutine get_block_values(n, blocks, bsize)
  integer, intent(in)                  :: n
  integer, allocatable, intent(inout)  :: blocks(:)
  integer, optional, intent(in)        :: bsize
  !local varaibles
  integer                              :: num_diff_blocks, i

  if (present(bsize)) then
    num_diff_blocks = 2
  else 
    num_diff_blocks = 1 + max(0, floor(log(n/4.0_wp)/log(2.0_wp)))
  end if

  allocate(blocks(num_diff_blocks))

  if (present(bsize)) then
    blocks(1) = 1
    blocks(2) = bsize
  else
    do i = 1, num_diff_blocks
      blocks(i) = 2**(i-1)
    end do
  end if
end subroutine get_block_values


!******************************************************************************** 
!
! Extract final std value from the array of stds for different block sizes
!
! Idea:
!    for perfectly correlated data sigma(i)/sigma(i-1) = sqrt(2)
!    defines threshold for stopping criterion:
!        sigma(i)/sigma(i-1) < thresh
!
! Good threshold value: (1 + sqrrt(2))/2
!
!******************************************************************************** 
function get_final_std(stds) result(std)
  real(wp), intent(in) :: stds(:)
  real(wp)             :: std
  !local variables
  real(wp), parameter  :: ratio = (1.0_wp+sqrt(2.0_wp))/2.0_wp
  integer              :: i

  do i = 2, size(stds)
    if (stds(i-1) == 0.0_wp) cycle
    std = stds(i)
    if (std/stds(i-1) < ratio) exit
  end do
end function get_final_std


!******************************************************************************** 
!
! Blocking analysis on array x
!
!******************************************************************************** 
subroutine blocking(x, std, corr_len, bsize, lwrite)
  real(wp), intent(in)          :: x(:)
  real(wp), intent(out)         :: std
  real(wp), intent(out)         :: corr_len
  integer, optional, intent(in) :: bsize
  logical, optional, intent(in) :: lwrite
  !local varaibles
  integer                       :: n, block_size, block_size_rest, num_blocks
  integer                       :: i, j, lower, upper
  real(wp), allocatable         :: means(:), stds(:), allmeans(:)
  integer, allocatable          :: blocks(:)
  
  n  = size(x)
  call get_block_values(n, blocks, bsize)

  allocate(allmeans(size(blocks)), stds(size(blocks)))
  allmeans = 0.0_wp
  stds = 0.0_wp

  do i = 1, size(blocks)
    block_size = blocks(i)
    num_blocks =  n / block_size
    block_size_rest = n - num_blocks * block_size
    
    if (block_size_rest == 0) then
      allocate(means(num_blocks))
    else
      allocate(means(num_blocks+1))
    end if
    
    do j = 1, num_blocks
      lower = block_size * (j-1) + 1
      upper = block_size * j
      means(j) = mean(x(lower:upper))
    end do
    
    if (block_size_rest /= 0) then
      means(num_blocks+1) = mean(x(upper+1:n))
    end if 
    
    allmeans(i) = mean(means)
    stds(i) = std_mean(means)
    
    if (allocated(means)) deallocate(means)
  end do
  
  std = get_final_std(stds)
  corr_len = (std/stds(1))**2
  
  if (present(lwrite))  call blocking_output(blocks, allmeans, stds)
end subroutine blocking

subroutine blocking_c(z, std, corr_len, bsize, lwrite)
  complex(wp), intent(in)       :: z(:)
  real(wp), intent(out)         :: std
  real(wp), intent(out)         :: corr_len
  integer, optional, intent(in) :: bsize
  logical, optional, intent(in) :: lwrite
  !local varaibles
  integer                       :: n, block_size, block_size_rest, num_blocks
  integer                       :: i, j, lower, upper
  complex(wp), allocatable      :: means(:), allmeans(:)
  real(wp), allocatable         :: stds(:)
  integer, allocatable          :: blocks(:)
  
  n  = size(z)
  call get_block_values(n, blocks, bsize)

  allocate(allmeans(size(blocks)), stds(size(blocks)))
  allmeans = 0.0_wp
  stds = 0.0_wp
  
  do i = 1, size(blocks)
    block_size = blocks(i)
    num_blocks =  n / block_size
    block_size_rest = n - num_blocks * block_size
    
    if (block_size_rest == 0) then
      allocate(means(num_blocks))
    else
      allocate(means(num_blocks+1))
    end if
    
    do j = 1, num_blocks
      lower = block_size * (j-1) + 1
      upper = block_size * j
      means(j) = mean(z(lower:upper))
    end do
    
    if (block_size_rest /= 0) then
      means(num_blocks+1) = mean(z(upper+1:n))
    end if 
    
    allmeans(i) = mean(means)
    stds(i) = std_mean(means)
    
    if (allocated(means)) deallocate(means)
  end do
  
  std = get_final_std(stds)
  corr_len = (std/stds(1))**2
  
  if (present(lwrite))  call blocking_output(blocks, allmeans, stds)
end subroutine blocking_c

subroutine blocking_array(x, std, corr_len, bsize, lwrite)
  real(wp), intent(in)          :: x(:,:)
  real(wp), intent(out)         :: std(:)
  real(wp), intent(out)         :: corr_len(:)
  integer, optional, intent(in) :: bsize
  logical, optional, intent(in) :: lwrite
  !local variables
  integer                       :: i
  
  do i = 1, size(x, 1)
    call blocking(x(i,:), std(i), corr_len(i), bsize, lwrite)
  end do
end subroutine blocking_array

subroutine blocking_array_c(z, std, corr_len, bsize, lwrite)
  complex(wp), intent(in)       :: z(:,:)
  real(wp), intent(out)         :: std(:)
  real(wp), intent(out)         :: corr_len(:)
  integer, optional, intent(in) :: bsize
  logical, optional, intent(in) :: lwrite
  !local variables
  integer                       :: i
  
  do i = 1, size(z, 1)
    call blocking(z(i,:), std(i), corr_len(i), bsize, lwrite)
  end do
end subroutine blocking_array_c
 

!******************************************************************************** 
!
! Blocking procedure suited for correlated quotient calculation
! of the form (\sum w_i E_i) / (\sum w_i)
!
!******************************************************************************** 
subroutine blocking_w(w, x, std, corr_len, bsize, lwrite)
  real(wp), intent(in)          :: w(:) 
  real(wp), intent(in)          :: x(:)
  real(wp), intent(out)         :: std
  real(wp), intent(out)         :: corr_len
  integer, optional, intent(in) :: bsize
  logical, optional, intent(in) :: lwrite
  !local varaibles
  integer                       :: n, block_size, block_size_rest, num_blocks
  integer                       :: i, j, lower, upper
  real(wp), allocatable         :: weights(:), means(:), stds(:), allmeans(:)
  integer, allocatable          :: blocks(:)
  
  n  = size(x)
  call get_block_values(n, blocks, bsize)

  allocate(allmeans(size(blocks)), stds(size(blocks)))
  allmeans = 0.0_wp
  stds = 0.0_wp
  
  do i = 1, size(blocks)
    block_size = blocks(i)
    num_blocks =  n / block_size
    block_size_rest = n - num_blocks * block_size
    
    if (block_size_rest == 0) then
      allocate(weights(num_blocks), means(num_blocks))
    else
      allocate(weights(num_blocks+1), means(num_blocks+1))
    end if
    
    do j = 1, num_blocks
      lower = block_size * (j-1) + 1
      upper = block_size * j
      weights(j) = sum(w(lower:upper))
      means(j)   = mean(w(lower:upper), x(lower:upper))
    end do
    
    if (block_size_rest /= 0) then
      weights(num_blocks+1) = sum(w(upper+1:n))
      means(num_blocks+1)   = mean(w(upper+1:n), x(upper+1:n))
    end if 
    
    allmeans(i) = mean(weights, means)
    stds(i) = std_mean(weights, means)
    
    if (allocated(means)) deallocate(means)
    if (allocated(weights)) deallocate(weights)
  end do
  
  std = get_final_std(stds)
  corr_len = (std/stds(1))**2
  
  if (present(lwrite))  call blocking_output(blocks, allmeans, stds)
end subroutine blocking_w

subroutine blocking_w_c(w, z, std, corr_len, bsize, lwrite)
  complex(wp), intent(in)       :: w(:) 
  complex(wp), intent(in)       :: z(:)
  real(wp), intent(out)         :: std
  real(wp), intent(out)         :: corr_len
  integer, optional, intent(in) :: bsize
  logical, optional, intent(in) :: lwrite
  !local varaibles
  integer                       :: n, block_size, block_size_rest, num_blocks
  integer                       :: i, j, lower, upper
  complex(wp), allocatable      :: weights(:), means(:), allmeans(:)
  real(wp), allocatable         :: stds(:)
  integer, allocatable          :: blocks(:)
  
  n  = size(z)
  call get_block_values(n, blocks, bsize)

  allocate(allmeans(size(blocks)), stds(size(blocks)))
  allmeans = 0.0_wp
  stds = 0.0_wp
  
  do i = 1, size(blocks)
    block_size = blocks(i) 
    num_blocks =  n / block_size
    block_size_rest = n - num_blocks * block_size
    
    if (block_size_rest == 0) then
      allocate(weights(num_blocks), means(num_blocks))
    else
      allocate(weights(num_blocks+1), means(num_blocks+1))
    end if
    
    do j = 1, num_blocks
      lower = block_size * (j-1) + 1
      upper = block_size * j
      weights(j) = sum(w(lower:upper))
      means(j)   = mean(w(lower:upper), z(lower:upper))
    end do
    
    if (block_size_rest /= 0) then
      weights(num_blocks+1) = sum(w(upper+1:n))
      means(num_blocks+1)   = mean(w(upper+1:n), z(upper+1:n))
    end if 
    
    allmeans(i) = mean(weights, means)
    stds(i) = real(std_mean(weights, means), wp)
    
    if (allocated(means)) deallocate(means)
    if (allocated(weights)) deallocate(weights)
  end do
  
  std = get_final_std(stds) 
  corr_len = (std/stds(1))**2

  if (present(lwrite))  call blocking_output(blocks, allmeans, stds)
end subroutine blocking_w_c

subroutine blocking_w_rc(w, z, std, corr_len, bsize, lwrite)
  real(wp), intent(in)          :: w(:) 
  complex(wp), intent(in)       :: z(:)
  real(wp), intent(out)         :: std
  real(wp), intent(out)         :: corr_len
  integer, optional, intent(in) :: bsize
  logical, optional, intent(in) :: lwrite
  !local varaibles
  integer                       :: n, block_size, block_size_rest, num_blocks
  integer                       :: i, j, lower, upper
  complex(wp), allocatable      :: means(:), allmeans(:)
  real(wp), allocatable         :: weights(:), stds(:)
  integer, allocatable          :: blocks(:)
  
  n  = size(z)
  call get_block_values(n, blocks, bsize)

  allocate(allmeans(size(blocks)), stds(size(blocks)))
  allmeans = 0.0_wp
  stds = 0.0_wp
  
  do i = 1, size(blocks)
    block_size = blocks(i)
    num_blocks =  n / block_size
    block_size_rest = n - num_blocks * block_size
    
    if (block_size_rest == 0) then
      allocate(weights(num_blocks), means(num_blocks))
    else
      allocate(weights(num_blocks+1), means(num_blocks+1))
    end if
    
    do j = 1, num_blocks
      lower = block_size * (j-1) + 1
      upper = block_size * j
      weights(j) = sum(w(lower:upper))
      means(j)   = mean(w(lower:upper), z(lower:upper))
    end do
    
    if (block_size_rest /= 0) then
      weights(num_blocks+1) = sum(w(upper+1:n))
      means(num_blocks+1)   = mean(w(upper+1:n), z(upper+1:n))
    end if 
    
    allmeans(i) = mean(weights, means)
    stds(i) = real(std_mean(weights, means), wp)
    
    if (allocated(means)) deallocate(means)
    if (allocated(weights)) deallocate(weights)
  end do
  
  std = get_final_std(stds)
  corr_len = (std/stds(1))**2
  
  if (present(lwrite))  call blocking_output(blocks, allmeans, stds)
end subroutine blocking_w_rc

subroutine blocking_w_array(w, x, std, corr_len, bsize, lwrite)
  real(wp), intent(in)          :: w(:) 
  real(wp), intent(in)          :: x(:,:)
  real(wp), intent(out)         :: std(:)
  real(wp), intent(out)         :: corr_len(:)
  integer, optional, intent(in) :: bsize
  logical, optional, intent(in) :: lwrite
  !local varaibles
  integer                       :: i
  
  do i = 1, size(x, 1)
    call blocking(w, x(i,:), std(i), corr_len(i), bsize, lwrite)
  end do
end subroutine blocking_w_array

subroutine blocking_w_array_c(w, z, std, corr_len, bsize, lwrite)
  complex(wp), intent(in)       :: w(:) 
  complex(wp), intent(in)       :: z(:,:)
  real(wp), intent(out)         :: std(:)
  real(wp), intent(out)         :: corr_len(:)
  integer, optional, intent(in) :: bsize
  logical, optional, intent(in) :: lwrite
  !local varaibles
  integer                       :: i
  
  do i = 1, size(z, 1)
    call blocking(w, z(i,:), std(i), corr_len(i), bsize, lwrite)
  end do
end subroutine blocking_w_array_c

subroutine blocking_w_array_rc(w, z, std, corr_len, bsize, lwrite)
  real(wp), intent(in)          :: w(:) 
  complex(wp), intent(in)       :: z(:,:)
  real(wp), intent(out)         :: std(:)
  real(wp), intent(out)         :: corr_len(:)
  integer, optional, intent(in) :: bsize
  logical, optional, intent(in) :: lwrite
  !local varaibles
  integer                       :: i
  
  do i = 1, size(z, 1)
    call blocking(w, z(i,:), std(i), corr_len(i), bsize, lwrite)
  end do
end subroutine blocking_w_array_rc


!******************************************************************************** 
!
!       Output of the blocking analysis
!
!******************************************************************************** 
subroutine blocking_output_real(blocks, means, stds)
  integer, intent(in)  :: blocks(:)
  real(wp), intent(in) :: means(:)
  real(wp), intent(in) :: stds(:)
  !local variables
  integer              :: i
  
  write(io%qmcfort_out%funit,100) "   block size  ", "      mean     ", "      std      "
  write(io%qmcfort_out%funit,100) "   ----------  ", "      ----     ", "      ---      "

  do i = 1, size(blocks)
    if (stds(i) == 0.0_wp ) cycle
    write(io%qmcfort_out%funit,101) blocks(i), means(i), stds(i)
  end do

  100 format(1X,T5,A,T20,A,T35,A)
  101 format(1X,T5,I6,T20,F12.6,T35,F12.6)
end subroutine blocking_output_real


subroutine blocking_output_cmplx(blocks, means, stds)
  integer, intent(in)     :: blocks(:)
  complex(wp), intent(in) :: means(:)
  real(wp), intent(in)    :: stds(:)
  !local variables
  integer                 :: i
  
  write(io%qmcfort_out%funit,100) "   block size  ", "    Re mean     ", "    Im mean     ", "      std      "
  write(io%qmcfort_out%funit,100) "   ----------  ", "    -------     ", "    -------     ", "      ---      "

  do i = 1, size(blocks)
    if (stds(i) == 0.0_wp ) cycle
    write(io%qmcfort_out%funit,101) blocks(i), means(i), stds(i)
  end do

  100 format(1X,T5,A,T20,A,T35,A,T50,A)
  101 format(1X,T5,I6,T20,F12.6,T35,F12.6,T50,F12.6)
end subroutine blocking_output_cmplx


!******************************************************************************** 
!
!   Chauvenet's Criterion to find outliers in Gaussian distributed random
!       values x_i
!
!           1. Calculate mean and sigma using N data points
!           2. Let z_i =  |x_i - mean| / sigma
!           3. If N erfc(z_i/sqrt(2)) < 0.5 reject the datapoint
!           4. Repeat process untill all points pass the chceck using only 
!               non-rejected values
!
!******************************************************************************** 
subroutine chauvenet_criterion(w, x, nout, sig, comm)
  real(wp), intent(inout)                        :: w(:)
  real(wp), intent(in)                           :: x(:)
  integer, intent(inout)                         :: nout
  real(wp), optional, intent(in)                 :: sig
  type(mpi_communicator), optional, intent(in)   :: comm
  !local variables
  real(wp)                                       :: m, s, sig_, z(size(w))
  logical                                        :: mask(size(w))
  integer                                        :: nout_old, nout_new, nout_diff
  !
  sig_ = 3.0_wp
  if (present(sig)) sig_ = sig
  !
  nout_diff = 1
  nout_new = 0
  !
  do while (nout_diff > 0)
    m = mean(w, x, comm=comm)
    s = sig_ * std(w, x, comm=comm)
    z = abs(x - m) / s
    mask = size(w)*erfc(z/sqrt(2.0_wp)) < 0.5_wp
    where (mask) w = 0.0_wp
    nout_old = nout_new
    nout_new = count(mask)
    nout_diff = nout_new - nout_old
    if (present(comm)) call comm%mpisum(nout_diff)
    nout = nout + nout_diff
  end do
end subroutine chauvenet_criterion


!******************************************************************************** 
!
!   Histogram function
!
!******************************************************************************** 
subroutine histogram_1d(x, hist, xhist, xmin, xmax, nbins, ldens, comm)
  real(wp), intent(in)                            :: x(:)
  real(wp), allocatable, intent(inout)            :: hist(:)
  real(wp), allocatable, optional, intent(out)    :: xhist(:)
  real(wp), optional, intent(in)                  :: xmin
  real(wp), optional, intent(in)                  :: xmax
  integer, optional, intent(in)                   :: nbins
  logical, optional, intent(in)                   :: ldens
  type(mpi_communicator), optional, intent(inout) :: comm
  !local variables
  real(wp)                                        :: xmin_, xmax_, dx, xval
  integer                                         :: i, nbins_, ind
  logical                                         :: ldens_
  
  if (present(xmin)) then
    xmin_ = xmin
  else
    xmin_ = minval(x)
  end if
  
  if (present(xmax)) then
    xmax_ = xmax
  else
    xmax_ = maxval(x)
  end if
  
  if (present(comm)) call comm%reduce(xmin_, op=mpi_min)
  if (present(comm)) call comm%reduce(xmax_, op=mpi_max)
  
  if (present(nbins)) then
    nbins_ = nbins
  else
    nbins_ = 1000
  end if

  if (present(ldens)) then
    ldens_ = ldens
  else
    ldens_ = .false.
  end if

  dx = (xmax_ - xmin_) / nbins_
  
  if (.not. allocated(hist)) allocate(hist(nbins_))
  hist = 0.0_wp
  
  do i = 1, size(x)
    ind = ceiling((x(i) - xmin_) / dx) 
    if (ind > nbins_) ind = nbins_
    if (ind < 1) ind = 1
    hist(ind) = hist(ind) + 1.0_wp
  end do
  
  if (present(comm)) call comm%mpisum(hist, 0)
  hist = hist / sum(hist)
  if (ldens_) hist = hist / dx
  
  if (present(xhist)) then
    if (allocated(xhist)) deallocate(xhist)
    allocate(xhist(nbins_))
    xval = xmin_ + dx/2.0_wp
    do i = 1, size(xhist)
      xhist(i) = xval
      xval = xval + dx
    end do
  end if
end subroutine histogram_1d

subroutine histogram_w_1d(w, x, hist, xhist, xmin, xmax, nbins, ldens, comm)
  real(wp), intent(in)                            :: w(:)
  real(wp), intent(in)                            :: x(:)
  real(wp), allocatable, intent(inout)            :: hist(:)
  real(wp), allocatable, optional, intent(out)    :: xhist(:)
  real(wp), optional, intent(in)                  :: xmin
  real(wp), optional, intent(in)                  :: xmax
  integer, optional, intent(in)                   :: nbins
  logical, optional, intent(in)                   :: ldens
  type(mpi_communicator), optional, intent(inout) :: comm
  !local variables
  real(wp)                                        :: xmin_, xmax_, dx, xval
  integer                                         :: i, nbins_, ind
  logical                                         :: ldens_
  
  if (present(xmin)) then
    xmin_ = xmin
  else
    xmin_ = minval(x)
  end if
  
  if (present(xmax)) then
    xmax_ = xmax
  else
    xmax_ = maxval(x)
  end if
  
  if (present(comm)) call comm%reduce(xmin_, op=mpi_min)
  if (present(comm)) call comm%reduce(xmax_, op=mpi_max)
  
  if (present(nbins)) then
    nbins_ = nbins
  else
    nbins_ = 1000
  end if
 
  if (present(ldens)) then
    ldens_ = ldens
  else
    ldens_ = .false.
  end if

  dx = (xmax_ - xmin_) / nbins_
  
  if (.not. allocated(hist)) allocate(hist(nbins_))
  hist = 0.0_wp
  
  do i = 1, size(x)
    ind = ceiling((x(i) - xmin_) / dx) 
    if (ind > nbins_) ind = nbins_
    if (ind < 1) ind = 1
    hist(ind) = hist(ind) + w(i)
  end do
  
  if (present(comm)) call comm%mpisum(hist, 0)
  hist = hist / sum(hist)
  if (ldens_) hist = hist / dx
  
  if (present(xhist)) then
    if (allocated(xhist)) deallocate(xhist)
    allocate(xhist(nbins_))
    xval = xmin_ + dx/2.0_wp
    do i = 1, size(xhist)
      xhist(i) = xval
      xval = xval + dx
    end do
  end if
end subroutine histogram_w_1d

subroutine histogram_2d(x, hist, xhist, xmin, xmax, nbins, ldens, comm)
  real(wp), intent(in)                            :: x(:,:)
  real(wp), allocatable, intent(inout)            :: hist(:)
  real(wp), allocatable, optional, intent(out)    :: xhist(:)
  real(wp), optional, intent(in)                  :: xmin
  real(wp), optional, intent(in)                  :: xmax
  integer, optional, intent(in)                   :: nbins
  logical, optional, intent(in)                   :: ldens
  type(mpi_communicator), optional, intent(inout) :: comm
  !local variables
  real(wp)                                        :: xmin_, xmax_, dx, xval
  integer                                         :: i, nbins_, ind, j
  logical                                         :: ldens_
  
  if (present(xmin)) then
    xmin_ = xmin
  else
    xmin_ = minval(x)
  end if
  
  if (present(xmax)) then
    xmax_ = xmax
  else
    xmax_ = maxval(x)
  end if
  
  if (present(comm)) call comm%reduce(xmin_, op=mpi_min)
  if (present(comm)) call comm%reduce(xmax_, op=mpi_max)
  
  if (present(nbins)) then
    nbins_ = nbins
  else 
    nbins_ = 1000
  end if
  
  if (present(ldens)) then
    ldens_ = ldens
  else
    ldens_ = .false.
  end if

  dx = (xmax_ - xmin_) / nbins_
  
  if (.not. allocated(hist)) allocate(hist(nbins_))
  hist = 0.0_wp
  
  do j = 1, size(x, 2)
    do i = 1, size(x, 1)
      ind = ceiling((x(i,j) - xmin_) / dx) 
      if (ind > nbins_) ind = nbins_
      if (ind < 1) ind = 1
      hist(ind) = hist(ind) + 1.0_wp
    end do
  end do
  
  if (present(comm)) call comm%mpisum(hist, 0)
  hist = hist / sum(hist)
  if (ldens_) hist = hist / dx
  
  if (present(xhist)) then
    if (allocated(xhist)) deallocate(xhist)
    allocate(xhist(nbins_))
    xval = xmin_ + dx/2.0_wp
    do i = 1, size(xhist)
      xhist(i) = xval
      xval = xval + dx
    end do
  end if
end subroutine histogram_2d

subroutine histogram_w_2d(w, x, hist, xhist, xmin, xmax, nbins, ldens, comm)
  real(wp), intent(in)                            :: w(:,:)
  real(wp), intent(in)                            :: x(:,:)
  real(wp), allocatable, intent(inout)            :: hist(:)
  real(wp), allocatable, optional, intent(out)    :: xhist(:)
  real(wp), optional, intent(in)                  :: xmin
  real(wp), optional, intent(in)                  :: xmax
  integer, optional, intent(in)                   :: nbins
  logical, optional, intent(in)                   :: ldens
  type(mpi_communicator), optional, intent(inout) :: comm
  !local variables
  real(wp)                                        :: xmin_, xmax_, dx, xval
  integer                                         :: i, j, nbins_, ind
  logical                                         :: ldens_
  
  if (present(xmin)) then
    xmin_ = xmin
  else
    xmin_ = minval(x)
  end if
  
  if (present(xmax)) then
    xmax_ = xmax
  else
    xmax_ = maxval(x)
  end if
  
  if (present(comm)) call comm%reduce(xmin_, op=mpi_min)
  if (present(comm)) call comm%reduce(xmax_, op=mpi_max)
  
  if (present(nbins))  then
    nbins_ = nbins
  else
    nbins_ = 1000
  end if

  if (present(ldens)) then
    ldens_ = ldens
  else
    ldens_ = .false.
  end if

  dx = (xmax_ - xmin_) / nbins_
  
  if (.not. allocated(hist)) allocate(hist(nbins_))
  hist = 0.0_wp
  
  do j = 1, size(x, 2)
    do i = 1, size(x, 1)
      ind = ceiling((x(i,j) - xmin_) / dx) 
      if (ind > nbins_) ind = nbins_
      if (ind < 1) ind = 1
      hist(ind) = hist(ind) + w(i,j)
    end do
  end do
  
  if (present(comm)) call comm%mpisum(hist, 0)
  hist = hist / sum(hist)
  if (ldens_) hist = hist / dx
  
  if (present(xhist)) then
    if (allocated(xhist)) deallocate(xhist)
    allocate(xhist(nbins_))
    xval = xmin_ + dx/2.0_wp
    do i = 1, size(xhist)
      xhist(i) = xval
      xval = xval + dx
    end do
  end if
end subroutine histogram_w_2d


!******************************************************************************** 
!
!       Dump histogram to the txt file
!
!******************************************************************************** 
subroutine dump_histogram(hist, xhist, fname)
  use file_handle, only: FileHandle
  real(wp), intent(in)                   :: hist(:)
  real(wp), optional, intent(in)         :: xhist(:)
  character(len=*), optional, intent(in) :: fname
  !local variables
  integer                                :: i
  character(len=:), allocatable          :: fname_
  type(FileHandle)                       :: hist_file
  
  if (comm_world%mpirank /= 0) return
  
  if (present(fname)) then
    fname_ = trim(fname)
  else
    fname_ = "histogram"
  end if
  
  hist_file = FileHandle(fname_)
  call hist_file%open(action="write", status="replace")
  
  do i = 1, size(hist)
    if (present(xhist)) then
      write(hist_file%funit, "(2f14.6)") xhist(i), hist(i)
    else
      write(hist_file%funit, "(i8,f14.6)") i, hist(i)
    end if
  end do
  
  call hist_file%close()
end subroutine dump_histogram

end module statistics
