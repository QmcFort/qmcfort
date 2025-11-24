! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module regression

#include "preproc.inc"
use constants
use mpi
use file_handle, only: FileHandle
use statistics
use lapack 

implicit none 

interface report_polyfit
  module procedure report_polyfit_r, report_polyfit_c
end interface report_polyfit

contains

subroutine linear_fit(x, y, coeff, dcoeff, rsq, dy)
  real(wp), intent(in)           :: x(:)
  real(wp), intent(in)           :: y(:)
  real(wp), intent(out)          :: coeff(:)
  real(wp), intent(out)          :: dcoeff(:)
  real(wp), intent(out)          :: rsq 
  real(wp), optional, intent(in) :: dy(:)
  !local variables
  integer                        :: n, i
  real(wp)                       :: wmx, wmy, wmdx, sigma
  real(wp)                       :: sse, sst
  real(wp), allocatable          :: w(:), deps(:), xmat(:,:), xtmat(:,:), xtxmat(:,:)

  n = size(x)
  allocate(w(n), deps(n), xmat(n,2), xtmat(2,n), xtxmat(2,2))

  if (present(dy)) then
    w = 1.0_wp / dy**2
    w = w / sum(w)
  else 
    w = 1.0_wp
  end if

  xmat(:,1) =  x 
  xmat(:,2) = 1.0_wp

  xtmat(1,:) = xmat(:,1) * w
  xtmat(2,:) = xmat(:,2) * w

  xtxmat = matmul(xtmat, xmat)
  call inverse(xtxmat)

  coeff = matmul(xtxmat, matmul(xtmat, y))
  
  deps = y - (coeff(1)*x + coeff(2))
  wmx = mean(w, x)
  wmy = mean(w, y)
  wmdx = sum((x - wmx)**2)
  sigma = std(w, deps)

  dcoeff(1) = sigma / sqrt(wmdx)
  dcoeff(2) = sigma * sqrt(1.0_wp/n + wmx**2/wmdx)

  sse = sum(w * deps**2)
  sst = sum(w * (y - wmy)**2)
  rsq = 1.0_wp - sse/sst
end subroutine linear_fit


subroutine polyfit(x, y, order, coeff, dcoeff, rsq, dy)
  real(wp), intent(in)           :: x(:)
  real(wp), intent(in)           :: y(:)
  integer, intent(in)            :: order
  real(wp), intent(out)          :: coeff(:)
  real(wp), intent(out)          :: dcoeff(:)
  real(wp), intent(out)          :: rsq 
  real(wp), optional, intent(in) :: dy(:)
  !local variables
  integer                        :: n, i
  real(wp)                       :: wmx, wmy, wmdx, sigma
  real(wp)                       :: sse, sst
  real(wp), allocatable          :: w(:), x_(:), deps(:), xmat(:,:), xtmat(:,:), xtxmat(:,:)

  n = size(x)
  allocate(w(n), x_(n), deps(n), xmat(n,2), xtmat(2,n), xtxmat(2,2))

  if (present(dy)) then
    w = 1.0_wp / dy**2
    w = w / sum(w)
  else 
    w = 1.0_wp
  end if

  x_ = x**(order)

  xmat(:,1) =  x_
  xmat(:,2) = 1.0_wp

  xtmat(1,:) = xmat(:,1) * w
  xtmat(2,:) = xmat(:,2) * w

  xtxmat = matmul(xtmat, xmat)
  call inverse(xtxmat)

  coeff = matmul(xtxmat, matmul(xtmat, y))
  
  deps = y - (coeff(1)*x_ + coeff(2))
  wmx = mean(w, x_)
  wmy = mean(w, y)
  wmdx = sum((x_ - wmx)**2)
  sigma = std(w, deps)

  dcoeff(1) = sigma / sqrt(wmdx)
  dcoeff(2) = sigma * sqrt(1.0_wp/n + wmx**2/wmdx)

  sse = sum(w * deps**2)
  sst = sum(w * (y - wmy)**2)
  rsq = 1.0_wp - sse/sst
end subroutine polyfit


subroutine report_linear_fit(fh, x, y, coeff, dcoeff, rsq, dy)
  type(FileHandle), intent(in)   :: fh
  real(wp), intent(in)           :: x(:)
  real(wp), intent(in)           :: y(:)
  real(wp), intent(in)           :: coeff(:)
  real(wp), intent(in)           :: dcoeff(:)
  real(wp), intent(in)           :: rsq 
  real(wp), optional, intent(in) :: dy(:)
  !local variables
  integer                        :: i, n

  n = size(x)

  write(fh%funit,*) "******************** Regression summary ********************"
  do i = 1, n
    if (present(dy)) then
      write(fh%funit,100) i, x(i), y(i), dy(i)
    else
      write(fh%funit,100) i, x(i), y(i)
    end if
  end do
  write(fh%funit,*) "******************************************" 
  write(fh%funit,101) "E(t-->0)  =      ", coeff(2), dcoeff(2)
  write(fh%funit,101) "slope     =      ", coeff(1), dcoeff(1)
  write(fh%funit,101) "r-squared:       ", rsq
  write(fh%funit,*) "*******************************************************************"

  100 format (1x,i6,t10,es12.4,t25,f12.6,t40,f12.6)
  101 format (1x,a,t20,f12.6,t35,f12.6)
end subroutine report_linear_fit


subroutine report_polyfit_r(fh, x, y, order, coeff, dcoeff, rsq, dy)
  type(FileHandle), intent(in)   :: fh
  real(wp), intent(in)           :: x(:)
  real(wp), intent(in)           :: y(:)
  integer, intent(in)            :: order
  real(wp), intent(in)           :: coeff(:)
  real(wp), intent(in)           :: dcoeff(:)
  real(wp), intent(in)           :: rsq 
  real(wp), optional, intent(in) :: dy(:)
  !local variables
  integer                        :: i, n
  character(len=charlen)         :: order_

  n = size(x)
  write(order_,"(i4)") order

  write(fh%funit,*) "*************** Regression summary k*x^{", trim(adjustl(order_)), "} ****************"
  do i = 1, n
    if (present(dy)) then
      write(fh%funit,100) i, x(i), y(i), dy(i)
    else
      write(fh%funit,100) i, x(i), y(i)
    end if
  end do
  write(fh%funit,*) "******************************************" 
  write(fh%funit,101) "E(t-->0)  =      ", coeff(2), dcoeff(2)
  write(fh%funit,101) "slope     =      ", coeff(1), dcoeff(1)
  write(fh%funit,101) "r-squared:       ", rsq
  write(fh%funit,*) "*******************************************************************"

  100 format (1x,i6,t10,es12.4,t25,f12.6,t40,f12.6)
  101 format (1x,a,t20,f12.6,t35,f12.6)
end subroutine report_polyfit_r

subroutine report_polyfit_c(fh, x, y, order, coeff, dcoeff, rsq, dy)
  type(FileHandle), intent(in)   :: fh
  real(wp), intent(in)           :: x(:)
  complex(wp), intent(in)        :: y(:)
  integer, intent(in)            :: order
  real(wp), intent(in)           :: coeff(:)
  real(wp), intent(in)           :: dcoeff(:)
  real(wp), intent(in)           :: rsq 
  real(wp), optional, intent(in) :: dy(:)
  !local variables
  integer                        :: i, n
  character(len=charlen)         :: order_

  n = size(x)
  write(order_,"(i4)") order

  write(fh%funit,*) "*************** Regression summary k*x^{", trim(adjustl(order_)), "} ****************"
  do i = 1, n
    if (present(dy)) then
      write(fh%funit,100) i, x(i), real(y(i),wp), aimag(y(i)), dy(i)
    else
      write(fh%funit,100) i, x(i), real(y(i),wp), aimag(y(i))
    end if
  end do
  write(fh%funit,*) "******************************************" 
  write(fh%funit,101) "E(t-->0)  =      ", coeff(2), dcoeff(2)
  write(fh%funit,101) "slope     =      ", coeff(1), dcoeff(1)
  write(fh%funit,101) "r-squared:       ", rsq
  write(fh%funit,*) "*******************************************************************"

  100 format (1x,i6,t10,es12.4,t25,f12.6,t40,f12.6,t55,f12.6)
  101 format (1x,a,t20,f12.6,t35,f12.6)
end subroutine report_polyfit_c

end module regression