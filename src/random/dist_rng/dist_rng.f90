! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module dist_rng

use constants

use basic_rng, only: BasicRng

implicit none

private
public :: DistRng, init_drng

type DistRng
  class(BasicRng), pointer :: brng
contains
  generic, public    :: uniform => uniform_r, uniform_r_1d, uniform_r_2d, &
                                   uniform_c, uniform_c_1d, uniform_c_2d
  procedure, public  :: uniform_r_1d !extended types can overwrite it
  procedure, private :: uniform_r, uniform_r_2d, uniform_c, uniform_c_1d, uniform_c_2d

  generic, public    :: normal => box_muller_r, box_muller_r_1d, box_muller_r_2d, &
                                  box_muller_c, box_muller_c_1d, box_muller_c_2d

  generic, public    :: box_muller => box_muller_r, box_muller_r_1d, box_muller_r_2d, &
                                      box_muller_c, box_muller_c_1d, box_muller_c_2d
  procedure, public  :: box_muller_r_1d !extended types can overwrite it
  procedure, private :: box_muller_r, box_muller_r_2d, box_muller_c, box_muller_c_1d, box_muller_c_2d

  generic, public    :: clt => clt_r, clt_r_1d, clt_r_2d, &
                               clt_c, clt_c_1d, clt_c_2d
  procedure, public  :: clt_r_1d !extended types can overwrite it
  procedure, private :: clt_r, clt_r_2d, clt_c, clt_c_1d, clt_c_2d

  generic, public    :: exp => exp_r, exp_r_1d, exp_r_2d
  procedure, public  :: exp_r_1d !extended types can overwrite it
  procedure, private :: exp_r, exp_r_2d

  generic, public    :: exp_range => exp_range_r, exp_range_r_1d, exp_range_r_2d
  procedure, public  :: exp_range_r_1d !extended types can overwrite it
  procedure, private :: exp_range_r, exp_range_r_2d

  generic, public    :: bernoulli => bernoulli_r, bernoulli_r_1d, bernoulli_r_2d
  procedure, public  :: bernoulli_r_1d !extended types can overwrite it
  procedure, private :: bernoulli_r, bernoulli_r_2d
end type DistRng

interface DistRng
  module procedure finit_drng
end interface DistRng

contains

!******************************************************************************** 
!
! Initialization of the DistRng object
!
!******************************************************************************** 
subroutine init_drng(brng, self)
  class(BasicRng), target, intent(in) :: brng
  type(DistRng), intent(out)          :: self

  self%brng => brng
end subroutine init_drng


!******************************************************************************** 
!
! DistRng constructor
!
!******************************************************************************** 
function finit_drng(brng) result(self)
  class(BasicRng), target, intent(in) :: brng
  type(DistRng)                       :: self

  call init_drng(brng, self)
end function finit_drng


!******************************************************************************** 
!
! Generate uniformly distributed random numbers on the interval [a,b]
!
! uniform_r_1d contains the actual implementation; other routines delegate to it
!
!******************************************************************************** 
subroutine uniform_r(self, sample, a, b)
  class(DistRng), intent(inout)  :: self
  real(wp), intent(out)          :: sample
  real(wp), optional, intent(in) :: a
  real(wp), optional, intent(in) :: b
  !local
  real(wp)                       :: sample_(1)

  call self%uniform(sample_, a, b)
  sample = sample_(1)
end subroutine uniform_r

subroutine uniform_r_1d(self, sample, a, b)
  class(DistRng), intent(inout)  :: self
  real(wp), intent(out)          :: sample(:)
  real(wp), optional, intent(in) :: a
  real(wp), optional, intent(in) :: b
  
  call self%brng%rand(sample)
  if (present(b)) sample = (b - a) * sample
  if (present(a)) sample = a + sample
end subroutine uniform_r_1d

subroutine uniform_r_2d(self, sample, a, b)
  class(DistRng), intent(inout)             :: self
  real(wp), target, contiguous, intent(out) :: sample(:,:)
  real(wp), optional, intent(in)            :: a
  real(wp), optional, intent(in)            :: b
  !local
  real(wp), pointer                         :: sample_(:)

  sample_(1:size(sample)) => sample
  call self%uniform(sample_, a, b)
end subroutine uniform_r_2d

subroutine uniform_c(self, sample, a, b)
  class(DistRng), intent(inout)  :: self
  complex(wp), intent(out)       :: sample
  real(wp), optional, intent(in) :: a
  real(wp), optional, intent(in) :: b
  !local
  real(wp)                       :: sample_(2)
  
  call self%uniform(sample_, a, b)
  sample = cmplx(sample_(1), sample_(2), kind=wp)
end subroutine uniform_c

subroutine uniform_c_1d(self, sample, a, b)
  class(DistRng), intent(inout)  :: self
  complex(wp), intent(out)       :: sample(:)
  real(wp), optional, intent(in) :: a
  real(wp), optional, intent(in) :: b
  !local
  real(wp), allocatable          :: sample_r(:), sample_i(:)
  
  allocate(sample_r(size(sample)), sample_i(size(sample)))
  
  call self%uniform(sample_r, a, b)
  call self%uniform(sample_i, a, b)
  sample = cmplx(sample_r, sample_i, kind=wp)
end subroutine uniform_c_1d

subroutine uniform_c_2d(self, sample, a, b)
  class(DistRng), intent(inout)                :: self
  complex(wp), target, contiguous, intent(out) :: sample(:,:)
  real(wp), optional, intent(in)               :: a
  real(wp), optional, intent(in)               :: b
  !local
  complex(wp), pointer                         :: sample_(:)
  
  sample_(1:size(sample)) => sample
  call self%uniform(sample_, a, b)
end subroutine uniform_c_2d


!******************************************************************************** 
!
! Generate normally distributed random numbers using Box-Muller method
!
!    r = sqrt(-2 * log(urand1))
!    phi = 2pi * urand2
!
!    x = r * cos(phi)
!    y = r * sin(phi)
!
! box_muller_r_1d contains the actual implementation; other routines delegate to it
!
!******************************************************************************** 
subroutine box_muller_r(self, sample, mean, std) 
  class(DistRng), intent(inout)  :: self
  real(wp), intent(out)          :: sample
  real(wp), optional, intent(in) :: mean
  real(wp), optional, intent(in) :: std
  !local
  real(wp)                       :: sample_(1)
  
  call self%box_muller(sample_, mean, std)
  sample = sample_(1)
end subroutine box_muller_r

subroutine box_muller_r_1d(self, sample, mean, std)
  class(DistRng), intent(inout)  :: self
  real(wp), intent(out)          :: sample(:)
  real(wp), optional, intent(in) :: mean
  real(wp), optional, intent(in) :: std
  !local
  integer                        :: m, n
  real(wp), allocatable          :: r(:), phi(:)
  
  n = size(sample)
  m = (n + 1) / 2
  allocate(r(m), phi(m))

  call self%uniform(r)
  call self%uniform(phi)

  !avoid problems with log
  where (r <= 0.0_wp)
    r = epsilon(0.0_wp)
  end where

  r = sqrt(-2.0_wp * log(r))
  phi = 2.0_wp * pi * phi

  sample(1:m)    = r * cos(phi)
  sample(n-m+1:) = r * sin(phi)

  if (present(std))  sample = sample * std
  if (present(mean)) sample = sample + mean
end subroutine box_muller_r_1d

subroutine box_muller_r_2d(self, sample, mean, std)
  class(DistRng), intent(inout)             :: self
  real(wp), target, contiguous, intent(out) :: sample(:,:)
  real(wp), optional, intent(in)            :: mean
  real(wp), optional, intent(in)            :: std
  !local
  real(wp), pointer                         :: sample_(:)
  
  sample_(1:size(sample)) => sample
  call self%box_muller(sample_, mean, std)
end subroutine box_muller_r_2d

subroutine box_muller_c(self, sample, mean, std) 
  class(DistRng), intent(inout)  :: self
  complex(wp), intent(out)       :: sample
  real(wp), optional, intent(in) :: mean
  real(wp), optional, intent(in) :: std
  !local
  real(wp)                       :: sample_(2)
  
  call self%box_muller(sample_, mean, std)
  sample = cmplx(sample_(1), sample_(2), kind=wp)
end subroutine box_muller_c

subroutine box_muller_c_1d(self, sample, mean, std) 
  class(DistRng), intent(inout)  :: self
  complex(wp), intent(out)       :: sample(:)
  real(wp), optional, intent(in) :: mean
  real(wp), optional, intent(in) :: std
  !local
  real(wp), allocatable          :: sample_r(:), sample_i(:)

  allocate(sample_r(size(sample)), sample_i(size(sample)))

  call self%box_muller(sample_r, mean, std)
  call self%box_muller(sample_i, mean, std)

  sample = cmplx(sample_r, sample_i, kind=wp)
end subroutine box_muller_c_1d

subroutine box_muller_c_2d(self, sample, mean, std) 
  class(DistRng), intent(inout)                :: self
  complex(wp), target, contiguous, intent(out) :: sample(:,:)
  real(wp), optional, intent(in)               :: mean
  real(wp), optional, intent(in)               :: std
  !local
  complex(wp), pointer                         :: sample_(:)

  sample_(1:size(sample)) => sample
  call self%box_muller(sample_, mean, std)
end subroutine box_muller_c_2d


!******************************************************************************** 
!
! Generate normally distributed random numbers using central limit theorem
!
!    x = clt_std * (\sum_i urand_i - clt_mean)
!    clt_mean = length / 2
!    clt_std = sqrt(12 / length)
!
! clt_r_1d contains the actual implementation; other routines delegate to it
!
!******************************************************************************** 
subroutine clt_r(self, sample, mean, std, length)
  class(DistRng), intent(inout)  :: self
  real(wp), intent(out)          :: sample
  real(wp), optional, intent(in) :: mean
  real(wp), optional, intent(in) :: std
  integer, optional, intent(in)  :: length
  !local
  real(wp)                       :: sample_(1)
  
  call self%clt(sample_, mean, std, length)
  sample = sample_(1)
end subroutine clt_r

subroutine clt_r_1d(self, sample, mean, std, length)
  class(DistRng), intent(inout)  :: self
  real(wp), intent(out)          :: sample(:)
  real(wp), optional, intent(in) :: mean
  real(wp), optional, intent(in) :: std
  integer, optional, intent(in)  :: length
  !local
  integer                        :: length_, i, lo, up
  real(wp)                       :: clt_mean, clt_std
  real(wp), allocatable          :: x(:)
  
  length_ = 12
  if (present(length)) length_ = length

  clt_std = sqrt(12.0_wp / real(length_, wp))
  clt_mean = real(length_, wp) / 2.0_wp

  allocate(x(length_ * size(sample)))
  call self%uniform(x)

  do i = 1, size(sample)
    lo = (i-1)*length_ + 1
    up = i * length_
    sample(i) = sum(x(lo:up)) 
  end do

  sample = clt_std * (sample - clt_mean)

  if (present(std))  sample = sample * std
  if (present(mean)) sample = sample + mean
end subroutine clt_r_1d

subroutine clt_r_2d(self, sample, mean, std, length)
  class(DistRng), intent(inout)             :: self
  real(wp), target, contiguous, intent(out) :: sample(:,:)
  real(wp), optional, intent(in)            :: mean
  real(wp), optional, intent(in)            :: std
  integer, intent(in), optional             :: length
  !local
  real(wp), pointer                         :: sample_(:)
  
  sample_(1:size(sample)) => sample
  call self%clt(sample_, mean, std, length)
end subroutine clt_r_2d

subroutine clt_c(self, sample, mean, std, length)
  class(DistRng), intent(inout)  :: self
  complex(wp), intent(out)       :: sample
  real(wp), optional, intent(in) :: mean
  real(wp), optional, intent(in) :: std
  integer, optional, intent(in)  :: length
  !local
  real(wp)                       :: sample_(2)
  
  call self%clt(sample_, mean, std, length)
  sample = cmplx(sample_(1), sample_(2), kind=wp)
end subroutine clt_c

subroutine clt_c_1d(self, sample, mean, std, length)
  class(DistRng), intent(inout)  :: self
  complex(wp), intent(out)       :: sample(:)
  real(wp), optional, intent(in) :: mean
  real(wp), optional, intent(in) :: std
  integer, optional, intent(in)  :: length
  !local
  real(wp), allocatable          :: sample_r(:), sample_i(:)
  
  allocate(sample_r(size(sample)), sample_i(size(sample)))
  call self%clt(sample_r, mean, std, length)
  call self%clt(sample_i, mean, std, length)
  sample = cmplx(sample_r, sample_i, kind=wp)
end subroutine clt_c_1d

subroutine clt_c_2d(self, sample, mean, std, length)
  class(DistRng), intent(inout)                :: self
  complex(wp), target, contiguous, intent(out) :: sample(:,:)
  real(wp), optional, intent(in)               :: mean
  real(wp), optional, intent(in)               :: std
  integer, optional, intent(in)                :: length
  !local
  complex(wp), pointer                         :: sample_(:)

  sample_(1:size(sample)) => sample
  call self%clt(sample_, mean, std, length)
end subroutine clt_c_2d


!******************************************************************************** 
!
! Random numbers from the exponential distribution p(x) = lambda*exp(-lambda*x)
!
!    x = - log(1-u) / lambda
!
! exp_r_1d contains the actual implementation; other routines delegate to it
!
!******************************************************************************** 
subroutine exp_r(self, sample, lambda)
  class(DistRng), intent(inout) :: self
  real(wp), intent(out)         :: sample
  real(wp), intent(in)          :: lambda
  !local
  real(wp)                      :: sample_(1)

  call self%exp(sample_, lambda)
  sample = sample_(1)
end subroutine exp_r

subroutine exp_r_1d(self, sample, lambda)
  class(DistRng), intent(inout) :: self
  real(wp), intent(out)         :: sample(:)
  real(wp), intent(in)          :: lambda

  call self%uniform(sample)
  sample = -log(sample) / lambda
end subroutine exp_r_1d

subroutine exp_r_2d(self, sample, lambda)
  class(DistRng), intent(inout)             :: self
  real(wp), target, contiguous, intent(out) :: sample(:,:)
  real(wp), intent(in)                      :: lambda
  !local
  real(wp), pointer                         :: sample_(:)

  sample_(1:size(sample)) => sample
  call self%exp(sample_, lambda)
end subroutine exp_r_2d


!******************************************************************************** 
!
! Random numbers from the exponential distribution on the range [x_min, x_max]
!
!    x = - log[e^{-lambda*x_min} - (e^{-lambda*x_min} - e^{-lambda*x_max}) u ]
!
! exp_range_r_1d contains the actual implementation; other routines delegate to it
!
!******************************************************************************** 
subroutine exp_range_r(self, sample, lambda, x_min, x_max)
  class(DistRng), intent(inout) :: self
  real(wp), intent(inout)       :: sample
  real(wp), intent(in)          :: lambda
  real(wp), intent(in)          :: x_min
  real(wp), intent(in)          :: x_max
  !local
  real(wp)                      :: sample_(1)

  call self%exp_range(sample_, lambda, x_min, x_max)
  sample = sample_(1)
end subroutine exp_range_r

subroutine exp_range_r_1d(self, sample, lambda, x_min, x_max)
  class(DistRng), intent(inout) :: self
  real(wp), intent(out)         :: sample(:)
  real(wp), intent(in)          :: lambda
  real(wp), intent(in)          :: x_min
  real(wp), intent(in)          :: x_max
  !local
  real(wp)                      :: c_min, c_max

  c_min = exp(-lambda*x_min)
  c_max = exp(-lambda*x_max)

  call self%uniform(sample)
  sample = - log(c_min - (c_min - c_max)*sample) / lambda
end subroutine exp_range_r_1d

subroutine exp_range_r_2d(self, sample, lambda, x_min, x_max)
  class(DistRng), intent(inout)             :: self
  real(wp), target, contiguous, intent(out) :: sample(:,:)
  real(wp), intent(in)                      :: lambda
  real(wp), intent(in)                      :: x_min
  real(wp), intent(in)                      :: x_max
  !local
  real(wp), pointer                         :: sample_(:)

  sample_(1:size(sample)) => sample
  call self%exp_range(sample_, lambda, x_min, x_max)
end subroutine exp_range_r_2d


!******************************************************************************** 
!
! Generate Bernoulli distributed random numbers
!
!                       { x,  with probability p
!    bernoulli(x,y,p) = {
!                       { y,  with probability 1-p
!
! bernoulli_r_1d contains the actual implementation; other routines delegate to it
!
!******************************************************************************** 
subroutine bernoulli_r(self, sample, x, y, p)
  class(DistRng), intent(inout) :: self
  real(wp), intent(out)         :: sample
  real(wp), intent(in)          :: x
  real(wp), intent(in)          :: y
  real(wp), intent(in)          :: p
  !local
  real(wp)                      :: sample_(1)

  call self%bernoulli(sample_, x, y, p)
  sample = sample_(1)
end subroutine bernoulli_r

subroutine bernoulli_r_1d(self, sample, x, y, p)
  class(DistRng), intent(inout) :: self
  real(wp), intent(out)         :: sample(:)
  real(wp), intent(in)          :: x
  real(wp), intent(in)          :: y
  real(wp), intent(in)          :: p

  call self%uniform(sample)

  where (sample <= p)
    sample = x
  else where
    sample = y
  end where
end subroutine bernoulli_r_1d

subroutine bernoulli_r_2d(self, sample, x, y, p)
  class(DistRng), intent(inout)             :: self
  real(wp), target, contiguous, intent(out) :: sample(:,:)
  real(wp), intent(in)                      :: x
  real(wp), intent(in)                      :: y
  real(wp), intent(in)                      :: p
  !local
  real(wp), pointer                         :: sample_(:)

  sample_(1:size(sample)) => sample
  call self%bernoulli(sample_, x, y, p)
end subroutine bernoulli_r_2d

end module dist_rng
