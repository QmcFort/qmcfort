! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module vsl_dist_rng

use constants
use mpi

use basic_rng, only: BasicRng
use dist_rng, only: DistRng

#ifdef MKL
use mkl_vsl
use mkl_vsl_type

use vsl_basic_rng, only: VslBasicRng

implicit none

private
public :: VslDistRng, init_vsl_drng

type, extends(DistRng) :: VslDistRng
contains
  !overwritten procedures
  procedure, public  :: uniform_r_1d
  procedure, public  :: box_muller_r_1d
  
  !new implementations only for VSL
  generic, public    :: icdf => icdf_r, icdf_r_1d, icdf_r_2d, &
                                icdf_c, icdf_c_1d, icdf_c_2d
  procedure, private :: icdf_r, icdf_r_1d, icdf_r_2d
  procedure, private :: icdf_c, icdf_c_1d, icdf_c_2d

  !actual interface to the VSL calls
  generic, private   :: vsl_uniform => vsl_uniform_r_1d_sp, vsl_uniform_r_1d_dp
  procedure, private :: vsl_uniform_r_1d_sp, vsl_uniform_r_1d_dp

  generic, private   :: vsl_gaussian => vsl_gaussian_r_1d_sp, vsl_gaussian_r_1d_dp
  procedure, private :: vsl_gaussian_r_1d_sp, vsl_gaussian_r_1d_dp
end type VslDistRng

interface VslDistRng
  module procedure finit_vsl_drng
end interface VslDistRng

contains

!******************************************************************************** 
!
! Initialization of the VslDistRng object
!
!******************************************************************************** 
subroutine init_vsl_drng(brng, self)
  class(BasicRng), target, intent(in) :: brng
  type(VslDistRng), intent(out)       :: self

  select type (brng)
    type is (VslBasicRng)
      self%brng => brng
    class default
!debug: rng: call bug here
  end select
end subroutine init_vsl_drng


!******************************************************************************** 
!
! VslDistRng constructor
!
!******************************************************************************** 
function finit_vsl_drng(brng) result(self)
  class(BasicRng), target, intent(in) :: brng
  type(VslDistRng)                    :: self

  call init_vsl_drng(brng, self)
end function finit_vsl_drng


!******************************************************************************** 
!
! Fill the comment later
! 
!******************************************************************************** 
subroutine uniform_r_1d(self, sample, a, b)
  class(VslDistRng), intent(inout) :: self
  real(wp), intent(out)            :: sample(:)
  real(wp), optional, intent(in)   :: a
  real(wp), optional, intent(in)   :: b

  call self%vsl_uniform(VSL_RNG_METHOD_UNIFORM_STD, sample, a, b)
end subroutine uniform_r_1d


!******************************************************************************** 
!
! Fill the comment later
!
!******************************************************************************** 
subroutine box_muller_r_1d(self, sample, mean, std)
  class(VslDistRng), intent(inout) :: self
  real(wp), intent(out)            :: sample(:)
  real(wp), intent(in), optional   :: mean
  real(wp), intent(in), optional   :: std
  
  call self%vsl_gaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2, sample, mean, std)
end subroutine box_muller_r_1d


!******************************************************************************** 
!
! Generate Gaussian random numbers using icdf method
!
!    icdf - inverse cumulativer distribution function
!
!******************************************************************************** 
subroutine icdf_r(self, sample, mean, std) 
  class(VslDistRng), intent(inout) :: self
  real(wp), intent(out)            :: sample
  real(wp), optional, intent(in)   :: mean
  real(wp), optional, intent(in)   :: std
  !local
  real(wp)                         :: sample_(1)

  call self%icdf(sample_, mean, std)
  sample = sample_(1)
end subroutine icdf_r

subroutine icdf_r_1d(self, sample, mean, std) 
  class(VslDistRng), intent(inout) :: self
  real(wp), intent(out)            :: sample(:)
  real(wp), optional, intent(in)   :: mean
  real(wp), optional, intent(in)   :: std

  call self%vsl_gaussian(vsl_rng_method_gaussian_icdf, sample, mean, std)
end subroutine icdf_r_1d

subroutine icdf_r_2d(self, sample, mean, std)
  class(VslDistRng), intent(inout)          :: self
  real(wp), target, contiguous, intent(out) :: sample(:,:)
  real(wp), optional, intent(in)            :: mean
  real(wp), optional, intent(in)            :: std
  !local
  real(wp), pointer                         :: sample_(:)

  sample_(1:size(sample)) => sample
  call self%icdf(sample_, mean, std)
end subroutine icdf_r_2d

subroutine icdf_c(self, sample, mean, std)
  class(VslDistRng), intent(inout) :: self
  complex(wp), intent(out)         :: sample
  real(wp), optional, intent(in)   :: mean
  real(wp), optional, intent(in)   :: std
  !local
  real(wp)                         :: sample_(2)

  call self%icdf(sample_, mean, std)
  sample = cmplx(sample_(1), sample_(2), kind=wp)
end subroutine icdf_c

subroutine icdf_c_1d(self, sample, mean, std)
  class(VslDistRng), intent(inout) :: self
  complex(wp), intent(out)         :: sample(:)
  real(wp), optional, intent(in)   :: mean
  real(wp), optional, intent(in)   :: std
  !local
  real(wp), allocatable            :: sample_r(:), sample_i(:)

  allocate(sample_r(size(sample)), sample_i(size(sample)))
  call self%icdf(sample_r, mean, std)
  call self%icdf(sample_i, mean, std)
  sample = cmplx(sample_r, sample_i, kind=wp)
end subroutine icdf_c_1d

subroutine icdf_c_2d(self, sample, mean, std)
  class(VslDistRng), intent(inout)             :: self
  complex(wp), target, contiguous, intent(out) :: sample(:,:)
  real(wp), optional, intent(in)               :: mean
  real(wp), optional, intent(in)               :: std
  !local
  complex(wp), pointer                         :: sample_(:)
  
  sample_(1:size(sample)) => sample
  call self%icdf(sample_, mean, std)
end subroutine icdf_c_2d


!******************************************************************************** 
!
! Direct calls to VSL uniform random numbers (needed to support both sp and dp)
!
!******************************************************************************** 
subroutine vsl_uniform_r_1d_sp(self, method, sample, a, b)
  class(VslDistRng), intent(inout) :: self
  integer, intent(in)              :: method
  real(sp), intent(out)            :: sample(:)
  real(sp), optional, intent(in)   :: a
  real(sp), optional, intent(in)   :: b
  !local
  integer                          :: ierr
  real(sp)                         :: a_, b_
  class(BasicRng), allocatable     :: brng_thread

  a_ = 0.0_sp
  if (present(a)) a_ = a

  b_ = 1.0_sp
  if (present(b)) b_ = b
  
  call self%brng%prepare_stream(size(sample), brng_thread)
  select type (brng_thread)
    type is (VslBasicRng)
      ierr = VsRngUniform(method, brng_thread%stream, size(sample), sample, a_, b_)
    class default
!debug: rng: bug here
  end select
  call self%brng%align_stream(size(sample))

  !ensure the range is (a,b] - originally, it is [a, b)
  sample = a_ + b_ - sample
end subroutine vsl_uniform_r_1d_sp

subroutine vsl_uniform_r_1d_dp(self, method, sample, a, b)
  class(VslDistRng), intent(inout) :: self
  integer, intent(in)              :: method
  real(dp), intent(out)            :: sample(:)
  real(dp), optional, intent(in)   :: a
  real(dp), optional, intent(in)   :: b
  !local
  integer                          :: ierr
  real(dp)                         :: a_, b_
  class(BasicRng), allocatable     :: brng_thread

  a_ = 0.0_dp
  if (present(a)) a_ = a

  b_ = 1.0_dp
  if (present(b)) b_ = b
  
  call self%brng%prepare_stream(size(sample), brng_thread)
  select type (brng_thread)
    type is (VslBasicRng)
      ierr = VdRngUniform(method, brng_thread%stream, size(sample), sample, a_, b_)
    class default
!debug: rng: bug here
  end select
  call self%brng%align_stream(size(sample))

  !ensure the range is (a,b] - originally, it is [a, b)
  sample = a_ + b_ - sample
end subroutine vsl_uniform_r_1d_dp


!******************************************************************************** 
!
! Direct calls to VSL Gaussian random numbers (needed to support both sp and dp)
!
!******************************************************************************** 
subroutine vsl_gaussian_r_1d_sp(self, method, sample, mean, std)
  class(VslDistRng), intent(inout) :: self
  integer, intent(in)              :: method
  real(sp), intent(out)            :: sample(:)
  real(sp), optional, intent(in)   :: mean
  real(sp), optional, intent(in)   :: std
  !local
  integer                          :: ierr
  real(sp)                         :: mean_, std_
  class(BasicRng), allocatable     :: brng_thread

  mean_ = 0.0_sp
  if (present(mean)) mean_ = mean

  std_ = 1.0_sp
  if (present(std)) std_ = std

  call self%brng%prepare_stream(size(sample), brng_thread)
  select type (brng_thread)
    type is (VslBasicRng)
      ierr = VsRngGaussian(method, brng_thread%stream, size(sample), sample, mean_, std_)
    class default
!debug: rng: bug here
  end select
  call self%brng%align_stream(size(sample))
end subroutine vsl_gaussian_r_1d_sp

subroutine vsl_gaussian_r_1d_dp(self, method, sample, mean, std)
  class(VslDistRng), intent(inout) :: self
  integer, intent(in)              :: method
  real(dp), intent(out)            :: sample(:)
  real(dp), optional, intent(in)   :: mean
  real(dp), optional, intent(in)   :: std
  !local
  integer                          :: ierr
  real(dp)                         :: mean_, std_
  class(BasicRng), allocatable     :: brng_thread

  mean_ = 0.0_dp
  if (present(mean)) mean_ = mean

  std_ = 1.0_dp
  if (present(std)) std_ = std
  
  call self%brng%prepare_stream(size(sample), brng_thread)
  select type (brng_thread)
    type is (VslBasicRng)
      ierr = VdRngGaussian(method, brng_thread%stream, size(sample), sample, mean_, std_)
    class default
!debug: rng: bug here
  end select
  call self%brng%align_stream(size(sample))
end subroutine vsl_gaussian_r_1d_dp

#endif
end module vsl_dist_rng