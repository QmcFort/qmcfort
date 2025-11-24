! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module drng_factories

use constants
use basic_rng, only: BasicRng
use dist_rng, only: DistRng, init_drng
#ifdef MKL
use vsl_basic_rng, only: VslBasicRng
use vsl_dist_rng, only: VslDistRng, init_vsl_drng
#endif

implicit none

public

!used to provide one drng for entire program
class(DistRng), allocatable :: drng_glb

contains

subroutine create_drng(brng, drng)
  class(BasicRng), intent(in)              :: brng
  class(DistRng), allocatable, intent(out) :: drng
  !local
  type(DistRng), allocatable               :: drng_t

  select type (brng)
#ifdef MKL
    type is (VslBasicRng)
      call create_vsl_drng(brng, drng)
#endif
    class default
      allocate(DistRng:: drng_t)     
      call init_drng(brng, drng_t)
      call move_alloc(drng_t, drng)
  end select
end subroutine create_drng

#ifdef MKL
subroutine create_vsl_drng(brng, drng)
  class(BasicRng), intent(in)              :: brng
  class(DistRng), allocatable, intent(out) :: drng
  !local
  type(VslDistRng), allocatable            :: vsl_drng

  allocate(VslDistRng:: vsl_drng)
  call init_vsl_drng(brng, vsl_drng)
  call move_alloc(vsl_drng, drng)
end subroutine create_vsl_drng
#endif

end module drng_factories