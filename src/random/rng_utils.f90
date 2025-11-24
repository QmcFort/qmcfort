! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module rng_utils

use constants
use mpi, only: comm_world
use qmcfort_in, only: add_input

implicit none

public

!default seed for each of the brngs
integer(i8), parameter :: def_brng_seed = 167901165_i8

!resolution for creating floating point random numbers
real(sp), parameter :: sp32_resolution = 2.0_sp**(-32)
real(dp), parameter :: dp32_resolution = 2.0_dp**(-32)
real(sp), parameter :: sp48_resolution = 2.0_sp**(-48)
real(dp), parameter :: dp48_resolution = 2.0_dp**(-48)
real(sp), parameter :: sp64_resolution = 2.0_sp**(-64)
real(dp), parameter :: dp64_resolution = 2.0_dp**(-64)

contains

subroutine get_random_seed(seed)
  integer(i8), intent(out) :: seed
  !local
  integer                  :: unit, ios
  
  if (comm_world%mpirank == 0) then
    inquire(iolength=unit) seed
    open(newunit=unit, file="/dev/urandom", access="stream", &
         form="unformatted", status="old", action="read", iostat=ios)
  
    if (ios == 0) then
      read(unit) seed
      close(unit, iostat=ios)
    else
      call get_seed_from_time(seed)
    end if
  end if

  call comm_world%bcast(seed, 0)
end subroutine get_random_seed


subroutine get_seed_from_time(seed)
  integer(i8), intent(out) :: seed
  !local
  integer                  :: i, clock, rate, max
  integer                  :: values(8)

  call date_and_time(values=values)
  call system_clock(count=clock, count_rate=rate, count_max=max)

  seed = int(clock + values(8) + i*37, kind=i8) 
  seed = ieor(seed, int(values(mod(i-1,8)+1), kind=i8) * i * 13)
end subroutine get_seed_from_time

end module rng_utils