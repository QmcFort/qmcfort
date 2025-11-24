! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module boys_functions

use constants

implicit none

contains

!******************************************************************************** 
!
!       Boys function for x=0
!
!******************************************************************************** 
real(wp) function boysfun_zero(n)
  integer, intent(in) :: n
  
  boysfun_zero = 1.0_wp / (2.0_wp*n + 1.0_wp)
end function boysfun_zero


!******************************************************************************** 
!
!       Boys function of the zero order
!
!******************************************************************************** 
real(wp) function boysfun_zero_order(x)
  real(wp), intent(in) :: x
  
  boysfun_zero_order = 0.5_wp * sqrt(pi/x) * erf(sqrt(x))
end function boysfun_zero_order


!******************************************************************************** 
!
!       Boys function for the small argument
!
!******************************************************************************** 
real(wp) function boysfun_small(x, n, kmax, tol)
  real(wp), intent(in)           :: x
  integer, intent(in)            :: n
  integer, optional, intent(in)  :: kmax
  real(wp), optional, intent(in) :: tol
  !local variables
  real(wp), parameter            :: tolx = 1.0e-12_wp
  integer                        :: kmax_
  real(wp)                       :: tol_, x_, fact, df
  integer                        :: k
  !
  if (present(kmax)) then
    kmax_ = kmax
  else
    kmax_ = 100
  end if
  !
  if (present(tol)) then
    tol_ = tol
  else
    tol_ = 1.0e-10_wp
  end if
  !
  boysfun_small =  1.0_wp / (2.0_wp*n + 1.0_wp) 
  if (x >= tolx) then
    k = 0
    x_ = 1.0_wp
    fact = 1.0_wp
    df = 1000.0_wp
    do while (k<kmax_ .and. abs(df)>tol_)
      k = k + 1
      x_ = -x_ * x 
      fact = fact * real(k, wp)
      df = x_ / (fact * (2.0_wp*n + 2.0_wp*k + 1.0_wp))
      boysfun_small = boysfun_small + df
    end do
  end if
end function boysfun_small


!******************************************************************************** 
!
!       Boys function for the large argument
!
!******************************************************************************** 
real(wp) function boysfun_large(x, n)
  use standalone, only: double_factorial_rec
  real(wp), intent(in) :: x
  integer, intent(in)  :: n
  !
  boysfun_large = double_factorial_rec(2*n-1) / (2.0_wp**(n+1)) * sqrt(pi/x**(2*n+1))
end function boysfun_large


!******************************************************************************** 
!
!       Analytic Boys function
!
!******************************************************************************** 
real(wp) function boysfun_an(x, n) 
  real(wp), intent(in) :: x
  integer, intent(in)  :: n
  !local variables
  integer(kind=4)      :: ifault
  real(dp)             :: n_
  real(dp)             :: gammad
  !
  n_ = real(n, dp) + 0.5_dp
  boysfun_an = gamma(n_) * gammad(real(x, dp), n_, ifault) / (2.0_wp*x**n_)
end function boysfun_an


!******************************************************************************** 
!
!       Using upward recursion relation to compute boys function
!
!               F_{n+1} = ((2n+1) F_n - exp(-x)) / 2x
!
!       Use carefully: stable only for large x arguments
!
!******************************************************************************** 
real(wp) function boysfun_up(x, n)
  real(wp), intent(in) :: x
  integer, intent(in)  :: n
  !local variables
  integer              :: i
  !
  boysfun_up = boysfun_zero_order(x)
  do i = 0, n-1
    boysfun_up = ((2.0_wp*i+1.0_wp)*boysfun_up - exp(-x)) / (2.0_wp * x)
  end do
end function boysfun_up


!******************************************************************************** 
!
!       Interface for Boys functions
!
!******************************************************************************** 
real(wp) function boysfun(x, n)
  real(wp), intent(in) :: x
  integer, intent(in)  :: n
  !
  if (x == 0.0_wp) then
    boysfun = boysfun_zero(n)
  else if (n == 0) then
    boysfun = boysfun_zero_order(x)
  else if (x <= 5.0_wp) then
    boysfun = boysfun_small(x, n, kmax=25)
  else if (x >= 35.0_wp) then
    boysfun = boysfun_large(x, n)
  else
    boysfun = boysfun_an(x, n)
    if (boysfun==0.0_wp .and. x<=10.0_wp) then
      boysfun = boysfun_small(x, n)
    else if (boysfun/=boysfun .and. x<=10.0_wp) then
      boysfun = boysfun_small(x, n)
    else if (boysfun==0.0_wp .and. x>10.0_wp) then
      boysfun = boysfun_up(x, n)
    else if (boysfun/=boysfun .and. x>10.0_wp) then
      boysfun = boysfun_up(x, n)
    end if
  end if
end function boysfun


!******************************************************************************** 
!
!       Boys functions table F_n(x),   n = min_ .... max_
!
!       Use recursion to relate different order Boys functions with same argument
!
!******************************************************************************** 
subroutine boysfun_table(max_, x, boys)
  integer, intent(in)                :: max_
  real(wp), intent(in)               :: x
  real(wp), allocatable, intent(out) :: boys(:)
  !local variables
  integer                            :: n
  real(wp), parameter                :: max_x = 20.0_wp
  !
  if (.not. allocated(boys)) allocate(boys(0:max_))
  if (x < max_x) then
    !Downward recursion
    boys(max_) = boysfun(x, max_)
    do n = max_-1, 0, -1
      boys(n) = (2.0_wp*x*boys(n+1) + exp(-x)) / (2.0_wp*n + 1.0_wp)
    end do
  else
    !Upward recursion
    boys(0) = boysfun_zero_order(x)
    do n = 0, max_-1
      boys(n+1) = ((2.0_wp*n + 1.0_wp)*boys(n) - exp(-x)) / (2.0_wp * x)
    end do
  end if
end subroutine boysfun_table

end module boys_functions
