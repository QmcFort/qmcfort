! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module standalone

#include "preproc.inc"
!$use omp_lib
use constants
use mpi
use file_handle, only: FileHandle

implicit none

interface heaviside
  module procedure heaviside_sp, heaviside_dp
end interface heaviside

interface stepfun
  module procedure stepfun_sp, stepfun_dp
end interface stepfun

interface imag
  module procedure imag_s, imag_d, imag_rs, imag_rd
  !module procedure imag_q, imag_rq
end interface imag

interface conjug
  module procedure conjg_r_s, conjg_r_d, conjg_s, conjg_d
  !module procedure conjg_q
end interface conjug

interface phase
  module procedure phase_r, phase_c
end interface phase

interface exp_
  module procedure exp_r, exp_c
end interface exp_

interface exp_euler
  module procedure exp_euler_r, exp_euler_c
end interface exp_euler

interface exp_cn
  module procedure exp_cn_r, exp_cn_c
end interface exp_cn

interface exp_taylor
  module procedure exp_taylor_r, exp_taylor_c
end interface exp_taylor

interface reverse
  module procedure reverse_1d_i, reverse_1d_r, reverse_2d_r
end interface reverse

interface eye
  module procedure eye_i, eye_r, eye_c
end interface eye

interface maximum
  procedure maximum_sp, maximum_dp, maximum_int
end interface maximum

interface trace
  module procedure trace_r_sp, trace_r_dp, trace_c_sp, trace_c_dp, trace_int
end interface trace

interface trace2
  module procedure trace2_r_sp, trace2_r_dp, trace2_c_sp, trace2_c_dp
  module procedure trace2_rc_sp, trace2_rc_dp, trace2_cr_sp, trace2_cr_dp, trace2_int 
end interface trace2

interface blocktranspose
  module procedure blocktranspose_r_sp, blocktranspose_r_dp
  module procedure blocktranspose_c_sp, blocktranspose_c_dp
end interface blocktranspose

interface extract_matrix_cols
  module procedure extract_matrix_cols_sp, extract_matrix_cols_dp
end interface extract_matrix_cols

interface dump_array
  module procedure dump_array_1d_int, dump_array_1d_r_sp, dump_array_1d_r_dp, dump_array_1d_c_sp, dump_array_1d_c_dp
  module procedure dump_array_2d_int, dump_array_2d_r_sp, dump_array_2d_r_dp, dump_array_2d_c_sp, dump_array_2d_c_dp
end interface dump_array

interface trans_symm
  module procedure trans_symm_i, trans_symm_r_sp, trans_symm_r_dp, trans_symm_c_sp, trans_symm_c_dp
  module procedure trans_symm_int_r_sp, trans_symm_int_r_dp
end interface trans_symm

interface report_energy_value
  module procedure report_real_energy, report_cmplx_energy, report_cmplx_std_energy
end interface report_energy_value

interface near_zero
  module procedure near_zero_r_sp, near_zero_r_dp, near_zero_c_sp, near_zero_c_dp
end interface near_zero

interface concatenate
  module procedure concatenate_i
end interface 

contains

!******************************************************************************** 
!
!       Get appropriate size for monitoring iterative procedure
!
!******************************************************************************** 
integer function get_print_size(size_, nprints)
  integer, intent(in) :: size_
  integer, intent(in) :: nprints

  get_print_size = max(1, size_/nprints)

  if (get_print_size <= 10) then
    get_print_size = 10
  else if (get_print_size <= 40) then
    get_print_size = 25
  else if (get_print_size <= 80) then
    get_print_size = 50
  else if (get_print_size <= 200) then
    get_print_size = 100
  else if (get_print_size <= 400) then
    get_print_size = 250
  else if (get_print_size <= 800) then
    get_print_size = 500
  else if (get_print_size <= 2000) then
    get_print_size = 1000
  else if (get_print_size <= 4000) then
    get_print_size = 2500
  else if (get_print_size <= 8000) then
    get_print_size = 5000
  else 
    get_print_size = 10000
  end if
end function get_print_size


!******************************************************************************** 
!       
!       factorial function - implementation by do_loop
!       
!******************************************************************************** 
elemental integer function factorial(n) 
  integer, intent(in) :: n
  integer :: i
  !
  factorial = 1
  !
  do i = 2, n
    factorial = factorial * i
  end do
end function factorial


!******************************************************************************** 
!       
!       double factorial function - implementation by do_loop
!       
!******************************************************************************** 
elemental integer function double_factorial(n)
  integer, intent(in) :: n
  !local variables
  integer             :: i
  !
  double_factorial = 1
  !
  do i = n, 1, -2
    double_factorial = double_factorial * i
  end do
end function double_factorial


!******************************************************************************** 
!       
!       factorial function - implementation by recursion
!       
!******************************************************************************** 
recursive function factorial_rec(n)  result(fact)
  integer, intent(in) :: n
  integer             :: fact
  !
  if (n == 0) then
    fact = 1
  else
    fact = n * factorial_rec(n-1)
  end if
end function factorial_rec


!******************************************************************************** 
!       
!       double factorial function - implementation by recursion
!       
!******************************************************************************** 
recursive function double_factorial_rec(n) result(dfact)
  integer, intent(in) :: n
  integer             :: dfact
  !
  if (n == -1) then
    dfact = 1
  else if (n == 0) then
    dfact = 1
  else if (n == 1) then
    dfact = 1
  else
    dfact = n * double_factorial_rec(n-2)
  end if    
end function double_factorial_rec


!******************************************************************************** 
!       
! Binom coefficient (n, k) = n! / (k! (n-k)!)
!       
!******************************************************************************** 
function binom(n, k) result(bin_coeff)
  integer, intent(in) :: n
  integer, intent(in) :: k
  integer             :: bin_coeff
  !local variables
  integer             :: i

  bin_coeff = 1
  do i = k+1, n
    bin_coeff = bin_coeff * i
  end do

  bin_coeff = bin_coeff / factorial(n-k)
end function binom


!******************************************************************************** 
!       Heaviside function implementation
!                           |   1       if      x >= 0
!               H(x) =     {
!                           |   0       if      x <  0
!******************************************************************************** 
elemental real(sp) function heaviside_sp(x)
  real(sp), intent(in) :: x
  !
  if (x >= 0.0_sp) then
    heaviside_sp = 1.0_sp
  else
    heaviside_sp = 0.0_sp
  end if
end function heaviside_sp

elemental real(dp) function heaviside_dp(x)
  real(dp), intent(in) :: x
  !
  if (x >= 0.0_dp) then
    heaviside_dp = 1.0_dp
  else
    heaviside_dp = 0.0_dp
  end if
end function heaviside_dp


!******************************************************************************** 
!       Step function implementation
!                           |   1       if      x >= 0
!               s(x) =     {
!                           |  -1       if      x <  0
!******************************************************************************** 
elemental real(sp) function stepfun_sp(x)
  real(sp), intent(in) :: x
  !
  if (x >= 0.0_sp) then
    stepfun_sp = 1.0_sp
  else
    stepfun_sp = -1.0_sp
  end if
end function stepfun_sp

elemental real(dp) function stepfun_dp(x)
  real(dp), intent(in) :: x
  !
  if (x >= 0.0_dp) then
    stepfun_dp = 1.0_dp
  else
    stepfun_dp = -1.0_dp
  end if
end function stepfun_dp


!******************************************************************************** 
!
!       Generic subroutine to determine imaginary part of the complex numbers
!
!******************************************************************************** 
elemental real(sp) function imag_s(z)
  complex(sp), intent(in) :: z
  !
  imag_s = aimag(z)
end function imag_s

elemental real(dp) function imag_d(z)
  complex(dp), intent(in) :: z
  !
  imag_d = aimag(z)
end function imag_d

!debug: qprec not accepted by nvidia compilers
!!elemental real(qprec) function imag_q(z)
!!  complex(qprec), intent(in) :: z
!!  !
!!  imag_q = aimag(z)
!!end function imag_q

elemental real(sp) function imag_rs(x)
  real(sp),intent(in) :: x
  !
  imag_rs = 0.0_sp
end function imag_rs

elemental real(dp) function imag_rd(x)
  real(dp), intent(in) :: x
  !
  imag_rd = 0.0_wp
end function imag_rd

!!elemental real(qprec) function imag_rq(x)
!!  real(qprec), intent(in) :: x
!!  !
!!  imag_rq = x
!!end function imag_rq

!******************************************************************************** 
!       Function which determines conjugate of the given complex number
!********************************************************************************
elemental real(sp) function conjg_r_s(z)
  real(sp),intent(in):: z

  conjg_r_s = z
end function conjg_r_s

elemental real(dp) function conjg_r_d(z)
  real(dp),intent(in):: z

  conjg_r_d = z
end function conjg_r_d
 
elemental complex(sp) function conjg_s(z)
  complex(sp),intent(in):: z
  
  conjg_s = conjg(z)
end function conjg_s

elemental complex(dp) function conjg_d(z)
  complex(dp),intent(in):: z
  
  conjg_d = conjg(z)
end function conjg_d

!!elemental complex(qprec) function conjg_q(z)
!!  complex(qprec),intent(in):: z
!!  
!!  conjg_q = conjg(z)
!!end function conjg_q


!******************************************************************************** 
!
! Obtain 1d index of the matrix element with index a,b 
!
! Matrix is of dimension n_a x n_b
! 
!******************************************************************************** 
function get_1d_indx(a, b, n_a) result(indx)
  integer, intent(in) :: a
  integer, intent(in) :: b
  integer, intent(in) :: n_a
  integer             :: indx

  indx = (b-1)*n_a + a
end function get_1d_indx


!******************************************************************************** 
!
! Obtain 2d matrix index from the composite index indx
!
! Matrix is of dimension n_a x n_b
! 
!******************************************************************************** 
subroutine get_2d_indx(indx, n_a, a, b)
  integer, intent(in)  :: indx
  integer, intent(in)  :: n_a
  integer, intent(out) :: a
  integer, intent(out) :: b

  a = 1 + mod(indx-1, n_a)
  b = 1 + (indx-1)/n_a
end subroutine get_2d_indx


!******************************************************************************** 
!
! Transform index of the supermatrix into indices of four dim array
!    Default mode: 
!        opt=1 :      [i,j] --> a
!        opt=2 :      a     --> [i,j]
!
!    or use following for n x m matrix and index idx we can calculate
!
!    idx --> i,j   as follows:
!         
!    i = 1 + mod(idx-1, n)
!    j = 1 + (idx-1)/n
!
!******************************************************************************** 
subroutine get_index(n, m, i, j, a, opt)
  integer, intent(in) :: n, m
  integer             :: i, j
  integer             :: a
  integer, optional   :: opt
  !local variables
  integer             :: opt_
  
  if (present(opt)) then
    opt_ = opt
  else
    opt_ = 1
  end if
  
  if (opt_ == 1) then                  ! [i,j] --> a
    a = n* (j-1) + i
  else if (opt_ == 2) then             ! a --> [i,j]
    if (modul(a,n)) then
      i = n
      j = a / n
    else
      j = a / n + 1
      i = a - m * (j-1)
    end if
  end if
end subroutine get_index


!******************************************************************************** 
!
! Calculate the phase of the complex number z
!
!    phase = atan2(z_r, z_im)
!
!******************************************************************************** 
elemental real(wp) function phase_c(z)
  complex(wp), intent(in) :: z
  
  !phase_c = imag(log(z))
  phase_c = atan2(aimag(z), real(z,wp))
end function phase_c

!******************************************************************************** 
!
! Calculate the phase of the real number x
!
! If the real number is input, returns 0 or pi      
!
!******************************************************************************** 
elemental real(wp) function phase_r(x)
  real(wp), intent(in) :: x
  
  phase_r = pi * heaviside(-x)
end function phase_r


!******************************************************************************** 
!
!       Determine number of blocks given total number of elements and blocksize
!
!******************************************************************************** 
function get_nblocks(size_, bsize) result(nblocks)
  integer, intent(in) :: size_
  integer, intent(in) :: bsize
  integer             :: nblocks

  nblocks = ceiling(real(size_, wp) / bsize)
end function get_nblocks


!******************************************************************************** 
!
! Fo integer array of length gives starting positions of each segment
!
! See example below
!    ----------------------
!    |    9    |   5  | 3 |       [9, 5, 3] --> [0, 9, 14]
!    ----------------------
!   0|        9|    14|
!
! Usually used to produce nell array from nel 
!
!******************************************************************************** 
function get_start_pos(x) result(y)
  integer, intent(in)  :: x(:)
  integer, allocatable :: y(:)
  !local variables     
  integer              :: i, sum_

  allocate(y, mold=x)
  y(1) = 0
  sum_ = 0

  do i = 2, size(y)
    sum_ = sum_ + x(i-1)
    y(i) = sum_
  end do
end function get_start_pos


!******************************************************************************** 
!
!       Exponential interface 
!
!******************************************************************************** 
elemental function exp_r(x, kmax, prop) result(y)
  real(wp), intent(in)                   :: x
  integer, optional, intent(in)          :: kmax
  character(len=*), optional, intent(in) :: prop
  real(wp)                               :: y
  !local variables
  character(len=charlen)                 :: prop_
  !
  if (present(prop)) then
    prop_ = prop
  else
    prop_ = "exact"
  end if
  !
  select case (trim(prop_))
    case ("s2", "split2", "split") 
      y = exp_taylor(x, kmax)
    case ("crank_nicolson", "cn")
      y = exp_cn(x)
    case ("euler", "e")
      y = exp_euler(x)
    case ("taylor", "t")
      y = exp_taylor(x, kmax)
    case ("exact")
      y = exp(x)
    case default
      y = exp(x)
  end select
end function exp_r

elemental function exp_c(x, kmax, prop) result(y)
  complex(wp), intent(in)                :: x
  integer, optional, intent(in)          :: kmax
  character(len=*), optional, intent(in) :: prop
  complex(wp)                            :: y
  !local variables
  character(len=charlen)                 :: prop_
  !
  if (present(prop)) then
    prop_ = prop
  else
    prop_ = "exact"
  end if
  !
  select case (trim(prop_))
    case ("s2", "split2", "split") 
      y = exp_taylor(x, kmax)
    case ("crank_nicolson", "cn")
      y = exp_cn(x)
    case ("euler", "e")
      y = exp_euler(x)
    case ("taylor", "t")
      y = exp_taylor(x, kmax)
    case ("exact")
      y = exp(x)
    case default
      y = exp(x)
  end select
end function exp_c


!******************************************************************************** 
!
!       Function to approximate exponential by Euler method
!
!******************************************************************************** 
elemental function exp_euler_r(x) result(y)
  real(wp), intent(in) :: x
  real(wp)             :: y        
  !
  y = oner + x
end function exp_euler_r

elemental function exp_euler_c(x) result(y)
  complex(wp), intent(in) :: x
  complex(wp)             :: y        
  !
  y = onec + real(x, wp)
  y = y * (cos(imag(x)) + im*sin(imag(x)))
end function exp_euler_c


!******************************************************************************** 
!
!       Function to approximate exponential by Crank-Nicolson method
!
!                       y = (1 + x/2)/(1 - x/2)
!
!******************************************************************************** 
elemental function exp_cn_r(x) result(y)
  real(wp), intent(in) :: x
  real(wp)             :: y
  !
  y = (oner + x/2.0_wp) / (oner - x/2.0_wp)
end function exp_cn_r

elemental function exp_cn_c(x) result(y)
  complex(wp), intent(in) :: x
  complex(wp)             :: y
  !
  y = (onec + real(x, wp)/2.0_wp) / (onec - real(x, dp)/2.0_wp)
  y = y * (cos(imag(x)) + im*sin(imag(x)))
end function exp_cn_c


!******************************************************************************** 
!
!       Function to approximate exponential by Taylor series
!
!               exp(x) = 1 + x + x^2/2 + x^3/6 + x^4/24 + ...
!
!******************************************************************************** 
elemental function exp_taylor_r(x, kmax) result(y)
  real(wp), intent(in)          :: x
  integer, intent(in), optional :: kmax
  real(wp)                      :: y        
  !local variables
  real(wp)                      :: dx
  integer                       :: k, kmax_
  !
  if (present(kmax)) then
    kmax_ = kmax
  else
    kmax_ = 6
  end if 
  !
  y = 1.0_wp
  dx = 1.0_wp
  !
  do k = 1, kmax_
    dx = dx * x / real(k, wp)
    y = y + dx
  end do
end function exp_taylor_r

elemental function exp_taylor_c(x, kmax) result(y)
  complex(wp), intent(in)       :: x
  integer, intent(in), optional :: kmax
  complex(wp)                   :: y
  !local variables
  complex(wp)                   :: dx
  integer                       :: k, kmax_
  !
  if (present(kmax)) then
    kmax_ = kmax
  else
    kmax_ = 6 
  end if 
  !
  y = onec
  dx = oner
  !
  do k = 1, kmax_
    dx = dx * x / real(k, wp)
    y = y + dx
  end do
  !
  y = y * (cos(imag(x)) + im*sin(imag(x)))
end function exp_taylor_c


!******************************************************************************** 
!
! Cap values outside the region  mean +- scale*sigma
!
!******************************************************************************** 
function cap_value(x, mean, sigma, scale) result(xx)
  real(wp), intent(in) :: x
  real(wp), intent(in) :: mean
  real(wp), intent(in) :: sigma
  real(wp), intent(in) :: scale
  real(wp)             :: xx

  if (x > mean + scale*sigma) then
    xx = mean + scale*sigma
  else if (x < mean - scale*sigma) then
    xx = mean - scale*sigma
  else 
    xx = x
  end if
end function cap_value


!******************************************************************************** 
!
! Tensor product of two 1d vectors
!
!******************************************************************************** 
function tensor_product(a, b) result(c)
  real(wp), intent(in)  :: a(:)
  real(wp), intent(in)  :: b(:)
  real(wp), allocatable :: c(:)
  !local variables
  integer               :: i, i1, i2

  allocate(c(size(a)*size(b)))

  do i = 1, size(b)
     i1 = (i-1)*size(a) + 1
     i2 = i * size(a)
     c(i1:i2) = a * b(i)
  end do
end function tensor_product


!******************************************************************************** 
!
!       Subroutine that returns identity matrix inplace
!
!               Called by interface eye
!
!******************************************************************************** 
subroutine eye_i(matrix)
  integer, intent(inout) :: matrix(:,:)
  !local
  integer                :: i, n, m
  
  n = size(matrix, 1)
  m = size(matrix, 2)

  matrix = 0

  do i = 1, min(n,m)
    matrix(i,i) = 1
  end do
end subroutine eye_i

subroutine eye_r(matrix)
  real(wp), intent(inout) :: matrix(:,:)
  !local
  integer                 :: i, n, m

  n = size(matrix, 1)
  m = size(matrix, 2)

  matrix = 0.0_wp

  do i = 1, min(n,m)
    matrix(i,i) = 1.0_wp
  end do
end subroutine eye_r

subroutine eye_c(matrix)
  complex(wp), intent(inout) :: matrix(:,:)
  !local
  integer                    :: i, n, m
  
  n = size(matrix, 1)
  m = size(matrix, 2)

  matrix = (0.0_wp, 0.0_wp)

  do i = 1, min(n,m)
    matrix(i,i) = (1.0_wp, 0.0_wp)
  end do
end subroutine eye_c


!******************************************************************************** 
!
!       Function that returns identity matrix inplace
!
!               Called by interface eye
!
!******************************************************************************** 
function eyef_i(n) result(eye)
  integer, intent(in) :: n
  integer             :: eye(n, n)
  !local variables
  integer             :: i 
  !
  eye = 0
  do i = 1, n
    eye(i,i) = 1
  end do
end function eyef_i

function eyef_r(n) result(eye)
  integer, intent(in) :: n
  real(wp)            :: eye(n, n)
  !local variables
  integer             :: i 
  !
  eye = 0.0_wp
  do i = 1, n
    eye(i,i) = 1.0_wp
  end do
end function eyef_r

function eyef_c(n) result(eye)
  integer, intent(in) :: n
  complex(wp)         :: eye(n, n)
  !local variables
  integer             :: i 
  !
  eye = (0.0_wp, 0.0_wp)
  do i = 1, n
    eye(i,i) = (1.0_wp, 0.0_wp)
  end do
end function eyef_c


!******************************************************************************** 
!
!       Reverse Array
!
!******************************************************************************** 
subroutine reverse_1d_i(x)
  integer, intent(inout) :: x(:)
  !local variables
  integer                :: i, n, mid
  integer                :: swap
  !
  n = size(x)
  mid = n / 2
  !
  do i = 1, mid
    swap = x(i)
    x(i) = x(n+1-i)
    x(n+1-i) = swap      
  end do
end subroutine reverse_1d_i

subroutine reverse_1d_r(x)
  real(wp), intent(inout) :: x(:)
  !local variables
  integer                 :: i, n, mid
  real(wp)                :: swap
  !
  n = size(x)
  mid = n / 2
  !
  do i = 1, mid
    swap = x(i)
    x(i) = x(n+1-i)
    x(n+1-i) = swap      
  end do
end subroutine reverse_1d_r

subroutine reverse_2d_r(x)
  real(wp), intent(inout) :: x(:,:)
  !local variables
  integer                 :: i, n, mid
  real(wp), allocatable   :: swap(:)
  !
  n = size(x, 2)
  mid = n / 2
  allocate(swap(n))
  !
  do i = 1, mid
    swap = x(:,i)
    x(:,i) = x(:,n+1-i)
    x(:,n+1-i) = swap
  end do
end subroutine reverse_2d_r


!******************************************************************************** 
!
!       Logical function that basically returns mod(a, b) == 0
!
!******************************************************************************** 
logical function modul(a, b)
  integer, intent(in) :: a, b
  !
  if (a == 0) modul = .false.
  !
  if (b <= 0) then
    modul = .false.
  else
    modul = mod(a,b) == 0
  end if  
end function modul


!******************************************************************************** 
!
!   Search for the maximum of the array
!
!******************************************************************************** 
function maximum_sp(N,vec) result(val)
  integer,intent(in):: N
  real(sp),dimension(N),intent(in):: vec
  real(sp):: val

  !local variables
  integer:: i
  
  val = vec(1)

  do i=2,N
    val = MAX(val,vec(i))    
  end do
end function maximum_sp

function maximum_dp(N,vec) result(val)
  integer,intent(in):: N
  real(dp),dimension(N),intent(in):: vec
  real(dp):: val

  !local variables
  integer:: i
  
  val = vec(1)

  do i=2,N
    val = MAX(val,vec(i))    
  end do
end function maximum_dp

function maximum_int(N,vec) result(val)
  integer,intent(in):: N
  integer,dimension(N),intent(in):: vec
  integer:: val

  !local variables
  integer:: i
  
  val = vec(1)

  do i=2,N
    val = MAX(val,vec(i))    
  end do
end function maximum_int


!******************************************************************************** 
!
!       trace - function that calculates Trace of the square matrix
!               inPUT:
!                       MAT - square matrix
!                       res - index that restricts summation on the diagonal
!               outPUT:
!                       TR  - Trace of the matrix
!
!
!       trace2 - function that caclulates Trace of the matrix product MAT1 * MAT2
!                Dimensions of the matrices MAT1 and MAT2 must match for the
!                matrix multiplicatio, MAT1 * MAT2 must be square matrix
!               inPUT:
!                       MAT1 - first matrix with dimensions N x M
!                       MAT2 - second matrix with dimension M X N
!                       res - index that restricts summation on diagonal of MAT1 * MAT2
!               outPUT: 
!                       TR - Trace of the matrix product
!
!******************************************************************************** 
function trace_r_dp(MAT,res) result(TR)
  real(dp),intent(in):: MAT(:,:)
  real(dp):: TR
  integer,optional:: res

  !local variables
  integer:: i
  integer:: N , M , myres, n1, m1

  N = size(MAT(:,1))
  M = size(MAT(1,:))

  if(N/=M) then
    write(*,*) "Problem in function trace()"
    write(*,*) N, N1
    write(*,*) M, M1
    write(*,*) "Inconsistent matrix dimensions - Trace is ill defined"
    return
  end if

  if(present(res)) then
    myres = res
  else
    myres = N
  end if

  TR = 0.0_dp

  do i=1,myres
    TR = TR + MAT(i,i) 
  end do
end function trace_r_dp

function trace_c_dp(MAT,res) result(TR)
  complex(dp),intent(in):: MAT(:,:)
  complex(dp):: TR
  integer,optional:: res

  !local variables
  integer:: i
  integer:: N , M , myres, n1, m1

  N = size(MAT(:,1))
  M = size(MAT(1,:))

  if(N/=M) then
    write(*,*) "Problem in function trace()"
    write(*,*) N, N1
    write(*,*) M, M1
    write(*,*) "Inconsistent matrix dimensions - Trace is ill defined"
    return
  end if

  if(present(res)) then
    myres = res
  else
    myres = N
  end if

  TR = (0.0_dp,0.0_dp)

  do i=1,myres
    TR = TR + MAT(i,i) 
  end do
end function trace_c_dp

function trace_r_sp(MAT,res) result(TR)
  real(sp),intent(in):: MAT(:,:)
  real(sp):: TR
  integer,optional:: res

  !local variables
  integer:: i
  integer:: N , M , myres, m1, n1

  N = size(MAT(:,1))
  M = size(MAT(1,:))

  if(N/=M) then
    write(*,*) "Problem in function trace()"
    write(*,*) N, N1
    write(*,*) M, M1
    write(*,*) "Inconsistent matrix dimensions - Trace is ill defined"
    return
  end if

  if(present(res))then
    myres = res
  else
    myres = N 
  end if

  TR = 0.0_sp

  do i=1,myres
    TR = TR + MAT(i,i) 
  end do
end function trace_r_sp

function trace_c_sp(MAT,res) result(TR)
  complex(sp),intent(in):: MAT(:,:)
  complex(sp):: TR
  integer,optional:: res

  !local variables
  integer:: i
  integer:: N , M , myres, n1, m1

  N = size(MAT(:,1))
  M = size(MAT(1,:))

  if(N/=M) then
    write(*,*) "Problem in function trace()"
    write(*,*) N, N1
    write(*,*) M, M1
    write(*,*) "Inconsistent matrix dimensions - Trace is ill defined"
    return
  end if

  if(present(res))then
    myres = res
  else
    myres = N 
  end if

  TR = (0.0_sp,0.0_sp)

  do i=1,myres
    TR = TR + MAT(i,i) 
  end do
end function trace_c_sp

function trace_int(MAT,res) result(TR)
  integer,intent(in):: MAT(:,:)
  integer:: TR
  integer,optional:: res

  !local variables
  integer:: i
  integer:: N , M , myres , n1, m1

  N = size(MAT(:,1))
  M = size(MAT(1,:))

  if(N/=M) then
    write(*,*) "Problem in function trace()"
    write(*,*) N, N1
    write(*,*) M, M1
    write(*,*) "Inconsistent matrix dimensions - Trace is ill defined"
    return
  end if

  if(present(res)) then
    myres = res
  else
    myres = N
  end if

  TR = 0

  do i=1,myres
    TR = TR + MAT(i,i) 
  end do
end function trace_int

function trace2_r_sp(a, b) result(tr)
  real(sp), intent(in) :: a(:,:)
  real(sp), intent(in) :: b(:,:)
  real(sp)             :: tr
  !local variables
  integer                 :: i, i1, i2, n1, m1, n2, m2
  
  n1 = size(a, 1)
  m1 = size(a, 2)
  n2 = size(b, 1)
  m2 = size(b, 2)

  if (mod(m2, m1)==0 .and. n1==n2) then
    if (n1 == m1) then
      tr = 0.0_sp
      do i = 1, m2/m1
        i1 = (i-1)*m1 + 1
        i2 = i * m1
        tr = tr + sum(a*transpose(b(:,i1:i2)))
      end do
    else
      tr = 0.0_sp
      do i = 1, m2/m1
        i1 = (i-1)*m1 + 1
        i2 = i * m1
        tr = tr + sum(a*b(:,i1:i2))
      end do
    end if
  else
    write(*,*) "Problem in function trace2 - inconsistent dimensions of factor matrices trace(A*B)"
    write(*,*) "First array dimensions  = ", n1, m1
    write(*,*) "Second array dimensions = ", n2, m2
    return
  end if
end function trace2_r_sp

function trace2_r_dp(a, b) result(tr)
  real(dp), intent(in) :: a(:,:)
  real(dp), intent(in) :: b(:,:)
  real(dp)             :: tr
  !local variables
  integer                 :: i, i1, i2, n1, m1, n2, m2
  
  n1 = size(a, 1)
  m1 = size(a, 2)
  n2 = size(b, 1)
  m2 = size(b, 2)

  if (mod(m2, m1)==0 .and. n1==n2) then
    if (n1 == m1) then
      tr = 0.0_dp
      do i = 1, m2/m1
        i1 = (i-1)*m1 + 1
        i2 = i * m1
        tr = tr + sum(a*transpose(b(:,i1:i2)))
      end do
    else 
      tr = 0.0_dp
      do i = 1, m2/m1
        i1 = (i-1)*m1 + 1
        i2 = i * m1
        tr = tr + sum(a*b(:,i1:i2))
      end do
    end if
  else
    write(*,*) "Problem in function trace2 - inconsistent dimensions of factor matrices trace(A*B)"
    write(*,*) "First array dimensions  = ", n1, m1
    write(*,*) "Second array dimensions = ", n2, m2
    return
  end if
end function trace2_r_dp

function trace2_c_sp(a, b) result(tr)
  complex(sp), intent(in) :: a(:,:)
  complex(sp), intent(in) :: b(:,:)
  complex(sp)             :: tr
  !local variables
  integer                    :: i, i1, i2, n1, m1, n2, m2
  
  n1 = size(a, 1)
  m1 = size(a, 2)
  n2 = size(b, 1)
  m2 = size(b, 2)

  if (mod(m2, m1)==0 .and. n1==n2) then
    if (n1 == m1) then
      tr = 0.0_sp
      do i = 1, m2/m1
        i1 = (i-1)*m1 + 1
        i2 = i * m1
        tr = tr + sum(a*transpose(b(:,i1:i2)))
      end do
    else 
      tr = 0.0_sp
      do i = 1, m2/m1
        i1 = (i-1)*m1 + 1
        i2 = i * m1
        tr = tr + sum(a*b(:,i1:i2))
      end do
    end if
  else
    write(*,*) "Problem in function trace2 - inconsistent dimensions of factor matrices trace(A*B)"
    write(*,*) "First array dimensions  = ", n1, m1
    write(*,*) "Second array dimensions = ", n2, m2
    return
  end if
end function trace2_c_sp

function trace2_c_dp(a, b) result(tr)
  complex(dp), intent(in) :: a(:,:)
  complex(dp), intent(in) :: b(:,:)
  complex(dp)             :: tr
  !local variables
  integer                    :: i, i1, i2, n1, m1, n2, m2
  
  n1 = size(a, 1)
  m1 = size(a, 2)
  n2 = size(b, 1)
  m2 = size(b, 2)

  if (mod(m2, m1)==0 .and. n1==n2) then
    if (n1 == m1) then
      tr = 0.0_dp
      do i = 1, m2/m1
        i1 = (i-1)*m1 + 1
        i2 = i * m1
        tr = tr + sum(a*transpose(b(:,i1:i2)))
      end do
    else 
      tr = 0.0_dp
      do i = 1, m2/m1
        i1 = (i-1)*m1 + 1
        i2 = i * m1
        tr = tr + sum(a*b(:,i1:i2))
      end do
    end if
  else
    write(*,*) "Problem in function trace2 - inconsistent dimensions of factor matrices trace(A*B)"
    write(*,*) "First array dimensions  = ", n1, m1
    write(*,*) "Second array dimensions = ", n2, m2
    return
  end if
end function trace2_c_dp

function trace2_rc_sp(a, b) result(tr)
  real(sp), intent(in)    :: a(:,:)
  complex(sp), intent(in) :: b(:,:)
  complex(sp)             :: tr
  !local variables
  integer                    :: i, i1, i2, n1, m1, n2, m2
  
  n1 = size(a, 1)
  m1 = size(a, 2)
  n2 = size(b, 1)
  m2 = size(b, 2)

  if (mod(m2, m1)==0 .and. n1==n2) then
    if (n1 == m1) then
      tr = 0.0_sp
      do i = 1, m2/m1
        i1 = (i-1)*m1 + 1
        i2 = i * m1
        tr = tr + sum(a*transpose(b(:,i1:i2)))
      end do
    else 
      tr = 0.0_sp
      do i = 1, m2/m1
        i1 = (i-1)*m1 + 1
        i2 = i * m1
        tr = tr + sum(a*b(:,i1:i2))
      end do
    end if
  else
    write(*,*) "Problem in function trace2 - inconsistent dimensions of factor matrices trace(A*B)"
    write(*,*) "First array dimensions  = ", n1, m1
    write(*,*) "Second array dimensions = ", n2, m2
    return
  end if
end function trace2_rc_sp

function trace2_rc_dp(a, b) result(tr)
  real(dp), intent(in)    :: a(:,:)
  complex(dp), intent(in) :: b(:,:)
  complex(dp)             :: tr
  !local variables
  integer                    :: i, i1, i2, n1, m1, n2, m2
  
  n1 = size(a, 1)
  m1 = size(a, 2)
  n2 = size(b, 1)
  m2 = size(b, 2)

  if (mod(m2, m1)==0 .and. n1==n2) then
    if (n1 == m1) then
      tr = 0.0_dp
      do i = 1, m2/m1
        i1 = (i-1)*m1 + 1
        i2 = i * m1
        tr = tr + sum(a*transpose(b(:,i1:i2)))
      end do
    else
      tr = 0.0_dp
      do i = 1, m2/m1
        i1 = (i-1)*m1 + 1
        i2 = i * m1
        tr = tr + sum(a*b(:,i1:i2))
      end do
    end if
  else
    write(*,*) "Problem in function trace2 - inconsistent dimensions of factor matrices trace(A*B)"
    write(*,*) "First array dimensions  = ", n1, m1
    write(*,*) "Second array dimensions = ", n2, m2
    return
  end if
end function trace2_rc_dp

function trace2_cr_sp(a, b) result(tr)
  complex(sp), intent(in) :: a(:,:)
  real(sp), intent(in)    :: b(:,:)
  complex(sp)             :: tr
  !local variables
  integer                    :: i, i1, i2, n1, m1, n2, m2
  
  n1 = size(a, 1)
  m1 = size(a, 2)
  n2 = size(b, 1)
  m2 = size(b, 2)

  if (mod(m2, m1)==0 .and. n1==n2) then
    if (n1 == m1) then
      tr = 0.0_sp
      do i = 1, m2/m1
        i1 = (i-1)*m1 + 1
        i2 = i * m1
        tr = tr + sum(a*transpose(b(:,i1:i2)))
      end do
    else 
      tr = 0.0_sp
      do i = 1, m2/m1
        i1 = (i-1)*m1 + 1
        i2 = i * m1
        tr = tr + sum(a*b(:,i1:i2))
      end do
    end if
  else
    write(*,*) "Problem in function trace2 - inconsistent dimensions of factor matrices trace(A*B)"
    write(*,*) "First array dimensions  = ", n1, m1
    write(*,*) "Second array dimensions = ", n2, m2
    return
  end if
end function trace2_cr_sp

function trace2_cr_dp(a, b) result(tr)
  complex(dp), intent(in) :: a(:,:)
  real(dp), intent(in)    :: b(:,:)
  complex(dp)             :: tr
  !local variables
  integer                    :: i, i1, i2, n1, m1, n2, m2
  
  n1 = size(a, 1)
  m1 = size(a, 2)
  n2 = size(b, 1)
  m2 = size(b, 2)

  if (mod(m2, m1)==0 .and. n1==n2) then
    if (n1 == m1) then
      tr = 0.0_dp
      do i = 1, m2/m1
        i1 = (i-1)*m1 + 1
        i2 = i * m1
        tr = tr + sum(a*transpose(b(:,i1:i2)))
      end do
    else
      tr = 0.0_dp
      do i = 1, m2/m1
        i1 = (i-1)*m1 + 1
        i2 = i * m1
        tr = tr + sum(a*b(:,i1:i2))
      end do
    end if
  else
    write(*,*) "Problem in function trace2 - inconsistent dimensions of factor matrices trace(A*B)"
    write(*,*) "First array dimensions  = ", n1, m1
    write(*,*) "Second array dimensions = ", n2, m2
    return
  end if
end function trace2_cr_dp

function trace2_int(a, b) result(tr)
  integer, intent(in) :: a(:,:)
  integer, intent(in) :: b(:,:)
  integer             :: tr
  !local variables
  integer             :: i, i1, i2, n1, m1, n2, m2
  
  n1 = size(a, 1)
  m1 = size(a, 2)
  n2 = size(b, 1)
  m2 = size(b, 2)

  if (mod(m2, m1)==0 .and. n1==n2) then
    tr = 0
    do i = 1, m2/m1
      i1 = (i-1)*m1 + 1
      i2 = i * m1
      tr = tr + sum(a*b(:,i1:i2))
    end do
  else
    write(*,*) "Problem in function trace2 - inconsistent dimensions of factor matrices trace(A*B)"
    write(*,*) "First array dimensions  = ", n1, m1
    write(*,*) "Second array dimensions = ", n2, m2
    return
  end if
end function trace2_int


!******************************************************************************** 
!
!       trans_symm routine:
!               Applies transposition symmetry on the matrix
!               per default - copies lower part to upper
!                               [i|j] = [j|i]
!
!******************************************************************************** 
subroutine trans_symm_i(a, low)
  use string, only: lowercase
  integer, intent(inout)                 :: a(:,:)
  character(len=*), optional, intent(in) :: low
  !local variables
  integer                                :: n, m
  integer                                :: i, j
  character(len=:), allocatable          :: low_
  !
  n = size(a, 1)
  m = size(a, 2) 
  if (n /= m) write(*,"(a,i6,i6)") "Transposal symmetrie applied to rectangular matrix: n, m =", n, m
  !
  if (present(low)) then
    low_ = low
  else 
    low_ = "low"
  end if
  !
  if (lowercase(low_(1:1)) == "u") then
    do i = 1, n
      do j = 1, i-1
        a(i,j) = a(j,i)
      end do
    end do
  else
    do i = 1, n
      do j = 1, i-1
        a(j,i) = a(i,j)
      end do
    end do
  end if
end subroutine trans_symm_i

subroutine trans_symm_r_sp(a, low)
  use string, only: lowercase
  real(sp), intent(inout)             :: a(:,:)
  character(len=*), optional, intent(in) :: low
  !local variables
  integer                                :: n, m
  integer                                :: i, j
  character(len=:), allocatable          :: low_
  !
  n = size(a, 1)
  m = size(a, 2) 
  if (n /= m) write(*,"(a,i6,i6)") "Transposal symmetrie applied to rectangular matrix: n, m =", n, m
  !
  if (present(low)) then
    low_ = low
  else 
    low_ = "low"
  end if
  !
  if (lowercase(low_(1:1)) == "u") then
    do i = 1, n
      do j = 1, i-1
        a(i,j) = a(j,i)
      end do
    end do
  else
    do i = 1, n
      do j = 1, i-1
        a(j,i) = a(i,j)
      end do
    end do
  end if
end subroutine trans_symm_r_sp

subroutine trans_symm_r_dp(a, low)
  use string, only: lowercase
  real(dp), intent(inout)             :: a(:,:)
  character(len=*), optional, intent(in) :: low
  !local variables
  integer                                :: n, m
  integer                                :: i, j
  character(len=:), allocatable          :: low_
  !
  n = size(a, 1)
  m = size(a, 2) 
  if (n /= m) write(*,"(a,i6,i6)") "Transposal symmetrie applied to rectangular matrix: n, m =", n, m
  !
  if (present(low)) then
    low_ = low
  else 
    low_ = "low"
  end if
  !
  if (lowercase(low_(1:1)) == "u") then
    do i = 1, n
      do j = 1, i-1
        a(i,j) = a(j,i)
      end do
    end do
  else
    do i = 1, n
      do j = 1, i-1
        a(j,i) = a(i,j)
      end do
    end do
  end if
end subroutine trans_symm_r_dp

subroutine trans_symm_c_sp(a, low)
  use string, only: lowercase
  complex(sp), intent(inout)          :: a(:,:)
  character(len=*), optional, intent(in) :: low
  !local variables
  integer                                :: n, m
  integer                                :: i, j
  character(len=:), allocatable          :: low_
  !
  n = size(a, 1)
  m = size(a, 2) 
  if (n /= m) write(*,"(a,i6,i6)") "Transposal symmetrie applied to rectangular matrix: n, m =", n, m
  !
  if (present(low)) then
    low_ = low
  else 
    low_ = "low"
  end if
  !
  if (lowercase(low_(1:1)) == "u") then
    do i = 1, n
      do j = 1, i-1
        a(i,j) = conjg(a(j,i))
      end do
    end do
  else
    do i = 1, n
      do j = 1, i-1
        a(j,i) = a(i,j)
      end do
    end do
  end if
end subroutine trans_symm_c_sp

subroutine trans_symm_c_dp(a, low)
  use string, only: lowercase
  complex(dp), intent(inout)          :: a(:,:)
  character(len=*), optional, intent(in) :: low
  !local variables
  integer                                :: n, m
  integer                                :: i, j
  character(len=:), allocatable          :: low_
  !
  n = size(a, 1)
  m = size(a, 2) 
  if (n /= m) write(*,"(a,i6,i6)") "Transposal symmetrie applied to rectangular matrix: n, m =", n, m
  !
  if (present(low)) then
    low_ = low
  else 
    low_ = "low"
  end if
  !
  if (lowercase(low_(1:1)) == "u") then
    do i = 1, n
      do j = 1, i-1
        a(i,j) = a(j,i)
      end do
    end do
  else
    do i = 1, n
      do j = 1, i-1
        a(j,i) = conjg(a(i,j))
      end do
    end do
  end if
end subroutine trans_symm_c_dp


!******************************************************************************** 
!
!       Applies symmetry on Molecular integrals
!           [pq|rs] = [qp|rs] = [pq|sr] = [qp|sr]
!         = [rs|pq] = [sr|pq] = [rs|qp] = [sr|qp]            
!
!******************************************************************************** 
subroutine trans_symm_int_r_sp(a)
  real(sp), intent(inout) :: a(:,:,:,:)
  !local variables
  integer                    :: n, p, q, r, s
  !
  n = size(a, 1)
  !$omp parallel do schedule(dynamic) private(p,q,r,s) 
  do s = 1, n
    do r = 1, s
      do q = 1, s
        do p = 1, q
          a(q,p,r,s) = a(p,q,r,s)
          a(p,q,s,r) = a(p,q,r,s)
          a(q,p,s,r) = a(p,q,r,s)
          a(r,s,p,q) = a(p,q,r,s)
          a(s,r,p,q) = a(p,q,r,s)
          a(r,s,q,p) = a(p,q,r,s)
          a(s,r,q,p) = a(p,q,r,s)
        end do
      end do
    end do
  end do
  !$omp end parallel do
end subroutine trans_symm_int_r_sp

subroutine trans_symm_int_r_dp(a)
  real(dp), intent(inout) :: a(:,:,:,:)
  !local variables
  integer                    :: n, p, q, r, s
  !
  n = size(a, 1)
  !$omp parallel do schedule(dynamic) private(p,q,r,s)
  do s = 1, n
    do r = 1, s
      do q = 1, s
        do p = 1, q
          a(q,p,r,s) = a(p,q,r,s)
          a(p,q,s,r) = a(p,q,r,s)
          a(q,p,s,r) = a(p,q,r,s)
          a(r,s,p,q) = a(p,q,r,s)
          a(s,r,p,q) = a(p,q,r,s)
          a(r,s,q,p) = a(p,q,r,s)
          a(s,r,q,p) = a(p,q,r,s)
        end do
      end do
    end do
  end do
  !$omp end parallel do
end subroutine trans_symm_int_r_dp


!******************************************************************************** 
!
!       Transpose a block of matrices strored in 2d array
!
!       Transposition done in dimension with larger size
!
!******************************************************************************** 
function blocktranspose_r_sp(a) result(b)
  real(sp), intent(in)  :: a(:,:)
  real(sp), allocatable :: b(:,:)
  !local variables
  integer                  :: n, m, i, w1, w2
  !
  allocate(b, mold=a)
  n = size(a, 1)
  m = size(a, 2)
  !
  if (n >= m) then
    do i = 1, n/m
      w1 = (i-1) * m + 1
      w2 = i * m
      b(w1:w2, 1:m) = transpose(a(w1:w2, 1:m))
    end do
  else
    do i = 1, m/n
      w1 = (i-1) * n + 1
      w2 = i * n
      b(1:n, w1:w2) = transpose(a(1:n, w1:w2))
    end do
  end if
end function blocktranspose_r_sp

function blocktranspose_r_dp(a) result(b)
  real(dp), intent(in)  :: a(:,:)
  real(dp), allocatable :: b(:,:)
  !local variables
  integer                  :: n, m, i, w1, w2
  !
  allocate(b, mold=a)
  n = size(a, 1)
  m = size(a, 2)
  !
  if (n >= m) then
    do i = 1, n/m
      w1 = (i-1) * m + 1
      w2 = i * m
      b(w1:w2, 1:m) = transpose(a(w1:w2, 1:m))
    end do
  else
    do i = 1, m/n
      w1 = (i-1) * n + 1
      w2 = i * n
      b(1:n, w1:w2) = transpose(a(1:n, w1:w2))
    end do
  end if
end function blocktranspose_r_dp

function blocktranspose_c_sp(a) result(b)
  complex(sp), intent(in)  :: a(:,:)
  complex(sp), allocatable :: b(:,:)
  !local variables
  integer                     :: n, m, i, w1, w2
  !
  allocate(b, mold=a)
  n = size(a, 1)
  m = size(a, 2)
  !
  if (n >= m) then
    do i = 1, n/m
      w1 = (i-1) * m + 1
      w2 = i * m
      b(w1:w2, 1:m) = transpose(a(w1:w2, 1:m))
    end do
  else
    do i = 1, m/n
      w1 = (i-1) * n + 1
      w2 = i * n
      b(1:n, w1:w2) = transpose(a(1:n, w1:w2))
    end do
  end if
end function blocktranspose_c_sp

function blocktranspose_c_dp(a) result(b)
  complex(dp), intent(in)  :: a(:,:)
  complex(dp), allocatable :: b(:,:)
  !local variables
  integer                     :: n, m, i, w1, w2
  !
  allocate(b, mold=a)
  n = size(a, 1)
  m = size(a, 2)
  !
  if (n >= m) then
    do i = 1, n/m
      w1 = (i-1) * m + 1
      w2 = i * m
      b(w1:w2, 1:m) = transpose(a(w1:w2, 1:m))
    end do
  else
    do i = 1, m/n
      w1 = (i-1) * n + 1
      w2 = i * n
      b(1:n, w1:w2) = transpose(a(1:n, w1:w2))
    end do
  end if
end function blocktranspose_c_dp

!******************************************************************************** 
!
!       Extract eigenvalues and eigenvectos of the matrix according to the mask
!
!       Input:
!               u       - eigenvectors
!               e       - eigenvalues
!               mask    - 1d mask to extract desired eigenvalues
!       Output:
!               u_      - extracted eigenvectors for which mask is true    
!               e_      - extracted eigenvalues for which mask is true
!
!******************************************************************************** 
subroutine extract_matrix_cols_sp(u, w, mask, u_, w_)
  real(sp), intent(in)               :: u(:,:)
  real(sp), intent(in)               :: w(:)
  logical, intent(in)                   :: mask(:)
  real(sp), allocatable, intent(out) :: u_(:,:)
  real(sp), allocatable, intent(out) :: w_(:)
  !local variables
  logical, allocatable                  :: mask_2d(:,:)
  real(sp), allocatable              :: u_pack(:)
  !
  w_ = pack(w, mask)
  mask_2d = spread(mask, 1, size(u,1))
  u_pack = pack(u, mask_2d)
  u_ = reshape(u_pack, [size(u,1), size(w_)])
end subroutine extract_matrix_cols_sp

subroutine extract_matrix_cols_dp(u, w, mask, u_, w_)
  real(dp), intent(in)               :: u(:,:)
  real(dp), intent(in)               :: w(:)
  logical, intent(in)                   :: mask(:)
  real(dp), allocatable, intent(out) :: u_(:,:)
  real(dp), allocatable, intent(out) :: w_(:)
  !local variables
  logical, allocatable                  :: mask_2d(:,:)
  real(dp), allocatable              :: u_pack(:)
  !
  w_ = pack(w, mask)
  mask_2d = spread(mask, 1, size(u,1))
  u_pack = pack(u, mask_2d)
  u_ = reshape(u_pack, [size(u,1), size(w_)])
end subroutine extract_matrix_cols_dp


!******************************************************************************** 
!
! Concatenate two 1d arrays
!
!******************************************************************************** 
function concatenate_i(array1, array2) result(array)
  integer, intent(in)  :: array1(:)
  integer, intent(in)  :: array2(:)
  integer, allocatable :: array(:)
  !local 
  integer              :: n1, n2

  n1 = size(array1)
  n2 = size(array2)

  allocate(array(n1+n2))
  array(1:n1) = array1
  array(n1+1:n1+n2) = array2
end function concatenate_i


!******************************************************************************** 
!
! Write structured array to the screen
!
! very useful for debugging purposes
!
! nrows/ncols specifies how many rows/columns to print
!    default is 10 elements
!    if negative integer supplied, all elements printed
!
!******************************************************************************** 
subroutine dump_array_1d_int(A, nrows)
  integer, intent(in)           :: A(:)
  integer, optional, intent(in) :: nrows
  !local variables
  integer                       :: nrows_, i
  integer, parameter            :: maxlen = 10
  
  if (present(nrows)) then
    if (nrows <= 0) then
      nrows_ = size(A)
    else
      nrows_ = min(nrows, size(A))
    end if
  else
    nrows_ = min(maxlen, size(A))
  end if

  do i = 1, nrows_
    write(*,"(I8)") A(i)
  end do

  write(*,*)
end subroutine dump_array_1d_int

subroutine dump_array_1d_r_sp(A, nrows)
  real(sp), intent(in)          :: A(:)
  integer, optional, intent(in) :: nrows
  !local variables
  integer                       :: nrows_, i
  integer, parameter            :: maxlen = 10
  
  if (present(nrows)) then
    if (nrows <= 0) then
      nrows_ = size(A)
    else
      nrows_ = min(nrows, size(A))
    end if
  else
    nrows_ = min(maxlen, size(A))
  end if

  do i = 1, nrows_
    write(*,"(f10.4)") A(i)
  end do

  write(*,*)
end subroutine dump_array_1d_r_sp

subroutine dump_array_1d_r_dp(A, nrows)
  real(dp), intent(in)          :: A(:)
  integer, optional, intent(in) :: nrows
  !local variables
  integer                       :: nrows_, i
  integer, parameter            :: maxlen = 10
  
  if (present(nrows)) then
    if (nrows <= 0) then
      nrows_ = size(A)
    else
      nrows_ = min(nrows, size(A))
    end if
  else
    nrows_ = min(maxlen, size(A))
  end if

  do i = 1, nrows_
    write(*,"(f10.4)") A(i)
  end do

  write(*,*)
end subroutine dump_array_1d_r_dp

subroutine dump_array_1d_c_sp(A, nrows)
  complex(sp), intent(in)       :: A(:)
  integer, optional, intent(in) :: nrows
  !local variables
  integer                       :: nrows_, i
  integer, parameter            :: maxlen = 10
  
  if (present(nrows)) then
    if (nrows <= 0) then
      nrows_ = size(A)
    else
      nrows_ = min(nrows, size(A))
    end if
  else
    nrows_ = min(maxlen, size(A))
  end if

  do i = 1, nrows_
    write(*,"(2f10.4)") A(i)
  end do

  write(*,*)
end subroutine dump_array_1d_c_sp

subroutine dump_array_1d_c_dp(A, nrows)
  complex(dp), intent(in)          :: A(:)
  integer, optional, intent(in) :: nrows
  !local variables
  integer                       :: nrows_, i
  integer, parameter            :: maxlen = 10
  
  if (present(nrows)) then
    if (nrows <= 0) then
      nrows_ = size(A)
    else
      nrows_ = min(nrows, size(A))
    end if
  else
    nrows_ = min(maxlen, size(A))
  end if

  do i = 1, nrows_
    write(*,"(2f10.4)") A(i)
  end do

  write(*,*)
end subroutine dump_array_1d_c_dp

subroutine dump_array_2d_int(A, ncols, nrows)
  integer, intent(in)           :: A(:,:)
  integer, optional, intent(in) :: ncols
  integer, optional, intent(in) :: nrows
  !local variables
  integer                       :: ncols_, nrows_, i, j
  integer, parameter            :: maxlen = 10
  
  if (present(ncols)) then
    if (ncols <= 0) then
      ncols_ = size(A,2)
    else
      ncols_ = min(ncols, size(A,2))
    end if
  else
    ncols_ = min(maxlen, size(A,2))
  end if

  if (present(nrows)) then
    if (nrows <= 0) then
      nrows_ = size(A,1)
    else
      nrows_ = min(nrows, size(A,1))
    end if
  else
    nrows_ = min(maxlen, size(A,1))
  end if

  do i = 1, nrows_
    write(*,"(100I8)") (A(i,j), j=1,ncols_)  
  end do

  write(*,*)
end subroutine dump_array_2d_int

subroutine dump_array_2d_r_sp(A, ncols, nrows)
  real(sp), intent(in)          :: A(:,:)
  integer, optional, intent(in) :: ncols
  integer, optional, intent(in) :: nrows
  !local variables
  integer                       :: ncols_, nrows_, i, j
  integer, parameter            :: maxlen = 10
  
  if (present(ncols)) then
    if (ncols <= 0) then
      ncols_ = size(A,2)
    else
      ncols_ = min(ncols, size(A,2))
    end if
  else
    ncols_ = min(maxlen, size(A,2))
  end if

  if (present(nrows)) then
    if (nrows <= 0) then
      nrows_ = size(A,1)
    else
      nrows_ = min(nrows, size(A,1))
    end if
  else
    nrows_ = min(maxlen, size(A,1))
  end if

  do i = 1, nrows_
    write(*,"(100F10.4)") (A(i,j), j=1,ncols_)  
  end do

  write(*,*)
end subroutine dump_array_2d_r_sp

subroutine dump_array_2d_r_dp(A, ncols, nrows)
  real(dp), intent(in)          :: A(:,:)
  integer, optional, intent(in) :: ncols
  integer, optional, intent(in) :: nrows
  !local variables
  integer                       :: ncols_, nrows_, i, j
  integer, parameter            :: maxlen = 10
  
  if (present(ncols)) then
    if (ncols <= 0) then
      ncols_ = size(A,2)
    else
      ncols_ = min(ncols, size(A,2))
    end if
  else
    ncols_ = min(maxlen, size(A,2))
  end if

  if (present(nrows)) then
    if (nrows <= 0) then
      nrows_ = size(A,1)
    else
      nrows_ = min(nrows, size(A,1))
    end if
  else
    nrows_ = min(maxlen, size(A,1))
  end if

  do i = 1, nrows_
    write(*,"(100F10.4)") (A(i,j), j=1,ncols_)  
  end do

  write(*,*)
end subroutine dump_array_2d_r_dp

subroutine dump_array_2d_c_sp(A, ncols, nrows)
  complex(sp), intent(in)       :: A(:,:)
  integer, optional, intent(in) :: ncols
  integer, optional, intent(in) :: nrows
  
  write(*,*) "Real part:"
  call dump_array(real(A, sp), ncols, nrows)
  write(*,*) "Imaginary part:"
  call dump_array(imag(A), ncols, nrows)
  write(*,*)
end subroutine dump_array_2d_c_sp

subroutine dump_array_2d_c_dp(A, ncols, nrows)
  complex(dp), intent(in)    :: A(:,:)
  integer, optional, intent(in) :: ncols
  integer, optional, intent(in) :: nrows
  
  write(*,*) "Real part:"
  call dump_array(real(A, dp), ncols, nrows)
  write(*,*) "Imaginary part:"
  call dump_array(imag(A), ncols, nrows)
  write(*,*)
end subroutine dump_array_2d_c_dp


!******************************************************************************** 
!
! Report real-valued deterministic energy
!
!******************************************************************************** 
subroutine report_real_energy(comm, fh, comment, energy, method)
  type(mpi_communicator), intent(in)     :: comm
  type(FileHandle), intent(in)           :: fh
  character(len=*), intent(in)           :: comment
  real(wp), intent(in)                   :: energy
  character(len=*), optional, intent(in) :: method
  !local variables
  character(len=:), allocatable          :: method_

  if (present(method)) then
    method_ = method // "_"
  else 
    method_ = ""
  end if

  if (comm%mpirank == 0) then
    write(fh%funit,*)
    write(fh%funit,*)   comment // ":"
    write(fh%funit,*)   "---------------------------------------------------------"
    write(fh%funit,100) "contribution", "     Re E      "
    write(fh%funit,100) "------------", "     ----      "
    write(fh%funit,101) method_ // "etot =", energy
  end if

  100 format (1x,a,t30,a)
  101 format (1x,a,t30,f14.8)
end subroutine report_real_energy

!******************************************************************************** 
!
! Report complex-valued deterministic energy
!
!******************************************************************************** 
subroutine report_cmplx_energy(comm, fh, comment, energy, method)
  type(mpi_communicator), intent(in)     :: comm
  type(FileHandle), intent(in)           :: fh
  character(len=*), intent(in)           :: comment
  complex(wp), intent(in)                :: energy
  character(len=*), optional, intent(in) :: method
  !local variables
  character(len=:), allocatable          :: method_

  if (present(method)) then
    method_ = method // "_"
  else 
    method_ = ""
  end if

  if (comm%mpirank == 0) then
    write(fh%funit,*)
    write(fh%funit,*)   comment // ":"
    write(fh%funit,*)   "---------------------------------------------------------"
    write(fh%funit,100) "contribution", "     Re E      ", "    Im E    "
    write(fh%funit,100) "------------", "     ----      ", "    ----    "
    write(fh%funit,101) method_ // "etot =", real(energy,wp), aimag(energy)
  end if

  100 format (1x,a,t30,2(a,2x))
  101 format (1x,a,t30,f14.8,2x,f12.8)
end subroutine report_cmplx_energy

!******************************************************************************** 
!
! Report complex-valued energy and its standard deviation
!
!******************************************************************************** 
subroutine report_cmplx_std_energy(comm, fh, comment, energy, stddev, corr_len, method)
  type(mpi_communicator), intent(in)     :: comm
  type(FileHandle), intent(in)           :: fh
  character(len=*), intent(in)           :: comment
  complex(wp), intent(in)                :: energy
  real(wp), intent(in)                   :: stddev
  real(wp), optional, intent(in)         :: corr_len  
  character(len=*), optional, intent(in) :: method
  !local variables
  character(len=:), allocatable          :: method_

  if (present(method)) then
    method_ = method // "_"
  else 
    method_ = ""
  end if

  if (comm%mpirank == 0) then
    write(fh%funit,*)
    write(fh%funit,*)   comment // ":"
    write(fh%funit,*)   "---------------------------------------------------------"

    if (present(corr_len)) then
      write(fh%funit,100) "contribution", "     Re E      ", "    Im E    ", "    stddev    ", "corr. length", " stddev final "
      write(fh%funit,100) "------------", "     ----      ", "    ----    ", "    ------    ", "------------", " ------------ "
      write(fh%funit,101) method_ // "etot =", real(energy,wp), aimag(energy), stddev, corr_len, stddev*sqrt(corr_len)
    else 
      write(fh%funit,102) "contribution", "     Re E      ", "    Im E    ", "    stddev    "
      write(fh%funit,102) "------------", "     ----      ", "    ----    ", "    ------    "
      write(fh%funit,103) method_ // "etot =", real(energy,wp), aimag(energy), stddev
    end if
  end if

  100 format (1x,a,t30,5(a,2x))
  101 format (1x,a,t30,f14.8,2x,f12.8,2x,es14.6,2x,2x,f8.2,2x,2x,es14.6)
  102 format (1x,a,t30,3(a,2x))
  103 format (1x,a,t30,f14.8,2x,f12.8,2x,es14.6)
end subroutine report_cmplx_std_energy


!******************************************************************************** 
!
! Robust equivalence for real/complex numbers
!
!******************************************************************************** 
elemental function near_zero_r_sp(x, y, tol) result(is_near)
  real(sp), intent(in)           :: x
  real(sp), intent(in)           :: y
  real(sp), optional, intent(in) :: tol
  logical                           :: is_near
  !local variables
  integer, parameter                :: fac = 5
  real(sp)                       :: tol_

  tol_ = fac * tiny(tol_)

  if (present(tol)) then
    if (tol > tol_) tol_ = tol
  end if

  is_near = abs(x-y) < tol_
end function near_zero_r_sp

elemental function near_zero_r_dp(x, y, tol) result(is_near)
  real(dp), intent(in)           :: x
  real(dp), intent(in)           :: y
  real(dp), optional, intent(in) :: tol
  logical                           :: is_near
  !local variables
  integer, parameter                :: fac = 5
  real(dp)                       :: tol_

  tol_ = fac * tiny(tol_)

  if (present(tol)) then
    if (tol > tol_) tol_ = tol
  end if

  is_near = abs(x-y) < tol_
end function near_zero_r_dp

elemental function near_zero_c_sp(x, y, tol) result(is_near)
  complex(sp), intent(in)        :: x
  complex(sp), intent(in)        :: y
  real(sp), optional, intent(in) :: tol
  logical                           :: is_near
  !local variables
  integer, parameter                :: fac = 5
  real(sp)                       :: tol_

  tol_ = fac * tiny(tol_)

  if (present(tol)) then
    if (tol > tol_) tol_ = tol
  end if

  is_near = abs(x-y) < tol_
end function near_zero_c_sp

elemental function near_zero_c_dp(x, y, tol) result(is_near)
  complex(dp), intent(in)        :: x
  complex(dp), intent(in)        :: y
  real(dp), optional, intent(in) :: tol
  logical                        :: is_near
  !local variables
  integer, parameter             :: fac = 5
  real(dp)                       :: tol_

  tol_ = fac * tiny(tol_)

  if (present(tol)) then
    if (tol > tol_) tol_ = tol
  end if

  is_near = abs(x-y) < tol_
end function near_zero_c_dp

end module standalone
