! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

!******************************************************************************** 
!
!       Sparse module - implements basic sparse matrix functions
!
!******************************************************************************** 
module sparse

use constants

implicit none

type sparse_matrix_base
  integer                      :: size_
  integer                      :: size_max
  integer                      :: width
  integer                      :: nrows
  integer                      :: ncols
  real(spprec)                 :: tol = 1.0e-09_spprec
  real(spprec)                 :: sparsity
  integer, allocatable         :: row(:)
  integer, allocatable         :: col(:)
  integer, allocatable         :: indx(:)
end type sparse_matrix_base

type, extends(sparse_matrix_base) :: sparse_matrix
  real(dp), allocatable :: val(:,:)
end type sparse_matrix

type, extends(sparse_matrix_base) :: sparse_matrix_c
  complex(dp), allocatable :: val(:,:)
end type sparse_matrix_c

interface sparse_matrix
  procedure :: init_sparse, init_sparse_space
end interface sparse_matrix

interface sparse_matrix_c
  procedure :: init_sparse_c, init_sparse_space_c
end interface sparse_matrix_c

interface operator(+)
  procedure :: add, add_c, add_rc, add_cr
end interface

interface operator(*)
  procedure :: left_scalar, right_scalar, left_array, right_array
  procedure :: left_scalar_c, right_scalar_c, left_array_c, right_array_c
  procedure :: left_scalar_cr, right_scalar_rc, left_array_cr, right_array_rc
end interface

interface assignment(=)
  procedure :: copy_, copy_c, copy_cr
end interface

interface alloc
  procedure :: allocate_, allocate_src, allocate_c, allocate_src_c, allocate_src_cr
end interface alloc

interface add_to
  procedure :: add_sparse_to_full, add_sparse_to_full_array, add_sparse_c_to_full, add_sparse_c_to_full_array
  procedure :: add_sparse_r_to_full, add_sparse_r_to_full_array
end interface add_to

interface copy
  procedure :: copy_, copytrim, copy_c, copytrim_c
end interface copy

interface get_full_matrix
  procedure :: get_full_matrix_, get_full_matrix_c, get_full_matrix_width, get_full_matrix_width_c
end interface get_full_matrix

contains

!******************************************************************************** 
!       
!       Initialize sparse matrix:        
!               init_sparse_sp, init_sparse_dp - 
!                       initialize sparse matrices by some full matrix
!               init_sapse_space_sp, init_sparse_space_dp
!                       initialize space for some sparse matrix
!               assign_sparse_sp, assign_sparse_dp
!                       implement a_sparse = b_sparse
!******************************************************************************** 
type(sparse_matrix) function init_sparse(matrix, tol, width)
  real(spprec), intent(in)            :: matrix(:,:)
  real(spprec), optional, intent(in)  :: tol
  integer, optional, intent(in)       :: width
  !local variables
  integer                             :: i, j, indx
  !
  init_sparse%width = 1
  if (present(width)) init_sparse%width = width
  if (present(tol)) init_sparse%tol = tol
  init_sparse%size_max = size(matrix)
  init_sparse%nrows = size(matrix, 1)
  init_sparse%ncols = size(matrix, 2)
  init_sparse%size_ = size(pack(matrix, abs(matrix)>=init_sparse%tol))
  call alloc(init_sparse)
  !
  indx = 0
  do j = 1, size(matrix, 2)
    do i = 1, size(matrix, 1)
      if (abs(matrix(i,j)) < init_sparse%tol) cycle
      indx = indx + 1
      init_sparse%row(indx) = i
      init_sparse%col(indx) = j
      init_sparse%indx(indx) = size(matrix, 1) * (j - 1) + i
      init_sparse%val(1, indx) = matrix(i,j)
    end do
  end do 
  !
  init_sparse%sparsity = real(init_sparse%size_, spprec) / real(init_sparse%size_max, spprec) 
end function init_sparse
 
type(sparse_matrix_c) function init_sparse_c(matrix, tol, width)
  complex(spprec), intent(in)           :: matrix(:,:)
  real(spprec), optional, intent(in)    :: tol
  integer, optional, intent(in)         :: width
  !local variables
  integer                               :: i, j, indx
  !
  init_sparse_c%width = 1
  if (present(width)) init_sparse_c%width = width
  if (present(tol)) init_sparse_c%tol = tol
  init_sparse_c%size_max = size(matrix)
  init_sparse_c%nrows = size(matrix, 1)
  init_sparse_c%ncols = size(matrix, 2)
  init_sparse_c%size_ = size(pack(matrix, abs(matrix)>=init_sparse_c%tol))
  call alloc(init_sparse_c)
  !
  indx = 0
  do j = 1, size(matrix, 2)
    do i = 1, size(matrix, 1)
      if (abs(matrix(i,j)) < init_sparse_c%tol) cycle
      indx = indx + 1
      init_sparse_c%row(indx) = i
      init_sparse_c%col(indx) = j
      init_sparse_c%indx(indx) = size(matrix, 1) * (j - 1) + i
      init_sparse_c%val(1, indx) = matrix(i,j)
    end do
  end do 
  !
  init_sparse_c%sparsity = real(init_sparse_c%size_, spprec) / real(init_sparse_c%size_max, spprec) 
end function init_sparse_c

type(sparse_matrix) function init_sparse_space(n, tol, width)
  integer, intent(in)           :: n
  real(spprec), intent(in)      :: tol
  integer, optional, intent(in) :: width
  !
  init_sparse_space%width = 1
  if (present(width)) init_sparse_space%width = width
  init_sparse_space%size_ = n
  init_sparse_space%size_max = n
  init_sparse_space%nrows = unset
  init_sparse_space%ncols = unset
  call alloc(init_sparse_space)
end function init_sparse_space

type(sparse_matrix_c) function init_sparse_space_c(n, tol, width)
  integer, intent(in)           :: n
  real(spprec), intent(in)      :: tol
  integer, optional, intent(in) :: width
  !
  init_sparse_space_c%width = 1
  if (present(width)) init_sparse_space_c%width = width
  init_sparse_space_c%size_ = n
  init_sparse_space_c%size_max = n
  init_sparse_space_c%nrows = unset
  init_sparse_space_c%ncols = unset
  call alloc(init_sparse_space_c)
end function init_sparse_space_c


!******************************************************************************** 
!
!       Get ful matrix from the sparse matrix
!
!******************************************************************************** 
subroutine get_full_matrix_(a_sp, a)
  type(sparse_matrix), intent(in)         :: a_sp
  real(spprec), allocatable, intent(out)  :: a(:,:)
  !local variables
  integer                                 :: i
  !
  if (.not. allocated(a)) allocate(a(a_sp%nrows, a_sp%ncols)) 
  if (size(a,1)/=a_sp%nrows .or. size(a,2)/=a_sp%ncols)  then
    write(*,*) "Sizes of the full and sparse matrices don't agree. Program stops now"
    call exit
  end if
  !
  a = 0.0_spprec
  do i = 1, a_sp%size_
    a(a_sp%row(i), a_sp%col(i)) = a_sp%val(1,i)
  end do
end subroutine get_full_matrix_

subroutine get_full_matrix_c(a_sp, a)
  type(sparse_matrix_c), intent(in)         :: a_sp
  complex(spprec), allocatable, intent(out) :: a(:,:)
  !local variables
  integer                                   :: i
  !
  if (.not. allocated(a)) allocate(a(a_sp%nrows, a_sp%ncols)) 
  if (size(a,1)/=a_sp%nrows .or. size(a,2)/=a_sp%ncols)  then
    write(*,*) "Sizes of the full and sparse matrices don't agree. Program stops now"
    call exit
  end if
  !
  a = (0.0_spprec, 0.0_spprec)
  do i = 1, a_sp%size_
    a(a_sp%row(i), a_sp%col(i)) = a_sp%val(1,i)
  end do
end subroutine get_full_matrix_c

subroutine get_full_matrix_width(a_sp, a)
  type(sparse_matrix), intent(in)        :: a_sp
  real(spprec), allocatable, intent(out) :: a(:,:,:)
  !local variables
  integer                                :: i
  !
  if (.not. allocated(a)) allocate(a(a_sp%nrows, a_sp%ncols, a_sp%width))
  if (size(a,1)/=a_sp%nrows .or. size(a,2)/=a_sp%ncols .or. size(a,3)/=a_sp%width)  then
    write(*,*) "Width of the sparse matrix doesn't match the third dimension of the target matrix. Program stops now."
    call exit
  end if
  !
  a = 0.0_spprec
  do i = 1, a_sp%size_
    a(a_sp%row(i),a_sp%col(i),:) = a_sp%val(:, i)
  end do
end subroutine get_full_matrix_width

subroutine get_full_matrix_width_c(a_sp, a)
  type(sparse_matrix_c), intent(in)         :: a_sp
  complex(spprec), allocatable, intent(out) :: a(:,:,:)
  !local variables
  integer                                   :: i
  !
  if (.not. allocated(a)) allocate(a(a_sp%nrows, a_sp%ncols, a_sp%width))
  if (size(a,1)/=a_sp%nrows .or. size(a,2)/=a_sp%ncols .or. size(a,3)/=a_sp%width)  then
    write(*,*) "Width of the sparse matrix doesn't match the third dimension of the target matrix. Program stops now."
    call exit
  end if
  !
  a = (0.0_spprec, 0.0_spprec)
  do i = 1, a_sp%size_
    a(a_sp%row(i),a_sp%col(i),:) = a_sp%val(:, i)
  end do
end subroutine get_full_matrix_width_c


!******************************************************************************** 
!
!       Copy of the sparse matrix, in other words:
!                       a = b
!
!******************************************************************************** 
subroutine copy_(a, b)
  type(sparse_matrix), intent(inout) :: a
  type(sparse_matrix), intent(in)    :: b
  !
  a%size_ = b%size_
  a%size_max = b%size_max
  a%width = b%width
  a%nrows = b%nrows
  a%ncols = b%ncols
  a%tol = b%tol
  a%sparsity = b%sparsity
  call alloc(a, b)
end subroutine copy_

subroutine copy_c(a, b)
  type(sparse_matrix_c), intent(inout) :: a
  type(sparse_matrix_c), intent(in)    :: b
  !
  a%size_ = b%size_
  a%size_max = b%size_max
  a%width = b%width
  a%nrows = b%nrows
  a%ncols = b%ncols
  a%tol = b%tol
  a%sparsity = b%sparsity
  call alloc(a, b)
end subroutine copy_c

subroutine copy_cr(a, b)
  type(sparse_matrix_c), intent(inout) :: a
  type(sparse_matrix), intent(in)      :: b
  !
  a%size_ = b%size_
  a%size_max = b%size_max
  a%width = b%width
  a%nrows = b%nrows
  a%ncols = b%ncols
  a%tol = b%tol
  a%sparsity = b%sparsity
  call alloc(a, b)
end subroutine copy_cr


!******************************************************************************** 
!
!       Copy part of the sparse matrix
!                       a = b(1:newsize)
!
!******************************************************************************** 
subroutine copytrim(a, b, newsize)
  type(sparse_matrix), intent(inout) :: a
  type(sparse_matrix), intent(in)    :: b
  integer, intent(in)                :: newsize
  !
  a%size_ = newsize
  a%size_max = b%size_max
  a%width = b%width
  a%nrows = b%nrows
  a%ncols = b%ncols
  a%tol = b%tol
  a%sparsity = real(a%size_, spprec) / real(a%size_max, spprec)
  call alloc(a, b)
end subroutine copytrim

subroutine copytrim_c(a, b, newsize)
  type(sparse_matrix_c), intent(inout) :: a
  type(sparse_matrix_c), intent(in)    :: b
  integer, intent(in)                  :: newsize
  !
  a%size_ = newsize
  a%size_max = b%size_max
  a%width = b%width
  a%nrows = b%nrows
  a%ncols = b%ncols
  a%tol = b%tol
  a%sparsity = real(a%size_, spprec) / real(a%size_max, spprec)
  call alloc(a, b)
end subroutine copytrim_c


!******************************************************************************** 
!
!       Allocate space for sparse matrix
!
!******************************************************************************** 
subroutine allocate_(a, size_)
  type(sparse_matrix), intent(inout) :: a
  integer, optional, intent(in)      :: size_
  !
  if (present(size_)) a%size_ = size_
  if (.not. allocated(a%indx)) allocate(a%indx(a%size_))
  if (.not. allocated(a%row)) allocate(a%row(a%size_))
  if (.not. allocated(a%col)) allocate(a%col(a%size_))
  if (.not. allocated(a%val)) allocate(a%val(a%width, a%size_))
end subroutine allocate_

subroutine allocate_c(a, size_)
  type(sparse_matrix_c), intent(inout) :: a
  integer, optional, intent(in)        :: size_
  !
  if (present(size_)) a%size_ = size_
  if (.not. allocated(a%indx)) allocate(a%indx(a%size_))
  if (.not. allocated(a%row)) allocate(a%row(a%size_))
  if (.not. allocated(a%col)) allocate(a%col(a%size_))
  if (.not. allocated(a%val)) allocate(a%val(a%width, a%size_))
end subroutine allocate_c

subroutine allocate_src(a, b, size_)
  type(sparse_matrix), intent(inout) :: a
  type(sparse_matrix), intent(in)    :: b
  integer, optional, intent(in)      :: size_
  !
  if (present(size_)) a%size_ = size_
  call alloc(a)
  a%indx = b%indx(1:a%size_)
  a%row = b%row(1:a%size_)
  a%col = b%col(1:a%size_)
  a%val = b%val(:,1:a%size_)
end subroutine allocate_src

subroutine allocate_src_c(a, b, size_)
  type(sparse_matrix_c), intent(inout) :: a
  type(sparse_matrix_c), intent(in)    :: b
  integer, optional, intent(in)        :: size_
  !
  if (present(size_)) a%size_ = size_
  call alloc(a)
  a%indx = b%indx(1:a%size_)
  a%row = b%row(1:a%size_)
  a%col = b%col(1:a%size_)
  a%val = b%val(:,1:a%size_)
end subroutine allocate_src_c

subroutine allocate_src_cr(a, b, size_)
  type(sparse_matrix_c), intent(inout) :: a
  type(sparse_matrix), intent(in)      :: b
  integer, optional, intent(in)        :: size_
  !
  if (present(size_)) a%size_ = size_
  call alloc(a)
  a%indx = b%indx(1:a%size_)
  a%row = b%row(1:a%size_)
  a%col = b%col(1:a%size_)
  a%val = b%val(:,1:a%size_)
end subroutine allocate_src_cr


!******************************************************************************** 
!
!       Addition of two sparse matrices
!               Overload operator (+) 
!               Usage: c= a + b
!
!******************************************************************************** 
type(sparse_matrix) function add(a, b)
  type(sparse_matrix), intent(in) :: a
  type(sparse_matrix), intent(in) :: b
  !local variables
  type(sparse_matrix)             :: temp
  integer                         :: i, j, k, rest
  real(spprec)                    :: diff = 1.0e-06_spprec
  !
  if (a%width/=b%width .or. a%nrows/=b%nrows .or. a%ncols/=b%ncols) then
    write(*,*) "Inconsistent sizes of sparse matrices"
    call exit
  end if
  if (abs(a%tol-b%tol) > diff) write(*,*) "Caution, two sparse matrices with different tolerances are added"
  !
  temp = sparse_matrix(a%size_+b%size_, min(a%tol,b%tol), a%width)
  i = 1
  j = 1 
  k = 1
  do while (i<=a%size_ .and. j<=b%size_)
    if (a%indx(i) == b%indx(j)) then
      temp%indx(k) = a%indx(i)
      temp%row(k)  = a%row(i)
      temp%col(k)  = a%col(i)
      temp%val(:,k)  = a%val(:,i) + b%val(:,j)
      i = i + 1
      j = j + 1
      k = k + 1 
    else if (a%indx(i) < b%indx(j)) then
      temp%indx(k) = a%indx(i)
      temp%row(k)  = a%row(i)
      temp%col(k)  = a%col(i)
      temp%val(:,k)  = a%val(:,i)
      i = i + 1
      k = k + 1
    else
      temp%indx(k) = b%indx(j)
      temp%row(k)  = b%row(j)
      temp%col(k)  = b%col(j)
      temp%val(:,k)  = b%val(:,j)
      j = j + 1
      k = k + 1
    end if
  end do 
  !Copy the rest of matrix a if not done
  do rest = i, a%size_
    temp%indx(k) = a%indx(rest)
    temp%row(k)  = a%row(rest)
    temp%col(k)  = a%col(rest)
    temp%val(:,k)  = a%val(:,rest)
    k = k + 1
  end do
  !Copy the rest of matrix b if not done
  do rest = j, b%size_
    temp%indx(k) = b%indx(rest)
    temp%row(k)  = b%row(rest)
    temp%col(k)  = b%col(rest)
    temp%val(:,k)  = b%val(:,rest)
    k = k + 1
  end do
  ! 
  temp%nrows = a%nrows
  temp%ncols = a%ncols
  call copy(add, temp, k-1)
end function add

type(sparse_matrix_c) function add_c(a, b)
  type(sparse_matrix_c), intent(in) :: a
  type(sparse_matrix_c), intent(in) :: b
  !local variables
  type(sparse_matrix_c)             :: temp
  integer                           :: i, j, k, rest
  real(spprec)                      :: diff = 1.0e-06_spprec
  !
  if (a%width/=b%width .or. a%nrows/=b%nrows .or. a%ncols/=b%ncols) then
    write(*,*) "Inconsistent sizes of sparse matrices"
    call exit
  end if
  if (abs(a%tol-b%tol) > diff) write(*,*) "Caution, two sparse matrices with different tolerances are added"
  !
  temp = sparse_matrix_c(a%size_+b%size_, min(a%tol,b%tol), a%width)
  i = 1
  j = 1 
  k = 1
  do while (i<=a%size_ .and. j<=b%size_)
    if (a%indx(i) == b%indx(j)) then
      temp%indx(k) = a%indx(i)
      temp%row(k)  = a%row(i)
      temp%col(k)  = a%col(i)
      temp%val(:,k)  = a%val(:,i) + b%val(:,j)
      i = i + 1
      j = j + 1
      k = k + 1 
    else if (a%indx(i) < b%indx(j)) then
      temp%indx(k) = a%indx(i)
      temp%row(k)  = a%row(i)
      temp%col(k)  = a%col(i)
      temp%val(:,k)  = a%val(:,i)
      i = i + 1
      k = k + 1
    else
      temp%indx(k) = b%indx(j)
      temp%row(k)  = b%row(j)
      temp%col(k)  = b%col(j)
      temp%val(:,k)  = b%val(:,j)
      j = j + 1
      k = k + 1
    end if
  end do 
  !Copy the rest of matrix a if not done
  do rest = i, a%size_
    temp%indx(k) = a%indx(rest)
    temp%row(k)  = a%row(rest)
    temp%col(k)  = a%col(rest)
    temp%val(:,k)  = a%val(:,rest)
    k = k + 1
  end do
  !Copy the rest of matrix b if not done
  do rest = j, b%size_
    temp%indx(k) = b%indx(rest)
    temp%row(k)  = b%row(rest)
    temp%col(k)  = b%col(rest)
    temp%val(:,k)  = b%val(:,rest)
    k = k + 1
  end do
  ! 
  temp%nrows = a%nrows
  temp%ncols = a%ncols
  call copy(add_c, temp, k-1)
end function add_c

type(sparse_matrix_c) function add_rc(a, b)
  type(sparse_matrix), intent(in)   :: a
  type(sparse_matrix_c), intent(in) :: b
  !local variables
  type(sparse_matrix_c)             :: temp
  integer                           :: i, j, k, rest
  real(spprec)                      :: diff = 1.0e-06_spprec
  !
  if (a%width/=b%width .or. a%nrows/=b%nrows .or. a%ncols/=b%ncols) then
    write(*,*) "Inconsistent sizes of sparse matrices"
    call exit
  end if
  if (abs(a%tol-b%tol) > diff) write(*,*) "Caution, two sparse matrices with different tolerances are added"
  !
  temp = sparse_matrix_c(a%size_+b%size_, min(a%tol,b%tol), a%width)
  i = 1
  j = 1 
  k = 1
  do while (i<=a%size_ .and. j<=b%size_)
    if (a%indx(i) == b%indx(j)) then
      temp%indx(k) = a%indx(i)
      temp%row(k)  = a%row(i)
      temp%col(k)  = a%col(i)
      temp%val(:,k)  = a%val(:,i) + b%val(:,j)
      i = i + 1
      j = j + 1
      k = k + 1 
    else if (a%indx(i) < b%indx(j)) then
      temp%indx(k) = a%indx(i)
      temp%row(k)  = a%row(i)
      temp%col(k)  = a%col(i)
      temp%val(:,k)  = a%val(:,i)
      i = i + 1
      k = k + 1
    else
      temp%indx(k) = b%indx(j)
      temp%row(k)  = b%row(j)
      temp%col(k)  = b%col(j)
      temp%val(:,k)  = b%val(:,j)
      j = j + 1
      k = k + 1
    end if
  end do 
  !Copy the rest of matrix a if not done
  do rest = i, a%size_
    temp%indx(k) = a%indx(rest)
    temp%row(k)  = a%row(rest)
    temp%col(k)  = a%col(rest)
    temp%val(:,k)  = a%val(:,rest)
    k = k + 1
  end do
  !Copy the rest of matrix b if not done
  do rest = j, b%size_
    temp%indx(k) = b%indx(rest)
    temp%row(k)  = b%row(rest)
    temp%col(k)  = b%col(rest)
    temp%val(:,k)  = b%val(:,rest)
    k = k + 1
  end do
  ! 
  temp%nrows = a%nrows
  temp%ncols = a%ncols
  call copy(add_rc, temp, k-1)
end function add_rc

type(sparse_matrix_c) function add_cr(a, b)
  type(sparse_matrix_c), intent(in) :: a
  type(sparse_matrix), intent(in)   :: b
  !local variables
  type(sparse_matrix_c)             :: temp
  integer                           :: i, j, k, rest
  real(spprec)                      :: diff = 1.0e-06_spprec
  !
  if (a%width/=b%width .or. a%nrows/=b%nrows .or. a%ncols/=b%ncols) then
    write(*,*) "Inconsistent sizes of sparse matrices"
    call exit
  end if
  if (abs(a%tol-b%tol) > diff) write(*,*) "Caution, two sparse matrices with different tolerances are added"
  !
  temp = sparse_matrix_c(a%size_+b%size_, min(a%tol,b%tol), a%width)
  i = 1
  j = 1 
  k = 1
  do while (i<=a%size_ .and. j<=b%size_)
    if (a%indx(i) == b%indx(j)) then
      temp%indx(k) = a%indx(i)
      temp%row(k)  = a%row(i)
      temp%col(k)  = a%col(i)
      temp%val(:,k)  = a%val(:,i) + b%val(:,j)
      i = i + 1
      j = j + 1
      k = k + 1 
    else if (a%indx(i) < b%indx(j)) then
      temp%indx(k) = a%indx(i)
      temp%row(k)  = a%row(i)
      temp%col(k)  = a%col(i)
      temp%val(:,k)  = a%val(:,i)
      i = i + 1
      k = k + 1
    else
      temp%indx(k) = b%indx(j)
      temp%row(k)  = b%row(j)
      temp%col(k)  = b%col(j)
      temp%val(:,k)  = b%val(:,j)
      j = j + 1
      k = k + 1
    end if
  end do 
  !Copy the rest of matrix a if not done
  do rest = i, a%size_
    temp%indx(k) = a%indx(rest)
    temp%row(k)  = a%row(rest)
    temp%col(k)  = a%col(rest)
    temp%val(:,k)  = a%val(:,rest)
    k = k + 1
  end do
  !Copy the rest of matrix b if not done
  do rest = j, b%size_
    temp%indx(k) = b%indx(rest)
    temp%row(k)  = b%row(rest)
    temp%col(k)  = b%col(rest)
    temp%val(:,k)  = b%val(:,rest)
    k = k + 1
  end do
  ! 
  temp%nrows = a%nrows
  temp%ncols = a%ncols
  call copy(add_cr, temp, k-1)
end function add_cr


!******************************************************************************** 
!
!       Multiplication by scalar:
!               Overload operator (*)
!                       b = lambda * a
!                       b = a * lambda 
!
!******************************************************************************** 
function left_scalar(x, a) result(b)
  real(spprec), intent(in)        :: x
  type(sparse_matrix), intent(in) :: a
  type(sparse_matrix)             :: b
  !local variables
  integer                         :: i
  !
  b = a
  do i = 1, a%size_
    b%val(:,i) = b%val(:,i) * x
  end do
end function left_scalar

function left_array(x, a) result(b)
  real(spprec), intent(in)        :: x(:)
  type(sparse_matrix), intent(in) :: a
  type(sparse_matrix)             :: b
  !local variables
  real(spprec), allocatable       :: val(:,:)
  integer                         :: i
  !
  if (a%width /= 1) then
    write(*,*) "Width of the sparse matrix in multiplication by array must be 1"
    call exit
  end if
  ! 
  b = a
  b%width = size(x)
  allocate(val(size(x), a%size_))
  do i = 1, a%size_
    val(:,i) = a%val(1,i) * x
  end do
  call move_alloc(val, b%val)
end function left_array

function right_scalar(a, x) result(b)
  type(sparse_matrix), intent(in) :: a
  real(spprec), intent(in)        :: x
  type(sparse_matrix)             :: b
  !local variables
  integer                         :: i
  !
  b = a
  do i = 1, a%size_
    b%val(:,i) = b%val(:,i) * x
  end do
end function right_scalar

function right_array(a, x) result(b)
  type(sparse_matrix), intent(in) :: a
  real(spprec), intent(in)        :: x(:)
  type(sparse_matrix)             :: b
  !local variables
  real(spprec), allocatable       :: val(:,:)
  integer                         :: i
  !
  if (a%width /= 1) then
    write(*,*) "Width of the sparse matrix in multiplication by array must be 1"
    call exit
  end if
  ! 
  b = a
  b%width = size(x)
  allocate(val(size(x), a%size_))
  do i = 1, a%size_
    val(:,i) = a%val(1,i) * x
  end do
  call move_alloc(val, b%val)
end function right_array

function left_scalar_c(x, a) result(b)
  complex(spprec), intent(in)       :: x
  type(sparse_matrix_c), intent(in) :: a
  type(sparse_matrix_c)             :: b
  !local variables
  integer                           :: i
  !
  b = a
  do i = 1, a%size_
    b%val(:,i) = b%val(:,i) * x
  end do
end function left_scalar_c

function left_array_c(x, a) result(b)
  complex(spprec), intent(in)       :: x(:)
  type(sparse_matrix_c), intent(in) :: a
  type(sparse_matrix_c)             :: b
  !local variables
  complex(spprec), allocatable      :: val(:,:)
  integer                           :: i
  !
  if (a%width /= 1) then
    write(*,*) "Width of the sparse matrix in multiplication by array must be 1"
    call exit
  end if
  ! 
  b = a
  b%width = size(x)
  allocate(val(size(x), a%size_))
  do i = 1, a%size_
    val(:,i) = a%val(1,i) * x
  end do
  call move_alloc(val, b%val)
end function left_array_c

function right_scalar_c(a, x) result(b)
  type(sparse_matrix_c), intent(in) :: a
  complex(spprec), intent(in)       :: x
  type(sparse_matrix_c)             :: b
  !local variables
  integer                           :: i
  !
  b = a
  do i = 1, a%size_
    b%val(:,i) = b%val(:,i) * x
  end do
end function right_scalar_c

function right_array_c(a, x) result(b)
  type(sparse_matrix_c), intent(in) :: a
  complex(spprec), intent(in)       :: x(:)
  type(sparse_matrix_c)             :: b
  !local variables
  complex(spprec), allocatable      :: val(:,:)
  integer                           :: i
  ! 
  if (a%width /= 1) then
    write(*,*) "Width of the sparse matrix in multiplication by array must be 1"
    call exit
  end if
  ! 
  b = a
  b%width = size(x)
  allocate(val(size(x), a%size_))
  do i = 1, a%size_
    val(:,i) = a%val(1,i) * x
  end do
  call move_alloc(val, b%val)
end function right_array_c

function left_scalar_cr(x, a) result(b)
  complex(spprec), intent(in)     :: x
  type(sparse_matrix), intent(in) :: a
  type(sparse_matrix_c)           :: b
  !local variables
  integer                         :: i
  !
  b = a
  do i = 1, a%size_
    b%val(:,i) = b%val(:,i) * x
  end do
end function left_scalar_cr

function left_array_cr(x, a) result(b)
  complex(spprec), intent(in)     :: x(:)
  type(sparse_matrix), intent(in) :: a
  type(sparse_matrix_c)           :: b
  !local variables
  complex(spprec), allocatable    :: val(:,:)
  integer                         :: i
  !
  if (a%width /= 1) then
    write(*,*) "Width of the sparse matrix in multiplication by array must be 1"
    call exit
  end if
  ! 
  b = a
  b%width = size(x)
  allocate(val(size(x), a%size_))
  do i = 1, a%size_
    val(:,i) = a%val(1,i) * x
  end do
  call move_alloc(val, b%val)
end function left_array_cr

function right_scalar_rc(a, x) result(b)
  type(sparse_matrix), intent(in) :: a
  complex(spprec), intent(in)     :: x
  type(sparse_matrix_c)           :: b
  !local variables
  integer                         :: i
  !
  b = a
  do i = 1, a%size_
    b%val(:,i) = b%val(:,i) * x
  end do
end function right_scalar_rc

function right_array_rc(a, x) result(b)
  type(sparse_matrix), intent(in) :: a
  complex(spprec), intent(in)     :: x(:)
  type(sparse_matrix_c)           :: b
  !local variables
  complex(spprec), allocatable    :: val(:,:)
  integer                         :: i
  ! 
  if (a%width /= 1) then
    write(*,*) "Width of the sparse matrix in multiplication by array must be 1"
    call exit
  end if
  ! 
  b = a
  b%width = size(x)
  allocate(val(size(x), a%size_))
  do i = 1, a%size_
    val(:,i) = a%val(1,i) * x
  end do
  call move_alloc(val, b%val)
end function right_array_rc


!******************************************************************************** 
!
!       Add sparse matrix to the full matrix
!
!******************************************************************************** 
subroutine add_sparse_to_full(a, a_, coeff)
  real(spprec), intent(inout)        :: a(:,:)
  type(sparse_matrix), intent(in)    :: a_
  real(spprec), optional, intent(in) :: coeff
  !local variables
  integer                            :: i
  real(spprec)                       :: coeff_
  !
  if (a_%width /= 1) then
    write(*,*) "Warning: add_to => add_sparse_to_full should use sparse matices with width 1 only"
    call exit
  end if
  !
  coeff_ = 1.0_spprec
  if (present(coeff)) coeff_ = coeff
  !
  do i = 1, a_%size_
    a(a_%row(i), a_%col(i)) = a(a_%row(i), a_%col(i)) + coeff_ * a_%val(1,i)
  end do
end subroutine add_sparse_to_full

subroutine add_sparse_to_full_array(a, a_, coeff)
  real(spprec), intent(inout)        :: a(:,:,:) 
  type(sparse_matrix), intent(in)    :: a_
  real(spprec), optional, intent(in) :: coeff(:)
  !local variables
  integer                            :: i
  real(spprec), allocatable          :: coeff_(:)
  !
  if (a_%width /= 1) then
    write(*,*) "Warning: add_to => add_sparse_to_full should use sparse matices with width 1 only"
    call exit
  end if
  !
  if (present(coeff)) then
    allocate(coeff_, source=coeff)
  else
    allocate(coeff_(size(a,1)))
    coeff_ = 1.0_spprec
  end if
  !
  do i = 1, a_%size_
    a(:, a_%row(i), a_%col(i)) = a(:, a_%row(i), a_%col(i)) + coeff_ * a_%val(1,i)
  end do
end subroutine add_sparse_to_full_array

subroutine add_sparse_c_to_full(a, a_, coeff)
  complex(spprec), intent(inout)        :: a(:,:)
  type(sparse_matrix_c), intent(in)     :: a_
  complex(spprec), optional, intent(in) :: coeff
  !local variables
  integer                               :: i
  complex(spprec)                       :: coeff_
  !
  if (a_%width /= 1) then
    write(*,*) "Warning: add_to => add_sparse_to_full should use sparse matices with width 1 only"
    call exit
  end if
  !
  coeff_ = (1.0_spprec, 0.0_spprec)
  if (present(coeff)) coeff_ = coeff
  !
  do i = 1, a_%size_
    a(a_%row(i), a_%col(i)) = a(a_%row(i), a_%col(i)) + coeff_ * a_%val(1,i)
  end do
end subroutine add_sparse_c_to_full

subroutine add_sparse_c_to_full_array(a, a_, coeff)
  complex(spprec), intent(inout)        :: a(:,:,:) 
  type(sparse_matrix_c), intent(in)     :: a_
  complex(spprec), optional, intent(in) :: coeff(:)
  !local variables
  integer                               :: i
  complex(spprec), allocatable          :: coeff_(:)
  !
  if (a_%width /= 1) then
    write(*,*) "Warning: add_to => add_sparse_to_full should use sparse matices with width 1 only"
    call exit
  end if
  !
  if (present(coeff)) then
    allocate(coeff_, source=coeff)
  else
    allocate(coeff_(size(a,1)))
    coeff_ = (1.0_spprec, 0.0_spprec)
  end if
  !
  do i = 1, a_%size_
    a(:, a_%row(i), a_%col(i)) = a(:, a_%row(i), a_%col(i)) + coeff_ * a_%val(1,i)
  end do
end subroutine add_sparse_c_to_full_array

subroutine add_sparse_r_to_full(a, a_, coeff)
  complex(spprec), intent(inout)        :: a(:,:)
  type(sparse_matrix), intent(in)       :: a_
  complex(spprec), optional, intent(in) :: coeff
  !local variables
  integer                               :: i
  complex(spprec)                       :: coeff_
  !
  if (a_%width /= 1) then
    write(*,*) "Warning: add_to => add_sparse_to_full should use sparse matices with width 1 only"
    call exit
  end if
  !
  coeff_ = (1.0_spprec, 0.0_spprec)
  if (present(coeff)) coeff_ = coeff
  !
  do i = 1, a_%size_
    a(a_%row(i), a_%col(i)) = a(a_%row(i), a_%col(i)) + coeff_ * a_%val(1,i)
  end do
end subroutine add_sparse_r_to_full

subroutine add_sparse_r_to_full_array(a, a_, coeff)
  complex(spprec), intent(inout)        :: a(:,:,:) 
  type(sparse_matrix), intent(in)       :: a_
  complex(spprec), optional, intent(in) :: coeff(:)
  !local variables
  integer                               :: i
  complex(spprec), allocatable          :: coeff_(:)
  !
  if (a_%width /= 1) then
    write(*,*) "Warning: add_to => add_sparse_to_full should use sparse matices with width 1 only"
    call exit
  end if
  !
  if (present(coeff)) then
    allocate(coeff_, source=coeff)
  else
    allocate(coeff_(size(a,1)))
    coeff_ = (1.0_spprec, 0.0_spprec)
  end if
  !
  do i = 1, a_%size_
    a(:, a_%row(i), a_%col(i)) = a(:, a_%row(i), a_%col(i)) + coeff_ * a_%val(1,i)
  end do
end subroutine add_sparse_r_to_full_array

end module sparse