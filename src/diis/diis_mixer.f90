! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module diis_mixer

use constants
use diis_config, only: DiisConfig

implicit none

private
public DiisMixer, DiisConfig

type DiisMixer
  integer               :: param_size 
  type(DiisConfig)      :: cfg

  integer               :: step
  integer               :: indx
  integer               :: diis_size

  real(wp), allocatable :: x(:,:)          !parameter vectors
  real(wp), allocatable :: e(:,:)          !error vectors
  real(wp), allocatable :: B(:,:)          !diis matrix
contains
  !generic               :: assignment(=) => assign_diis_mixer
  !procedure, private    :: assign_diis_mixer

  procedure, private    :: alloc => allocate_diis_arrays

  procedure             :: is_ready => is_diis_ready
  procedure             :: get_diis_size
  procedure, private    :: get_diis_indx

  procedure             :: update => update_diis
  procedure, private    :: update_diis_indx
  procedure, private    :: update_diis_space
  procedure, private    :: update_diis_matrix

  procedure             :: solve => solve_diis
  procedure             :: rms_error => calc_diis_rms_error
!debug: diis: to be implemented
  !procedure             :: sizeof => sizeof_diis_mixer
  !procedure             :: report => report_diis_mixer
end type DiisMixer

interface DiisMixer
  module procedure finit_diis_mixer
end interface DiisMixer

contains

!******************************************************************************** 
!
! Initialization of the DiisMixer object
! 
!******************************************************************************** 
subroutine init_diis_mixer(self, param_size, cfg)
  type(DiisMixer), intent(out) :: self
  integer, intent(in)          :: param_size
  type(DiisConfig), intent(in) :: cfg

  self%param_size = param_size
  self%cfg = cfg

  !keeps track of the SCF steps
  self%step = 0

  !indx shows at which row/column to store new vectors
  self%indx = 0

  !actual number of stored DIIS parameter and error vectors
  self%diis_size = 0

  !allocate DIIS arrays
  call self%alloc()
end subroutine init_diis_mixer


!******************************************************************************** 
!
! DiisMixer constructor
!
!******************************************************************************** 
function finit_diis_mixer(param_size, cfg) result(self)
  integer, intent(in)          :: param_size
  type(DiisConfig), intent(in) :: cfg
  type(DiisMixer)              :: self

  call init_diis_mixer(self, param_size, cfg)
end function finit_diis_mixer


!******************************************************************************** 
!
! Allocate arrays needed to store DIIS parameter and error vectors
!
!******************************************************************************** 
subroutine allocate_diis_arrays(self)
  class(DiisMixer), intent(inout) :: self
  !local
  integer                         :: param_size, max_diis_size

  param_size = self%param_size
  max_diis_size = self%cfg%max_diis_size

  !parameter vectors
  if (allocated(self%x)) deallocate(self%x)
  allocate(self%x(param_size, max_diis_size))
  self%x = 0.0_wp

  !error vectors
  if (allocated(self%e)) deallocate(self%e)
  allocate(self%e(param_size, max_diis_size))
  self%e = 0.0_wp

  !DIIS matrix
  if (allocated(self%B)) deallocate(self%B)
  allocate(self%B(max_diis_size, max_diis_size))
  self%B = 0.0_wp
end subroutine allocate_diis_arrays


!******************************************************************************** 
!
! Assignment operator for the DiisMixer object
!
!******************************************************************************** 
!subroutine assign_diis_mixer(to, from)
!  class(DiisMixer), intent(inout) :: to
!  type(DiisMixer), intent(in)     :: from
!
!  to%max_diis_size = from%max_diis_size
!  to%param_size = from%param_size
!  to%indx = from%indx
!  to%diis_size = from%diis_size
!
!  if (allocated(from%x)) then
!    if (allocated(to%x)) deallocate(to%x)
!    allocate(to%x, source=from%x)
!  end if
!
!  if (allocated(from%e)) then
!    if (allocated(to%e)) deallocate(to%e)
!    allocate(to%e, source=from%e)
!  end if
!
!  if (allocated(from%B)) then
!    if (allocated(to%B)) deallocate(to%B)
!    allocate(to%B, source=from%B)
!  end if
!end subroutine assign_diis_mixer


!******************************************************************************** 
!
! Returns the column index in the Diis buffer where to store new x, e vectors
!
!******************************************************************************** 
function get_diis_indx(self) result(indx)
  class(DiisMixer), intent(in) :: self
  integer                      :: indx
  
  indx = self%indx
end function get_diis_indx


!******************************************************************************** 
!
! Returns the current number of stored DIIS vectors
!
!******************************************************************************** 
function get_diis_size(self) result(diis_size)
  class(DiisMixer), intent(in) :: self
  integer                      :: diis_size
  
  if (self%cfg%enabled) then
    diis_size = self%diis_size
  else
    diis_size = 0
  end if
end function get_diis_size


!******************************************************************************** 
!
! Check is DIIS ready for extrapolating
!
!******************************************************************************** 
function is_diis_ready(self) result(diis_ready)
  class(DiisMixer), intent(in) :: self
  logical                      :: diis_ready

  diis_ready = self%cfg%enabled .and. (self%step >= self%cfg%start_solving_step)
end function is_diis_ready


!******************************************************************************** 
!
! Main routine in DIIS module:
!
!    1) updates step variable
!    2) Store new parameter and error vectors in DiisMixer object
!    3) Update DIIS matrix B
!
! Note: at the moment soft reset is implemented using the ring buffer:
!       the oldest parameter vector is replaced with the new one
!
!******************************************************************************** 
subroutine update_diis(self, x, e)
  class(DiisMixer), intent(inout) :: self
  real(wp), intent(in)            :: x(:)
  real(wp), intent(in)            :: e(:)

  self%step = self%step + 1

  !start storing at step start_sampling_step
  if (self%step >= self%cfg%start_sampling_step) then
    call self%update_diis_indx()
    call self%update_diis_space(x, e)
    call self%update_diis_matrix()
  end if
end subroutine update_diis


!******************************************************************************** 
!
! Determine 
!    1) at whihc index to store DIIS parameter and error vectors
!                                 &
!    2) actual size of the DIIS space
!
!******************************************************************************** 
subroutine update_diis_indx(self)
  class(DiisMixer), intent(inout) :: self

  self%indx = mod(self%indx+1, self%cfg%max_diis_size)
  if (self%indx == 0) self%indx = self%cfg%max_diis_size

  self%diis_size = self%diis_size + 1
  if (self%diis_size > self%cfg%max_diis_size) self%diis_size = self%cfg%max_diis_size
end subroutine update_diis_indx


!******************************************************************************** 
!
! Store new parameter and error vector into the DiisMixer object
!
!******************************************************************************** 
subroutine update_diis_space(self, x, e)
  class(DiisMixer), intent(inout) :: self
  real(wp), intent(in)            :: x(:)
  real(wp), intent(in)            :: e(:)
  !local
  integer                         :: indx

  !determine the index of oldest DIIS vectors
  !and overwrite them with the new ones
  indx = self%get_diis_indx()

  self%x(:,indx) = x
  self%e(:,indx) = e
end subroutine update_diis_space


!******************************************************************************** 
!
! Update diis matrix 
!
!    B_{ij} = <e_i | e_j>  
!
! A new column and row are added at the position 
! stored in DiisMixer%get_diis_indx()
!
!******************************************************************************** 
subroutine update_diis_matrix(self)
  class(DiisMixer), intent(inout) :: self
  !local
  integer                         :: i, indx, diis_size
  
  indx = self%get_diis_indx()
  diis_size = self%get_diis_size()

  do i = 1, diis_size
    self%B(i, indx) = dot_product(self%e(:,i), self%e(:,indx))
    self%B(indx, i) = self%B(i, indx)
  end do
end subroutine update_diis_matrix


!******************************************************************************** 
!
! Solve DIIS equation
!
!    B c = g
!
!
!        |              1  |       | c_1 |       | 0  |
!        |    B_ij      1  |       | c_2 |       | 0  |
!    B = |             ... |,  c = | ... |,  g = | ...|
!        |             ... |       | c_n |       | ...| 
!        |1   1  ...    0  |       | lam |       | 1  |
!
! and construct the new optimal parameter vector
!
!    x = \sum_i c_i x_i
!
!******************************************************************************** 
subroutine solve_diis(self, x)
  use lapack, only: inverse
  class(DiisMixer), intent(in) :: self
  real(wp), intent(out)        :: x(:)
  !local
  integer                      :: i, diis_size
  real(wp), parameter          :: reg=1.0E-10_wp
  real(wp), allocatable        :: B(:,:)
  real(wp), allocatable        :: g(:), c(:)

  !early return if still not ready for solving DIIS
  if (.not. self%is_ready()) return

  diis_size = self%get_diis_size()

  !setup full B matrix
  allocate(B(diis_size+1, diis_size+1))
  B(1:diis_size,1:diis_size) = self%B(1:diis_size, 1:diis_size)
  B(diis_size+1,1:diis_size) = 1.0_wp
  B(1:diis_size,diis_size+1) = 1.0_wp
  B(diis_size+1,diis_size+1) = 0.0_wp

  !regularize B matrix to avoid singularities
  do i = 1, diis_size
    B(i,i) = B(i,i) + reg
  end do

  !setup g vector
  allocate(g(diis_size+1))
  g = 0.0_wp
  g(diis_size+1) = 1.0_wp

  allocate(c(diis_size+1))
  c = 0.0_wp

  !solve Bc = g
  call inverse(B)
  c = matmul(B, g)

  x = 0.0_wp
  do i = 1, diis_size
    x = x + c(i) * self%x(:,i)
  end do
end subroutine solve_diis


!******************************************************************************** 
!
! Calculate RMS Error of the most recent DIIS error vector
!
!******************************************************************************** 
function calc_diis_rms_error(self) result(diis_err)
  class(DiisMixer), intent(in) :: self
  real(wp)                     :: diis_err
  !local
  integer                      :: indx

  indx = self%get_diis_indx()

  if (indx == 0) then
    !still no error vectors in the buffer
    diis_err = 0.0_wp
  else
    diis_err = norm2(self%e(:,indx)) / sqrt(real(size(self%e(:,indx)), wp))
  end if
end function calc_diis_rms_error

end module diis_mixer