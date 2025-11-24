! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module stat_analysis

  use constants
  use file_handle, only: FileHandle
  use mpi
  
  implicit none
    
  type stat_var
    real(wp) :: cnt
    real(wp) :: max_ 
    real(wp) :: min_
    real(wp) :: mean
    real(wp) :: std
    integer  :: min_index
    integer  :: max_index
  contains
    procedure :: init => init_stat_var
    generic :: update => update_stat_var, update_stat_var_1d
    procedure :: finalize => finalize_stat_var
    procedure :: print_ => print_stat_var
    procedure, private :: update_stat_var, update_stat_var_1d
  end type stat_var

  type stat_var_c
    real(wp)    :: cnt
    complex(wp) :: max_ 
    complex(wp) :: min_
    complex(wp) :: mean
    real(wp)    :: std
    integer     :: min_index
    integer     :: max_index
  contains
    procedure :: init => init_stat_var_c
    generic :: update => update_stat_var_c, update_stat_var_1d_c
    procedure :: finalize => finalize_stat_var_c
    procedure :: print_ => print_stat_var_c
    procedure, private :: update_stat_var_c, update_stat_var_1d_c
  end type stat_var_c
contains

!******************************************************************************** 
!
!       Initialize stat_var
!
!******************************************************************************** 
subroutine init_stat_var(self)
  class(stat_var), intent(inout) :: self

  self%cnt = 0.0_wp
  self%max_ = 0.0_wp
  self%min_ = 0.0_wp
  self%mean = 0.0_wp
  self%std = 0.0_wp
  self%min_index = 0
  self%max_index = 0
end subroutine init_stat_var    


!******************************************************************************** 
!
!       Update stat_var
!
!******************************************************************************** 
subroutine update_stat_var(self, x, step)
  class(stat_var), intent(inout) :: self
  real(wp), intent(in)           :: x
  integer, optional, intent(in)  :: step

  if (self%cnt < 0.5_wp) then
    self%cnt = 1.0_wp
    self%max_ = x
    self%max_index = 1
    self%min_ = x 
    self%min_index = 1
    self%mean = x
    self%std = x**2
    return
  end if

  self%cnt = self%cnt + 1.0_wp

  if (x > self%max_) then
    self%max_ = x
    if (present(step)) then
      self%max_index = step
    else
      self%max_index = nint(self%cnt)
    end if
  end if

  if (x < self%min_) then
    self%min_ = x
    if (present(step)) then
      self%min_index = step
    else
      self%min_index = self%cnt
    end if
  end if

  self%mean = self%mean + x 
  self%std = self%std + x**2
end subroutine update_stat_var

subroutine update_stat_var_1d(self, x, step)
  class(stat_var), intent(inout) :: self
  real(wp), intent(in)           :: x(:)
  integer, optional, intent(in)  :: step
  !local variables
  real(wp)                       :: xmax, xmin
  integer                        :: xmaxloc, xminloc, pos

  if (self%cnt < 0.5_wp) then
    self%cnt = real(size(x), wp)
    pos = maxloc(x, 1)
    self%max_ = x(pos)
    if (present(step)) then
      self%max_index = step
    else
      self%max_index = pos
    end if
    pos = minloc(x, 1)
    self%min_ = x(pos)
    if (present(step)) then
      self%min_index = step
    else
      self%min_index = pos
    end if 
    self%mean = sum(x)
    self%std = sum(x**2)
    return
  end if

  xmaxloc = maxloc(x, 1)
  xmax = x(xmaxloc)

  if (xmax > self%max_) then
    self%max_ = xmax
    if (present(step)) then
      self%max_index = step
    else
      self%max_index = self%cnt + xmaxloc
    end if
  end if

  xminloc = minloc(x, 1)
  xmin = x(xminloc)

  if (xmin < self%min_) then
    self%min_ = xmin
    if (present(step)) then
      self%min_index = step
    else
      self%min_index = self%cnt + xminloc
    end if
  end if

  self%cnt = self%cnt + size(x)
  self%mean = self%mean + sum(x)
  self%std = self%std + sum(x**2)
end subroutine update_stat_var_1d


!******************************************************************************** 
!
!       Finalize stat_var
!
!******************************************************************************** 
subroutine finalize_stat_var(self, comm)
  class(stat_var), intent(inout)               :: self
  type(mpi_communicator), optional, intent(in) :: comm 
  !local variables
  type(mpi_communicator)                       :: comm_
  real(wp), allocatable                        :: temp(:)
  integer, allocatable                         :: temp_ind(:)
  integer                                      :: pos

  if (present(comm)) then
    comm_ = comm
  else
    comm_ = comm_world
  end if

  call comm_%mpisum(self%cnt)

  call comm_%mpisum(self%mean)
  call comm_%mpisum(self%std)
  self%mean = self%mean / self%cnt
  self%std = sqrt(self%std/self%cnt - self%mean**2)

  call comm_%gather(self%max_, temp)
  call comm_%gather(self%max_index, temp_ind)
  pos = maxloc(temp, 1)
  self%max_ = temp(pos)
  self%max_index = temp_ind(pos)

  call comm_%gather(self%min_, temp)
  call comm_%gather(self%min_index, temp_ind)
  pos = minloc(temp, 1)
  self%min_ = temp(pos)
  self%min_index = temp_ind(pos)
end subroutine finalize_stat_var


!******************************************************************************** 
!
!       Print stat_var to the fh file
!
!******************************************************************************** 
subroutine print_stat_var(self, name, fh, comm)
  class(stat_var), intent(inout)               :: self
  character(len=*), intent(in)                 :: name
  type(FileHandle), intent(in)                 :: fh
  type(mpi_communicator), optional, intent(in) :: comm 

  write(fh%funit, *) 
  write(fh%funit, *)   "variable name: " // trim(name)
  write(fh%funit, *)   "  counter      : ", self%cnt
  write(fh%funit, *) "  mean value   : ", self%mean
  write(fh%funit, *) "  stddev       : ", self%std
  write(fh%funit, *)   "  min value    : ", self%min_
  write(fh%funit, *)   "  min index    : ", self%min_index
  write(fh%funit, *)   "  max value    : ", self%max_
  write(fh%funit, *)   "  max index    : ", self%max_index

  100 format (1x,a,f12.6)
  101 format (1x,a,i12)
end subroutine print_stat_var


!******************************************************************************** 
!
!       Initialize stat_var_c
!
!******************************************************************************** 
subroutine init_stat_var_c(self)
  class(stat_var_c), intent(inout) :: self
  !local variable
  real(wp)                         :: x

  self%cnt = 0.0_wp
  self%max_ = 0.0_wp
  self%min_ = 0.0_wp
  self%mean = 0.0_wp
  self%std = 0.0_wp
  self%min_index = 0
  self%max_index = 0
end subroutine init_stat_var_c


!******************************************************************************** 
!
!       Update stat_var_c
!
!******************************************************************************** 
subroutine update_stat_var_c(self, x, step)
  class(stat_var_c), intent(inout) :: self
  complex(wp), intent(in)          :: x
  integer, optional, intent(in)    :: step

  if (self%cnt < 0.5_wp) then
    self%cnt = 1.0_wp
    self%max_ = x
    self%max_index = 1
    self%min_ = x 
    self%min_index = 1
    self%mean = x
    self%std = abs(x)**2
    return
  end if

  self%cnt = self%cnt + 1.0_wp

  if (real(x, wp) > real(self%max_, wp)) then
    self%max_ = x
    if (present(step)) then
      self%max_index = step
    else
      self%max_index = self%cnt
    end if
  end if

  if (real(x, wp) < real(self%min_, wp)) then
    self%min_ = x
    if (present(step)) then
      self%min_index = step
    else
      self%min_index = self%cnt
    end if
  end if

  self%mean = self%mean + x
  self%std = self%std + abs(x)**2
end subroutine update_stat_var_c

subroutine update_stat_var_1d_c(self, x, step)
  class(stat_var_c), intent(inout) :: self
  complex(wp), intent(in)          :: x(:)
  integer, optional, intent(in)    :: step
  !local variables
  real(wp)                         :: xmax, xmin
  integer                          :: xmaxloc, xminloc, pos

  if (self%cnt < 0.5_wp) then
    self%cnt = real(size(x), wp)
    pos = maxloc(real(x, wp), 1)
    self%max_ = x(pos)
    if (present(step)) then
      self%max_index = step
    else
      self%max_index = pos
    end if
    pos = minloc(real(x, wp), 1)
    self%min_ = x(pos)
    if (present(step)) then
      self%min_index = step
    else
      self%min_index = pos
    end if
    self%mean = sum(x)
    self%std = sum(abs(x)**2)
    return
  end if

  xmax = maxval(real(x, wp))
  xmaxloc = maxloc(real(x, wp), 1)

  if (xmax > real(self%max_,wp)) then
    self%max_ = x(xmaxloc)
    if (present(step)) then
      self%max_index = step
    else
      self%max_index = self%cnt + xmaxloc
    end if
  end if

  xmin = minval(real(x, wp))
  xminloc = minloc(real(x, wp), 1)

  if (xmin < real(self%min_,wp)) then
    self%min_ = x(xminloc)
    if (present(step)) then
      self%min_index = step
    else
      self%min_index = self%cnt + xminloc
    end if
  end if

  self%cnt = self%cnt + real(size(x), wp)
  self%mean = self%mean + sum(x)
  self%std = self%std + sum(abs(x)**2)
end subroutine update_stat_var_1d_c


!******************************************************************************** 
!
!       Finalize stat_var
!
!******************************************************************************** 
subroutine finalize_stat_var_c(self, comm)
  class(stat_var_c), intent(inout)             :: self
  type(mpi_communicator), optional, intent(in) :: comm 
  !local variables
  type(mpi_communicator)                       :: comm_
  complex(wp), allocatable                     :: temp(:)
  integer, allocatable                         :: temp_ind(:)
  integer                                      :: pos

  if (present(comm)) then
    comm_ = comm
  else
    comm_ = comm_world
  end if

  call comm_%mpisum(self%cnt)

  call comm_%mpisum(self%mean)
  call comm_%mpisum(self%std)
  self%mean = self%mean / self%cnt
  self%std = sqrt(self%std/self%cnt - abs(self%mean)**2)

  call comm_%gather(self%max_, temp)
  call comm_%gather(self%max_index, temp_ind)
  pos = maxloc(real(temp, wp), 1)
  self%max_ = temp(pos)
  self%max_index = temp_ind(pos)

  call comm_%gather(self%min_, temp)
  call comm_%gather(self%min_index, temp_ind)
  pos = minloc(real(temp, wp), 1)
  self%min_ = temp(pos)
  self%min_index = temp_ind(pos)
end subroutine finalize_stat_var_c


!******************************************************************************** 
!
!       Print stat_var_c to the fh file
!
!******************************************************************************** 
subroutine print_stat_var_c(self, name, fh, comm)
  class(stat_var_c), intent(inout)             :: self
  character(len=*), intent(in)                 :: name
  type(FileHandle), intent(in)                 :: fh
  type(mpi_communicator), optional, intent(in) :: comm 

  write(fh%funit, *) 
  write(fh%funit, *)   "variable name: " // trim(name)
  write(fh%funit, *)   "  counter      : ", self%cnt
  write(fh%funit, *) "  mean value   : ", self%mean
  write(fh%funit, *) "  stddev       : ", self%std
  write(fh%funit, *)   "  min value    : ", self%min_
  write(fh%funit, *)   "  min index    : ", self%min_index
  write(fh%funit, *)   "  max value    : ", self%max_
  write(fh%funit, *)   "  max index    : ", self%max_index

  100 format (1x,a,2f12.6)
  101 format (1x,a,i12)
end subroutine print_stat_var_c

end module stat_analysis
