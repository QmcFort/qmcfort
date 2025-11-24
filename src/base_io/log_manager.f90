! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module log_manager_mod

use constants
use mpi
use file_handle, only: FileHandle

implicit none

private
public :: LogManager, logging

enum, bind(c)
  enumerator :: STOP_MESSAGE=1, ALLOC_ERROR 
end enum

type LogManager
  type(FileHandle), allocatable :: fhs(:)
contains
  generic, public    :: info => log_info_msg, log_info_msg_id
  generic, public    :: warning => log_warning_msg, log_warning_msg_id
  generic, public    :: error => log_error_msg, log_error_msg_id
  generic, public    :: bug => log_bug_msg, log_bug_msg_id

  procedure, private :: get_msg => get_typed_msg
  procedure, private :: send_msg_to_io_rank
  procedure, private :: write_msg_to_files
  procedure, private :: log_info_msg, log_info_msg_id
  procedure, private :: log_warning_msg, log_warning_msg_id
  procedure, private :: log_error_msg, log_error_msg_id
  procedure, private :: log_bug_msg, log_bug_msg_id
end type LogManager

interface LogManager
  procedure finit_log_manager
end interface LogManager

type LogFormatter
  integer :: text_width = 80
  integer :: pos_msg = 15
contains
  procedure          :: format_msg
  procedure          :: generate_header
  procedure          :: generate_footer
  procedure, private :: get_next_line
end type LogFormatter

type(LogManager) :: logging

contains

!******************************************************************************** 
!
! Initialization of the LogManager object
!
!******************************************************************************** 
subroutine init_log_manager(fhs, self)
  type(FileHandle), intent(in)  :: fhs(:)
  type(LogManager), intent(out) :: self
  !local
  integer                       :: i

  allocate(self%fhs(size(fhs)))

  do i = 1, size(fhs)
    self%fhs(i) = fhs(i)
  end do
end subroutine init_log_manager


!******************************************************************************** 
!
! LogManager constructor
!
!******************************************************************************** 
function finit_log_manager(fhs) result(self)
  type(FileHandle), intent(in) :: fhs(:)
  type(LogManager)             :: self

  call init_log_manager(fhs, self)
end function finit_log_manager


!******************************************************************************** 
!
! Format the raw INFO message and log it to the files
!
!******************************************************************************** 
subroutine log_info_msg(self, msg)
  class(LogManager), intent(in) :: self
  character(len=*), intent(in)  :: msg
  !local
  character(len=*), parameter   :: tag = "LOG: INFO: "
  character(len=:), allocatable :: formatted_msg
  type(LogFormatter)            :: formatter

  formatted_msg = formatter%format_msg(tag, msg)
  if (comm_world%mpirank == 0) call self%write_msg_to_files(formatted_msg)
end subroutine log_info_msg


!******************************************************************************** 
!
! Extract the raw INFO message, and call the info
!
!******************************************************************************** 
subroutine log_info_msg_id(self, msg_id)
  class(LogManager), intent(in) :: self
  integer, intent(in)           :: msg_id
  !local
  character(len=:), allocatable :: msg

  msg = self%get_msg(msg_id)
  call self%info(msg)
end subroutine log_info_msg_id


!******************************************************************************** 
!
! Format the raw WARNING message and log it to the files
!
!******************************************************************************** 
subroutine log_warning_msg(self, msg)
  class(LogManager), intent(in) :: self
  character(len=*), intent(in)  :: msg
  !local
  character(len=*), parameter   :: tag = "LOG: WARNING: "
  character(len=:), allocatable :: formatted_msg
  type(LogFormatter)            :: formatter

  formatted_msg = formatter%format_msg(tag, msg)
  if (comm_world%mpirank == 0) call self%write_msg_to_files(formatted_msg)
end subroutine log_warning_msg


!******************************************************************************** 
!
! Extract the raw WARNING message from msg_id and issue the warning
!
!******************************************************************************** 
subroutine log_warning_msg_id(self, msg_id)
  class(LogManager), intent(in) :: self
  integer, intent(in)           :: msg_id
  !local
  character(len=:), allocatable :: msg

  msg = self%get_msg(msg_id)
  call self%warning(msg)
end subroutine log_warning_msg_id


!******************************************************************************** 
!
! Format the raw ERROR message, log it to the files and stop the program
!
!******************************************************************************** 
subroutine log_error_msg(self, is_trigger, msg, file, line)
  class(LogManager), intent(in) :: self
  logical, intent(in)           :: is_trigger
  character(len=*), intent(in)  :: msg
  character(len=*), intent(in)  :: file
  integer, intent(in)           :: line
  !local
  character(len=*), parameter   :: tag = "LOG: ERROR: "
  character(len=:), allocatable :: header, formatted_msg, footer, full_msg, stop_msg
  type(LogFormatter)            :: formatter
  
  stop_msg = self%get_msg(STOP_MESSAGE)
  header = formatter%generate_header(tag, file, line, comm_world%mpirank)
  formatted_msg = formatter%format_msg(tag, msg)
  footer = formatter%generate_footer(tag, stop_msg)
  full_msg = header // formatted_msg // footer

  call self%send_msg_to_io_rank(is_trigger, full_msg)

  if (comm_world%mpirank == 0) call self%write_msg_to_files(full_msg)
  call stop_program()
end subroutine log_error_msg


!******************************************************************************** 
!
! Extract the raw ERROR message from msg_id, and raise an error
!
!******************************************************************************** 
subroutine log_error_msg_id(self, is_trigger, msg_id, file, line)
  class(LogManager), intent(in) :: self
  logical, intent(in)           :: is_trigger
  integer, intent(in)           :: msg_id
  character(len=*), intent(in)  :: file
  integer, intent(in)           :: line
  !local
  character(len=:), allocatable :: msg

  msg = self%get_msg(msg_id)
  call self%error(is_trigger, msg, file, line)
end subroutine log_error_msg_id


!******************************************************************************** 
!
! Format the raw BUG message, log it to the files and stop the program
!
!******************************************************************************** 
subroutine log_bug_msg(self, is_trigger, msg, file, line)
  class(LogManager), intent(in) :: self
  logical, intent(in)           :: is_trigger
  character(len=*), intent(in)  :: msg
  character(len=*), intent(in)  :: file
  integer, intent(in)           :: line
  !local
  character(len=*), parameter   :: tag = "LOG: BUG: "
  character(len=:), allocatable :: header, formatted_msg, footer, full_msg, stop_msg
  type(LogFormatter)            :: formatter
  
  stop_msg = self%get_msg(STOP_MESSAGE)
  header = formatter%generate_header(tag, file, line, comm_world%mpirank)
  formatted_msg = formatter%format_msg(tag, msg)
  footer = formatter%generate_footer(tag, stop_msg)
  full_msg = header // formatted_msg // footer

  call self%send_msg_to_io_rank(is_trigger, full_msg)

  if (comm_world%mpirank == 0) call self%write_msg_to_files(full_msg)
  call stop_program()
end subroutine log_bug_msg


!******************************************************************************** 
! 
! Extract the raw BUG message from msg_id, and raise a BUG
!
!******************************************************************************** 
subroutine log_bug_msg_id(self, is_trigger, msg_id, file, line)
  class(LogManager), intent(in) :: self
  logical, intent(in)           :: is_trigger
  integer, intent(in)           :: msg_id
  character(len=*), intent(in)  :: file
  integer, intent(in)           :: line
  !local
  character(len=:), allocatable :: msg

  msg = self%get_msg(msg_id)
  call self%bug(is_trigger, msg, file, line)
end subroutine log_bug_msg_id


!******************************************************************************** 
!
! Determine the rank whose message should be written, and send the message
! to the IO rank
!
!******************************************************************************** 
subroutine send_msg_to_io_rank(self, is_trigger, msg)
  class(LogManager), intent(in)                :: self
  logical, intent(in)                          :: is_trigger
  character(len=:), allocatable, intent(inout) :: msg
  !local
  integer                                      :: trigger_rank, msg_len

  if (is_trigger) then
    trigger_rank = comm_world%mpirank
  else
    trigger_rank = comm_world%mpisize
  end if

  call comm_world%reduce(trigger_rank, op=MPI_Min)

  if (comm_world%mpirank == trigger_rank) then
    call comm_world%send(len(msg), 0, tag=0)
    call comm_world%send(msg, 0, tag=0)
  else if (comm_world%mpirank == 0) then
    call comm_world%recv(msg_len, trigger_rank, tag=0)
    if (allocated(msg)) deallocate(msg)
    allocate(character(len=msg_len) :: msg)
    call comm_world%recv(msg, trigger_rank, tag=0)
  end if
end subroutine send_msg_to_io_rank


!******************************************************************************** 
!
! Write formatted message to log files
!
!******************************************************************************** 
subroutine write_msg_to_files(self, msg)
  class(LogManager), intent(in) :: self
  character(len=*), intent(in)  :: msg
  !local
  integer                       :: i

  do i = 1, size(self%fhs)
    write(self%fhs(i)%funit, "(a)") msg 
  end do
end subroutine write_msg_to_files


!******************************************************************************** 
!
! Retrieve the log message from the message id
!
!******************************************************************************** 
function get_typed_msg(self, msg_id) result(msg)
  class(LogManager), intent(in) :: self
  integer, intent(in)           :: msg_id
  character(len=:), allocatable :: msg

  select case (msg_id)
    case (ALLOC_ERROR)
      msg = "Allocation of the array failed!"
    case (STOP_MESSAGE)
      msg = "Program stops now!"
    case default
      msg = "Unknown message id requested!"
  end select
end function get_typed_msg


!******************************************************************************** 
!
! Generate a formatted message header for ERROR and BUG messages.
! The format is:
!    [Rank n] (file: fname.f90, line: m)
!
!******************************************************************************** 
function generate_header(self, tag, file, line, rank) result(header)
  class(LogFormatter), intent(in) :: self
  character(len=*), intent(in)    :: tag
  character(len=*), intent(in)    :: file
  integer, intent(in)             :: line
  integer, intent(in)             :: rank
  character(len=:), allocatable   :: header

  header = tag // "[Rank " // i_to_a(rank) // "] (file: " // file // ", line: " // i_to_a(line) // ")" // new_line("a") 
end function generate_header


!******************************************************************************** 
!
! Generate a formatted message footer for ERROR and BUG messages.
!
!******************************************************************************** 
function generate_footer(self, tag, msg) result(footer)
  class(LogFormatter), intent(in) :: self
  character(len=*), intent(in)    :: tag
  character(len=*), intent(in)    :: msg
  character(len=:), allocatable   :: footer
  !local
  integer                         :: indent

  indent = (self%text_width - len(msg)) / 2
  footer = tag // repeat(" ", indent) // msg // new_line("a")
end function generate_footer


!******************************************************************************** 
!
! Format a raw log message for writing to log files.
!
!******************************************************************************** 
function format_msg(self, tag, msg) result(formatted_msg)
  class(LogFormatter), intent(in) :: self
  character(len=*), intent(in)    :: tag
  character(len=*), intent(in)    :: msg
  character(len=:), allocatable   :: formatted_msg
  !local
  integer                         :: indent
  character(len=:), allocatable   :: next_line, residual_msg 

  indent = self%pos_msg - len(tag)

  formatted_msg = ""
  residual_msg = msg

  do 
    if (len_trim(residual_msg) == 0)  then
      formatted_msg = formatted_msg // tag // new_line("a")
      exit
    end if
    call self%get_next_line(residual_msg, next_line)
    formatted_msg = formatted_msg // tag // repeat(" ", indent) // next_line // new_line("a")
  end do
end function format_msg


!******************************************************************************** 
!
! Extract a text line of at most text_width characters from the residual message,
! and update the residual message by removing the extracted portion
!
!******************************************************************************** 
subroutine get_next_line(self, residual_msg, next_line)       
  class(LogFormatter), intent(in)              :: self
  character(len=:), allocatable, intent(inout) :: residual_msg
  character(len=:), allocatable, intent(out)   :: next_line
  !local
  integer                                      :: left, right, break_pos
                                                              
  if (len(residual_msg) <= self%text_width) then
    next_line = residual_msg
    residual_msg = ""
  else if (scan(residual_msg, " ") == 0) then
    !check for completeness — though this should already be covered by the first test.
    !it would be unusual for the residual_msg to exceed the text_width and still not contain a space.
    next_line = residual_msg
    residual_msg = ""
  else
    !decide where to break the line
    if (residual_msg(self%text_width:self%text_width) == " ") then
      break_pos = self%text_width
    else 
      !find the closest space to residual_msg(text_width)
      left = scan(residual_msg(1:self%text_width), " ", back=.true.)
      right = scan(residual_msg(self%text_width+1:), " ")

      if (self%text_width-left <= right) then
        break_pos = left
      else
        break_pos = self%text_width + right
      end if
    end if

    next_line = residual_msg(1:break_pos-1)
    residual_msg = residual_msg(break_pos+1:)
  end if
end subroutine get_next_line


!******************************************************************************** 
!
! Convert integer to the string. 
! Length of the string matches the number of digits in the integer
!
! DEBUG: should go to the string module on the long run
!
!******************************************************************************** 
function i_to_a(i) result(s)
  integer, intent(in)           :: i
  character(len=:), allocatable :: s
  !local
  character(len=charlen)        :: buffer

  write(buffer, '(I0)') i
  s = trim(buffer)
end function i_to_a


!******************************************************************************** 
!
! Terminate the program using MPI_Abort.
! If MPI_Abort fails, fall back to Fortran's stop statement.
!
!******************************************************************************** 
subroutine stop_program(error_code)
  integer, optional, intent(in) :: error_code
  !local
  integer                       :: error_code_

  call comm_world%barrier()

  error_code_ = 1
  if (present(error_code)) error_code_ = error_code

  call comm_world%abort(error_code_)

  !fallback in case MPI_Abort did not work
  stop error_code_
end subroutine stop_program

end module log_manager_mod    