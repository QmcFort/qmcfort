! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module qmcfort_io

use iso_fortran_env, only: output_unit

use constants
use file_handle, only: FileHandle
use mpi

implicit none

interface write_perc
  module procedure write_perc_sp, write_perc_dp
end interface write_perc

type io_var
  type(FileHandle) :: screen
  type(FileHandle) :: qmcfort_log
  type(FileHandle) :: qmcfort_in
  type(FileHandle) :: qmcfort_out
  type(FileHandle) :: basis_set
  type(FileHandle) :: timing
end type io_var

type(io_var) :: io

contains

!********************************************************************************
!
!         Touch OUTPUT file and write first words into it
!
!********************************************************************************
subroutine init_io
  io%screen = FileHandle("screen", funit=output_unit)
  io%qmcfort_log = FileHandle("qmcfort_log")
  io%qmcfort_in = FileHandle("qmcfort_in")
  io%qmcfort_out = FileHandle("qmcfort_out")
  io%timing = FileHandle("timing")
  
  if (comm_world%mpirank == 0) then
    call io%qmcfort_log%open(status="replace", action="write")
    call io%qmcfort_out%open(status="replace", action="write")
    call dump_header()
  end if
  
  call comm_world%print_info(io%screen%funit)
  call comm_world%print_info(io%qmcfort_out%funit)
  call comm_world%print_info(io%qmcfort_log%funit)
end subroutine init_io


subroutine finalize_io
  if (comm_world%mpirank == 0) then
    call dump_footer()
    call io%qmcfort_log%close()
    call io%qmcfort_out%close()
  end if
end subroutine finalize_io


!******************************************************************************** 
!
! Write headers to the output files
!
!******************************************************************************** 
subroutine dump_header
  character(len=charlen) :: date, date_
  character(len=charlen) :: time, time_

  call date_and_time(date=date_, time=time_)
  date = get_month(date_(5:6)) // " " // date_(7:8) // " " // date_(1:4)
  time = time_(1:2) // ":" // time_(3:4) // ":" // time_(5:6)

  write(io%screen%funit,100)"***************************************************************************"
  write(io%screen%funit,100)"***                                                                     ***"
  write(io%screen%funit,100)"***                            Q m c F o r t                            ***"
  write(io%screen%funit,100)"***                   Quantum Monte Carlo with Fortran                  ***"
  write(io%screen%funit,100)"***                                                                     ***" 
  write(io%screen%funit,101)"***   started at: ", trim(date), "  ", trim(time),                     "***" 
  write(io%screen%funit,101)"***   built at  : ", __DATE__, "  ", __TIME__,                         "***" 
  write(io%screen%funit,102)"***   version   : ", VERSION,                                          "***" 
  write(io%screen%funit,100)"***                                                                     ***" 
  write(io%screen%funit,100)"***************************************************************************"
  
  write(io%qmcfort_log%funit,100)"***************************************************************************"
  write(io%qmcfort_log%funit,100)"***                                                                     ***"
  write(io%qmcfort_log%funit,100)"***                            Q m c F o r t                            ***"
  write(io%qmcfort_log%funit,100)"***                   Quantum Monte Carlo with Fortran                  ***"
  write(io%qmcfort_log%funit,100)"***                                                                     ***" 
  write(io%qmcfort_log%funit,101)"***   started at: ", trim(date), "  ", trim(time),                     "***" 
  write(io%qmcfort_log%funit,101)"***   built at  : ", __DATE__, "  ", __TIME__,                         "***" 
  write(io%qmcfort_log%funit,102)"***   version   : ", VERSION,                                          "***" 
  write(io%qmcfort_log%funit,100)"***                                                                     ***" 
  write(io%qmcfort_log%funit,100)"***************************************************************************"

  write(io%qmcfort_out%funit,100)"***************************************************************************"
  write(io%qmcfort_out%funit,100)"***                                                                     ***"
  write(io%qmcfort_out%funit,100)"***                            Q m c F o r t                            ***"
  write(io%qmcfort_out%funit,100)"***                   Quantum Monte Carlo with Fortran                  ***"
  write(io%qmcfort_out%funit,100)"***                                                                     ***" 
  write(io%qmcfort_out%funit,101)"***   started at: ", trim(date), "  ", trim(time),                     "***" 
  write(io%qmcfort_out%funit,101)"***   built at  : ", __DATE__, "  ", __TIME__,                         "***" 
  write(io%qmcfort_out%funit,102)"***   version   : ", VERSION,                                          "***" 
  write(io%qmcfort_out%funit,100)"***                                                                     ***" 
  write(io%qmcfort_out%funit,100)"***************************************************************************"

  100 format(1x,a)
  101 format(1x,a,a,a,a,t74,a)
  102 format(1x,a,a,t74,a)
end subroutine dump_header


!******************************************************************************** 
!
! Write footers to the output files
!
!******************************************************************************** 
subroutine dump_footer
  character(len=charlen) :: date, date_
  character(len=charlen) :: time, time_

  call date_and_time(date=date_, time=time_)
  date = get_month(date_(5:6)) // " " // date_(7:8) // " " // date_(1:4)
  time = time_(1:2) // ":" // time_(3:4) // ":" // time_(5:6)

  write(io%screen%funit,*)
  write(io%screen%funit,101)"     finished at: ", trim(date), "  ", trim(time)
  write(io%screen%funit,*)
  write(io%screen%funit,100)"***************************************************************************"
  
  write(io%qmcfort_log%funit,*)
  write(io%qmcfort_log%funit,101)"     finished at: ", trim(date), "  ", trim(time)
  write(io%qmcfort_log%funit,*)
  write(io%qmcfort_log%funit,100)"***************************************************************************"

  write(io%qmcfort_out%funit,*)
  write(io%qmcfort_out%funit,101)"     finished at: ", trim(date), "  ", trim(time)
  write(io%qmcfort_out%funit,*)
  write(io%qmcfort_out%funit,100)"***************************************************************************"

  100 format(1x,a)
  101 format(1x,a,a,a,a)
end subroutine dump_footer


function get_month(mm) result(mm_)
  character(len=*), intent(in)  :: mm
  character(len=:), allocatable :: mm_

  select case (trim(mm))
    case ("1", "01")
      mm_ = "Jan"
    case ("2", "02")
      mm_ = "Feb"
    case ("3", "03")
      mm_ = "Mar"
    case ("4", "04")
      mm_ = "Apr"
    case ("5", "05")
      mm_ = "May"
    case ("6", "06")
      mm_ = "Jun"
    case ("7", "07")
      mm_ = "Jul"
    case ("8", "08")
      mm_ = "Aug"
    case ("9", "09")
      mm_ = "Sep"
    case ("10")
      mm_ = "Oct"
    case ("11")
      mm_ = "Nov"
    case ("12")
      mm_ = "Dec"
    case default
      mm_ = "XXX"
  end select
end function get_month


!******************************************************************************** 
!
! IO function for writing structured timings:
!
!    X s --> H hours M minutes S seconds       
!
!******************************************************************************** 
function write_time(time) result(io_time)
  real(dp), intent(in)    :: time
  character(len=charlen)  :: io_time
  !local variables
  integer                 :: hours, mins, secs
  integer                 :: leng
  real(dp)                :: rest

  !determine how many hours
  hours = floor(time / 3600)
  rest = time - hours * 3600

  !determine how many minutes
  mins = floor(rest / 60)
  rest = rest - mins * 60

  !remained are seconds
  secs = floor(rest)
  
  leng = -1

  if (hours==0 .and. mins==0 .and. secs==0) then
    write(io_time,1003) rest
    1003 format(f8.6)
    io_time = trim(io_time) // "s"
  else if (hours==0 .and. mins==0) then
    write(io_time,1004) rest
    1004 format(f6.2)
    io_time = trim(io_time) // "s"
  else
    if (hours > 0) then
      write(io_time,"(i0)") hours
      io_time = trim(io_time) //  "h "
      leng = len(trim(io_time))
    end if
    
    if (mins > 0) then
      write(io_time(leng+2:),"(i0)") mins
      io_time = trim(io_time) // "min "
      leng = len(trim(io_time))
    end if
    
    write(io_time(leng+2:),"(i0)") secs
    io_time = trim(io_time) // "s"
  end if
end function write_time


!******************************************************************************** 
!
! IO function for writing number of FLOPS with appropriate unit
!
!******************************************************************************** 
function write_flops(flop, time) result(io_flops)
  real(dp), intent(in)   :: flop
  real(dp), intent(in)   :: time
  character(len=charlen) :: io_flops
  !local variables 
  real(dp)               :: flops
  integer                :: expon
 
  if (flop < 0.1_dp) then
    io_flops = ""
  else
    flops = flop / time
    expon = floor(log10(flops))
    expon = expon - mod(expon, 3)
    flops = flops / 10**(expon)
    write(io_flops,"(f6.1)") flops
    io_flops = trim(io_flops) // get_flops_prefix(expon)
  end if
contains
  function get_flops_prefix(expon) result(flops_unit)
    integer, intent(in)           :: expon
    character(len=:), allocatable :: flops_unit
    
    select case (expon)
      case (0)
        flops_unit = "FLOPS"
      case (3)
        flops_unit = "KFLOPS"
      case (6)
        flops_unit = "MFLOPS"
      case (9) 
        flops_unit = "GFLOPS"
      case (12)
        flops_unit = "TFLOPS"
      case (15)
        flops_unit = "PFLOPS"
      case default
        flops_unit = "XFLOPS"
    end select
  end function get_flops_prefix
end function write_flops


!******************************************************************************** 
!
! IO function for writing FLOP number with appropriate unit
!
!******************************************************************************** 
function write_flop(flop) result(io_flop)
  real(dp), intent(in)   :: flop
  character(len=charlen) :: io_flop
  !local variables
  real(dp)               :: flop_
  integer                :: expon
 
  if (flop < 1.0_dp) then
    io_flop = "0FLOP"
  else
    expon = floor(log10(flop))
    expon = expon - mod(expon, 3)
    flop_ = flop / (10.0_dp)**(expon)
    write(io_flop,"(f6.1)") flop_
    io_flop = trim(io_flop) // get_flop_prefix(expon)
  end if
contains
  function get_flop_prefix(expon) result(flop_unit)
    integer, intent(in)           :: expon
    character(len=:), allocatable :: flop_unit
    
    select case (expon)
      case (0)
        flop_unit = "FLOP"
      case (3)
        flop_unit = "kFLOP"
      case (6)
        flop_unit = "MFLOP"
      case (9) 
        flop_unit = "GFLOP"
      case (12)
        flop_unit = "TFLOP"
      case (15)
        flop_unit = "PFLOP"
      case (18)
        flop_unit = "EFLOP"
      case default
        flop_unit = "xFLOP"
    end select
  end function get_flop_prefix
end function write_flop


!******************************************************************************** 
!
! IO function for writing number of bytes with appropriate unit
!
!******************************************************************************** 
function write_bytes(bytes) result(io_bytes)
  real(dp), intent(in)    :: bytes
  character(len=charlen)  :: io_bytes
  !local variables
  real(dp)                :: bytes_
  integer                 :: expon
 
  expon = floor(log10(bytes))
  expon = expon - mod(expon, 3)
  bytes_ = bytes / 10**(expon)
  write(io_bytes,"(f6.1)") bytes_
  io_bytes = trim(io_bytes) // get_bytes_prefix(expon)
contains
  function get_bytes_prefix(expon) result(bytes_unit)
    integer, intent(in)           :: expon
    character(len=:), allocatable :: bytes_unit
    
    select case (expon)
      case (0)
        bytes_unit = "B"
      case (3)
        bytes_unit = "KB"
      case (6)
        bytes_unit = "MB"
      case (9) 
        bytes_unit = "GB"
      case (12)
        bytes_unit = "TB"
      case (15)
        bytes_unit = "PB"
      case default
        bytes_unit = "XB"
    end select
  end function get_bytes_prefix
end function write_bytes


!******************************************************************************** 
!
! IO function for writing structured percentages
!
!******************************************************************************** 
function write_perc_sp(x) result(io_perc)
  real(sp), intent(in)   :: x
  character(len=charlen) :: io_perc
  !local
  real(sp)               :: p
  
  p = 100.0_sp * x
  write(io_perc,100) p
  io_perc = trim(io_perc) // "%"
  100 format(f6.2)
end function write_perc_sp

function write_perc_dp(x) result(io_perc)
  real(dp), intent(in)   :: x
  character(len=charlen) :: io_perc
  !local
  real(dp)               :: p
  
  p = 100.0_dp * x
  write(io_perc,100) p
  io_perc = trim(io_perc) // "%"
  100 format(f6.2)
end function write_perc_dp


!******************************************************************************** 
!
! IO function for writing structured percentages
!
!******************************************************************************** 
function write_prom(x) result(io_prom)
  real(wp), intent(in)          :: x
  character(len=:), allocatable :: prom_str
  !local
  character(len=charlen)        :: io_prom
  real(wp)                      :: p
  
  p = 1000.0_wp * x
  write(io_prom,100) p
  io_prom = trim(io_prom) // "%%"
  prom_str = trim(io_prom)

  100 format(f7.2)
end function write_prom

end module qmcfort_io
