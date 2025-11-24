! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module mpi

#include "preproc.inc"
#ifdef MPI
use mpi_f08
#endif
!$ use omp_lib

use constants

implicit none

public 

type mpi_communicator
#ifdef MPI
  type(MPI_Comm)     :: mpi_comm
#endif
  integer            :: mpirank
  integer            :: mpinode
  integer            :: mpisize
  
  logical            :: lomp = .false.
  logical            :: initialized = .false.
  integer            :: ompsize
  integer            :: blasthr
  integer            :: thread_levels
  
  integer            :: totsize
contains
  procedure, public  :: init       => init_mpi            
  procedure, public  :: pseudoinit => pseudoinit_mpi
  procedure, public  :: finalize   => finalize_mpi
  procedure, public  :: barrier    => barrier
  procedure, public  :: abort      => abort
  procedure, public  :: print_info => print_comm_info
  generic, public    :: send       => send_char, send_log, send_int, send_int_1d, send_int8, send_int8_1d, &
                                      send_r_sp, send_r1_sp, send_r2_sp, send_r3_sp, send_r4_sp, &
                                      send_r_dp, send_r1_dp, send_r2_dp, send_r3_dp, send_r4_dp, &
                                      send_c_sp, send_c1_sp, send_c2_sp, send_c3_sp, send_c4_sp, &
                                      send_c_dp, send_c1_dp, send_c2_dp, send_c3_dp, send_c4_dp 
  generic, public    :: recv       => recv_char, recv_log, recv_int, recv_int_1d, recv_int8, recv_int8_1d, &
                                      recv_r_sp, recv_r1_sp, recv_r2_sp, recv_r3_sp, recv_r4_sp, &
                                      recv_r_dp, recv_r1_dp, recv_r2_dp, recv_r3_dp, recv_r4_dp, &
                                      recv_c_sp, recv_c1_sp, recv_c2_sp, recv_c3_sp, recv_c4_sp, &
                                      recv_c_dp, recv_c1_dp, recv_c2_dp, recv_c3_dp, recv_c4_dp 
  generic, public    :: bcast      => bcast_char, bcast_char_1d, bcast_log, &
                                      bcast_int, bcast_int_1d, bcast_int_2d, &
                                      bcast_int8, bcast_int8_1d, &
                                      bcast_r_sp, bcast_r1_sp, bcast_r2_sp, bcast_r3_sp, bcast_r4_sp, bcast_r5_sp, &
                                      bcast_r_dp, bcast_r1_dp, bcast_r2_dp, bcast_r3_dp, bcast_r4_dp, bcast_r5_dp, &
                                      bcast_c_sp, bcast_c1_sp, bcast_c2_sp, bcast_c3_sp, bcast_c4_sp, bcast_c5_sp, &
                                      bcast_c_dp, bcast_c1_dp, bcast_c2_dp, bcast_c3_dp, bcast_c4_dp, bcast_c5_dp
  generic, public    :: gather     => gather_int, &
                                      gather_r_sp, gather_r1_sp, gather_r2_sp, &
                                      gather_r_dp, gather_r1_dp, gather_r2_dp, &
                                      gather_c_sp, gather_c1_sp, gather_c2_sp, &
                                      gather_c_dp, gather_c1_dp, gather_c2_dp
  generic, public    :: mpisum     => mpi_sum_r_sp, mpi_sum_r1_sp, mpi_sum_r2_sp, mpi_sum_c_sp, mpi_sum_c1_sp, mpi_sum_c2_sp, &
                                      mpi_sum_r_dp, mpi_sum_r1_dp, mpi_sum_r2_dp, mpi_sum_c_dp, mpi_sum_c1_dp, mpi_sum_c2_dp, &
                                      mpi_sum_i, mpi_sum_i1
  procedure, private :: send_char, send_log, send_int, send_int_1d, send_int8, send_int8_1d
  procedure, private :: send_r_sp, send_r1_sp, send_r2_sp, send_r3_sp, send_r4_sp
  procedure, private :: send_r_dp, send_r1_dp, send_r2_dp, send_r3_dp, send_r4_dp
  procedure, private :: send_c_sp, send_c1_sp, send_c2_sp, send_c3_sp, send_c4_sp
  procedure, private :: send_c_dp, send_c1_dp, send_c2_dp, send_c3_dp, send_c4_dp 
  procedure, private :: recv_char, recv_log, recv_int, recv_int_1d, recv_int8, recv_int8_1d
  procedure, private :: recv_r_sp, recv_r1_sp, recv_r2_sp, recv_r3_sp, recv_r4_sp
  procedure, private :: recv_r_dp, recv_r1_dp, recv_r2_dp, recv_r3_dp, recv_r4_dp
  procedure, private :: recv_c_sp, recv_c1_sp, recv_c2_sp, recv_c3_sp, recv_c4_sp
  procedure, private :: recv_c_dp, recv_c1_dp, recv_c2_dp, recv_c3_dp, recv_c4_dp 
  procedure, private :: bcast_char, bcast_char_1d, bcast_log, bcast_int, bcast_int8, bcast_int_1d, bcast_int8_1d, bcast_int_2d
  procedure, private :: bcast_r_sp, bcast_r1_sp, bcast_r2_sp, bcast_r3_sp, bcast_r4_sp, bcast_r5_sp
  procedure, private :: bcast_r_dp, bcast_r1_dp, bcast_r2_dp, bcast_r3_dp, bcast_r4_dp, bcast_r5_dp
  procedure, private :: bcast_c_sp, bcast_c1_sp, bcast_c2_sp, bcast_c3_sp, bcast_c4_sp, bcast_c5_sp
  procedure, private :: bcast_c_dp, bcast_c1_dp, bcast_c2_dp, bcast_c3_dp, bcast_c4_dp, bcast_c5_dp
  procedure, private :: gather_int 
  procedure, private :: gather_r_sp, gather_r1_sp, gather_r2_sp
  procedure, private :: gather_r_dp, gather_r1_dp, gather_r2_dp
  procedure, private :: gather_c_sp, gather_c1_sp, gather_c2_sp
  procedure, private :: gather_c_dp, gather_c1_dp, gather_c2_dp
  procedure, private :: mpi_sum_r_sp, mpi_sum_r1_sp, mpi_sum_r2_sp, mpi_sum_c_sp, mpi_sum_c1_sp, mpi_sum_c2_sp
  procedure, private :: mpi_sum_r_dp, mpi_sum_r1_dp, mpi_sum_r2_dp, mpi_sum_c_dp, mpi_sum_c1_dp, mpi_sum_c2_dp
  procedure, private :: mpi_sum_i, mpi_sum_i1
#ifdef MPI
  generic, public    :: reduce  => reduce_i, reduce_log, &
                                   reduce_r_sp, reduce_r1_sp, reduce_r2_sp, reduce_c_sp, reduce_c1_sp, reduce_c2_sp, &
                                   reduce_r_dp, reduce_r1_dp, reduce_r2_dp, reduce_c_dp, reduce_c1_dp, reduce_c2_dp
  procedure, private :: reduce_i, reduce_log
  procedure, private :: reduce_r_sp, reduce_r1_sp, reduce_r2_sp, reduce_c_sp, reduce_c1_sp, reduce_c2_sp
  procedure, private :: reduce_r_dp, reduce_r1_dp, reduce_r2_dp, reduce_c_dp, reduce_c1_dp, reduce_c2_dp
#endif
end type mpi_communicator

type(mpi_communicator) :: comm_world

!debug:
#ifdef MKL
interface 
integer function mkl_get_max_threads()
end function mkl_get_max_threads
end interface
#endif

contains

subroutine init_mpi(comm)
  class(mpi_communicator), intent(out) :: comm
  !local variables
  integer                              :: ierror, request
  
  comm%thread_levels = 2

#ifdef MPI
  !$ comm%lomp = .true.
  if (comm%lomp) then
    !call MPI_Init_Thread(MPI_Thread_Funneled, request, ierror)
    !if (MPI_Thread_Funneled /= request) then
    call MPI_Init_Thread(MPI_Thread_Serialized, request, ierror)
    if (MPI_Thread_Serialized /= request) then
      write(*,*) "MPI thread initialization not satisfied - Program stops now!"
      stop
    end if
    call omp_set_max_active_levels(comm%thread_levels)
  else 
    call MPI_Init(ierror)
    comm%ompsize = 1
  end if
  
  comm%mpi_comm = MPI_Comm_World
  call MPI_Comm_Rank(comm%mpi_comm, comm%mpirank, ierror)
  call MPI_Comm_Size(comm%mpi_comm, comm%mpisize, ierror)
  comm%mpinode = comm%mpirank + 1
  
  comm%blasthr = 1

  if (comm%lomp) then
    !$omp parallel
#ifdef _OPENMP
    comm%ompsize = omp_get_num_threads()
#ifdef MKL
    comm%blasthr = mkl_get_max_threads()
#endif
#endif
    !$omp end parallel
  else
    comm%ompsize = 1
#ifdef MKL
    comm%blasthr = mkl_get_max_threads()
#endif
  end if
  !
  if (ierror /= 0) then
    write(*,*) "MPI Communicator can not be initialized - Program stops now!" 
    stop
  end if
#else
  comm%mpirank = 0
  comm%mpisize = 1
  comm%mpinode = 1
  !
  if (comm%lomp) then
    !$omp parallel
#ifdef _OPENMP
    comm%ompsize = omp_get_num_threads()
#endif
    !$omp end parallel
  else
    comm%ompsize = 1
  end if
#endif
  comm%totsize = comm%mpisize * comm%ompsize 
  comm%initialized = .true.
end subroutine init_mpi


!******************************************************************************** 
!
!       Initializes mpi_communicator variable without calling any MPI routine
!
!       Useful for unit tests
!
!******************************************************************************** 
subroutine pseudoinit_mpi(comm)
  class(mpi_communicator), intent(out) :: comm
  !
  !$ comm%lomp = .true.
  comm%mpirank = 0
  comm%mpisize = 1
  comm%mpinode = 1
  !
  if (comm%lomp) then
    !$omp parallel
#ifdef _OPENMP
    comm%ompsize = omp_get_num_threads()
#endif
    !$omp end parallel
  else
    comm%ompsize = 1
  end if
  comm%totsize = comm%mpisize * comm%ompsize 
  comm%initialized = .true.
end subroutine pseudoinit_mpi


!******************************************************************************** 
!
! Routine to finalize MPI module
!
!******************************************************************************** 
subroutine finalize_mpi(comm)
  class(mpi_communicator), intent(in) :: comm
  !local variables
  integer                             :: ierror
#ifdef MPI
  call MPI_Finalize(ierror)
#endif
end subroutine finalize_mpi


!******************************************************************************** 
!
! Print basic info about the MPI communicator including OpenMP & BLAS
!
!******************************************************************************** 
subroutine print_comm_info(comm, funit)
  class(mpi_communicator), intent(in) :: comm
  integer, intent(in)                 :: funit

  if (comm%mpirank == 0) then
    write(funit,100) "MPI Communicator information:"
    write(funit,100) "-----------------------------" 
    write(funit,101) "  Number of MPI processes      = ", comm%mpisize
    write(funit,*)
    write(funit,102) "  OpenMP allowed?              = ", comm%lomp
    write(funit,101) "  Allowed no. of OpenMP levels = ", comm%thread_levels
    write(funit,101) "  Number of OpenMP threads     = ", comm%ompsize
    write(funit,101) "  Number of BLAS threads       = ", comm%blasthr
    write(funit,*)
    write(funit,101) "  Total number of comp. units  = ", comm%mpisize * comm%ompsize * comm%blasthr
    write(funit,*) starlong
  end if

  100 format(1x,a)
  101 format(1x,a,i6)
  102 format(1x,a,l4)
end subroutine print_comm_info


!******************************************************************************** 
!
! Interface routine for MPI_Barrier
!
!******************************************************************************** 
subroutine barrier(comm, timeout)
  class(mpi_communicator), intent(in) :: comm
  real(dp), optional, intent(in)      :: timeout
  !local variables
  integer                             :: ierror
  real(dp)                            :: tstart, tend
  !
#ifdef MPI
  if (present(timeout)) then
    tstart = MPI_Wtime()
    tend = MPI_Wtime()
    do while (tend-tstart < timeout)
      tend = MPI_Wtime()
    end do
  end if
  call MPI_Barrier(comm%mpi_comm, ierror)
#endif
end subroutine barrier


!******************************************************************************** 
!
! Interface to MPI_Abort
!
!******************************************************************************** 
subroutine abort(comm, error_code)
  class(mpi_communicator), intent(in) :: comm
  integer, optional, intent(in)       :: error_code
  !local
  integer                             :: ierror, error_code_

  error_code_ = 1
  if (present(error_code)) error_code_ = error_code
  
#ifdef MPI
  call MPI_Abort(comm%mpi_comm, error_code_, ierror)
#endif
end subroutine abort


!******************************************************************************** 
!
! Interface routines for MPI_Send
!
!******************************************************************************** 
subroutine send_char(comm, buffer, dest, tag)
  class(mpi_communicator), intent(in) :: comm
  character(len=*), intent(in)        :: buffer
  integer, intent(in)                 :: dest
  integer, intent(in)                 :: tag

#ifdef MPI
  call MPI_Send(buffer, len(buffer), MPI_Character, dest, tag, comm%mpi_comm)
#endif
end subroutine send_char

subroutine send_log(comm, buffer, dest, tag)
  class(mpi_communicator), intent(in) :: comm
  logical, intent(in)                 :: buffer
  integer, intent(in)                 :: dest
  integer, intent(in)                 :: tag

#ifdef MPI
  call MPI_Send(buffer, 1, MPI_Logical, dest, tag, comm%mpi_comm)
#endif
end subroutine send_log

subroutine send_int(comm, buffer, dest, tag)
  class(mpi_communicator), intent(in) :: comm
  integer, intent(in)                 :: buffer
  integer, intent(in)                 :: dest
  integer, intent(in)                 :: tag

#ifdef MPI
  call MPI_Send(buffer, 1, MPI_Integer, dest, tag, comm%mpi_comm)
#endif
end subroutine send_int

subroutine send_int_1d(comm, buffer, dest, tag)
  class(mpi_communicator), intent(in) :: comm
  integer, intent(in)                 :: buffer(:)
  integer, intent(in)                 :: dest
  integer, intent(in)                 :: tag

#ifdef MPI
  call MPI_Send(buffer, size(buffer), MPI_Integer, dest, tag, comm%mpi_comm)
#endif
end subroutine send_int_1d

subroutine send_int8(comm, buffer, dest, tag)
  class(mpi_communicator), intent(in) :: comm
  integer(i8), intent(in)             :: buffer
  integer, intent(in)                 :: dest
  integer, intent(in)                 :: tag

#ifdef MPI
  call MPI_Send(buffer, 1, MPI_Integer8, dest, tag, comm%mpi_comm)
#endif
end subroutine send_int8

subroutine send_int8_1d(comm, buffer, dest, tag)
  class(mpi_communicator), intent(in) :: comm
  integer(i8), intent(in)             :: buffer(:)
  integer, intent(in)                 :: dest
  integer, intent(in)                 :: tag

#ifdef MPI
  call MPI_Send(buffer, size(buffer), MPI_Integer8, dest, tag, comm%mpi_comm)
#endif
end subroutine send_int8_1d

subroutine send_r_sp(comm, buffer, dest, tag)
  class(mpi_communicator), intent(in) :: comm
  real(sp), intent(in)                :: buffer
  integer, intent(in)                 :: dest
  integer, intent(in)                 :: tag

#ifdef MPI
  call MPI_Send(buffer, 1, MPI_Real, dest, tag, comm%mpi_comm)
#endif
end subroutine send_r_sp

subroutine send_r1_sp(comm, buffer, dest, tag)
  class(mpi_communicator), intent(in) :: comm
  real(sp), intent(in)                :: buffer(:)
  integer, intent(in)                 :: dest
  integer, intent(in)                 :: tag

#ifdef MPI
  call MPI_Send(buffer, size(buffer), MPI_Real, dest, tag, comm%mpi_comm)
#endif
end subroutine send_r1_sp

subroutine send_r2_sp(comm, buffer, dest, tag)
  class(mpi_communicator), intent(in) :: comm
  real(sp), intent(in)                :: buffer(:,:)
  integer, intent(in)                 :: dest
  integer, intent(in)                 :: tag

#ifdef MPI
  call MPI_Send(buffer, size(buffer), MPI_Real, dest, tag, comm%mpi_comm)
#endif
end subroutine send_r2_sp

subroutine send_r3_sp(comm, buffer, dest, tag)
  class(mpi_communicator), intent(in) :: comm
  real(sp), intent(in)                :: buffer(:,:,:)
  integer, intent(in)                 :: dest
  integer, intent(in)                 :: tag

#ifdef MPI
  call MPI_Send(buffer, size(buffer), MPI_Real, dest, tag, comm%mpi_comm)
#endif
end subroutine send_r3_sp

subroutine send_r4_sp(comm, buffer, dest, tag)
  class(mpi_communicator), intent(in) :: comm
  real(sp), intent(in)                :: buffer(:,:,:,:)
  integer, intent(in)                 :: dest
  integer, intent(in)                 :: tag

#ifdef MPI
  call MPI_Send(buffer, size(buffer), MPI_Real, dest, tag, comm%mpi_comm)
#endif
end subroutine send_r4_sp

subroutine send_r_dp(comm, buffer, dest, tag)
  class(mpi_communicator), intent(in) :: comm
  real(dp), intent(in)                :: buffer
  integer, intent(in)                 :: dest
  integer, intent(in)                 :: tag

#ifdef MPI
  call MPI_Send(buffer, 1, MPI_Double_Precision, dest, tag, comm%mpi_comm)
#endif
end subroutine send_r_dp

subroutine send_r1_dp(comm, buffer, dest, tag)
  class(mpi_communicator), intent(in) :: comm
  real(dp), intent(in)                :: buffer(:)
  integer, intent(in)                 :: dest
  integer, intent(in)                 :: tag

#ifdef MPI
  call MPI_Send(buffer, size(buffer), MPI_Double_Precision, dest, tag, comm%mpi_comm)
#endif
end subroutine send_r1_dp

subroutine send_r2_dp(comm, buffer, dest, tag)
  class(mpi_communicator), intent(in) :: comm
  real(dp), intent(in)                :: buffer(:,:)
  integer, intent(in)                 :: dest
  integer, intent(in)                 :: tag

#ifdef MPI
  call MPI_Send(buffer, size(buffer), MPI_Double_Precision, dest, tag, comm%mpi_comm)
#endif
end subroutine send_r2_dp

subroutine send_r3_dp(comm, buffer, dest, tag)
  class(mpi_communicator), intent(in) :: comm
  real(dp), intent(in)                :: buffer(:,:,:)
  integer, intent(in)                 :: dest
  integer, intent(in)                 :: tag

#ifdef MPI
  call MPI_Send(buffer, size(buffer), MPI_Double_Precision, dest, tag, comm%mpi_comm)
#endif
end subroutine send_r3_dp

subroutine send_r4_dp(comm, buffer, dest, tag)
  class(mpi_communicator), intent(in) :: comm
  real(dp), intent(in)                :: buffer(:,:,:,:)
  integer, intent(in)                 :: dest
  integer, intent(in)                 :: tag

#ifdef MPI
  call MPI_Send(buffer, size(buffer), MPI_Double_Precision, dest, tag, comm%mpi_comm)
#endif
end subroutine send_r4_dp

subroutine send_c_sp(comm, buffer, dest, tag)
  class(mpi_communicator), intent(in) :: comm
  complex(sp), intent(in)             :: buffer
  integer, intent(in)                 :: dest
  integer, intent(in)                 :: tag

#ifdef MPI
  call MPI_Send(buffer, 1, MPI_Complex, dest, tag, comm%mpi_comm)
#endif
end subroutine send_c_sp

subroutine send_c1_sp(comm, buffer, dest, tag)
  class(mpi_communicator), intent(in) :: comm
  complex(sp), intent(in)             :: buffer(:)
  integer, intent(in)                 :: dest
  integer, intent(in)                 :: tag

#ifdef MPI
  call MPI_Send(buffer, size(buffer), MPI_Complex, dest, tag, comm%mpi_comm)
#endif
end subroutine send_c1_sp

subroutine send_c2_sp(comm, buffer, dest, tag)
  class(mpi_communicator), intent(in) :: comm
  complex(sp), intent(in)             :: buffer(:,:)
  integer, intent(in)                 :: dest
  integer, intent(in)                 :: tag

#ifdef MPI
  call MPI_Send(buffer, size(buffer), MPI_Complex, dest, tag, comm%mpi_comm)
#endif
end subroutine send_c2_sp

subroutine send_c3_sp(comm, buffer, dest, tag)
  class(mpi_communicator), intent(in) :: comm
  complex(sp), intent(in)             :: buffer(:,:,:)
  integer, intent(in)                 :: dest
  integer, intent(in)                 :: tag

#ifdef MPI
  call MPI_Send(buffer, size(buffer), MPI_Complex, dest, tag, comm%mpi_comm)
#endif
end subroutine send_c3_sp

subroutine send_c4_sp(comm, buffer, dest, tag)
  class(mpi_communicator), intent(in) :: comm
  complex(sp), intent(in)             :: buffer(:,:,:,:)
  integer, intent(in)                 :: dest
  integer, intent(in)                 :: tag

#ifdef MPI
  call MPI_Send(buffer, size(buffer), MPI_Complex, dest, tag, comm%mpi_comm)
#endif
end subroutine send_c4_sp

subroutine send_c_dp(comm, buffer, dest, tag)
  class(mpi_communicator), intent(in) :: comm
  complex(dp), intent(in)             :: buffer
  integer, intent(in)                 :: dest
  integer, intent(in)                 :: tag

#ifdef MPI
  call MPI_Send(buffer, 1, MPI_Double_Complex, dest, tag, comm%mpi_comm)
#endif
end subroutine send_c_dp

subroutine send_c1_dp(comm, buffer, dest, tag)
  class(mpi_communicator), intent(in) :: comm
  complex(dp), intent(in)             :: buffer(:)
  integer, intent(in)                 :: dest
  integer, intent(in)                 :: tag

#ifdef MPI
  call MPI_Send(buffer, size(buffer), MPI_Double_Complex, dest, tag, comm%mpi_comm)
#endif
end subroutine send_c1_dp

subroutine send_c2_dp(comm, buffer, dest, tag)
  class(mpi_communicator), intent(in) :: comm
  complex(dp), intent(in)             :: buffer(:,:)
  integer, intent(in)                 :: dest
  integer, intent(in)                 :: tag

#ifdef MPI
  call MPI_Send(buffer, size(buffer), MPI_Double_Complex, dest, tag, comm%mpi_comm)
#endif
end subroutine send_c2_dp

subroutine send_c3_dp(comm, buffer, dest, tag)
  class(mpi_communicator), intent(in) :: comm
  complex(dp), intent(in)             :: buffer(:,:,:)
  integer, intent(in)                 :: dest
  integer, intent(in)                 :: tag

#ifdef MPI
  call MPI_Send(buffer, size(buffer), MPI_Double_Complex, dest, tag, comm%mpi_comm)
#endif
end subroutine send_c3_dp

subroutine send_c4_dp(comm, buffer, dest, tag)
  class(mpi_communicator), intent(in) :: comm
  complex(dp), intent(in)             :: buffer(:,:,:,:)
  integer, intent(in)                 :: dest
  integer, intent(in)                 :: tag

#ifdef MPI
  call MPI_Send(buffer, size(buffer), MPI_Double_Complex, dest, tag, comm%mpi_comm)
#endif
end subroutine send_c4_dp


!******************************************************************************** 
!
! Interface routines for MPI_Recv
!
!******************************************************************************** 
subroutine recv_char(comm, buffer, source, tag)
  class(mpi_communicator), intent(in) :: comm
  character(len=*), intent(in)        :: buffer
  integer, intent(in)                 :: source
  integer, intent(in)                 :: tag
  !local
#ifdef MPI
  type(MPI_Status)                    :: status

  call MPI_Recv(buffer, len(buffer), MPI_Character, source, tag, comm%mpi_comm, status)
#endif
end subroutine recv_char

subroutine recv_log(comm, buffer, source, tag)
  class(mpi_communicator), intent(in) :: comm
  logical, intent(inout)              :: buffer
  integer, intent(in)                 :: source
  integer, intent(in)                 :: tag
  !local
#ifdef MPI
  type(MPI_Status)                    :: status

  call MPI_Recv(buffer, 1, MPI_Logical, source, tag, comm%mpi_comm, status)
#endif
end subroutine recv_log

subroutine recv_int(comm, buffer, source, tag)
  class(mpi_communicator), intent(in) :: comm
  integer, intent(inout)              :: buffer
  integer, intent(in)                 :: source
  integer, intent(in)                 :: tag
  !local
#ifdef MPI
  type(MPI_Status)                    :: status

  call MPI_Recv(buffer, 1, MPI_Integer, source, tag, comm%mpi_comm, status)
#endif
end subroutine recv_int

subroutine recv_int_1d(comm, buffer, source, tag)
  class(mpi_communicator), intent(in) :: comm
  integer, intent(inout)              :: buffer(:)
  integer, intent(in)                 :: source
  integer, intent(in)                 :: tag
  !local
#ifdef MPI
  type(MPI_Status)                    :: status

  call MPI_Recv(buffer, size(buffer), MPI_Integer, source, tag, comm%mpi_comm, status)
#endif
end subroutine recv_int_1d

subroutine recv_int8(comm, buffer, source, tag)
  class(mpi_communicator), intent(in) :: comm
  integer(i8), intent(inout)          :: buffer
  integer, intent(in)                 :: source
  integer, intent(in)                 :: tag
  !local
#ifdef MPI
  type(MPI_Status)                    :: status

  call MPI_Recv(buffer, 1, MPI_Integer8, source, tag, comm%mpi_comm, status)
#endif
end subroutine recv_int8

subroutine recv_int8_1d(comm, buffer, source, tag)
  class(mpi_communicator), intent(in) :: comm
  integer(i8), intent(inout)          :: buffer(:)
  integer, intent(in)                 :: source
  integer, intent(in)                 :: tag
  !local
#ifdef MPI
  type(MPI_Status)                    :: status

  call MPI_Recv(buffer, size(buffer), MPI_Integer8, source, tag, comm%mpi_comm, status)
#endif
end subroutine recv_int8_1d

subroutine recv_r_sp(comm, buffer, source, tag)
  class(mpi_communicator), intent(in) :: comm
  real(sp), intent(inout)             :: buffer
  integer, intent(in)                 :: source
  integer, intent(in)                 :: tag
  !local
#ifdef MPI
  type(MPI_Status)                    :: status

  call MPI_Recv(buffer, 1, MPI_Real, source, tag, comm%mpi_comm, status)
#endif
end subroutine recv_r_sp

subroutine recv_r1_sp(comm, buffer, source, tag)
  class(mpi_communicator), intent(in) :: comm
  real(sp), intent(inout)             :: buffer(:)
  integer, intent(in)                 :: source
  integer, intent(in)                 :: tag
  !local
#ifdef MPI
  type(MPI_Status)                    :: status

  call MPI_Recv(buffer, size(buffer), MPI_Real, source, tag, comm%mpi_comm, status)
#endif
end subroutine recv_r1_sp

subroutine recv_r2_sp(comm, buffer, source, tag)
  class(mpi_communicator), intent(in) :: comm
  real(sp), intent(inout)             :: buffer(:,:)
  integer, intent(in)                 :: source
  integer, intent(in)                 :: tag
  !local
#ifdef MPI
  type(MPI_Status)                    :: status

  call MPI_Recv(buffer, size(buffer), MPI_Real, source, tag, comm%mpi_comm, status)
#endif
end subroutine recv_r2_sp

subroutine recv_r3_sp(comm, buffer, source, tag)
  class(mpi_communicator), intent(in) :: comm
  real(sp), intent(inout)             :: buffer(:,:,:)
  integer, intent(in)                 :: source
  integer, intent(in)                 :: tag
  !local
#ifdef MPI
  type(MPI_Status)                    :: status

  call MPI_Recv(buffer, size(buffer), MPI_Real, source, tag, comm%mpi_comm, status)
#endif
end subroutine recv_r3_sp

subroutine recv_r4_sp(comm, buffer, source, tag)
  class(mpi_communicator), intent(in) :: comm
  real(sp), intent(inout)             :: buffer(:,:,:,:)
  integer, intent(in)                 :: source
  integer, intent(in)                 :: tag
  !local
#ifdef MPI
  type(MPI_Status)                    :: status

  call MPI_Recv(buffer, size(buffer), MPI_Real, source, tag, comm%mpi_comm, status)
#endif
end subroutine recv_r4_sp

subroutine recv_r_dp(comm, buffer, source, tag)
  class(mpi_communicator), intent(in) :: comm
  real(dp), intent(inout)             :: buffer
  integer, intent(in)                 :: source
  integer, intent(in)                 :: tag
  !local
#ifdef MPI
  type(MPI_Status)                    :: status

  call MPI_Recv(buffer, 1, MPI_Double_Precision, source, tag, comm%mpi_comm, status)
#endif
end subroutine recv_r_dp

subroutine recv_r1_dp(comm, buffer, source, tag)
  class(mpi_communicator), intent(in) :: comm
  real(dp), intent(inout)             :: buffer(:)
  integer, intent(in)                 :: source
  integer, intent(in)                 :: tag
  !local
#ifdef MPI
  type(MPI_Status)                    :: status

  call MPI_Recv(buffer, size(buffer), MPI_Double_Precision, source, tag, comm%mpi_comm, status)
#endif
end subroutine recv_r1_dp

subroutine recv_r2_dp(comm, buffer, source, tag)
  class(mpi_communicator), intent(in) :: comm
  real(dp), intent(inout)             :: buffer(:,:)
  integer, intent(in)                 :: source
  integer, intent(in)                 :: tag
  !local
#ifdef MPI
  type(MPI_Status)                    :: status

  call MPI_Recv(buffer, size(buffer), MPI_Double_Precision, source, tag, comm%mpi_comm, status)
#endif
end subroutine recv_r2_dp

subroutine recv_r3_dp(comm, buffer, source, tag)
  class(mpi_communicator), intent(in) :: comm
  real(dp), intent(inout)             :: buffer(:,:,:)
  integer, intent(in)                 :: source
  integer, intent(in)                 :: tag
  !local
#ifdef MPI
  type(MPI_Status)                    :: status

  call MPI_Recv(buffer, size(buffer), MPI_Double_Precision, source, tag, comm%mpi_comm, status)
#endif
end subroutine recv_r3_dp

subroutine recv_r4_dp(comm, buffer, source, tag)
  class(mpi_communicator), intent(in) :: comm
  real(dp), intent(inout)             :: buffer(:,:,:,:)
  integer, intent(in)                 :: source
  integer, intent(in)                 :: tag
  !local
#ifdef MPI
  type(MPI_Status)                    :: status

  call MPI_Recv(buffer, size(buffer), MPI_Double_Precision, source, tag, comm%mpi_comm, status)
#endif
end subroutine recv_r4_dp

subroutine recv_c_sp(comm, buffer, source, tag)
  class(mpi_communicator), intent(in) :: comm
  complex(sp), intent(inout)          :: buffer
  integer, intent(in)                 :: source
  integer, intent(in)                 :: tag
  !local
#ifdef MPI
  type(MPI_Status)                    :: status

  call MPI_Recv(buffer, 1, MPI_Complex, source, tag, comm%mpi_comm, status)
#endif
end subroutine recv_c_sp

subroutine recv_c1_sp(comm, buffer, source, tag)
  class(mpi_communicator), intent(in) :: comm
  complex(sp), intent(inout)          :: buffer(:)
  integer, intent(in)                 :: source
  integer, intent(in)                 :: tag
  !local
#ifdef MPI
  type(MPI_Status)                    :: status

  call MPI_Recv(buffer, size(buffer), MPI_Complex, source, tag, comm%mpi_comm, status)
#endif
end subroutine recv_c1_sp

subroutine recv_c2_sp(comm, buffer, source, tag)
  class(mpi_communicator), intent(in) :: comm
  complex(sp), intent(inout)          :: buffer(:,:)
  integer, intent(in)                 :: source
  integer, intent(in)                 :: tag
  !local
#ifdef MPI
  type(MPI_Status)                    :: status

  call MPI_Recv(buffer, size(buffer), MPI_Complex, source, tag, comm%mpi_comm, status)
#endif
end subroutine recv_c2_sp

subroutine recv_c3_sp(comm, buffer, source, tag)
  class(mpi_communicator), intent(in) :: comm
  complex(sp), intent(inout)          :: buffer(:,:,:)
  integer, intent(in)                 :: source
  integer, intent(in)                 :: tag
  !local
#ifdef MPI
  type(MPI_Status)                    :: status

  call MPI_Recv(buffer, size(buffer), MPI_Complex, source, tag, comm%mpi_comm, status)
#endif
end subroutine recv_c3_sp

subroutine recv_c4_sp(comm, buffer, source, tag)
  class(mpi_communicator), intent(in) :: comm
  complex(sp), intent(inout)          :: buffer(:,:,:,:)
  integer, intent(in)                 :: source
  integer, intent(in)                 :: tag
  !local
#ifdef MPI
  type(MPI_Status)                    :: status

  call MPI_Recv(buffer, size(buffer), MPI_Complex, source, tag, comm%mpi_comm, status)
#endif
end subroutine recv_c4_sp

subroutine recv_c_dp(comm, buffer, source, tag)
  class(mpi_communicator), intent(in) :: comm
  complex(dp), intent(inout)          :: buffer
  integer, intent(in)                 :: source
  integer, intent(in)                 :: tag
  !local
#ifdef MPI
  type(MPI_Status)                    :: status

  call MPI_Recv(buffer, 1, MPI_Double_Complex, source, tag, comm%mpi_comm, status)
#endif
end subroutine recv_c_dp

subroutine recv_c1_dp(comm, buffer, source, tag)
  class(mpi_communicator), intent(in) :: comm
  complex(dp), intent(inout)          :: buffer(:)
  integer, intent(in)                 :: source
  integer, intent(in)                 :: tag
  !local
#ifdef MPI
  type(MPI_Status)                    :: status

  call MPI_Recv(buffer, size(buffer), MPI_Double_Complex, source, tag, comm%mpi_comm, status)
#endif
end subroutine recv_c1_dp

subroutine recv_c2_dp(comm, buffer, source, tag)
  class(mpi_communicator), intent(in) :: comm
  complex(dp), intent(inout)          :: buffer(:,:)
  integer, intent(in)                 :: source
  integer, intent(in)                 :: tag
  !local
#ifdef MPI
  type(MPI_Status)                    :: status

  call MPI_Recv(buffer, size(buffer), MPI_Double_Complex, source, tag, comm%mpi_comm, status)
#endif
end subroutine recv_c2_dp

subroutine recv_c3_dp(comm, buffer, source, tag)
  class(mpi_communicator), intent(in) :: comm
  complex(dp), intent(inout)          :: buffer(:,:,:)
  integer, intent(in)                 :: source
  integer, intent(in)                 :: tag
  !local
#ifdef MPI
  type(MPI_Status)                    :: status

  call MPI_Recv(buffer, size(buffer), MPI_Double_Complex, source, tag, comm%mpi_comm, status)
#endif
end subroutine recv_c3_dp

subroutine recv_c4_dp(comm, buffer, source, tag)
  class(mpi_communicator), intent(in) :: comm
  complex(dp), intent(inout)          :: buffer(:,:,:,:)
  integer, intent(in)                 :: source
  integer, intent(in)                 :: tag
  !local
#ifdef MPI
  type(MPI_Status)                    :: status

  call MPI_Recv(buffer, size(buffer), MPI_Double_Complex, source, tag, comm%mpi_comm, status)
#endif
end subroutine recv_c4_dp


!******************************************************************************** 
!
! Interface routines for MPI_Bcast
!
!******************************************************************************** 
subroutine bcast_char(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  character(len=*), intent(inout)     :: buffer
  integer, intent(in)                 :: root
  !local
  integer                             :: ierror
  
#ifdef MPI
  call MPI_Bcast(buffer, len(buffer), MPI_Character, root, comm%mpi_comm, ierror)
#endif
end subroutine bcast_char

subroutine bcast_char_1d(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  character(len=*), intent(inout)     :: buffer(:)
  integer, intent(in)                 :: root
  !local
  integer                             :: ierror, size_

  if (size(buffer) > 0) then
    size_ = len(buffer(1)) * size(buffer)
  else
    size_ = 0
  end if
  
#ifdef MPI
  call MPI_Bcast(buffer, size_, MPI_Character, root, comm%mpi_comm, ierror)
#endif
end subroutine bcast_char_1d

subroutine bcast_log(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  logical, intent(inout)              :: buffer
  integer, intent(in)                 :: root
  !local
  integer                             :: ierror
  
#ifdef MPI
  call MPI_Bcast(buffer, 1, MPI_Logical, root, comm%mpi_comm, ierror)
#endif
end subroutine bcast_log

subroutine bcast_int(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  integer, intent(inout)              :: buffer
  integer, intent(in)                 :: root
  !local
  integer                             :: ierror
  
#ifdef MPI
  call MPI_Bcast(buffer, 1, MPI_Integer, root, comm%mpi_comm, ierror)
#endif
end subroutine bcast_int

subroutine bcast_int_1d(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  integer, intent(inout)              :: buffer(:)
  integer, intent(in)                 :: root
  !local
  integer                             :: ierror
  
#ifdef MPI
  call MPI_Bcast(buffer, size(buffer), MPI_Integer, root, comm%mpi_comm, ierror)
#endif
end subroutine bcast_int_1d

subroutine bcast_int8(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  integer(i8), intent(inout)          :: buffer
  integer, intent(in)                 :: root
  !local
  integer                             :: ierror
  
#ifdef MPI
  call MPI_Bcast(buffer, 1, MPI_Integer8, root, comm%mpi_comm, ierror)
#endif
end subroutine bcast_int8

subroutine bcast_int8_1d(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  integer(i8), intent(inout)          :: buffer(:)
  integer, intent(in)                 :: root
  !local
  integer                             :: ierror
  
#ifdef MPI
  call MPI_Bcast(buffer, size(buffer), MPI_Integer8, root, comm%mpi_comm, ierror)
#endif
end subroutine bcast_int8_1d

subroutine bcast_int_2d(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  integer, intent(inout)              :: buffer(:,:)
  integer, intent(in)                 :: root
  !local
  integer                             :: ierror
  
#ifdef MPI
  call MPI_Bcast(buffer, size(buffer), MPI_Integer, root, comm%mpi_comm, ierror)
#endif
end subroutine bcast_int_2d

subroutine bcast_r_sp(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  real(sp), intent(inout)             :: buffer
  integer, intent(in)                 :: root
  !local
  integer                             :: ierror
  
#ifdef MPI
  call MPI_Bcast(buffer, 1, MPI_Real, root, comm%mpi_comm, ierror)
#endif
end subroutine bcast_r_sp

subroutine bcast_r1_sp(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  real(sp), intent(inout)             :: buffer(:)
  integer, intent(in)                 :: root
  !local
  integer                             :: ierror
  
#ifdef MPI
  call MPI_Bcast(buffer, size(buffer), MPI_Real, root, comm%mpi_comm, ierror)
#endif
end subroutine bcast_r1_sp

subroutine bcast_r2_sp(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  real(sp), intent(inout)             :: buffer(:,:)
  integer, intent(in)                 :: root
  !local
  integer                             :: ierror
  
#ifdef MPI
  call MPI_Bcast(buffer, size(buffer), MPI_Real, root, comm%mpi_comm, ierror)
#endif
end subroutine bcast_r2_sp

subroutine bcast_r3_sp(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  real(sp), intent(inout)             :: buffer(:,:,:)
  integer, intent(in)                 :: root
  !local
  integer                             :: ierror
  
#ifdef MPI
  call MPI_Bcast(buffer, size(buffer), MPI_real, root, comm%mpi_comm, ierror)
#endif
end subroutine bcast_r3_sp

subroutine bcast_r4_sp(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  real(sp), intent(inout)             :: buffer(:,:,:,:)
  integer, intent(in)                 :: root
  !local
  integer                             :: ierror
  
#ifdef MPI
  call MPI_Bcast(buffer, size(buffer), MPI_Real, root, comm%mpi_comm, ierror)
#endif
end subroutine bcast_r4_sp

subroutine bcast_r5_sp(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  real(sp), intent(inout)             :: buffer(:,:,:,:,:)
  integer, intent(in)                 :: root
  !local
  integer                             :: ierror
  
#ifdef MPI
  call MPI_Bcast(buffer, size(buffer), MPI_Real, root, comm%mpi_comm, ierror)
#endif
end subroutine bcast_r5_sp

subroutine bcast_r_dp(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  real(dp), intent(inout)             :: buffer
  integer, intent(in)                 :: root
  !local
  integer                             :: ierror
  
#ifdef MPI
  call MPI_Bcast(buffer, 1, MPI_Double_Precision, root, comm%mpi_comm, ierror)
#endif
end subroutine bcast_r_dp

subroutine bcast_r1_dp(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  real(dp), intent(inout)             :: buffer(:)
  integer, intent(in)                 :: root
  !local
  integer                             :: ierror
  
#ifdef MPI
  call MPI_Bcast(buffer, size(buffer), MPI_Double_Precision, root, comm%mpi_comm, ierror)
#endif
end subroutine bcast_r1_dp

subroutine bcast_r2_dp(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  real(dp), intent(inout)             :: buffer(:,:)
  integer, intent(in)                 :: root
  !local
  integer                             :: ierror
  
#ifdef MPI
  call MPI_Bcast(buffer, size(buffer), MPI_Double_Precision, root, comm%mpi_comm, ierror)
#endif
end subroutine bcast_r2_dp

subroutine bcast_r3_dp(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  real(dp), intent(inout)             :: buffer(:,:,:)
  integer, intent(in)                 :: root
  !local
  integer                             :: ierror
  
#ifdef MPI
  call MPI_Bcast(buffer, size(buffer), MPI_Double_Precision, root, comm%mpi_comm, ierror)
#endif
end subroutine bcast_r3_dp

subroutine bcast_r4_dp(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  real(dp), intent(inout)             :: buffer(:,:,:,:)
  integer, intent(in)                 :: root
  !local
  integer                             :: ierror
  
#ifdef MPI
  call MPI_Bcast(buffer, size(buffer), MPI_Double_Precision, root, comm%mpi_comm, ierror)
#endif
end subroutine bcast_r4_dp

subroutine bcast_r5_dp(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  real(dp), intent(inout)             :: buffer(:,:,:,:,:)
  integer, intent(in)                 :: root
  !local
  integer                             :: ierror
  integer                             :: i, j, k
  
#ifdef MPI
  do i = 1, size(buffer, 5)
    do j = 1, size(buffer, 4)
      do k = 1, size(buffer, 3)
        call MPI_Bcast(buffer(:,:,k,j,i), size(buffer(:,:,k,j,i)), MPI_Double_Precision, root, comm%mpi_comm, ierror)
      end do
    end do
  end do
#endif
end subroutine bcast_r5_dp

subroutine bcast_c_sp(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  complex(sp), intent(inout)          :: buffer
  integer, intent(in)                 :: root
  !local
  integer                             :: ierror
  
#ifdef MPI
  call MPI_Bcast(buffer, 1, MPI_Complex, root, comm%mpi_comm, ierror)
#endif
end subroutine bcast_c_sp

subroutine bcast_c1_sp(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  complex(sp), intent(inout)          :: buffer(:)
  integer, intent(in)                 :: root
  !local
  integer                             :: ierror
  
#ifdef MPI
  call MPI_Bcast(buffer, size(buffer), MPI_Complex, root, comm%mpi_comm, ierror)
#endif
end subroutine bcast_c1_sp

subroutine bcast_c2_sp(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  complex(sp), intent(inout)          :: buffer(:,:)
  integer, intent(in)                 :: root
  !local
  integer                             :: ierror
  
#ifdef MPI
  call MPI_Bcast(buffer, size(buffer), MPI_Complex, root, comm%mpi_comm, ierror)
#endif
end subroutine bcast_c2_sp

subroutine bcast_c3_sp(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  complex(sp), intent(inout)          :: buffer(:,:,:)
  integer, intent(in)                 :: root
  !local
  integer                             :: ierror
  
#ifdef MPI
  call MPI_Bcast(buffer, size(buffer), MPI_Complex, root, comm%mpi_comm, ierror)
#endif
end subroutine bcast_c3_sp

subroutine bcast_c4_sp(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  complex(sp), intent(inout)          :: buffer(:,:,:,:)
  integer, intent(in)                 :: root
  !local
  integer                             :: ierror
  
#ifdef MPI
  call MPI_Bcast(buffer, size(buffer), MPI_Complex, root, comm%mpi_comm, ierror)
#endif
end subroutine bcast_c4_sp

subroutine bcast_c5_sp(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  complex(sp), intent(inout)          :: buffer(:,:,:,:,:)
  integer, intent(in)                 :: root
  !local
  integer                             :: ierror
  integer                             :: i, j, k
  
#ifdef MPI
  do i = 1, size(buffer, 5)
    do j = 1, size(buffer, 4)
      do k = 1, size(buffer, 3)
        call MPI_Bcast(buffer(:,:,k,j,i), size(buffer(:,:,k,j,i)), MPI_Complex, root, comm%mpi_comm, ierror)
      end do
    end do
  end do
#endif
end subroutine bcast_c5_sp

subroutine bcast_c_dp(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  complex(dp), intent(inout)          :: buffer
  integer, intent(in)                 :: root
  !local
  integer                             :: ierror
  
#ifdef MPI
  call MPI_Bcast(buffer, 1, MPI_Double_Complex, root, comm%mpi_comm, ierror)
#endif
end subroutine bcast_c_dp

subroutine bcast_c1_dp(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  complex(dp), intent(inout)          :: buffer(:)
  integer, intent(in)                 :: root
  !local
  integer                             :: ierror
  
#ifdef MPI
  call MPI_Bcast(buffer, size(buffer), MPI_Double_Complex, root, comm%mpi_comm, ierror)
#endif
end subroutine bcast_c1_dp

subroutine bcast_c2_dp(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  complex(dp), intent(inout)          :: buffer(:,:)
  integer, intent(in)                 :: root
  !local
  integer                             :: ierror
  
#ifdef MPI
  call MPI_Bcast(buffer, size(buffer), MPI_Double_Complex, root, comm%mpi_comm, ierror)
#endif
end subroutine bcast_c2_dp

subroutine bcast_c3_dp(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  complex(dp), intent(inout)          :: buffer(:,:,:)
  integer, intent(in)                 :: root
  !local
  integer                             :: ierror
  
#ifdef MPI
  call MPI_Bcast(buffer, size(buffer), MPI_Double_Complex, root, comm%mpi_comm, ierror)
#endif
end subroutine bcast_c3_dp

subroutine bcast_c4_dp(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  complex(dp), intent(inout)          :: buffer(:,:,:,:)
  integer, intent(in)                 :: root
  !local
  integer                             :: ierror
  
#ifdef MPI
  call MPI_Bcast(buffer, size(buffer), MPI_Double_Complex, root, comm%mpi_comm, ierror)
#endif
end subroutine bcast_c4_dp

subroutine bcast_c5_dp(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  complex(dp), intent(inout)          :: buffer(:,:,:,:,:)
  integer, intent(in)                 :: root
  !local
  integer                             :: ierror
  integer                             :: i, j, k
  
#ifdef MPI
  do i = 1, size(buffer, 5)
    do j = 1, size(buffer, 4)
      do k = 1, size(buffer, 3)
        call MPI_Bcast(buffer(:,:,k,j,i), size(buffer(:,:,k,j,i)), MPI_Double_Complex, root, comm%mpi_comm, ierror)
      end do
    end do
  end do
#endif
end subroutine bcast_c5_dp


!******************************************************************************** 
!
! Interface for MPI_Gather, MPI_Allgather
!
!******************************************************************************** 
subroutine gather_int(comm, buffer, buffer_, root)
  class(mpi_communicator), intent(in) :: comm
  integer, intent(in)                 :: buffer
  integer, allocatable, intent(out)   :: buffer_(:)
  integer, optional, intent(in)       :: root

#ifdef MPI
  if (present(root)) then
    if (comm%mpirank == root) then 
      allocate(buffer_(comm%mpisize))
    else
      allocate(buffer_(1))
    end if
    call MPI_Gather(buffer, 1, mpi_integer, buffer_, 1, mpi_integer, root, comm%mpi_comm)
  else
    allocate(buffer_(comm%mpisize))
    call MPI_Allgather(buffer, 1, mpi_integer, buffer_, 1, mpi_integer, comm%mpi_comm)
  end if
#else
  allocate(buffer_(1))
  buffer_(1) = buffer
#endif
end subroutine gather_int

subroutine gather_r_sp(comm, buffer, buffer_, root)
  class(mpi_communicator), intent(in) :: comm
  real(sp), intent(in)                :: buffer
  real(sp), allocatable, intent(out)  :: buffer_(:)
  integer, optional, intent(in)       :: root

#ifdef MPI
  if (present(root)) then
    if (comm%mpirank == root) then 
      allocate(buffer_(comm%mpisize))
    else
      allocate(buffer_(1))
    end if
    call MPI_Gather(buffer, 1, MPI_real, buffer_, 1, MPI_real, root, comm%mpi_comm)
  else
    allocate(buffer_(comm%mpisize))
    call MPI_Allgather(buffer, 1, MPI_real, buffer_, 1, MPI_real, comm%mpi_comm)
  end if
#else
  allocate(buffer_(1))
  buffer_(1) = buffer
#endif
end subroutine gather_r_sp

subroutine gather_r1_sp(comm, buffer, buffer_, root)
  class(mpi_communicator), intent(in) :: comm
  real(sp), allocatable, intent(in)   :: buffer(:)
  real(sp), allocatable, intent(out)  :: buffer_(:)
  integer, optional, intent(in)       :: root

#ifdef MPI
  if (present(root)) then
    if (comm%mpirank == root) then
       allocate(buffer_(size(buffer)*comm%mpisize))
    else
      allocate(buffer_(1))
    end if
    call MPI_Gather(buffer, size(buffer), MPI_Real, buffer_, size(buffer), MPI_Real, root, comm%mpi_comm)
  else
    allocate(buffer_(size(buffer)*comm%mpisize))
    call MPI_Allgather(buffer, size(buffer), MPI_Real, buffer_, size(buffer), MPI_Real, comm%mpi_comm)
  end if
#else
  allocate(buffer_(size(buffer)))
  buffer_ = buffer
#endif
end subroutine gather_r1_sp

subroutine gather_r2_sp(comm, buffer, buffer_, root)
  class(mpi_communicator), intent(in) :: comm
  real(sp), allocatable, intent(in)   :: buffer(:,:)
  real(sp), allocatable, intent(out)  :: buffer_(:,:)
  integer, optional, intent(in)       :: root

#ifdef MPI
  if (present(root)) then
    if (comm%mpirank == root) then 
      allocate(buffer_(size(buffer,1), size(buffer,2)*comm%mpisize))
    else
      allocate(buffer_(1,1))
    end if
    call MPI_Gather(buffer, size(buffer), MPI_Real, buffer_, size(buffer), MPI_Real, root, comm%mpi_comm)
  else
    allocate(buffer_(size(buffer,1), size(buffer,2)*comm%mpisize))
    call MPI_Allgather(buffer, size(buffer), MPI_Real, buffer_, size(buffer), MPI_Real, comm%mpi_comm)
  end if
#else
  allocate(buffer_, source=buffer)
#endif
end subroutine gather_r2_sp

subroutine gather_r_dp(comm, buffer, buffer_, root)
  class(mpi_communicator), intent(in) :: comm
  real(dp), intent(in)                :: buffer
  real(dp), allocatable, intent(out)  :: buffer_(:)
  integer, optional, intent(in)       :: root

#ifdef MPI
  if (present(root)) then
    if (comm%mpirank == root) then
      allocate(buffer_(comm%mpisize))
    else
      allocate(buffer_(1))
    end if
    call MPI_Gather(buffer, 1, MPI_Double_Precision, buffer_, 1, MPI_Double_Precision, root, comm%mpi_comm)
  else
    allocate(buffer_(comm%mpisize))
    call MPI_Allgather(buffer, 1, MPI_Double_Precision, buffer_, 1, MPI_Double_Precision, comm%mpi_comm)
  end if
#else
  allocate(buffer_(1))
  buffer_(1) = buffer
#endif
end subroutine gather_r_dp

subroutine gather_r1_dp(comm, buffer, buffer_, root)
  class(mpi_communicator), intent(in) :: comm
  real(dp), allocatable, intent(in)   :: buffer(:)
  real(dp), allocatable, intent(out)  :: buffer_(:)
  integer, optional, intent(in)       :: root

#ifdef MPI
  if (present(root)) then
    if (comm%mpirank == root) then
      allocate(buffer_(size(buffer)*comm%mpisize))
    else
      allocate(buffer_(1))
    end if
    call MPI_Gather(buffer, size(buffer), MPI_Double_Precision, buffer_, size(buffer), MPI_Double_Precision, root, comm%mpi_comm)
  else
    allocate(buffer_(size(buffer)*comm%mpisize))
    call MPI_Allgather(buffer, size(buffer), MPI_Double_Precision, buffer_, size(buffer), MPI_Double_Precision, comm%mpi_comm)
  end if
#else
  allocate(buffer_(size(buffer)))
  buffer_ = buffer
#endif
end subroutine gather_r1_dp

subroutine gather_r2_dp(comm, buffer, buffer_, root)
  class(mpi_communicator), intent(in) :: comm
  real(dp), allocatable, intent(in)   :: buffer(:,:)
  real(dp), allocatable, intent(out)  :: buffer_(:,:)
  integer, optional, intent(in)       :: root

#ifdef MPI
  if (present(root)) then
    if (comm%mpirank == root) then 
      allocate(buffer_(size(buffer,1), size(buffer,2)*comm%mpisize))
    else
      allocate(buffer_(1,1))
    end if
    call MPI_Gather(buffer, size(buffer), MPI_Double_Precision, buffer_, size(buffer), MPI_Double_Precision, root, comm%mpi_comm)
  else
    allocate(buffer_(size(buffer,1), size(buffer,2)*comm%mpisize))
    call MPI_Allgather(buffer, size(buffer), MPI_Double_Precision, buffer_, size(buffer), MPI_Double_Precision, comm%mpi_comm)
  end if
#else
  allocate(buffer_, source=buffer)
#endif
end subroutine gather_r2_dp

subroutine gather_c_sp(comm, buffer, buffer_, root)
  class(mpi_communicator), intent(in)   :: comm
  complex(sp), intent(in)               :: buffer
  complex(sp), allocatable, intent(out) :: buffer_(:)
  integer, optional, intent(in)         :: root

#ifdef MPI
  if (present(root)) then
    if (comm%mpirank == root) then 
      allocate(buffer_(comm%mpisize))
    else
      allocate(buffer_(1))
    end if
    call MPI_Gather(buffer, 1, MPI_Complex , buffer_, 1, MPI_Complex, root, comm%mpi_comm)
  else
    allocate(buffer_(comm%mpisize))
    call MPI_Allgather(buffer, 1, MPI_Complex, buffer_, 1, MPI_Complex, comm%mpi_comm)
  end if
#else
  allocate(buffer_(1))
  buffer_(1) = buffer
#endif
end subroutine gather_c_sp

subroutine gather_c1_sp(comm, buffer, buffer_, root)
  class(mpi_communicator), intent(in)   :: comm
  complex(sp), allocatable, intent(in)  :: buffer(:)
  complex(sp), allocatable, intent(out) :: buffer_(:)
  integer, optional, intent(in)         :: root

#ifdef MPI
  if (present(root)) then
    if (comm%mpirank == root) then 
      allocate(buffer_(size(buffer)*comm%mpisize))
    else
      allocate(buffer_(1))
    end if
    call MPI_Gather(buffer, size(buffer), MPI_Complex, buffer_, size(buffer), MPI_Complex, root, comm%mpi_comm)
  else
    allocate(buffer_(size(buffer)*comm%mpisize))
    call MPI_Allgather(buffer, size(buffer), MPI_Complex, buffer_, size(buffer), MPI_Complex, comm%mpi_comm)
  end if
#else
  allocate(buffer_(size(buffer)))
  buffer_ = buffer
#endif
end subroutine gather_c1_sp

subroutine gather_c2_sp(comm, buffer, buffer_, root)
  class(mpi_communicator), intent(in)   :: comm
  complex(sp), allocatable, intent(in)  :: buffer(:,:)
  complex(sp), allocatable, intent(out) :: buffer_(:,:)
  integer, optional, intent(in)         :: root

#ifdef MPI
  if (present(root)) then
    if (comm%mpirank == root) then 
      allocate(buffer_(size(buffer,1), size(buffer,2)*comm%mpisize))
    else
      allocate(buffer_(1,1))
    end if
    call MPI_Gather(buffer, size(buffer), MPI_Complex, buffer_, size(buffer), MPI_Complex, root, comm%mpi_comm)
  else
    allocate(buffer_(size(buffer,1), size(buffer,2)*comm%mpisize))
    call MPI_Allgather(buffer, size(buffer), MPI_Complex, buffer_, size(buffer), MPI_Complex, comm%mpi_comm)
  end if
#else
  allocate(buffer_, source=buffer)
#endif
end subroutine gather_c2_sp

subroutine gather_c_dp(comm, buffer, buffer_, root)
  class(mpi_communicator), intent(in)   :: comm
  complex(dp), intent(in)               :: buffer
  complex(dp), allocatable, intent(out) :: buffer_(:)
  integer, optional, intent(in)         :: root

#ifdef MPI
  if (present(root)) then
    if (comm%mpirank == root) then
      allocate(buffer_(comm%mpisize))
    else
      allocate(buffer_(1))
    end if
    call MPI_Gather(buffer, 1, MPI_Double_Complex , buffer_, 1, MPI_Double_Complex, root, comm%mpi_comm)
  else
    allocate(buffer_(comm%mpisize))
    call MPI_Allgather(buffer, 1, MPI_Double_Complex, buffer_, 1, MPI_Double_Complex, comm%mpi_comm)
  end if
#else
  allocate(buffer_(1))
  buffer_(1) = buffer
#endif
end subroutine gather_c_dp

subroutine gather_c1_dp(comm, buffer, buffer_, root)
  class(mpi_communicator), intent(in)   :: comm
  complex(dp), allocatable, intent(in)  :: buffer(:)
  complex(dp), allocatable, intent(out) :: buffer_(:)
  integer, optional, intent(in)         :: root

#ifdef MPI
  if (present(root)) then
    if (comm%mpirank == root) then 
      allocate(buffer_(size(buffer)*comm%mpisize))
    else
      allocate(buffer_(1))
    end if
    call MPI_Gather(buffer, size(buffer), MPI_Double_Complex, buffer_, size(buffer), MPI_Double_Complex, root, comm%mpi_comm)
  else
    allocate(buffer_(size(buffer)*comm%mpisize))
    call MPI_Allgather(buffer, size(buffer), MPI_Double_Complex, buffer_, size(buffer), MPI_Double_Complex, comm%mpi_comm)
  end if
#else
  allocate(buffer_(size(buffer)))
  buffer_ = buffer
#endif
end subroutine gather_c1_dp

subroutine gather_c2_dp(comm, buffer, buffer_, root)
  class(mpi_communicator), intent(in)   :: comm
  complex(dp), allocatable, intent(in)  :: buffer(:,:)
  complex(dp), allocatable, intent(out) :: buffer_(:,:)
  integer, optional, intent(in)         :: root

#ifdef MPI
  if (present(root)) then
    if (comm%mpirank == root) then
      allocate(buffer_(size(buffer,1), size(buffer,2)*comm%mpisize))
    else 
      allocate(buffer_(1,1))
    end if
    call MPI_Gather(buffer, size(buffer), MPI_Double_Complex, buffer_, size(buffer), MPI_Double_Complex, root, comm%mpi_comm)
  else
    allocate(buffer_(size(buffer,1), size(buffer,2)*comm%mpisize))
    call MPI_Allgather(buffer, size(buffer), MPI_Double_Complex, buffer_, size(buffer), MPI_Double_Complex, comm%mpi_comm)
  end if
#else
  allocate(buffer_, source=buffer)
#endif
end subroutine gather_c2_dp


!******************************************************************************** 
!
! Interface routine for MPI_Allreduce with mpi_sum
!
!******************************************************************************** 
subroutine mpi_sum_i(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  integer, intent(inout)              :: buffer
  integer, optional, intent(in)       :: root
  !
#ifdef MPI
  if (present(root)) then
    if (comm%mpirank == root) then
      call MPI_Reduce(mpi_in_place, buffer, 1, mpi_integer, mpi_sum, root, comm%mpi_comm)
    else
      call MPI_Reduce(buffer, buffer, 1, mpi_integer, mpi_sum, root, comm%mpi_comm)
    end if
  else
    call MPI_Allreduce(mpi_in_place, buffer, 1, mpi_integer, mpi_sum, comm%mpi_comm)
  end if
#endif  
end subroutine mpi_sum_i

subroutine mpi_sum_i1(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  integer, intent(inout)              :: buffer(:)
  integer, optional, intent(in)       :: root
  
#ifdef MPI
  if (present(root)) then
    if (comm%mpirank == root) then
      call MPI_Reduce(mpi_in_place, buffer, size(buffer), mpi_integer, mpi_sum, root, comm%mpi_comm)
    else
      call MPI_Reduce(buffer, buffer, size(buffer), mpi_integer, mpi_sum, root, comm%mpi_comm)
    end if
  else
    call MPI_Allreduce(mpi_in_place, buffer, size(buffer), mpi_integer, mpi_sum, comm%mpi_comm)
  end if
#endif  
end subroutine mpi_sum_i1

subroutine mpi_sum_r_sp(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  real(sp), intent(inout)             :: buffer
  integer, optional, intent(in)       :: root
  
#ifdef MPI
  if (present(root)) then
    if (comm%mpirank == root) then
      call MPI_Reduce(mpi_in_place, buffer, 1, mpi_real, mpi_sum, root, comm%mpi_comm)
    else
      call MPI_Reduce(buffer, buffer, 1, mpi_real, mpi_sum, root, comm%mpi_comm)
    end if
  else
    call MPI_Allreduce(mpi_in_place, buffer, 1, mpi_real, mpi_sum, comm%mpi_comm)
  end if
#endif  
end subroutine mpi_sum_r_sp

subroutine mpi_sum_r1_sp(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  real(sp), intent(inout)             :: buffer(:)
  integer, optional, intent(in)       :: root
  
#ifdef MPI
  if (present(root)) then
    if (comm%mpirank == root) then
      call MPI_Reduce(mpi_in_place, buffer, size(buffer), mpi_real, mpi_sum, root, comm%mpi_comm)
    else
      call MPI_Reduce(buffer, buffer, size(buffer), mpi_real, mpi_sum, root, comm%mpi_comm)
    end if
  else
    call MPI_Allreduce(mpi_in_place, buffer, size(buffer), mpi_real, mpi_sum, comm%mpi_comm)
  end if
#endif  
end subroutine mpi_sum_r1_sp

subroutine mpi_sum_r2_sp(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  real(sp), intent(inout)             :: buffer(:,:)
  integer, optional, intent(in)       :: root
  
#ifdef MPI
  if (present(root)) then
    if (comm%mpirank == root) then
      call MPI_Reduce(mpi_in_place, buffer, size(buffer), mpi_real, mpi_sum, root, comm%mpi_comm)
    else
      call MPI_Reduce(buffer, buffer, size(buffer), mpi_real, mpi_sum, root, comm%mpi_comm)
    end if
  else
    call MPI_Allreduce(mpi_in_place, buffer, size(buffer), mpi_real, mpi_sum, comm%mpi_comm)
  end if
#endif  
end subroutine mpi_sum_r2_sp

subroutine mpi_sum_c_sp(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  complex(sp), intent(inout)          :: buffer
  integer, optional, intent(in)       :: root
  
#ifdef MPI
  if (present(root)) then
    if (comm%mpirank == root) then
      call MPI_Reduce(mpi_in_place, buffer, 1, mpi_complex, mpi_sum, root, comm%mpi_comm)
    else
      call MPI_Reduce(buffer, buffer, 1, mpi_complex, mpi_sum, root, comm%mpi_comm)
    end if
  else
    call MPI_Allreduce(mpi_in_place, buffer, 1, mpi_complex, mpi_sum, comm%mpi_comm)
  end if
#endif  
end subroutine mpi_sum_c_sp

subroutine mpi_sum_c1_sp(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  complex(sp), intent(inout)          :: buffer(:)
  integer, optional, intent(in)       :: root
  
#ifdef MPI
  if (present(root)) then
    if (comm%mpirank == root) then
      call MPI_Reduce(mpi_in_place, buffer, size(buffer), mpi_complex, mpi_sum, root, comm%mpi_comm)
    else
      call MPI_Reduce(buffer, buffer, size(buffer), mpi_complex, mpi_sum, root, comm%mpi_comm)
    end if
  else
    call MPI_Allreduce(mpi_in_place, buffer, size(buffer), mpi_complex, mpi_sum, comm%mpi_comm)
  end if
#endif  
end subroutine mpi_sum_c1_sp

subroutine mpi_sum_c2_sp(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  complex(sp), intent(inout)          :: buffer(:,:)
  integer, optional, intent(in)       :: root
  
#ifdef MPI
  if (present(root)) then
    if (comm%mpirank == root) then
      call MPI_Reduce(mpi_in_place, buffer, size(buffer), mpi_complex, mpi_sum, root, comm%mpi_comm)
    else
      call MPI_Reduce(buffer, buffer, size(buffer), mpi_complex, mpi_sum, root, comm%mpi_comm)
    end if
  else
    call MPI_Allreduce(mpi_in_place, buffer, size(buffer), mpi_complex, mpi_sum, comm%mpi_comm)
  end if
#endif  
end subroutine mpi_sum_c2_sp

subroutine mpi_sum_r_dp(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  real(dp), intent(inout)             :: buffer
  integer, optional, intent(in)       :: root
  
#ifdef MPI
  if (present(root)) then
    if (comm%mpirank == root) then
      call MPI_Reduce(mpi_in_place, buffer, 1, mpi_double_precision, mpi_sum, root, comm%mpi_comm)
    else
      call MPI_Reduce(buffer, buffer, 1, mpi_double_precision, mpi_sum, root, comm%mpi_comm)
    end if
  else
    call MPI_Allreduce(mpi_in_place, buffer, 1, mpi_double_precision, mpi_sum, comm%mpi_comm)
  end if
#endif  
end subroutine mpi_sum_r_dp

subroutine mpi_sum_r1_dp(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  real(dp), intent(inout)             :: buffer(:)
  integer, optional, intent(in)       :: root
  
#ifdef MPI
  if (present(root)) then
    if (comm%mpirank == root) then
      call MPI_Reduce(mpi_in_place, buffer, size(buffer), mpi_double_precision, mpi_sum, root, comm%mpi_comm)
    else
      call MPI_Reduce(buffer, buffer, size(buffer), mpi_double_precision, mpi_sum, root, comm%mpi_comm)
    end if
  else
    call MPI_Allreduce(mpi_in_place, buffer, size(buffer), mpi_double_precision, mpi_sum, comm%mpi_comm)
  end if
#endif  
end subroutine mpi_sum_r1_dp

subroutine mpi_sum_r2_dp(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  real(dp), intent(inout)             :: buffer(:,:)
  integer, optional, intent(in)       :: root
  
#ifdef MPI
  if (present(root)) then
    if (comm%mpirank == root) then
      call MPI_Reduce(mpi_in_place, buffer, size(buffer), mpi_double_precision, mpi_sum, root, comm%mpi_comm)
    else
      call MPI_Reduce(buffer, buffer, size(buffer), mpi_double_precision, mpi_sum, root, comm%mpi_comm)
    end if
  else
    call MPI_Allreduce(mpi_in_place, buffer, size(buffer), mpi_double_precision, mpi_sum, comm%mpi_comm)
  end if
#endif  
end subroutine mpi_sum_r2_dp

subroutine mpi_sum_c_dp(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  complex(dp), intent(inout)          :: buffer
  integer, optional, intent(in)       :: root
  
#ifdef MPI
  if (present(root)) then
    if (comm%mpirank == root) then
      call MPI_Reduce(mpi_in_place, buffer, 1, mpi_double_complex, mpi_sum, root, comm%mpi_comm)
    else
      call MPI_Reduce(buffer, buffer, 1, mpi_double_complex, mpi_sum, root, comm%mpi_comm)
    end if
  else
    call MPI_Allreduce(mpi_in_place, buffer, 1, mpi_double_complex, mpi_sum, comm%mpi_comm)
  end if
#endif  
end subroutine mpi_sum_c_dp

subroutine mpi_sum_c1_dp(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  complex(dp), intent(inout)          :: buffer(:)
  integer, optional, intent(in)       :: root
  
#ifdef MPI
  if (present(root)) then
    if (comm%mpirank == root) then
      call MPI_Reduce(mpi_in_place, buffer, size(buffer), mpi_double_complex, mpi_sum, root, comm%mpi_comm)
    else
      call MPI_Reduce(buffer, buffer, size(buffer), mpi_double_complex, mpi_sum, root, comm%mpi_comm)
    end if
  else
    call MPI_Allreduce(mpi_in_place, buffer, size(buffer), mpi_double_complex, mpi_sum, comm%mpi_comm)
  end if
#endif  
end subroutine mpi_sum_c1_dp

subroutine mpi_sum_c2_dp(comm, buffer, root)
  class(mpi_communicator), intent(in) :: comm
  complex(dp), intent(inout)          :: buffer(:,:)
  integer, optional, intent(in)       :: root
  
#ifdef MPI
  if (present(root)) then
    if (comm%mpirank == root) then
      call MPI_Reduce(mpi_in_place, buffer, size(buffer), mpi_double_complex, mpi_sum, root, comm%mpi_comm)
    else
      call MPI_Reduce(buffer, buffer, size(buffer), mpi_double_complex, mpi_sum, root, comm%mpi_comm)
    end if
  else
    call MPI_Allreduce(mpi_in_place, buffer, size(buffer), mpi_double_complex, mpi_sum, comm%mpi_comm)
  end if
#endif  
end subroutine mpi_sum_c2_dp


!******************************************************************************** 
!
! Interface routines for MPI_Reduce MPI_Allreduce
!
!******************************************************************************** 
#ifdef MPI
subroutine reduce_i(comm, buffer, root, op)
  class(mpi_communicator), intent(in) :: comm
  integer, intent(inout)              :: buffer
  integer, optional, intent(in)       :: root
  type(mpi_op), optional, intent(in)  :: op
  !local
  type(mpi_op)                        :: op_
  
  op_ = mpi_sum
  if (present(op)) op_ = op
  
  if (present(root)) then
    if (comm%mpirank == root) then
      call MPI_Reduce(mpi_in_place, buffer, 1, mpi_integer, op_, root, comm%mpi_comm)
    else
      call MPI_Reduce(buffer, buffer, 1, mpi_integer, op_, root, comm%mpi_comm)
    end if
  else
    call MPI_Allreduce(mpi_in_place, buffer, 1, mpi_integer, op_, comm%mpi_comm)
  end if
end subroutine reduce_i
#endif

#ifdef MPI
subroutine reduce_log(comm, buffer, root, op)
  class(mpi_communicator), intent(in) :: comm
  logical, intent(inout)              :: buffer
  integer, optional, intent(in)       :: root
  type(mpi_op), optional, intent(in)  :: op
  !local
  type(mpi_op)                        :: op_
  
  op_ = mpi_sum
  if (present(op)) op_ = op
  
  if (present(root)) then
    if (comm%mpirank == root) then
      call MPI_Reduce(mpi_in_place, buffer, 1, MPI_Logical, op_, root, comm%mpi_comm)
    else
      call MPI_Reduce(buffer, buffer, 1, MPI_Logical, op_, root, comm%mpi_comm)
    end if
  else
    call MPI_Allreduce(mpi_in_place, buffer, 1, MPI_Logical, op_, comm%mpi_comm)
  end if
end subroutine reduce_log
#endif

#ifdef MPI
subroutine reduce_r_sp(comm, buffer, root, op)
  class(mpi_communicator), intent(in) :: comm
  real(sp), intent(inout)             :: buffer
  integer, optional, intent(in)       :: root
  type(mpi_op), optional, intent(in)  :: op
  !local
  type(mpi_op)                        :: op_
  
  op_ = mpi_sum
  if (present(op)) op_ = op
  
  if (present(root)) then
    if (comm%mpirank == root) then
      call MPI_Reduce(mpi_in_place, buffer, 1, mpi_real, op_, root, comm%mpi_comm)
    else
      call MPI_Reduce(buffer, buffer, 1, mpi_real, op_, root, comm%mpi_comm)
    end if
  else
    call MPI_Allreduce(mpi_in_place, buffer, 1, mpi_real, op_, comm%mpi_comm)
  end if
end subroutine reduce_r_sp
#endif

#ifdef MPI
subroutine reduce_r1_sp(comm, buffer, root, op)
  class(mpi_communicator), intent(in) :: comm
  real(sp), intent(inout)             :: buffer(:)
  integer, optional, intent(in)       :: root
  type(mpi_op), optional, intent(in)  :: op
  !local
  type(mpi_op)                        :: op_
  
  op_ = mpi_sum
  if (present(op)) op_ = op
  
  if (present(root)) then
    if (comm%mpirank == root) then
      call MPI_Reduce(mpi_in_place, buffer, size(buffer), mpi_real, op_, root, comm%mpi_comm)
    else
      call MPI_Reduce(buffer, buffer, size(buffer), mpi_real, op_, root, comm%mpi_comm)
    end if
  else
    call MPI_Allreduce(mpi_in_place, buffer, size(buffer), mpi_real, op_, comm%mpi_comm)
  end if
end subroutine reduce_r1_sp
#endif

#ifdef MPI
subroutine reduce_r2_sp(comm, buffer, root, op)
  class(mpi_communicator), intent(in) :: comm
  real(sp), intent(inout)             :: buffer(:,:)
  integer, optional, intent(in)       :: root
  type(mpi_op), optional, intent(in)  :: op
  !local
  type(mpi_op)                        :: op_
  
  op_ = mpi_sum
  if (present(op)) op_ = op
  
  if (present(root)) then
    if (comm%mpirank == root) then
      call MPI_Reduce(mpi_in_place, buffer, size(buffer), mpi_real, op_, root, comm%mpi_comm)
    else
      call MPI_Reduce(buffer, buffer, size(buffer), mpi_real, op_, root, comm%mpi_comm)
    end if
  else
    call MPI_Allreduce(mpi_in_place, buffer, size(buffer), mpi_real, op_, comm%mpi_comm)
  end if
end subroutine reduce_r2_sp
#endif

#ifdef MPI
subroutine reduce_c_sp(comm, buffer, root, op)
  class(mpi_communicator), intent(in) :: comm
  complex(sp), intent(inout)          :: buffer
  integer, optional, intent(in)       :: root
  type(mpi_op), optional, intent(in)  :: op
  !local
  type(mpi_op)                        :: op_
  
  op_ = mpi_sum
  if (present(op)) op_ = op
  
  if (present(root)) then
    if (comm%mpirank == root) then
      call MPI_Reduce(mpi_in_place, buffer, 1, mpi_complex, op_, root, comm%mpi_comm)
    else
      call MPI_Reduce(buffer, buffer, 1, mpi_complex, op_, root, comm%mpi_comm)
    end if
  else
    call MPI_Allreduce(mpi_in_place, buffer, 1, mpi_complex, op_, comm%mpi_comm)
  end if
end subroutine reduce_c_sp
#endif

#ifdef MPI
subroutine reduce_c1_sp(comm, buffer, root, op)
  class(mpi_communicator), intent(in) :: comm
  complex(sp), intent(inout)          :: buffer(:)
  integer, optional, intent(in)       :: root
  type(mpi_op), optional, intent(in)  :: op
  !local
  type(mpi_op)                        :: op_
  
  op_ = mpi_sum
  if (present(op)) op_ = op
  
  if (present(root)) then
    if (comm%mpirank == root) then
      call MPI_Reduce(mpi_in_place, buffer, size(buffer), mpi_complex, op_, root, comm%mpi_comm)
    else
      call MPI_Reduce(buffer, buffer, size(buffer), mpi_complex, op_, root, comm%mpi_comm)
    end if
  else
    call MPI_Allreduce(mpi_in_place, buffer, size(buffer), mpi_complex, op_, comm%mpi_comm)
  end if
end subroutine reduce_c1_sp
#endif

#ifdef MPI
subroutine reduce_c2_sp(comm, buffer, root, op)
  class(mpi_communicator), intent(in) :: comm
  complex(sp), intent(inout)          :: buffer(:,:)
  integer, optional, intent(in)       :: root
  type(mpi_op), optional, intent(in)  :: op
  !local
  type(mpi_op)                        :: op_
  
  op_ = mpi_sum
  if (present(op)) op_ = op
  
  if (present(root)) then
    if (comm%mpirank == root) then
      call MPI_Reduce(mpi_in_place, buffer, size(buffer), mpi_complex, op_, root, comm%mpi_comm)
    else
      call MPI_Reduce(buffer, buffer, size(buffer), mpi_complex, op_, root, comm%mpi_comm)
    end if
  else
    call MPI_Allreduce(mpi_in_place, buffer, size(buffer), mpi_complex, op_, comm%mpi_comm)
  end if
end subroutine reduce_c2_sp
#endif

#ifdef MPI
subroutine reduce_r_dp(comm, buffer, root, op)
  class(mpi_communicator), intent(in) :: comm
  real(dp), intent(inout)             :: buffer
  integer, optional, intent(in)       :: root
  type(mpi_op), optional, intent(in)  :: op
  !local
  type(mpi_op)                        :: op_
  
  op_ = mpi_sum
  if (present(op)) op_ = op
  
  if (present(root)) then
    if (comm%mpirank == root) then
      call MPI_Reduce(mpi_in_place, buffer, 1, mpi_double_precision, op_, root, comm%mpi_comm)
    else
      call MPI_Reduce(buffer, buffer, 1, mpi_double_precision, op_, root, comm%mpi_comm)
    end if
  else
    call MPI_Allreduce(mpi_in_place, buffer, 1, mpi_double_precision, op_, comm%mpi_comm)
  end if
end subroutine reduce_r_dp
#endif

#ifdef MPI
subroutine reduce_r1_dp(comm, buffer, root, op)
  class(mpi_communicator), intent(in) :: comm
  real(dp), intent(inout)             :: buffer(:)
  integer, optional, intent(in)       :: root
  type(mpi_op), optional, intent(in)  :: op
  !local
  type(mpi_op)                        :: op_
  
  op_ = mpi_sum
  if (present(op)) op_ = op
  
  if (present(root)) then
    if (comm%mpirank == root) then
      call MPI_Reduce(mpi_in_place, buffer, size(buffer), mpi_double_precision, op_, root, comm%mpi_comm)
    else
      call MPI_Reduce(buffer, buffer, size(buffer), mpi_double_precision, op_, root, comm%mpi_comm)
    end if
  else
    call MPI_Allreduce(mpi_in_place, buffer, size(buffer), mpi_double_precision, op_, comm%mpi_comm)
  end if
end subroutine reduce_r1_dp
#endif

#ifdef MPI
subroutine reduce_r2_dp(comm, buffer, root, op)
  class(mpi_communicator), intent(in) :: comm
  real(dp), intent(inout)             :: buffer(:,:)
  integer, optional, intent(in)       :: root
  type(mpi_op), optional, intent(in)  :: op
  !local
  type(mpi_op)                        :: op_
  
  op_ = mpi_sum
  if (present(op)) op_ = op
  
  if (present(root)) then
    if (comm%mpirank == root) then
      call MPI_Reduce(mpi_in_place, buffer, size(buffer), mpi_double_precision, op_, root, comm%mpi_comm)
    else
      call MPI_Reduce(buffer, buffer, size(buffer), mpi_double_precision, op_, root, comm%mpi_comm)
    end if
  else
    call MPI_Allreduce(mpi_in_place, buffer, size(buffer), mpi_double_precision, op_, comm%mpi_comm)
  end if
end subroutine reduce_r2_dp
#endif

#ifdef MPI
subroutine reduce_c_dp(comm, buffer, root, op)
  class(mpi_communicator), intent(in) :: comm
  complex(dp), intent(inout)          :: buffer
  integer, optional, intent(in)       :: root
  type(mpi_op), optional, intent(in)  :: op
  !local
  type(mpi_op)                        :: op_
  
  op_ = mpi_sum
  if (present(op)) op_ = op
  
  if (present(root)) then
    if (comm%mpirank == root) then
      call MPI_Reduce(mpi_in_place, buffer, 1, mpi_double_complex, op_, root, comm%mpi_comm)
    else
      call MPI_Reduce(buffer, buffer, 1, mpi_double_complex, op_, root, comm%mpi_comm)
    end if
  else
    call MPI_Allreduce(mpi_in_place, buffer, 1, mpi_double_complex, op_, comm%mpi_comm)
  end if
end subroutine reduce_c_dp
#endif

#ifdef MPI
subroutine reduce_c1_dp(comm, buffer, root, op)
  class(mpi_communicator), intent(in) :: comm
  complex(dp), intent(inout)          :: buffer(:)
  integer, optional, intent(in)       :: root
  type(mpi_op), optional, intent(in)  :: op
  !local
  type(mpi_op)                        :: op_
  
  op_ = mpi_sum
  if (present(op)) op_ = op
  
  if (present(root)) then
    if (comm%mpirank == root) then
      call MPI_Reduce(mpi_in_place, buffer, size(buffer), mpi_double_complex, op_, root, comm%mpi_comm)
    else
      call MPI_Reduce(buffer, buffer, size(buffer), mpi_double_complex, op_, root, comm%mpi_comm)
    end if
  else
    call MPI_Allreduce(mpi_in_place, buffer, size(buffer), mpi_double_complex, op_, comm%mpi_comm)
  end if
end subroutine reduce_c1_dp
#endif

#ifdef MPI
subroutine reduce_c2_dp(comm, buffer, root, op)
  class(mpi_communicator), intent(in) :: comm
  complex(dp), intent(inout)          :: buffer(:,:)
  integer, optional, intent(in)       :: root
  type(mpi_op), optional, intent(in)  :: op
  !local
  type(mpi_op)                        :: op_
  
  op_ = mpi_sum
  if (present(op)) op_ = op
  
  if (present(root)) then
    if (comm%mpirank == root) then
      call MPI_Reduce(mpi_in_place, buffer, size(buffer), mpi_double_complex, op_, root, comm%mpi_comm)
    else
      call MPI_Reduce(buffer, buffer, size(buffer), mpi_double_complex, op_, root, comm%mpi_comm)
    end if
  else
    call MPI_Allreduce(mpi_in_place, buffer, size(buffer), mpi_double_complex, op_, comm%mpi_comm)
  end if
end subroutine reduce_c2_dp
#endif

end module  mpi