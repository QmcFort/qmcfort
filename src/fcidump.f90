! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module fcidump_io

#include "preproc.inc"
use constants
use file_handle, only: FileHandle
use method_base
use mpi
use profiling
use hamilton_vars
use hamilton_type

implicit none

contains

!******************************************************************************** 
!
!       Creates FCIDUMP file format
!
!******************************************************************************** 
subroutine fcidump(ham)
  type(Hamilton), intent(inout) :: ham
  !local
  integer                       :: p, q, r, s
  integer                       :: p_, q_, r_, s_
  integer                       :: spin, spin1, spin2
  real(wp)                      :: e0
  type(FileHandle)              :: fh
  type(QmcFortMethod)           :: method
  real(wp), parameter           :: tol=1.0e-08_wp
  character(len=*), parameter   :: proc_name = "fcidump"
  character(len=*), parameter   :: fname = "FCIDUMP"

  method = QmcFortMethod(method="fcidump", def_active=.false., basis="mo", integral_mode="eri")
  if (.not. method%is_active()) return
  call ham%meet_method_req(method)

  if (use_profiling) call start_profiling(proc_name)

  if (comm_world%mpirank /= 0) then
    if (use_profiling) call end_profiling(proc_name)
    return
  end if

  fh = FileHandle(fname)
  call fh%open(status="replace", action="write")

  call print_fcidump(ham%des)

  call fcidump_header(fh, ham%des%ispin*ham%des%n, ham%des%nel)

  !write ERIs
  do spin1 = 1, ham%des%ispin
    do spin2 = 1, ham%des%ispin
      do s = 1, ham%des%n
        do r = 1, ham%des%n
          do q = 1, ham%des%n
            do p = 1, ham%des%n
              spin =  ham%des%get_spin2(spin1, spin2)
              if (abs(ham%h2(p,q,r,s,spin)) < tol) cycle
              p_ = (spin1-1)*ham%des%n + p
              q_ = (spin1-1)*ham%des%n + q
              r_ = (spin2-1)*ham%des%n + r
              s_ = (spin2-1)*ham%des%n + s
              write(fh%funit, 100) ham%h2(p,q,r,s,spin), p_, q_, r_, s_
            end do
          end do
        end do
      end do
    end do
  end do

  !write one-core hamiltonian
  do spin1 = 1, ham%des%ispin
    do p = 1, ham%des%n
      do q = 1, ham%des%n
        spin =  ham%des%get_spin1(spin1)
        if (abs(ham%h1(p,q,spin)) < tol) cycle
        p_ = (spin1-1)*ham%des%n + p
        q_ = (spin1-1)*ham%des%n + q
        write (fh%funit, 100) ham%h1(p,q,spin), p_, q_, 0, 0
      end do
    end do
  end do

  !write constant energy
  e0 = ham%enuc + ham%h0
  write(fh%funit, 100) e0, 0, 0, 0, 0

  100 format (f18.10, 4i6)

  if (use_profiling) call end_profiling(proc_name)
end subroutine fcidump


!******************************************************************************** 
!
!       Write header of the FCIDUMP file
!
!******************************************************************************** 
subroutine fcidump_header(fh, n, nel)
  type(FileHandle), intent(in)  :: fh
  integer, intent(in)           :: n
  integer, intent(in)           :: nel(:)
  !local variables
  integer                       :: nelect, spin

  nelect = sum(nel)
  spin = nel(1) - nel(2)

  write(fh%funit,100) "&FCI  NORB=", n, ", NELEC=", nelect, ", MS=", spin, ","
  write(fh%funit,*) "ORBSYM=", symm_string(n)
  write(fh%funit,*) "ISYM=1,"
  write(fh%funit,*) "&END"
  100 format (a,i6,a,i6,a,i6,a)
end subroutine fcidump_header


!******************************************************************************** 
!
!       Small function to create symmetry string in FCIDUMP file
!
!******************************************************************************** 
function symm_string(n) result(symm_str)
  integer, intent(in)           :: n
  character(len=:), allocatable :: symm_str
  !local variables
  integer                       :: i

  symm_str = ""

  do i = 1, n
    symm_str = symm_str // "1,"
  end do
end function symm_string


!******************************************************************************** 
!
!       Print header for the fcidump
!
!******************************************************************************** 
subroutine print_fcidump(ham_des)
  type(hamil_descriptor), intent(in) :: ham_des

  if (comm_world%mpirank == 0) then
    write(io%screen%funit,*) "FCIDUMP"
    write(io%screen%funit,*) starshort
    write(io%screen%funit,100) "spin polarization      =", ham_des%ispin
    write(io%screen%funit,100) "number of orbitals     =", ham_des%n
    write(io%screen%funit,100) "number of up electrons =", ham_des%nel(1)
    write(io%screen%funit,100) "number of dn electrons =", ham_des%nel(2)
    write(io%screen%funit,*) starlong

    write(io%qmcfort_out%funit,*) "FCIDUMP"
    write(io%qmcfort_out%funit,*) starshort
    write(io%qmcfort_out%funit,100) "spin polarization      =", ham_des%ispin
    write(io%qmcfort_out%funit,100) "number of orbitals     =", ham_des%n
    write(io%qmcfort_out%funit,100) "number of up electrons =", ham_des%nel(1)
    write(io%qmcfort_out%funit,100) "number of dn electrons =", ham_des%nel(2)
    write(io%qmcfort_out%funit,*) starlong
  end if

  100 format (1x,a,i6)
end subroutine print_fcidump

end module fcidump_io