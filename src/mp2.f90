! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module mp2

#include "preproc.inc"
use constants
use mpi, only: comm_world
use profiling
use standalone
use qmcfort_pos
use method_base
use hamilton_vars, only: hamil_descriptor
use hamilton_type, only: Hamilton
use energy_types, only: Energy
use orth, only: get_integral

implicit none 

private
public do_mp2

contains

subroutine do_mp2(ham, struc)
  use hfproc
  type(hamilton), intent(inout) :: ham
  type(structure), intent(in)   :: struc
  !local variables
  character(len=*), parameter   :: proc_name = "do_mp2"
  type(Energy)                  :: e
  type(QmcFortMethod)           :: method
  
  method = QmcFortMethod(method="mp2", def_active=.false., basis="mo", integral_mode="any")
  if (.not. method%is_active()) return
  call ham%meet_method_req(method)

  if (use_profiling) call start_profiling(proc_name)

  call print_mp2_header(ham%des)

  e%enuc = struc%calc_enuc()

  call calculate_gse_full(ham, e)

  call mp2_energy(ham, e)

  call print_mp2_footer(e)

  if (use_profiling) call end_profiling(proc_name)
end subroutine do_mp2


!******************************************************************************** 
!
!       Prints header for the MP2 calculation
!
!******************************************************************************** 
subroutine print_mp2_header(ham_des)
  use qmcfort_io
  type(hamil_descriptor), intent(in) :: ham_des
  
  if (comm_world%mpirank == 0) then
    write(io%screen%funit,*) "MP2 CALCULATION"
    write(io%screen%funit,*) starshort
    
    write(io%qmcfort_log%funit,*) "MP2 CALCULATION"
    write(io%qmcfort_log%funit,*) starshort
    
    write(io%qmcfort_out%funit,*) "MP2 CALCULATION"
    write(io%qmcfort_out%funit,*) starshort
    write(io%qmcfort_out%funit,100) "Number of ocupied orbitals per spin = ", ham_des%nel
    write(io%qmcfort_out%funit,100) "Number of virtual orbitals per spin = ", ham_des%n - ham_des%nel
  end if
  100 format(1x,a,2i6)
end subroutine print_mp2_header


!******************************************************************************** 
!
!       Prints header for the MP2 calculation
!
!******************************************************************************** 
subroutine print_mp2_footer(e)
  use qmcfort_io
  type(Energy), intent(in)           :: e
  
  if (comm_world%mpirank == 0) then
    call e%report(io%screen, "mp2")
    write(io%screen%funit,*) starlong
   
    call e%report(io%qmcfort_log, "mp2")
    write(io%qmcfort_log%funit,*) starlong
   
    write(io%qmcfort_out%funit,*)
    call e%report(io%qmcfort_out, "mp2")
    write(io%qmcfort_out%funit,*) starlong
  end if
end subroutine print_mp2_footer


!******************************************************************************** 
!
!       Calculates mp2 energy
!
!******************************************************************************** 
subroutine mp2_energy(ham, e)
  type(hamilton), intent(in)  :: ham
  type(Energy), intent(inout) :: e
  !local variables
  integer                     :: spin1, spin2, spin1_f, spin2_f, sp
  integer                     :: i, j, a, b, ind1, ind2, ind3, ind4
  real(wp)                    :: deig, integral1, integral2
  character(len=*), parameter :: proc_name = "mp2_energy"
  
  if (use_profiling) call start_profiling(proc_name)

  do spin1_f = 1, ham%des%ispin_fock
    spin1 = spin1_f
    if (ham%des%ispin == 1) spin1  = 1
    do spin2_f = 1, ham%des%ispin_fock
      spin2 = spin2_f
      if (ham%des%ispin == 1) spin2  = 1
      sp = ham%des%get_spin2(spin1, spin2)
      do b = ham%des%nel(spin2_f)+1, ham%des%n
        do j = 1, ham%des%nel(spin2_f)
          do a = ham%des%nel(spin1_f)+1, ham%des%n
            do i = 1, ham%des%nel(spin1_f)
              integral1 = 0.0_wp
              integral2 = 0.0_wp
              if (trim(ham%des%integral_mode) == "eri") then
                if (trim(ham%des%basis) == "mo") then
                  integral1 = ham%h2(i,a,j,b,sp)
                  if (spin1_f == spin2_f) integral2 = ham%h2(i,b,j,a,sp)
                else
                  integral1 = get_integral(ham%h2, ham%phi(:,:,spin1), ham%phi(:,:,spin2), i, a, j, b)
                  if (spin1_f == spin2_f) integral2 = get_integral(ham%h2, ham%phi(:,:,spin1), ham%phi(:,:,spin2), i, b, j, a)
                end if
              else
                  call get_index(ham%des%n, ham%des%n, i, a, ind1, 1)
                  call get_index(ham%des%n, ham%des%n, j, b, ind2, 1)
                  call get_index(ham%des%n, ham%des%n, i, b, ind3, 1)
                  call get_index(ham%des%n, ham%des%n, j, a, ind4, 1)
                if (trim(ham%des%basis) == "mo") then
                  integral1 = get_integral(ham%h2_gres, ind1, ind2, spin1, spin2)
                  if (spin1_f == spin2_f) integral2 = get_integral(ham%h2_gres, ind3, ind4, spin1, spin2)
                else
                  integral1 = get_integral(ham%h2_gres, ham%phi(:,:,spin1), ham%phi(:,:,spin2), ind1, ind2)
                  if (spin1_f == spin2_f) integral2 = get_integral(ham%h2_gres, ham%phi(:,:,spin1), ham%phi(:,:,spin2), ind3, ind4)
                end if
              end if
              
              deig = ham%eigenval(i,spin1) + ham%eigenval(j,spin2) - ham%eigenval(a,spin1) - ham%eigenval(b,spin2)
              e%eh = e%eh + 0.5_wp * abs(integral1)**2 * ham%des%rspin**2 / deig
              e%ex = e%ex - 0.5_wp * integral1 * integral2 * ham%des%rspin / deig
            end do
          end do
        end do
      end do
    end do
  end do

  if (use_profiling) call end_profiling(proc_name)
end subroutine mp2_energy

end module mp2
