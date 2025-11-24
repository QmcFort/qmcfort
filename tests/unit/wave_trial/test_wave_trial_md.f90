! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module test_wave_trial_md

use fruit
use file_handle, only: FileHandle
!use wave_trial_md

implicit none

contains

!!!subroutine test_wave_trial_md_reader
!!!  type(wf_trial_md) :: phi_t
!!!  type(file_handle) :: fh
!!!
!!!  !phi_t%hdes%n = 14
!!!  !phi_t%hdes%nel = [7, 5]
!!!
!!!  fh = file_handle("cas_space", "txt")
!!!  call fh%open_(action_="write", status_="replace")
!!!  write(fh%funit,*) "# of dets   # of active orbitals   # of active electrons"
!!!  write(fh%funit,*) "5     4      2  2"
!!!  write(fh%funit,*) "alpha string          beta string            ci_coeff"
!!!  write(fh%funit,*) "[ 0 1 ]  [ 0 1 ]  0.9986" 
!!!  write(fh%funit,*) "[ 0 2 ]  [ 0 2 ]  0.0829"
!!!  write(fh%funit,*) "[ 0 3 ]  [ 0 3 ]  0.0223"
!!!  write(fh%funit,*) "[ 1 2 ]  [ 1 2 ]  0.1086"
!!!  write(fh%funit,*) "[ 0 1 ]  [ 0 2 ]  0.1024"
!!!  call fh%close_
!!!
!!!  call phi_t%read_cas_space()  
!!!end subroutine test_wave_trial_md_reader

end module test_wave_trial_md