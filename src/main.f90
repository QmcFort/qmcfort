! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

program main

#include "preproc.inc"

use constants
use mpi
use log_manager_mod
use profiling
use standalone
use lapack
use method_base, only: QmcFortMethod
use qmcfort_io
use qmcfort_pos
use brng_factories, only: brng_glb, create_qmcfort_brng
use drng_factories, only: drng_glb, create_drng
use hamilton_type, only: Hamilton, read_files
use frozen_core_approx, only: frozen_core
use hfproc, only: do_hartree_fock
use afqmc, only: do_afqmc
use noci_afqmc_selection, only: do_noci_afqmc_selection
use mp2
use ci
use fciqmc
use fci_wrap
use fcidump_io

implicit none

!********************************************************************************
! Data declaration    
!********************************************************************************
type(Structure)         :: struc                !Geometry of the molecule
type(Hamilton), target  :: ham                  !Hamiltonain definition

type(QmcFortMethod)     :: fci_method

!Symmetry group informations
!TYPE(SYMM_GROUP):: SYMMS
!INTEGER:: NUM_SEA
!TYPE(SEA_SET),DIMENSION(:),ALLOCATABLE:: SEA
!local variables
integer                :: new_orb_num
logical                :: retour

!******************************************************************************** 
! Initialize MPI
!******************************************************************************** 
call comm_world%init

!******************************************************************************** 
! Initialize timers for PROFILING if defined
!******************************************************************************** 
call initialize_timers()

!********************************************************************************
! Start program and touch the files which will be used later   
!********************************************************************************
call init_io()

!********************************************************************************
! Create logging object and initialize qmcfort_log file  
!********************************************************************************
logging = LogManager([io%qmcfort_log, io%screen])

!********************************************************************************
! Read Molecular geometry if present
!********************************************************************************
call read_qmcfort_structure(struc)

!********************************************************************************
! Initialize a global random number generator
!********************************************************************************
call create_qmcfort_brng(brng_glb)
call create_drng(brng_glb, drng_glb)

!******************************************************************************** 
! Main File reader: qmcfort_in basis_set ham files
!********************************************************************************  
call read_files(struc, info, ham, sym)

!******************************************************************************** 
! Point group exploration of the system
!******************************************************************************** 
!call symm_mol(struc, symms)

!********************************************************************************
! Hartree-Fock calculation   
!********************************************************************************
call do_hartree_fock(ham, struc)

!********************************************************************************
! Change initial Gaussian basis to some orthogonal basis
!********************************************************************************
call ham%transform_basis()

!******************************************************************************** 
! Here comes the frozen core approximation
!******************************************************************************** 
call frozen_core(struc, ham, add_pot=.true.)

!******************************************************************************** 
! Write Hamiltonian files
!******************************************************************************** 
call ham%write_()

!******************************************************************************** 
! Create FCIDUMP file 
!******************************************************************************** 
!call fcidump(ham)

!********************************************************************************
! Calculate MP2 correlation energy
!********************************************************************************
call do_mp2(ham, struc)

!********************************************************************************
! START CI       
!********************************************************************************
fci_method = QmcFortMethod(method="fci", def_active=.false., basis="mo", integral_mode="eri")
if (fci_method%is_active()) then
  call ham%meet_method_req(fci_method)
  call do_fci(info, struc, ham%des, sym, ham, en, ham%h1, ham%h2)
end if

!******************************************************************************** 
! AFQMC Calculation
!******************************************************************************** 
call do_afqmc(ham, struc)

!******************************************************************************** 
! NOCI Selection based on the AFQMC walk
!******************************************************************************** 
call do_noci_afqmc_selection(ham, struc)

!********************************************************************************
! Write time profiling to file TIMING        
!********************************************************************************
call finalize_timers

!********************************************************************************
! Dump footers and close output files
!********************************************************************************
call finalize_io()

!******************************************************************************** 
! Finalize MPI
!******************************************************************************** 
call comm_world%finalize

end program main