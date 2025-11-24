! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module constants

#include "preproc.inc"

implicit none 

!******************************************************************************** 
! Precision of  floating point numbers and integers  
!******************************************************************************** 
integer, parameter  :: intkind  = 8                             !for strings in (QM)CI
integer, parameter  :: intkind2 = 4                             !for walkers in  QMCI
integer, parameter  :: i4       = intkind2                      !4 byte integer
integer, parameter  :: i8       = intkind                       !8 byte integer
integer, parameter  :: sp       = selected_real_kind(6,37)      !single precision
integer, parameter  :: dp       = selected_real_kind(15,307)    !double precision
!debug: qprec not accepted by nvidia compilers
!integer, parameter  :: qprec    = selected_real_kind(33,4931)  !quadruple precision
integer, parameter  :: wp       = dp   
integer, parameter  :: spprec   = dp                            !precision for the sparse matrices


!******************************************************************************** 
! Constants and parameters
!******************************************************************************** 
integer, parameter  :: i4size       = 32
integer, parameter  :: i8size       = 64
integer, parameter  :: unset        = -1000
real(sp), parameter :: zeror_sp     = 0.0_sp
real(dp), parameter :: zeror_dp     = 0.0_dp
real(sp), parameter :: oner_sp      = 1.0_sp
real(dp), parameter :: oner_dp      = 1.0_dp
real(wp), parameter :: pi           = acos(-1.0_wp)
complex(sp), parameter :: zeroc_sp  = (0.0_sp, 0.0_sp)
complex(dp), parameter :: zeroc_dp  = (0.0_dp, 0.0_dp)
complex(sp), parameter :: onec_sp   = (1.0_sp, 0.0_sp)
complex(dp), parameter :: onec_dp   = (1.0_dp, 0.0_dp)


!******************************************************************************** 
! Character sizes and some default characters
!******************************************************************************** 
integer, parameter  :: charlen      = 60
integer, parameter  :: linec        = 256
integer, parameter  :: longc        = 2500
integer, parameter  :: atomc        = 4
character(len=60)   :: starshort = "************************************************************"
character(len=75)   :: starlong  = "***************************************************************************"


!**************************************************
! Physical constants and conversion factors          
!**************************************************
real(wp), parameter :: h_eV = 27.211386_wp
real(wp), parameter :: eV_h = 0.0367493225_wp
real(wp), parameter :: angstrom_to_bohr = 1.8897259886_wp             ! R[B] = angstrom_to_bohr * R[A]


!**************************************************
! Control variables                     
!**************************************************
TYPE CONTROLL
  !Sysyem controll
  LOGICAL:: LDMAT = .FALSE.             !Calculation of FCI density matrices
  LOGICAL:: LDMET = .FALSE.             !Calculating energy based one DMET 
  LOGICAL:: LCI_VEC = .TRUE.            !Vectorized CI algorithm
  LOGICAL:: LWRITE_CI = .FALSE.         !Write CI vector to file
  LOGICAL:: LRANDOMIZE = .FALSE.        !Random seed for RNG
  LOGICAL,DIMENSION(:),ALLOCATABLE:: ORB_MASK
  LOGICAL:: LMASK = .FALSE.             !Determines whether ORB_MASK is used
  !FCI input parameters
  INTEGER::CI_MAXITER = 10              !Maximal number of iterations in FCI iterative diag procedure
  REAL(wp):: CI_EPS = 0.0001_wp         !FCI convergence criterion
  LOGICAL:: LSPIN_RES=.TRUE.            !True iif only one Sz sector should be calc
  INTEGER:: CI_MODE = 1                 !Mode of CI calculation
  INTEGER:: CI_DIAG_MODE = 2            !Diagonalization method in CI
  INTEGER:: NUM_EIGS = 1                !Number of desired CI eigenvalues
  INTEGER:: CI_STEPS                    !Number of CI steps performed
  INTEGER:: LDMET_OPT = 1               !Connected with FCI Density matrix embedding 
  !
  INTEGER                :: write_files   = 1           !Controls amount of informations written to the files
  !Spin polarization variables
  INTEGER                :: ispin         = 1           !1 for restricted spin and 2 for unrestricted spin
  INTEGER                :: rspin         = 2           !Spin multiplicity  
  INTEGER                :: ispin1        = 1           !Spin dimension of one-electron Hamiltonian matrix elements
  INTEGER                :: ispin2        = 1           !Spin dimension of two-electron Hamiltonian matrix elements
  INTEGER                :: spin          = 0           !S_z component in the system
  !Number of electrons variables
  INTEGER                :: nelect        = 0           !Total number of electrons
  INTEGER                :: nocc          = 0           !Number of occupied states        
  INTEGER                :: nelect_a      = 0           !Number of spin up electrons
  INTEGER                :: nelect_b      = 0           !Number of spin down electrons
  INTEGER                :: nel(2)        = 0           !Number of electrons per spin channel (nelect_a, nelect_b)
  INTEGER                :: nell(3)       = 0           !Start index for spin loops (0, nelect_a, nelect_a+nelect_b)
  !Dimensions of arrays
  INTEGER                :: n             = 0           !Number of basis functions
  INTEGER                :: nbtot         = 0           !Total number of bands - max(nbands_up, nbands_down)
  INTEGER                :: ngtot         = 0           !Total number of auxiliary fields
  INTEGER, ALLOCATABLE   :: ng(:)                       !Number of auxiliary fields per spin channel - dim to bie ispin2
  !FCIQMC input parameters
  INTEGER:: EQ_STEPS     = 10000        !Number of steps in S-mode FCIQMC
  INTEGER:: WARM_UP      = 0            !Warm up phase without sampling
  INTEGER:: NUM_MC_STEPS = 0            !Number of steps in Nw-mode FCIQMC
  INTEGER:: PRINT_STEP   = 1000         !Print control variables to file
  INTEGER:: A_per        = 1            !Update E_SHIFT eac A_per steps
  INTEGER:: W_start      = 1            !Initial number of walkers
  REAL(wp):: DAMP        = 0.01_wp      !Damp parameter for change of E_SHIFT
  REAL(wp):: MC_TSTEP    = 0.001_wp     !time step in FCIQMC algorithm
  REAL(wp):: f_c         = 1.0_wp       !Nw/Ndet

  CHARACTER(len=10):: UNITS = "h"       !Unit used in simulation
END TYPE CONTROLL


!**************************************************
!       Energy variables used in methods
!**************************************************
TYPE fci_ENERGY
  REAL(wp):: ENUC
  REAL(wp):: EHF
  REAL(wp):: EH
  REAL(wp):: EF
  REAL(wp):: EHF_TOT
  REAL(wp):: EMP2
  REAL(wp):: EMP2_TOT
  REAL(wp):: EMP2_CORR
  REAL(wp):: EFCI
  REAL(wp):: EFCI_TOT
  REAL(wp):: EFCI_CORR
  REAL(wp):: EAFQMC
  REAL(wp):: EAFQMC_TOT
  REAL(wp):: EAFQMC_CORR
  REAL(wp):: E_SHIFT = 0.0_wp
  REAL(wp),DIMENSION(:),ALLOCATABLE:: CI_EIG
END TYPE fci_ENERGY

type(fci_energy) :: en

!**************************************************
!      Symmetry informations about system
!**************************************************
TYPE SYMMETRY
  LOGICAL:: CLSH = .TRUE.
  INTEGER:: SYM = 0
  INTEGER:: SPIN = 0
END TYPE SYMMETRY


!**************************************************
!  Block Hamiltonian for exact CI diagonalization
!**************************************************
TYPE BLOCKING_HAM
  INTEGER:: dimen
  INTEGER:: dimen_a , dimen_b
  INTEGER:: SPIN
  REAL(wp),DIMENSION(:,:),ALLOCATABLE:: CI_EV
  REAL(wp),DIMENSION(:),ALLOCATABLE:: CI_EW
  INTEGER(KIND=intkind),DIMENSION(:),ALLOCATABLE:: a_str_arr
  INTEGER(KIND=intkind),DIMENSION(:),ALLOCATABLE:: b_str_arr
END TYPE BLOCKING_HAM


!**************************************************
!       Representation of a walker in QMCI        
!**************************************************
TYPE WALKER_c
  INTEGER(KIND=intkind2):: walk
  INTEGER:: a_adr
  INTEGER:: b_adr
  TYPE(WALKER_c),POINTER:: ptr
END TYPE WALKER_c


!**************************************************
!    Data structure for exapnsion space and for   
!          representation of CI_vectors           
!      if more than 1 eigenvalues are desired     
!**************************************************
TYPE :: EXP_SPACE
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: exp_vec
END TYPE EXP_SPACE


!**************************************************
!debug: obsolete file units
!**************************************************
INTEGER, PARAMETER :: out = 74450
INTEGER, PARAMETER :: vect = 700
INTEGER, PARAMETER :: ci_eigenval = 1700

!**************************************************
!    global values - probabilities to choose      
!             excitations in QMCI                
!**************************************************
REAL(wp) :: P_s             !prob for single exc.
REAL(wp) :: P_d             !prob for double exc.
REAL(wp) :: P_s_mode        !prob to choose singl a/b exc.
REAL(wp) :: P_I_a , P_I_b   !total prob for single a/b exc
REAL(wp) :: P_II_aa , P_II_ab , P_II_bb !taotal prob for double excitations
REAL(wp) :: P_d_mode_aa , P_d_mode_bb , P_d_mode_ab !probability to choose aa/bb/ab

!control variables from control file INPUT
TYPE(CONTROLL):: INFO

!Symmetry inforamtions
TYPE(SYMMETRY):: SYM

!**************************************************
!       Array for new spawned walkers if          
!         arrays are used - version 2             
!**************************************************
INTEGER,PARAMETER:: MAXWALK=1000000
TYPE(WALKER_c),DIMENSION(MAXWALK):: sp_array
INTEGER:: spawned=0

end module constants