! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

MODULE FCI_WRAP
!********************************************************************************
!
!  Collects all (F)CI routines and contains a wrap subroutine to call one of 
!  those routines. Additionally, this file contains routines for calculation
!          of one-electron and two-electron reduced density matrices         
!
!********************************************************************************

#include "preproc.inc"
USE CONSTANTS
use mpi, only: comm_world
use qmcfort_pos
USE STANDALONE
USE LAPACK
USE ORTH
USE CI
USE FCIQMC
use hamilton_type, only: Hamilton
use hamilton_vars, only: hamil_descriptor

implicit none 

INTERFACE CI_ASSEMBLER
  MODULE PROCEDURE CI_ASSEMBLER_str
  MODULE PROCEDURE CI_ASSEMBLER_det
END INTERFACE CI_ASSEMBLER

INTERFACE WRITE_BLOCK_EW
  MODULE PROCEDURE WRITE_BLOCK_EW_arr
  MODULE PROCEDURE WRITE_BLOCK_EW_one
END INTERFACE WRITE_BLOCK_EW

CONTAINS


SUBROUTINE CI_ASSEMBLER_str(ORB_NUM,NELECT_a,NELECT_b,dimen_a,dimen_b, &
                        W_wa,W_wb,Y_wa,Y_wb,a_str_arr,b_str_arr)
  INTEGER,INTENT(IN):: ORB_NUM, NELECT_a, NELECT_b
  INTEGER,INTENT(OUT):: dimen_a, dimen_b
  INTEGER(KIND=intkind),DIMENSION(:),ALLOCATABLE,INTENT(OUT) :: a_str_arr , b_str_arr
  INTEGER,DIMENSION(0:ORB_NUM,0:NELECT_a),INTENT(OUT):: W_wa
  INTEGER,DIMENSION(0:ORB_NUM,0:NELECT_b),INTENT(OUT):: W_wb
  INTEGER,DIMENSION(0:1,ORB_NUM,NELECT_a),INTENT(OUT):: Y_wa
  INTEGER,DIMENSION(0:1,ORB_NUM,NELECT_b),INTENT(OUT):: Y_wb

  !Graphical representation for alpha stings
  IF(NELECT_a==0) THEN 
    W_wa = 1
    Y_wa = 0
  
    dimen_a = 1
 
    IF(.not. ALLOCATED(a_str_arr)) THEN
      ALLOCATE(a_str_arr(1))
    ENDIF

    a_str_arr(1) = 0
  ELSE
    CALL CALC_VERTEX_WEIGHTS(W_wa,ORB_NUM,NELECT_a)
    CALL CALC_ARC_WEIGHTS(Y_wa,W_wa,ORB_NUM,NELECT_a)
    CALL STRING_ASSEMBLER(ORB_NUM,NELECT_a,dimen_a,Y_wa,a_str_arr) !within routines
  ENDIF

  !Graphical representation for beta stings
  IF(NELECT_b==0) THEN
    W_wb = 1
    Y_wb = 0
   
    dimen_b = 1

    IF(.not. ALLOCATED(b_str_arr)) THEN
      ALLOCATE(b_str_arr(1))
    ENDIF

    b_str_arr(1) = 0
  ELSE
    CALL CALC_VERTEX_WEIGHTS(W_wb,ORB_NUM,NELECT_b)
    CALL CALC_ARC_WEIGHTS(Y_wb,W_wb,ORB_NUM,NELECT_b)
    CALL STRING_ASSEMBLER(ORB_NUM,NELECT_b,dimen_b,Y_wb,b_str_arr) !String list allocated
  ENDIF

  IF(comm_world%mpirank==0) THEN
    write(*,*)"Number of a strings = " , dimen_a
    write(*,*)"Number of b strings = " , dimen_b
  ENDIF

END SUBROUTINE CI_ASSEMBLER_str


SUBROUTINE CI_ASSEMBLER_det(ORB_NUM,NELECT,dimen,W_w,Y_w,str_arr)
  INTEGER,INTENT(IN):: ORB_NUM, NELECT
  INTEGER,INTENT(OUT):: dimen
  INTEGER(KIND=intkind),DIMENSION(:),ALLOCATABLE,INTENT(OUT) :: str_arr 
  INTEGER,DIMENSION(0:ORB_NUM,0:NELECT),INTENT(OUT):: W_w
  INTEGER,DIMENSION(0:1,ORB_NUM,NELECT),INTENT(OUT):: Y_w

  IF(NELECT==0) THEN
    W_w = 1 
    Y_w = 0

    dimen = 1

    IF(.not. ALLOCATED(str_arr)) THEN
      ALLOCATE(str_arr(1))
    ENDIF 

    str_arr(1) = 0
  ELSE  
    CALL CALC_VERTEX_WEIGHTS(W_w,ORB_NUM,NELECT)
    CALL CALC_ARC_WEIGHTS(Y_w,W_w,ORB_NUM,NELECT)
    CALL STRING_ASSEMBLER(ORB_NUM,NELECT,dimen,Y_w,str_arr) !String list allocated
  ENDIF

  IF(comm_world%mpirank==0) THEN
    WRITE(*,*) "Number of determinants = ", dimen
  ENDIF
END SUBROUTINE CI_ASSEMBLER_det


SUBROUTINE DO_FCI(INFO,STRUC,hdes,SYM,ham,EN,HAM_1,INT_2)
  TYPE(CONTROLL),INTENT(INOUT):: INFO
  TYPE(STRUCTURE),INTENT(INOUT):: STRUC
  type(hamil_descriptor), intent(inout) :: hdes
  TYPE(SYMMETRY),INTENT(INOUT):: SYM
  type(Hamilton), intent(inout) :: ham
  TYPE(fci_ENERGY),INTENT(INOUT):: EN
  REAL(wp), INTENT(INOUT):: HAM_1(:,:,:)
  REAL(wp),INTENT(INOUT):: INT_2(:,:,:,:,:)

  IF(comm_world%mpirank==0) THEN
    WRITE(*,*) "CI CALCULATION"
  ENDIF

  CI_MODE: IF(INFO%CI_MODE==0) THEN
    IF(INFO%LSPIN_RES) THEN
      CALL DO_EX_FCI_SP(INFO,STRUC,hdes,HAM_1,INT_2)
    ELSE
      CALL DO_EX_FCI_FULL(INFO,STRUC,hdes,HAM_1,INT_2)
    ENDIF
  ELSEIF(INFO%CI_MODE==1) THEN
    CALL DO_DIR_FCI(INFO,STRUC,hdes,SYM,ham,EN,HAM_1,INT_2)
    
    EN%EFCI = EN%CI_EIG(1)
    EN%EFCI_TOT = EN%EFCI + EN%ENUC
    EN%EFCI_CORR = EN%EFCI - EN%EHF
  ELSEIF (INFO%CI_MODE==2) THEN
    CALL DO_FCIQMC(INFO,STRUC,hdes,SYM,EN,HAM_1,INT_2)
    
    EN%EFCI = EN%EHF + EN%EFCI_CORR
    EN%EFCI_TOT = EN%EFCI + EN%ENUC
  ENDIF CI_MODE

END SUBROUTINE DO_FCI


SUBROUTINE DO_EX_FCI_SP(INFO,STRUC,hdes,HAM_1,INT_2)
  TYPE(CONTROLL),INTENT(INOUT):: INFO
  TYPE(STRUCTURE),INTENT(INOUT):: STRUC
  type(hamil_descriptor), intent(inout) :: hdes
  REAL(wp),DIMENSION(hdes%n,hdes%n,INFO%ISPIN1),INTENT(INOUT):: HAM_1
  REAL(wp),DIMENSION(hdes%n,hdes%n,hdes%n,hdes%n,INFO%ISPIN2),INTENT(INOUT):: INT_2

  !local variables 
  TYPE(BLOCKING_HAM):: BLOCK_HAM

  CALL DO_EX_FCI(INFO,STRUC,hdes,HAM_1,INT_2,&
                      BLOCK_HAM%CI_EV,BLOCK_HAM%CI_EW,&
                      BLOCK_HAM%a_str_arr,BLOCK_HAM%b_str_arr)

  BLOCK_HAM%SPIN = hdes%nel(1) - hdes%nel(2)
  BLOCK_HAM%dimen = SIZE(BLOCK_HAM%CI_EW)
  BLOCK_HAM%dimen_a = SIZE(BLOCK_HAM%a_str_arr)
  BLOCK_HAM%dimen_b = SIZE(BLOCK_HAM%b_str_arr)

  IF(comm_world%mpirank==0) THEN
    CALL WRITE_BLOCK_EW(BLOCK_HAM)
  ENDIF

END SUBROUTINE DO_EX_FCI_SP


SUBROUTINE DO_EX_FCI_FULL(INFO,STRUC,hdes,HAM_1,INT_2)
  TYPE(CONTROLL),INTENT(INOUT):: INFO
  TYPE(STRUCTURE),INTENT(INOUT):: STRUC
  type(hamil_descriptor), intent(inout) :: hdes
  REAL(wp),DIMENSION(hdes%n,hdes%n,INFO%ISPIN1),INTENT(INOUT):: HAM_1
  REAL(wp),DIMENSION(hdes%n,hdes%n,hdes%n,hdes%n,INFO%ISPIN2),INTENT(INOUT):: INT_2

  !local variables
  TYPE(BLOCKING_HAM),DIMENSION(:),ALLOCATABLE:: BLOCK_HAM
  INTEGER:: SPIN, MAX_SPIN , MIN_SPIN , NSZ, i, nel

  MAX_SPIN = sum(hdes%nel)
  IF(mod(sum(hdes%nel),2)==0) THEN
    MIN_SPIN = 0
  ELSE
    MIN_SPIN = 1
  ENDIF

  NSZ = (MAX_SPIN - MIN_SPIN)/2 + 1

  IF(.not. ALLOCATED(BLOCK_HAM)) THEN
    ALLOCATE(BLOCK_HAM(NSZ))
  ENDIF

  SPIN = MAX_SPIN + 2

  DO i=1,NSZ
    SPIN = SPIN - 2

    nel = sum(hdes%nel)

    hdes%nel(1) = (nel + SPIN)/2
    hdes%nel(2) = nel - hdes%nel(1)

    IF(comm_world%mpirank==0) THEN
      WRITE(*,*) "SPIN = " , SPIN
    ENDIF

    CALL DO_EX_FCI(INFO,STRUC,hdes,HAM_1,INT_2,&
                      BLOCK_HAM(i)%CI_EV,BLOCK_HAM(i)%CI_EW,&
                      BLOCK_HAM(i)%a_str_arr,BLOCK_HAM(i)%b_str_arr)

    BLOCK_HAM(i)%SPIN = SPIN
    BLOCK_HAM(i)%dimen = SIZE(BLOCK_HAM(i)%CI_EW)
    BLOCK_HAM(i)%dimen_a = SIZE(BLOCK_HAM(i)%a_str_arr)
    BLOCK_HAM(i)%dimen_b = SIZE(BLOCK_HAM(i)%b_str_arr)
  ENDDO

  IF(comm_world%mpirank==0) THEN
    CALL WRITE_BLOCK_EW(NSZ,BLOCK_HAM)
  ENDIF

END SUBROUTINE DO_EX_FCI_FULL


SUBROUTINE DO_EX_FCI(INFO,STRUC,hdes,HAM_1,INT_2,CI_EV,CI_EW,a_str_arr,b_str_arr)
  TYPE(CONTROLL),INTENT(INOUT):: INFO
  TYPE(STRUCTURE),INTENT(INOUT):: STRUC
  type(hamil_descriptor), intent(inout) :: hdes
  REAL(wp),DIMENSION(hdes%n,hdes%n,INFO%ISPIN1),INTENT(INOUT):: HAM_1
  REAL(wp),DIMENSION(hdes%n,hdes%n,hdes%n,hdes%n,INFO%ISPIN2),INTENT(INOUT):: INT_2
  REAL(wp),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: CI_EV
  REAL(wp),DIMENSION(:),ALLOCATABLE,INTENT(OUT):: CI_EW
  INTEGER(KIND=intkind),DIMENSION(:),ALLOCATABLE,INTENT(OUT):: a_str_arr , b_str_arr

  !local variables
  INTEGER:: dimen_a , dimen_b , dimen
  INTEGER,DIMENSION(0:hdes%n,0:hdes%nel(1)):: W_wa
  INTEGER,DIMENSION(0:hdes%n,0:hdes%nel(2)):: W_wb
  INTEGER,DIMENSION(0:1,hdes%n,hdes%nel(1)):: Y_wa
  INTEGER,DIMENSION(0:1,hdes%n,hdes%nel(2)):: Y_wb
  REAL(wp),DIMENSION(:,:),ALLOCATABLE:: CI_MAT

  !Set Slater Determinants
  CALL CI_ASSEMBLER(hdes%n,hdes%nel(1),hdes%nel(2),dimen_a,dimen_b,W_wa,W_wb, &
                    Y_wa,Y_wb,a_str_arr,b_str_arr)

  dimen = dimen_a * dimen_b

  IF(.not. ALLOCATED(CI_MAT)) THEN
    ALLOCATE(CI_MAT(dimen,dimen))
  ENDIF
  IF(.not. ALLOCATED(CI_EW)) THEN
    ALLOCATE(CI_EW(dimen))
  ENDIF
  IF(.not. ALLOCATED(CI_EV)) THEN
    ALLOCATE(CI_EV(dimen,dimen))
  ENDIF
  
  !Calculate CI Matrix
  CALL CALC_CI_MATRIX(INFO,STRUC,hdes,dimen_a,dimen_b,a_str_arr,b_str_arr,Y_wa,Y_wb,HAM_1,INT_2,CI_MAT)

  !Diagonalize CI_MATRIX 
  CALL SYM_EIGV(dimen,CI_MAT,CI_EW,CI_EV) 

END SUBROUTINE DO_EX_FCI


SUBROUTINE DO_DIR_FCI(INFO,STRUC,hdes,SYM,ham,EN,HAM_1,INT_2)
  TYPE(CONTROLL),INTENT(INOUT):: INFO
  TYPE(STRUCTURE),INTENT(INOUT):: STRUC
  type(hamil_descriptor), intent(inout) :: hdes
  TYPE(SYMMETRY),INTENT(INOUT):: SYM
  type(Hamilton), intent(inout) :: ham
  TYPE(fci_ENERGY),INTENT(INOUT):: EN
  REAL(wp),DIMENSION(hdes%n,hdes%n,INFO%ISPIN1),INTENT(INOUT):: HAM_1
  REAL(wp),DIMENSION(hdes%n,hdes%n,hdes%n,hdes%n,INFO%ISPIN2),INTENT(INOUT):: INT_2
  
  !local variables
  INTEGER:: dimen_a,dimen_b
  INTEGER(KIND=intkind),DIMENSION(:),ALLOCATABLE:: a_str_arr , b_str_arr
  INTEGER,DIMENSION(0:hdes%n,0:hdes%nel(1)):: W_wa
  INTEGER,DIMENSION(0:hdes%n,0:hdes%nel(2)):: W_wb
  INTEGER,DIMENSION(0:1,hdes%n,hdes%nel(1)):: Y_wa
  INTEGER,DIMENSION(0:1,hdes%n,hdes%nel(2)):: Y_wb
  REAL(wp) :: ONE_DMAT(hdes%n, hdes%n, INFO%ISPIN)
  REAL(wp),DIMENSION(hdes%n,hdes%n,hdes%n,hdes%n):: TWO_DMAT
  REAL(wp),DIMENSION(hdes%n,hdes%n,hdes%n,hdes%n,INFO%ISPIN1):: HAM_1T
  REAL(wp),DIMENSION(hdes%n, info%ispin):: OCCUP
  REAL(wp),DIMENSION(hdes%n,hdes%n, info%ispin):: NO
  REAL(wp),DIMENSION(:,:),ALLOCATABLE:: CI_H_DIAG , PROJ
  TYPE(EXP_SPACE),DIMENSION(:),ALLOCATABLE:: CI_vector
  REAL(wp):: alpha=1.0_wp
  INTEGER:: istat  
  REAL(wp):: DMATEN

  CALL CI_ASSEMBLER(hdes%n,hdes%nel(1),hdes%nel(2),dimen_a,dimen_b, &
                    W_wa,W_wb,Y_wa,Y_wb,a_str_arr,b_str_arr)

  IF(.NOT. ALLOCATED(CI_H_DIAG)) THEN
    ALLOCATE(CI_H_DIAG(dimen_a,dimen_b),STAT=istat)
    IF(istat/=0) THEN
      IF(comm_world%mpirank==0) THEN
        WRITE(*,*) "Allocation of CI_H_DIAG failed"
      ENDIF
    ENDIF
  ENDIF
 
  CALL CALC_CI_DIAG_HAM(a_str_arr,b_str_arr,dimen_a,dimen_b,HAM_1,INT_2, &
                        hdes%n,hdes%nel(1),hdes%nel(2),Y_wa,Y_wb,SYM%SYM,CI_H_DIAG )
  WRITE(*,*) "Energy of HF det = ", CI_H_DIAG(1,1)

  IF(.NOT. ALLOCATED(EN%CI_EIG)) THEN
    ALLOCATE(EN%CI_EIG(INFO%NUM_EIGS),STAT=istat)
    IF(istat/=0) THEN
      IF(comm_world%mpirank==0) THEN
        WRITE(*,*) "Allocation of CI_EIG fails"
      ENDIF
      CALL EXIT
    ENDIF
  ENDIF
 
  !Initialize trial vectors for iterative diagonalization
  !and allocate CI_vector
  CALL INITIALIZE_CI_VECTOR(dimen_a,dimen_b,CI_H_DIAG,CI_vector,INFO%NUM_EIGS)

  !Transform one-electron integrals before direct CI calculation
  CALL TRAFO_HAM1_CI(hdes%n,HAM_1,INT_2,HAM_1T)

  NUM_EIGS: IF(INFO%NUM_EIGS>1) THEN
    CALL CI_DAVIDSON_LIU_LOOP(hdes%n,hdes%nel(1),hdes%nel(2),dimen_a,dimen_b,HAM_1T,INT_2, &
                Y_wa,Y_wb,a_str_arr,b_str_arr,CI_H_DIAG,INFO%NUM_EIGS,INFO%LCI_VEC,CI_vector,EN%CI_EIG, &
                INFO%CI_STEPS,INFO%CI_MAXITER,INFO%CI_EPS,SYM%SYM,SYM%SPIN)
  ELSE
    DIAG_MODE: IF(INFO%CI_DIAG_MODE==1) THEN
      CALL CI_STEEPEST_DESCENT_LOOP(hdes%n,hdes%nel(1),hdes%nel(2),dimen_a,dimen_b,HAM_1T,INT_2, &
                 Y_wa,Y_wb,a_str_arr,b_str_arr,CI_H_DIAG,INFO%LCI_VEC,CI_vector(1)%exp_vec,EN%CI_EIG(1),   &
                 INFO%CI_STEPS,INFO%CI_MAXITER,INFO%CI_EPS,SYM%SYM,SYM%SPIN,alpha )
    ELSEIF(INFO%CI_DIAG_MODE==2) THEN
      CALL CI_DAVIDSON_LOOP(hdes%n,hdes%nel(1),hdes%nel(2),dimen_a,dimen_b,HAM_1T,INT_2, &
                 Y_wa,Y_wb,a_str_arr,b_str_arr,CI_H_DIAG,INFO%LCI_VEC,CI_vector(1)%exp_vec,EN%CI_EIG,INFO%CI_STEPS,&
                 INFO%CI_MAXITER,INFO%CI_EPS,SYM%SYM,SYM%SPIN)

    ELSEIF(INFO%CI_DIAG_MODE==3) THEN
      CALL CI_DAVIDSON_LIU_LOOP(hdes%n,hdes%nel(1),hdes%nel(2),dimen_a,dimen_b,HAM_1T,INT_2,&
                  Y_wa,Y_wb,a_str_arr,b_str_arr,CI_H_DIAG,INFO%NUM_EIGS,INFO%LCI_VEC,CI_vector,EN%CI_EIG,INFO%CI_STEPS, &
                  INFO%CI_MAXITER,INFO%CI_EPS,SYM%SYM,SYM%SPIN)
    ENDIF DIAG_MODE
  ENDIF NUM_EIGS

write(*,*) "ci davidson step done   "

  IF(comm_world%mpirank==0) THEN
    WRITE_CI: IF(INFO%LWRITE_CI) THEN
      CALL WRITE_CI_VECTOR(CI_vector,dimen_a,dimen_b,INFO%NUM_EIGS)
    ENDIF WRITE_CI
  ENDIF

  IF(comm_world%mpirank==0) THEN
    WRITE(*,*) "Weight of the HF determinant in FCI WF = ", CI_vector(1)%exp_vec(1,1)
  ENDIF

  !At this moment, the CI_vector is calculated 
  LDMAT: IF(INFO%LDMAT)  THEN
    CALL CALC_ONE_DMAT(STRUC,hdes,SYM,dimen_a,dimen_b,a_str_arr,b_str_arr,Y_wa,Y_wb,CI_vector(1)%exp_vec,ONE_DMAT)
    CALL CALC_TWO_DMAT(STRUC,hdes,SYM,dimen_a,dimen_b,a_str_arr,b_str_arr,Y_wa,Y_wb,CI_vector(1)%exp_vec,ONE_DMAT,TWO_DMAT)
    
    LDMET: IF(INFO%LDMET) THEN
      !First read projection matrix from whole canonical basis to local basis
      ALLOCATE(PROJ(1,hdes%n),STAT=istat)
      IF(istat/=0) THEN
        IF(comm_world%mpirank==0) THEN
          WRITE(*,*) "Allocation of PROJ fails"
        ENDIF
        CALL EXIT
      ENDIF
      block
        use array_numpy_file_io, only: ArrayNumpyFileIO
        type(ArrayNumpyFileIO) :: numpy_io
        call numpy_io%load(ham%h_io%orbitals%fname, proj)
      end block
      IF(INFO%LDMET_opt==1) THEN
        CALL GEN_TRAFO_INT_1(hdes%n,1,PROJ,HAM_1,INT_2)
        CALL GEN_TRAFO_INT_1(hdes%n,1,PROJ,ONE_DMAT(:,:,1),TWO_DMAT)
        CALL CALC_EN_DMAT_EMB(hdes%n,1,HAM_1,INT_2,ONE_DMAT(:,:,1),TWO_DMAT,DMATEN)
      ELSE
        CALL CALC_EN_DMAT_EMB2(hdes%n,1,PROJ,HAM_1,INT_2,ONE_DMAT(:,:,1),TWO_DMAT,DMATEN)
      ENDIF
    ELSE
      CALL CALC_EN_DMAT(hdes%n,HAM_1,INT_2,ONE_DMAT,TWO_DMAT,DMATEN)    
    ENDIF LDMET

    IF(comm_world%mpirank==0) THEN
      WRITE(*,*) "Energy from density matrices= ", DMATEN
    ENDIF

    CALL CALC_NO(ham, hdes%n,ONE_DMAT,OCCUP,NO)

  ENDIF LDMAT

END SUBROUTINE DO_DIR_FCI


SUBROUTINE DO_FCIQMC(INFO,STRUC,hdes,SYM,EN,HAM_1,INT_2)
  TYPE(CONTROLL),INTENT(INOUT):: INFO
  TYPE(STRUCTURE),INTENT(INOUT):: STRUC
  type(hamil_descriptor), intent(inout) :: hdes
  TYPE(SYMMETRY),INTENT(INOUT):: SYM
  TYPE(fci_ENERGY),INTENT(INOUT):: EN
  REAL(wp),DIMENSION(hdes%n,hdes%n,INFO%ISPIN1),INTENT(INOUT):: HAM_1
  REAL(wp),DIMENSION(hdes%n,hdes%n,hdes%n,hdes%n,INFO%ISPIN2),INTENT(INOUT):: INT_2

  !local variables
  INTEGER:: dimen_a,dimen_b
  INTEGER(KIND=intkind),DIMENSION(:),ALLOCATABLE:: a_str_arr , b_str_arr
  INTEGER,DIMENSION(0:hdes%n,0:hdes%nel(1)):: W_wa
  INTEGER,DIMENSION(0:hdes%n,0:hdes%nel(2)):: W_wb
  INTEGER,DIMENSION(0:1,hdes%n,hdes%nel(1)):: Y_wa
  INTEGER,DIMENSION(0:1,hdes%n,hdes%nel(2)):: Y_wb
  REAL(wp),DIMENSION(:,:),ALLOCATABLE:: CI_H_DIAG
  INTEGER(KIND=intkind2),DIMENSION(:,:),ALLOCATABLE:: W_list 
  INTEGER:: istat  

  CALL CI_ASSEMBLER(hdes%n,hdes%nel(1),hdes%nel(2),dimen_a,dimen_b, &
                    W_wa,W_wb,Y_wa,Y_wb,a_str_arr,b_str_arr)

  IF(.NOT. ALLOCATED(CI_H_DIAG)) THEN
    ALLOCATE(CI_H_DIAG(dimen_a,dimen_b),STAT=istat)
    IF(istat/=0) THEN
      IF(comm_world%mpirank==0) THEN
        WRITE(*,*) "Allocation of CI_H_DIAG failed"
      ENDIF
    ENDIF
  ENDIF
  
  CALL CALC_CI_DIAG_HAM(a_str_arr,b_str_arr,dimen_a,dimen_b,HAM_1,INT_2, &
                        hdes%n,hdes%nel(1),hdes%nel(2),Y_wa,Y_wb,SYM%SYM,CI_H_DIAG)
  IF(comm_world%mpirank==0) THEN
    WRITE(*,*) "Diagonal CI matrix elements calculated; Energy of HF det = ", CI_H_DIAG(1,1)
  ENDIF

  ALLOCATE(W_list(dimen_a,dimen_b),STAT=istat)
  IF(istat/=0) THEN
    IF(comm_world%mpirank==0) THEN
      WRITE(*,*) "Allocation of W_list failed"
    ENDIF
  ENDIF

  W_list = 0 
  W_list(1,1) = INFO%W_start

  CALL QMCI_PROC(STRUC,hdes,INFO,SYM,EN,HAM_1,INT_2,CI_H_DIAG,W_list,dimen_a,dimen_b, &
                 a_str_arr,b_str_arr,Y_wa,Y_wb)

END SUBROUTINE DO_FCIQMC


SUBROUTINE CALC_ONE_DMAT(STRUC,hdes,SYM,dimen_a,dimen_b,a_str_arr,b_str_arr,Y_wa,Y_wb,CI_vector,ONE_DMAT)
  TYPE(Structure),INTENT(IN):: STRUC
  type(hamil_descriptor), intent(in) :: hdes
  TYPE(Symmetry),INTENT(IN):: SYM
  INTEGER,INTENT(IN):: dimen_a, dimen_b
  INTEGER(KIND=intkind),DIMENSION(dimen_a),INTENT(IN):: a_str_arr
  INTEGER(KIND=intkind),DIMENSION(dimen_b),INTENT(IN):: b_str_arr
  INTEGER,DIMENSION(0:1,hdes%n,hdes%nel(1)),INTENT(IN):: Y_wa
  INTEGER,DIMENSION(0:1,hdes%n,hdes%nel(2)),INTENT(IN):: Y_wb
  REAL(wp),DIMENSION(dimen_a,dimen_b),INTENT(IN):: CI_vector
  REAL(wp),INTENT(OUT):: ONE_DMAT(:,:,:)

  !local variables
  INTEGER:: i
  INTEGER:: p, q
  INTEGER:: counter
  INTEGER:: adrl , adrr
  INTEGER,DIMENSION(dimen_a):: ALIST_L, ALIST_R, APHASE
  INTEGER,DIMENSION(dimen_b):: BLIST_L, BLIST_R, BPHASE
 
  ONE_DMAT = 0.0_wp
 
  IF(SYM%SYM==0) THEN
    ! For open-shell calculations and Ms /= 0
    DO q=1,hdes%n
      DO p=q,hdes%n   

        !Part for alpha spin 
        !Get a list of excitations made with E_{pq}^{alpha}
        CALL CREATE_S_EXC_LISTpq(p,q,hdes%n,hdes%nel(1),dimen_a,a_str_arr,Y_wa,counter,ALIST_L,ALIST_R,APHASE) 
        DO i=1,counter
          adrl = ALIST_L(i)
          adrr = ALIST_R(i)

          ONE_DMAT(p,q,1) = ONE_DMAT(p,q,1) + DOT_PRODUCT(CI_vector(adrl,:),CI_vector(adrr,:))*APHASE(i)
        ENDDO

        ONE_DMAT(q,p,1) = ONE_DMAT(p,q,1)

        !Part for beta spin
        !Get a list of beta excitations with E_{pq}^{beta}
        CALL CREATE_S_EXC_LISTpq(p,q,hdes%n,hdes%nel(2),dimen_b,b_str_arr,Y_wb,counter,BLIST_L,BLIST_R,BPHASE)
        DO i=1,counter
          adrl = BLIST_L(i)
          adrr = BLIST_R(i)

          ONE_DMAT(p,q,2) = ONE_DMAT(p,q,2) + DOT_PRODUCT(CI_vector(:,adrl),CI_vector(:,adrr))*BPHASE(i)
        ENDDO

        ONE_DMAT(q,p,2) = ONE_DMAT(p,q,2)

      ENDDO
    ENDDO

  ELSE
  !For closed-shell calculations with  Ms = 0
    DO q=1,hdes%n
      DO p=q,hdes%n   
       
 
        CALL CREATE_S_EXC_LISTpq(p,q,hdes%n,hdes%nel(1),dimen_a,a_str_arr,Y_wa,counter,ALIST_L,ALIST_R,APHASE) 
        DO i=1,counter
          adrl = ALIST_L(i)
          adrr = ALIST_R(i)

          ONE_DMAT(p,q,1) = ONE_DMAT(p,q,1) + DOT_PRODUCT(CI_vector(adrl,:),CI_vector(adrr,:))*APHASE(i)
        ENDDO
 
        ONE_DMAT(q,p,1) = ONE_DMAT(p,q,1)
      
      ENDDO
    ENDDO

    ONE_DMAT = 2.0_wp * ONE_DMAT

  ENDIF
END SUBROUTINE CALC_ONE_DMAT


SUBROUTINE CALC_TWO_DMAT(STRUC,hdes,SYM,dimen_a,dimen_b,a_str_arr,b_str_arr,Y_wa,Y_wb,CI_vector,ONE_DMAT,TWO_DMAT)
  TYPE(Structure),INTENT(IN):: STRUC
  type(hamil_descriptor), intent(in) :: hdes
  TYPE(Symmetry),INTENT(IN):: SYM
  INTEGER,INTENT(IN):: dimen_a, dimen_b
  INTEGER(KIND=intkind),DIMENSION(dimen_a),INTENT(IN):: a_str_arr
  INTEGER(KIND=intkind),DIMENSION(dimen_b),INTENT(IN):: b_str_arr
  INTEGER,DIMENSION(0:1,hdes%n,hdes%nel(1)),INTENT(IN):: Y_wa
  INTEGER,DIMENSION(0:1,hdes%n,hdes%nel(2)),INTENT(IN):: Y_wb
  REAL(wp),DIMENSION(dimen_a,dimen_b),INTENT(IN):: CI_vector
  REAL(wp),INTENT(IN):: ONE_DMAT(:,:,:)
  REAL(wp),DIMENSION(hdes%n,hdes%n,hdes%n,hdes%n),INTENT(OUT):: TWO_DMAT 

  !local variables
  INTEGER:: i, j
  INTEGER:: p, q, r, s
  INTEGER:: counter, counter2
  INTEGER:: adrl , adrr, bdrl, bdrr
  INTEGER,DIMENSION(dimen_a):: ALIST_L, ALIST_R, APHASE
  INTEGER,DIMENSION(dimen_b):: BLIST_L, BLIST_R, BPHASE
  
  TWO_DMAT = 0.0_wp
  
  IF(SYM%SYM==1) THEN
    DO s=1,hdes%n
      DO r=1,hdes%n
        DO q=1,hdes%n
          DO p=1,hdes%n

            CALL CREATE_D_EXC_LISTpqrs(p,q,r,s,hdes%n,hdes%nel(1),a_str_arr,dimen_a,Y_wa,counter,ALIST_L,ALIST_R,APHASE)
            DO i=1,counter
              adrl = ALIST_L(i)
              adrr = ALIST_R(i)

              TWO_DMAT(p,q,r,s) = TWO_DMAT(p,q,r,s) + DOT_PRODUCT(CI_vector(adrl,:),CI_vector(adrr,:)) * APHASE(i)
            ENDDO
                       
          ENDDO
        ENDDO
      ENDDO
    ENDDO
   
    DO s=1,hdes%n
      DO r=1,hdes%n

        CALL CREATE_S_EXC_LISTpq(r,s,hdes%n,hdes%nel(1),dimen_a,a_str_arr,Y_wa,counter,ALIST_L,ALIST_R,APHASE)

        DO q=1,hdes%n
          DO p=1,hdes%n

            CALL CREATE_S_EXC_LISTpq(p,q,hdes%n,hdes%nel(2),dimen_b,b_str_arr,Y_wb,counter2,BLIST_L,BLIST_R,BPHASE) 
            
            DO i=1,counter
              adrl = ALIST_L(i)
              adrr = ALIST_R(i)
              DO j=1,counter2
                bdrl = BLIST_L(j)
                bdrr = BLIST_R(j)


                TWO_DMAT(p,q,r,s) = TWO_DMAT(p,q,r,s) + APHASE(i)*BPHASE(j)*CI_vector(adrl,bdrl) * CI_vector(adrr,bdrr)
              ENDDO
            ENDDO
                       
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    TWO_DMAT = TWO_DMAT * 2.0_wp

  ELSE
    DO s=1,hdes%n
      DO r=1,hdes%n
        DO q=1,hdes%n
          DO p=1,hdes%n

            CALL CREATE_D_EXC_LISTpqrs(p,q,r,s,hdes%n,hdes%nel(1),a_str_arr,dimen_a,Y_wa,counter,ALIST_L,ALIST_R,APHASE)
            DO i=1,counter
              adrl = ALIST_L(i)
              adrr = ALIST_R(i)

              TWO_DMAT(p,q,r,s) = TWO_DMAT(p,q,r,s) + DOT_PRODUCT(CI_vector(adrl,:),CI_vector(adrr,:)) * APHASE(i)
            ENDDO
                       
            CALL CREATE_D_EXC_LISTpqrs(p,q,r,s,hdes%n,hdes%nel(2),b_str_arr,dimen_b,Y_wb,counter,BLIST_L,BLIST_R,BPHASE)
            DO i=1,counter
              bdrl = BLIST_L(i)
              bdrr = BLIST_R(i)

              TWO_DMAT(p,q,r,s) = TWO_DMAT(p,q,r,s) + DOT_PRODUCT(CI_vector(:,bdrl),CI_vector(:,bdrr)) * BPHASE(i)
            ENDDO
                       
          ENDDO
        ENDDO
      ENDDO
    ENDDO
   

    DO s=1,hdes%n
      DO r=1,hdes%n

        CALL CREATE_S_EXC_LISTpq(r,s,hdes%n,hdes%nel(1),dimen_a,a_str_arr,Y_wa,counter,ALIST_L,ALIST_R,APHASE)

        DO q=1,hdes%n
          DO p=1,hdes%n

            CALL CREATE_S_EXC_LISTpq(p,q,hdes%n,hdes%nel(2),dimen_b,b_str_arr,Y_wb,counter2,BLIST_L,BLIST_R,BPHASE) 
            
            DO i=1,counter
              adrl = ALIST_L(i)
              adrr = ALIST_R(i)
              DO j=1,counter2
                bdrl = BLIST_L(j)
                bdrr = BLIST_R(j)


                TWO_DMAT(p,q,r,s) = TWO_DMAT(p,q,r,s) + APHASE(i)*BPHASE(j)*CI_vector(adrl,bdrl) * CI_vector(adrr,bdrr)
              ENDDO
            ENDDO
                       
          ENDDO
        ENDDO
      ENDDO
    ENDDO


    DO s=1,hdes%n
      DO r=1,hdes%n

        CALL CREATE_S_EXC_LISTpq(r,s,hdes%n,hdes%nel(2),dimen_b,b_str_arr,Y_wb,counter,BLIST_L,BLIST_R,BPHASE)

        DO q=1,hdes%n
          DO p=1,hdes%n

            CALL CREATE_S_EXC_LISTpq(p,q,hdes%n,hdes%nel(1),dimen_a,a_str_arr,Y_wa,counter2,ALIST_L,ALIST_R,APHASE) 
            
            DO i=1,counter
              bdrl = BLIST_L(i)
              bdrr = BLIST_R(i)
              DO j=1,counter2
                adrl = ALIST_L(j)
                adrr = ALIST_R(j)


                TWO_DMAT(p,q,r,s) = TWO_DMAT(p,q,r,s) + APHASE(j)*BPHASE(i)*CI_vector(adrl,bdrl) * CI_vector(adrr,bdrr)
              ENDDO
            ENDDO
                       
          ENDDO
        ENDDO
      ENDDO
    ENDDO

  ENDIF
 
 
  DO i=1,hdes%n
    TWO_DMAT(:,i,i,:) = TWO_DMAT(:,i,i,:) -ONE_DMAT(:,:,1)
  ENDDO

END SUBROUTINE CALC_TWO_DMAT

SUBROUTINE CREATE_S_EXC_LISTpq(p,q,ORB_NUM,NELECT,dimen,str_arr,Y_w,counter,LIST_L,LIST_R,LPHASE)
  INTEGER,INTENT(IN):: p, q
  INTEGER,INTENT(IN):: ORB_NUM, NELECT
  INTEGER,INTENT(IN):: dimen
  INTEGER(KIND=intkind),DIMENSION(dimen):: str_arr
  INTEGER,DIMENSION(0:1,ORB_NUM,NELECT):: Y_w
  INTEGER,INTENT(OUT):: counter
  INTEGER,DIMENSION(dimen),INTENT(OUT):: LIST_L , LIST_R , LPHASE

  !lcoal varaibles
  INTEGER:: i,j
  INTEGER(KIND=intkind) :: str , tempstr

  counter = 0

  DO i=1,dimen
    str = str_arr(i)
    IF(btest(str,q-1)) THEN
      tempstr = ibclr(str,q-1)
      IF(.not. btest(tempstr,p-1)) THEN
        tempstr = ibset(tempstr,p-1)
        counter = counter + 1 
        j = CALC_ADRESS_BINARY(ORB_NUM,NELECT,tempstr,Y_w)
        LIST_L(counter) = i
        LIST_R(counter) = j
        LPHASE(counter) = CI_PHASE(str,p,q)
      ENDIF
    ENDIF
  ENDDO
END SUBROUTINE CREATE_S_EXC_LISTpq

SUBROUTINE CREATE_D_EXC_LISTpqrs(p,q,r,s,ORB_NUM,NELECT,str_arr,dimen,Y_w,counter,LIST_L,LIST_R,LPHASE)
  INTEGER,INTENT(IN):: p, q, r, s
  INTEGER,INTENT(IN):: ORB_NUM, NELECT
  INTEGER,INTENT(IN):: dimen
  INTEGER(KIND=intkind),DIMENSION(dimen):: str_arr
  INTEGER,DIMENSION(0:1,ORB_NUM,NELECT):: Y_w
  INTEGER,INTENT(OUT):: counter
  INTEGER,DIMENSION(dimen),INTENT(OUT):: LIST_L , LIST_R , LPHASE

  !lcoal varaibles
  INTEGER:: i,j
  INTEGER(KIND=intkind) :: str , tempstr , tempstr2

  counter = 0

  DO i=1,dimen
    str = str_arr(i)
    IF(btest(str,s-1)) THEN
      tempstr = ibclr(str,s-1)
      IF(.not. btest(tempstr,r-1)) THEN
        tempstr = ibset(tempstr,r-1)
        IF(btest(tempstr,q-1)) THEN
          tempstr2 = ibclr(tempstr,q-1)
          IF(.not. btest(tempstr2,p-1)) THEN
            tempstr2 = ibset(tempstr2,p-1)    
            counter = counter + 1
            j = CALC_ADRESS_BINARY(ORB_NUM,NELECT,tempstr2,Y_w) 
            LIST_L(counter) = i
            LIST_R(counter) = j   
            LPHASE(counter) = CI_PHASE(str,r,s) * CI_PHASE(tempstr,p,q)
          ENDIF
        ENDIF
      ENDIF
    ENDIF
  ENDDO

END SUBROUTINE CREATE_D_EXC_LISTpqrs


!********************************************************************************
!
!  This Subroutine calculates the total FCI energy with reduced one-particle
!                          and two-particle denisties
!
!********************************************************************************
SUBROUTINE CALC_EN_DMAT(ORB_NUM,HAM_1,INT_2,ONE_DMAT,TWO_DMAT,ENDMAT)
  INTEGER,INTENT(IN):: ORB_NUM
  REAL(wp),DIMENSION(ORB_NUM,ORB_NUM),INTENT(IN):: HAM_1
  real(wp), intent(in) :: ONE_DMAT(:,:,:)
  REAL(wp),DIMENSION(ORB_NUM,ORB_NUM,ORB_NUM,ORB_NUM),INTENT(IN):: INT_2, TWO_DMAT
  REAL(wp),INTENT(OUT):: ENDMAT
  !local variables

  ENDMAT = 0.0_wp

  ENDMAT= SUM(HAM_1*ONE_DMAT(:,:,1)) + 0.5_wp*SUM(INT_2*TWO_DMAT)

END SUBROUTINE CALC_EN_DMAT


!********************************************************************************
!
!  This Subroutine calculates the energy of a fragment in DMET
!        difference to CALC_EN_DMAT is the fact that one      
!        index is always restricted to the local orbitals
!
!********************************************************************************
SUBROUTINE CALC_EN_DMAT_EMB(ORB_NUM,ORB_LOC,HAM_1,INT_2,ONE_DMAT,TWO_DMAT,ENDMAT)
  INTEGER,INTENT(IN):: ORB_NUM, ORB_LOC
  REAL(wp),DIMENSION(ORB_NUM,ORB_NUM),INTENT(IN):: HAM_1, ONE_DMAT
  REAL(wp),DIMENSION(ORB_NUM,ORB_NUM,ORB_NUM,ORB_NUM),INTENT(IN):: INT_2, TWO_DMAT
  REAL(wp),INTENT(OUT):: ENDMAT
  !local variables
  INTEGER:: p, q, r, s

  ENDMAT = 0.0_wp

  DO s=1,ORB_LOC
    DO r=1, ORB_NUM
      ENDMAT = ENDMAT + HAM_1(r,s) * ONE_DMAT(r,s)
      DO q=1,ORB_NUM
        DO p=1,ORB_NUM
          ENDMAT = ENDMAT + 0.5_wp*INT_2(p,q,r,s)*TWO_DMAT(p,q,r,s)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
END SUBROUTINE CALC_EN_DMAT_EMB


!********************************************************************************
!
!  This Subroutine calculates the energy of a fragment in DMET
!        difference to CALC_EN_DMAT is the fact that one      
!        index is always restricted to the local orbitals
!
!********************************************************************************
SUBROUTINE CALC_EN_DMAT_EMB2(ORB_NUM,ORB_LOC,PROJ,HAM_1,INT_2,ONE_DMAT,TWO_DMAT,ENDMAT)
  INTEGER,INTENT(IN):: ORB_NUM, ORB_LOC
  REAL(wp),DIMENSION(ORB_NUM,ORB_LOC),INTENT(IN):: PROJ
  REAL(wp),DIMENSION(ORB_NUM,ORB_NUM),INTENT(IN):: HAM_1, ONE_DMAT
  REAL(wp),DIMENSION(ORB_NUM,ORB_NUM,ORB_NUM,ORB_NUM),INTENT(IN):: INT_2, TWO_DMAT
  REAL(wp),INTENT(OUT):: ENDMAT
  !local variables
  INTEGER:: p, q
  REAL(wp),DIMENSION(ORB_NUM,ORB_NUM):: RHOP
  REAL(wp),DIMENSION(ORB_NUM,ORB_NUM):: MAT

  CALL GEMM("N","T",ORB_NUM,ORB_NUM,ORB_LOC,1.0_wp,PROJ(:,1:ORB_LOC),ORB_NUM,PROJ(:,1:ORB_LOC),ORB_NUM,0.0_wp,RHOP,ORB_NUM)

  ENDMAT = 0.0_wp

  DO p=1,ORB_NUM
    DO q=1,ORB_NUM
      MAT(p,q) = SUM(INT_2(p,:,:,:)*TWO_DMAT(q,:,:,:)) + SUM(HAM_1(p,:)*ONE_DMAT(q,:))
    ENDDO
  ENDDO

  ENDMAT = SUM(RHOP*MAT)
END SUBROUTINE CALC_EN_DMAT_EMB2


SUBROUTINE CALC_NO(ham, ORB_NUM,ONE_DMAT,OCCUP,NO)
  type(Hamilton), intent(inout) :: ham
  INTEGER,INTENT(IN):: ORB_NUM
  REAL(wp), INTENT(IN):: ONE_DMAT(:,:,:)
  REAL(wp),INTENT(OUT):: OCCUP(:,:)
  REAL(wp),INTENT(OUT):: NO(:,:,:)
  !local varaibles
  INTEGER:: spin
  !
  do spin = 1, size(one_dmat, 3)
    CALL SYM_EIGV(ORB_NUM, one_dmat(:,:,spin), occup(:,spin), no(:,:,spin))
    CALL reverse(no(:,:,spin))
    CALL reverse(OCCUP(:,spin))
  end do
  !
  IF(comm_world%mpirank==0) THEN
    block
      use array_numpy_file_io, only: ArrayNumpyFileIO
      type(ArrayNumpyFileIO) :: numpy_io
      call numpy_io%save(ham%h_io%orbitals%fname, no)
      call numpy_io%save(ham%h_io%occup%fname, occup)
    end block
  ENDIF
END SUBROUTINE CALC_NO

SUBROUTINE WRITE_BLOCK_EW_arr(ndim,B_H)
  INTEGER,INTENT(IN):: ndim
  TYPE(BLOCKING_HAM),DIMENSION(ndim):: B_H

  !local variables
  INTEGER:: i , j , k , modul
  INTEGER:: a , b

  OPEN(UNIT=ci_eigenval,FILE="CI_EIGENVAL",STATUS="replace",ACTION="write")
  
  WRITE(ci_eigenval,*) "Number of Sz Sectors = " , ndim
  WRITE(ci_eigenval,*)
  WRITE(ci_eigenval,*)
  
  DO i=1,ndim
    WRITE(ci_eigenval,*) "Sector number " , i , "SPIN = " , B_H(i)%SPIN

    DO j=1,B_H(i)%dimen
      WRITE(ci_eigenval,*) "  EIG No. ", j , B_H(i)%CI_EW(j)
      DO k=1,B_H(i)%dimen

        modul = mod(k,B_H(i)%dimen_b)
        IF(modul==0) THEN
          a = k / B_H(i)%dimen_b
          b = B_H(i)%dimen_b
        ELSE
          a = k / B_H(i)%dimen_b + 1
          b = modul
        ENDIF
        WRITE(ci_eigenval,*) "    W:", B_H(i)%a_str_arr(a), B_H(i)%b_str_arr(b) , B_H(i)%CI_EV(k,j)
      ENDDO 
    ENDDO

  ENDDO
  CLOSE(UNIT=ci_eigenval)
END SUBROUTINE WRITE_BLOCK_EW_arr


SUBROUTINE WRITE_BLOCK_EW_one(B_H)
  TYPE(BLOCKING_HAM):: B_H

  !local variables
  INTEGER:: i , j , modul
  INTEGER:: a , b

  OPEN(UNIT=ci_eigenval,FILE="CI_EIGENVAL",STATUS="replace",ACTION="write")
  
  WRITE(ci_eigenval,*) "Number of Sz Sectors = " , 1
  WRITE(ci_eigenval,*)
  WRITE(ci_eigenval,*)
  
  WRITE(ci_eigenval,*) "Sector number " , 1 , "SPIN = " , B_H%SPIN

  DO i=1,B_H%dimen
    WRITE(ci_eigenval,*) "  EIG No. ", i , B_H%CI_EW(i)
    DO j=1,B_H%dimen

      modul = mod(j,B_H%dimen_b)
      IF(modul==0) THEN
        a = j / B_H%dimen_b
        b = B_H%dimen_b
      ELSE
        a = j / B_H%dimen_b + 1
        b = modul
      ENDIF
      WRITE(ci_eigenval,*) "    W:", B_H%a_str_arr(a), B_H%b_str_arr(b) , B_H%CI_EV(j,i)
    ENDDO 
  ENDDO

  CLOSE(UNIT=ci_eigenval)
END SUBROUTINE WRITE_BLOCK_EW_one


END MODULE FCI_WRAP
