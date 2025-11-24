! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

MODULE FCIQMC
!*******************************************************************************
!             Module that contains subroutines that perform                    *
!                  Quantum Monte Carlo configuration                           *
!                       interaction calculations                               *
!                         Author: Zoran Sukurma                                *
!*******************************************************************************  

#include "preproc.inc"
USE CONSTANTS
use mpi, only: comm_world
use qmcfort_pos
use standalone
use statistics
USE CI
use hamilton_vars, only: hamil_descriptor

IMPLICIT NONE

INTEGER:: doab=0 ,doaabb=0 , sab=0
INTEGER:: doabcnt=0 , doaabbcnt=0 , sabcnt=0
INTEGER:: n_died=0 , n_died_p=0 , n_died_m=0
INTEGER:: n_sp=0 
!INTEGER:: n_sp_p=0 , n_sp_m=0
INTEGER:: sp_contr=0 
!INTEGER:: sp_contr_p=0 ,sp_contr_m=0
INTEGER:: n_cl=0 
!INTEGER:: n_cl_p=0 , n_cl_m=0
INTEGER:: n_ann=0

INTEGER(KIND=8):: try_s=0 , try_dab=0 , try_daabb=0
INTEGER(KIND=8):: try_death=0 , try_clone=0

CONTAINS

!**************************************************
!      Calculate probabilities for all kind       *
!                of excitations                   *
!**************************************************
SUBROUTINE CALC_PROBS(ORB_NUM,NELECT_a,NELECT_b)
  INTEGER,INTENT(IN):: ORB_NUM,NELECT_a,NELECT_b
  
  !local variables
  INTEGER:: N_I_a , N_I_b
  INTEGER:: N_II_aa , N_II_bb , N_II_ab
  INTEGER:: N_I , N_II
  INTEGER:: NELECT 

  NELECT = NELECT_a + NELECT_b

  N_I_a = NELECT_a*(ORB_NUM-NELECT_a)
  N_I_b = NELECT_b*(ORB_NUM-NELECT_b)
  N_I = N_I_a + N_I_b

  N_II_aa = NELECT_a*(NELECT_a-1)*(ORB_NUM-NELECT_a)*(ORB_NUM-NELECT_a-1)/4
  N_II_bb = NELECT_B*(NELECT_b-1)*(ORB_NUM-NELECT_b)*(ORB_NUM-NELECT_b-1)/4
  N_II_ab = NELECT_a*NELECT_b*(ORB_NUM-NELECT_a)*(ORB_NUM-NELECT_b)
  N_II = N_II_aa + N_II_bb + N_II_ab

  IF(P_s==0.0_wp) THEN
    P_s = REAL(N_I)/REAL(N_I+N_II)
  ENDIF
  P_d = 1.0_wp - P_s

  P_s_mode = REAL(NELECT_a)/REAL(NELECT)

  P_d_mode_aa = REAL(NELECT_a*(NELECT_a-1))/REAL(NELECT*(NELECT-1))
  P_d_mode_bb = REAL(NELECT_b*(NELECT_b-1))/REAL(NELECT*(NELECT-1))
  P_d_mode_ab = 2.0_wp*REAL(NELECT_a*NELECT_b)/REAL(NELECT*(NELECT-1))

  P_I_a = P_s*P_s_mode/REAL(NELECT_a*(ORB_NUM-NELECT_a))
  P_I_b = P_s*(1.0_wp-P_s_mode)/REAL(NELECT_b*(ORB_NUM-NELECT_b))

  IF(NELECT_a>1 .and. (ORB_NUM-NELECT_a)>1) THEN
    P_II_aa = P_d*P_d_mode_aa*4.0_wp/REAL(NELECT_a*(NELECT_a-1)* &
             (ORB_NUM-NELECT_a)*(ORB_NUM-NELECT_a-1))
  ELSE
    P_II_aa = 0.0_wp
  ENDIF

  IF(NELECT_b>1 .and. (ORB_NUM-NELECT_b)>1) THEN
    P_II_bb = P_d*P_d_mode_bb*4.0_wp/REAL(NELECT_b*(NELECT_b-1)* &
              (ORB_NUM-NELECT_b)*(ORB_NUM-NELECT_b-1))
  ELSE
    P_II_bb = 0.0_wp
  ENDIF

  P_II_ab = 1.0_wp*P_d*P_d_mode_ab/REAL(NELECT_a*NELECT_b*(ORB_NUM-NELECT_a)*(ORB_NUM-NELECT_b))

  IF(comm_world%mpirank==0) THEN
    WRITE(*,*)"probabilities"
    write(*,*)"P_s" , P_s
    write(*,*)"P_d" , P_d
    write(*,*)"P_s_mode" , P_s_mode
    write(*,*)"P_d_mode_aa" , P_d_mode_aa
    write(*,*)"P_d_mode_ab" , P_d_mode_ab
    write(*,*)"P_d_mode_bb" , P_d_mode_bb
    write(*,*)"P_I_a" , P_I_a
    write(*,*)"P_I_b" , P_I_b
    write(*,*)"P_II_aa" , P_II_aa
    write(*,*)"P_II_ab" , P_II_ab
    write(*,*)"P_II_bb" , P_II_bb
    write(*,*)"Normalization" , P_I_a*N_I_a + P_I_b*N_I_b + &
            P_II_aa*N_II_aa + P_II_ab*N_II_ab + P_II_bb*N_II_bb
  ENDIF
END SUBROUTINE CALC_PROBS

SUBROUTINE CREATE_LISTS(str,ORB_NUM,NELECT,occ_list,unocc_list)
!**************************************************
!    Creates List of occupied and unoccupied      *
!          orbitals in a given string             *
!**************************************************
  INTEGER(KIND=intkind),INTENT(IN):: str
  INTEGER,INTENT(IN):: ORB_NUM,NELECT
  INTEGER,DIMENSION(NELECT),INTENT(OUT):: occ_list 
  INTEGER,DIMENSION(ORB_NUM-NELECT),INTENT(OUT):: unocc_list
  
  !local variables
  INTEGER:: i
  INTEGER:: size1 , size2
  size1 = 0 
  size2 = 0
  DO i=1,ORB_NUM
    IF(btest(str,i-1)) THEN
      size1 = size1 + 1
      occ_list(size1) = i
    ELSE
      size2 = size2 + 1 
      unocc_list(size2) = i
    ENDIF
  ENDDO
END SUBROUTINE CREATE_LISTS

SUBROUTINE DO_SINGLE_EXC(str,str_new,occ1,unocc1,ORB_NUM,NELECT)
!**************************************************
!    Make a single excitation of a given string   *
!**************************************************
  INTEGER,INTENT(IN):: ORB_NUM , NELECT
  INTEGER(KIND=intkind),INTENT(IN):: str
  INTEGER(KIND=intkind),INTENT(OUT):: str_new
  INTEGER,INTENT(out):: occ1 , unocc1
 
  !local variables
  INTEGER,DIMENSION(NELECT):: occ_list
  INTEGER,DIMENSION(ORB_NUM-NELECT):: unocc_list
  REAL(wp):: rand
  INTEGER:: ind

  CALL CREATE_LISTS(str,ORB_NUM,NELECT,occ_list,unocc_list)

  !pick up an occupied orbital from str
  CALL RANDOM_NUMBER(rand)
  ind = 1 + FLOOR(rand*NELECT)
  occ1 = occ_list(ind)

  !pick up an uoccupied orbital from str
  CALL RANDOM_NUMBER(rand)
  ind = 1 + FLOOR(rand*(ORB_NUM-NELECT))
  unocc1 = unocc_list(ind)

  str_new = ibset(ibclr(str,occ1-1),unocc1-1)

END SUBROUTINE DO_SINGLE_EXC

SUBROUTINE DO_SINGLE_EXC_1(str,str_new,ORB_NUM,i,j)
!**************************************************
!       Pick up one single excitation             *
!**************************************************
  INTEGER,INTENT(IN):: ORB_NUM
  INTEGER(KIND=intkind),INTENT(IN):: str
  INTEGER(KIND=intkind),INTENT(OUT):: str_new
  INTEGER,INTENT(OUT):: i,j

  !local variables
  INTEGER:: ind
  REAL(wp):: rand

  !pick up occupied orbital
  DO
    CALL RANDOM_NUMBER(rand)
    ind = FLOOR(rand*ORB_NUM)
    IF (btest(str,ind)) THEN
      i = ind
      EXIT
    ENDIF
  ENDDO

  !pick up unoccupied orbital
  DO
    CALL RANDOM_NUMBER(rand)
    ind = FLOOR(rand*ORB_NUM)
    IF(.NOT. btest(str,ind)) THEN
      j = ind
      EXIT
    ENDIF
  ENDDO

  str_new = ibset(ibclr(str,i),j)

END SUBROUTINE DO_SINGLE_EXC_1

SUBROUTINE DO_DOUBLE_EXC_aa(str,str_new,occ1,occ2,unocc1,unocc2,ORB_NUM,NELECT)
!**************************************************
!      Make a double alpha-alpha or beta-beta     *
!           excitation of a given string          * 
!**************************************************
  INTEGER,INTENT(IN):: ORB_NUM , NELECT
  INTEGER(KIND=intkind),INTENT(IN):: str
  INTEGER(KIND=intkind),INTENT(OUT):: str_new
  INTEGER,INTENT(out):: occ1 , occ2
  INTEGER,INTENT(out):: unocc1 , unocc2

  !local variables
  INTEGER,DIMENSION(NELECT):: occ_list
  INTEGER,DIMENSION(ORB_NUM-NELECT):: unocc_list
  REAL(wp):: rand
  INTEGER:: ind1 , ind2

  CALL CREATE_LISTS(str,ORB_NUM,NELECT,occ_list,unocc_list)

  !pick up two occupied orbitals from str
  DO
    CALL RANDOM_NUMBER(rand)
    ind1 = 1 + FLOOR(rand*NELECT)
  
    CALL RANDOM_NUMBER(rand)
    ind2 = 1 + FLOOR(rand*NELECT)
 
    IF(ind2/=ind1) EXIT
  ENDDO
  
  occ1 = occ_list(ind1)
  occ2 = occ_list(ind2)

  !pick up two unoccupied orbitals from str
  DO
    CALL RANDOM_NUMBER(rand)
    ind1 = 1 + FLOOR(rand*(ORB_NUM-NELECT))

    CALL RANDOM_NUMBER(rand)
    ind2 = 1 + FLOOR(rand*(ORB_NUM-NELECT))

    IF(ind1/=ind2) EXIT
  ENDDO
 
  unocc1 = unocc_list(ind1)
  unocc2 = unocc_list(ind2)

  str_new = ibset(ibset(ibclr(ibclr(str,occ1-1),occ2-1),unocc1-1),unocc2-1)

END SUBROUTINE DO_DOUBLE_EXC_aa

SUBROUTINE DO_DOUBLE_EXC_aa_1(str,str_new,ORB_NUM,i,j,k,l)
!**************************************************
!    Pick up a-a or b-b double excitation         *
!**************************************************
  INTEGER(KIND=intkind),INTENT(IN)::str
  INTEGER(KIND=intkind),INTENT(OUT):: str_new
  INTEGER,INTENT(IN):: ORB_NUM
  INTEGER,INTENT(OUT):: i,j    !occupied orbitals in str
  INTEGER,INTENT(OUT):: k,l    !unoccupied orbitals in str

  !local variables
  REAL(wp):: rand1,rand2
  INTEGER:: ind1,ind2
  INTEGER:: temp

  !pick up two occupied orbitals
  DO
    CALL RANDOM_NUMBER(rand1)
    ind1 = FLOOR(rand1*ORB_NUM)

    DO
      CALL RANDOM_NUMBER(rand2)
      ind2 = FLOOR(rand2*ORB_NUM)
      IF(ind1/=ind2) THEN
        EXIT
      ENDIF
    ENDDO

    IF(btest(str,ind1) .AND. btest(str,ind2)) THEN
      i = ind1
      j = ind2
      EXIT
    ENDIF
  ENDDO

  IF(i>j) THEN
    temp = j
    j = i
    i = temp
  ENDIF

  !pick up two unoccupied orbitals
  DO
    CALL RANDOM_NUMBER(rand1)
    ind1 = FLOOR(rand1*ORB_NUM)

    DO
      CALL RANDOM_NUMBER(rand2)
      ind2 = FLOOR(rand2*ORB_NUM)
      IF(ind1/=ind2) THEN
        EXIT
      ENDIF
    ENDDO

    IF((.NOT. btest(str,ind1)) .AND. (.NOT. btest(str,ind2))) THEN
      k = ind1
      l = ind2
      EXIT
    ENDIF
  ENDDO

  IF(k>l) THEN
    temp = l
    l = k
    k = temp
  ENDIF

  str_new = ibset(ibset(ibclr(ibclr(str,i),j),k),l)

END SUBROUTINE DO_DOUBLE_EXC_aa_1

SUBROUTINE DO_DOUBLE_EXC_ab(a_str,b_str,a_str_new,b_str_new,occ1,occ2,unocc1,unocc2, &
                            ORB_NUM,NELECT_a,NELECT_b)
!**************************************************
!     Make a double alpha-beta excitation of      *
!                  a given string                 *
!**************************************************
  INTEGER,INTENT(IN):: ORB_NUM , NELECT_a , NELECT_b
  INTEGER(KIND=intkind),INTENT(IN):: a_str , b_str
  INTEGER(KIND=intkind),INTENT(OUT):: a_str_new , b_str_new
  INTEGER,INTENT(out):: occ1 , occ2
  INTEGER,INTENT(out):: unocc1 , unocc2

  !local variables
  INTEGER,DIMENSION(NELECT_a):: occ_list_a
  INTEGER,DIMENSION(ORB_NUM-NELECT_a):: unocc_list_a
  INTEGER,DIMENSION(NELECT_b):: occ_list_b 
  INTEGER,DIMENSION(ORB_NUM-NELECT_b):: unocc_list_b
  REAL(wp):: rand
  INTEGER:: ind

  CALL CREATE_LISTS(a_str,ORB_NUM,NELECT_a,occ_list_a,unocc_list_a)
  CALL CREATE_LISTS(b_str,ORB_NUM,NELECT_b,occ_list_b,unocc_list_b)
 
  !pick up an occupied orbital from a_str
  CALL RANDOM_NUMBER(rand)
  ind = 1 + FLOOR(rand*NELECT_a)
  occ1 = occ_list_a(ind)

  !pick up an occupied orbital from b_str
  CALL RANDOM_NUMBER(rand)
  ind = 1 + FLOOR(rand*NELECT_b)
  occ2 = occ_list_b(ind)

  !pick up an unoccupied orbital from a_str
  CALL RANDOM_NUMBER(rand)
  ind = 1 + FLOOR(rand*(ORB_NUM-NELECT_a))
  unocc1 = unocc_list_a(ind)

  !pick up an unoccupied orbital from b_str
  CALL RANDOM_NUMBER(rand)
  ind = 1 + FLOOR(rand*(ORB_NUM-NELECT_b))
  unocc2 = unocc_list_b(ind)

  a_str_new = ibset(ibclr(a_str,occ1-1),unocc1-1)
  b_str_new = ibset(ibclr(b_str,occ2-1),unocc2-1)

END SUBROUTINE DO_DOUBLE_EXC_ab

SUBROUTINE DO_DOUBLE_EXC_ab_1(a_str,b_str,a_str_new,b_str_new,ORB_NUM,i,j,k,l)
!**************************************************
!         Pick up  a-b double excitation          *
!**************************************************
  INTEGER(KIND=intkind),INTENT(IN):: a_str , b_str
  INTEGER(KIND=intkind),INTENT(OUT):: a_str_new , b_str_new
  INTEGER,INTENT(IN):: ORB_NUM
  INTEGER,INTENT(OUT):: i,j     !occupied orbitals in both strings
  INTEGER,INTENT(OUT):: k,l     !unoccupied orbitals in both strings

  !local variables
  REAL(wp):: rand1 , rand2
  INTEGER:: ind1 , ind2

  !pick up two occupied orbitals
  DO
    CALL RANDOM_NUMBER(rand1)
    CALL RANDOM_NUMBER(rand2)
    ind1 = FLOOR(rand1*ORB_NUM)
    ind2 = FLOOR(rand2*ORB_NUM)

    IF(btest(a_str,ind1) .AND. btest(b_str,ind2)) THEN
      i = ind1
      j = ind2
      EXIT
    ENDIF
  ENDDO

  !pick up two unoccupied orbitals
  DO
    CALL RANDOM_NUMBER(rand1)
    CALL RANDOM_NUMBER(rand2)
    ind1 = FLOOR(rand1*ORB_NUM)
    ind2 = FLOOR(rand2*ORB_NUM)

    IF((.NOT. btest(a_str,ind1)) .AND. (.NOT. btest(b_str,ind2))) THEN
      k = ind1
      l = ind2
      EXIT
    ENDIF
  ENDDO

  a_str_new = ibset(ibclr(a_str,i),k)
  b_str_new = ibset(ibclr(b_str,j),l)

END SUBROUTINE DO_DOUBLE_EXC_ab_1


SUBROUTINE SPAWNING_STEP(INFO,STRUC,hdes,SYM,dimen_a,dimen_b,a_str_arr,b_str_arr,Y_wa,Y_wb, &
                         W_list,HAM_DIAG,HAM_1,INT_2,E_SHIFT)
!**************************************************
!   Spawning step procedure of  quantum MC CI     *
!     At the end the diagonal death/cloning       *
!             process is done, too                *
!**************************************************
  TYPE(CONTROLL),INTENT(INOUT):: INFO
  TYPE(STRUCTURE),INTENT(INOUT):: STRUC
  type(hamil_descriptor), intent(inout) :: hdes
  TYPE(SYMMETRY),INTENT(INOUT):: SYM 
  INTEGER,INTENT(IN):: dimen_a, dimen_b
  INTEGER(KIND=intkind),DIMENSION(dimen_a),INTENT(IN):: a_str_arr
  INTEGER(KIND=intkind),DIMENSION(dimen_b),INTENT(IN):: b_str_arr
  INTEGER,DIMENSION(0:1,hdes%n,hdes%nel(1)),INTENT(IN):: Y_wa
  INTEGER,DIMENSION(0:1,hdes%n,hdes%nel(2)),INTENT(IN):: Y_wb
  INTEGER(KIND=intkind2),DIMENSION(dimen_a,dimen_b),INTENT(OUT):: W_list
  REAL(wp),DIMENSION(hdes%n,hdes%n,INFO%ISPIN1),INTENT(IN):: HAM_1
  REAL(wp),DIMENSION(hdes%n,hdes%n,hdes%n,hdes%n,INFO%ISPIN2),INTENT(IN):: INT_2
  REAL(wp),DIMENSION(dimen_a,dimen_b),INTENT(IN):: HAM_DIAG
  REAL(wp),INTENT(INOUT):: E_SHIFT 

  !local variables
  INTEGER:: p,q
  INTEGER:: i
  INTEGER(KIND=intkind):: a_str , b_str , a_str_new , b_str_new
  REAL(wp):: rand
  REAL(wp):: EN
  INTEGER:: occ1 , occ2 , unocc1 , unocc2
  REAL(wp):: prob
  INTEGER(KIND=intkind2) :: walk
  INTEGER:: lowind
  INTEGER:: inda , indb , indab

  IF(INFO%ISPIN1==1) THEN
    inda = 1 
    indb = 1
    indab = 1 
  ELSEIF(INFO%ISPIN1==2) THEN
    inda = 1 
    indb = 2 
    indab = 3
    SYM%SYM = 0
  ELSE
    IF(comm_world%mpirank==0) THEN
      WRITE(*,*) "Unrecognized ISPIN value"
    ENDIF
    CALL EXIT
  ENDIF

  spawned = 0
  sp_contr = 0
  !sp_contr_p = 0
  !sp_contr_m = 0

  DO q=1,dimen_b
    b_str = b_str_arr(q)

    lowind = 1 
    IF(SYM%SYM==1 .and. SYM%SPIN==0) THEN
      lowind = q
    ELSE
      lowind = 1
    ENDIF

    DO p=lowind,dimen_a
      a_str = a_str_arr(p)
      walk = W_list(p,q)
      IF(walk/=0) THEN
        DO i=1,ABS(walk)

          !choose single/double excitation
          CALL RANDOM_NUMBER(rand)
          IF(rand<=P_s) THEN
            !PERFORM SINGLE EXCITAION
            !choose which single excitation
            CALL RANDOM_NUMBER(rand)
            IF(rand<=P_s_mode) THEN
              !perform alpha excitation
              CALL DO_SINGLE_EXC(a_str,a_str_new,occ1,unocc1,hdes%n,hdes%nel(1))
              !CALL DO_SINGLE_EXC_1(a_str,a_str_new,ORB_NUM,occ1,unocc1)
              b_str_new = b_str
              CALL SINGLE_EXC_EN(hdes%n,HAM_1(:,:,inda),INT_2,a_str,b_str,occ1,unocc1,EN)
              CALL SPAWN(P_I_a)
              try_s = try_s + 1 
            ELSE
              !perform beta excitation
              CALL DO_SINGLE_EXC(b_str,b_str_new,occ1,unocc1,hdes%n,hdes%nel(2))
              !CALL DO_SINGLE_EXC_1(b_str,b_str_new,ORB_NUM,occ1,unocc1)
              a_str_new = a_str
              CALL SINGLE_EXC_EN(hdes%n,HAM_1(:,:,indb),INT_2,b_str,a_str,occ1,unocc1,EN)
              CALL SPAWN(P_I_b)
              try_s = try_s + 1
            ENDIF
          ELSE
            !PERFORM DOUBLE EXCITATION
            CALL RANDOM_NUMBER(rand)
            !WRITE(*,*) "random for mode" , rand
            IF(rand<=P_d_mode_aa) THEN
              !perform alpha-alpha excitation
              CALL DO_DOUBLE_EXC_aa(a_str,a_str_new,occ1,occ2,unocc1,unocc2,hdes%n,hdes%nel(1))
              !CALL DO_DOUBLE_EXC_aa_1(a_str,a_str_new,ORB_NUM,occ1,occ2,unocc1,unocc2)
              b_str_new = b_str
              CALL DOUBLE_EXC_EN(hdes%n,INT_2(:,:,:,:,inda),a_str,b_str,occ1,occ2,unocc1,unocc2,EN,1)
              CALL SPAWN(P_II_aa)
              try_daabb = try_daabb + 1
            ELSEIF(rand>P_d_mode_aa .and. rand<=(P_d_mode_aa+P_d_mode_bb)) THEN
              !perform beta-beta excitation
              CALL DO_DOUBLE_EXC_aa(b_str,b_str_new,occ1,occ2,unocc1,unocc2,hdes%n,hdes%nel(2))
              !CALL DO_DOUBLE_EXC_aa_1(b_str,b_str_new,ORB_NUM,occ1,occ2,unocc1,unocc2)
              a_str_new = a_str
              CALL DOUBLE_EXC_EN(hdes%n,INT_2(:,:,:,:,indb),b_str,a_str,occ1,occ2,unocc1,unocc2,EN,1)
              CALL SPAWN(P_II_bb)
              try_daabb = try_daabb + 1
            ELSE
              !perform alpha-beta excitation
              CALL DO_DOUBLE_EXC_ab(a_str,b_str,a_str_new,b_str_new,occ1,occ2, &
                                  unocc1,unocc2,hdes%n,hdes%nel(1),hdes%nel(2))
              !CALL DO_DOUBLE_EXC_ab_1(a_str,b_str,a_str_new,b_str_new,ORB_NUM, &
              !                        occ1,occ2,unocc1,unocc2)
              CALL DOUBLE_EXC_EN(hdes%n,INT_2(:,:,:,:,indab),a_str,b_str,occ1,occ2,unocc1,unocc2,EN,2)
              CALL SPAWN(P_II_ab)
              try_dab = try_dab + 1
            ENDIF
          ENDIF

          !DO DIAGONAL DEATH/CLONING PROCESS - Version 1
          prob = INFO%MC_tstep*(HAM_DIAG(p,q)-E_SHIFT)
          CALL RANDOM_NUMBER(rand)
          IF(prob>=0.0_wp) THEN
            !death

            try_death = try_death + 1

            IF(rand<=prob) THEN
              n_died =  n_died + 1
              IF(walk>0) THEN
                !n_died_p = n_died_p + 1
                W_list(p,q) = W_list(p,q) - 1
              ELSEIF(walk<0) THEN
                !n_died_m = n_died_m + 1
                W_list(p,q) = W_list(p,q) + 1
              ENDIF
            ENDIF
          ELSE
            !cloning

            try_clone = try_clone + 1

            IF(rand<=ABS(prob)) THEN
              n_cl = n_cl + 1
              IF(walk>0) THEN
                !n_cl_p = n_cl_p + 1
                W_list(p,q) = W_list(p,q) + 1
              ELSEIF(walk<0) THEN
                !n_cl_m = n_cl_m + 1
                W_list(p,q) = W_list(p,q) - 1
              ENDIF
            ENDIF
          ENDIF

        
        ENDDO
      ENDIF
    ENDDO
  ENDDO

CONTAINS

SUBROUTINE SPAWN(PROB)
  REAL(wp),INTENT(IN):: PROB
    
  !local probability
  REAL(wp):: TOT_PROB
  REAL(wp):: rand
  INTEGER:: a_adr , b_adr , swap_adr
  INTEGER(KIND=intkind2):: contrib

  TOT_PROB  = INFO%MC_tstep*ABS(EN)/PROB
  contrib = FLOOR(TOT_PROB)

  CALL RANDOM_NUMBER(rand)
  IF(rand<=TOT_PROB-FLOOR(TOT_PROB)) THEN
    contrib = contrib + 1
  ENDIF   

  !spawn   
  IF(contrib>0) THEN
    a_adr = CALC_ADRESS_BINARY(hdes%n,hdes%nel(1),a_str_new,Y_wa)
    b_adr = CALC_ADRESS_BINARY(hdes%n,hdes%nel(2),b_str_new,Y_wb)

    !SYM=1 SPIN=0
    IF(SYM%SYM==1 .and. SYM%SPIN==0) THEN
      IF(b_adr > a_adr) THEN
        swap_adr = a_adr
        a_adr = b_adr
        b_adr = swap_adr
      ENDIF
    ENDIF

    IF(EN<0) THEN
      contrib = contrib * (walk/ABS(walk))
    ELSE
      contrib = -contrib * (walk/ABS(walk))
    ENDIF

    spawned = spawned +  1
    sp_array(spawned)%a_adr = a_adr
    sp_array(spawned)%b_adr = b_adr
    sp_array(spawned)%walk = contrib

    sp_contr = sp_contr + ABS(contrib)
    !IF(contrib>0) THEN
    !  sp_contr_p = sp_contr_p + contrib
    !ELSE
    !  sp_contr_m = sp_contr_m + ABS(contrib)
    !ENDIF

    IF(PROB==P_I_a) THEN
      sab=sab+ABS(contrib)
      sabcnt = sabcnt + 1
    ELSEIF(PROB==P_II_aa) THEN
      doaabb=doaabb+ABS(contrib)
      doaabbcnt = doaabbcnt + 1
    ELSE
      doab=doab+ABS(contrib)
      doabcnt = doabcnt + 1
    ENDIF
  ELSEIF(contrib<0) THEN
    WRITE(*,*)"WARNING, SOMETHING IS NOT OK, contrib < 0 "
  ENDIF

END SUBROUTINE SPAWN

END SUBROUTINE SPAWNING_STEP

SUBROUTINE MERGE_ANNIHILATE(dimen_a,dimen_b,W_list)
!**************************************************
!   Merge/annihilate new list of spawned walkers  *
!           with main list of walkers             *
!**************************************************
  INTEGER,INTENT(IN):: dimen_a , dimen_b 
  INTEGER(KIND=intkind2),DIMENSION(dimen_a,dimen_b),INTENT(INOUT):: W_list
  TYPE(WALKER_c):: s_walker

  !local variables
  INTEGER:: a , b 
  INTEGER:: i

  DO i=1,spawned
    s_walker = sp_array(i)
    a=s_walker%a_adr
    b=s_walker%b_adr
    IF(W_list(a,b)*s_walker%walk<0) THEN
      IF(ABS(W_list(a,b))>=ABS(s_walker%walk)) THEN
        n_ann = n_ann + ABS(s_walker%walk)*2
      ELSE
        n_ann = n_ann + ABS(W_list(a,b))*2
      ENDIF
    ENDIF
 
    W_list(a,b) = W_list(a,b) + s_walker%walk
  ENDDO

END SUBROUTINE MERGE_ANNIHILATE

SUBROUTINE CALC_HF_EXC_LIST(ORB_NUM,NELECT_a,NELECT_b,HF_str_a,HF_str_b, &
                      dim_list,HAM_1,INT_2,Y_w_a,Y_w_b,en_list,a_adr_list,b_adr_list, &
                      SYM,SPIN)
!**************************************************
!  Creates list of single and double excitations  *
!        of HF determinant and store the          *
!            energies  <HF|H|EXC>                 *
!**************************************************
  INTEGER,INTENT(IN):: ORB_NUM , NELECT_a , NELECT_b 
  INTEGER(KIND=intkind),INTENT(IN):: HF_str_a , HF_str_b
  INTEGER,INTENT(OUT):: dim_list
  INTEGER,DIMENSION(0:1,ORB_NUM,NELECT_a):: Y_w_a
  INTEGER,DIMENSION(0:1,ORB_NUM,NELECT_b):: Y_w_b
  REAL(wp),DIMENSION(ORB_NUM,ORB_NUM,INFO%ISPIN1),INTENT(IN):: HAM_1
  REAL(wp),DIMENSION(ORB_NUM,ORB_NUM,ORB_NUM,ORB_NUM,INFO%ISPIN2),INTENT(IN):: INT_2
  REAL(wp),DIMENSION(:),ALLOCATABLE,INTENT(OUT):: en_list
  INTEGER,DIMENSION(:),ALLOCATABLE,INTENT(OUT):: a_adr_list , b_adr_list
  INTEGER,INTENT(INOUT):: SYM , SPIN

  !local variables
  INTEGER:: Nsa , Nsb, Ndaa , Ndbb, Ndab
  INTEGER:: cnt
  INTEGER,DIMENSION(NELECT_a):: occ_list_a
  INTEGER,DIMENSION(ORB_NUM-NELECT_a):: unocc_list_a
  INTEGER,DIMENSION(NELECT_b):: occ_list_b
  INTEGER,DIMENSION(ORB_NUM-NELECT_b):: unocc_list_b
  INTEGER:: i,j,k,l
  INTEGER:: homeadr_a , homeadr_b
  INTEGER:: ind1 , ind2 , ind3 , ind4
  INTEGER(KIND=intkind):: tempstr_a  , tempstr_b
  REAL(wp):: EN
  INTEGER:: inda , indb , indab

  IF(INFO%ISPIN1==1) THEN
    inda = 1 
    indb = 1
    indab = 1
  ELSEIF(INFO%ISPIN1==2) THEN
    inda = 1
    indb = 2
    indab = 3 
    SYM = 0
  ELSE
    IF(comm_world%mpirank==0) THEN
      WRITE(*,*) "Unrecognized ISPIN value"
    ENDIF
    CALL EXIT
  ENDIF

  Nsa = NELECT_a*(ORB_NUM-NELECT_a)
  Nsb = NELECT_b*(ORB_NUM-NELECT_b)
  IF(NELECT_a>1) THEN
    Ndaa = NELECT_a*(NELECT_a-1)*(ORB_NUM-NELECT_a)*(ORB_NUM-NELECT_a-1)/4
  ELSE
    Ndaa = 0
  ENDIF
  IF(NELECT_b>1) THEN
    Ndbb = NELECT_b*(NELECT_b-1)*(ORB_NUM-NELECT_b)*(ORB_NUM-NELECT_b-1)/4
  ELSE
    Ndbb = 0
  ENDIF
  Ndab = NELECT_a*NELECT_b*(ORB_NUM-NELECT_a)*(ORB_NUM-NELECT_b)

  dim_list = Nsa + Nsb + Ndaa + Ndbb + Ndab

  IF(comm_world%mpirank==0) THEN
    WRITE(*,*) "Number of excitations of HF determinant= ", dim_list
  ENDIF

  ALLOCATE(en_list(dim_list))
  ALLOCATE(a_adr_list(dim_list))
  ALLOCATE(b_adr_list(dim_list))


  CALL CREATE_LISTS(HF_str_a,ORB_NUM,NELECT_a,occ_list_a,unocc_list_a)  
  CALL CREATE_LISTS(HF_str_b,ORB_NUM,NELECT_b,occ_list_b,unocc_list_b)
  
  homeadr_a = CALC_ADRESS_BINARY(ORB_NUM,NELECT_a,HF_str_a,Y_w_a)
  homeadr_b = CALC_ADRESS_BINARY(ORB_NUM,NELECT_b,HF_str_b,Y_w_b)

  !pick up single alpha excitations
  cnt=0

  DO i=1,NELECT_a
     ind1 = occ_list_a(i)
    DO j=1,ORB_NUM-NELECT_a
      ind2= unocc_list_a(j)
      tempstr_a = ibset(ibclr(HF_str_a,ind1-1),ind2-1)
     
      cnt = cnt + 1 
 
      a_adr_list(cnt) = CALC_ADRESS_BINARY(ORB_NUM,NELECT_a,tempstr_a,Y_w_a)
      b_adr_list(cnt) = homeadr_b
      
      CALL SINGLE_EXC_EN(ORB_NUM,HAM_1(:,:,inda),INT_2,HF_str_a,HF_str_b,ind1,ind2,EN)
      en_list(cnt) = EN
    ENDDO
  ENDDO
  
  !do the same for single beta excitations
  IF(SYM==1 .and. SPIN==0) THEN
    a_adr_list(cnt+1:2*cnt) = b_adr_list(1:cnt)
    b_adr_list(cnt+1:2*cnt) = a_adr_list(1:cnt)
    en_list(cnt+1:2*cnt) = en_list(1:cnt)

    cnt = 2*cnt
  ELSE
    DO i=1,NELECT_b
      ind1 = occ_list_b(i)
      DO j=1,ORB_NUM-NELECT_b
        ind2= unocc_list_b(j)
        tempstr_b = ibset(ibclr(HF_str_b,ind1-1),ind2-1)

        cnt = cnt + 1 

        b_adr_list(cnt) = CALC_ADRESS_BINARY(ORB_NUM,NELECT_b,tempstr_b,Y_w_b)
        a_adr_list(cnt) = homeadr_a

        CALL SINGLE_EXC_EN(ORB_NUM,HAM_1(:,:,indb),INT_2,HF_str_b,HF_str_a,ind1,ind2,EN)
        en_list(cnt) = EN
      ENDDO
    ENDDO
  ENDIF

  !do same for alpha-alpha excitations
  IF(NELECT_a>1) THEN
    DO i=1,NELECT_a
      ind1 = occ_list_a(i)
      DO j=i+1,NELECT_a
        ind2 = occ_list_a(j)
        DO k=1,ORB_NUM-NELECT_a
          ind3 = unocc_list_a(k)
          DO l=k+1,ORB_NUM-NELECT_a
            ind4 = unocc_list_a(l)
            tempstr_a = ibset(ibset(ibclr(ibclr(HF_str_a,ind1-1),ind2-1),ind3-1),ind4-1)
         
            cnt = cnt + 1
 
            a_adr_list(cnt) = CALC_ADRESS_BINARY(ORB_NUM,NELECT_a,tempstr_a,Y_w_a)
            b_adr_list(cnt) = homeadr_b

            CALL DOUBLE_EXC_EN(ORB_NUM,INT_2(:,:,:,:,inda),HF_str_a,HF_str_b,ind1,ind2,ind3,ind4,EN,1) 
            en_list(cnt) = EN
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  IF(NELECT_b>1) THEN
    !same for beta-beta excitations
    IF(SYM==1 .and. SPIN==0) THEN
      a_adr_list(cnt+1:cnt+Ndbb) = b_adr_list(cnt-Ndaa+1:cnt)
      b_adr_list(cnt+1:cnt+Ndbb) = a_adr_list(cnt-Ndaa+1:cnt)
      en_list(cnt+1:cnt+Ndbb) = en_list(cnt-Ndaa+1:cnt)
      cnt = cnt + Ndbb
    ELSE
      DO i=1,NELECT_b
        ind1 = occ_list_b(i)
        DO j=i+1,NELECT_b
          ind2 = occ_list_b(j)
          DO k=1,ORB_NUM-NELECT_b
            ind3 = unocc_list_b(k)
            DO l=k+1,ORB_NUM-NELECT_b
              ind4 = unocc_list_b(l)
              cnt = cnt + 1

              tempstr_b = ibset(ibset(ibclr(ibclr(HF_str_b,ind1-1),ind2-1),ind3-1),ind4-1)

              b_adr_list(cnt) = CALC_ADRESS_BINARY(ORB_NUM,NELECT_b,tempstr_b,Y_w_b)
              a_adr_list(cnt) = homeadr_a

              CALL DOUBLE_EXC_EN(ORB_NUM,INT_2(:,:,:,:,indb),HF_str_b,HF_str_a,ind1,ind2,ind3,ind4,EN,1)
              en_list(cnt) = EN
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDIF
  ENDIF

  !and at the end for alpha-beta excitations
  DO i=1,NELECT_a
    ind1 = occ_list_a(i)
    DO j=1,NELECT_b
      ind2 = occ_list_b(j)
      DO k=1,ORB_NUM-NELECT_a
        ind3 = unocc_list_a(k)
        DO l=1,ORB_NUM-NELECT_b
          ind4 = unocc_list_b(l)
 
          tempstr_a = ibset(ibclr(HF_str_a,ind1-1),ind3-1)
          tempstr_b = ibset(ibclr(HF_str_b,ind2-1),ind4-1)

          cnt = cnt + 1

          a_adr_list(cnt) = CALC_ADRESS_BINARY(ORB_NUM,NELECT_a,tempstr_a,Y_w_a)
          b_adr_list(cnt) = CALC_ADRESS_BINARY(ORB_NUM,NELECT_b,tempstr_b,Y_w_b)

          CALL DOUBLE_EXC_EN(ORB_NUM,INT_2(:,:,:,:,indab),HF_str_a,HF_str_b,ind1,ind2,ind3,ind4,EN,2)
          en_list(cnt) = EN
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  IF(comm_world%mpirank==0) THEN
    WRITE(*,*) "Number of picked up excitations = ", cnt
  ENDIF

END SUBROUTINE CALC_HF_EXC_LIST


FUNCTION CALC_CORR(dimen_a,dimen_b,dim_list,W_list,en_list,a_adr_list,b_adr_list) RESULT(E_corr)
!**************************************************
!     Calculates correlation energy using         *
!      created lists in CALC_HF_EXC_LIST          *
!**************************************************
  INTEGER,INTENT(IN):: dimen_a , dimen_b , dim_list
  INTEGER(KIND=intkind2),DIMENSION(dimen_a,dimen_b):: W_list
  INTEGER,DIMENSION(dim_list),INTENT(IN):: a_adr_list , b_adr_list
  REAL(wp),DIMENSION(dim_list),INTENT(IN):: en_list
  REAL(wp):: E_corr

  !local variables
  INTEGER:: i
  INTEGER:: a_adr , b_adr

  E_corr  = 0.0_wp

  DO i=1,dim_list
    a_adr = a_adr_list(i)
    b_adr = b_adr_list(i)
    !IF(SYM==1 .and. SPIN==0) THEN
    !  IF(a_adr>b_adr) THEN
    !     E_corr = E_corr + en_list(i) * W_list(a_adr,b_adr) * 2.0_wp 
    !  ELSEIF(a_adr==b_adr) THEN
    !     E_corr = E_corr + en_list(i) * W_list(a_adr,b_adr) 
    !  ENDIF
    !ELSE
      E_corr = E_corr + en_list(i) * W_list(a_adr,b_adr)
    !ENDIF
  ENDDO
 
  E_corr = E_corr/REAL(W_list(1,1),wp)

END FUNCTION CALC_CORR 

FUNCTION NUM_WALKERS(dimen_a,dimen_b,SYM,SPIN,W_list) RESULT(N_WALK)
!**************************************************
!       Calculates total number of walkers        *
!**************************************************
  INTEGER,INTENT(IN):: dimen_a , dimen_b
  INTEGER(KIND=intkind2),DIMENSION(dimen_a,dimen_b):: W_list
  INTEGER:: N_WALK
  INTEGER,INTENT(IN):: SYM, SPIN

  !local variables
  INTEGER:: p,q
  
  N_WALK = 0
  IF(SYM==1 .and. SPIN==0) THEN
    DO q=1,dimen_b
      DO p=q+1,dimen_a
        N_WALK = N_WALK + ABS(W_list(p,q))
      ENDDO
    ENDDO

    DO q=1,dimen_b
      N_WALK = N_WALK + ABS(W_list(q,q))
    ENDDO
  ELSE
    DO q=1,dimen_b
      DO p=1,dimen_a
        N_WALK = N_WALK + ABS(W_list(p,q))
      ENDDO
    ENDDO
  ENDIF

END FUNCTION

SUBROUTINE QMCI_PROC(STRUC,hdes,INFO,SYM,EN,HAM_1,INT_2,HAM_DIAG,W_list,dimen_a,dimen_b, &
                     a_str_arr,b_str_arr,Y_wa,Y_wb)
!**************************************************
!     Main procedure to perfrom QMCI:             *
!      1. Spawning for each walker                *
!      2. Diagonal cloning/death for each walker  *
!      3. Annihilation/Merge lists of walkers     *
!                                                 *
!     Three modes of simulations :                *
!      1. Constant S - exponential growth         *
!      2. Constant Nw simulation                  *
!**************************************************
  TYPE(CONTROLL),INTENT(INOUT):: INFO
  TYPE(STRUCTURE),INTENT(INOUT):: STRUC
  type(hamil_descriptor), intent(inout) :: hdes
  TYPE(fci_ENERGY),INTENT(INOUT):: EN
  TYPE(SYMMETRY),INTENT(INOUT):: SYM
  INTEGER,INTENT(IN):: dimen_a, dimen_b
  INTEGER(KIND=intkind),DIMENSION(dimen_a),INTENT(IN):: a_str_arr
  INTEGER(KIND=intkind),DIMENSION(dimen_b),INTENT(IN):: b_str_arr
  INTEGER,DIMENSION(0:1,hdes%n,hdes%nel(1)),INTENT(IN):: Y_wa
  INTEGER,DIMENSION(0:1,hdes%n,hdes%nel(2)),INTENT(IN):: Y_wb
  REAL(wp),DIMENSION(hdes%n,hdes%n,INFO%ISPIN1),INTENT(IN):: HAM_1
  REAL(wp),DIMENSION(hdes%n,hdes%n,hdes%n,hdes%n,INFO%ISPIN2),INTENT(IN):: INT_2
  REAL(wp),DIMENSION(dimen_a,dimen_b),INTENT(INOUT):: HAM_DIAG
  INTEGER(KIND=intkind2),DIMENSION(dimen_a,dimen_b),INTENT(INOUT):: W_list
   
  !local varaibles
  LOGICAL:: STOPPED
  INTEGER:: STOP_ST, Ncrit
  INTEGER,DIMENSION(:),ALLOCATABLE:: a_adr_list,b_adr_list
  REAL(wp),DIMENSION(:),ALLOCATABLE:: en_list
  INTEGER:: N_WALK , N_WALK_OLD
  INTEGER::  count 
  REAL(wp):: EHF
  INTEGER(KIND=intkind):: HF_str1 , HF_str2
  INTEGER::dim_list
  INTEGER::i , le, leng

  REAL(wp):: mean_1, mean_2
  REAL(wp):: var_1, var_2
  REAL(wp),DIMENSION(:),ALLOCATABLE:: EN_EST_1, EN_EST_2, EN_EST_3, EN_EST_4
  REAL(wp):: Norma

  Ncrit = INT(INFO%f_c*dimen_a*dimen_b)

  IF(INFO%LRANDOMIZE) THEN
    CALL RANDOM_SEED()
  ENDIF


  EHF = HAM_DIAG(1,1)
  HAM_DIAG = HAM_DIAG - EHF

  HF_str1 = a_str_arr(1)
  HF_str2 = b_str_arr(1)

  CALL CALC_PROBS(hdes%n,hdes%nel(1),hdes%nel(2))

  CALL CALC_HF_EXC_LIST(hdes%n,hdes%nel(1),hdes%nel(2),HF_str1,HF_str2, &
                      dim_list,HAM_1,INT_2,Y_wa,Y_wb,en_list,a_adr_list,b_adr_list,SYM%SYM,SYM%SPIN)

   
  !open OUTPUT FILE
  IF(comm_world%mpirank==0) THEN
    OPEN(UNIT=out,FILE="qmcfort_out",STATUS="old",POSITION="append",ACTION="write")
   
    WRITE(out,*) "START QUANTUM MONTE CARLO SIMULATION"
    WRITE(out,*) starshort

    ! Equilibration Phase
    WRITE(out,*) "EQUILIBRATION PHASE - Desired number of walkers:", Ncrit
  ENDIF

  DO i=1,INFO%EQ_STEPS

    CALL SPAWNING_STEP(INFO,STRUC,hdes,SYM,dimen_a,dimen_b,a_str_arr,b_str_arr,Y_wa,Y_wb, &
                       W_list,HAM_DIAG,HAM_1,INT_2,EN%E_SHIFT )

    n_sp = n_sp + sp_contr
    !n_sp_p = n_sp_p + sp_contr_p
    !n_sp_m = n_sp_m + sp_contr_m

    CALL MERGE_ANNIHILATE(dimen_a,dimen_b,W_list)

    IF(i==INFO%PRINT_STEP) THEN
      N_WALK = NUM_WALKERS(dimen_a,dimen_b,SYM%SYM,SYM%SPIN,W_list)
      IF(comm_world%mpirank==0) THEN
        WRITE(out,*) "STEP = " , i , "NW =" , N_WALK , "W1 =" , W_list(1,1)
      ENDIF

    ELSEIF(mod(i,INFO%PRINT_STEP)==0) THEN
      N_WALK = NUM_WALKERS(dimen_a,dimen_b,SYM%SYM,SYM%SPIN,W_list)
      IF(comm_world%mpirank==0) THEN
        WRITE(out,*) "STEP = " , i , "NW =" , N_WALK , "W1 =" , W_list(1,1)
      ENDIF

      IF(N_WALK>=Ncrit) THEN
        IF(comm_world%mpirank==0) THEN
          WRITE(out,*) "Number of walkers reached critical value"
        ENDIF
        STOPPED = .TRUE.
        STOP_ST = i
        EXIT
      ENDIF
    ENDIF
  ENDDO

  !WARM UP PHASE 
  IF(STOPPED) THEN
    count = STOP_ST
  ELSE
    count = INFO%EQ_STEPS
  ENDIF

  IF(INFO%WARM_UP==0) THEN
    INFO%WARM_UP = INT(0.1*count)
    IF(mod(INFO%WARM_UP,INFO%PRINT_STEP)==0) THEN
      INFO%WARM_UP = (INFO%WARM_UP/INFO%PRINT_STEP) * INFO%PRINT_STEP
    ELSE
      INFO%WARM_UP = (INFO%WARM_UP/INFO%PRINT_STEP+1) * INFO%PRINT_STEP
    ENDIF
  ENDIF

  IF(comm_world%mpirank==0) THEN
    WRITE(out,*) "CONSTANT NW SIMULATION - WARM UP; LENGTH:", INFO%WARM_UP
    WRITE(out,*)
  ENDIF

  DO i = count+1,count+INFO%WARM_UP

    CALL SPAWNING_STEP(INFO,STRUC,hdes,SYM,dimen_a,dimen_b,a_str_arr,b_str_arr,Y_wa,Y_wb, &
                       W_list,HAM_DIAG,HAM_1,INT_2,EN%E_SHIFT )

    n_sp = n_sp + sp_contr
    !n_sp_p = n_sp_p + sp_contr_p
    !n_sp_m = n_sp_m + sp_contr_m

    CALL MERGE_ANNIHILATE(dimen_a,dimen_b,W_list)

    IF(mod(i,INFO%A_per)==0) THEN
      N_WALK_old = N_WALK
      N_WALK = NUM_WALKERS(dimen_a,dimen_b,SYM%SYM,SYM%SPIN,W_list)
      EN%E_SHIFT = EN%E_SHIFT - INFO%DAMP*LOG(REAL(N_WALK,wp)/REAL(N_WALK_old,wp))/(INFO%A_per*INFO%MC_tstep)
    ENDIF

    IF(mod(i,INFO%PRINT_STEP)==0) THEN
      EN%EFCI_CORR = CALC_CORR(dimen_a,dimen_b,dim_list,W_list,en_list,a_adr_list,b_adr_list)
      IF(comm_world%mpirank==0) THEN
        WRITE(out,*) "STEP = ", i , "NW =" , N_WALK , "W1 =" , W_list(1,1) 
      ENDIF
    ENDIF
  ENDDO 

  !CONSTANTN NW SIMULATION - SAMPLING
  count = count + INFO%WARM_UP

  !determine length of constant Nw simulation mode
  IF(INFO%NUM_MC_STEPS==0) THEN
    IF(STOPPED) THEN
      INFO%NUM_MC_STEPS = ((INT(1.5*STOP_ST))/INFO%PRINT_STEP)*INFO%PRINT_STEP
    ELSE
      INFO%NUM_MC_STEPS = ((INT(1.5*INFO%EQ_STEPS))/INFO%PRINT_STEP)*INFO%PRINT_STEP
    ENDIF
  ELSE
    INFO%NUM_MC_STEPS = (INFO%NUM_MC_STEPS/INFO%PRINT_STEP)*INFO%PRINT_STEP
  ENDIF

  leng = INFO%NUM_MC_STEPS/INFO%A_per
  IF(mod(INFO%NUM_MC_STEPS,INFO%A_per)==0) THEN
    ALLOCATE(EN_EST_1(leng),EN_EST_2(leng),EN_EST_3(leng),EN_EST_4(leng))
  ELSE
    leng = leng +1
    ALLOCATE(EN_EST_1(leng),EN_EST_2(leng),EN_EST_3(leng),EN_EST_4(leng))
  ENDIF

  EN_EST_1 = 0.0_wp
  EN_EST_2 = 0.0_wp
  EN_EST_3 = 0.0_wp
  EN_EST_4 = 0
  le = 0
  
  IF(comm_world%mpirank==0) THEN
    WRITE(out,*) "CONSTANT NW SIMULATION; Length:", INFO%NUM_MC_STEPS
    WRITE(out,*)
  ENDIF

  DO i=count+1,count+INFO%NUM_MC_STEPS

    CALL SPAWNING_STEP(INFO,STRUC,hdes,SYM,dimen_a,dimen_b,a_str_arr,b_str_arr,Y_wa,Y_wb, &
                       W_list,HAM_DIAG,HAM_1,INT_2,EN%E_SHIFT )
    
    n_sp = n_sp + sp_contr
    !n_sp_p = n_sp_p + sp_contr_p
    !n_sp_m = n_sp_m + sp_contr_m

    CALL MERGE_ANNIHILATE(dimen_a,dimen_b,W_list)

    IF(mod(i,INFO%A_per)==0) THEN
      N_WALK_old = N_WALK
      N_WALK = NUM_WALKERS(dimen_a,dimen_b,SYM%SYM,SYM%SPIN,W_list)
      EN%E_SHIFT = EN%E_SHIFT - INFO%DAMP*LOG(REAL(N_WALK,wp)/REAL(N_WALK_old,wp))/(REAL(INFO%A_per,wp)*INFO%MC_tstep)

      !SAMPLING
      le = le + 1
      EN%EFCI_CORR = CALC_CORR(dimen_a,dimen_b,dim_list,W_list,en_list,a_adr_list,b_adr_list)
      EN_EST_1(le) = EN%E_SHIFT
      EN_EST_2(le) = EN%EFCI_CORR
    ENDIF

    IF(mod(i,INFO%PRINT_STEP)==0) THEN
      !CALL FLYVBJERG(le,EN_EST_1(1:le),mean_1,var_1)
      !CALL FLYVBJERG(le,EN_EST_2(1:le),mean_2,var_2)
    
      IF(comm_world%mpirank==0) THEN
        WRITE(out,*) "STEP = ", i , "NW =" , N_WALK, "W0 =", W_list(1,1)  
        WRITE(out,*) "  <S> = ", mean_1, "+-", var_1 
        WRITE(out,*) "  <E> = ", mean_2, "+-", var_2
        WRITE(out,*) "  S =", EN%E_SHIFT, "E =", EN%EFCI_CORR

        IF(ABS(var_1/mean_1)<0.001.and.ABS(var_2/mean_2)<0.0001.and.ABS(mean_1-mean_2)/mean_2<0.005) THEN
          WRITE(*,*) "Monte Carlo sampling converged"
          WRITE(out,*) "Monte Carlo sampling converged at the step: " , i
          EXIT
        ENDIF
      ENDIF
    ENDIF
  ENDDO 

  EN%EFCI_CORR = mean_2
  Norma = SQRT(REAL(SUM((W_list**2)),wp))
  N_WALK = NUM_WALKERS(dimen_a,dimen_b,SYM%SYM,SYM%SPIN,W_list)

  IF(comm_world%mpirank==0) THEN
    WRITE(out,*) 
    WRITE(out,*) "LEN_TOT  = " ,  count+INFO%NUM_MC_STEPS
    WRITE(out,*) "LEN_NW  = " , INFO%NUM_MC_STEPS
    WRITE(out,*) starshort
    WRITE(out,*) "<S>  =" , mean_1, "+-", var_1, ";", ABS(var_1/mean_1)*100.0
    WRITE(out,*) "<E>  =" , mean_2, "+-", var_2, ";", ABS(var_2/mean_2)*100.0
    WRITE(out,*) "Nw  =", N_WALK
    WRITE(out,*) "W0  =" , W_list(1,1) 
    WRITE(out,*) starshort
    WRITE(out,*) "Ratio spawned/tried:"
    WRITE(out,*) "Singles", REAL(sab)/REAL(sabcnt)
    WRITE(out,*) "Doubles ab", REAL(doab)/REAL(doabcnt)
    WRITE(out,*) "Doubles aa/bb", REAL(doaabb)/REAL(doaabbcnt)
    WRITE(out,*) starshort
    WRITE(out,*) "Acceptance Test:"
    WRITE(out,*) "Single acceptance:", (REAL(sab)/REAL(try_s))*100.0 , "%"
    WRITE(out,*) "Double ab acceptance:", (REAL(doab)/REAL(try_dab))*100.0 , "%"
    WRITE(out,*) "Double aa/bb acceptande:", (REAL(doaabb)/REAL(try_daabb))*100.0 , "%"
    WRITE(out,*) "Double acceptance:", (REAL(doab+doaabb)/REAL(try_dab+try_daabb))*100.0, "%"
    WRITE(out,*) starshort
    WRITE(out,*)"SPAWNED  =" , n_sp
    !WRITE(out,*)"+SPAWNED  =", n_sp_p
    !WRITE(out,*)"-SPAWNED  =", n_sp_m
    WRITE(out,*)"DIED  =" , n_died
    !WRITE(out,*)"+DIED  =", n_died_p
    !WRITE(out,*)"-DIED  =", n_died_m
    WRITE(out,*)"CLONED  =" , n_cl
    !WRITE(out,*)"+CLONED  =", n_cl_p
    !WRITE(out,*)"-CLONED  =", n_cl_m
    WRITE(out,*)"ANNIHILATED  =" , n_ann
    WRITE(out,*)starshort
    WRITE(out,*)"C(1,1)  =" , REAL(W_list(1,1),wp)/Norma
    WRITE(out,*)"C(2,1)  =" , REAL(W_list(2,1),wp)/Norma
    WRITE(out,*)"C(1,2)  =" , REAL(W_list(1,2),wp)/Norma
    WRITE(out,*)"C(2,2)  =" , REAL(W_list(2,2),wp)/Norma
    WRITE(out,*)"C(1,3)  =" , REAL(W_list(1,3),wp)/Norma
    WRITE(out,*)"C(3,1)  =" , REAL(W_list(3,1),wp)/Norma
    WRITE(out,*)"C(2,3)  =" , REAL(W_list(2,3),wp)/Norma
    WRITE(out,*)"C(3,2)  =" , REAL(W_list(3,2),wp)/Norma
    WRITE(out,*)"C(3,3)  =" , REAL(W_list(3,3),wp)/Norma
    CLOSE(UNIT=out)
  ENDIF

END SUBROUTINE QMCI_PROC

END MODULE FCIQMC
