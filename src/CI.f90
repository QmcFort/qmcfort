! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

MODULE CI
!********************************************************************************
!                                                                               
!  Module that contains subroutines and functions to perform CI calculation     
!     on HF groundstate. Aditionally the routines for transformation of         
!         molecular integrals from basis functions to spin orbitals             
!                        are insidide of this module                            
!                           Author: Zoran Sukurma                               
!                                                                                
!********************************************************************************

#include "preproc.inc"
USE CONSTANTS
use mpi, only: comm_world
use qmcfort_pos
use hamilton_vars, only: hamil_descriptor
USE STANDALONE
USE LAPACK

IMPLICIT NONE

CONTAINS 

SUBROUTINE TRAFO_HAM1_CI(ORB_NUM,HAM_1,INT_2,HAM_1T)
  INTEGER,INTENT(IN):: ORB_NUM
  REAL(wp),DIMENSION(ORB_NUM,ORB_NUM,INFO%ISPIN1),INTENT(IN):: HAM_1
  REAL(wp),DIMENSION(ORB_NUM,ORB_NUM,ORB_NUM,ORB_NUM,INFO%ISPIN2),INTENT(IN):: INT_2
  REAL(wp),DIMENSION(ORB_NUM,ORB_NUM,INFO%ISPIN1),INTENT(OUT)::  HAM_1T

  !local variables
  INTEGER:: i , loop

  HAM_1T = 0.0_wp
  
  DO loop=1,INFO%ISPIN1
    DO i=1,ORB_NUM
      HAM_1T(:,:,loop) = HAM_1T(:,:,loop) + INT_2(:,i,i,:,loop)
    ENDDO
  ENDDO

  HAM_1T = -0.5_wp * HAM_1T
  HAM_1T = HAM_1T + HAM_1

END SUBROUTINE TRAFO_HAM1_CI


!********************************************************************************
!       This recursive subroutine calculates the node weights of a graph        
!          that represents alpha (i.e. beta) strings. The arc weights           
!                 will be calculated from this weights later                    
!********************************************************************************
RECURSIVE FUNCTION W_WEIGHT(k,m,ORB_NUM,NELECT) RESULT(w)
  INTEGER:: k                                  !ordering number of electron
  INTEGER:: m                                  !number of electron up to k orbital
  INTEGER:: w                                  !weight
  INTEGER:: ORB_NUM                            !number od orbitals
  INTEGER:: NELECT                             !number of alpha(beta) electrons

  IF(k==0 .and. m==0) THEN
    w=1
  ELSEIF(k-m<0 .or. m<0 .or. k-m>ORB_NUM-NELECT .or. m>NELECT) THEN
    w=0
  ELSE
    w=W_WEIGHT(k-1,m,ORB_NUM,NELECT) + W_WEIGHT(k-1,m-1,ORB_NUM,NELECT)
  ENDIF
END FUNCTION W_WEIGHT


!********************************************************************************
!        Store the wertex weight in a matrix that represents a graph            
!********************************************************************************
SUBROUTINE CALC_VERTEX_WEIGHTS(W_w,ORB_NUM,NELECT)
  !Stores vertex weights for alpha(beta) stings into a matrix
  INTEGER:: ORB_NUM , NELECT                   !number of orbitals and alpha(beta) electrons
  INTEGER,DIMENSION(0:ORB_NUM,0:NELECT) :: W_w

  !local varaibles
  INTEGER:: k,m                                !loop indices

  DO m=0,NELECT
    DO k=0,ORB_NUM
      W_w(k,m) = W_WEIGHT(k,m,ORB_NUM,NELECT)
    ENDDO
  ENDDO
END SUBROUTINE CALC_VERTEX_WEIGHTS


!********************************************************************************
!       Suborutine that calculates arc weights of each arc in a graph           
!       For vertiacl arcs (no electron in that orbital) - weight is 0           
!      For horizontal arcs - sum of W_Weights from 2 neighboring nodes          
!********************************************************************************
SUBROUTINE CALC_ARC_WEIGHTS(Y_w,W_w,ORB_NUM,NELECT)
  INTEGER,INTENT(IN):: ORB_NUM
  INTEGER,INTENT(IN):: NELECT
  INTEGER,DIMENSION(0:ORB_NUM,0:NELECT),INTENT(IN):: W_w
  INTEGER,DIMENSION(0:1,1:ORB_NUM,1:NELECT),INTENT(OUT):: Y_w

  Y_w(0,:,:) = 0
  Y_w(1,:,:) = W_w(0:ORB_NUM-1,1:NELECT)
END SUBROUTINE CALC_ARC_WEIGHTS


FUNCTION CALC_ADRESS_naive(ORB_NUM,NELECT,OCCUP,Y_w) RESULT(ADRESS)
  !calculate adress of the string
  !each string is represented with an array of integers &
  !of length ORB_NUM with the m-th electron in k-th orbital
  INTEGER,INTENT(IN):: ORB_NUM , NELECT
  INTEGER,DIMENSION(ORB_NUM),INTENT(IN):: OCCUP
  INTEGER,DIMENSION(0:1,ORB_NUM,NELECT),INTENT(IN):: Y_w
  INTEGER:: ADRESS

  !local variables
  INTEGER:: k

  ADRESS = 1 + Y_w(OCCUP(1),1,OCCUP(1)) 
  DO k=2,ORB_NUM
    ADRESS = ADRESS + Y_w(OCCUP(k)-OCCUP(k-1),k,OCCUP(k))
  ENDDO
END FUNCTION CALC_ADRESS_naive


!********************************************************************************
!       Sum up all arc weights of a path in graph and determines the            
!                 ordering number of a alpha (beta) string                      
!********************************************************************************
FUNCTION CALC_ADRESS_BINARY(ORB_NUM,NELECT,OCCUP,Y_w) RESULT(ADRESS)
  !calculate adress of the string 
  !the string is represented as integer(kind=8)
  INTEGER,INTENT(IN):: ORB_NUM , NELECT
  INTEGER(KIND=intkind):: OCCUP
  INTEGER,DIMENSION(0:1,ORB_NUM,NELECT),INTENT(IN):: Y_w
  INTEGER:: ADRESS

  !local variables
  INTEGER:: k
  INTEGER:: m  

  ADRESS = 1 
  m=0
  DO k=1,ORB_NUM
    IF (btest(OCCUP,k-1)) THEN
      m = m + 1
      ADRESS = ADRESS + Y_w(1,k,m)
    ENDIF
  ENDDO
END FUNCTION CALC_ADRESS_BINARY


!********************************************************************************
!             Initialize string genrator - array of strings                     
!********************************************************************************
SUBROUTINE SET_STR_GEN(str_arr,myint,ORB_NUM,NELECT,starter,ender,ind,counter,dimen)
  INTEGER,INTENT(IN)::ORB_NUM,NELECT
  INTEGER,INTENT(OUT)::starter,ender,ind,counter,dimen
  INTEGER(KIND=intkind),ALLOCATABLE,DIMENSION(:),INTENT(OUT):: str_arr
  INTEGER(KIND=intkind),DIMENSION(NELECT),INTENT(OUT):: myint

  !local variables
  INTEGER:: i
  REAL(wp):: temp

  starter = 1
  ender  = ORB_NUM
  ind = 1
  counter = 0

  !dimen = 1
  !DO i=ORB_NUM,ORB_NUM-NELECT+1,-1
  !  dimen = dimen * i
  !ENDDO
  !
  !DO i=1,NELECT
  !  dimen = dimen / i
  !ENDDO
  
  temp=1.0_wp
  DO i=1,NELECT
    temp=temp *real(ORB_NUM-NELECT+i) / real(NELECT+1-i)
  ENDDO
  dimen=int(temp)
  myint = 0
  ALLOCATE(str_arr(dimen))
END SUBROUTINE SET_STR_GEN


!********************************************************************************
!        Recursive subroutine that generate all alpha(beta) strings in           
!                            reverse lexical ordering                           
!             Each sting is represented as integer of kind intkind              
!********************************************************************************
RECURSIVE SUBROUTINE GEN_STRINGS(str_arr,myint,NELECT,starter,ender,ind,counter,dimen)
  INTEGER:: NELECT,starter,ender,ind,counter,dimen
  INTEGER(KIND=intkind),DIMENSION(dimen):: str_arr
  INTEGER(KIND=intkind),DIMENSION(NELECT):: myint

  !local variables
  INTEGER:: i  
  INTEGER(KIND=intkind)::nullint = 0

  !conditon for return from the subroutine
  IF(ind==NELECT+1) THEN
    counter = counter + 1 
    str_arr(counter) = sum(myint)
    return
  ENDIF 
  
  i = starter
  DO WHILE (i<=ender .and. ender-i>=NELECT-ind)
    myint(ind) = ibset(nullint,i-1)
    CALL GEN_STRINGS(str_arr,myint,NELECT,starter+1,ender,ind+1,counter,dimen)
    i=i+1
    starter = starter + 1
  ENDDO 
END SUBROUTINE GEN_STRINGS


!********************************************************************************
!        String assembeler = initialize generator and generate stings           
!********************************************************************************
SUBROUTINE STRING_ASSEMBLER(ORB_NUM,NELECT,dimen,Y_w,str_arr)
  INTEGER:: ORB_NUM , NELECT
  INTEGER:: dimen
  INTEGER,DIMENSION(0:1,ORB_NUM,NELECT):: Y_w
  INTEGER(KIND=intkind),DIMENSION(:),ALLOCATABLE:: temp_str_arr , str_arr

  !local variables
  INTEGER(KIND=intkind),DIMENSION(NELECT):: myint
  INTEGER:: starter,ender,ind,counter
  INTEGER:: i 
  INTEGER(KIND=intkind):: str
  INTEGER:: adr

  CALL SET_STR_GEN(temp_str_arr,myint,ORB_NUM,NELECT,starter,ender,ind,counter,dimen)

  CALL GEN_STRINGS(temp_str_arr,myint,NELECT,starter,ender,ind,counter,dimen)

  ALLOCATE(str_arr(dimen))

  DO i=1,dimen
    str = temp_str_arr(i)
    adr = CALC_ADRESS_BINARY(ORB_NUM,NELECT,str,Y_w)
    str_arr(adr) = str
  ENDDO  

END SUBROUTINE


!********************************************************************************
!         Get occupation number of an orbital in one alpha (beta) string        
!                                  either 0 or 1                                
!********************************************************************************
FUNCTION GET_OCC(bit_test) RESULT(occ)
  LOGICAL :: bit_test
  INTEGER(KIND=1):: occ
  
  IF (bit_test) THEN
    occ = 1
  ELSE 
    occ = 0
  ENDIF
END FUNCTION GET_OCC


!********************************************************************************
!         Calculates diagonal  elements of a CI Matrix  <Phi|H|Phi>             
!               and stores elements in a matrix D(alpha,beta)                   
!     If alpha nd beta strings are indentical, the matrix is symmetric          
!********************************************************************************
SUBROUTINE CALC_CI_DIAG_HAM(a_str_arr,b_str_arr,dimen_a,dimen_b,HAM_1,INT_2,ORB_NUM, &
                            NELECT_a,NELECT_b,Y_wa,Y_wb,SYM,CI_H_DIAG)
  INTEGER,INTENT(IN):: dimen_a , dimen_b    !number of alpha and beta strings
  INTEGER(KIND=intkind),DIMENSION(dimen_a),INTENT(IN):: a_str_arr
  INTEGER(KIND=intkind),DIMENSION(dimen_b),INTENT(IN):: b_str_arr
  INTEGER,INTENT(IN):: ORB_NUM , NELECT_a , NELECT_b
  INTEGER,DIMENSION(0:1,ORB_NUM,NELECT_a),INTENT(IN):: Y_wa 
  INTEGER,DIMENSION(0:1,ORB_NUM,NELECT_b),INTENT(IN):: Y_wb
  INTEGER,INTENT(INOUT):: SYM              !sym=0 - diff strings ; sym=1 - eqv strings
  REAL(wp),DIMENSION(ORB_NUM,ORB_NUM,INFO%ISPIN1),INTENT(IN) :: HAM_1
  REAL(wp),DIMENSION(ORB_NUM,ORB_NUM,ORB_NUM,ORB_NUM,INFO%ISPIN2),INTENT(IN) :: INT_2
  REAL(wp),DIMENSION(dimen_a,dimen_b),INTENT(OUT):: CI_H_DIAG

  !local variables
  INTEGER(KIND=intkind):: a_str , b_str
  INTEGER:: p,q,i,j
  INTEGER:: a_adr , b_adr
  REAL(wp):: tempval
  INTEGER(KIND=1):: nia , nib , nja , njb
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
      WRITE(*,*) "Unrecognised ISPIN value"
    ENDIF
    CALL EXIT
  ENDIF

  IF (SYM==0) THEN    ! alpha and beta strings are not indentical
    !loop over beta strings
    DO p=1,dimen_b
      b_str = b_str_arr(p)
      b_adr = CALC_ADRESS_BINARY(ORB_NUM,NELECT_b,b_str,Y_wb)
      !loop over alpha strings
      DO q=1,dimen_a
        a_str = a_str_arr(q)
        a_adr = CALC_ADRESS_BINARY(ORB_NUM,NELECT_a,a_str,Y_wa)
        tempval = 0.0_wp
        DO j=1,ORB_NUM
          nja = GET_OCC(btest(a_str,j-1))
          njb = GET_OCC(btest(b_str,j-1))
          tempval = tempval + nja*HAM_1(j,j,inda) + njb*HAM_1(j,j,indb)
          tempval = tempval + nja*njb*INT_2(j,j,j,j,indab)
          DO i=j+1,ORB_NUM
            nia = GET_OCC(btest(a_str,i-1))
            nib = GET_OCC(btest(b_str,i-1))

            tempval = tempval + nia*nja*INT_2(i,i,j,j,inda) + nib*njb*INT_2(i,i,j,j,indb) &
                    + (nia*njb + nib*nja)*INT_2(i,i,j,j,indab) & 
                    - nia*nja*INT_2(i,j,j,i,inda) - nib*njb*INT_2(i,j,j,i,indb)                    
          ENDDO
        ENDDO

        CI_H_DIAG(a_adr,b_adr) = tempval
      ENDDO
    ENDDO

  ELSEIF (SYM==1) THEN  !alpha and beta strings are indentical
    !loop over beta strings
    DO p=1,dimen_b
      b_str = b_str_arr(p)
      b_adr = CALC_ADRESS_BINARY(ORB_NUM,NELECT_b,b_str,Y_wb)
      !loop over alpha strings
      DO q=p,dimen_a
        a_str = a_str_arr(q)
        a_adr = CALC_ADRESS_BINARY(ORB_NUM,NELECT_a,a_str,Y_wa)
        tempval = 0.0_wp
        DO j=1,ORB_NUM
          nja = GET_OCC(btest(a_str,j-1))
          njb = GET_OCC(btest(b_str,j-1))
          tempval = tempval + nja*HAM_1(j,j,inda) + njb*HAM_1(j,j,indb)
          tempval = tempval + nja*njb*INT_2(j,j,j,j,indab)
          DO i=j+1,ORB_NUM
            nia = GET_OCC(btest(a_str,i-1))
            nib = GET_OCC(btest(b_str,i-1))

            tempval = tempval + nia*nja*INT_2(i,i,j,j,inda) + nib*njb*INT_2(i,i,j,j,indb) &
                    + (nia*njb + nib*nja)*INT_2(i,i,j,j,indab) & 
                    - nia*nja*INT_2(i,j,j,i,inda) - nib*njb*INT_2(i,j,j,i,indb)                    
          ENDDO
        ENDDO

        CI_H_DIAG(a_adr,b_adr) = tempval
        CI_H_DIAG(b_adr,a_adr) = tempval
      ENDDO
    ENDDO
  ELSE
    IF(comm_world%mpirank==0) THEN
      WRITE(*,*) "An Error is occured with the value of var SYM"
    ENDIF
  ENDIF

END SUBROUTINE CALC_CI_DIAG_HAM


!********************************************************************************
!      This subroutine updates a CI vector in a Steepest descent method         
!                       Corresponds to CI_DIAG_MODE = 1                        
!********************************************************************************
SUBROUTINE STEEPEST_DESCENT_STEP(dimen_a,dimen_b,CI_vector,sigma,CI_H_DIAG,CI_eigval,CI_norm,alpha)
  INTEGER,INTENT(IN):: dimen_a , dimen_b
  REAL(wp),DIMENSION(dimen_a,dimen_b),INTENT(IN):: sigma 
  REAL(wp),DIMENSION(dimen_a,dimen_b),INTENT(IN):: CI_H_DIAG 
  REAL(wp),DIMENSION(dimen_a,dimen_b),INTENT(INOUT):: CI_VECTOR
  REAL(wp),INTENT(OUT):: CI_eigval , CI_norm
  REAL(wp),INTENT(IN):: alpha
  
  !local variables
  REAL(wp),DIMENSION(dimen_a,dimen_b) :: correction

  CI_eigval = SUM(CI_vector*sigma)
  correction = sigma - CI_eigval*CI_vector
  correction = -correction/(CI_H_DIAG-CI_eigval-2.0_wp*CI_vector*correction+alpha*CI_vector**2)
  
  !IF(ABS(CI_EIGVAL-CI_H_DIAG(1,1))<1.e-08_wp) THEN
  !  correction(1,1) = 0.0_wp
  !ENDIF 
  !WRITE(*,*)"correction " , correction
  CI_norm = NORM2(correction)
  CI_vector = CI_vector + correction
  CI_vector = CI_vector/NORM2(CI_vector)

END SUBROUTINE STEEPEST_DESCENT_STEP


!here come routines for the step sigma=HC
!********************************************************************************
!         Function that calculates the phase factor <string|E(j,i)|string>       
!                  by an excitation, where E(j,i)=a_cr(j)a_in(i)                
!********************************************************************************
FUNCTION CI_PHASE(string,i,j) RESULT(ph)
  INTEGER(KIND=intkind)::string
  INTEGER:: i,j

  !local variables
  INTEGER:: counter
  INTEGER:: ph
  INTEGER:: k

  counter=0
  IF(j>i) THEN
    DO k=i+1,j-1
      IF(btest(string,k-1)) THEN
        counter = counter + 1
      ENDIF
    ENDDO
  ELSEIF(j<i) THEN
    DO k=j+1,i-1
      IF(btest(string,k-1)) THEN
        counter = counter + 1
      ENDIF
    ENDDO
  ELSE
    counter = 0
  ENDIF

  ph = (-1)**counter

END FUNCTION CI_PHASE


!********************************************************************************
!         Function that calculates the phase factor <string|E(j,i)|string>       
!                  by an excitation, where E(j,i)=a_cr(j)a_in(i)                
!********************************************************************************
FUNCTION CI_PHASE_N(string1,string2,i,j) RESULT(ph)
  INTEGER(KIND=intkind)::string1 , string2
  INTEGER:: i,j

  !local variables
  INTEGER:: counter
  INTEGER:: ph
  INTEGER:: k

  counter=0
  IF(j>i) THEN
    DO k=i+1,j-1
      IF(btest(string1,k-1) .and. btest(string2,k-1)) THEN
        counter = counter + 1
      ENDIF
    ENDDO
  ELSEIF(j<i) THEN
    DO k=j+1,i-1
      IF(btest(string1,k-1) .and. btest(string2,k-1)) THEN
        counter = counter + 1
      ENDIF
    ENDDO
  ELSE
    counter = 0
  ENDIF

  ph = (-1)**counter

END FUNCTION CI_PHASE_N


SUBROUTINE SINGLE_EXC_EN(ORB_NUM,HAM_1,INT_2,a_str,b_str,occ1,unocc1,EN)
!**************************************************
!    Calculates energy of a single excitation     *
!      Order of a_str and b_str is important      *
!**************************************************
  INTEGER,INTENT(IN):: ORB_NUM
  REAL(wp),DIMENSION(ORB_NUM,ORB_NUM):: HAM_1
  REAL(wp),DIMENSION(ORB_NUM,ORB_NUM,ORB_NUM,ORB_NUM,INFO%ISPIN2):: INT_2
  INTEGER(KIND=intkind),INTENT(IN):: a_str,b_str
  INTEGER,INTENT(IN):: occ1 , unocc1
  REAL(wp),INTENT(OUT):: EN

  !local variables
  INTEGER:: i
  INTEGER:: na , nb
  INTEGER:: ph
  INTEGER:: inda , indb , indab

  IF(INFO%ISPIN1==1) THEN
    inda = 1 
    indb = 1 
    indab = 1
  ELSEIF(INFO%ISPIN1==2) THEN
    inda = 1 
    indb = 2
    indab = 3
  ELSE
    IF(comm_world%mpirank==0) THEN
      WRITE(*,*) "Unrecognized ISPIN value"
    ENDIF
    CALL EXIT
  ENDIF

  EN = HAM_1(unocc1,occ1)

  DO i=1,ORB_NUM
    na = GET_OCC(btest(a_str,i-1))
    nb = GET_OCC(btest(b_str,i-1))

    EN = EN + na*INT_2(occ1,unocc1,i,i,inda) + nb*INT_2(occ1,unocc1,i,i,indb) &
            - na*INT_2(occ1,i,i,unocc1,inda)
  ENDDO

  ph = CI_PHASE(a_str,occ1,unocc1)

  EN = EN * ph

END SUBROUTINE SINGLE_EXC_EN


SUBROUTINE DOUBLE_EXC_EN(ORB_NUM,INT_2,a_str,b_str,occ1,occ2,unocc1,unocc2,EN,mode)
!**************************************************
!    Calculates energy of a double excitation     *
!**************************************************
  INTEGER,INTENT(IN):: ORB_NUM
  INTEGER,INTENT(IN):: occ1 , occ2
  INTEGER(KIND=intkind),INTENT(IN):: a_str ,  b_str
  INTEGER,INTENT(IN):: unocc1 , unocc2
  REAL(wp),DIMENSION(ORB_NUM,ORB_NUM,ORB_NUM,ORB_NUM),INTENT(IN):: INT_2
  REAL(wp),INTENT(OUT):: EN
  INTEGER,INTENT(IN):: mode

  !local variables
  INTEGER:: ph
  INTEGER(KIND=intkind):: tempstr

  IF(mode==1) THEN 
  !ALPHA-ALPHA or BETA-BETA EXCITATION
    EN = INT_2(occ1,unocc1,occ2,unocc2) - INT_2(occ1,unocc2,occ2,unocc1)
    
    tempstr = ibset(ibclr(a_str,occ1-1),unocc1-1)

    ph = CI_PHASE(a_str,occ1,unocc1) * CI_PHASE(tempstr,occ2,unocc2)

    EN = EN * ph     
  ELSEIF(mode==2) THEN
  !ALPHA-BETA EXCITATION
    EN = INT_2(occ1,unocc1,occ2,unocc2)

    ph = CI_PHASE(a_str,occ1,unocc1) * CI_PHASE(b_str,occ2,unocc2)

    EN = EN * ph

  ENDIF

END SUBROUTINE DOUBLE_EXC_EN


SUBROUTINE MAKE_EXC_LIST(str,dimen,ORB_NUM,dimen1,dimen2,exc1_list,exc2_list,ind1_list,ind2_list)
  INTEGER,INTENT(IN):: dimen,ORB_NUM
  INTEGER(KIND=intkind):: str
  INTEGER,INTENT(OUT):: dimen1,dimen2
  INTEGER(KIND=intkind),DIMENSION(dimen),INTENT(INOUT):: exc1_list
  INTEGER(KIND=intkind),DIMENSION(dimen**2),INTENT(INOUT):: exc2_list
  INTEGER,DIMENSION(2,dimen),INTENT(INOUT):: ind1_list
  INTEGER,DIMENSION(4,dimen**2),INTENT(INOUT):: ind2_list

  !local variables
  INTEGER:: i,j
  INTEGER:: k,l
  INTEGER(KIND=intkind)::tempstr , temp2str , temp3str , temp4str

  dimen1 = 0
  dimen2 = 0

  DO l=1,ORB_NUM
    IF(btest(str,l-1)) THEN
      tempstr = ibclr(str,l-1)
      DO k=1,ORB_NUM
        IF(.not. btest(str,k-1)) THEN
          dimen1  = dimen1 + 1
          temp2str =  ibset(tempstr,k-1)
          exc1_list(dimen1) = temp2str
          ind1_list(1,dimen1) = l
          ind1_list(2,dimen1) = k
          DO j=l+1,ORB_NUM
            IF(btest(str,j-1)) THEN
              IF(j/=k) THEN
                temp3str=ibclr(temp2str,j-1)

                DO i=k+1,ORB_NUM
                  IF(i/=l .and. i/=j) THEN
                    dimen2 = dimen2 + 1
                    temp4str = ibset(temp3str,i-1)
                    exc2_list(dimen2) = temp4str
                    ind2_list(1,dimen2) = l
                    ind2_list(2,dimen2) = j
                    ind2_list(3,dimen2) = k
                    ind2_list(4,dimen2) = i
                  ENDIF
                ENDDO
              ENDIF
            ENDIF
          ENDDO
        ENDIF
      ENDDO
    ENDIF
  ENDDO

END SUBROUTINE MAKE_EXC_LIST


SUBROUTINE CALC_CI_MATRIX(INFO,STRUC,hdes,dimen_a,dimen_b,a_str_arr,b_str_arr,Y_wa,Y_wb,HAM_1,INT_2,CI_MAT)
  TYPE(CONTROLL):: INFO
  TYPE(STRUCTURE):: STRUC
  type(hamil_descriptor) :: hdes
  INTEGER,INTENT(IN):: dimen_a,dimen_b
  INTEGER(KIND=intkind),DIMENSION(dimen_a),INTENT(IN):: a_str_arr
  INTEGER(KIND=intkind),DIMENSION(dimen_b),INTENT(IN):: b_str_arr
  INTEGER,DIMENSION(0:1,hdes%n,hdes%nel(1)),INTENT(IN):: Y_wa  
  INTEGER,DIMENSION(0:1,hdes%n,hdes%nel(2)),INTENT(IN):: Y_wb
  REAL(wp),DIMENSION(hdes%n,hdes%n,INFO%ISPIN1),INTENT(IN):: HAM_1
  REAL(wp),DIMENSION(hdes%n,hdes%n,hdes%n,hdes%n,INFO%ISPIN2),INTENT(IN):: INT_2
  REAL(wp),DIMENSION(dimen_a*dimen_b,dimen_a*dimen_b),INTENT(OUT):: CI_MAT

  !local variables
  INTEGER:: a , b , i , j
  INTEGER(KIND=intkind):: a_str , b_str
  INTEGER:: dim_a1 , dim_a2 , dim_b1 , dim_b2
  INTEGER(KIND=intkind),DIMENSION(dimen_a):: exc_a1
  INTEGER(KIND=intkind),DIMENSION(dimen_b):: exc_b1
  INTEGER(KIND=intkind),DIMENSION(dimen_a**2):: exc_a2
  INTEGER(KIND=intkind),DIMENSION(dimen_b**2):: exc_b2
  INTEGER,DIMENSION(2,dimen_a):: indl1a
  INTEGER,DIMENSION(2,dimen_b):: indl1b
  INTEGER,DIMENSION(4,dimen_a**2):: indl2a
  INTEGER,DIMENSION(4,dimen_b**2):: indl2b

  INTEGER(KIND=intkind):: ea1 , ea2 , eb1 , eb2
  INTEGER:: adr1 , adr2 , eadr1 , eadr2, ebdr1 , ebdr2
  INTEGER:: ind1 , ind2 , ind3 , ind4
  REAL(wp),DIMENSION(dimen_a,dimen_b):: CI_H_DIAG
  REAL(wp):: matel
  INTEGER:: inda, indb , indab , SYM

  IF(INFO%ISPIN1==1) THEN
    inda = 1 
    indb = 1 
    indab = 1
    IF(hdes%nel(1)==hdes%nel(2)) THEN
      SYM = 1
    ELSE
      SYM = 0
    ENDIF
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

  !Initialize CI Matrix
  CI_MAT = 0.0_wp
 
  !CALCULATION OF OFF-DIAGONAL CI MATRIX ELEMENTS
  
  DO a=1,dimen_a
    a_str = a_str_arr(a)
    CALL MAKE_EXC_LIST(a_str,dimen_a,hdes%n,dim_a1,dim_a2,exc_a1,exc_a2,indl1a,indl2a)
 
    !1. Calculate <Phib|<Phia|H|Phia+>|Phib>
    DO i=1,dim_a1
      ea1 = exc_a1(i)
      eadr1 = CALC_ADRESS_BINARY(hdes%n,hdes%nel(1),ea1,Y_wa)
      DO b=1,dimen_b
        b_str = b_str_arr(b)

        adr1 = (a-1)*dimen_b + b
        adr2 = (eadr1-1)*dimen_b + b 

        IF(adr1>adr2) THEN
          ind1 = indl1a(1,i)
          ind2 = indl1a(2,i)

          CALL SINGLE_EXC_EN(hdes%n,HAM_1(:,:,inda),INT_2,a_str,b_str,ind1,ind2,matel)
          
          CI_MAT(adr1,adr2) = matel
        ENDIF
      ENDDO
    ENDDO

    !2. Calculate <Phib|<Phia|H|Phia++>|Phib>
    DO i=1,dim_a2
      ea2 = exc_a2(i)
      eadr2 = CALC_ADRESS_BINARY(hdes%n,hdes%nel(1),ea2,Y_wa)
      DO b=1,dimen_b
        b_str = b_str_arr(b)

        adr1 = (a-1)*dimen_b + b
        adr2 = (eadr2-1)*dimen_b + b 
        
        IF(adr1>adr2) THEN
          ind1 = indl2a(1,i)
          ind2 = indl2a(2,i)
          ind3 = indl2a(3,i)
          ind4 = indl2a(4,i)

          CALL DOUBLE_EXC_EN(hdes%n,INT_2(:,:,:,:,inda),a_str,b_str,ind1,ind2,ind3,ind4,matel,1)

          CI_MAT(adr1,adr2) = matel
        ENDIF
      ENDDO
    ENDDO    
  ENDDO


  DO b=1,dimen_b
    b_str = b_str_arr(b)
    CALL MAKE_EXC_LIST(b_str,dimen_b,hdes%n,dim_b1,dim_b2,exc_b1,exc_b2,indl1b,indl2b)
    
    !3. Calculate <Phib|<Phia|H|Phia>|Phib+>
    DO i=1,dim_b1
      eb1 = exc_b1(i)
      ebdr1 = CALC_ADRESS_BINARY(hdes%n,hdes%nel(2),eb1,Y_wb)
      DO a=1,dimen_a
        a_str= a_str_arr(a)

        adr1 = (a-1)*dimen_b + b
        adr2 = (a-1)*dimen_b + ebdr1

        IF(adr1>adr2) THEN
          ind1 = indl1b(1,i)
          ind2 = indl1b(2,i)

          CALL SINGLE_EXC_EN(hdes%n,HAM_1(:,:,indb),INT_2,b_str,a_str,ind1,ind2,matel)
          
          CI_MAT(adr1,adr2) = matel
        ENDIF
      ENDDO
    ENDDO

    !4. Calculate <Phib|<Phia|H|Phia>|Phib++>
    DO i=1,dim_b2
      eb2 = exc_b2(i)
      ebdr2 = CALC_ADRESS_BINARY(hdes%n,hdes%nel(2),eb2,Y_wb)
      DO a=1,dimen_a
        a_str = a_str_arr(a)

        adr1 = (a-1)*dimen_b + b
        adr2 = (a-1)*dimen_b + ebdr2

        IF(adr1>adr2) THEN
          ind1 = indl2b(1,i)
          ind2 = indl2b(2,i)
          ind3 = indl2b(3,i)
          ind4 = indl2b(4,i)

          CALL DOUBLE_EXC_EN(hdes%n,INT_2(:,:,:,:,indb),b_str,a_str,ind1,ind2,ind3,ind4,matel,1)

          CI_MAT(adr1,adr2) = matel
        ENDIF
      ENDDO
    ENDDO    
  ENDDO

  !5. Calculation of <Phib|<Phia|H|Phia+>|Phib+>
  DO a=1,dimen_a
    a_str = a_str_arr(a)
    CALL MAKE_EXC_LIST(a_str,dimen_a,hdes%n,dim_a1,dim_a2,exc_a1,exc_a2,indl1a,indl2a)
    DO i=1,dim_a1
      ea1 = exc_a1(i)
      eadr1 = CALC_ADRESS_BINARY(hdes%n,hdes%nel(1),ea1,Y_wa)
      DO b=1,dimen_b
        b_str = b_str_arr(b)
        CALL MAKE_EXC_LIST(b_str,dimen_b,hdes%n,dim_b1,dim_b2,exc_b1,exc_b2,indl1b,indl2b)
        DO j=1,dim_b1
          eb1 = exc_b1(j)
          ebdr1 = CALC_ADRESS_BINARY(hdes%n,hdes%nel(2),eb1,Y_wb) 

          adr1 = (a-1)*dimen_b + b
          adr2 = (eadr1-1)*dimen_b + ebdr1

          IF(adr1>adr2) THEN
            ind1 = indl1a(1,i)
            ind2 = indl1b(1,j)
            ind3 = indl1a(2,i)
            ind4 = indl1b(2,j)

            CALL DOUBLE_EXC_EN(hdes%n,INT_2(:,:,:,:,indab),a_str,b_str,ind1,ind2,ind3,ind4,matel,2)

            CI_MAT(adr1,adr2) = matel
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  CI_MAT = CI_MAT + TRANSPOSE(CI_MAT)

  !Calcuate diagonal elements of CI Matrix
  CALL CALC_CI_DIAG_HAM(a_str_arr,b_str_arr,dimen_a,dimen_b,HAM_1,INT_2,hdes%n,&
                        hdes%nel(1),hdes%nel(2),Y_wa,Y_wb,SYM,CI_H_DIAG)

  !Reassemble CI_H_DIAG into CI_MAT
  DO a=1,dimen_a
    DO b=1,dimen_b
      adr1 = (a-1)*dimen_b + b
      !WRITE(*,*) "  " , a , b , CI_H_DIAG(a,b)
      CI_MAT(adr1,adr1) = CI_H_DIAG(a,b)
    ENDDO
  ENDDO 

END SUBROUTINE CALC_CI_MATRIX


SUBROUTINE BETA_BETA_ROUTINE_VEC(ORB_NUM,NELECT_b,HAM_1,INT_2,dimen_a,dimen_b,Y_wb,b_str_arr,CI_vector,sigma)
  INTEGER,INTENT(IN):: ORB_NUM,NELECT_b,dimen_a,dimen_b
  INTEGER,DIMENSION(2,ORB_NUM,NELECT_b):: Y_wb
  INTEGER(KIND=intkind),DIMENSION(dimen_b),INTENT(IN):: b_str_arr
  REAL(wp),DIMENSION(dimen_a,dimen_b),INTENT(IN):: CI_vector 
  REAL(wp),DIMENSION(ORB_NUM,ORB_NUM),INTENT(IN):: HAM_1
  REAL(wp),DIMENSION(ORB_NUM,ORB_NUM,ORB_NUM,ORB_NUM),INTENT(IN):: INT_2  
  REAL(wp),DIMENSION(dimen_a,dimen_b),INTENT(INOUT):: sigma

  !local variables
  INTEGER:: p
  INTEGER:: i,j
  INTEGER:: k,l
  INTEGER:: adr1 , adr2
  INTEGER:: phase1 , phase2
  INTEGER(KIND=intkind):: str, tempstr , temp2str , temp3str , temp4str
  REAL(wp),DIMENSION(dimen_b):: temparray

  DO p=1,dimen_b                  !loop over all beta(alpha) strings

    str=b_str_arr(p)
    temparray = 0.0_wp

    DO l=1,ORB_NUM              !looop over all --
      IF(btest(str,l-1)) THEN
        tempstr = ibclr(str,l-1)
        DO k=1,ORB_NUM            !single excitations of str
          IF(.not. btest(tempstr,k-1)) THEN
            temp2str =  ibset(tempstr,k-1)
            phase1 = CI_PHASE(str,k,l)         
            !phase1 = CI_PHASE_N(str,temp2str,k,l)
            adr1 = CALC_ADRESS_BINARY(ORB_NUM,NELECT_b,temp2str,Y_wb) 
            temparray(adr1) = temparray(adr1) + phase1*HAM_1(k,l)
            DO j=1,ORB_NUM      !loop over all --
              IF (btest(temp2str,j-1)) THEN
                temp3str = ibclr(temp2str,j-1)
                DO i=1,ORB_NUM    !double excitations of str
                  IF(.not. btest(temp3str,i-1)) THEN
                      temp4str = ibset(temp3str,i-1)
                      phase2 = phase1 * CI_PHASE(temp2str,i,j)
                      !phase2 = phase1 * CI_PHASE_N(temp2str,temp4str,i,j)
                      adr2 = CALC_ADRESS_BINARY(ORB_NUM,NELECT_b,temp4str,Y_wb)
                      temparray(adr2) = temparray(adr2) + 0.5_wp*phase2*INT_2(i,j,k,l)
                  ENDIF
                ENDDO
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDIF
    ENDDO
 
    CALL GEMV("N",dimen_a,dimen_b,1.0_wp,CI_vector,dimen_a,temparray,1,1.0_wp,sigma(:,p),1)

  ENDDO

END SUBROUTINE BETA_BETA_ROUTINE_VEC


SUBROUTINE BETA_BETA_ROUTINE(ORB_NUM,NELECT_b,HAM_1,INT_2,dimen_a,dimen_b,Y_wb,b_str_arr,CI_vector,sigma)
  INTEGER,INTENT(IN):: ORB_NUM,NELECT_b,dimen_a,dimen_b
  INTEGER,DIMENSION(2,ORB_NUM,NELECT_b):: Y_wb
  INTEGER(KIND=intkind),DIMENSION(dimen_b),INTENT(IN):: b_str_arr
  REAL(wp),DIMENSION(dimen_a,dimen_b),INTENT(IN):: CI_vector 
  REAL(wp),DIMENSION(ORB_NUM,ORB_NUM),INTENT(IN):: HAM_1
  REAL(wp),DIMENSION(ORB_NUM,ORB_NUM,ORB_NUM,ORB_NUM),INTENT(IN):: INT_2  
  REAL(wp),DIMENSION(dimen_a,dimen_b),INTENT(INOUT):: sigma
 
  !local variables
  INTEGER:: p,q
  INTEGER:: i,j
  INTEGER:: k,l
  INTEGER:: adr1 , adr2
  INTEGER:: phase1 , phase2
  INTEGER(KIND=intkind):: str, tempstr , temp2str , temp3str , temp4str


  DO p=1,dimen_b                  !loop over all beta(alpha) strings

    str=b_str_arr(p)

    DO l=1,ORB_NUM              !looop over all --
      IF(btest(str,l-1)) THEN
        tempstr = ibclr(str,l-1)
        DO k=1,ORB_NUM            !single excitations of str
          IF(.not. btest(tempstr,k-1)) THEN
            temp2str =  ibset(tempstr,k-1)
            phase1 = CI_PHASE(str,k,l)
            !phase1 = CI_PHASE_N(str,temp2str,k,l)         
            adr1 = CALC_ADRESS_BINARY(ORB_NUM,NELECT_b,temp2str,Y_wb)
            
            DO q=1,dimen_a 
                sigma(q,p) = sigma(q,p) + phase1*HAM_1(k,l)*CI_vector(q,adr1)
            ENDDO

            DO j=1,ORB_NUM      !loop over all --
              IF (btest(temp2str,j-1)) THEN
                temp3str = ibclr(temp2str,j-1)
                DO i=1,ORB_NUM    !double excitations of str
                  IF(.not. btest(temp3str,i-1)) THEN
                    temp4str = ibset(temp3str,i-1)
                    phase2 = phase1 * CI_PHASE(temp2str,i,j)
                    !phase2 = phase1 * CI_PHASE_N(temp2str,temp4str,i,j)
                    adr2 = CALC_ADRESS_BINARY(ORB_NUM,NELECT_b,temp4str,Y_wb)
                   
                    DO q=1,dimen_a
                        sigma(q,p) = sigma(q,p) + 0.5_wp*phase2*INT_2(i,j,k,l)* &
                                        CI_vector(q,adr2)
                    ENDDO

                  ENDIF
                ENDDO
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDIF
    ENDDO
 
  ENDDO
END SUBROUTINE BETA_BETA_ROUTINE


SUBROUTINE ALPHA_ALPHA_ROUTINE(ORB_NUM,NELECT_a,HAM_1,INT_2,dimen_a,dimen_b,Y_wa,a_str_arr,CI_vector,sigma)
  INTEGER,INTENT(IN):: ORB_NUM,NELECT_a,dimen_a,dimen_b
  INTEGER,DIMENSION(2,ORB_NUM,NELECT_a):: Y_wa
  INTEGER(KIND=intkind),DIMENSION(dimen_a),INTENT(IN):: a_str_arr
  REAL(wp),DIMENSION(dimen_a,dimen_b),INTENT(IN):: CI_vector 
  REAL(wp),DIMENSION(ORB_NUM,ORB_NUM),INTENT(IN):: HAM_1
  REAL(wp),DIMENSION(ORB_NUM,ORB_NUM,ORB_NUM,ORB_NUM),INTENT(IN):: INT_2  
  REAL(wp),DIMENSION(dimen_a,dimen_b),INTENT(INOUT):: sigma

  !local variables
  INTEGER:: p,q
  INTEGER:: i,j
  INTEGER:: k,l
  INTEGER:: adr1 , adr2
  INTEGER:: phase1 , phase2
  INTEGER(KIND=intkind):: str, tempstr , temp2str , temp3str , temp4str

  DO p=1,dimen_a                  !loop over all beta(alpha) strings

    str=a_str_arr(p)

    DO l=1,ORB_NUM              !looop over all --
      IF(btest(str,l-1)) THEN
        tempstr = ibclr(str,l-1)
        DO k=1,ORB_NUM            !single excitations of str
          IF(.not. btest(tempstr,k-1)) THEN
            temp2str =  ibset(tempstr,k-1)
            phase1 = CI_PHASE(str,k,l)    
            !phase1 = CI_PHASE_N(str,temp2str,k,l)     
            adr1 = CALC_ADRESS_BINARY(ORB_NUM,NELECT_a,temp2str,Y_wa)
            
            DO q=1,dimen_b
                sigma(p,q) = sigma(p,q) + phase1*HAM_1(k,l)*CI_vector(adr1,q)
            ENDDO

            DO j=1,ORB_NUM      !loop over all --
              IF (btest(temp2str,j-1)) THEN
                temp3str = ibclr(temp2str,j-1)
                DO i=1,ORB_NUM    !double excitations of str
                  IF(.not. btest(temp3str,i-1)) THEN
                    temp4str = ibset(temp3str,i-1)
                    phase2 = phase1 * CI_PHASE(temp2str,i,j)
                    !phase2 = phase1 * CI_PHASE_N(temp2str,temp4str,i,j)
                    adr2 = CALC_ADRESS_BINARY(ORB_NUM,NELECT_a,temp4str,Y_wa)
                 
                    DO q=1,dimen_b
                        sigma(p,q) = sigma(p,q) + 0.5_wp*phase2*INT_2(i,j,k,l)* &
                                        CI_vector(adr2,q)
                    ENDDO

                  ENDIF
                ENDDO
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDIF
    ENDDO
 
  ENDDO

END SUBROUTINE ALPHA_ALPHA_ROUTINE


SUBROUTINE ALPHA_ALPHA_ROUTINE_VEC(ORB_NUM,NELECT_a,HAM_1,INT_2,dimen_a,dimen_b, &
                               Y_wa,a_str_arr,CI_vector,sigma)
  INTEGER,INTENT(IN):: ORB_NUM,NELECT_a,dimen_a,dimen_b
  INTEGER,DIMENSION(2,ORB_NUM,NELECT_a):: Y_wa
  INTEGER(KIND=intkind),DIMENSION(dimen_a),INTENT(IN):: a_str_arr
  REAL(wp),DIMENSION(dimen_a,dimen_b),INTENT(IN):: CI_vector
  REAL(wp),DIMENSION(ORB_NUM,ORB_NUM),INTENT(IN):: HAM_1
  REAL(wp),DIMENSION(ORB_NUM,ORB_NUM,ORB_NUM,ORB_NUM),INTENT(IN):: INT_2
  REAL(wp),DIMENSION(dimen_a,dimen_b),INTENT(INOUT):: sigma

  !local variables
  INTEGER:: p
  INTEGER:: i,j
  INTEGER:: k,l
  INTEGER:: adr1 , adr2
  INTEGER:: phase1 , phase2
  INTEGER(KIND=intkind):: str, tempstr , temp2str , temp3str , temp4str
  REAL(wp),DIMENSION(dimen_a):: temparray

  DO p=1,dimen_a                  !loop over all beta(alpha) strings

    str=a_str_arr(p)
    temparray = 0.0_wp

    DO l=1,ORB_NUM              !looop over all --
      IF(btest(str,l-1)) THEN   !all single excitations of str
        tempstr = ibclr(str,l-1)
        DO k=1,ORB_NUM          
          IF(.not. btest(tempstr,k-1)) THEN
            temp2str =  ibset(tempstr,k-1)
            phase1 = CI_PHASE(str,k,l)
            !phase1 = CI_PHASE_N(str,temp2str,k,l)
            adr1 = CALC_ADRESS_BINARY(ORB_NUM,NELECT_a,temp2str,Y_wa)
            temparray(adr1) = temparray(adr1) + phase1*HAM_1(k,l)
            DO j=1,ORB_NUM                 !loop over all --
              IF (btest(temp2str,j-1)) THEN       !double excitations of str
                temp3str = ibclr(temp2str,j-1)
                DO i=1,ORB_NUM  
                  IF(.not. btest(temp3str,i-1)) THEN
                      temp4str = ibset(temp3str,i-1)
                      phase2 = phase1 * CI_PHASE(temp2str,i,j)
                      !phase2 = phase1 * CI_PHASE_N(temp2str,temp4str,i,j)
                      adr2 = CALC_ADRESS_BINARY(ORB_NUM,NELECT_a,temp4str,Y_wa)
                      temparray(adr2) = temparray(adr2) + 0.5_wp*phase2*INT_2(i,j,k,l)
                  ENDIF
                ENDDO
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDIF
    ENDDO

    CALL GEMV("T",dimen_a,dimen_b,1.0_wp,CI_vector,dimen_a,temparray,1,1.0_wp,sigma(p,:),1)

  ENDDO

END SUBROUTINE ALPHA_ALPHA_ROUTINE_VEC


SUBROUTINE ALPHA_BETA_ROUTINE(ORB_NUM,NELECT_a,NELECT_b,INT_2,dimen_a,dimen_b, &
                              Y_wa,Y_wb,a_str_arr,b_str_arr,CI_vector,sigma,SYM,SPIN)

  INTEGER,INTENT(IN):: ORB_NUM,NELECT_a,NELECT_b,dimen_a,dimen_b
  INTEGER(KIND=intkind),DIMENSION(dimen_b),INTENT(IN)::a_str_arr,b_str_arr
  INTEGER,DIMENSION(2,ORB_NUM,NELECT_a):: Y_wa
  INTEGER,DIMENSION(2,ORB_NUM,NELECT_b):: Y_wb
  REAL(wp),DIMENSION(dimen_a,dimen_b),INTENT(IN):: CI_vector
  REAL(wp),DIMENSION(ORB_NUM,ORB_NUM,ORB_NUM,ORB_NUM),INTENT(IN):: INT_2
  REAL(wp),DIMENSION(dimen_a,dimen_b),INTENT(INOUT):: sigma
  INTEGER,INTENT(IN):: SYM , SPIN

  !local variables
  INTEGER:: p,q
  INTEGER:: i,j
  INTEGER:: k,l
  INTEGER:: a_adr1 , b_adr1
  INTEGER:: phase1 , phase2
  INTEGER(KIND=intkind):: a_str, b_str , a_tempstr , b_tempstr , a_temp2str , b_temp2str
  REAL(wp),DIMENSION(dimen_a,dimen_b):: tempsigma

  tempsigma = 0.0_wp

  IF(SYM==1) THEN

    DO q=1,dimen_b

      b_str = b_str_arr(q)

      DO l=1,ORB_NUM
        IF(btest(b_str,l-1)) THEN
          b_tempstr = ibclr(b_str,l-1)
          DO k=1,ORB_NUM
            IF(.not. btest(b_tempstr,k-1)) THEN
              b_temp2str = ibset(b_tempstr,k-1)
              b_adr1 = CALC_ADRESS_BINARY(ORB_NUM,NELECT_b,b_temp2str,Y_wb)
              phase1 = CI_PHASE(b_str,k,l)
              !phase1 = CI_PHASE_N(b_str,b_temp2str,k,l)
 
              DO p=q+1,dimen_a

                a_str = a_str_arr(p)

                DO j=1,ORB_NUM
                  IF(btest(a_str,j-1)) THEN
                    a_tempstr = ibclr(a_str,j-1)
                    DO i=1,ORB_NUM
                      IF(.NOT. btest(a_tempstr,i-1)) THEN
                        a_temp2str = ibset(a_tempstr,i-1)
                        a_adr1 = CALC_ADRESS_BINARY(ORB_NUM,NELECT_a,a_temp2str,Y_wa)
                        phase2 = phase1 *  CI_PHASE(a_str,i,j)
                        !phase2 = phase1 * CI_PHASE_N(a_str,a_temp2str,i,j)
                        
                        tempsigma(p,q) = tempsigma(p,q) + &
                        phase2*INT_2(i,j,k,l)*CI_vector(a_adr1,b_adr1)

                      ENDIF
                    ENDDO
                  ENDIF
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDIF
      ENDDO
    ENDDO

    tempsigma = tempsigma + (-1)**SPIN*TRANSPOSE(tempsigma)

  
    !add diagonal elements
    DO q=1,dimen_b
      b_str = b_str_arr(q)
 
      DO l=1,ORB_NUM
        IF(btest(b_str,l-1)) THEN
          b_tempstr = ibclr(b_str,l-1)
          DO k=1,ORB_NUM
            IF(.not. btest(b_tempstr,k-1)) THEN
              b_temp2str = ibset(b_tempstr,k-1)
              b_adr1 = CALC_ADRESS_BINARY(ORB_NUM,NELECT_b,b_temp2str,Y_wb)
              phase1 = CI_PHASE(b_str,k,l)
              !phase1 = CI_PHASE_N(b_str,b_temp2str,k,l)

              DO j=1,ORB_NUM
                IF(btest(b_str,j-1)) THEN
                  a_tempstr = ibclr(b_str,j-1)
                  DO i=1,ORB_NUM
                    IF(.NOT. btest(a_tempstr,i-1)) THEN
                      a_temp2str = ibset(a_tempstr,i-1)
                      a_adr1 =CALC_ADRESS_BINARY(ORB_NUM,NELECT_a,a_temp2str,Y_wa)
                      phase2 = phase1 *  CI_PHASE(b_str,i,j)
                      !phase2 = phase1 * CI_PHASE_N(b_str,a_temp2str,i,j)

                      tempsigma(q,q) = tempsigma(q,q) + &
                      phase2*INT_2(i,j,k,l)*CI_vector(a_adr1,b_adr1)

                    ENDIF
                  ENDDO
                ENDIF
              ENDDO
            ENDIF
          ENDDO
        ENDIF
      ENDDO
    ENDDO
 
    sigma = sigma + tempsigma
  
  ELSE 

    DO q=1,dimen_b

      b_str = b_str_arr(q)
      
      DO l=1,ORB_NUM
        IF(btest(b_str,l-1)) THEN
          b_tempstr = ibclr(b_str,l-1)
          DO k=1,ORB_NUM
            IF(.not. btest(b_tempstr,k-1)) THEN
              b_temp2str = ibset(b_tempstr,k-1)
              b_adr1 = CALC_ADRESS_BINARY(ORB_NUM,NELECT_b,b_temp2str,Y_wb)
              phase1 = CI_PHASE(b_str,k,l)
              !phase1 = CI_PHASE_N(b_str,b_temp2str,k,l)

              DO p=1,dimen_a

                a_str = a_str_arr(p)

                DO j=1,ORB_NUM
                  IF(btest(a_str,j-1)) THEN
                    a_tempstr = ibclr(a_str,j-1)
                    DO i=1,ORB_NUM
                      IF(.NOT. btest(a_tempstr,i-1)) THEN
                        a_temp2str = ibset(a_tempstr,i-1)
                        a_adr1 = CALC_ADRESS_BINARY(ORB_NUM,NELECT_a,a_temp2str,Y_wa)
                        phase2 = phase1 *  CI_PHASE(a_str,i,j)
                        !phase2 = phase1 * CI_PHASE_N(a_str,a_temp2str,i,j)

                        sigma(p,q) = sigma(p,q) + &
                        phase2*INT_2(i,j,k,l)*CI_vector(a_adr1,b_adr1)

                      ENDIF
                    ENDDO
                  ENDIF
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDIF
      ENDDO
    ENDDO

  ENDIF

END SUBROUTINE ALPHA_BETA_ROUTINE

SUBROUTINE ALPHA_BETA_ROUTINE_VEC(ORB_NUM,NELECT_a,NELECT_b,INT_2,dimen_a,dimen_b, &
                              Y_wa,Y_wb,a_str_arr,b_str_arr,CI_vector,sigma,SYM)

  INTEGER,INTENT(IN):: ORB_NUM,NELECT_a,NELECT_b,dimen_a,dimen_b
  INTEGER(KIND=intkind),DIMENSION(dimen_b),INTENT(IN)::a_str_arr,b_str_arr
  INTEGER,DIMENSION(2,ORB_NUM,NELECT_a):: Y_wa
  INTEGER,DIMENSION(2,ORB_NUM,NELECT_b):: Y_wb
  REAL(wp),DIMENSION(dimen_a,dimen_b),INTENT(IN):: CI_vector
  REAL(wp),DIMENSION(ORB_NUM,ORB_NUM,ORB_NUM,ORB_NUM),INTENT(IN):: INT_2
  REAL(wp),DIMENSION(dimen_a,dimen_b),INTENT(INOUT):: sigma
  INTEGER,INTENT(IN):: SYM 

  !local variables
  INTEGER:: p,q
  INTEGER:: i,j , s
  INTEGER:: k,l
  INTEGER:: b_adr
  INTEGER:: phase
  INTEGER,DIMENSION(2*dimen_a):: L_adr_list , R_adr_list
  REAL,DIMENSION(2*dimen_a) :: SIGNUM

  INTEGER(KIND=intkind):: a_str, b_str , a_tempstr , a_tempstr2 , b_tempstr ,  b_tempstr2
  INTEGER:: list_size , ind1
  REAL(wp),DIMENSION(2*dimen_a,dimen_b):: CI_vector2
  REAL(wp),DIMENSION(dimen_b):: F_array
  REAL(wp),DIMENSION(2*dimen_a):: V_array
  REAL(sp):: presign

  IF(SYM==0) THEN
    DO l=1,ORB_NUM
      DO k=l,ORB_NUM
         list_size = 0
         !Set up lists L_adr_list, R_adr_list , SIGNUM
         IF(k==l) THEN
           presign = 0.5
         ELSE
           presign = 1.0
         ENDIF

         DO p=1,dimen_a
           a_str= a_str_arr(p)

           !E_{kl} exciaation
           IF(btest(a_str,l-1)) THEN
             a_tempstr = ibclr(a_str,l-1)
             IF(.NOT. btest(a_tempstr,k-1)) THEN
               a_tempstr2 = ibset(a_tempstr,k-1)
               list_size = list_size + 1 
               R_adr_list(list_size) = p
               L_adr_list(list_size) = CALC_ADRESS_BINARY(ORB_NUM,NELECT_a,a_tempstr2,Y_wa)
               SIGNUM(list_size) = REAL(CI_PHASE(a_str,k,l)) * presign
             ENDIF 
           ENDIF

           !E_{lk} excitation
           IF(btest(a_str,k-1)) THEN
             a_tempstr = ibclr(a_str,k-1)
             IF(.NOT. btest(a_tempstr,l-1)) THEN
               a_tempstr2 = ibset(a_tempstr,l-1)
               list_size = list_size + 1 
               R_adr_list(list_size) = p
               L_adr_list(list_size) = CALC_ADRESS_BINARY(ORB_NUM,NELECT_a,a_tempstr2,Y_wa)
               SIGNUM(list_size) = REAL(CI_PHASE(a_str,l,k)) * presign
             ENDIF
           ENDIF

         ENDDO !p
    

         !Form C'
         DO s=1,list_size
          ind1 = L_adr_list(s)
          CI_vector2(s,:) = CI_vector(ind1,:) * SIGNUM(s)
         ENDDO

         !loop over beta strings
         DO q=1,dimen_b
           !Set up F array
           F_array = 0.0_wp
           b_str = b_str_arr(q)
           DO j=1,ORB_NUM
             IF(btest(b_str,j-1)) THEN
               b_tempstr = ibclr(b_str,j-1)
               DO i=1,ORB_NUM
                 IF(.NOT. btest(b_tempstr,i-1)) THEN
                   b_tempstr2 = ibset(b_tempstr,i-1)
                   b_adr = CALC_ADRESS_BINARY(ORB_NUM,NELECT_b,b_tempstr2,Y_wb)
                   phase = CI_PHASE(b_str,i,j)
                   F_array(b_adr) = F_array(b_adr) + phase*INT_2(i,j,k,l)
                 ENDIF
               ENDDO !i         
             ENDIF
           ENDDO !j
           CALL MULT_CF(list_size,dimen_b,CI_vector2(1:list_size,:),F_array,V_array(1:list_size))
           DO s=1,list_size
             ind1 = R_adr_list(s)
             sigma(ind1,q) = sigma(ind1,q) + V_array(s)
           ENDDO
         ENDDO !q
      ENDDO !k
    ENDDO !l
  ELSE !SYM=1

    DO l=1,ORB_NUM
      DO k=1,ORB_NUM
         list_size = 0
         
         !Set up lists L_adr_list, R_adr_list , SIGNUM
         DO p=1,dimen_a
           a_str= a_str_arr(p)

           !E_{kl} Excitation
           IF(btest(a_str,l-1)) THEN
             a_tempstr = ibclr(a_str,l-1)
             IF(.NOT. btest(a_tempstr,k-1)) THEN
               a_tempstr2 = ibset(a_tempstr,k-1)
               list_size = list_size + 1 
               R_adr_list(list_size) = p
               L_adr_list(list_size) = CALC_ADRESS_BINARY(ORB_NUM,NELECT_a,a_tempstr2,Y_wa)
               SIGNUM(list_size) = REAL(CI_PHASE(a_str,k,l))
             ENDIF 
           ENDIF

         ENDDO !p

         !Form C'
         DO s=1,list_size
          ind1 = L_adr_list(s)
          CI_vector2(s,:) = CI_vector(ind1,:) * SIGNUM(s)
         ENDDO

         !loop over beta strings
         !Set up F array
         DO q=1,dimen_b
           F_array = 0.0_wp
           b_str = b_str_arr(q)
           DO j=l+1,ORB_NUM  !(ij)>=(kl)
             IF(btest(b_str,j-1)) THEN
               b_tempstr = ibclr(b_str,j-1)
               DO i=k+1,ORB_NUM !(ij)>=(kl)
                 IF(.NOT. btest(b_tempstr,i-1)) THEN
                   b_tempstr2 = ibset(b_tempstr,i-1)
                   b_adr = CALC_ADRESS_BINARY(ORB_NUM,NELECT_b,b_tempstr2,Y_wb)
                   phase = CI_PHASE(b_str,i,j)
                   F_array(b_adr) = F_array(b_adr) + phase*INT_2(i,j,k,l)
                 ENDIF
               ENDDO !i         
             ENDIF
           ENDDO !j

           !Add i>=k , j<l
           !DO j=1,l-1
           DO j=1,l-1
             IF(btest(b_str,j-1)) THEN
               b_tempstr = ibclr(b_str,j-1)
               DO i=k+1,ORB_NUM !(ij)>=(kl)
                 IF(.NOT. btest(b_tempstr,i-1)) THEN
                   b_tempstr2 = ibset(b_tempstr,i-1)
                   b_adr = CALC_ADRESS_BINARY(ORB_NUM,NELECT_b,b_tempstr2,Y_wb)
                   phase = CI_PHASE(b_str,i,j)
                   F_array(b_adr) = F_array(b_adr) + phase*INT_2(i,j,k,l)
                 ENDIF
               ENDDO !i
             ENDIF
           ENDDO !j

           !Add i=k , j>l
           DO j=l+1,ORB_NUM
             IF(btest(b_str,j-1)) THEN
               b_tempstr = ibclr(b_str,j-1)
               IF(.NOT. btest(b_tempstr,k-1)) THEN
                 b_tempstr2 = ibset(b_tempstr,k-1)
                 b_adr = CALC_ADRESS_BINARY(ORB_NUM,NELECT_b,b_tempstr2,Y_wb)
                 phase = CI_PHASE(b_str,k,j)
                 F_array(b_adr) = F_array(b_adr) + phase*INT_2(k,j,k,l)
               ENDIF
             ENDIF
           ENDDO !j

           !Add j=l , i>k
           IF(btest(b_str,l-1)) THEN
             b_tempstr = ibclr(b_str,l-1)
             DO i=k+1,ORB_NUM
               IF(.NOT. btest(b_tempstr,i-1)) THEN
                 b_tempstr2 = ibset(b_tempstr,i-1)
                 b_adr = CALC_ADRESS_BINARY(ORB_NUM,NELECT_b,b_tempstr2,Y_wb)
                 phase = CI_PHASE(b_str,i,l)
                 F_array(b_adr) = F_array(b_adr) + phase*INT_2(i,l,k,l)
               ENDIF
             ENDDO
           ENDIF
             
           !Add i=k , j=l
           IF(btest(b_str,l-1)) THEN
             b_tempstr = ibclr(b_str,l-1)
             IF(.NOT. btest(b_tempstr,k-1)) THEN
               b_tempstr2 = ibset(b_tempstr,k-1)
               b_adr = CALC_ADRESS_BINARY(ORB_NUM,NELECT_b,b_tempstr2,Y_wb)
               phase = CI_PHASE(b_str,k,l) 
               F_array(b_adr) = F_array(b_adr) + 0.5_wp*phase*INT_2(k,l,k,l)
             ENDIF
           ENDIF          

           CALL MULT_CF(list_size,dimen_b,CI_vector2(1:list_size,:),F_array,V_array(1:list_size))

           DO s=1,list_size
             ind1 = R_adr_list(s)
             sigma(ind1,q) = sigma(ind1,q) + V_array(s)
           ENDDO

          ENDDO !q
      ENDDO !k
    ENDDO !l
  ENDIF !SYM

  CONTAINS

  SUBROUTINE MULT_CF(N,M,A,x,y)
    INTEGER,INTENT(IN):: N,M
    REAL(wp),DIMENSION(N,M),INTENT(IN):: A
    REAL(wp),DIMENSION(M),INTENT(IN):: x
    REAL(wp),DIMENSION(N),INTENT(OUT):: y
   
    !local variables
    REAL(wp):: alpha , beta
    alpha = 1.0_wp
    beta = 0.0_wp

    IF(wp==sp) THEN
      CALL SGEMV("N",N,M,alpha,A,N,x,1,beta,y,1)
    ELSEIF(wp==dp) THEN
      CALL DGEMV("N",N,M,alpha,A,N,x,1,beta,y,1)
    ELSE
      IF(comm_world%mpirank==0) THEN
        WRITE(*,*) "There is no specific LAPACK routine for this precision: " , wp
      ENDIF
    ENDIF
  END SUBROUTINE MULT_CF

END SUBROUTINE ALPHA_BETA_ROUTINE_VEC


SUBROUTINE SIGMA_STEP(ORB_NUM,NELECT_a,NELECT_b,HAM_1,INT_2,dimen_a,dimen_b, &
                      Y_wa,Y_wb,a_str_arr,b_str_arr,LCI_VEC,CI_vector,sigma,SYM,SPIN)
  INTEGER,INTENT(IN):: ORB_NUM,NELECT_a,NELECT_b,dimen_a,dimen_b
  LOGICAL,INTENT(IN):: LCI_VEC
  INTEGER(KIND=intkind),DIMENSION(dimen_a),INTENT(IN)::a_str_arr
  INTEGER(KIND=intkind),DIMENSION(dimen_b),INTENT(IN)::b_str_arr
  INTEGER,DIMENSION(2,ORB_NUM,NELECT_a):: Y_wa
  INTEGER,DIMENSION(2,ORB_NUM,NELECT_b):: Y_wb
  REAL(wp),DIMENSION(dimen_a,dimen_b),INTENT(IN):: CI_vector
  REAL(wp),DIMENSION(ORB_NUM,ORB_NUM,INFO%ISPIN1),INTENT(IN):: HAM_1
  REAL(wp),DIMENSION(ORB_NUM,ORB_NUM,ORB_NUM,ORB_NUM,INFO%ISPIN2),INTENT(IN):: INT_2
  REAL(wp),DIMENSION(dimen_a,dimen_b),INTENT(OUT):: sigma
  INTEGER,INTENT(INOUT):: SYM , SPIN

  !local variables
  INTEGER:: inda , indb , indab

  sigma = 0.0_wp

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

  IF(LCI_VEC) THEN !Vectorized algorithm
    CALL BETA_BETA_ROUTINE_VEC(ORB_NUM,NELECT_b,HAM_1(:,:,indb),INT_2(:,:,:,:,indb), &
                               dimen_a,dimen_b,Y_wb,b_str_arr,CI_vector,sigma)
    IF(SYM==1) THEN
      CALL ALPHA_BETA_ROUTINE_VEC(ORB_NUM,NELECT_a,NELECT_b,INT_2(:,:,:,:,indab),dimen_a,dimen_b, &
                                  Y_wa,Y_wb,a_str_arr,b_str_arr,CI_vector,sigma,SYM)
      sigma = sigma + (-1)**SPIN*TRANSPOSE(sigma)
    ELSE
      CALL ALPHA_ALPHA_ROUTINE_VEC(ORB_NUM,NELECT_a,HAM_1(:,:,inda),INT_2(:,:,:,:,inda), &
                                   dimen_a,dimen_b,Y_wa,a_str_arr,CI_vector,sigma)
      CALL ALPHA_BETA_ROUTINE_VEC(ORB_NUM,NELECT_a,NELECT_b,INT_2(:,:,:,:,indab),dimen_a,dimen_b, &
                                   Y_wa,Y_wb,a_str_arr,b_str_arr,CI_vector,sigma,SYM) 
    ENDIF
  ELSE  !Element-Wise algorithm
    CALL BETA_BETA_ROUTINE(ORB_NUM,NELECT_b,HAM_1(:,:,indb),INT_2(:,:,:,:,indb), &
                           dimen_a,dimen_b,Y_wb,b_str_arr,CI_vector,sigma)
    IF (SYM==1) THEN
      sigma = sigma + (-1)**SPIN*TRANSPOSE(sigma)
    ELSE
      CALL ALPHA_ALPHA_ROUTINE(ORB_NUM,NELECT_a,HAM_1(:,:,inda),INT_2(:,:,:,:,inda), &
                               dimen_a,dimen_b,Y_wa,a_str_arr,CI_vector,sigma)
    ENDIF
      CALL ALPHA_BETA_ROUTINE(ORB_NUM,NELECT_a,NELECT_b,INT_2(:,:,:,:,indab),dimen_a,dimen_b, &
                              Y_wa,Y_wb,a_str_arr,b_str_arr,CI_vector,sigma,SYM,SPIN)
  ENDIF

END SUBROUTINE SIGMA_STEP


SUBROUTINE INITIALIZE_CI_VECTOR(dimen_a,dimen_b,CI_H_DIAG,CI_vector,NUM_EIGS)
  INTEGER,INTENT(IN):: dimen_a , dimen_b 
  INTEGER,INTENT(INOUT):: NUM_EIGS
  REAL(wp),DIMENSION(dimen_a,dimen_b),INTENT(IN):: CI_H_DIAG
  TYPE(EXP_SPACE),DIMENSION(:),ALLOCATABLE,INTENT(OUT):: CI_vector

  !local variables
  INTEGER:: adr_a , adr_b
  INTEGER:: maxdim
  REAL(wp),DIMENSION(NUM_EIGS):: minvals
  INTEGER(KIND=8),DIMENSION(NUM_EIGS):: ind
  REAL(wp),DIMENSION(:),ALLOCATABLE:: tempvec 
  REAL(wp):: letter , former
  
  INTEGER(KIND=8):: i
  INTEGER:: tempmaxind
  REAL(wp):: maxv , tempval

  ALLOCATE(CI_vector(NUM_EIGS))
  DO i=1,NUM_EIGS
    ALLOCATE(CI_vector(i)%exp_vec(dimen_a,dimen_b))
  ENDDO

  IF (NUM_EIGS == 1) THEN
    CI_vector(1)%exp_vec = 0.0_wp
    CI_vector(1)%exp_vec(1,1) = 1.0_wp
  ELSE
    maxdim = dimen_a*dimen_b
    ALLOCATE(tempvec(maxdim))
    tempvec = RESHAPE(CI_H_DIAG,(/maxdim/))
    
    IF(NUM_EIGS>maxdim) THEN
      IF(comm_world%mpirank==0) THEN
        WRITE(*,*) "Number of desired eigenpairs exceed dimension of CI Matrix"
        WRITE(*,*) "NUM_EIGS set to max. possible value"
      ENDIF
      NUM_EIGS = maxdim
    ENDIF

    !Fill the minvals with NUM_EIGS first values from tempvec
    DO i=1,NUM_EIGS
      minvals(i) = tempvec(i)
      ind(i) = i
    ENDDO

    CALL FIND_MAX(maxv,tempmaxind)

    !Find the minimal NUM_EIGS elements from CI_H_DIAG and store them in vector
    !and their indices in integer array ind - use tempvec(dimen_a*dimen_b)
    DO i=NUM_EIGS+1,maxdim
      tempval = tempvec(i)
      IF(tempval<maxv) THEN
        minvals(tempmaxind) = tempval
        ind(tempmaxind) = i
        CALL FIND_MAX(maxv,tempmaxind)
      ENDIF
    ENDDO
   
    !Sort the list of minimal NUM_EIGS elements
    CALL SORT_LIST()

    !Store NUM_EIGS guess CI vectors approximating CI Matrix with
    !its diagonal matrix
    IF(dimen_a==dimen_b) THEN
      IF(comm_world%mpirank==0) THEN
        WRITE(*,*) "minvals",  minvals
      ENDIF
      letter = 0.0_wp
      former = 0.0_wp

      DO i=1,NUM_EIGS
        adr_b = 1 + ind(i)/dimen_a
        adr_a = ind(i) - dimen_a*(adr_b-1)
 
        former = letter
        letter = minvals(i)

        CI_vector(i)%exp_vec = 0.0_wp
      
        IF(comm_world%mpirank==0) THEN
          WRITE(*,*) "Guess", adr_a , adr_b , letter
        ENDIF

        IF(adr_a==adr_b) THEN
          CI_vector(i)%exp_vec(adr_a,adr_b) = 1.0_wp
        ELSE
          IF(ABS(former-letter)<1.0e-6_wp) THEN
            CI_vector(i)%exp_vec(adr_a,adr_b) = 1.0_wp/SQRT(2.0_wp)
            CI_vector(i)%exp_vec(adr_b,adr_a) = 1.0_wp/SQRT(2.0_wp)
          ELSE
            CI_vector(i)%exp_vec(adr_a,adr_b) = -1.0_wp/SQRT(2.0_wp)
            CI_vector(i)%exp_vec(adr_b,adr_a) = 1.0_wp/SQRT(2.0_wp)
          ENDIF
        ENDIF

      ENDDO
    ELSE
      DO i=1,NUM_EIGS
        adr_b = 1 + ind(i)/dimen_a
        adr_a = ind(i) - dimen_a*(adr_b-1)
        
        CI_vector(i)%exp_vec = 0.0_wp
        IF(comm_world%mpirank==0) THEN
          WRITE(*,*) "Guess", adr_a , adr_b , letter
        ENDIF
        CI_vector(i)%exp_vec(adr_a,adr_b) = 1.0_wp
      ENDDO
    ENDIF
  ENDIF

  CONTAINS

  SUBROUTINE FIND_MAX(maxv,tempmaxind)
    REAL(wp),INTENT(OUT):: maxv
    INTEGER,INTENT(OUT):: tempmaxind
 
    !local variables
    INTEGER:: i 
 
    maxv = minvals(1)
    tempmaxind = 1
    DO i=2,NUM_EIGS
      IF(minvals(i)>maxv) THEN
        maxv = minvals(i)
        tempmaxind = i
      ENDIF
    ENDDO
  END SUBROUTINE

  SUBROUTINE SORT_LIST()
    INTEGER:: i , j
    INTEGER(KIND=8):: tempind
    REAL(wp):: tempval

    DO i =1,NUM_EIGS-1
      DO j=i+1,NUM_EIGS
       IF(minvals(j)<minvals(i)) THEN
         tempval = minvals(i)
         minvals(i) = minvals(j)
         minvals(j) = tempval

         tempind = ind(i)
         ind(i) = ind(j)
         ind(j) = tempind
       ENDIF  
      ENDDO
    ENDDO 
  END SUBROUTINE

END SUBROUTINE INITIALIZE_CI_VECTOR

SUBROUTINE CI_STEEPEST_DESCENT_LOOP(ORB_NUM,NELECT_a,NELECT_b,dimen_a,dimen_b,HAM_1,INT_2, &
                   Y_wa,Y_wb,a_str_arr,b_str_arr,CI_H_DIAG,LCI_VEC,CI_vector,CI_EIG,counter, &
                   MAXITER,EPS,SYM,SPIN,alpha)
  INTEGER,INTENT(IN):: ORB_NUM , NELECT_a , NELECT_b
  INTEGER,INTENT(IN):: dimen_a , dimen_b
  LOGICAL,INTENT(IN):: LCI_VEC
  REAL(wp),DIMENSION(ORB_NUM,ORB_NUM,INFO%ISPIN1),INTENT(IN):: HAM_1
  REAL(wp),DIMENSION(ORB_NUM,ORB_NUM,ORB_NUM,ORB_NUM,INFO%ISPIN2),INTENT(IN):: INT_2
  INTEGER,DIMENSION(2,ORB_NUM,NELECT_a):: Y_wa
  INTEGER,DIMENSION(2,ORB_NUM,NELECT_b):: Y_wb
  INTEGER(KIND=intkind),DIMENSION(dimen_a),INTENT(IN)::a_str_arr
  INTEGER(KIND=intkind),DIMENSION(dimen_b),INTENT(IN)::b_str_arr
  REAL(wp),DIMENSION(dimen_a,dimen_b),INTENT(IN):: CI_H_DIAG
  REAL(wp),DIMENSION(dimen_a,dimen_b),INTENT(INOUT):: CI_vector 
  REAL(wp),INTENT(OUT):: CI_EIG
  INTEGER,INTENT(INOUT):: SYM , SPIN
  INTEGER,INTENT(IN):: MAXITER
  REAL(wp),INTENT(IN):: EPS 
  REAL(wp),INTENT(IN):: alpha
  INTEGER,INTENT(OUT):: counter

  !local variables
  REAL(wp):: CI_NORM , CI_NORM_OLD
  REAL(wp):: CI_EIG_OLD
  REAL(wp),DIMENSION(dimen_a,dimen_b):: sigma



  CI_NORM = 1000.0_wp
  CI_EIG = 0.0_wp
  counter = 0

  IF(comm_world%mpirank==0) THEN
    !open files where you want to write results
    OPEN(UNIT=out,FILE="qmcfort_out",STATUS="old",POSITION="append",ACTION="write")

    WRITE(out,*) starshort
    WRITE(out,*) "ENTERING CI IERATIVE DIAGONALIZATION PROCEDURE:"
    WRITE(out,109) "CI Step" , "CI Norm " , "CI EV"
    109 FORMAT(1X,T5,A,T15,A,T35,A) 
  ENDIF

  DO WHILE (CI_NORM>EPS .and. counter<MAXITER)
    CI_NORM_OLD = CI_NORM
    CI_EIG_OLD = CI_EIG
    counter = counter + 1
    sigma=0.0_wp

    CALL SIGMA_STEP(ORB_NUM,NELECT_a,NELECT_b,HAM_1,INT_2,dimen_a,dimen_b, & 
                    Y_wa,Y_wb,a_str_arr,b_str_arr,LCI_VEC,CI_vector,sigma,SYM,SPIN)

    CALL STEEPEST_DESCENT_STEP(dimen_a,dimen_b,CI_vector,sigma,CI_H_DIAG,CI_EIG,CI_NORM,alpha)

    IF(comm_world%mpirank==0) THEN
      WRITE(out,110) counter, CI_NORM , CI_EIG
      110 FORMAT (1X,T5,I3,T15,F10.6,T35,F10.6)
    ENDIF

    IF(CI_NORM>CI_NORM_OLD) THEN
      IF(comm_world%mpirank==0) THEN
        WRITE(out,*) "  CI diagonalization procedure failed - do not converge"
      ENDIF
      CI_NORM = CI_NORM_OLD
      CI_EIG = CI_EIG_OLD
      EXIT
    ENDIF
  ENDDO

  IF(comm_world%mpirank==0) THEN
    WRITE(out,120) counter
    120 FORMAT (1X,"CI Diagonalization procedure converged after ",I3, " steps")
    WRITE(out,130) CI_EIG
    130 FORMAT (1X,"E_CI",T20,"= ",F10.6)
    WRITE(out,*) starshort

    !close open files
    CLOSE(UNIT=out)
  ENDIF

END SUBROUTINE CI_STEEPEST_DESCENT_LOOP

SUBROUTINE CI_DAVIDSON_LOOP(ORB_NUM,NELECT_a,NELECT_b,dimen_a,dimen_b,HAM_1,INT_2, &
                   Y_wa,Y_wb,a_str_arr,b_str_arr,CI_H_DIAG,LCI_VEC,CI_vector,CI_EIG,  &
                   counter,MAXITER,EPS,SYM,SPIN)
  INTEGER,INTENT(IN):: ORB_NUM , NELECT_a , NELECT_b
  INTEGER,INTENT(IN):: dimen_a , dimen_b
  LOGICAL,INTENT(IN):: LCI_VEC
  REAL(wp),DIMENSION(ORB_NUM,ORB_NUM,INFO%ISPIN1),INTENT(IN):: HAM_1
  REAL(wp),DIMENSION(ORB_NUM,ORB_NUM,ORB_NUM,ORB_NUM,INFO%ISPIN2),INTENT(IN):: INT_2
  INTEGER,DIMENSION(2,ORB_NUM,NELECT_a):: Y_wa
  INTEGER,DIMENSION(2,ORB_NUM,NELECT_b):: Y_wb
  INTEGER(KIND=intkind),DIMENSION(dimen_a),INTENT(IN)::a_str_arr
  INTEGER(KIND=intkind),DIMENSION(dimen_b),INTENT(IN)::b_str_arr
  REAL(wp),DIMENSION(dimen_a,dimen_b),INTENT(IN):: CI_H_DIAG
  REAL(wp),DIMENSION(dimen_a,dimen_b),INTENT(INOUT):: CI_vector
  REAL(wp),DIMENSION(1),INTENT(OUT):: CI_EIG
  INTEGER,INTENT(INOUT):: SYM , SPIN
  INTEGER,INTENT(IN):: MAXITER
  REAL(wp),INTENT(IN):: EPS
  INTEGER,INTENT(OUT):: counter

  !local variables
  REAL(wp):: CI_NORM_OLD , CI_NORM , RITZ_corr
  REAL(wp),DIMENSION(dimen_a,dimen_b):: sigma
  INTEGER:: j , k 

  TYPE(EXP_SPACE),DIMENSION(MAXITER):: DAV_SPACE
  TYPE(EXP_SPACE),DIMENSION(MAXITER):: SIGMA_SPACE
  REAL(wp),DIMENSION(MAXITER,MAXITER):: B , B_aux
  REAL(wp),DIMENSION(MAXITER,1):: B_vec 
  REAL(wp),DIMENSION(MAXITER):: s_vec

  B = 0.0_wp
  B_vec = 0.0_wp
  s_vec = 0.0_wp

  IF(comm_world%mpirank==0) THEN
    !open files where you want to write results
    OPEN(UNIT=out,FILE="qmcfort_out",STATUS="old",POSITION="append",ACTION="write")
    WRITE(out,*) "CI ITERATIVE DIAGONALIZATION PROCEDURE:"
    WRITE(out,159) "CI Step" , "CI Norm" , "CI EV"
    159 FORMAT (1X,T5,A,T15,A,T35,A)
  ENDIF

  DO j=1,MAXITER
    counter = j 
   
    !Allocate j-th expansion vector
    ALLOCATE(DAV_SPACE(j)%exp_vec(dimen_a,dimen_b))
    ALLOCATE(SIGMA_SPACE(j)%exp_vec(dimen_a,dimen_b))
     
    IF (j==1) THEN
      DAV_SPACE(j)%exp_vec = CI_vector
    ELSE
      DAV_SPACE(j)%exp_vec = sigma
    ENDIF

    CALL SIGMA_STEP(ORB_NUM,NELECT_a,NELECT_b,HAM_1,INT_2,dimen_a,dimen_b, &
                   Y_wa,Y_wb,a_str_arr,b_str_arr,LCI_VEC,DAV_SPACE(j)%exp_vec,SIGMA_SPACE(j)%exp_vec, &
                   SYM,SPIN)

    sigma = SIGMA_SPACE(j)%exp_vec

    DO k=1,j-1
      B(k,j) = SUM(DAV_SPACE(k)%exp_vec*sigma)
    ENDDO    
 
    B(j,j) = SUM(DAV_SPACE(j)%exp_vec*sigma)
   
    B_aux = B 
    !diagonalize B and find the lowest eigenvalue and corresponding eigenvector
    CI_NORM_OLD = CI_NORM
    CALL DIAG_MATRIX(B_aux(1:j,1:j),B_vec(1:j,1),CI_EIG,j,1)
    !residuum vector
    CI_vector = SPACE_DOT_PROD(DAV_SPACE,B_vec,MAXITER,j,dimen_a,dimen_b)
    sigma = SPACE_DOT_PROD(SIGMA_SPACE,B_vec,MAXITER,j,dimen_a,dimen_b)  - &
            CI_EIG(1)*CI_vector
    CI_NORM = norm2(sigma)


    IF(comm_world%mpirank==0) THEN
      WRITE(out,160) counter, CI_NORM , CI_EIG
      160 FORMAT (1X,T5,I2,T15,F10.6,T35,F10.6)
    ENDIF

    RITZ_corr = SUM(CI_vector*sigma/(CI_H_DIAG-CI_EIG(1)-2.0_wp*CI_vector*sigma+CI_vector**2)) / &
                SUM(CI_vector**2/(CI_H_DIAG-CI_EIG(1)-2.0_wp*CI_vector*sigma+CI_vector**2))

    !Preconditioning
    sigma = (-sigma + RITZ_corr*CI_vector)/(CI_H_DIAG-CI_EIG(1)-2.0_wp*CI_vector*sigma+CI_vector**2)

    sigma = sigma / CI_NORM

    !Orthogonalization
    s_vec = SPACE_PROD(DAV_SPACE,sigma,MAXITER,j,dimen_a,dimen_b)
    sigma = sigma - SPACE_DOT_PROD(DAV_SPACE,s_vec,MAXITER,j,dimen_a,dimen_b)
    
    !Normalisation
    sigma = sigma/NORM2(sigma) 

    IF(CI_NORM<EPS) THEN
      EXIT
    ENDIF

  ENDDO

  IF(comm_world%mpirank==0) THEN
    IF (counter==MAXITER) THEN
      WRITE(out,*) "CI diagonalization procedure does not converge"
      WRITE(out,161) MAXITER
      161 FORMAT (1X,"Maximal number of iterations = ",I3," is reached")
      WRITE(out,162) CI_EIG(1)
      162 FORMAT (1X,"Current eigenvalue = ",F10.6)
    ELSE
      WRITE(out,164) counter
      164 FORMAT (1X,"CI diagonalization procedure converged after ",I3," steps")
      WRITE(out,163) CI_EIG(1)
      163 FORMAT (1X,"E_CI",T20,"= ",F10.6)
    ENDIF
    
    WRITE(out,*) starshort
    !close open files
    CLOSE(UNIT=out)
  ENDIF
END SUBROUTINE CI_DAVIDSON_LOOP

SUBROUTINE CI_DAVIDSON_LIU_LOOP(ORB_NUM,NELECT_a,NELECT_b,dimen_a,dimen_b,HAM_1,INT_2, &
                   Y_wa,Y_wb,a_str_arr,b_str_arr,CI_H_DIAG,NUM_EIGS,LCI_VEC,CI_vector,CI_EIG,  &
                   counter,MAXITER,EPS,SYM,SPIN)
  INTEGER,INTENT(IN):: ORB_NUM , NELECT_a , NELECT_b , NUM_EIGS
  INTEGER,INTENT(IN):: dimen_a , dimen_b
  LOGICAL,INTENT(IN):: LCI_VEC
  REAL(wp),DIMENSION(ORB_NUM,ORB_NUM,INFO%ISPIN1),INTENT(IN):: HAM_1
  REAL(wp),DIMENSION(ORB_NUM,ORB_NUM,ORB_NUM,ORB_NUM,INFO%ISPIN2),INTENT(IN):: INT_2
  INTEGER,DIMENSION(2,ORB_NUM,NELECT_a):: Y_wa
  INTEGER,DIMENSION(2,ORB_NUM,NELECT_b):: Y_wb
  INTEGER(KIND=intkind),DIMENSION(dimen_a),INTENT(IN)::a_str_arr
  INTEGER(KIND=intkind),DIMENSION(dimen_b),INTENT(IN)::b_str_arr
  REAL(wp),DIMENSION(dimen_a,dimen_b),INTENT(IN):: CI_H_DIAG
  TYPE(EXP_SPACE),DIMENSION(NUM_EIGS),INTENT(INOUT):: CI_vector
  REAL(wp),DIMENSION(NUM_EIGS),INTENT(OUT):: CI_EIG
  INTEGER,INTENT(INOUT):: SYM , SPIN
  INTEGER,INTENT(IN):: MAXITER
  REAL(wp),INTENT(IN):: EPS 
  INTEGER,INTENT(OUT):: counter
  
  !local variables
  REAL(wp):: CI_NORM , RITZ_corr
  REAL(wp),DIMENSION(dimen_a,dimen_b):: sigma
  INTEGER:: i, j , k 
  INTEGER:: vec_conv ,vec_contr
  INTEGER:: nb , nb_old , dnb

  TYPE(EXP_SPACE),DIMENSION(NUM_EIGS*MAXITER):: DAV_SPACE
  TYPE(EXP_SPACE),DIMENSION(NUM_EIGS*MAXITER):: SIGMA_SPACE
  REAL(wp),DIMENSION(NUM_EIGS*MAXITER,NUM_EIGS*MAXITER):: B , B_aux
  REAL(wp),DIMENSION(NUM_EIGS*MAXITER,NUM_EIGS):: B_vec
  REAL(wp),DIMENSION(NUM_EIGS):: B_EIG
  REAL(wp),DIMENSION(NUM_EIGS*MAXITER):: s_vec

  B = 0.0_wp
  B_vec = 0.0_wp
  s_vec = 0.0_wp
  B_EIG = 0.0_wp
  CI_EIG = 0.0_wp

  IF(comm_world%mpirank==0) THEN
    !open files where you want to write results
    OPEN(UNIT=out,FILE="qmcfort_out",STATUS="old",POSITION="append",ACTION="write")
    WRITE(out,*) "Entering CI diagonalization procedure"
    WRITE(out,*) starshort
    WRITE(out,170) "CI Step","CI Norm","CI EV"
    170 FORMAT(1X,T5,A,T15,A,T35,A) 
  ENDIF

  !Initialize first NUM_EIGS vectors of subspace
  nb_old = 0
  nb = NUM_EIGS
  dnb = 0
  vec_conv = 0
  vec_contr = 0

  !Start a loop now!  
  DO i=1,MAXITER-1
    counter = i

    !Calculate new W vectors w=Hv and new subspace matrix elements
    DO j=nb_old+1,nb

      IF(i==1) THEN
          ALLOCATE(DAV_SPACE(j)%exp_vec(dimen_a,dimen_b))
          ALLOCATE(SIGMA_SPACE(j)%exp_vec(dimen_a,dimen_b))
      
          DAV_SPACE(j)%exp_vec = CI_vector(j)%exp_vec
      ENDIF

      CALL SIGMA_STEP(ORB_NUM,NELECT_a,NELECT_b,HAM_1,INT_2,dimen_a,dimen_b, &
                     Y_wa,Y_wb,a_str_arr,b_str_arr,LCI_VEC,DAV_SPACE(j)%exp_vec,SIGMA_SPACE(j)%exp_vec, &
                     SYM,SPIN)

      DO k=1,j-1
        B(k,j) = SUM(DAV_SPACE(k)%exp_vec*SIGMA_SPACE(j)%exp_vec)
      ENDDO    
 
      B(j,j) = SUM(DAV_SPACE(j)%exp_vec*SIGMA_SPACE(j)%exp_vec)
    ENDDO

    B_aux = B
    CALL DIAG_MATRIX(B_aux(1:nb,1:nb),B_vec(1:nb,1:NUM_EIGS),B_EIG,nb,NUM_EIGS)

    DO j=vec_conv+1,NUM_EIGS
      CI_vector(j)%exp_vec = SPACE_DOT_PROD(DAV_SPACE,B_vec(:,j),MAXITER*NUM_EIGS, &
                                            nb,dimen_a,dimen_b)
      sigma = SIGMA_SPACE(j)%exp_vec

      !calculate residue vector
      sigma = SPACE_DOT_PROD(SIGMA_SPACE,B_vec(:,j),MAXITER*NUM_EIGS,nb,dimen_a,dimen_b) &
              - B_EIG(j)*CI_vector(j)%exp_vec
      CI_NORM = NORM2(sigma)

      IF(comm_world%mpirank==0) THEN
        IF(j==1) THEN
          WRITE(out,180)i, j, CI_NORM, B_EIG(j)
          180 FORMAT (1X,T5,I2,T15,I2,'. ',F10.6,T35,F10.6)
        ELSE
          WRITE(out,181) j, CI_NORM, B_EIG(j)
          181 FORMAT (1X,T15,I2,'. ',F10.6,T35,F10.6)
        ENDIF
      ENDIF

      IF(CI_NORM<=EPS) THEN
        vec_contr = vec_contr + 1
        CI_EIG(vec_conv+vec_contr) = B_EIG(vec_conv+vec_contr)
        CYCLE
      ENDIF

     !Rayleigh-Ritz correction
      RITZ_corr = SUM(CI_vector(j)%exp_vec*sigma/(CI_H_DIAG-B_EIG(j) - &
                      2.0_wp*CI_vector(j)%exp_vec*sigma+CI_vector(j)%exp_vec**2)) /&
                  SUM(CI_vector(j)%exp_vec**2/(CI_H_DIAG-B_EIG(j) - &
                      2.0_wp*CI_vector(j)%exp_vec*sigma+CI_vector(j)%exp_vec**2))

      !Preconditioner
      sigma = (-sigma + RITZ_corr*CI_vector(j)%exp_vec)/(CI_H_DIAG-B_EIG(j) - &
                   2.0_wp*CI_vector(j)%exp_vec*sigma+CI_vector(j)%exp_vec**2)

      !sigma = -sigma/(CI_H_DIAG-B_EIG(j))

      sigma = sigma/NORM2(sigma)

      !GS orthogonalization with DAV Space

      s_vec = SPACE_PROD(DAV_SPACE,sigma,MAXITER*NUM_EIGS,nb+dnb,dimen_a,dimen_b)
      sigma = sigma - SPACE_DOT_PROD(DAV_SPACE,s_vec,MAXITER*NUM_EIGS,nb+dnb,dimen_a,dimen_b)

      !GS orthogonalization with converged eigenvectors

      s_vec = SPACE_PROD(CI_vector,sigma,MAXITER*NUM_EIGS,vec_conv,dimen_a,dimen_b)
      sigma = sigma - SPACE_DOT_PROD(CI_vector,s_vec,MAXITER*NUM_EIGS,vec_conv,dimen_a,dimen_b)

      CI_NORM = NORM2(sigma)

      IF (CI_NORM > 0.001_wp) THEN
        dnb = dnb + 1
        IF(.NOT. ALLOCATED(DAV_SPACE(nb+dnb)%exp_vec)) THEN
          ALLOCATE(DAV_SPACE(nb+dnb)%exp_vec(dimen_a,dimen_b))
          ALLOCATE(SIGMA_SPACE(nb+dnb)%exp_vec(dimen_a,dimen_b))
        ENDIF

        DAV_SPACE(nb+dnb)%exp_vec = sigma/CI_NORM
      ENDIF
    ENDDO

    
    vec_conv = vec_conv + vec_contr
   
   ! Deflate and restart DAV space
   ! IF(dnb==0 .or. vec_contr/=0) THEN
   !   !Restart Dav Space
   !   DAV_SPACE(1:NUM_EIGS-vec_conv) = CI_vector(1:vec_conv) !CI_vector(vec_conv+1:NUM_EIGS)
   !   DAV_SPACE(NUM_EIGS-vec_conv+1:NUM_EIGS-vec_conv+dnb) = DAV_SPACE(nb_old+1:nb_old+dnb)

   !   tempnb = nb
   !   nb_old = 0
   !   nb = dnb + NUM_EIGS - vec_conv
   !   dnb = 0
   !   WRITE(out,*) "Deflation from dim" , tempnb , "to" , nb 
   ! ELSE
      !expand Dav Space
    nb_old = nb
    nb = nb + dnb
    dnb = 0
    IF(comm_world%mpirank==0) THEN
      WRITE(*,*)"it", i , "exp space dim:" , nb
    ENDIF
   ! ENDIF
    vec_contr = 0

    IF(vec_conv==NUM_EIGS) THEN
      EXIT
    ENDIF

  ENDDO

  IF(comm_world%mpirank==0) THEN
    WRITE(out,*) starshort

    IF (counter==MAXITER-1) THEN
      WRITE(out,261) MAXITER
      261 FORMAT (1X,"Maximal number of iterations = ",I3," is reached")
      WRITE(out,*) "Not all eigenvalues converge"
      WRITE(out,262) vec_conv
      262 FORMAT ("Number of converged eigenvalues= ",I3)
      DO i=1,NUM_EIGS
        WRITE(out,263) i, CI_EIG(i)
        263 FORMAT (1X,"Current eigenvalue No. ",I3," = ",F10.6)
      ENDDO
    ELSE
      WRITE(out,264) counter
      264 FORMAT (1X,"CI diagonalization procedure converged after ",I3," steps")
      DO i=1,NUM_EIGS
        WRITE(out,265) i, CI_EIG(i)
        265 FORMAT (1X,"Current eigenvalue No. ",I3," = ",F10.6)
      ENDDO
    ENDIF
    WRITE(out,*)

    !close open files
    CLOSE(UNIT=out)
  ENDIF

END SUBROUTINE CI_DAVIDSON_LIU_LOOP

FUNCTION SPACE_DOT_PROD(U,s,MAXITER,N,dimen_a,dimen_b) RESULT(sigma)
  INTEGER,INTENT(IN):: N , dimen_a , dimen_b , MAXITER
  TYPE(EXP_SPACE),DIMENSION(MAXITER),INTENT(IN):: U
  REAL(wp),DIMENSION(MAXITER),INTENT(IN):: s
  REAL(wp),DIMENSION(dimen_a,dimen_b):: sigma

  !local variable
  INTEGER:: i 

  sigma = 0.0_wp

  DO i=1,N
    sigma = sigma  + U(i)%exp_vec*s(i)
  ENDDO 

END FUNCTION SPACE_DOT_PROD

FUNCTION SPACE_PROD(U,sigma,MAXITER,N,dimen_a,dimen_b) RESULT(s)
  INTEGER,INTENT(IN):: dimen_a , dimen_b , N , MAXITER
  TYPE(EXP_SPACE),DIMENSION(MAXITER),INTENT(IN):: U
  REAL(wp),DIMENSION(dimen_a,dimen_b),INTENT(IN):: sigma
  REAL(wp),DIMENSION(MAXITER):: s
  
  !local variable
  INTEGER:: i

  DO i=1,N
    s(i) = SUM(U(i)%exp_vec*sigma)
  ENDDO
END FUNCTION SPACE_PROD

SUBROUTINE DIAG_MATRIX(B,B_vec,B_eig,N,m)
  INTEGER,INTENT(IN):: N , m
  REAL(wp),DIMENSION(N,N),INTENT(INOUT):: B
  REAL(wp),DIMENSION(N,m),INTENT(OUT):: B_vec
  REAL(wp),DIMENSION(m),INTENT(OUT)::B_eig

  !local variables
  REAL(wp),DIMENSION(N):: EIGS

  CALL SYEV(B,EIGS)

  B_eig = EIGS(1:m)  
  B_vec = B(:,1:m)
END SUBROUTINE DIAG_MATRIX

SUBROUTINE WRITE_CI_VECTOR(CI_vector,dimen_a,dimen_b,NUM_EIGS)
  INTEGER,INTENT(IN):: NUM_EIGS
  TYPE(EXP_SPACE),DIMENSION(NUM_EIGS),INTENT(IN):: CI_vector
  INTEGER,INTENT(IN):: dimen_a , dimen_b

  !local variables
  INTEGER:: i,j,n

  OPEN(UNIT=vect,FILE="CI_VECTOR",STATUS="REPLACE",ACTION="WRITE")
  
  DO n=1,NUM_EIGS
    WRITE(vect,*) "CI vector numb." , n
    DO j=1,dimen_b
      DO i=1,dimen_a
        WRITE(vect,*) i , j , CI_vector(n)%exp_vec(i,j)
      ENDDO
    ENDDO
  ENDDO

  CLOSE(UNIT=vect)
    
END SUBROUTINE WRITE_CI_VECTOR

END MODULE CI
