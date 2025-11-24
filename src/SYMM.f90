! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

MODULE SYMM

#include "preproc.inc"
use constants
use qmcfort_pos
use lapack

IMPLICIT NONE

!THRESHOLD FOR VECTOR EQUALITY
REAL(wp),PARAMETER:: dvec=0.0001_wp !dvec=0.0001_wp
REAL(wp),PARAMETER:: delem=0.01_wp

!********************************************************************************
!  Symmmetrically equivalent atoms
!********************************************************************************
TYPE SEA_SET
  INTEGER:: TIP
  INTEGER:: SEALEN                              !Number of atoms in a set
  INTEGER:: distc                               !Corresponding column in DIST_MAT
  REAL(wp),DIMENSION(:,:),ALLOCATABLE:: POS     !Position of each atom 
  REAL(wp),DIMENSION(:),ALLOCATABLE:: M         !Mass for each atom
  REAL(wp),DIMENSION(3):: CM                    !Center of Mass of SEA
  REAL(wp),DIMENSION(3):: IMOM                  !Moment of inertia of SEA
  REAL(wp),DIMENSION(3,3):: IAX                 !Principial axis of SEA
END TYPE SEA_SET

!********************************************************************************
!  Symmetry Element + corresponding Symmetry Operation
!********************************************************************************
TYPE SYMM_EL
  CHARACTER(LEN=5):: NAM
  REAL(wp),DIMENSION(3):: AXIS = 0.0_wp
END TYPE SYMM_EL

TYPE SYMM_GROUP
  CHARACTER(LEN=5):: NAM
  TYPE(SYMM_EL):: ID
  TYPE(SYMM_EL):: INV
  TYPE(SYMM_EL):: MAINROT
  INTEGER:: N_ROT
  INTEGER:: N_ROT2
  INTEGER:: N_REF_ROT
  INTEGER:: N_REF
  TYPE(SYMM_EL),DIMENSION(:),ALLOCATABLE:: ROT
  TYPE(SYMM_EL),DIMENSION(:),ALLOCATABLE:: ROT2
  TYPE(SYMM_EL),DIMENSION(:),ALLOCATABLE:: REF_ROT
  TYPE(SYMM_EL),DIMENSION(:),ALLOCATABLE:: REF
END TYPE SYMM_GROUP

CONTAINS

!********************************************************************************
!  Routine calculates center of Mass of entire molecule or SEA set
!********************************************************************************
SUBROUTINE CALC_CM(N,POS,MASS,CM)
  INTEGER,INTENT(IN):: N
  REAL(wp),DIMENSION(3,N),INTENT(IN):: POS
  REAL(wp),DIMENSION(N),INTENT(IN):: MASS
  REAL(wp),DIMENSION(3),INTENT(OUT):: CM

  !local variables
  INTEGER:: i

  CM = 0.0_wp
  DO i=1,N
    CM = CM + MASS(i)*POS(:,i)
  ENDDO
  
  CM = CM/SUM(MASS)
END SUBROUTINE CALC_CM


!********************************************************************************
!  Calculates the Moment of inertia of entire molecule, i.e. of SEA
!********************************************************************************
SUBROUTINE CALC_IMOM(N,POS,CM,MASS,IMOM,IAX)
  INTEGER,INTENT(IN):: N
  REAL(wp),DIMENSION(3,N),INTENT(IN):: POS
  REAL(wp),DIMENSION(3),INTENT(IN):: CM
  REAL(wp),DIMENSION(N),INTENT(IN):: MASS
  REAL(wp),DIMENSION(3),INTENT(OUT):: IMOM
  REAL(wp),DIMENSION(3,3),INTENT(OUT):: IAX

  !local variables
  INTEGER:: i 
  REAL(wp),DIMENSION(3,3):: ITEN
  REAL(wp),DIMENSION(3,N):: POS2

  DO i=1,N
    POS2(:,i) = POS(:,i) - CM
  ENDDO

  ITEN = 0.0_wp

  DO i=1,N
    ITEN = ITEN + MASS(i)*(DOT_PRODUCT(POS2(:,i),POS2(:,i))*IDENT(3)-TENPROD(3,POS2(:,i),POS2(:,i)))
  ENDDO

  !Daigonalize Moment of Inertia
  CALL SYM_EIGV(3,ITEN,IMOM,IAX)
END SUBROUTINE CALC_IMOM


!********************************************************************************
!  Routine that construct SEA sets from all atoms in molecules
!********************************************************************************
SUBROUTINE FIND_SEA_SETS(STRUC,NUM_SEA,SEA)
  TYPE(STRUCTURE),INTENT(INOUT):: STRUC
  INTEGER,INTENT(OUT):: NUM_SEA
  TYPE(SEA_SET),DIMENSION(:),ALLOCATABLE:: SEA
  
  !local variables
  INTEGER:: i, j, k
  TYPE(SEA_SET),DIMENSION(STRUC%natoms):: TEMP_SEA
  REAL(wp),DIMENSION(STRUC%natoms,STRUC%natoms):: DIST_MAT
  INTEGER:: upk
  LOGICAL:: logres



  !Allocate Temporal and maximal SEA set
  DO i=1,STRUC%natoms
    ALLOCATE(TEMP_SEA(i)%POS(3,STRUC%natoms))
    ALLOCATE(TEMP_SEA(i)%M(STRUC%natoms))
    TEMP_SEA(i)%SEALEN = 0
  ENDDO

  !Calculate distance matrix D(i,j) = |r(i) - r(j)|
  CALL CALC_DIST_MAT()

  !Divide atoms in SEA SETS according to DIST_MAT
  TEMP_SEA(1)%SEALEN = TEMP_SEA(1)%SEALEN + 1
  TEMP_SEA(1)%POS(:,TEMP_SEA(1)%SEALEN) = STRUC%POS(:,1)
  TEMP_SEA(1)%M(TEMP_SEA(1)%SEALEN) = STRUC%MASS(1)
  TEMP_SEA(1)%distc = 1 
  NUM_SEA  = 1

  DO i=2,STRUC%natoms
    upk = NUM_SEA
    DO k=1,upk
      j=TEMP_SEA(k)%distc
      CALL COMPARE(i,j,logres)
  
      !Add atom to existing SEA SET
      IF(logres) THEN
        TEMP_SEA(k)%SEALEN = TEMP_SEA(k)%SEALEN + 1
        TEMP_SEA(k)%POS(:,TEMP_SEA(k)%SEALEN) = STRUC%POS(:,i)
        TEMP_SEA(k)%M(TEMP_SEA(k)%SEALEN) = STRUC%MASS(i)

        EXIT
      ENDIF

      !open new SEA SET
      IF(k==upk) THEN
        NUM_SEA = NUM_SEA + 1
        TEMP_SEA(NUM_SEA)%SEALEN = TEMP_SEA(NUM_SEA)%SEALEN + 1
        TEMP_SEA(NUM_SEA)%POS(:,TEMP_SEA(NUM_SEA)%SEALEN) = STRUC%POS(:,i)
        TEMP_SEA(NUM_SEA)%M(TEMP_SEA(NUM_SEA)%SEALEN) = STRUC%MASS(i)
        TEMP_SEA(NUM_SEA)%distc = i
      ENDIF
    ENDDO
  ENDDO

  !Redefine SEA_SET according to informations contained in TEMP_SEA_SET
  CALL REDEFINE_SEA(NUM_SEA,SEA,TEMP_SEA)

  WRITE(*,*) "Number of SEA SETS", NUM_SEA 
  WRITE(*,*) starshort
  DO i=1,NUM_SEA
    WRITE(*,*)"  Length of SEA SET", SEA(i)%SEALEN
    WRITE(*,*)"  SEA TYPE" , SEA(i)%TIP
    WRITE(*,*)
  ENDDO 
  WRITE(*,*) starshort
CONTAINS

!********************************************************************************
!  Calculates distance matrix D(i,j)=|r(i)-r(j)|
!********************************************************************************
  SUBROUTINE CALC_DIST_MAT()
    INTEGER:: i, j
  
    DIST_MAT = 0.0_wp

    DO j=1,STRUC%natoms
      DO i=1,j-1
        DIST_MAT(i,j) = NORM2(STRUC%POS(:,i)-STRUC%POS(:,j))
      ENDDO
    ENDDO

    DIST_MAT = DIST_MAT + TRANSPOSE(DIST_MAT)
  END SUBROUTINE CALC_DIST_MAT

!********************************************************************************
!  Compare two columns of distance matrix 
!       INPUT:
!               i , j  - indices of columns 
!       OUTPUT:
!               logres - TRUE if columns are same
!********************************************************************************
  SUBROUTINE COMPARE(i,j,logres)
    INTEGER,INTENT(IN):: i, j 
    LOGICAL,INTENT(OUT):: logres

    !local variables
    REAL(wp),DIMENSION(STRUC%natoms):: vec1 , vec2
    REAL(wp):: swap
    INTEGER:: match, match_old
    INTEGER:: k,l
 
    vec1 = DIST_MAT(:,i)
    vec2 = DIST_MAT(:,j)

    match = 0
    DO k=1,STRUC%natoms
      match_old = match
      DO l=match_old+1,STRUC%natoms
        IF(ABS(vec1(k)-vec2(l))<dvec) THEN
          match = match + 1 
          swap = vec2(l)
          vec2(l) = vec2(match)
          vec2(match) = swap
          EXIT
        ENDIF
      ENDDO

      IF(match_old==match) THEN
        logres = .FALSE.
        RETURN
      ENDIF
    ENDDO

    IF(match==STRUC%natoms) THEN
      logres = .TRUE.
    ELSE
      WRITE(*,*)"Something strange in COMPARE(i,j)- logres set to false"
      logres = .FALSE.
    ENDIF
  END SUBROUTINE COMPARE

!********************************************************************************
!  Assign final SEA set from temporal one TEMP_SEA
!  Calculates also moment of inertia and center of 
!  mass for each SEA_SET and for whole molecule
!********************************************************************************
  SUBROUTINE REDEFINE_SEA(NUM_SEA,SEA,TEMP_SEA)
    INTEGER,INTENT(IN):: NUM_SEA
    TYPE(SEA_SET),DIMENSION(:),ALLOCATABLE,INTENT(OUT):: SEA
    TYPE(SEA_SET),DIMENSION(STRUC%natoms):: TEMP_SEA

    !local variables
    INTEGER:: i 

    ALLOCATE(SEA(NUM_SEA))
    DO i=1,NUM_SEA
      SEA(i)%SEALEN = TEMP_SEA(i)%SEALEN
      SEA(i)%distc = TEMP_SEA(i)%distc
      ALLOCATE(SEA(i)%POS(3,SEA(i)%SEALEN))
      ALLOCATE(SEA(i)%M(SEA(i)%SEALEN))
      
      SEA(i)%POS = TEMP_SEA(i)%POS(:,1:TEMP_SEA(i)%SEALEN)
      SEA(i)%M = TEMP_SEA(i)%M(1:TEMP_SEA(i)%SEALEN)

      !Calculate CM and Moment of Inertia for each SEATIP 
      CALL CALC_CM(SEA(i)%SEALEN,SEA(i)%POS,SEA(i)%M,SEA(i)%CM)
      CALL CALC_IMOM(SEA(i)%SEALEN,SEA(i)%POS,SEA(i)%CM,SEA(i)%M,SEA(i)%IMOM,SEA(i)%IAX)

      !Determine 
      CALL DET_SEA_TYPE(SEA(i)%SEALEN,SEA(i)%IMOM,SEA(i)%TIP)
    ENDDO
  END SUBROUTINE REDEFINE_SEA

!********************************************************************************
!  Determine Type of SEA Set
!       INPUT:
!               SEALEN - length of SEA Set
!               IMOM - moment of inertia
!       OUTPUT:
!               SEATIP :
!                 1  - Single Atom
!                       No informations about symmetry
!                 2  - Linear arrangemenet
!                       any Cn(Ic) and perpendicular C2
!                 31 - Regular polygonal arrang
!                 32 - Irregular polygonal arrang
!                 41 - Prolate symmetric polyhedron
!                 42 - Oblate symmetric polyhedron
!                 43 - primitive tetragonal cell
!                 0  - ERROR
!********************************************************************************
  SUBROUTINE DET_SEA_TYPE(SEALEN,IMOM,SEATIP)
    INTEGER,INTENT(IN):: SEALEN
    REAL(wp),DIMENSION(3),INTENT(IN):: IMOM
    INTEGER,INTENT(OUT):: SEATIP
 
    !local variables

    IF(SEALEN==1) THEN
      !Single atom
      SEATIP = 1
    ELSEIF(SEALEN==2) THEN
      !Linear arrangement
      IF(IMOM(1)<dvec .and. ABS(IMOM(2)-IMOM(3))<dvec) THEN
        SEATIP = 2 
      ELSE
        SEATIP = 2 
        WRITE(*,*) "Warning: untypical situation - analyse of I failed"
      ENDIF
    ELSE
      IF(ABS(IMOM(1)+IMOM(2)-IMOM(3))<delem) THEN
        !Polygonal arrangement
        IF(ABS(IMOM(1)-IMOM(2))<delem) THEN
          !Regular polygon
          SEATIP = 31
        ELSE
          !Irregular polygon
          SEATIP = 32
        ENDIF
      ELSEIF(ABS(IMOM(1)-IMOM(2))>delem .and. ABS(IMOM(2)-IMOM(3))>delem .and. ABS(IMOM(1)-IMOM(3))>delem) THEN
        !Polyhedrar Arrangement
        IF(ABS(IMOM(1)-IMOM(2))>delem .and. ABS(IMOM(2)-IMOM(3))<delem) THEN
          !Prolate Symmetric Top
          SEATIP = 41
        ELSEIF(ABS(IMOM(1)-IMOM(2))<delem .and. ABS(IMOM(2)-IMOM(3))>delem) THEN
          !Oblate Symmetric Top
          SEATIP = 42
        ELSE
          !primitive orthorombic cell
          SEATIP = 43
        ENDIF
      ELSE
        WRITE(*,*) "Error can be expected - the SEA TIP not recognised "
        SEATIP = 0
      ENDIF
    ENDIF
  END SUBROUTINE DET_SEA_TYPE
END SUBROUTINE FIND_SEA_SETS


!********************************************************************************
!  First main subroutine in this module - determines symmetry point group of a 
!                  molecule and all symmetry elements
!********************************************************************************
SUBROUTINE SYMM_MOL(STRUC,SYMMS)
  TYPE(STRUCTURE),INTENT(INOUT):: STRUC
  TYPE(SYMM_GROUP),INTENT(OUT):: SYMMS

  !local varaibles
  REAL(wp),DIMENSION(3,STRUC%natoms):: TEMP_POS
  TYPE(SYMM_EL),DIMENSION(100):: TEMP_EL
  TYPE(SEA_SET),DIMENSION(:),ALLOCATABLE:: SEA
  TYPE(SYMM_EL):: EL
  INTEGER:: NUM_SEA ,  NUM_SYMMS = 0 , NUM_FIX_PTS = 0
  REAL(wp),DIMENSION(:,:),ALLOCATABLE:: FIX_PTS
  INTEGER:: i   
  REAL(wp),DIMENSION(3):: pt
  LOGICAL:: res , PLAN

  !Assign masses of each atom in STRUC according to its Z number
  call struc%set_atomic_mass
 
  CALL CALC_CM(STRUC%natoms,STRUC%POS,STRUC%MASS,STRUC%CM)

  !Redefine positions in respect to Center of mass: CM =(0 0 0)
  DO i=1,STRUC%natoms
    STRUC%POS(:,i) = STRUC%POS(:,i) - STRUC%CM
  ENDDO
  
  CALL CALC_CM(STRUC%natoms,STRUC%POS,STRUC%MASS,STRUC%CM)

  TEMP_POS = STRUC%POS

  CALL CALC_IMOM(STRUC%natoms,TEMP_POS,STRUC%CM,STRUC%MASS,STRUC%IMOM,STRUC%IAX)
  
  CALL DET_ROT_TYPE(STRUC%natoms,STRUC%IMOM,STRUC%ROT_TIP,PLAN)

  CALL SORT_POS(TEMP_POS,STRUC%natoms)
  
  IF(STRUC%ROT_TIP==-1) THEN
    !Single Atom
    WRITE(*,*) "The structure corresponds to single atom. Symmetries are not found"
    RETURN
  ELSEIF(STRUC%ROT_TIP==0) THEN
    !Error
    WRITE(*,*) "Rotor type of molecule is not recognised "
    RETURN
  ELSEIF(STRUC%ROT_TIP==1) THEN
    !Linear Molecule
    WRITE(*,*) "The molecule is linear"
    EL%NAM = 'cb'
    EL%AXIS = STRUC%IAX(:,1)
    CALL ADD_EL(EL)

    SYMMS%N_ROT = NUM_SYMMS
    ALLOCATE(SYMMS%ROT(SYMMS%N_ROT))
    SYMMS%ROT(1) = TEMP_EL(1)
    SYMMS%MAINROT = TEMP_EL(1)
    !Update Num_Symms
    NUM_SYMMS = 0
  ELSEIF(STRUC%ROT_TIP==2) THEN
    !Spherical Molecule
    WRITE(*,*) "Still not implemented for spherical molecules"
    RETURN
  ELSEIF(STRUC%ROT_TIP==3 .or. STRUC%ROT_TIP==4) THEN
    !Symmetric or Asymmetric molecule
  
    !Divide atom positions to SEA sets
    CALL FIND_SEA_SETS(STRUC,NUM_SEA,SEA)
 
    !Sort positions in SEA sets - to be able to find symm elements
    DO i=1,NUM_SEA
      CALL SORT_POS(SEA(i)%POS,SEA(i)%SEALEN)
    ENDDO

    ALLOCATE(FIX_PTS(3,NUM_SEA))
    !Add center of mass of entire molecule to fix point
    NUM_FIX_PTS = NUM_FIX_PTS + 1
    FIX_PTS(:,NUM_FIX_PTS) = STRUC%CM  
    !Do only once, and not for all symmetry operations
    DO i=1,NUM_SEA
      IF(SEA(i)%SEALEN==1) THEN
        pt = SEA(i)%CM
        CALL COMPARE_PT_SET(pt,FIX_PTS,NUM_FIX_PTS,res)
        IF(.not. res) THEN
          NUM_FIX_PTS = NUM_FIX_PTS + 1 
          FIX_PTS(:,NUM_FIX_PTS) = pt
        ENDIF
      ENDIF  
    ENDDO
 
    CALL FIND_ROT()
  
    CALL FIND_C2_ROT()

    CALL FIND_REF_ROT()

    CALL FIND_REF()
   
  ENDIF

  !Do for all ROTOR TYPES: inversion(i) and identity(E)
  !Check for inversion element
  CALL FIND_INV(res)
  SYMMS%INV%NAM = 'i'
  IF(res) THEN
    SYMMS%INV%AXIS = 1.0_wp
  ELSE
    SYMMS%INV%AXIS = 0.0_wp
  ENDIF

  !Add identity 
  SYMMS%ID%NAM = 'e'
  SYMMS%ID%AXIS = 1.0_wp

  WRITE(*,*) starshort
  WRITE(*,*) "Identity               " ,  .TRUE.
  WRITE(*,*) "Inversion              " ,  res
  WRITE(*,*) "Number of Rotations    " ,  SYMMS%N_ROT
  WRITE(*,*) "Number of C2 Rotations " ,  SYMMS%N_ROT2
  WRITE(*,*) "Number of Improper Rot " ,  SYMMS%N_REF_ROT
  WRITE(*,*) "Number of Reflections  " ,  SYMMS%N_REF
  CALL DET_PT_GR(STRUC%ROT_TIP,SYMMS)
  WRITE(*,*)
  WRITE(*,*) "MAINROT " , SYMMS%MAINROT%NAM , " :" , SYMMS%MAINROT%AXIS
  WRITE(*,*)
  WRITE(*,*) "Point group is " , SYMMS%NAM
  WRITE(*,*) starshort

CONTAINS

!********************************************************************************
!  Find all possible rotations (along moment of inertia) add them 
!           to SYMMS%ROT  and restart counter NUM_SYMMS
!********************************************************************************
  SUBROUTINE FIND_ROT()
    LOGICAL:: res 
    INTEGER:: i , j  
    INTEGER:: nmax , m
    INTEGER:: cmax = 1 , cmaxind
    REAL(wp),DIMENSION(3,3):: MAT
    REAL(wp),DIMENSION(3):: ax
    CHARACTER(len=10):: numb

    DO i=1,NUM_SEA
      IF(SEA(i)%TIP==1) THEN
        !No new information about symmetry
        CYCLE
      ELSEIF(SEA(i)%TIP==2) THEN
        !Add Cb rottaion along Ia axis and perpendicular C2 rotations
        !Should not be imortant
        CYCLE
        EL%NAM = "cb"
        EL%AXIS = SEA(i)%IAX(:,1)
        CALL ADD_EL(EL)      
        cmax = 2 
        !In this Moment no additional informations about C2 axes
      ELSEIF(SEA(i)%TIP==31) THEN
        !regular polygon
        nmax = SEA(i)%SEALEN
        ax = SEA(i)%IAX(:,3)
        DO j=1,nmax-1
          IF(MOD(nmax,j)==0) THEN
            EL%NAM = "c" 
            m = nmax/j
            MAT = ROT_MAT(m,ax)
            CALL SYMM_ACT(MAT,TEMP_POS,STRUC%natoms,res)
            IF(res) THEN
              EL%AXIS = ax
              IF(m>9) THEN
                WRITE(numb,'(I2)') m
              ELSE
                WRITE(numb,'(I1)') m 
              ENDIF
              !numb = adjustl(numb)
              EL%NAM = trim(EL%NAM) // trim(numb)                
              CALL ADD_EL(EL)
              IF(m > cmax) THEN
                cmax = m 
                cmaxind = NUM_SYMMS
              ENDIF
            ENDIF
          ENDIF
        ENDDO
      ELSEIF(SEA(i)%TIP==32) THEN
        !iregular polygon 
        nmax = SEA(i)%SEALEN
        ax = SEA(i)%IAX(:,3)
        DO j=2,nmax-1
          IF(MOD(nmax,j)==0) THEN
            EL%NAM = "c"
            m = nmax/j
            MAT = ROT_MAT(m,ax)
            CALL SYMM_ACT(MAT,TEMP_POS,STRUC%natoms,res)
            IF(res) THEN
              EL%AXIS = ax
              IF(m>9) THEN
                WRITE(numb,'(I2)') m
              ELSE
                WRITE(numb,'(I1)') m 
              ENDIF
              numb = adjustl(numb)
              EL%NAM = trim(EL%NAM) // trim(numb)
              CALL ADD_EL(EL)
              IF(m > cmax) THEN
                cmax = m 
                cmaxind = NUM_SYMMS
              ENDIF
            ENDIF
          ENDIF
        ENDDO
      ELSEIF(SEA(i)%TIP==41) THEN
        !Prolate Symmetric Top
        IF(MOD(SEA(i)%SEALEN,2)==0) THEN
          nmax = SEA(i)%SEALEN/2
        ELSE
          WRITE(*,*) "POLYHEDRON HAS ODD NUMBER OF VERTICES - unrecognised"
          CYCLE
        ENDIF

        ax = SEA(i)%IAX(:,1)
        DO j=1,nmax-1
          IF(MOD(nmax,j)==0) THEN
            EL%NAM = "c"
            m = nmax/j
            MAT = ROT_MAT(m,ax)
            CALL SYMM_ACT(MAT,TEMP_POS,STRUC%natoms,res)
            IF(res) THEN
              EL%AXIS = ax
              IF(m>9) THEN
                WRITE(numb,'(I2)') m
              ELSE
                WRITE(numb,'(I1)') m 
              ENDIF
              numb = adjustl(numb)
              EL%NAM = trim(EL%NAM) // trim(numb)
              CALL ADD_EL(EL)
              IF(m > cmax) THEN
                cmax = m 
                cmaxind = NUM_SYMMS
              ENDIF
            ENDIF
          ENDIF
        ENDDO
      ELSEIF(SEA(i)%TIP==42) THEN
        !Oblate Symmetric Top
        IF(MOD(SEA(i)%SEALEN,2)==0) THEN
          nmax = SEA(i)%SEALEN/2
        ELSE
          WRITE(*,*) "POLYHEDRON HAS ODD NUMBER OF VERTICES - unrecognised!"
          WRITE(*,*) "Switch to new SEA set!"
          CYCLE
        ENDIF

        ax = SEA(i)%IAX(:,3)
        DO j=1,nmax-1
          IF(MOD(nmax,j)==0) THEN
            EL%NAM = "c"
            m = nmax/j
            MAT = ROT_MAT(m,ax)
            CALL SYMM_ACT(MAT,TEMP_POS,STRUC%natoms,res)
            IF(res) THEN
              EL%AXIS = ax
              IF(m>9) THEN
                WRITE(numb,'(I2)') m
              ELSE
                WRITE(numb,'(I1)') m 
              ENDIF
              numb = adjustl(numb)
              EL%NAM = trim(EL%NAM) // trim(numb)
              CALL ADD_EL(EL)
              IF(m > cmax) THEN
                cmax = m 
                cmaxind = NUM_SYMMS
              ENDIF
            ENDIF
          ENDIF
        ENDDO
      ELSEIF(SEA(i)%TIP==43) THEN
        EL%NAM = 'c2'
        DO j=1,3
          ax = SEA(i)%IAX(:,j)
          MAT = ROT_MAT(2,ax)
          CALL SYMM_ACT(MAT,TEMP_POS,STRUC%natoms,res)
          IF(res) THEN
            EL%AXIS = ax
            CALL ADD_EL(EL)
            IF(m > cmax) THEN
              cmax = m 
              cmaxind = NUM_SYMMS
            ENDIF
          ENDIF
        ENDDO
      ELSE
        WRITE(*,*) "Type of SEA set unrecognised" 
      ENDIF
    ENDDO

    SYMMS%N_ROT = NUM_SYMMS
    IF(NUM_SYMMS>0) THEN
      ALLOCATE(SYMMS%ROT(SYMMS%N_ROT))
      DO i=1,SYMMS%N_ROT
        SYMMS%ROT(i) = TEMP_EL(i)
      ENDDO
      SYMMS%MAINROT = TEMP_EL(cmaxind)
    ELSE
      SYMMS%MAINROT%NAM = "c1"
      SYMMS%MAINROT%AXIS = 1.0_wp
    ENDIF

    !Re-update TEMP_EL
    NUM_SYMMS = 0 
  END SUBROUTINE FIND_ROT

!********************************************************************************
!  Find all planes of reflection and add them to 
!     SYMMS%REF and restart counter NUM_SYMMS  
!********************************************************************************
  SUBROUTINE FIND_REF()
    !local varaibles
    INTEGER:: i , j , k
    REAL(wp),DIMENSION(3):: ax , main_rot_ax, vec1 , vec2
    REAL(wp),DIMENSION(3,3):: MAT
    LOGICAL:: res2
    CHARACTER(LEN=1):: ext

    main_rot_ax = SYMMS%MAINROT%AXIS
    
    DO i=1,NUM_SEA
      IF(SEA(i)%SEALEN>1) THEN
        DO j=1,SEA(i)%SEALEN
          DO k=j+1,SEA(i)%SEALEN
            vec1 = SEA(i)%POS(:,j)
            vec2 = SEA(i)%POS(:,k)

            EL%NAM = 'sig' 
            ax = (vec2-vec1)/NORM2(vec2-vec1)
            MAT = REF_MAT(ax)
        
            CALL SYMM_ACT(MAT,TEMP_POS,STRUC%natoms,res)
          
            IF(res) THEN
              CALL CHECK_NAME(ax,main_rot_ax,ext)
              EL%NAM = trim(EL%NAM) // trim(ext)
              EL%AXIS = ax
              CALL CHECK_SYMM(EL,NUM_SYMMS,TEMP_EL,res2)
              IF(.not. res2) THEN
                NUM_SYMMS = NUM_SYMMS + 1 
                TEMP_EL(NUM_SYMMS)%NAM = EL%NAM
                TEMP_EL(NUM_SYMMS)%AXIS = EL%AXIS
              ENDIF
            ENDIF

          ENDDO
        ENDDO
      ENDIF
    ENDDO

    IF(PLAN) THEN
      EL%NAM = 'sig'
      DO i=2,STRUC%natoms
        ax = CROSSPROD(STRUC%POS(:,1),STRUC%POS(:,i))
        IF(NORM2(ax)>dvec) THEN
          ax = ax/NORM2(ax)
          EXIT  
        ENDIF
      ENDDO
      MAT = REF_MAT(ax)
      
      CALL SYMM_ACT(MAT,TEMP_POS,STRUC%natoms,res)
      
      IF(res) THEN
        CALL CHECK_NAME(ax,main_rot_ax,ext)
        EL%NAM = trim(EL%NAM) // trim(ext)
        EL%AXIS = ax
        CALL CHECK_SYMM(EL,NUM_SYMMS,TEMP_EL,res2)
        IF(.not. res2) THEN
          NUM_SYMMS = NUM_SYMMS + 1 
          TEMP_EL(NUM_SYMMS)%NAM = EL%NAM
          TEMP_EL(NUM_SYMMS)%AXIS = EL%AXIS
        ENDIF
      ENDIF
    ENDIF
 
    SYMMS%N_REF = NUM_SYMMS
    IF(NUM_SYMMS>0) THEN
      ALLOCATE(SYMMS%REF(SYMMS%N_REF))
      DO i=1,SYMMS%N_REF
        SYMMS%REF(i) = TEMP_EL(i)
      ENDDO
    ENDIF
    
    NUM_SYMMS = 0
    
  END SUBROUTINE FIND_REF

!********************************************************************************
!  Check name of a symmetry plane, used in FIND_REF()
!       h - horizontal
!       v - vertical
!       d - dihedral
!********************************************************************************
  SUBROUTINE CHECK_NAME(ax,mainax,ext)
    REAL(wp),DIMENSION(3):: ax , mainax
    CHARACTER(LEN=1):: ext
    INTEGER:: i, j
    !local variables
    REAL(wp):: prod
    REAL(wp),DIMENSION(3):: helpax

    prod = DOT_PRODUCT(ax,mainax)

    IF(ABS(prod)<dvec) THEN
      ext = 'v'
      outer: DO i=1,SYMMS%N_ROT2-1
        inner: DO j=i+1,SYMMS%N_ROT2
          helpax = SYMMS%ROT2(i)%AXIS - SYMMS%ROT2(j)%AXIS
          helpax = helpax/NORM2(helpax)
          prod = DOT_PRODUCT(helpax,ax)
          IF(ABS(prod)<dvec) THEN
            ext = 'd'
            EXIT outer
          ELSEIF(ABS(ABS(prod)-1.0_wp)<dvec) THEN
            ext = 'd'
            EXIT outer
          ENDIF
        ENDDO inner
      ENDDO outer        
    ELSEIF(ABS(ABS(prod)-1.0_wp)<dvec) THEN
      ext = 'h'
    ELSE
      ext = 'x'
    ENDIF
  END SUBROUTINE CHECK_NAME

!********************************************************************************
!  Find all possible C2 rotations (they cant be found from moment of inertia)
!       They are added to SYMMS%ROT2 and NUM_SYMMS is restarted
!********************************************************************************
  SUBROUTINE FIND_C2_ROT()
    !local variables
    INTEGER:: i , j , k
    REAL(wp),DIMENSION(3):: ax , vec1 , vec2
    REAL(wp),DIMENSION(3,3):: MAT
    LOGICAL:: res, ress, res2 , ress1 , ress2
    INTEGER:: ind1 , ind2 ,  DI_numb
    INTEGER,DIMENSION(:),ALLOCATABLE:: DI_indices
 
    EL%NAM = 'c2'
    WRITE(*,*) "NUM FIX PTS= " , NUM_FIX_PTS
    IF(NUM_FIX_PTS<=1) THEN
      !Middlepoint of two atoms in SEA SET criterion
      DO i=1,NUM_SEA
        DO j=1,SEA(i)%SEALEN
          DO k=1,j-1
            ax = 0.5_wp*(SEA(i)%POS(:,j)+SEA(i)%POS(:,k))
            ax = ax/NORM2(ax)
            CALL COMPARE_VECS(ax,STRUC%CM,res)
            CALL COMPARE_VECS(ax,SYMMS%MAINROT%AXIS,ress1)
            CALL COMPARE_VECS(ax,-SYMMS%MAINROT%AXIS,ress2)
            res = .not. res
            ress1 = .not. ress1
            ress2 = .not. ress2
            ress = ress1 .and. ress2
            IF(res .and. ress) THEN
              MAT = ROT_MAT(2,ax)
              CALL SYMM_ACT(MAT,TEMP_POS,STRUC%natoms,res2)
              IF(res2) THEN
                EL%AXIS = ax
                CALL ADD_EL(EL)
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      !An atom in SEA set criterion
      DO i=1,NUM_SEA
        DO j=1,SEA(i)%SEALEN
          ax = SEA(i)%POS(:,j)
          ax = ax/NORM2(ax)
          CALL COMPARE_VECS(ax,STRUC%CM,res)
          CALL COMPARE_VECS(ax,SYMMS%MAINROT%AXIS,ress1)
          CALL COMPARE_VECS(ax,-SYMMS%MAINROT%AXIS,ress2)
          res = .not. res
          ress1 = .not. ress1
          ress2 = .not. ress2
          ress = ress1 .and. ress2
          IF(res .and. ress) THEN
            MAT = ROT_MAT(2,ax)
            CALL SYMM_ACT(MAT,TEMP_POS,STRUC%natoms,res2)
            IF(res2) THEN
              EL%AXIS = ax
              CALL ADD_EL(EL)
            ENDIF
          ENDIF
        ENDDO
      ENDDO

      !More diatomic SEA sets criterion
      CALL FIND_DI_SEA(DI_numb,DI_indices)
      IF(DI_numb>1) THEN
        DO i=1,DI_numb
          DO j=1,i-1
            ind1 = DI_indices(i)
            ind2 = DI_indices(j)
            vec1 = SEA(ind1)%POS(:,1) - SEA(ind1)%POS(:,2)
            vec2 = SEA(ind2)%POS(:,1) - SEA(ind2)%POS(:,2)

            ax = CROSSPROD(vec1,vec2)
            ax = ax/NORM2(ax)
            MAT = ROT_MAT(2,ax)
            CALL SYMM_ACT(MAT,TEMP_POS,STRUC%natoms,res)
            IF(res) THEN
              EL%AXIS = ax
              CALL ADD_EL(EL)
            ENDIF

          ENDDO
        ENDDO
      ENDIF
    ELSE
      ax = FIX_PTS(:,2) - FIX_PTS(:,1)
      ax = ax/NORM2(ax)
      MAT = ROT_MAT(2,ax)
      CALL SYMM_ACT(MAT,TEMP_POS,STRUC%natoms,res)
      IF(res) THEN
        EL%AXIS = ax
        CALL ADD_EL(EL)
      ENDIF
    ENDIF

    SYMMS%N_ROT2 = NUM_SYMMS
    IF(NUM_SYMMS>0) THEN
      ALLOCATE(SYMMS%ROT2(SYMMS%N_ROT2))
      DO i=1,SYMMS%N_ROT2
        SYMMS%ROT2(i) = TEMP_EL(i)
      ENDDO
    ENDIF

    IF(ABS(SUM(SYMMS%MAINROT%AXIS)-3.0_wp)<dvec) THEN
      IF(SYMMS%N_ROT2==1) THEN
       SYMMS%MAINROT = SYMMS%ROT2(1)
      ENDIF
    ENDIF

    NUM_SYMMS = 0

  END SUBROUTINE FIND_C2_ROT

!********************************************************************************
!  Calculates number and indices of diatomic SEA sets, used in FINC_C2_ROT()
!********************************************************************************
  SUBROUTINE FIND_DI_SEA(numb,indices)
     INTEGER,INTENT(OUT):: numb
     INTEGER,DIMENSION(:),ALLOCATABLE,INTENT(OUT):: indices

     !local variables
     INTEGER:: i 
     INTEGER:: tempnumb=0
     INTEGER,DIMENSION(NUM_SEA):: tempindices

     DO i=1,NUM_SEA
       IF(SEA(i)%SEALEN==2) THEN
         tempnumb = tempnumb + 1
         tempindices(tempnumb) = i
       ENDIF
     ENDDO

     numb = tempnumb
     ALLOCATE(indices(numb))

     DO i=1,numb
       indices(i) = tempindices(i)
     ENDDO
  END SUBROUTINE FIND_DI_SEA

!********************************************************************************
!  Find all improper axes of rotations and store them to SYMMS%REF_ROT
!********************************************************************************
  SUBROUTINE FIND_REF_ROT()
    INTEGER:: i , j , m , nm
    REAL(wp),DIMENSION(3):: ax
    REAL(wp),DIMENSION(3,3):: MAT
    CHARACTER(LEN=10):: numb
    LOGICAL:: res
    
    DO i=1,SYMMS%N_ROT
      EL = SYMMS%ROT(i)
      ax = EL%AXIS
      READ(EL%NAM(2:),*) m
      DO j=1,2
        nm = m * j
        MAT = REF_ROT_MAT(nm,ax)
        CALL SYMM_ACT(MAT,TEMP_POS,STRUC%natoms,res)
        IF(res) THEN
          IF(nm>9) THEN
            WRITE(numb,'(I2)') nm
          ELSE
            WRITE(numb,'(I1)') nm 
          ENDIF
          numb = adjustl(numb)
          EL%NAM  = 's' // trim(numb)
          CALL ADD_EL(EL)
        ENDIF
      ENDDO      
    ENDDO

    SYMMS%N_REF_ROT = NUM_SYMMS
    IF(NUM_SYMMS>0) THEN
      ALLOCATE(SYMMS%REF_ROT(SYMMS%N_REF_ROT))
      DO i=1,SYMMS%N_REF_ROT
        SYMMS%REF_ROT(i) = TEMP_EL(i)
      ENDDO
    ENDIF
    
    NUM_SYMMS = 0    
  END SUBROUTINE FIND_REF_ROT

!********************************************************************************
!  Check if the molecule has center of inversion
!********************************************************************************
  SUBROUTINE FIND_INV(res)
    LOGICAL,INTENT(OUT):: res

    !local variables
    REAL(wp):: crit

    crit = SUM(TEMP_POS)
    IF(ABS(crit)<dvec) THEN
      res = .TRUE.
    ELSE
      res = .FALSE.
    ENDIF
  END SUBROUTINE FIND_INV

!********************************************************************************
!  Add Symmetry element to TEMP_EL, if it is not already in TEMP_EL
!********************************************************************************
  SUBROUTINE ADD_EL(EL)
    TYPE(SYMM_EL),INTENT(IN):: EL
    LOGICAL:: res

    CALL CHECK_SYMM(EL,NUM_SYMMS,TEMP_EL,res)
    IF(.not. res) THEN
      NUM_SYMMS = NUM_SYMMS + 1 
      TEMP_EL(NUM_SYMMS)%NAM = EL%NAM
      TEMP_EL(NUM_SYMMS)%AXIS = EL%AXIS
    ENDIF
  END SUBROUTINE ADD_EL
END SUBROUTINE SYMM_MOL

!********************************************************************************
!  Determines the symmetry point group of molecule
!********************************************************************************
SUBROUTINE DET_PT_GR(ROT_TIP,SYMMS)
  INTEGER,INTENT(IN):: ROT_TIP
  TYPE(SYMM_GROUP),INTENT(INOUT):: SYMMS

  !local variables
  INTEGER:: N_ROT5
  LOGICAL:: highC , sigh , s2n
  INTEGER:: maxC , nc2 , nsigd , nsigv
  INTEGER:: i , tempint
  REAL(wp),DIMENSION(3):: ax1 , ax2
  REAL(wp):: prod
  CHARACTER(LEN=2) :: ext , ext2
  !Find maxC from MAINROT
  IF(TRIM(SYMMS%MAINROT%NAM)/='cb') THEN
    READ(SYMMS%MAINROT%NAM(2:),*) maxC
    IF(maxC>1) THEN
      highC = .TRUE.
    ELSE
      highC = .FALSE.
    ENDIF
    ax1 = SYMMS%MAINROT%AXIS
    ext = SYMMS%MAINROT%NAM(2:3)

    !Find number of perpendicular C2 to mainrot
    nc2 = 0
    DO i=1,SYMMS%N_ROT2
      ax2 = SYMMS%ROT2(i)%AXIS
      prod = DOT_PRODUCT(ax1,ax2)
      IF(abs(prod)<dvec) THEN
        nc2 = nc2 + 1
      ENDIF
    ENDDO

    !Find if horizontal plane exists and how many dihedral plane exist
    sigh = .FALSE.
    nsigd = 0
    nsigv = 0
    DO i=1,SYMMS%N_REF
      IF(SYMMS%REF(i)%NAM=='sigh') THEN
        sigh = .TRUE.
      ENDIF

      IF(SYMMS%REF(i)%NAM=='sigd') THEN
        nsigd = nsigd + 1
      ENDIF

      IF(SYMMS%REF(i)%NAM=='sigv') THEN
        nsigv = nsigv + 1
      ENDIF
    ENDDO

    !Does S2n exist?
    s2n=.FALSE.
    DO i=1,SYMMS%N_REF_ROT
      READ(SYMMS%REF_ROT(i)%NAM(2:),*) tempint
      IF(tempint==2*maxC) THEN
        s2n = .TRUE.
        ext2 = SYMMS%REF_ROT(i)%NAM(2:3)
      ENDIF
    ENDDO
  ENDIF

  IF(ROT_TIP==1) THEN
    IF(ABS(SUM(SYMMS%INV%AXIS)-3.0_wp)<dvec) THEN
      SYMMS%NAM = 'dinfh'
    ELSE
      SYMMS%NAM = 'cinfv'
    ENDIF
  ELSEIF(ROT_TIP==2) THEN
    IF(ABS(SUM(SYMMS%INV%AXIS)-3.0_wp)<dvec) THEN
      IF(N_ROT5>1) THEN
        SYMMS%NAM = 'ih'
      ELSE
        SYMMS%NAM = 'oh'
      ENDIF
    ELSE
      SYMMS%NAM = 'td'
    ENDIF
  ELSEIF(ROT_TIP==3 .or. ROT_TIP==4) THEN
    IF(highC) THEN
      IF(maxC==nc2) THEN
        IF(sigh) THEN
          SYMMS%NAM = 'd' // TRIM(ext) // 'h'
        ELSE
          IF(nsigd==maxC) THEN
            SYMMS%NAM = 'd' // TRIM(ext) // 'd'
          ELSE
            SYMMS%NAM = 'd' // TRIM(ext)
          ENDIF
        ENDIF
      ELSE
        IF(sigh) THEN
          SYMMS%NAM = 'c' // TRIM(ext) // 'h'
        ELSE
          IF(nsigv==maxc) THEN
            SYMMS%NAM = 'c' // TRIM(ext) // 'v'
          ELSE
            IF(s2n) THEN
              SYMMS%NAM = 's' // ext2
            ELSE
              SYMMS%NAM = 'c' // ext
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ELSE
      IF(SYMMS%N_REF>=1) THEN
        SYMMS%NAM = 'cs'
      ELSE
        IF(ABS(SUM(SYMMS%INV%AXIS)-3.0_wp)<dvec) THEN
          SYMMS%NAM = 'ci'
        ELSE
          SYMMS%NAM = 'c1'
        ENDIF
      ENDIF
    ENDIF
  ELSE
    WRITE(*,*) "Point group can not be determined!"
  ENDIF 
END SUBROUTINE DET_PT_GR


!********************************************************************************
!  Check if the symmetry EL is contained in SYM_SET
!       OUTPUT:
!         res: TRUE OR FALSE
!         ind: if TRUE  -> index of EL in SYM_SET
!              else     -> SYM_SIZE + 1 
!********************************************************************************
SUBROUTINE CHECK_SYMM(EL,SYM_SIZE,SYM_SET,res)
  TYPE(SYMM_EL),INTENT(IN):: EL
  INTEGER,INTENT(IN):: SYM_SIZE
  TYPE(SYMM_EL),DIMENSION(SYM_SIZE),INTENT(IN):: SYM_SET 
  LOGICAL,INTENT(OUT):: res

  !logical varaibles
  INTEGER:: i 
  LOGICAL:: res2 , res3
  
  res = .FALSE.

  DO i=1,SYM_SIZE
    IF(EL%NAM==SYM_SET(i)%NAM) THEN
      CALL COMPARE_VECS(EL%AXIS,SYM_SET(i)%AXIS,res2)
      CALL COMPARE_VECS(-EL%AXIS,SYM_SET(i)%AXIS,res3)
      IF(res2 .or. res3) THEN
        res = .TRUE.
        EXIT
      ENDIF
    ENDIF
  ENDDO

END SUBROUTINE CHECK_SYMM


!********************************************************************************
!  Compare two vectors 
!       PREC: dvec
!       IF they are equivalent -> TRUE , otherwise FALSE
!********************************************************************************
SUBROUTINE COMPARE_VECS(vec1,vec2,res)
  REAL(wp),DIMENSION(3),INTENT(IN):: vec1, vec2
  LOGICAL,INTENT(OUT):: res
  
  !local variables
  REAL(wp):: dist
  dist = SUM((vec1-vec2)**2)
  IF (dist<=dvec) THEN
    res = .TRUE.
  ELSE
    res = .FALSE.
  ENDIF
END SUBROUTINE COMPARE_VECS


!********************************************************************************
!  Compare if the vector vec is already contained in a set of vectors set
!       PREC: SEE -> COMPARE_VECS -> dvec
!********************************************************************************
SUBROUTINE COMPARE_PT_SET(vec,set,set_len,res)
  REAL(wp),DIMENSION(3),INTENT(IN):: vec 
  INTEGER,INTENT(IN):: set_len
  REAL(wp),DIMENSION(3,set_len):: set
  LOGICAL,INTENT(OUT):: res

  !local variables
  INTEGER:: i
  LOGICAL:: tempres
  
  res = .TRUE.
  
  DO i=1,set_len
    CALL COMPARE_VECS(vec,set(:,i),tempres)
    res = res .and. tempres
  ENDDO
END SUBROUTINE COMPARE_PT_SET


!********************************************************************************
!  Compares two set of sorted positions (SORT_POS()) and check if the 
!                 symmetry operations is present or not
!********************************************************************************
SUBROUTINE COMPARE_POS(POS1,POS2,NAT,res)
  INTEGER,INTENT(IN):: NAT
  REAL(wp),DIMENSION(3,NAT),INTENT(IN):: POS1 , POS2
  LOGICAL,INTENT(OUT):: res

  !local variables
  INTEGER:: i
  LOGICAL:: locres
  DO i=1,NAT
    CALL COMPARE_VECS(POS1(:,i),POS2(:,i),locres)
    IF(.not. locres) THEN
      res = .FALSE.
      RETURN
    ENDIF
  ENDDO
 
  res = .TRUE.
END SUBROUTINE COMPARE_POS


!********************************************************************************
!  Sorts position of SEA set to be able to compare action of a Symm. operation
!********************************************************************************
SUBROUTINE SORT_POS(POS,NAT)
  INTEGER,INTENT(IN):: NAT
  REAL(wp),DIMENSION(3,NAT),INTENT(INOUT):: POS
  
  !local variables
  INTEGER:: i, j, k, l 
  INTEGER:: vec_len
  INTEGER:: ng=0 , ng_old=0
  INTEGER,DIMENSION(2,NAT):: bounds 
  REAL(wp),DIMENSION(3):: swap
  REAL(wp):: swapx
  REAL(wp),DIMENSION(NAT):: vec

  bounds = 0 

  ng_old = 0
  ng = 1 
  bounds(1,1) = 1
  bounds(2,1) = NAT

  IF(NAT==1) THEN
    RETURN
  ENDIF

  DO i=1,3  !iterate along coordinates
    DO j=ng_old+1,ng
      vec_len = bounds(2,j)-bounds(1,j) + 1
      vec(1:vec_len) = POS(i,bounds(1,j):bounds(2,j))
      !SORT particular group of positions with same coordinates i
      IF(vec_len>1) THEN
        DO k=1,vec_len-1
          DO l=k+1,vec_len
             IF(vec(l)<vec(k)) THEN
               swap = POS(:,bounds(1,j)+l-1)
               POS(:,bounds(1,j)+l-1) = POS(:,bounds(1,j)+k-1)
               POS(:,bounds(1,j)+k-1) = swap

               swapx = vec(l)
               vec(l) = vec(k)
               vec(k) = swapx  
             ENDIF
          ENDDO
        ENDDO
      ENDIF

    
    ENDDO
     
    !Scan same coordinates for different atoms in a particular group
    ng_old = ng         !and update ng in a loop below 
    j = 1
    DO WHILE(j<NAT)
      k = j + 1
      DO WHILE(ABS(POS(i,k)-POS(i,j))<dvec .and. k<NAT)
          k = k + 1
      ENDDO

      IF((k-j)/=1) THEN
        ng = ng + 1
        bounds(1,ng) = j
        bounds(2,ng) = k-1
      ENDIF
    
      j = k 
    ENDDO

    IF(ng_old==ng) THEN
      EXIT
    ENDIF

  ENDDO

END SUBROUTINE SORT_POS


!********************************************************************************
!  Symmetry action on set of POSITIONS
!       INPUT:
!         SYMM_MAT(3,3) - SYMMETRY MATRIX
!         POS(N,3)      - ORIGINAL POSITIONS OF ATOMS IN MOL/SEA Set
!       OUTPUT:
!         res           - TRUE, IF new positions are only interchanged
!                         FALSE, otherwise
!********************************************************************************
SUBROUTINE SYMM_ACT(MAT,POS,N,res)
  REAL(wp),DIMENSION(3,3),INTENT(IN):: MAT
  INTEGER,INTENT(IN):: N
  REAL(wp),DIMENSION(3,N),INTENT(IN):: POS
  LOGICAL,INTENT(OUT):: res

  !local varaibles
  INTEGER:: i 
  REAL(wp),DIMENSION(3,N):: POS2
  
  DO i=1,N
    POS2(:,i) = MATMUL(MAT,POS(:,i))
  ENDDO  

  CALL SORT_POS(POS2,N)

  CALL COMPARE_POS(POS,POS2,N,res)
END SUBROUTINE SYMM_ACT


!********************************************************************************
!  Determine Rotor type of entire molecule
!       INPUT:
!               NAT - number of atoms
!               IMOM - moment of inertia
!       OUTPUT:
!               ROTTIP                  Possible Molecular Point Groups :
!                 1 - Linear            Cbv , Dbh
!                 2 - Spherical         Td , Th , Oh , Ih , K
!                 3 - Symmetric         Cn , Cnh , Cnv , Dn , Dnh , Dnd , Sn
!                 4 - Asymmetric        C1 , Cs , Ci , Cn , Cnh , Cnv , Dn , Dnh
!                                       Dnd , Sn
!                 0 - ERROR
!                -1 - Atom
!********************************************************************************
SUBROUTINE DET_ROT_TYPE(NAT,IMOM,ROTTIP,PLAN)
  INTEGER,INTENT(IN):: NAT
  REAL(wp),DIMENSION(3),INTENT(IN):: IMOM
  INTEGER,INTENT(OUT):: ROTTIP
  LOGICAL,INTENT(OUT):: PLAN

  !local variables
  IF(NAT==1) THEN
    ROTTIP = -1
    RETURN
  ENDIF

  IF(ABS(IMOM(1))<delem .and. ABS(IMOM(2)-IMOM(3))<delem) THEN
    !Linear Rotor Type
    ROTTIP = 1
  ELSEIF(ABS(IMOM(1)-IMOM(2))<delem .and. ABS(IMOM(2)-IMOM(3))<delem) THEN
    !Spherical Rotor Type
    ROTTIP = 2
  ELSEIF(ABS(IMOM(1)-IMOM(2))<delem .and. ABS(IMOM(2)-IMOM(3))>delem) THEN
    !Symmetric Rotor Type
    ROTTIP = 3
  ELSEIF(ABS(IMOM(1)-IMOM(2))>delem .and. ABS(IMOM(2)-IMOM(3))<delem) THEN
    !Symmetric Rotor Type
    ROTTIP = 3
  ELSEIF(ABS(IMOM(1)-IMOM(2))>delem .and. ABS(IMOM(2)-IMOM(3))>delem) THEN
    !Astmmetric Rotor Type
    ROTTIP = 4
  ELSE
    WRITE(*,*) "Error can be expected - the Rotor type of molecule not recognised "
    ROTTIP = 0
  ENDIF

  IF(ABS(IMOM(1)+IMOM(2)-IMOM(3))<delem) THEN
    PLAN = .TRUE.
  ELSE
    PLAN = .FALSE.
  ENDIF
END SUBROUTINE DET_ROT_TYPE


!********************************************************************************
!  Generates Identity Matrix of dimension N
!********************************************************************************
FUNCTION IDENT(N) RESULT(E)
  INTEGER:: N
  REAL(wp),DIMENSION(N,N):: E

  !local variables
  INTEGER:: i

  E = 0.0_wp

  DO i=1,N
    E(i,i) = 1.0_wp
  ENDDO
END FUNCTION IDENT


!********************************************************************************
!  Calculate a cross product between two vectors 
!********************************************************************************
FUNCTION CROSSPROD(v1,v2) RESULT(v)
  REAL(wp),DIMENSION(3):: v1, v2, v

  !local variables

  v(1) = v1(2)*v2(3) - v1(3)*v2(2)
  v(2) = v1(3)*v2(1) - v1(1)*v2(3)
  v(3) = v1(1)*v2(2) - v1(2)*v2(1)
END FUNCTION CROSSPROD


!********************************************************************************
!  Generates a tensor vector product of two vectors of length N
!              Output matrix M is always symmetric
!********************************************************************************
FUNCTION TENPROD(N,v1,v2) RESULT(MAT)
  INTEGER:: N
  REAL(wp),DIMENSION(N):: v1, v2
  REAL(wp),DIMENSION(N,N):: MAT

  INTEGER:: i,j

  MAT = 0.0_wp

  !Calculate off-diagonal elements
  DO i=1,N 
    DO j=i+1,N
      MAT(j,i) = v1(j) * v2(i)
    ENDDO
  ENDDO

  MAT = MAT + TRANSPOSE(MAT)

  DO i=1,N
    MAT(i,i) = v1(i)*v2(i)
  ENDDO
END FUNCTION TENPROD


!********************************************************************************
!  Cenerates Cross product of matrix for vector n with |n|=1
!********************************************************************************
FUNCTION CROSSMATPROD(nvec) RESULT(MAT)
  REAL(wp),DIMENSION(3):: nvec
  REAL(wp),DIMENSION(3,3):: MAT

  MAT(1,1) = 0.0_wp
  MAT(2,1) = nvec(3)
  MAT(3,1) = -nvec(2)
  MAT(1,2) = -nvec(3)
  MAT(2,2) = 0.0_wp
  MAT(3,2) = nvec(1)
  MAT(1,3) = nvec(2)
  MAT(2,3) = -nvec(1)
  MAT(3,3) = 0.0_wp
END FUNCTION CROSSMATPROD


!********************************************************************************
!  Generates Rotation Matrix from SO3 Group Cn = C(n,nvec
!      Rotation about axis nvec by angle phi = 2pi/n
!********************************************************************************
FUNCTION ROT_MAT(n,nvec) RESULT(Cn)
  INTEGER:: n
  REAL(wp),DIMENSION(3):: nvec
  REAL(wp),DIMENSION(3,3):: Cn

  !local variables
  REAL(wp):: phi

  nvec = nvec/NORM2(nvec)

  !Rotation angle
  phi = 2.0_wp*pi/REAL(n)

  Cn = COS(phi)*IDENT(3) + (1.0_wp-COS(phi))*TENPROD(3,nvec,nvec) + SIN(phi)*CROSSMATPROD(nvec)
END FUNCTION ROT_MAT


!********************************************************************************
!  Generates Reflectin Matrix with normal vector nvec, |nvec|=1
!********************************************************************************
FUNCTION REF_MAT(nvec) RESULT(Sigma)
  REAL(wp),DIMENSION(3):: nvec
  REAL(wp),DIMENSION(3,3):: Sigma

  nvec = nvec/NORM2(nvec)
  Sigma = IDENT(3) - 2.0_wp*TENPROD(3,nvec,nvec)
END FUNCTION REF_MAT


!********************************************************************************
!  Generates Rotation-Reflection Matrix about axis nvec by angle phi = 2pi/n
!********************************************************************************
FUNCTION REF_ROT_MAT(n,nvec) RESULT(Sn)
  INTEGER:: n
  REAL(wp),DIMENSION(3):: nvec
  REAL(wp),DIMENSION(3,3):: Sn  

  Sn = MATMUL(REF_MAT(nvec),ROT_MAT(n,nvec))
END FUNCTION REF_ROT_MAT

END MODULE SYMM
