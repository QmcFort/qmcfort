! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

MODULE LAPACK

#include "preproc.inc"
use constants
use string, only: lowercase
use standalone

implicit none

INTERFACE GEMM                  !Matrix matrix multiplication
  MODULE PROCEDURE GEMM_r_sp
  MODULE PROCEDURE GEMM_r_dp
  MODULE PROCEDURE GEMM_c_sp
  MODULE PROCEDURE GEMM_c_dp
  MODULE PROCEDURE GEMM_rc_sp
  MODULE PROCEDURE GEMM_cr_sp
  MODULE PROCEDURE GEMM_rc_dp
  MODULE PROCEDURE GEMM_cr_dp
END INTERFACE GEMM

interface trmm                  !Matrix L(U) Matrix multiplication
  module procedure trmm_r_sp, trmm_r_dp, trmm_c_sp, trmm_c_dp
end interface trmm

INTERFACE GEMV                  !Matrix vector multiplication
  MODULE PROCEDURE GEMV_r_sp
  MODULE PROCEDURE GEMV_r_dp
  MODULE PROCEDURE GEMV_c_sp
  MODULE PROCEDURE GEMV_c_dp
  module procedure gemv_rc_sp, gemv_rc_dp, gemv_cr_sp, gemv_cr_dp
END INTERFACE GEMV

INTERFACE SYEV                  !Symmetric eigenvalue decomposition
  MODULE PROCEDURE SYEV_r_sp
  MODULE PROCEDURE SYEV_r_dp
  MODULE PROCEDURE SYEV_c_sp
  MODULE PROCEDURE SYEV_c_dp
END INTERFACE SYEV

interface gsyev                 !Generalized symmetric eigenvalue decomposition
  module procedure gsyev_r_sp, gsyev_r_dp, gsyev_c_sp, gsyev_c_dp
end interface gsyev

interface geev 
  module procedure geev_r_sp, geev_r_dp, geev_c_sp, geev_c_dp
end interface geev

interface svd                   !Singular value decomposition
  module procedure svd_all_r_sp, svd_all_r_dp, svd_all_c_sp, svd_all_c_dp
  module procedure svd_none_r_sp, svd_none_r_dp, svd_none_c_sp, svd_none_c_dp
end interface

INTERFACE GETRF                 !PLU decomposition
  MODULE PROCEDURE GETRF_r_sp
  MODULE PROCEDURE GETRF_r_dp
  MODULE PROCEDURE GETRF_c_sp
  MODULE PROCEDURE GETRF_c_dp
END INTERFACE GETRF

interface getqr                 !qr decomposition
  module procedure getqr_r_sp, getqr_r_dp, getqr_c_sp, getqr_c_dp
end interface getqr

interface getqr_det             !qr decomposition with determinant
  module procedure getqr_det_r_sp, getqr_det_r_dp, getqr_det_c_sp, getqr_det_c_dp
end interface getqr_det

interface potrf                 !cholesky decomposition
  module procedure potrf_r_sp, potrf_r_dp, potrf_c_sp, potrf_c_dp
end interface potrf

INTERFACE INVERSE                   !Matrix inversion; done with help of LU dec
  MODULE PROCEDURE INVERSE_r_sp
  MODULE PROCEDURE INVERSE_r_dp
  MODULE PROCEDURE INVERSE_c_sp
  MODULE PROCEDURE INVERSE_c_dp
END INTERFACE INVERSE

interface matrix_power
  module procedure matrix_power_c_dp
end interface matrix_power

interface pseudoinverse
  module procedure pseudoinverse_r_dp, pseudoinverse_c_dp
end interface pseudoinverse

interface pseudoinverse_tikhonov
  module procedure pseudoinverse_tikhonov_c_dp
end interface pseudoinverse_tikhonov

INTERFACE INVERSE_TRI                !Inversion of the triangular matrix
  MODULE PROCEDURE INVERSE_TRI_r_sp
  MODULE PROCEDURE INVERSE_TRI_r_dp
  MODULE PROCEDURE INVERSE_TRI_c_sp
  MODULE PROCEDURE INVERSE_TRI_c_dp
END INTERFACE INVERSE_TRI

INTERFACE DETERMINANT           !Computation of the determinant via PLU factorization
  MODULE PROCEDURE DETERMINANT_r_sp
  MODULE PROCEDURE DETERMINANT_r_dp
  MODULE PROCEDURE DETERMINANT_c_sp
  MODULE PROCEDURE DETERMINANT_c_dp
END INTERFACE DETERMINANT

INTERFACE SOLVE_LS
  MODULE PROCEDURE SOLVE_LS_r_sp
  MODULE PROCEDURE SOLVE_LS_r_dp
  MODULE PROCEDURE SOLVE_LS_c_sp
  MODULE PROCEDURE SOLVE_LS_c_dp
END INTERFACE SOLVE_LS

!********************************************************************************
!       Module for usage of LAPACK routinen
!********************************************************************************

!********************************************************************************
!       Routines:
!********************************************************************************
!       MAT_MAT_PROD            (N,A,B,C,mode_A,mode_B)
!       GEN_MAT_MAT_PROD        (M,N,P,A,B,C)
!       GEN_MAT_VEC_PROD        (N,M,A,X,Y)
!       SYM_MAT_VEC_PROD        (N,A,X,Y)
!       SYM_EIG                 (N,A,EW)
!       SYM_EIGV                (N,A,EW,EV)
!       GEN_EIGV                (N,A,EW,EV)
!       SIGV                    (M,N,A,U,S,VT)
!********************************************************************************
CONTAINS

!******************************************************************************** 
!               GEMM Routine - Matrix multiplication 
!******************************************************************************** 
SUBROUTINE GEMM_r_sp(TA,TB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
  CHARACTER:: TA, TB
  INTEGER:: M, N, K, LDA, LDB, LDC 
  REAL(sp):: ALPHA, BETA
  REAL(sp):: A(:,:) , B(:,:), C(:,:)

  CALL SGEMM(TA,TB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
END SUBROUTINE GEMM_r_sp

SUBROUTINE GEMM_r_dp(TA,TB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
  CHARACTER:: TA, TB
  INTEGER:: M, N, K, LDA, LDB, LDC 
  REAL(dp):: ALPHA, BETA
  REAL(dp):: A(:,:) , B(:,:), C(:,:)

  CALL DGEMM(TA,TB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
END SUBROUTINE GEMM_r_dp

SUBROUTINE GEMM_c_sp(TA,TB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
  CHARACTER:: TA, TB
  INTEGER:: M, N, K, LDA, LDB, LDC 
  COMPLEX(sp):: ALPHA, BETA
  COMPLEX(sp):: A(:,:) , B(:,:), C(:,:)

  CALL CGEMM(TA,TB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
END SUBROUTINE GEMM_c_sp

SUBROUTINE GEMM_c_dp(TA,TB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
  CHARACTER:: TA, TB
  INTEGER:: M, N, K, LDA, LDB, LDC 
  COMPLEX(dp):: ALPHA, BETA
  COMPLEX(dp):: A(:,:) , B(:,:), C(:,:)

  CALL ZGEMM(TA,TB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
END SUBROUTINE GEMM_c_dp

subroutine gemm_rc_sp(ta,tb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
  character, intent(in)      :: ta, tb
  integer, intent(in)        :: m, n, k, lda, ldb, ldc
  complex(sp), intent(in)    :: alpha, beta
  real(sp), intent(in)       :: a(:,:)
  complex(sp), intent(in)    :: b(:,:)
  complex(sp), intent(inout) :: c(:,:)
  !local variables
  character                     :: tb_
  real(sp), parameter        :: tol=1.0e-06_sp
  real(sp)                   :: fact, alpha_r, beta_r
  real(sp), allocatable      :: b_(:,:)
  real(sp), allocatable      :: c_(:,:)
  
  if (lowercase(tb) =="c") then
    tb_ = 't'
    fact = -1.0_sp
  else
    tb_ = tb
    fact = 1.0_sp
  end if
  
  allocate(b_(k,2*n), c_(m,2*n))

  if (tb_ == "t") then
    b_(:,1:n) = transpose(real(b, kind=sp))
    b_(:,n+1:2*n) = transpose(imag(b))
  else 
    b_(:,1:n) = real(b, kind=sp)
    b_(:,n+1:2*n) = imag(b)
  end if

  call sgemm(ta,"n",m,2*n,k,1.0_sp,a,lda,b_,k,0.0_sp,c_,ldc)
  
  if (abs(beta) < tol) then
    c = alpha*cmplx(c_(:,1:n),fact*c_(:,n+1:2*n),kind=sp)
  else
    c = alpha*cmplx(c_(:,1:n),fact*c_(:,n+1:2*n),kind=sp) + beta*c
  end if
end subroutine gemm_rc_sp

subroutine gemm_rc_dp(ta,tb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
  character, intent(in)      :: ta, tb
  integer, intent(in)        :: m, n, k, lda, ldb, ldc
  complex(dp), intent(in)    :: alpha, beta
  real(dp), intent(in)       :: a(:,:)
  complex(dp), intent(in)    :: b(:,:)
  complex(dp), intent(inout) :: c(:,:)
  !local variables
  character                     :: tb_
  real(dp), parameter        :: tol=1.0e-06_dp
  real(dp)                   :: fact, alpha_r, beta_r
  real(dp), allocatable      :: b_(:,:)
  real(dp), allocatable      :: c_(:,:)
  
  if (lowercase(tb) =="c") then
    tb_ = 't'
    fact = -1.0_dp
  else
    tb_ = tb
    fact = 1.0_dp
  end if
  
  allocate(b_(k,2*n), c_(m,2*n))

  if (tb_ == "t") then
    b_(:,1:n) = transpose(real(b, kind=dp))
    b_(:,n+1:2*n) = transpose(imag(b))
  else 
    b_(:,1:n) = real(b, kind=dp)
    b_(:,n+1:2*n) = imag(b)
  end if

  call dgemm(ta,"n",m,2*n,k,1.0_dp,a,lda,b_,k,0.0_dp,c_,ldc)
  
  if (abs(beta) < tol) then
    c = alpha*cmplx(c_(:,1:n),fact*c_(:,n+1:2*n),kind=dp)
  else
    c = alpha*cmplx(c_(:,1:n),fact*c_(:,n+1:2*n),kind=dp) + beta*c
  end if
end subroutine gemm_rc_dp

subroutine gemm_cr_sp(ta,tb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
  character, intent(in)         :: ta, tb
  integer, intent(in)           :: m, n, k, lda, ldb, ldc
  complex(sp), intent(in)    :: alpha, beta
  complex(sp), intent(in)    :: a(:,:)
  real(sp), intent(in)    :: b(:,:)
  complex(sp), intent(inout) :: c(:,:)
  !local variables
  character                     :: ta_
  real(sp), parameter        :: tol=1.0e-06_sp
  real(sp)                   :: fact, alpha_r, beta_r
  real(sp), allocatable      :: a_(:,:)
  real(sp), allocatable      :: c_(:,:)
  
  if (lowercase(ta) =="c") then
    ta_ = 't'
    fact = -1.0_sp
  else
    ta_ = ta
    fact = 1.0_sp
  end if
  
  allocate(a_(2*m,k), c_(2*m,n))

  if (ta_ == "t") then
    a_(1:m,:) = transpose(real(a, kind=sp))
    a_(m+1:2*m,:) = transpose(imag(a))
  else 
    a_(1:m,:) = real(a, kind=sp)
    a_(m+1:2*m,:) = imag(a)
  end if

  call sgemm("n",tb,2*m,n,k,1.0_sp,a_,2*m,b,ldb,0.0_sp,c_,2*ldc)
  
  if (abs(beta) < tol) then
    c = alpha*cmplx(c_(1:m,:),fact*c_(m+1:2*m,:),kind=sp)
  else
    c = alpha*cmplx(c_(1:m,:),fact*c_(m+1:2*m,:),kind=sp) + beta*c
  end if
end subroutine gemm_cr_sp

subroutine gemm_cr_dp(ta,tb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
  character, intent(in)      :: ta, tb
  integer, intent(in)        :: m, n, k, lda, ldb, ldc
  complex(dp), intent(in)    :: alpha, beta
  complex(dp), intent(in)    :: a(:,:)
  real(dp), intent(in)       :: b(:,:)
  complex(dp), intent(inout) :: c(:,:)
  !local variables
  character                  :: ta_
  real(dp), parameter        :: tol=1.0e-06_dp
  real(dp)                   :: fact, alpha_r, beta_r
  real(dp), allocatable      :: a_(:,:)
  real(dp), allocatable      :: c_(:,:)
  
  if (lowercase(ta) =="c") then
    ta_ = 't'
    fact = -1.0_dp
  else
    ta_ = ta
    fact = 1.0_dp
  end if
  
  allocate(a_(2*m,k), c_(2*m,n))

  if (ta_ == "t") then
    a_(1:m,:) = transpose(real(a, kind=dp))
    a_(m+1:2*m,:) = transpose(imag(a))
  else 
    a_(1:m,:) = real(a, kind=dp)
    a_(m+1:2*m,:) = imag(a)
  end if

  call dgemm("n",tb,2*m,n,k,1.0_dp,a_,2*m,b,ldb,0.0_dp,c_,2*ldc)
  
  if (abs(beta) < tol) then
    c = alpha*cmplx(c_(1:m,:),fact*c_(m+1:2*m,:),kind=dp)
  else
    c = alpha*cmplx(c_(1:m,:),fact*c_(m+1:2*m,:),kind=dp) + beta*c
  end if
end subroutine gemm_cr_dp

!******************************************************************************** 
!
!       trmm Routine - Multiply matrix by lower (upper) triangular matrix
!
!******************************************************************************** 
subroutine trmm_r_sp(side, uplo, ta, diag, m, n, alpha, a, lda, b, ldb)
  character, intent(in)      :: side, uplo, ta, diag
  integer, intent(in)        :: m, n, lda, ldb
  real(sp), intent(in)    :: alpha
  real(sp), intent(in)    :: a(:,:)
  real(sp), intent(inout) :: b(:,:)

  call strmm(side, uplo, ta, diag, m, n, alpha, a, lda, b, ldb)
end subroutine trmm_r_sp

subroutine trmm_r_dp(side, uplo, ta, diag, m, n, alpha, a, lda, b, ldb)
  character, intent(in)      :: side, uplo, ta, diag
  integer, intent(in)        :: m, n, lda, ldb
  real(dp), intent(in)    :: alpha
  real(dp), intent(in)    :: a(:,:)
  real(dp), intent(inout) :: b(:,:)

  call dtrmm(side, uplo, ta, diag, m, n, alpha, a, lda, b, ldb)
end subroutine trmm_r_dp

subroutine trmm_c_sp(side, uplo, ta, diag, m, n, alpha, a, lda, b, ldb)
  character, intent(in)         :: side, uplo, ta, diag
  integer, intent(in)           :: m, n, lda, ldb
  complex(sp), intent(in)    :: alpha
  complex(sp), intent(in)    :: a(:,:)
  complex(sp), intent(inout) :: b(:,:)

  call ctrmm(side, uplo, ta, diag, m, n, alpha, a, lda, b, ldb)
end subroutine trmm_c_sp

subroutine trmm_c_dp(side, uplo, ta, diag, m, n, alpha, a, lda, b, ldb)
  character, intent(in)         :: side, uplo, ta, diag
  integer, intent(in)           :: m, n, lda, ldb
  complex(dp), intent(in)    :: alpha
  complex(dp), intent(in)    :: a(:,:)
  complex(dp), intent(inout) :: b(:,:)

  call ztrmm(side, uplo, ta, diag, m, n, alpha, a, lda, b, ldb)
end subroutine trmm_c_dp


!******************************************************************************** 
!       GEMV Routine - Matrix vector multiplication
!******************************************************************************** 
SUBROUTINE GEMV_r_sp(T,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
  CHARACTER:: T
  INTEGER:: M, N, LDA, INCX, INCY
  REAL(sp):: ALPHA, BETA
  REAL(sp):: A(:,:)
  REAL(sp):: X(:), Y(:)

  CALL SGEMV(T,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
END SUBROUTINE GEMV_r_sp

SUBROUTINE GEMV_r_dp(T,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
  CHARACTER:: T
  INTEGER:: M, N, LDA, INCX, INCY
  REAL(dp):: ALPHA, BETA
  REAL(dp):: A(:,:)
  REAL(dp):: X(:), Y(:)

  CALL DGEMV(T,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
END SUBROUTINE GEMV_r_dp

SUBROUTINE GEMV_c_sp(T,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
  CHARACTER:: T
  INTEGER:: M, N, LDA, INCX, INCY
  COMPLEX(sp):: ALPHA, BETA
  COMPLEX(sp):: A(:,:)
  COMPLEX(sp):: X(:), Y(:)

  CALL CGEMV(T,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
END SUBROUTINE GEMV_c_sp

SUBROUTINE GEMV_c_dp(T,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
  CHARACTER:: T
  INTEGER:: M, N, LDA, INCX, INCY
  COMPLEX(dp):: ALPHA, BETA
  COMPLEX(dp):: A(:,:)
  COMPLEX(dp):: X(:), Y(:)

  CALL ZGEMV(T,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
END SUBROUTINE GEMV_c_dp

subroutine gemv_rc_sp(t, m, n, alpha, a, lda, x, incx, beta, y, incy)
  character                :: t
  integer                  :: m, n, lda, incx, incy
  complex(sp)           :: alpha, beta
  real(sp)              :: a(:,:)
  complex(sp)           :: x(:), y(:)
  !local variables
  real(sp), allocatable :: yr(:), yi(:)

  allocate(yr(size(y)), yi(size(y)))
  yr = 0.0_sp
  yi = 0.0_sp

  call sgemv(t, m, n, 1.0_sp, a, lda, real(x, sp), incx, 0.0_sp, yr, incy)
  call sgemv(t, m, n, 1.0_sp, a, lda, imag(x), incx, 0.0_sp, yi, incy)

  y = alpha*cmplx(yr, yi, kind=sp) + beta*y
end subroutine gemv_rc_sp

subroutine gemv_rc_dp(t, m, n, alpha, a, lda, x, incx, beta, y, incy)
  character                :: t
  integer                  :: m, n, lda, incx, incy
  complex(dp)           :: alpha, beta
  real(dp)              :: a(:,:)
  complex(dp)           :: x(:), y(:)
  !local variables
  real(dp), allocatable :: yr(:), yi(:)

  allocate(yr(size(y)), yi(size(y)))
  yr = 0.0_dp
  yi = 0.0_dp

  call dgemv(t, m, n, 1.0_dp, a, lda, real(x, dp), incx, 0.0_dp, yr, incy)
  call dgemv(t, m, n, 1.0_dp, a, lda, imag(x), incx, 0.0_dp, yi, incy)

  y = alpha*cmplx(yr, yi, kind=dp) + beta*y
end subroutine gemv_rc_dp

subroutine gemv_cr_sp(t, m, n, alpha, a, lda, x, incx, beta, y, incy)
  character      :: t
  integer        :: m, n, lda, incx, incy
  complex(sp) :: alpha, beta
  complex(sp) :: a(:,:), y(:)
  real(sp)    :: x(:)
  !local variables
  character                :: t_
  real(sp)              :: fact
  real(sp), allocatable :: yr(:), yi(:)

  if (t=='c' .or. t=='C') then
    t_ = 'T'
    fact = -1.0_sp
  else
    t_ = t
    fact = 1.0_sp
  end if

  allocate(yr(size(y)), yi(size(y)))

  call sgemv(t_, m, n, 1.0_sp, real(a, sp), lda, x, incx, 0.0_sp, yr, incy)
  call sgemv(t_, m, n, 1.0_sp, imag(a), lda, x, incx, 0.0_sp, yi, incy)

  y = alpha*cmplx(yr, fact*yi, kind=sp) + beta*y
end subroutine gemv_cr_sp

subroutine gemv_cr_dp(t, m, n, alpha, a, lda, x, incx, beta, y, incy)
  character      :: t
  integer        :: m, n, lda, incx, incy
  complex(dp) :: alpha, beta
  complex(dp) :: a(:,:), y(:)
  real(dp)    :: x(:)
  !local variables
  character                :: t_
  real(sp)              :: fact
  real(dp), allocatable :: yr(:), yi(:)

  if (t=='c' .or. t=='C') then
    t_ = 'T'
    fact = -1.0_dp
  else
    t_ = t
    fact = 1.0_dp
  end if

  allocate(yr(size(y)), yi(size(y)))

  call dgemv(t_, m, n, 1.0_dp, real(a, dp), lda, x, incx, 0.0_dp, yr, incy)
  call dgemv(t_, m, n, 1.0_dp, imag(a), lda, x, incx, 0.0_dp, yi, incy)

  y = alpha*cmplx(yr, fact*yi, kind=dp) + beta*y
end subroutine gemv_cr_dp


!********************************************************************************
!       SYEV Routine - Symmetric eigenvalue problem
!               INPUT:
!                       A       - Matrix that have to be diagonalized
!                       JOBZ    - V if eigenvectors are needed too
!                       UPLO    - U if the matrix A is stored in upper triangular
!                                 part and L for lower triangular part           
!               OUTPUT:
!                       W       - Array of eigenvalues
!                       A       - each column is eigenvector    
!******************************************************************************** 
SUBROUTINE SYEV_r_sp(A,W,JOBZ,UPLO)
  REAL(sp):: A(:,:), W(:)
  CHARACTER,OPTIONAL:: JOBZ, UPLO

  !local variables
  INTEGER:: N, M, LDA, LWORK, INFO
  CHARACTER:: myjob, myup
  REAL(sp),ALLOCATABLE:: WORK(:)
  REAL(sp):: PREWORK(1)

  IF(PRESENT(JOBZ)) THEN
    myjob = JOBZ
  ELSE
    myjob = "V"
  ENDIF

  IF(PRESENT(UPLO)) THEN
    myup = UPLO
  ELSE
    myup = "U"
  ENDIF

  N = SIZE(A(:,1))
  M = SIZE(A(1,:))

  IF(N/=M) THEN
    STOP "Matrix A is not square matrix - SSYEV will fail"
  ENDIF
 
  LDA = N

  !Determine size of the WORK array
  CALL SSYEV(myjob,myup,N,A,LDA,W,PREWORK,-1,INFO)

  LWORK = INT(PREWORK(1))
  ALLOCATE(WORK(LWORK))

  !Do diagonalization
  CALL SSYEV(myjob,myup,N,A,LDA,W,WORK,LWORK,INFO)

  IF(INFO/=0) THEN
    STOP "Problem in SSYEV routine"
  ENDIF

END SUBROUTINE SYEV_r_sp

SUBROUTINE SYEV_r_dp(A,W,JOBZ,UPLO)
  REAL(dp):: A(:,:), W(:)
  CHARACTER,OPTIONAL:: JOBZ, UPLO

  !local variables
  INTEGER:: N, M, LDA, LWORK, INFO
  CHARACTER:: myjob, myup
  REAL(dp),ALLOCATABLE:: WORK(:)
  REAL(dp):: PREWORK(1)

  IF(PRESENT(JOBZ)) THEN
    myjob = JOBZ
  ELSE
    myjob = "V"
  ENDIF

  IF(PRESENT(UPLO)) THEN
    myup = UPLO
  ELSE
    myup = "U"
  ENDIF

  N = SIZE(A(:,1))
  M = SIZE(A(1,:))

  IF(N/=M) THEN
    STOP "Matrix A is not square matrix - DSYEV will fail"
  ENDIF
 
  LDA = N

  !Determine size of the WORK array
  CALL DSYEV(myjob,myup,N,A,LDA,W,PREWORK,-1,INFO)

  LWORK = INT(PREWORK(1))
  ALLOCATE(WORK(LWORK))

  !Do diagonalization
  CALL DSYEV(myjob,myup,N,A,LDA,W,WORK,LWORK,INFO)

  IF(INFO/=0) THEN
    STOP "Problem in DSYEV routine"
  ENDIF

END SUBROUTINE SYEV_r_dp

SUBROUTINE SYEV_c_sp(A,W,JOBZ,UPLO)
  COMPLEX(sp):: A(:,:)
  REAL(sp):: W(:)
  CHARACTER,OPTIONAL:: JOBZ, UPLO

  !local variables
  INTEGER:: N, M, LDA, LWORK, LWORK2, INFO
  CHARACTER:: myjob, myup
  COMPLEX(sp),ALLOCATABLE:: WORK(:)
  COMPLEX(sp):: PREWORK(1)
  REAL(sp),ALLOCATABLE:: RWORK(:)

  IF(PRESENT(JOBZ)) THEN
    myjob = JOBZ
  ELSE
    myjob = "V"
  ENDIF

  IF(PRESENT(UPLO)) THEN
    myup = UPLO
  ELSE
    myup = "U"
  ENDIF

  N = SIZE(A(:,1))
  M = SIZE(A(1,:))

  IF(N/=M) THEN
    STOP "Matrix A is not square matrix - CHEEV will fail"
  ENDIF
 
  LDA = N
 
  !Determine size of RWORK array
  LWORK2 = MAX(1,3*N-2)
  ALLOCATE(RWORK(LWORK2))

  !Determine size of the WORK array
  CALL CHEEV(myjob,myup,N,A,LDA,W,PREWORK,-1,RWORK,INFO)

  LWORK = INT(PREWORK(1))
  ALLOCATE(WORK(LWORK))

  !Do diagonalization
  CALL CHEEV(myjob,myup,N,A,LDA,W,WORK,LWORK,RWORK,INFO)

  IF(INFO/=0) THEN
    STOP "Problem in CHEEV routine"
  ENDIF

END SUBROUTINE SYEV_c_sp

SUBROUTINE SYEV_c_dp(A,W,JOBZ,UPLO)
  COMPLEX(dp):: A(:,:)
  REAL(dp):: W(:)
  CHARACTER,OPTIONAL:: JOBZ, UPLO

  !local variables
  INTEGER:: N, M, LDA, LWORK, LWORK2, INFO
  CHARACTER:: myjob, myup
  COMPLEX(dp),ALLOCATABLE:: WORK(:)
  COMPLEX(dp):: PREWORK(1)
  REAL(dp),ALLOCATABLE:: RWORK(:)

  IF(PRESENT(JOBZ)) THEN
    myjob = JOBZ
  ELSE
    myjob = "V"
  ENDIF

  IF(PRESENT(UPLO)) THEN
    myup = UPLO
  ELSE
    myup = "U"
  ENDIF

  N = SIZE(A(:,1))
  M = SIZE(A(1,:))

  IF(N/=M) THEN
    STOP "Matrix A is not square matrix - ZHEEV will fail"
  ENDIF
 
  LDA = N
 
  !Determine size of RWORK array
  LWORK2 = MAX(1,3*N-2)
  ALLOCATE(RWORK(LWORK2))

  !Determine size of the WORK array
  CALL ZHEEV(myjob,myup,N,A,LDA,W,PREWORK,-1,RWORK,INFO)

  LWORK = INT(PREWORK(1))
  ALLOCATE(WORK(LWORK))

  !Do diagonalization
  CALL ZHEEV(myjob,myup,N,A,LDA,W,WORK,LWORK,RWORK,INFO)

  IF(INFO/=0) THEN
    STOP "Problem in ZHEEV routine"
  ENDIF

END SUBROUTINE SYEV_c_dp


!********************************************************************************
!
! GSYEV Routine - Generalized symmetric eigenvalue decomposition
!
!            A x = w B x
!
!    input:
!            A       - Matrix that have to be diagonalized
!            jobz    - V if eigenvectors are needed too
!            uplo    - U if the matrix A is stored in upper triangular
!                      part and L for lower triangular part           
!    output:
!            w       - Array of eigenvalues
!            A       - each column is eigenvector    
!
!******************************************************************************** 
subroutine gsyev_r_sp(A, B, w, jobz, uplo)
  real(sp)              :: A(:,:)
  real(sp)              :: B(:,:)
  real(sp)              :: w(:)
  character, optional   :: jobz
  character, optional   :: uplo
  !local
  integer               :: n, m, lwork, info
  character             :: jobz_, uplo_
  real(sp)              :: prework(1)
  real(sp), allocatable :: work(:)

  if (present(jobz)) then
    jobz_ = jobz
  else
    jobz_ = "V"
  end if

  if (present(uplo)) then
    uplo_ = uplo
  else
    uplo_ = "U"
  end if

  n = size(A, 1)
  m = size(A, 2)

  if (n /= m) stop "Matrix A is not square matrix - gsyev_r_sp will fail"

  !Determine size of the WORK array
  call ssygv(1, jobz_, uplo_, n, A, n, B, n, w, prework, -1, info) 

  lwork = int(prework(1))
  allocate(work(lwork))

  !actual lapack call
  call ssygv(1, jobz_, uplo_, n, A, n, B, n, w, work, lwork, info) 

  if (info /= 0) stop "Problem in gsyev_r_sp routine"
end subroutine gsyev_r_sp

subroutine gsyev_r_dp(A, B, w, jobz, uplo)
  real(dp)              :: A(:,:)
  real(dp)              :: B(:,:)
  real(dp)              :: w(:)
  character, optional   :: jobz
  character, optional   :: uplo
  !local
  integer               :: n, m, lwork, info
  character             :: jobz_, uplo_
  real(dp)              :: prework(1)
  real(dp), allocatable :: work(:)

  if (present(jobz)) then
    jobz_ = jobz
  else
    jobz_ = "V"
  end if

  if (present(uplo)) then
    uplo_ = uplo
  else
    uplo_ = "U"
  end if

  n = size(A, 1)
  m = size(A, 2)

  if (n /= m) stop "Matrix A is not square matrix - gsyev_r_dp will fail"

  !Determine size of the WORK array
  call dsygv(1, jobz_, uplo_, n, A, n, B, n, w, prework, -1, info) 

  lwork = int(prework(1))
  allocate(work(lwork))

  !actual lapack call
  call dsygv(1, jobz_, uplo_, n, A, n, B, n, w, work, lwork, info) 

  if (info /= 0) stop "Problem in gsyev_r_dp routine"
end subroutine gsyev_r_dp

subroutine gsyev_c_sp(A, B, w, jobz, uplo)
  complex(sp)              :: A(:,:)
  complex(sp)              :: B(:,:)
  real(sp)                 :: w(:)
  character, optional      :: jobz
  character, optional      :: uplo
  !local
  integer                  :: n, m, lwork, info
  character                :: jobz_, uplo_
  complex(sp)              :: prework(1)
  real(sp), allocatable    :: rwork(:)
  complex(sp), allocatable :: work(:)

  if (present(jobz)) then
    jobz_ = jobz
  else
    jobz_ = "V"
  end if

  if (present(uplo)) then
    uplo_ = uplo
  else
    uplo_ = "U"
  end if

  n = size(A, 1)
  m = size(A, 2)

  if (n /= m) stop "Matrix A is not square matrix - gsyev_c_sp will fail"

  allocate(rwork(max(1, 3*n-2)))

  !Determine size of the WORK array
  call chegv(1, jobz_, uplo_, n, A, n, B, n, w, prework, -1, rwork, info) 

  lwork = int(prework(1))
  allocate(work(lwork))

  !actual lapack call
  call chegv(1, jobz_, uplo_, n, A, n, B, n, w, work, lwork, rwork, info) 

  if (info /= 0) stop "Problem in gsyev_c_sp routine"
end subroutine gsyev_c_sp

subroutine gsyev_c_dp(A, B, w, jobz, uplo)
  complex(dp)              :: A(:,:)
  complex(dp)              :: B(:,:)
  real(dp)                 :: w(:)
  character, optional      :: jobz
  character, optional      :: uplo
  !local
  integer                  :: n, m, lwork, info
  character                :: jobz_, uplo_
  complex(dp)              :: prework(1)
  real(dp), allocatable    :: rwork(:)
  complex(dp), allocatable :: work(:)

  if (present(jobz)) then
    jobz_ = jobz
  else
    jobz_ = "V"
  end if

  if (present(uplo)) then
    uplo_ = uplo
  else
    uplo_ = "U"
  end if

  n = size(A, 1)
  m = size(A, 2)

  if (n /= m) stop "Matrix A is not square matrix - gsyev_c_dp will fail"

  allocate(rwork(max(1, 3*n-2)))

  !Determine size of the WORK array
  call zhegv(1, jobz_, uplo_, n, A, n, B, n, w, prework, -1, rwork, info) 

  lwork = int(prework(1))
  allocate(work(lwork))

  !actual lapack call
  call zhegv(1, jobz_, uplo_, n, A, n, B, n, w, work, lwork, rwork, info) 

  if (info /= 0) stop "Problem in gsyev_c_dp routine"
end subroutine gsyev_c_dp


!********************************************************************************
!
!       geev routine - general matrices eigenvalue problem
!
!               INPUT:
!                       A       - Matrix that have to be diagonalized
!                       JOBZ    - V if eigenvectors are needed too
!                       UPLO    - U if the matrix A is stored in upper triangular
!                                 part and L for lower triangular part           
!               OUTPUT:
!                       W       - Array of eigenvalues
!                       A       - each column is eigenvector    
!******************************************************************************** 
subroutine geev_r_sp(A, ev, U, jobz)
  real(sp), intent(in)         :: A(:,:)
  complex(sp), intent(out)     :: ev(:)
  complex(sp), intent(out)     :: U(:,:)
  character, optional, intent(in) :: jobz
  !local variables
  integer                         :: i, n, m, lda, lwork, info
  character                       :: jobz_
  real(sp), allocatable        :: work(:), evr(:), evi(:), A_(:,:), Ul(:,:), Ur(:,:)
  real(sp)                     :: prework(1)
  real(sp), parameter          :: tol = 1.0e-06_sp

  if (present(jobz)) then
    jobz_ = lowercase(jobz)
  else
    jobz_ = "v"
  end if

  n = size(A, 1)
  m = size(A, 2)

  if (n /= m) then
    stop "matrix a is not square matrix - ssyev will fail"
  end if
 
  lda = n

  allocate(evr(n), evi(n), Ul(n,n), Ur(n,n))
  allocate(A_, source=A)

  !determine size of the work array
  call sgeev("n", jobz_, n, A_, lda, evr, evi, Ul, lda, Ur, lda, prework, -1, info)

  lwork = int(prework(1))
  allocate(work(lwork))

  !do diagonalization
  call sgeev("n", jobz_, n, A_, lda, evr, evi, Ul, lda, Ur, lda, work, lwork, info)

  if (info /= 0) then
    stop "problem in sgeev routine"
  end if

  ev = cmplx(evr, evi, kind=sp)
  
  if (jobz_ /= "v") return

  do i = 1, n
    if (abs(evi(i)) < tol) then
      U(:,i) = Ur(:,i)
    else 
      if (abs(evr(i)-evr(i+1))<tol .and. abs(evi(i)-evi(i+1))<tol) then
        U(:,i) = cmplx(Ur(:,i), Ur(:,i+1), kind=sp)
        U(:,i+1) = cmplx(Ur(:,i), -Ur(:,i+1), kind=sp)
      else 
        cycle
      end if
    end if
  end do
end subroutine geev_r_sp

subroutine geev_r_dp(A, ev, U, jobz)
  real(dp), intent(in)         :: A(:,:)
  complex(dp), intent(out)     :: ev(:)
  complex(dp), intent(out)     :: U(:,:)
  character, optional, intent(in) :: jobz
  !local variables
  integer                         :: i, n, m, lda, lwork, info
  character                       :: jobz_
  real(dp), allocatable        :: work(:), evr(:), evi(:), A_(:,:), Ul(:,:), Ur(:,:)
  real(dp)                     :: prework(1)
  real(dp), parameter          :: tol = 1.0e-06_dp

  if (present(jobz)) then
    jobz_ = lowercase(jobz)
  else
    jobz_ = "v"
  end if

  n = size(A, 1)
  m = size(A, 2)

  if (n /= m) then
    stop "matrix a is not square matrix - ssyev will fail"
  end if
 
  lda = n

  allocate(evr(n), evi(n), Ul(n,n), Ur(n,n))
  allocate(A_, source=A)

  !determine size of the work array
  call dgeev("n", jobz_, n, A_, lda, evr, evi, Ul, lda, Ur, lda, prework, -1, info)

  lwork = int(prework(1))
  allocate(work(lwork))

  !do diagonalization
  call dgeev("n", jobz_, n, A_, lda, evr, evi, Ul, lda, Ur, lda, work, lwork, info)

  if (info /= 0) then
    stop "problem in dgeev routine"
  end if

  ev = cmplx(evr, evi, kind=dp)
  
  if (jobz_ /= "v") return

  do i = 1, n
    if (abs(evi(i)) < tol) then
      U(:,i) = Ur(:,i)
    else 
      if (abs(evr(i)-evr(i+1))<tol .and. abs(evi(i)-evi(i+1))<tol) then
        U(:,i) = cmplx(Ur(:,i), Ur(:,i+1), kind=dp)
        U(:,i+1) = cmplx(Ur(:,i), -Ur(:,i+1), kind=dp)
      else 
        cycle
      end if
    end if
  end do
end subroutine geev_r_dp

subroutine geev_c_sp(A, ev, U, jobz)
  complex(sp), intent(in)      :: A(:,:)
  complex(sp), intent(out)     :: ev(:)
  complex(sp), intent(out)     :: U(:,:)
  character, optional, intent(in) :: jobz
  !local variables
  integer                         :: n, m, lda, lwork, info
  character                       :: jobz_
  real(sp), allocatable        :: rwork(:)
  complex(sp), allocatable     :: work(:), A_(:,:), Ul(:,:)
  complex(sp)                  :: prework(1)
  complex(sp), parameter       :: tol = 1.0e-06_sp

  if (present(jobz)) then
    jobz_ = lowercase(jobz)
  else
    jobz_ = "v"
  end if

  n = size(A, 1)
  m = size(A, 2)

  if (n /= m) then
    stop "matrix a is not square matrix - ssyev will fail"
  end if
 
  lda = n

  allocate(rwork(2*n), Ul(n,n))
  allocate(A_, source=A)

  !determine size of the work array
  call cgeev("n", jobz_, n, A_, lda, ev, Ul, lda, U, lda, prework, -1, rwork, info)

  lwork = int(prework(1))
  allocate(work(lwork))

  !do diagonalization
  call cgeev("n", jobz_, n, A_, lda, ev, Ul, lda, U, lda, work, lwork, rwork, info)

  if (info /= 0) then
    stop "problem in cgeev routine"
  end if
end subroutine geev_c_sp

subroutine geev_c_dp(A, ev, U, jobz)
  complex(dp), intent(in)      :: A(:,:)
  complex(dp), intent(out)     :: ev(:)
  complex(dp), intent(out)     :: U(:,:)
  character, optional, intent(in) :: jobz
  !local variables
  integer                         :: n, m, lda, lwork, info
  character                       :: jobz_
  real(dp), allocatable        :: rwork(:)
  complex(dp), allocatable     :: work(:), A_(:,:), Ul(:,:)
  complex(dp)                  :: prework(1)
  complex(dp), parameter       :: tol = 1.0e-06_dp

  if (present(jobz)) then
    jobz_ = lowercase(jobz)
  else
    jobz_ = "v"
  end if

  n = size(A, 1)
  m = size(A, 2)

  if (n /= m) then
    stop "matrix a is not square matrix - ssyev will fail"
  end if
 
  lda = n

  allocate(rwork(2*n), Ul(n,n))
  allocate(A_, source=A)

  !determine size of the work array
  call zgeev("n", jobz_, n, A_, lda, ev, Ul, lda, U, lda, prework, -1, rwork, info)

  lwork = int(prework(1))
  allocate(work(lwork))

  !do diagonalization
  call zgeev("n", jobz_, n, A_, lda, ev, Ul, lda, U, lda, work, lwork, rwork, info)

  if (info /= 0) then
    stop "problem in cgeev routine"
  end if
end subroutine geev_c_dp


!******************************************************************************** 
!
!       Singular Value decomposition
!       
!       Default:
!               both left and right singulat vectos calculated jobu = jobvt = "a" 
!
!               To compute only singular values use jobu = jobvt = "n"
!
!               To overwrite original matrix use jobu = "o" or jobvt = "o"
!
!******************************************************************************** 
subroutine svd_all_r_sp(A, sigma, U, VT)
  real(sp), intent(inout) :: A(:,:)
  real(sp), intent(out)   :: sigma(:)
  real(sp), intent(out)   :: U(:,:)
  real(sp), intent(out)   :: VT(:,:)
  !local variables
  integer                    :: m, n, lwork, info
  real(sp)                :: prework(1)
  real(sp), allocatable   :: work(:)
  !
  m = size(A, 1)
  n = size(A, 2)
  !
  call sgesvd("a", "a", m, n, A, m, sigma, U, m, VT, n, prework, -1, info)
  !
  if (info /= 0) then
    write(*,*) "problem in svd - program stops now"
    call exit
  end if
  !
  lwork = int(prework(1))
  allocate(work(lwork))
  !
  call sgesvd("a", "a", m, n, A, m, sigma, U, m, VT, n, work, lwork, info)
  !
  if (info /= 0) then
    write(*,*) "problem in svd - program stops now"
    call exit
  end if
end subroutine svd_all_r_sp

subroutine svd_all_r_dp(A, sigma, U, VT)
  real(dp), intent(inout) :: A(:,:)
  real(dp), intent(out)   :: sigma(:)
  real(dp), intent(out)   :: U(:,:)
  real(dp), intent(out)   :: VT(:,:)
  !local variables
  integer                    :: m, n, lwork, info
  real(dp)                :: prework(1)
  real(dp), allocatable   :: work(:)
  !
  m = size(A, 1)
  n = size(A, 2)
  !
  call dgesvd("a", "a", m, n, A, m, sigma, U, m, VT, n, prework, -1, info)
  !
  if (info /= 0) then
    write(*,*) "problem in svd - program stops now"
    call exit
  end if
  !
  lwork = int(prework(1))
  allocate(work(lwork))
  !
  call dgesvd("a", "a", m, n, A, m, sigma, U, m, VT, n, work, lwork, info)
  !
  if (info /= 0) then
    write(*,*) "problem in svd - program stops now"
    call exit
  end if
end subroutine svd_all_r_dp

subroutine svd_all_c_sp(A, sigma, U, VT)
  complex(sp), intent(inout) :: A(:,:)
  real(sp), intent(out)      :: sigma(:)
  complex(sp), intent(out)   :: U(:,:)
  complex(sp), intent(out)   :: VT(:,:)
  !local variables
  integer                    :: m, n, lwork, info
  real(sp), allocatable      :: rwork(:)
  complex(sp)                :: prework(1)
  complex(sp), allocatable   :: work(:)
  
  m = size(A, 1)
  n = size(A, 2)
  
  allocate(rwork(5*min(n,m)))

  call cgesvd("a", "a", m, n, A, m, sigma, U, m, VT, n, prework, -1, rwork, info)
  
  if (info /= 0) then
    write(*,*) "problem in cgesvd routine"
    call exit
  end if
  
  lwork = int(prework(1))
  allocate(work(lwork))
  
  call cgesvd("a", "a", m, n, A, m, sigma, U, m, VT, n, work, lwork, rwork, info)
  
  if (info /= 0) then
    write(*,*) "problem in cgesvd routine"
    call exit
  end if
end subroutine svd_all_c_sp

subroutine svd_all_c_dp(A, sigma, U, VT)
  complex(dp), intent(inout) :: A(:,:)
  real(dp), intent(out)      :: sigma(:)
  complex(dp), intent(out)   :: U(:,:)
  complex(dp), intent(out)   :: VT(:,:)
  !local variables
  integer                    :: m, n, lwork, info
  real(dp), allocatable      :: rwork(:)
  complex(dp)                :: prework(1)
  complex(dp), allocatable   :: work(:)
  
  m = size(A, 1)
  n = size(A, 2)
  
  allocate(rwork(5*min(n,m)))

  call zgesvd("a", "a", m, n, A, m, sigma, U, m, VT, n, prework, -1, rwork, info)
  
  if (info /= 0) then
    write(*,*) "problem in zgesvd routine"
    call exit
  end if
  
  lwork = int(prework(1))
  allocate(work(lwork))
  
  call zgesvd("a", "a", m, n, A, m, sigma, U, m, VT, n, work, lwork, rwork, info)
  
  if (info /= 0) then
    write(*,*) "problem in zgesvd routine"
    call exit
  end if
end subroutine svd_all_c_dp

subroutine svd_none_r_sp(A, sigma)
  real(sp), intent(inout)      :: A(:,:)
  real(sp), intent(out)        :: sigma(:)
  !local variables
  integer                         :: m, n, lwork, info
  real(sp)                     :: prework(1), U(1,1), VT(1,1)
  real(sp), allocatable        :: work(:)
  !
  m = size(A, 1)
  n = size(A, 2)
  !
  call sgesvd("o", "n", m, n, A, m, sigma, U, m, VT, n, prework, -1, info)
  !
  if (info /= 0) then
    write(*,*) "problem in svd - program stops now"
    call exit
  end if
  !
  lwork = int(prework(1))
  allocate(work(lwork))
  !
  call sgesvd("o", "n", m, n, A, m, sigma, U, m, VT, n, work, lwork, info)
  !
  if (info /= 0) then
    write(*,*) "problem in svd - program stops now"
    call exit
  end if
end subroutine svd_none_r_sp

subroutine svd_none_r_dp(A, sigma)
  real(dp), intent(inout)      :: A(:,:)
  real(dp), intent(out)        :: sigma(:)
  !local variables
  integer                         :: m, n, lwork, info
  real(dp)                     :: prework(1), U(1,1), VT(1,1)
  real(dp), allocatable        :: work(:)
  !
  m = size(A, 1)
  n = size(A, 2)
  !
  call dgesvd("o", "n", m, n, A, m, sigma, U, m, VT, n, prework, -1, info)
  !
  if (info /= 0) then
    write(*,*) "problem in svd - program stops now"
    call exit
  end if
  !
  lwork = int(prework(1))
  allocate(work(lwork))
  !
  call dgesvd("o", "n", m, n, A, m, sigma, U, m, VT, n, work, lwork, info)
  !
  if (info /= 0) then
    write(*,*) "problem in svd - program stops now"
    call exit
  end if
end subroutine svd_none_r_dp

subroutine svd_none_c_sp(A, sigma)
  complex(sp), intent(inout) :: A(:,:)
  real(sp), intent(out)      :: sigma(:)
  !local variables
  integer                    :: m, n, lwork, info
  real(sp), allocatable      :: rwork(:)
  complex(sp)                :: prework(1), U(1,1), VT(1,1)
  complex(sp), allocatable   :: work(:)
  
  m = size(A, 1)
  n = size(A, 2)

  allocate(rwork(5*min(n,m)))
  
  call cgesvd("o", "n", m, n, A, m, sigma, U, m, VT, n, prework, -1, rwork, info)
  
  if (info /= 0) then
    write(*,*) "problem in cgesvd routine"
    call exit
  end if
  
  lwork = int(prework(1))
  allocate(work(lwork))
  
  call cgesvd("o", "n", m, n, A, m, sigma, U, m, VT, n, work, lwork, rwork, info)

  if (info /= 0) then
    write(*,*) "problem in cgesvd routine"
    call exit
  end if
end subroutine svd_none_c_sp

subroutine svd_none_c_dp(A, sigma)
  complex(dp), intent(inout) :: A(:,:)
  real(dp), intent(out)      :: sigma(:)
  !local variables
  integer                    :: m, n, lwork, info
  real(dp), allocatable      :: rwork(:)
  complex(dp)                :: prework(1), U(1,1), VT(1,1)
  complex(dp), allocatable   :: work(:)
  
  m = size(A, 1)
  n = size(A, 2)

  allocate(rwork(5*min(n,m)))
  
  call zgesvd("o", "n", m, n, A, m, sigma, U, m, VT, n, prework, -1, rwork, info)
  
  if (info /= 0) then
    write(*,*) "problem in zgesvd routine"
    call exit
  end if
  
  lwork = int(prework(1))
  allocate(work(lwork))
  
  call zgesvd("o", "n", m, n, A, m, sigma, U, m, VT, n, work, lwork, rwork, info)

  if (info /= 0) then
    write(*,*) "problem in zgesvd routine"
    call exit
  end if
end subroutine svd_none_c_dp


!******************************************************************************** 
!       GETRF Routine - PLU Facorization 
!               INPUT:
!                       A       - matrix that need to be factorized
!               OUTPUT: 
!                       A       - strict lower triangular part of A represents L
!                                 L has unit diagonal elements and upper triangular 
!                                 part of A corresponds to U
!                       IPIV    - permuation array: row i of the matrix A was 
!                                 interchanged with the row IPIV(i)                
!******************************************************************************** 
SUBROUTINE GETRF_r_sp(A,IPIV)
  INTEGER:: IPIV(:)
  REAL(sp):: A(:,:)

  !local variables
  INTEGER:: M, N, LDA, INFO 

  M = SIZE(A(:,1))
  N = SIZE(A(1,:))
  LDA = M

  CALL SGETRF(M,N,A,LDA,IPIV,INFO)

  IF(INFO < 0) THEN
    STOP "Problem in SGETRF routine"
  ENDIF

END SUBROUTINE GETRF_r_sp

SUBROUTINE GETRF_r_dp(A,IPIV)
  INTEGER:: IPIV(:)
  REAL(dp):: A(:,:)

  !local variables
  INTEGER:: M, N, LDA, INFO

  M = SIZE(A(:,1))
  N = SIZE(A(1,:))
  LDA = M

  CALL DGETRF(M,N,A,LDA,IPIV,INFO)

  IF(INFO < 0) THEN
    STOP "Problem in DGETRF routine"
  ENDIF

END SUBROUTINE GETRF_r_dp

SUBROUTINE GETRF_c_sp(A,IPIV)
  INTEGER:: IPIV(:)
  COMPLEX(sp):: A(:,:)

  !local variables
  INTEGER:: M, N, LDA, INFO

  M = SIZE(A(:,1))
  N = SIZE(A(1,:))
  LDA = M

  CALL CGETRF(M,N,A,LDA,IPIV,INFO)

  IF(INFO < 0) THEN
    STOP "Problem in CGETRF routine"
  ENDIF

END SUBROUTINE GETRF_c_sp

SUBROUTINE GETRF_c_dp(A,IPIV)
  INTEGER:: IPIV(:)
  COMPLEX(dp):: A(:,:)

  !local variables
  INTEGER:: M, N, LDA, INFO

  M = SIZE(A(:,1))
  N = SIZE(A(1,:))
  LDA = M

  CALL ZGETRF(M,N,A,LDA,IPIV,INFO)

  IF(INFO < 0) THEN
    STOP "Problem in ZGETRF routine"
  ENDIF

END SUBROUTINE GETRF_c_dp


!******************************************************************************** 
!       getqr routine - qr factorization
!               input:
!                       a       - matrix that should be facorized
!               output: 
!                       a       - orthogonal matrix q computed by _geqrf and _orgqr
!                       det     - determinant of the matrix r 
!******************************************************************************** 
subroutine getqr_r_sp(A, R)
  real(sp), intent(inout) :: A(:,:)
  real(sp), intent(inout) :: R(:,:)
  !local variables
  integer                    :: i, j, m, n, lda, lwork, info, taudim
  real(sp)                :: prework(1)
  real(sp), allocatable   :: work(:), tau(:), work2(:)

  !determine dimensions
  m = size(A(:,1))
  n = size(A(1,:))
  lda = m
  
  R = 0.0_sp

  !determine dimension of tau
  taudim = min(m, n)
  allocate(tau(taudim))

  !determine optimal value of the array work
  call sgeqrf(m, n, A, lda, tau, prework, -1, info)

  !allocate work array
  lwork = int(prework(1))
  allocate(work(lwork))

  !call lapack routine now
  call sgeqrf(m, n, A, lda, tau, work, lwork, info)

  if (info /= 0) then
    stop "problem with sgeqrf routine"
  end if

  do i = 1, taudim
    do j = 1, i
      R(j,i) = A(j,i)
    end do
  end do
   
  !determine optimal value of the array work2
  call sorgqr(m, n, taudim, A, lda, tau, prework, -1, info)

  !allocate work array for sorgqr
  lwork = int(prework(1))
  allocate(work2(lwork))

  !call lapack routine
  call sorgqr(m, n, taudim, A, lda, tau, work2, lwork, info)

  if (info /= 0) then
    stop "problem in sorgqr routine"
  end if

end subroutine getqr_r_sp

subroutine getqr_r_dp(A, R)
  real(dp), intent(inout) :: A(:,:)
  real(dp), intent(inout) :: R(:,:)
  !local variables
  integer                    :: i, j, m, n, lda, lwork, info, taudim
  real(dp)                :: prework(1)
  real(dp), allocatable   :: work(:), tau(:), work2(:)

  !determine dimensions
  m = size(A(:,1))
  n = size(A(1,:))
  lda = m

  R = 0.0_dp

  !determine dimension of tau
  taudim = min(m, n)
  allocate(tau(taudim))

  !determine optimal value of the array work
  call dgeqrf(m, n, A, lda, tau, prework, -1, info)

  !allocate work array
  lwork = int(prework(1))
  allocate(work(lwork))

  !call lapack routine now
  call dgeqrf(m, n, A, lda, tau, work, lwork, info)

  if (info /= 0) then
    stop "problem with dgeqrf routine"
  end if

  do i = 1, taudim
    do j = 1, i
      R(j,i) = A(j,i)
    end do
  end do
   
  !determine optimal value of the array work2
  call dorgqr(m, n, taudim, A, lda, tau, prework, -1, info)

  !allocate work array for sorgqr
  lwork = int(prework(1))
  allocate(work2(lwork))

  !call lapack routine
  call dorgqr(m, n, taudim, A, lda, tau, work2, lwork, info)

  if (info /= 0) then
    stop "problem in dorgqr routine"
  end if

end subroutine getqr_r_dp

subroutine getqr_c_sp(A, R)
  complex(sp), intent(inout) :: A(:,:)
  complex(sp), intent(inout) :: R(:,:)
  !local variables
  integer                       :: i, j, m, n, lda, lwork, info, taudim
  complex(sp)                :: prework(1)
  complex(sp), allocatable   :: work(:), tau(:), work2(:)

  !determine dimensions
  m = size(a(:,1))
  n = size(a(1,:))
  lda = m

  R = (0.0_sp, 0.0_sp)

  !determine dimension of tau
  taudim = min(m, n)
  allocate(tau(taudim))

  !determine optimal value of the array work
  call cgeqrf(m, n, A, lda, tau, prework, -1, info)

  !allocate work array
  lwork = int(prework(1))
  allocate(work(lwork))

  !call lapack routine now
  call cgeqrf(m, n, A, lda, tau, work, lwork, info)

  if (info /= 0) then
    stop "problem with cgeqrf routine"
  end if

  do i = 1,taudim
    do j = 1, i
      R(j,i) = A(j,i)
    end do
  end do
   
  !determine optimal value of the array work2
  call cungqr(m, n, taudim, A, lda, tau, prework, -1, info)

  !allocate work array for sorgqr
  lwork = int(prework(1))
  allocate(work2(lwork))

  !call lapack routine
  call cungqr(m, n, taudim, A, lda, tau, work2, lwork, info)

  if (info /= 0) then
    stop "problem in cungqr routine"
  end if

end subroutine getqr_c_sp

subroutine getqr_c_dp(A, R)
  complex(dp), intent(inout) :: A(:,:)
  complex(dp), intent(inout) :: R(:,:)
  !local variables
  integer                       :: i, j, m, n, lda, lwork, info, taudim
  complex(dp)                :: prework(1)
  complex(dp), allocatable   :: work(:), tau(:), work2(:)

  !determine dimensions
  m = size(a(:,1))
  n = size(a(1,:))
  lda = m

  R = (0.0_dp, 0.0_dp)

  !determine dimension of tau
  taudim = min(m, n)
  allocate(tau(taudim))

  !determine optimal value of the array work
  call zgeqrf(m, n, A, lda, tau, prework, -1, info)

  !allocate work array
  lwork = int(prework(1))
  allocate(work(lwork))

  !call lapack routine now
  call zgeqrf(m, n, A, lda, tau, work, lwork, info)

  if (info /= 0) then
    stop "problem with zgeqrf routine"
  end if

  do i = 1, taudim
    do j = 1, i
      R(j,i) = A(j,i)
    end do
  end do
   
  !determine optimal value of the array work2
  call zungqr(m, n, taudim, A, lda, tau, prework, -1, info)

  !allocate work array for sorgqr
  lwork = int(prework(1))
  allocate(work2(lwork))

  !call lapack routine
  call zungqr(m, n, taudim, A, lda, tau, work2, lwork, info)

  if (info /= 0) then
    stop "problem in zungqr routine"
  end if

end subroutine getqr_c_dp


!******************************************************************************** 
!       GETQR Routine - QR factorization
!               INPUT:
!                       A       - Matrix that should be facorized
!               OUTPUT: 
!                       A       - Orthogonal matrix Q computed by _GEQRF and _ORGQR
!                       DET     - Determinant of the matrix R 
!******************************************************************************** 
SUBROUTINE GETQR_det_r_sp(A,det)
  REAL(sp):: A(:,:)
  REAL(sp):: det
  
  !local variables
  INTEGER:: i, M, N, LDA, LWORK, INFO, TAUDIM
  REAL(sp):: PREWORK(1)
  REAL(sp),ALLOCATABLE:: WORK(:), TAU(:), WORK2(:)

  det = 1.0_sp

  !Determine dimensions
  M = SIZE(A(:,1))
  N = SIZE(A(1,:))
  LDA = M

  !Determine dimension of TAU
  TAUDIM = MIN(M,N)
  ALLOCATE(TAU(TAUDIM))

  !Determine optimal value of the array WORK
  CALL SGEQRF(M,N,A,LDA,TAU,PREWORK,-1,INFO)

  !Allocate WORK array
  LWORK = INT(PREWORK(1))
  ALLOCATE(WORK(LWORK))

  !Call LAPACK routine now
  CALL SGEQRF(M,N,A,LDA,TAU,WORK,LWORK,INFO)

  IF(INFO/=0) THEN
    STOP "Problem with SGEQRF routine"
  ENDIF

  DO i=1,TAUDIM
    det = det * A(i,i) 
  ENDDO
   
  !Determine optimal value of the array WORK2
  CALL SORGQR(M,N,TAUDIM,A,LDA,TAU,PREWORK,-1,INFO)

  !Allocate WORK array for SORGQR
  LWORK = INT(PREWORK(1))
  ALLOCATE(WORK2(LWORK))

  !Call LAPACK routine
  CALL SORGQR(M,N,TAUDIM,A,LDA,TAU,WORK2,LWORK,INFO)

  IF(INFO/=0) THEN
    STOP "Problem in SORGQR routine"
  ENDIF

END SUBROUTINE GETQR_det_r_sp

SUBROUTINE GETQR_det_r_dp(A,det)
  REAL(dp):: A(:,:)
  REAL(dp):: det
  
  !local variables
  INTEGER:: i, M, N, LDA, LWORK, INFO, TAUDIM
  REAL(dp):: PREWORK(1)
  REAL(dp),ALLOCATABLE:: WORK(:), TAU(:), WORK2(:)

  det = 1.0_dp

  !Determine dimensions
  M = SIZE(A(:,1))
  N = SIZE(A(1,:))
  LDA = M

  !Determine dimension of TAU
  TAUDIM = MIN(M,N)
  ALLOCATE(TAU(TAUDIM))

  !Determine optimal value of the array WORK
  CALL DGEQRF(M,N,A,LDA,TAU,PREWORK,-1,INFO)

  !Allocate WORK array
  LWORK = INT(PREWORK(1))
  ALLOCATE(WORK(LWORK))

  !Call LAPACK routine now
  CALL DGEQRF(M,N,A,LDA,TAU,WORK,LWORK,INFO)

  IF(INFO/=0) THEN
    STOP "Problem in DGEQRF routine"
  ENDIF

  DO i=1,TAUDIM
    det = det * A(i,i) 
  ENDDO
   
  !Determine optimal value of the array WORK2
  CALL DORGQR(M,N,TAUDIM,A,LDA,TAU,PREWORK,-1,INFO)

  !Allocate WORK array for SORGQR
  LWORK = INT(PREWORK(1))
  ALLOCATE(WORK2(LWORK))

  !Call LAPACK routine
  CALL DORGQR(M,N,TAUDIM,A,LDA,TAU,WORK2,LWORK,INFO)

  IF(INFO/=0) THEN
    STOP "Problem in DORGQR routine"
  ENDIF

END SUBROUTINE GETQR_det_r_dp

SUBROUTINE GETQR_det_c_sp(A,det)
  COMPLEX(sp):: A(:,:)
  COMPLEX(sp):: det
  
  !local variables
  INTEGER:: i, M, N, LDA, LWORK, INFO, TAUDIM
  COMPLEX(sp):: PREWORK(1)
  COMPLEX(sp),ALLOCATABLE:: WORK(:), TAU(:), WORK2(:)

  det = (1.0_sp,0.0_sp)

  !Determine dimensions
  M = SIZE(A(:,1))
  N = SIZE(A(1,:))
  LDA = M

  !Determine dimension of TAU
  TAUDIM = MIN(M,N)
  ALLOCATE(TAU(TAUDIM))

  !Determine optimal value of the array WORK
  CALL CGEQRF(M,N,A,LDA,TAU,PREWORK,-1,INFO)

  !Allocate WORK array
  LWORK = INT(PREWORK(1))
  ALLOCATE(WORK(LWORK))

  !Call LAPACK routine now
  CALL CGEQRF(M,N,A,LDA,TAU,WORK,LWORK,INFO)

  IF(INFO/=0) THEN
    STOP "Problem with CGEQRF routine"
  ENDIF

  DO i=1,TAUDIM
    det = det * A(i,i) 
  ENDDO
   
  !Determine optimal value of the array WORK2
  CALL CUNGQR(M,N,TAUDIM,A,LDA,TAU,PREWORK,-1,INFO)

  !Allocate WORK array for SORGQR
  LWORK = INT(PREWORK(1))
  ALLOCATE(WORK2(LWORK))

  !Call LAPACK routine
  CALL CUNGQR(M,N,TAUDIM,A,LDA,TAU,WORK2,LWORK,INFO)

  IF(INFO/=0) THEN
    STOP "Problem in CUNGQR routine"
  ENDIF

END SUBROUTINE GETQR_det_c_sp

SUBROUTINE GETQR_det_c_dp(A,det)
  COMPLEX(dp):: A(:,:)
  COMPLEX(dp):: det
  
  !local variables
  INTEGER:: i, M, N, LDA, LWORK, INFO, TAUDIM
  COMPLEX(dp):: PREWORK(1)
  COMPLEX(dp),ALLOCATABLE:: WORK(:), TAU(:), WORK2(:)

  det = (1.0_dp,0.0_dp)

  !Determine dimensions
  M = SIZE(A(:,1))
  N = SIZE(A(1,:))
  LDA = M

  !Determine dimension of TAU
  TAUDIM = MIN(M,N)
  ALLOCATE(TAU(TAUDIM))

  !Determine optimal value of the array WORK
  CALL ZGEQRF(M,N,A,LDA,TAU,PREWORK,-1,INFO)

  !Allocate WORK array
  LWORK = INT(PREWORK(1))
  ALLOCATE(WORK(LWORK))

  !Call LAPACK routine now
  CALL ZGEQRF(M,N,A,LDA,TAU,WORK,LWORK,INFO)

  IF(INFO/=0) THEN
    STOP "Problem with ZGEQRF routine"
  ENDIF

  DO i=1,TAUDIM
    det = det * A(i,i) 
  ENDDO
   
  !Determine optimal value of the array WORK2
  CALL ZUNGQR(M,N,TAUDIM,A,LDA,TAU,PREWORK,-1,INFO)

  !Allocate WORK array for SORGQR
  LWORK = INT(PREWORK(1))
  ALLOCATE(WORK2(LWORK))

  !Call LAPACK routine
  CALL ZUNGQR(M,N,TAUDIM,A,LDA,TAU,WORK2,LWORK,INFO)

  IF(INFO/=0) THEN
    STOP "Problem in ZUNGQR routine"
  ENDIF

END SUBROUTINE GETQR_det_c_dp


!******************************************************************************** 
!       Cholesky decomposition 
!               INPUT:
!                       A       - matrix that need to be Cholesky decomposed
!                       UPLO    - L for LL^+ version and U for U^+U version
!               OUTPUT:
!                       A       - IF UPLO = L lower triangular part of A is L and
!                                 if UPLO = U upper triangular part of A is U.  
!******************************************************************************** 
subroutine potrf_r_sp(a, uplo, det)
  real(sp), intent(inout)         :: a(:,:)
  character, optional, intent(in)    :: uplo
  real(sp), optional, intent(out) :: det
  !local variables
  character                          :: uplo_
  integer                            :: i, n, m, info
  
  n = size(a, 1)
  m = size(a, 2)

  if (n /= m) then
    write(*,*) "non-square matrix can't be cholesky decomposed"
    call exit
  end if

  if (present(uplo)) then
    uplo_ = uplo
  else
    uplo_ = "u"
  end if

  call spotrf(uplo_, n, a, n, info)

  if (info /= 0) then
    stop "problem in spotrf routine"
  end if

  if (present(det)) then
    det = 1.0_sp
    do i = 1, n
      det = det * a(i,i)
    end do
  end if
end subroutine potrf_r_sp

subroutine potrf_r_dp(a, uplo, det)
  real(dp), intent(inout)         :: a(:,:)
  character, optional, intent(in)    :: uplo
  real(dp), optional, intent(out) :: det
  !local variables
  character                          :: uplo_
  integer                            :: i, n, m, info
  
  n = size(a, 1)
  m = size(a, 2)

  if(n /= m) then
    write(*,*) "non-square matrix can't be cholesky decomposed"
    call exit
  end if

  if (present(uplo)) then
    uplo_ = uplo
  else
    uplo_ = "u"
  end if

  call dpotrf(uplo_, n, a, n, info)

  if (info /= 0) then
    stop "problem in spotrf routine"
  end if

  if (present(det)) then
    det = 1.0_dp
    do i = 1, n
      det = det * a(i,i)
    end do
  end if
end subroutine potrf_r_dp

subroutine potrf_c_sp(a, uplo, det)
  complex(sp), intent(inout)         :: a(:,:)
  character, optional, intent(in)       :: uplo
  complex(sp), optional, intent(out) :: det
  !local variables
  character                             :: uplo_
  integer                               :: i, n, m, info
  
  n = size(a, 1)
  m = size(a, 2)

  if(n /= m) then
    write(*,*) "non-square matrix can't be cholesky decomposed"
    call exit
  end if

  if (present(uplo)) then
    uplo_ = uplo
  else
    uplo_ = "u"
  end if

  call cpotrf(uplo_, n, a, n, info)

  if (info /= 0) then
    stop "problem in spotrf routine"
  end if

  if (present(det)) then
    det = (1.0_sp, 0.0_sp)
    do i = 1, n
      det = det * a(i,i)
    end do
  end if
end subroutine potrf_c_sp

subroutine potrf_c_dp(a, uplo, det)
  complex(dp), intent(inout)         :: a(:,:)
  character, optional, intent(in)       :: uplo
  complex(dp), optional, intent(out) :: det
  !local variables
  character                             :: uplo_
  integer                               :: i, n, m, info
  
  n = size(a, 1)
  m = size(a, 2)

  if (n /= m) then
    write(*,*) "non-square matrix can't be cholesky decomposed"
    call exit
  end if

  if (present(uplo)) then
    uplo_ = uplo
  else
    uplo_ = "u"
  end if

  call zpotrf(uplo_, n, a, n, info)

  if (info /= 0) then
    stop "problem in spotrf routine"
  end if

  if (present(det)) then
    det = (1.0_dp, 0.0_dp)
    do i = 1, n
      det = det * a(i,i)
    end do
  end if
end subroutine potrf_c_dp


!******************************************************************************** 
!       INVERSE Routine - Uses Lapack routines to inverse a Matrix A
!               INPUT:
!                       A       - Matrix that need to be inverted
!               OUTPUT:
!                       A       - inverted matrix 
!******************************************************************************** 
SUBROUTINE INVERSE_r_sp(A)
  REAL(sp):: A(:,:)
  !local variables
  INTEGER:: N, M, INFORM, LWORK
  INTEGER,ALLOCATABLE:: IPIV(:)
  REAL(sp),ALLOCATABLE:: WORK(:)
  REAL(sp):: PREWORK(1)

  N = SIZE(A(:,1))
  M = SIZE(A(1,:))

  IF(N/=M) THEN
    STOP "Problem in INV: Matrix is not square matrix"
  ENDIF

  !Allocate IPIV for PLU decomposition
  ALLOCATE(IPIV(N))

  !Call PLU decomposition
  CALL SGETRF(N,N,A,N,IPIV,INFORM)
  
  IF(INFORM/=0) THEN
    STOP "LU factorization fails - problem in SGETRF routine"
  ENDIF

  !Determine optimal size of WORK array
  CALL SGETRI(N,A,N,IPIV,PREWORK,-1,INFORM)

  !Allocate WORK array
  LWORK = INT(PREWORK(1))
  ALLOCATE(WORK(LWORK))  

  CALL SGETRI(N,A,N,IPIV,WORK,LWORK,INFORM)

  IF(INFORM/=0) THEN
    STOP "Inversion of the matrix fails - problem in SGETRI routine"
  ENDIF

END SUBROUTINE INVERSE_r_sp

SUBROUTINE INVERSE_r_dp(A)
  REAL(dp):: A(:,:)
  !local variables
  INTEGER:: N, M, INFORM, LWORK
  INTEGER,ALLOCATABLE:: IPIV(:)
  REAL(dp),ALLOCATABLE:: WORK(:)
  REAL(dp):: PREWORK(1)

  N = SIZE(A(:,1))
  M = SIZE(A(1,:))

  IF(N/=M) THEN
    STOP "Problem in INV: Matrix is not square matrix"
  ENDIF

  !Allocate IPIV for PLU decomposition
  ALLOCATE(IPIV(N))

  !Call PLU decomposition
  CALL DGETRF(N,N,A,N,IPIV,INFORM)
  
  IF(INFORM/=0) THEN
    STOP "LU factorization fails - problem in DGETRF routine"
  ENDIF

  !Determine optimal size of WORK array
  CALL DGETRI(N,A,N,IPIV,PREWORK,-1,INFORM)

  !Allocate WORK array
  LWORK = INT(PREWORK(1))
  ALLOCATE(WORK(LWORK))  

  CALL DGETRI(N,A,N,IPIV,WORK,LWORK,INFORM)

  IF(INFORM/=0) THEN
    STOP "Inversion of the matrix fails - problem in DGETRI routine"
  ENDIF

END SUBROUTINE INVERSE_r_dp

SUBROUTINE INVERSE_c_sp(A)
  COMPLEX(sp):: A(:,:)
  !local variables
  INTEGER:: N, M, INFORM, LWORK
  INTEGER,ALLOCATABLE:: IPIV(:)
  COMPLEX(sp),ALLOCATABLE:: WORK(:)
  COMPLEX(sp):: PREWORK(1)

  N = SIZE(A(:,1))
  M = SIZE(A(1,:))

  IF(N/=M) THEN
    STOP "Problem in INV: Matrix is not square matrix"
  ENDIF

  !Allocate IPIV for PLU decomposition
  ALLOCATE(IPIV(N))

  !Call PLU decomposition
  CALL CGETRF(N,N,A,N,IPIV,INFORM)
  
  IF(INFORM/=0) THEN
    STOP "LU factorization fails - problem in CGETRF routine"
  ENDIF

  !Determine optimal size of WORK array
  CALL CGETRI(N,A,N,IPIV,PREWORK,-1,INFORM)

  !Allocate WORK array
  LWORK = INT(PREWORK(1))
  ALLOCATE(WORK(LWORK))  

  CALL CGETRI(N,A,N,IPIV,WORK,LWORK,INFORM)

  IF(INFORM/=0) THEN
    STOP "Inversion of the matrix fails - problem in CGETRI routine"
  ENDIF

END SUBROUTINE INVERSE_c_sp

SUBROUTINE INVERSE_c_dp(A)
  COMPLEX(dp):: A(:,:)
  !local variables
  INTEGER:: N, M, INFORM, LWORK
  INTEGER,ALLOCATABLE:: IPIV(:)
  COMPLEX(dp),ALLOCATABLE:: WORK(:)
  COMPLEX(dp):: PREWORK(1)

  N = SIZE(A(:,1))
  M = SIZE(A(1,:))

  IF(N/=M) THEN
    STOP "Problem in INV: Matrix is not square matrix"
  ENDIF

  !Allocate IPIV for PLU decomposition
  ALLOCATE(IPIV(N))

  !Call PLU decomposition
  CALL ZGETRF(N,N,A,N,IPIV,INFORM)
  
  IF(INFORM/=0) THEN
    STOP "LU factorization fails - problem in ZGETRF routine"
  ENDIF

  !Determine optimal size of WORK array
  CALL ZGETRI(N,A,N,IPIV,PREWORK,-1,INFORM)

  !Allocate WORK array
  LWORK = INT(PREWORK(1))
  ALLOCATE(WORK(LWORK))  

  CALL ZGETRI(N,A,N,IPIV,WORK,LWORK,INFORM)

  IF(INFORM/=0) THEN
    STOP "Inversion of the matrix fails - problem in ZGETRI routine"
  ENDIF

END SUBROUTINE INVERSE_c_dp


!******************************************************************************** 
!
! Calculate matrix to the power p, for any real exponent p
!
!    A --> A^{p} = U \Sigma^{p} U^{\dagger}
!
!******************************************************************************** 
subroutine matrix_power_c_dp(A, p)
  complex(dp), intent(inout) :: A(:,:)
  real(dp), intent(in)       :: p
  !local
  integer                    :: n, m, j
  real(dp), allocatable      :: w(:)
  complex(dp), allocatable   :: U(:,:), Uw(:,:)

  n = size(A, 1)
  m = size(A, 2)

  if (n /= m) stop "matrix A in matrix_power_c_dp must be a square matrix"

  allocate(w(n))
  call syev(A, w)

  w = w**p

  allocate(U, source=A)
  allocate(Uw, mold=A)

  do j = 1, n
    Uw(:,j) = U(:,j) * w(j)
  end do

  call gemm("n", "c", n, n, n, onec_dp, Uw, n, U, n, zeroc_dp, A, n)
end subroutine matrix_power_c_dp


!******************************************************************************** 
!
! Calculate Moore-Penrose inverse of the matrix A
!
!    A = U Sigma V^{\dagger}
!
!    A' = V Sigma' U^{\dagger}
!
!    Sigma'_i  = 1 / Sigma_i if Sigma_i > thresh
!    Sigma'_i  = 0 if Sigma_i <= thresh
!
!******************************************************************************** 
subroutine pseudoinverse_r_dp(A, thresh)
  real(wp), intent(inout)    :: A(:,:)
  real(wp), intent(in)       :: thresh
  !local
  integer                    :: i, n
  real(wp), allocatable      :: sigma(:)
  real(wp), allocatable      :: U(:,:), VT(:,:)

  n = size(A, 1)
  allocate(sigma(n), U(n,n), VT(n,n))

  call svd(A, sigma, U, VT)

  where (sigma < thresh) 
    sigma = 0.0_wp
  else where
    sigma = 1.0_wp / sigma
  end where

  do i = 1, size(VT, 2)
    VT(:,i) = sigma * VT(:,i)
  end do

  call gemm("t", "t", n, n, n, oner, VT, n, U, n, zeror, A, n)
end subroutine pseudoinverse_r_dp

subroutine pseudoinverse_c_dp(A, thresh)
  complex(wp), intent(inout) :: A(:,:)
  real(wp), intent(in)       :: thresh
  !local
  integer                    :: i, n
  real(wp), allocatable      :: sigma(:)
  complex(wp), allocatable   :: U(:,:), VT(:,:)

  n = size(A, 1)
  allocate(sigma(n), U(n,n), VT(n,n))

  call svd(A, sigma, U, VT)

  where (sigma < thresh) 
    sigma = 0.0_wp
  else where
    sigma = 1.0_wp / sigma
  end where

  do i = 1, size(VT, 2)
    VT(:,i) = sigma * VT(:,i)
  end do

  call gemm("c", "c", n, n, n, onec, VT, n, U, n, zeroc, A, n)
end subroutine pseudoinverse_c_dp


!******************************************************************************** 
!
! Calculate Moore-Penrose inverse of the matrix A using Tikhonov regularization
!
!    A = U Sigma V^{\dagger}
!
!    A' = V Sigma' U^{\dagger}
!
!    Sigma'_i  = sigma_i / (sigma_i^2 + \lambda^2)
!
!******************************************************************************** 
subroutine pseudoinverse_tikhonov_c_dp(A, lambda)
  complex(wp), intent(inout) :: A(:,:)
  real(wp), intent(in)       :: lambda
  !local
  integer                    :: i, n
  real(wp), allocatable      :: sigma(:)
  complex(wp), allocatable   :: U(:,:), VT(:,:)

  n = size(A, 1)
  allocate(sigma(n), U(n,n), VT(n,n))

  call svd(A, sigma, U, VT)

  sigma = sigma / (sigma**2 + lambda**2)

  do i = 1, size(VT, 2)
    VT(:,i) = sigma * VT(:,i)
  end do

  call gemm("c", "c", n, n, n, onec, VT, n, U, n, zeroc, A, n)
end subroutine pseudoinverse_tikhonov_c_dp


!******************************************************************************** 
!       Inversion of triangular matrices
!               INPUT:
!                       A       - upper or lower triangular matrix
!                       UPLO    - U for upper and L  for lower triangular matrix
!                       DIAG    - N for non-unit and U for unit triangular matrix
!               OUTPUT:
!                       A       - inverse of the matrix in same format as input
!******************************************************************************** 
SUBROUTINE INVERSE_TRI_r_sp(A,UPLO,DIAG)
  REAL(sp):: A(:,:)
  CHARACTER,OPTIONAL:: UPLO, DIAG
  !local variables
  INTEGER:: N, M
  INTEGER:: INFO
  CHARACTER:: myuplo, mydiag
 
  N = SIZE(A(:,1))
  M = SIZE(A(1,:))

  IF(N/=M) THEN
    STOP "Input matrix is not square matrix"
  ENDIF

  IF(PRESENT(UPLO)) THEN
    myuplo = UPLO
  ELSE
    myuplo = "U"
  ENDIF

  IF(PRESENT(DIAG)) THEN
    mydiag = DIAG
  ELSE
    mydiag = "N"
  ENDIF

  CALL STRTRI(myuplo,mydiag,N,A,N,INFO)
    
  IF(INFO/=0) THEN
    STOP "Inverion of triangular matrix fails - problem in STRTRI routine"
  ENDIF

END SUBROUTINE INVERSE_TRI_r_sp

SUBROUTINE INVERSE_TRI_r_dp(A,UPLO,DIAG)
  REAL(dp):: A(:,:)
  CHARACTER,OPTIONAL:: UPLO, DIAG
  !local variables
  INTEGER:: N, M
  INTEGER:: INFO
  CHARACTER:: myuplo, mydiag
 
  N = SIZE(A(:,1))
  M = SIZE(A(1,:))

  IF(N/=M) THEN
    STOP "Input matrix is not square matrix"
  ENDIF

  IF(PRESENT(UPLO)) THEN
    myuplo = UPLO
  ELSE
    myuplo = "U"
  ENDIF

  IF(PRESENT(DIAG)) THEN
    mydiag = DIAG
  ELSE
    mydiag = "N"
  ENDIF

  CALL DTRTRI(myuplo,mydiag,N,A,N,INFO)
    
  IF(INFO/=0) THEN
    STOP "Inverion of triangular matrix fails - problem in DTRTRI routine"
  ENDIF

END SUBROUTINE INVERSE_TRI_r_dp

SUBROUTINE INVERSE_TRI_c_sp(A,UPLO,DIAG)
  COMPLEX(sp):: A(:,:)
  CHARACTER,OPTIONAL:: UPLO, DIAG
  !local variables
  INTEGER:: N, M
  INTEGER:: INFO
  CHARACTER:: myuplo, mydiag
 
  N = SIZE(A(:,1))
  M = SIZE(A(1,:))

  IF(N/=M) THEN
    STOP "Input matrix is not square matrix"
  ENDIF

  IF(PRESENT(UPLO)) THEN
    myuplo = UPLO
  ELSE
    myuplo = "U"
  ENDIF

  IF(PRESENT(DIAG)) THEN
    mydiag = DIAG
  ELSE
    mydiag = "N"
  ENDIF

  CALL CTRTRI(myuplo,mydiag,N,A,N,INFO)
    
  IF(INFO/=0) THEN
    STOP "Inverion of triangular matrix fails - problem in CTRTRI routine"
  ENDIF

END SUBROUTINE INVERSE_TRI_c_sp

SUBROUTINE INVERSE_TRI_c_dp(A,UPLO,DIAG)
  COMPLEX(dp):: A(:,:)
  CHARACTER,OPTIONAL:: UPLO, DIAG
  !local variables
  INTEGER:: N, M
  INTEGER:: INFO
  CHARACTER:: myuplo, mydiag
 
  N = SIZE(A(:,1))
  M = SIZE(A(1,:))

  IF(N/=M) THEN
    STOP "Input matrix is not square matrix"
  ENDIF

  IF(PRESENT(UPLO)) THEN
    myuplo = UPLO
  ELSE
    myuplo = "U"
  ENDIF

  IF(PRESENT(DIAG)) THEN
    mydiag = DIAG
  ELSE
    mydiag = "N"
  ENDIF

  CALL ZTRTRI(myuplo,mydiag,N,A,N,INFO)
    
  IF(INFO/=0) THEN
    STOP "Inverion of triangular matrix fails - problem in ZTRTRI routine"
  ENDIF

END SUBROUTINE INVERSE_TRI_c_dp


!******************************************************************************** 
!       This subroutine calculates determinant of the square matrix
!               via PLU factorization
!               INPUT:
!                       A       - matrix which determinant is needed
!               OUTPUT:
!                       DET     - determinant of the matrix A
!******************************************************************************** 
SUBROUTINE DETERMINANT_r_sp(A,DET)
  REAL(sp):: A(:,:)
  REAL(sp):: DET

  !local variables
  INTEGER,ALLOCATABLE:: PIV(:)
  INTEGER:: N, M, i
  REAL(sp):: fact
  INTEGER:: Numb

  N = SIZE(A(:,1))
  M = SIZE(A(1,:))
 
  IF(N/=M) THEN
    STOP "DETERMINAT will fail: matrix is not square"
  ENDIF

  !Allocate PIV array
  ALLOCATE(PIV(N))
 
  !Do PLU factorization
  CALL GETRF_r_sp(A,PIV)

  DET = 1.0_sp
 
  DO i=1,N
    DET = DET * A(i,i)
  ENDDO
  
  Numb = PERMUTE(PIV)

  IF(mod(Numb,2)==0) THEN
    fact = 1.0_sp
  ELSE
    fact = -1.0_sp
  ENDIF

  DET = DET * fact

END SUBROUTINE DETERMINANT_r_sp

SUBROUTINE DETERMINANT_r_dp(A,DET)
  REAL(dp):: A(:,:)
  REAL(dp):: DET

  !local variables
  INTEGER,ALLOCATABLE:: PIV(:)
  INTEGER:: N, M, i
  REAL(dp):: fact
  INTEGER:: Numb

  N = SIZE(A(:,1))
  M = SIZE(A(1,:))
 
  IF(N/=M) THEN
    STOP "DETERMINAT will fail: matrix is not square"
  ENDIF

  !Allocate PIV array
  ALLOCATE(PIV(N))
 
  !Do PLU factorization
  CALL GETRF_r_dp(A,PIV)

  DET = 1.0_dp
 
  DO i=1,N
    DET = DET * A(i,i)
  ENDDO
  
  Numb = PERMUTE(PIV)
  IF(mod(Numb,2)==0) THEN
    fact = 1.0_dp
  ELSE
    fact = -1.0_dp
  ENDIF

  DET = DET * fact

END SUBROUTINE DETERMINANT_r_dp

SUBROUTINE DETERMINANT_c_sp(A,DET)
  COMPLEX(sp):: A(:,:)
  COMPLEX(sp):: DET

  !local variables
  INTEGER,ALLOCATABLE:: PIV(:)
  INTEGER:: N, M, i
  COMPLEX(sp):: fact
  INTEGER:: Numb

  N = SIZE(A(:,1))
  M = SIZE(A(1,:))
 
  IF(N/=M) THEN
    STOP "DETERMINAT will fail: matrix is not square"
  ENDIF

  !Allocate PIV array
  ALLOCATE(PIV(N))
 
  !Do PLU factorization
  CALL GETRF_c_sp(A,PIV)

  DET = (1.0_sp,0.0_sp)
 
  DO i=1,N
    DET = DET * A(i,i)
  ENDDO
  
  Numb = PERMUTE(PIV)
  IF(mod(Numb,2)==0) THEN
    fact = (1.0_sp,0.0_sp)
  ELSE
    fact = (-1.0_sp,0.0_sp)
  ENDIF

  DET = DET * fact

END SUBROUTINE DETERMINANT_c_sp

SUBROUTINE DETERMINANT_c_dp(A,DET)
  COMPLEX(dp):: A(:,:)
  COMPLEX(dp):: DET

  !local variables
  INTEGER,ALLOCATABLE:: PIV(:)
  INTEGER:: N, M, i
  COMPLEX(dp):: fact
  INTEGER:: Numb
 
  N = SIZE(A(:,1))
  M = SIZE(A(1,:))
 
  IF(N/=M) THEN
    STOP "DETERMINAT will fail: matrix is not square"
  ENDIF

  !Allocate PIV array
  ALLOCATE(PIV(N))
 
  !Do PLU factorization
  CALL GETRF_c_dp(A,PIV)

  DET = (1.0_dp,0.0_dp)
 
  DO i=1,N
    DET = DET * A(i,i)
  ENDDO
  
  Numb = PERMUTE(PIV)
  IF(mod(Numb,2)==0) THEN
    fact = (1.0_dp,0.0_dp)
  ELSE
    fact = (-1.0_dp,0.0_dp)
  ENDIF

  DET = DET * fact

END SUBROUTINE DETERMINANT_c_dp


!******************************************************************************** 
!       This subroutine solves the system of linear equations A X = B
!               INPUT:
!                       A       -       square matrix
!               INOUT:  
!                       X       -       On entry res.vecto; on exit solution
!******************************************************************************** 
SUBROUTINE SOLVE_LS_r_sp(A,X)
  REAL(sp),INTENT(IN):: A(:,:)
  REAL(sp),INTENT(INOUT):: X(:,:)

  !local variables
  INTEGER:: N,M,NB
  INTEGER:: INFO, LWORK
  INTEGER,ALLOCATABLE:: IPIV(:)
  REAL(sp),ALLOCATABLE:: WORK(:)
  REAL(sp):: PREWORK(1)  

  N = SIZE(A(:,1))
  M = SIZE(A(1,:))

  IF(N/=M) THEN
    STOP "MATRIX A IN LINEAR SOLVER IS NOT DIAGONAL"
  ENDIF
  
  NB = SIZE(X(1,:))

  !Allovate pivot array
  ALLOCATE(IPIV(N))

  !Determine optimal size of the array WORK
  CALL SSYSV("U",N,NB,A,N,IPIV,X,N,PREWORK,-1,INFO)

  LWORK = INT(PREWORK(1))
  ALLOCATE(WORK(LWORK))

  !Solve linear system of equations
  CALL SSYSV("U",N,NB,A,N,IPIV,X,N,WORK,LWORK,INFO)

  IF(INFO/=0) THEN
    STOP "SOLVE LINEAR SYSTEM DOES NOT WORK"
  ENDIF  
END SUBROUTINE SOLVE_LS_r_sp

SUBROUTINE SOLVE_LS_r_dp(A,X)
  REAL(dp),INTENT(IN):: A(:,:)
  REAL(dp),INTENT(INOUT):: X(:,:)

  !local variables
  INTEGER:: N,M,NB
  INTEGER:: INFO, LWORK
  INTEGER,ALLOCATABLE:: IPIV(:)
  REAL(dp),ALLOCATABLE:: WORK(:)
  REAL(dp):: PREWORK(1)  

  N = SIZE(A(:,1))
  M = SIZE(A(1,:))

  IF(N/=M) THEN
    STOP "MATRIX A IN LINEAR SOLVER IS NOT DIAGONAL"
  ENDIF
  
  NB = SIZE(X(1,:))

  !Allovate pivot array
  ALLOCATE(IPIV(N))

  !Determine optimal size of the array WORK
  CALL DSYSV("U",N,NB,A,N,IPIV,X,N,PREWORK,-1,INFO)

  LWORK = INT(PREWORK(1))
  ALLOCATE(WORK(LWORK))

  !Solve linear system of equations
  CALL DSYSV("U",N,NB,A,N,IPIV,X,N,WORK,LWORK,INFO)

  IF(INFO/=0) THEN
    STOP "SOLVE LINEAR SYSTEM DOES NOT WORK"
  ENDIF  
END SUBROUTINE SOLVE_LS_r_dp

SUBROUTINE SOLVE_LS_c_sp(A,X)
  COMPLEX(sp),INTENT(IN):: A(:,:)
  COMPLEX(sp),INTENT(INOUT):: X(:,:)

  !local variables
  INTEGER:: N,M,NB
  INTEGER:: INFO, LWORK
  INTEGER,ALLOCATABLE:: IPIV(:)
  COMPLEX(sp),ALLOCATABLE:: WORK(:)
  COMPLEX(sp):: PREWORK(1)  

  N = SIZE(A(:,1))
  M = SIZE(A(1,:))

  IF(N/=M) THEN
    STOP "MATRIX A IN LINEAR SOLVER IS NOT DIAGONAL"
  ENDIF
  
  NB = SIZE(X(1,:))

  !Allovate pivot array
  ALLOCATE(IPIV(N))

  !Determine optimal size of the array WORK
  !CALL CHESV("U",N,NB,A,N,IPIV,X,N,PREWORK,-1,INFO)
  CALL CSYSV("U",N,NB,A,N,IPIV,X,N,PREWORK,-1,INFO)

  LWORK = INT(PREWORK(1))
  ALLOCATE(WORK(LWORK))

  !Solve linear system of equations
  !CALL CHESV("U",N,NB,A,N,IPIV,X,N,WORK,LWORK,INFO)
  CALL CSYSV("U",N,NB,A,N,IPIV,X,N,WORK,LWORK,INFO)

  IF(INFO/=0) THEN
    STOP "SOLVE LINEAR SYSTEM DOES NOT WORK"
  ENDIF  
END SUBROUTINE SOLVE_LS_c_sp

SUBROUTINE SOLVE_LS_c_dp(A,X)
  COMPLEX(dp),INTENT(IN):: A(:,:)
  COMPLEX(dp),INTENT(INOUT):: X(:,:)

  !local variables
  INTEGER:: N,M,NB
  INTEGER:: INFO, LWORK
  INTEGER,ALLOCATABLE:: IPIV(:)
  COMPLEX(dp),ALLOCATABLE:: WORK(:)
  COMPLEX(dp):: PREWORK(1)  

  N = SIZE(A(:,1))
  M = SIZE(A(1,:))

  IF(N/=M) THEN
    STOP "MATRIX A IN LINEAR SOLVER IS NOT DIAGONAL"
  ENDIF
  
  NB = SIZE(X(1,:))

  !Allovate pivot array
  ALLOCATE(IPIV(N))

  !Determine optimal size of the array WORK
  !CALL ZHESV("U",N,NB,A,N,IPIV,X,N,PREWORK,-1,INFO)
  CALL ZSYSV("U",N,NB,A,N,IPIV,X,N,PREWORK,-1,INFO)

  LWORK = INT(PREWORK(1))
  ALLOCATE(WORK(LWORK))

  !Solve linear system of equations
  !CALL ZHESV("U",N,NB,A,N,IPIV,X,N,WORK,LWORK,INFO)
  CALL ZSYSV("U",N,NB,A,N,IPIV,X,N,WORK,LWORK,INFO)

  IF(INFO/=0) THEN
    STOP "SOLVE LINEAR SYSTEM DOES NOT WORK"
  ENDIF  
END SUBROUTINE SOLVE_LS_c_dp


!******************************************************************************** 
!       Small function that calculates order of the  permutation 
!       array; Used in DETERMINANT subroutine
!       Permutation as it is defined in LAPACK
!******************************************************************************** 
FUNCTION PERMUTE(IPIV) RESULT(N)
  INTEGER,INTENT(IN):: IPIV(:)
  INTEGER:: N

  !local variables
  INTEGER:: M, i

  M = SIZE(IPIV)
  N = 0 

  DO i=1,M
    IF(i/=IPIV(i)) THEN
      N = N + 1 
    ENDIF
  ENDDO

END FUNCTION PERMUTE


!********************************************************************************
!       Matrix product A*B  for quadratic matrices
!               Input:
!                       A(N,N)
!                       B(N,N)
!                       mode_opt = 1/2/3
!               Output:
!                       C(N,N)=A^T*B
!********************************************************************************
SUBROUTINE MAT_MAT_PROD(N,A,B,C,mode_opt_A,mode_opt_B)
  INTEGER,INTENT(IN):: N 
  REAL(wp),DIMENSION(N,N),INTENT(IN):: A
  REAL(wp),DIMENSION(N,N),INTENT(IN):: B
  REAL(wp),DIMENSION(N,N),INTENT(OUT):: C
  CHARACTER(LEN=1),INTENT(IN),OPTIONAL:: mode_opt_A, mode_opt_B

  !local variables
  REAL(wp):: alpha=1.0_wp , beta=0.0_wp
  CHARACTER(LEN=1):: mode_A , mode_B

  IF(PRESENT(mode_opt_A)) THEN
    mode_A = mode_opt_A
  ELSE
    mode_A = 'N'
  ENDIF

  IF(PRESENT(mode_opt_B)) THEN
    mode_B = mode_opt_B
  ELSE
    mode_B = 'N'
  ENDIF

  IF(wp==sp) THEN
    CALL SGEMM(mode_A , mode_B , N , N , N , alpha , A , N , B , N , beta , C , N)
  ELSEIF(wp==dp) THEN
    CALL DGEMM(mode_A , mode_B , N , N , N , alpha , A , N , B , N , beta , C , N)
  ELSE
    WRITE(*,*) "There is no specific routine for this precision: " , wp
  ENDIF
END SUBROUTINE MAT_MAT_PROD


!********************************************************************************
!       Matrix product A*B  for general matrices
!               Input:
!                       A(M1,N1)
!                       B(M2,N2)
!                       Take care that dimensiona do agree
!               Output:
!                       C - must be allocatable array
!********************************************************************************
SUBROUTINE GEN_MAT_MAT_PROD(M1,N1,A,M2,N2,B,C,mode_opt_A,mode_opt_B)
  INTEGER,INTENT(IN):: M1, N1, M2 , N2 
  REAL(wp),DIMENSION(M1,N1),INTENT(IN):: A
  REAL(wp),DIMENSION(M2,N2),INTENT(IN):: B
  REAL(wp),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: C
  CHARACTER(LEN=1),INTENT(IN),OPTIONAL:: mode_opt_A,mode_opt_B

  !local variables
  REAL(wp):: alpha=1.0_wp , beta=0.0_wp
  CHARACTER(LEN=1):: mode_A , mode_B

  IF(PRESENT(mode_opt_A)) THEN
    mode_A = mode_opt_A
  ELSE
    mode_A = 'N'
  ENDIF

  IF(PRESENT(mode_opt_B)) THEN
    mode_B = mode_opt_B
  ELSE
    mode_B = 'N'
  ENDIF

  IF((mode_A=="N") .and. (mode_B=="N")) THEN
    IF(N1==M2) THEN
      ALLOCATE(C(M1,N2))
      IF(wp==sp) THEN
        CALL SGEMM( mode_A , mode_B , M1 , N2 , N1 , alpha , A , M1 , B , M2 , beta , C , M1)
      ELSEIF(wp==dp) THEN
        CALL DGEMM( mode_A , mode_B , M1 , N2 , N1 , alpha , A , M1 , B , M2 , beta , C , M1)
      ELSE
        WRITE(*,*) "There is no specific LAPACK routine for this precision: " ,wp
      ENDIF
    ELSE
      WRITE(*,*) "Dimension of input matrices do not match"
      CALL EXIT
    ENDIF
  ELSEIF((mode_A=="N") .and. (mode_B=="T")) THEN
    IF(N1==N2) THEN
      ALLOCATE(C(M1,M2))
      IF(wp==sp) THEN
        CALL SGEMM( mode_A , mode_B , M1 , M2 , N1 , alpha , A , M1 , B , M2 , beta , C , M1)
      ELSEIF(wp==dp) THEN
        CALL DGEMM( mode_A , mode_B , M1 , M2 , N1 , alpha , A , M1 , B , M2 , beta , C , M1)
      ELSE
        WRITE(*,*) "There is no specific LAPACK routine for this precision: " ,wp
      ENDIF
    ELSE
      WRITE(*,*) "Dimension of input matrices do not match"
      CALL EXIT
    ENDIF
  ELSEIF((mode_A=="T") .and. (mode_B=="N")) THEN
    IF(M1==M2) THEN
      ALLOCATE(C(N1,N2))
      IF(wp==sp) THEN
        CALL SGEMM( mode_A , mode_B , N1 , N2 , M1 , alpha , A , M1 , B , M2 , beta , C , N1)
      ELSEIF(wp==dp) THEN
        CALL DGEMM( mode_A , mode_B , N1 , N2 , M1 , alpha , A , M1 , B , M2 , beta , C , N1)
      ELSE
        WRITE(*,*) "There is no specific LAPACK routine for this precision: " ,wp
      ENDIF
    ELSE
      WRITE(*,*) "Dimension of input matrices do not match"
      CALL EXIT
    ENDIF
  ELSEIF((mode_A=="T") .and. (mode_B=="T")) THEN
    IF(M1==N2) THEN
      ALLOCATE(C(N1,M2))
      IF(wp==sp) THEN
        CALL SGEMM( mode_A , mode_B , N1 , M2 , M1 , alpha , A , M1 , B , M2 , beta , C , N1)
      ELSEIF(wp==dp) THEN 
        CALL DGEMM( mode_A , mode_B , N1 , M2 , M1 , alpha , A , M1 , B , M2 , beta , C , N1)
      ELSE
        WRITE(*,*) "There is no specific LAPACK routine for this precision: " ,wp
      ENDIF
    ELSE
      WRITE(*,*) "Dimension of input matrices do not match"
      CALL EXIT
    ENDIF
  ELSE
     WRITE(*,*) "Uncrecognised case in GEN_MAT_MAT_PROD"
     CALL EXIT
  ENDIF

END SUBROUTINE GEN_MAT_MAT_PROD


!********************************************************************************
!       Matrix vector multiplication for general matrix A
!               INPUT:
!                       A(N,M)
!                       X(M)
!               OUTPUT:
!                       Y(N)
!********************************************************************************
SUBROUTINE GEN_MAT_VEC_PROD(N,M,A,X,Y)
  INTEGER,INTENT(IN):: N
  INTEGER,INTENT(IN):: M
  REAL(wp),DIMENSION(N,M),INTENT(IN):: A
  REAL(wp),DIMENSION(M),INTENT(IN):: X
  REAL(wp),DIMENSION(N),INTENT(OUT):: Y

  !local variables  
  REAL(wp):: alpha=1.0_wp , beta=0.0_wp

  IF(wp==sp) THEN
    CALL SGEMV('N' , N , M , alpha , A , N , X , 1 , beta , Y , 1)
  ELSEIF(wp==dp) THEN
    CALL DGEMV('N' , N , M , alpha , A , N , X , 1 , beta , Y , 1)
  ELSE
    WRITE(*,*) "There is no specific LAPACK routine for this precision: " ,wp
  ENDIF

END SUBROUTINE GEN_MAT_VEC_PROD


!********************************************************************************
!       Matrix vector multiplication for symmetric matrix A
!               INPUT:
!                       A(N,M)
!                       X(M)
!               OUTPUT:
!                       Y(N)
!********************************************************************************
SUBROUTINE SYM_MAT_VEC_PROD(N,A,X,Y)
  INTEGER,INTENT(IN):: N
  REAL(wp),DIMENSION(N,N),INTENT(IN):: A
  REAL(wp),DIMENSION(N),INTENT(IN):: X
  REAL(wp),DIMENSION(N),INTENT(OUT):: Y

  !local variables  
  REAL(wp):: alpha=1.0_wp , beta=0.0_wp
  
  IF(wp==sp) THEN
    CALL SSYMV('U' , N , alpha , A , N , X , 1 , beta , Y , 1)
  ELSEIF(wp==dp) THEN
    CALL DSYMV('U' , N , alpha , A , N , X , 1 , beta , Y , 1)
  ELSE
    WRITE(*,*) "There is no specific LAPACK routine for this precision: " ,wp
  ENDIF
END SUBROUTINE SYM_MAT_VEC_PROD


!********************************************************************************
!       Eigenvalue solver of symmetric real matrix A - returns only eigenvalues
!               INPUT:
!                       A(N,N)
!               OUTPUT:
!                       EW(N)
!********************************************************************************
SUBROUTINE SYM_EIG(N,A,EW)
  INTEGER,INTENT(IN):: N
  REAL(wp),DIMENSION(N,N),INTENT(IN):: A
  REAL(wp),DIMENSION(N),INTENT(OUT):: EW

  !local variables
  INTEGER:: LWORK
  REAL(wp),DIMENSION(:),ALLOCATABLE:: WORK
  REAL(wp),DIMENSION(N,N):: aux
  INTEGER:: INFO
  
  LWORK = 3*N
  ALLOCATE(WORK(LWORK))

  aux = A

  IF(wp==sp) THEN
    CALL SSYEV('N' , 'U' , N , aux , N , EW , WORK , LWORK , info)
  ELSEIF(wp==dp) THEN
    CALL DSYEV('N' , 'U' , N , aux , N , EW , WORK , LWORK , info)
  ELSE
    WRITE(*,*) "There is no specific LAPACK routine for this precision: " ,wp
  ENDIF

END SUBROUTINE SYM_EIG  


!********************************************************************************
!       Eigenvalue solver of symmetric real matrix A - returns eigenvalues and EV
!               INPUT:
!                       A(N,N)
!               OUTPUT:
!                       EW(N)
!                       EV(N,N)
!********************************************************************************
SUBROUTINE SYM_EIGV(N,A,EW,EV)
  INTEGER,INTENT(IN):: N
  REAL(wp),DIMENSION(N,N),INTENT(IN):: A
  REAL(wp),DIMENSION(N),INTENT(OUT):: EW
  REAL(wp),DIMENSION(N,N),INTENT(OUT):: EV

  !local variables
  INTEGER:: LWORK
  REAL(wp),DIMENSION(:),ALLOCATABLE:: WORK
  INTEGER:: INFO
  
  LWORK = 3*N
  ALLOCATE(WORK(LWORK))

  EV = A
  IF(wp==sp) THEN
    CALL SSYEV('V' , 'U' , N , EV , N , EW , WORK , LWORK , INFO)
  ELSEIF(wp==dp) THEN
    CALL DSYEV('V' , 'U' , N , EV , N , EW , WORK , LWORK , INFO)
  ELSE
    WRITE(*,*) "There is no specific LAPACK routine for this precision: " ,wp
  ENDIF
  
  IF(INFO /= 0 ) THEN 
    WRITE(*,*) 'Termination of program, error in DSYEV LAPACK routine'
    CALL EXIT
  ENDIF
END SUBROUTINE SYM_EIGV  


!********************************************************************************
!       Eigenvalue solver of general real matrix A - returns eigenvalues and 
!       right eigenvectors 
!               INPUT:
!                       A(N,N)
!               OUTPUT:
!                       EW(N)
!                       EV(N,N) - right eigenvectors of A
!                         v(j)  = EV(:,j) + iEV(:,j+1)
!                         v(j+1)= EV(:,j) - iEV(:,j+1)                          
!********************************************************************************
SUBROUTINE GEN_EIGV(N,A,EW,EV)
  INTEGER,INTENT(IN):: N
  REAL(wp),DIMENSION(N,N),INTENT(IN):: A
  REAL(wp),DIMENSION(N,2),INTENT(OUT):: EW
  REAL(wp),DIMENSION(N,N),INTENT(OUT):: EV

  !local variables
  INTEGER:: LWORK
  REAL(wp),DIMENSION(:),ALLOCATABLE:: WORK
  INTEGER:: INFO
 
  LWORK = 4*N
  ALLOCATE(WORK(LWORK))

  IF(wp==sp) THEN
    CALL SGEEV('N' , 'V' , N , A , N , EW(:,1) , EW(:,2), EV , N , EV , N , WORK , LWORK , INFO)
  ELSEIF(wp==dp) THEN
    CALL DGEEV('N' , 'V' , N , A , N , EW(:,1) , EW(:,2), EV , N , EV , N , WORK , LWORK , INFO)
  ELSE
    WRITE(*,*) "There is no specific LAPACK routine for this precision: " ,wp
  ENDIF

  IF(INFO /= 0 ) THEN 
    WRITE(*,*) 'Termination of program, error in DGEEV LAPACK routine'
    CALL EXIT
  ENDIF
END SUBROUTINE GEN_EIGV 


!********************************************************************************
!       Solve linear system of equations Ax = b, where A is N by N symm matrix
!               INPUT:
!                        A(N,N)
!                        b(N)
!               OUTPUT:
!                        x(N)
!********************************************************************************
!SUBROUTINE SOLVE_X(N,A,b,x)
!  INTEGER,INTENT(IN):: N
!  REAL(wp),DIMENSION(N,N),INTENT(IN):: A
!  REAL(wp),DIMENSION(N),INTENT(IN):: b
!  REAL(wp),DIMENSION(N),INTENT(OUT):: x
!
!  !local variables
!  INTEGER:: INFO
!  INTEGER,PARAMETER:: LWORK = N
!  REAL(wp),DIMENSION(LWORK):: WORK
!  REAL(wp),DIMENSION(N,N):: aux
!  REAL(wp),DIMENSION(N):: IPIV
!
!  aux = A
!  CALL DSYSV('U' , N , 1 , aux , N , IPIV , b , N , WORK , LWORK , INFO)
!  IF(INFO /= 0 ) THEN 
!    WRITE(*,*) 'Termination of program, error in DSYSV LAPACK routine'
!    CALL EXIT
!  ENDIF
!END SUBROUTINE SOLVE_X


!********************************************************************************
!       Solve Singular value problem of matrix A   
!                       A = U S V**T                
!               INPUT:                           
!                       A(M,N)   
!               OUTPUT:
!                       S(M,N)   - singular values
!                       U(M,M)   - first N columns of U and VT are 
!                       V(N,N)     left and right singular vectors
!********************************************************************************
SUBROUTINE SIGV(M,N,A,U,S,VT) 
  INTEGER,INTENT(IN):: M
  INTEGER,INTENT(IN):: N
  REAL(wp),DIMENSION(M,N):: A
  REAL(wp),DIMENSION(M,N),INTENT(OUT):: S
  REAL(wp),DIMENSION(M,M),INTENT(OUT):: U
  REAL(wp),DIMENSION(N,N),INTENT(OUT):: VT

  !local variables
  INTEGER:: INFO
  INTEGER:: LWORK
  INTEGER:: IWORK
  REAL(wp),DIMENSION(:),ALLOCATABLE:: WORK
  REAL(wp),DIMENSION(M,N):: aux

  LWORK = 4*(M*N)**2+6*M*N
  IF(M>N) THEN
    LWORK = LWORK + M
    IWORK = 8*N
  ELSE
    LWORK = LWORK + N
    IWORK = 8*M
  ENDIF

  ALLOCATE(WORK(LWORK))

  aux = A

  IF(wp==sp) THEN
    CALL SGESDD('A' , M , N , aux , M , S , U , M , VT , N , WORK , LWORK , IWORK, INFO)
  ELSEIF(wp==dp) THEN
    CALL DGESDD('A' , M , N , aux , M , S , U , M , VT , N , WORK , LWORK , IWORK, INFO)
  ELSE
    WRITE(*,*) "There is no specific LAPACK routine for this precision: " ,wp
  ENDIF

  IF(INFO /= 0 ) THEN 
    WRITE(*,*) 'Termination of program, error in DGESDD LAPACK routine'
    CALL EXIT
  ENDIF
END SUBROUTINE SIGV

END MODULE LAPACK
