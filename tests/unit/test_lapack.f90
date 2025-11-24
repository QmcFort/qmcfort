! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module test_lapack

use mpi

use fruit, only: assert_equals, run_test_case, test_set_summary
use lcg48_basic_rng, only: Lcg48BasicRng

implicit none

private 
public :: test_lapack_set

contains

subroutine test_lapack_set
  call run_test_case(test_gemv_rc, "real complex matrix vector multiplication check")
  call run_test_case(test_svd_none_trivial, "singular values of the small symmetric matrix")
  call run_test_case(test_svd_all_consist, "check consistency of the svd decomposition a = u*sigma*vt")
  call run_test_case(test_svd_all_none, "check consistency of the svd_all and svd_none")
  call run_test_case(test_gemm_rc_nn_20, "gemm_rc test with ta=n, tb=n, alpha=2 beta=0")
  call run_test_case(test_gemm_rc_nc_11, "gemm_rc test with ta=n, tb=c, alpha=1 beta=1")
  call run_test_case(test_gemm_cr_cons, "gemm_cr test with ta=c, tb=t, alpha=1 beta=1")
  call run_test_case(test_cgeev, "test cgeev routine comparing to cheev")
  call run_test_case(test_zgeev_trivial, "trivial test of the zgeev routine")
  call run_test_case(test_zgeev_cons, "test consistency of the zgeev routine")
  call run_test_case(test_qr_cons, "test consistency of qr decomposition")
  call run_test_case(test_qr_det_cons, "test consistency of qr_det routine")

  call test_set_summary("src/lapack.f90")
end subroutine test_lapack_set

subroutine test_gemv_rc
  use constants, only: wp
  use lapack, only: gemv
  use standalone, only: imag
  integer, parameter            :: n=50, prec=wp
  real(prec), parameter         :: tol = 1.0e-04_prec
  real(prec), allocatable       :: a(:,:), br(:), bi(:), cr(:), ci(:)
  complex(prec)                 :: alpha, beta
  complex(prec), allocatable    :: b(:), c(:), c1(:)
  type(Lcg48BasicRng)           :: brng

  if (.not. comm_world%initialized) call comm_world%init()
  brng = Lcg48BasicRng()

  alpha = (2.0_prec, 0.0_prec)
  beta = (0.0_prec, 1.0_prec)

  allocate(a(n,n), br(n), bi(n), b(n), cr(n), ci(n), c(n), c1(n))
  call brng%rand(a)
  call brng%rand(br)
  call brng%rand(bi)
  call brng%rand(cr)
  call brng%rand(ci)
  b = cmplx(br, bi, kind=prec)
  c = cmplx(cr, ci, kind=prec)
  c1 = c

  c = alpha*matmul(a, b) + beta*c
  call gemv("n", n, n, alpha, A, n, b, 1, beta, c1, 1)

  call assert_equals(c, c1, n, tol)
end subroutine test_gemv_rc

subroutine test_svd_none_trivial
  use lapack, only: svd
  use constants, only: sp
  integer, parameter     :: n = 3
  real(sp), parameter :: tol = 1.0e-05_sp
  real(sp)            :: matrix(n,n)
  real(sp)            :: expected(n)
  real(sp)            :: result_(n)
  
  expected = [2.66907909, 2.14510269, 0.5239764]
  matrix = reshape([1.0, -1.0, 0.0, -1.0, 0.0, 2.0, 0.0, 2.0, -1.0], shape=[n, n])
  
  call svd(matrix, result_)
  call assert_equals(expected, result_, n, tol)
end subroutine test_svd_none_trivial

subroutine test_svd_all_consist
  use lapack, only: svd
  use constants, only: wp
  integer, parameter            :: prec=wp, n=20, m=10
  real(prec), parameter         :: tol = 1.0e-04_prec
  integer                       :: i
  real(prec)                    :: a(n,m), u(n,n), vt(m,m), sigma(n), sigma_(n,m)
  real(prec)                    :: expected(n,m), result_(n,m)
  type(Lcg48BasicRng)           :: brng

  if (.not. comm_world%initialized) call comm_world%init()
  brng = Lcg48BasicRng()

  call brng%rand(a)
  expected = a 
  
  call svd(a, sigma, u, vt)
  
  sigma_ = 0.0_prec
  do i = 1, min(n, m)
    sigma_(i,i) = sigma(i) 
  end do

  result_ = matmul(matmul(u,sigma_), vt)

  call assert_equals(expected, result_, n, m, tol)
end subroutine test_svd_all_consist

subroutine test_svd_all_none
  use lapack, only: svd
  use constants, only: wp
  integer, parameter            :: prec=wp, n=20, m=10
  real(prec), parameter         :: tol = 1.0e-04_prec
  integer                       :: i
  real(prec)                    :: a(n,m), a2(n,m), vt(m,m), sigma(n)
  real(prec)                    :: expected(n,m), result_(n,n)
  type(Lcg48BasicRng)           :: brng

  if (.not. comm_world%initialized) call comm_world%init()
  brng = Lcg48BasicRng()

  if (.not. comm_world%initialized) call comm_world%init
  
  call brng%rand(a)
  expected = a 
  
  call svd(a, sigma, result_, vt)
  call svd(expected, sigma)
  
  call assert_equals(abs(expected), abs(result_(:,1:m)), n, m, tol)
end subroutine test_svd_all_none

subroutine test_gemm_rc_nn_20
  use constants, only: dp
  use lapack, only: gemm
  real(dp), parameter      :: tol=1.0e-04_dp
  integer                  :: m=2, n=2, k=3 
  complex(dp)              :: alpha=(2.0_dp, 0.0_dp), beta=(0.0_dp, 0.0_dp)
  real(dp), allocatable    :: a(:,:), b_(:,:), d(:,:)
  complex(dp), allocatable :: b(:,:), c(:,:), c_(:,:)
  
  allocate(a(m,k), b(k,n), c(m,n), c_(m,n), b_(2*n,k), d(2*n,m))
  
  a = 0.0_dp
  a(:,1) = [1.0_dp, 3.0_dp]
  a(:,2) = [0.0_dp, -1.0_dp]
  a(:,3) = [1.0_dp, 0.0_dp]
  
  b = (0.0_dp, 0.0_dp)
  b(:,1) = [(1.0_dp, 1.0_dp), (1.0_dp, -1.0_dp), (0.0_dp, 0.0_dp)]
  b(:,2) = [(0.0_dp, 0.0_dp), (1.0_dp, 1.0_dp), (2.0_dp, -1.0_dp)]
  
  c = (0.0_dp, 0.0_dp)
  c_ = (0.0_dp, 0.0_dp)
  c_(:,1) = [(2.0_dp, 2.0_dp), (4.0_dp, 8.0_dp)]
  c_(:,2) = [(4.0_dp, -2.0_dp), (-2.0_dp, -2.0_dp)]
  
  call gemm("n", "n", m, n, k, alpha, a, m ,b, k, beta, c, m)
  
  call assert_equals(c_, c, m, n, tol)
end subroutine test_gemm_rc_nn_20

subroutine test_gemm_rc_nc_11
  use constants, only: wp
  use lapack, only: gemm
  integer, parameter            :: prec = wp
  real(prec), parameter         :: tol=1.0e-04_prec
  integer                       :: m=10, n=6, k=8 
  complex(prec)                 :: alpha=(1.0_prec, 0.0_prec), beta=(1.0_prec, 0.0_prec)
  real(prec), allocatable       :: a(:,:), b_r(:,:), b_i(:,:), c_r(:,:), c_i(:,:)
  complex(prec), allocatable    :: b(:,:), c(:,:), c_(:,:)
  type(Lcg48BasicRng)           :: brng

  if (.not. comm_world%initialized) call comm_world%init()
  brng = Lcg48BasicRng()

  allocate(a(m,k), b(n,k), c(m,n), c_(m,n))
  allocate(b_r(n,k), b_i(n,k), c_r(m,n), c_i(m,n))
  
  call brng%rand(a)
  call brng%rand(b_r)
  call brng%rand(b_i)
  call brng%rand(c_r)
  call brng%rand(c_i)

  b = cmplx(b_r, b_i, kind=prec)
  c = cmplx(c_r, c_i, kind=prec)
  c_ = c
  
  c_ = beta*c_ + alpha*matmul(a, transpose(b))
  call gemm("n", "t", m, n, k, alpha, a, m ,b, n, beta, c, m)
  
  call assert_equals(c_, c, m, n, tol)
end subroutine test_gemm_rc_nc_11

subroutine test_gemm_cr_cons
  use constants, only: wp
  use lapack, only: gemm
  integer, parameter            :: prec=wp
  real(prec), parameter         :: tol=1.0e-04_prec
  integer                       :: m=10, n=8, k=12
  complex(prec)                 :: alpha=(1.0_prec, 0.0_prec), beta=(1.0_prec, 0.0_prec)
  real(prec), allocatable       :: a_r(:,:), a_i(:,:), b(:,:), c_r(:,:), c_i(:,:)
  complex(prec), allocatable    :: a(:,:), c(:,:), c_(:,:)
  type(Lcg48BasicRng)           :: brng

  if (.not. comm_world%initialized) call comm_world%init()
  brng = Lcg48BasicRng()

  allocate(a(k,m), b(n,k), c(m,n), c_(m,n))
  allocate(a_r(k,m), a_i(k,m), c_r(m,n), c_i(m,n))
  
  call brng%rand(a_r)
  call brng%rand(a_i)
  call brng%rand(b)
  call brng%rand(c_r)
  call brng%rand(c_i)

  a = cmplx(a_r, a_i, kind=prec)
  c = cmplx(c_r, c_i, kind=prec)
  c_ = c
  
  c_ = beta*c_ + alpha*matmul(conjg(transpose(a)), transpose(b))
  call gemm("c", "t", m, n, k, alpha, a, k, b, n, beta, c, m)
  
  call assert_equals(c_, c, m, n, tol)
end subroutine test_gemm_cr_cons

subroutine test_cgeev
  use constants, only: wp
  use lapack, only: syev, geev
  integer, parameter            :: n=5, prec=wp
  real(prec), parameter         :: tol=1.0e-04_prec
  integer                       :: i, j
  complex(prec), allocatable    :: a(:,:), a_(:,:)
  real(prec), allocatable       :: ev_(:), evr(:)
  complex(prec), allocatable    :: ev__(:), U_(:,:), U__(:,:)
  real(prec), allocatable       :: ar(:,:), ai(:,:)
  real(prec)                    :: swap
  type(Lcg48BasicRng)           :: brng

  if (.not. comm_world%initialized) call comm_world%init()
  brng = Lcg48BasicRng()

  allocate(a(n,n), ar(n,n), ai(n,n))
  allocate(ev_(n), ev__(n), U_(n,n), U__(n,n))
  call brng%rand(ar)
  call brng%rand(ai)

  ar = (ar + transpose(ar)) / 2.0_sp
  ai = (ai - transpose(ai)) / 2.0_sp
  a = cmplx(ar, ai, kind=prec)

  allocate(a_, source=a)
  
  call syev(a_, ev_)
  call geev(a, ev__, U__)

  evr = real(ev__, prec)
  do i = 1, size(evr)
    do j = i+1, size(evr)
      if (evr(i) > evr(j)) then
        swap = evr(i)
        evr(i) = evr(j)
        evr(j) = swap
      end if
    end do

  end do

  call assert_equals(ev_, evr, n, tol)
end subroutine test_cgeev

subroutine test_zgeev_trivial
  use constants, only: wp
  use lapack, only: geev, inverse
  integer, parameter         :: prec=wp, n=3
  real(prec), parameter      :: tol=1.0e-04_prec
  integer                    :: i
  complex(prec), allocatable :: a(:,:)
  complex(prec), allocatable :: ev(:), ev_(:)
  complex(prec), allocatable :: U(:,:), U_(:,:), Dev(:,:)

  allocate(a(n,n), ev(n), ev_(n), U(n,n), U_(n,n))

  a(:,1) = [(1.0_prec, 1.0_prec), (2.0_prec, -1.0_prec), (0.0_prec, 1.0_prec)]
  a(:,2) = [(1.0_prec, -1.0_prec), (2.0_prec, 0.0_prec), (0.0_prec, 0.0_prec)]
  a(:,3) = [(0.0_prec, 0.0_prec), (1.0_prec, -1.0_prec), (3.0_prec, -1.0_prec)]


  call geev(a, ev, U, jobz="V")
  U_ = U
  call inverse(U_)
  Dev = matmul(U_, matmul(a, U))
  
  do i = 1, n
    ev_(i) = Dev(i,i)
  end do

  call assert_equals(ev, ev_, n, tol)
end subroutine test_zgeev_trivial

subroutine test_zgeev_cons
  use constants, only: wp
  use lapack, only: geev, inverse
  integer, parameter            :: prec=wp, n=10
  real(prec), parameter         :: tol=1.0e-04_prec
  integer                       :: i
  complex(prec), allocatable    :: a(:,:)
  complex(prec), allocatable    :: ev(:), ev_(:)
  complex(prec), allocatable    :: U(:,:), U_(:,:), Dev(:,:)
  real(prec), allocatable       :: ar(:,:), ai(:,:)
  type(Lcg48BasicRng)           :: brng

  if (.not. comm_world%initialized) call comm_world%init()
  brng = Lcg48BasicRng()

  allocate(a(n,n), ar(n,n), ai(n,n))
  allocate(ev(n), ev_(n), U(n,n), U_(n,n), Dev(n,n))

  call brng%rand(ar)
  call brng%rand(ai)
  a = cmplx(ar, ai, kind=dp)

  call geev(a, ev, U, jobz="V")
  U_ = U
  call inverse(U_)

  Dev = matmul(U_, matmul(a, U))
  
  do i = 1, n
    ev_(i) = Dev(i,i)
  end do

  call assert_equals(ev, ev_, n, tol)
end subroutine test_zgeev_cons

subroutine test_qr_cons
  use constants, only: wp
  use lapack, only: getqr
  integer, parameter            :: n=10, m=4, prec=wp
  real(prec), parameter         :: tol=1.0e-04_prec
  integer                       :: i
  complex(prec), allocatable    :: A(:,:), A_(:,:), Q(:,:), R(:,:)
  real(prec), allocatable       :: Ar(:,:), Ai(:,:)
  type(Lcg48BasicRng)           :: brng

  if (.not. comm_world%initialized) call comm_world%init()
  brng = Lcg48BasicRng()

  allocate(A(n,m), A_(n,m), Q(n,m), R(m,m))
  allocate(Ar(n,m), Ai(n,m))

  call brng%rand(Ar)
  call brng%rand(Ai)
  A = cmplx(Ar, Ai, kind=prec)
  A = 0.005_prec * A
  do i = 1, m
    A(i,i) = (1.0_prec, 0.0_prec)
  end do

  Q = A
  call getqr(Q, R)
  A_ = matmul(Q, R)

  call assert_equals(A, A_, n, m, tol)
end subroutine test_qr_cons

subroutine test_qr_det_cons
  use constants, only: wp
  use lapack, only: getqr, getqr_det
  integer, parameter            :: n=10, m=4, prec=wp
  real(prec), parameter         :: tol=1.0e-04_prec
  integer                       :: i
  real(prec), allocatable       :: A(:,:), Q(:,:), R(:,:)
  real(prec)                    :: det, det_
  type(Lcg48BasicRng)           :: brng

  if (.not. comm_world%initialized) call comm_world%init()
  brng = Lcg48BasicRng()

  allocate(A(n,m), Q(n,m), R(m,m))

  call brng%rand(A)
  A = 0.005_prec * A
  do i = 1, m
    A(i,i) = 1.0_prec
  end do

  Q = A
  call getqr(Q, R)
  det = 1.0_prec
  do i = 1, m
    det = det * R(i,i)
  end do

  Q = A
  call getqr_det(Q, det_)

  call assert_equals(det, det_, tol)
end subroutine test_qr_det_cons

end module test_lapack