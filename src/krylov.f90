! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module krylov

#include "preproc.inc"
use constants
use profiling
use lapack

implicit none

interface arnoldi
  module procedure arnoldi_r_1d, arnoldi_c_1d, arnoldi_rc_1d, arnoldi_cr_1d
  module procedure arnoldi_r_2d, arnoldi_c_2d, arnoldi_rc_2d, arnoldi_cr_2d
end interface arnoldi

interface block_arnoldi
  module procedure block_arnoldi_r, block_arnoldi_c, block_arnoldi_rc, block_arnoldi_cr
end interface block_arnoldi

contains

!******************************************************************************** 
!
! Arnoldi process to create orthonormal basis B_m = [b1, b2... bm]
! for a Krylov subspace K_m(A, v) = span(v, Av, A^2v ... A^{m-1}v)
!
!     AB_m = B_mH_m + h_{m+1,m}b_{m+1}e_m^{T}
! 
!     H_m - upper Hessenberg matrix 
!
!******************************************************************************** 
subroutine arnoldi_r_1d(A, v, k, B, H)
  real(wp), intent(in)                 :: A(:,:)
  real(wp), intent(in)                 :: v(:)
  integer, intent(in)                  :: k
  real(wp), allocatable, intent(inout) :: B(:,:)
  real(wp), allocatable, intent(inout) :: H(:,:)
  !local variables
  integer                              :: kk, i, j, n, flops
  real(wp), allocatable                :: Av(:)
  character(len=*), parameter          :: proc_name = "arnoldi"

  if (profile_code) call start_profiling(proc_name) 

  n = size(v)
  kk = min(n, k+1)

  if (allocated(B)) deallocate(B)
  if (allocated(H)) deallocate(H)
  allocate(B(n,kk), H(kk,kk), Av(n))

  B = zeror
  H = zeror
  B(:,1) = v / norm2(v)

  do j = 1, kk
    call gemv("n", n, n, oner, A, n, B(:,j), 1, zeror, Av, 1)
    do i = 1, j
       H(i,j) = dot_product(B(:,i), Av)
       Av = Av - H(i,j) * B(:,i)
    end do
    if (j /= kk) then
      H(j+1,j) = norm2(Av)
      B(:,j+1) = Av / H(j+1,j)
    end if
  end do

  flops = 2*n*kk*(n+ (kk+2)/2)
  if (profile_code) call end_profiling(proc_name, m=flops) 
end subroutine arnoldi_r_1d

subroutine arnoldi_c_1d(A, v, k, B, H)
  complex(wp), intent(in)                 :: A(:,:)
  complex(wp), intent(in)                 :: v(:)
  integer, intent(in)                     :: k
  complex(wp), allocatable, intent(inout) :: B(:,:)
  complex(wp), allocatable, intent(inout) :: H(:,:)
  !local variables
  integer                                 :: kk, i, j, n, flops
  complex(wp), allocatable                :: Av(:)
  character(len=*), parameter             :: proc_name = "arnoldi"

  if (profile_code) call start_profiling(proc_name)

  n = size(v)
  kk = min(n, k+1)

  if (allocated(B)) deallocate(B)
  if (allocated(H)) deallocate(H)
  allocate(B(n,kk), H(kk,kk), Av(n))

  B = zeroc
  H = zeroc
  B(:,1) = v / sqrt(real(dot_product(v,v), kind=wp))

  do j = 1, kk
    call gemv("n", n, n, onec, A, n, B(:,j), 1, zeroc, Av, 1)
    do i = 1, j
       H(i,j) = dot_product(B(:,i), Av)
       Av = Av - H(i,j) * B(:,i)
    end do
    if (j /= kk) then
      H(j+1,j) = sqrt(real(dot_product(Av, Av), kind=wp))
      B(:,j+1) = Av / H(j+1,j)
    end if
  end do

  flops = 8*n*kk*(n+ (kk+2)/2)
  if (profile_code) call end_profiling(proc_name, m=flops)
end subroutine arnoldi_c_1d

subroutine arnoldi_rc_1d(A, v, k, B, H)
  real(wp), intent(in)                    :: A(:,:)
  complex(wp), intent(in)                 :: v(:)
  integer, intent(in)                     :: k
  complex(wp), allocatable, intent(inout) :: B(:,:)
  complex(wp), allocatable, intent(inout) :: H(:,:)
  !local variables
  integer                                 :: kk, i, j, n, flops
  complex(wp), allocatable                :: Av(:)
  character(len=*), parameter             :: proc_name = "arnoldi"

  if (profile_code) call start_profiling(proc_name)

  n = size(v)
  kk = min(n, k+1)

  if (allocated(B)) deallocate(B)
  if (allocated(H)) deallocate(H)
  allocate(B(n,kk), H(kk,kk), Av(n))

  B = zeroc
  H = zeroc
  B(:,1) = v / sqrt(real(dot_product(v,v), kind=wp))

  do j = 1, kk
    call gemv("n", n, n, onec, A, n, B(:,j), 1, zeroc, Av, 1)
    do i = 1, j
       H(i,j) = dot_product(B(:,i), Av)
       Av = Av - H(i,j) * B(:,i)
    end do
    if (j /= kk) then
      H(j+1,j) = sqrt(real(dot_product(Av, Av), kind=wp))
      B(:,j+1) = Av / H(j+1,j)
    end if
  end do

  flops = 4*n*n*kk + 4*n*kk*(kk+2)
  if (profile_code) call end_profiling(proc_name)
end subroutine arnoldi_rc_1d

subroutine arnoldi_cr_1d(A, v, k, B, H)
  complex(wp), intent(in)                 :: A(:,:)
  real(wp), intent(in)                    :: v(:)
  integer, intent(in)                     :: k
  complex(wp), allocatable, intent(inout) :: B(:,:)
  complex(wp), allocatable, intent(inout) :: H(:,:)
  !local variables
  complex(wp), allocatable                :: v_(:)

  allocate(v_(size(v)))
  v_ = cmplx(v, 0.0_wp, kind=wp)

  call arnoldi(A, v_, k, B, H)
end subroutine arnoldi_cr_1d

subroutine arnoldi_r_2d(A, V, k, B, H)
  real(wp), intent(in)                 :: A(:,:)
  real(wp), intent(in)                 :: V(:,:)
  integer, intent(in)                  :: k
  real(wp), allocatable, intent(inout) :: B(:,:,:)
  real(wp), allocatable, intent(inout) :: H(:,:,:)
  !local variables
  integer                              :: kk, i, j, n, m, p, flops
  real(wp), allocatable                :: Av(:,:)
  character(len=*), parameter          :: proc_name = "arnoldi"

  if (profile_code) call start_profiling(proc_name)

  n = size(V, 1)
  m = size(V, 2)
  kk = min(n, k+1)

  if (allocated(B)) deallocate(B)
  if (allocated(H)) deallocate(H)
  allocate(B(n,m,kk), H(m,kk,kk), Av(n,m))

  B = zeror
  H = zeror

  do i = 1, m
    B(:,i,1) = V(:,i) / norm2(V(:,i))
  end do

  do j = 1, kk
    if (profile_code) call start_profiling("arnoldi_gemm")
    call gemm("n", "n", n, m, n, oner, A, n, B(:,:,j), n, zeror, Av, n)
    if (profile_code) call end_profiling("arnoldi_gemm", n, m, n, "r")
    do i = 1, j
      do p = 1, m
        H(p,i,j) = dot_product(B(:,p,i), Av(:,p))
        Av(:,p) = Av(:,p) - H(p,i,j) * B(:,p,i)
      end do
    end do
    if (j /= kk) then
      do p = 1, m
        H(p,j+1,j) = norm2(Av(:,p))
        B(:,p,j+1) = Av(:,p) / H(p,j+1,j)
      end do
    end if
  end do

  flops = 2*n*m*kk*(n+ (kk+2)/2)
  if (profile_code) call end_profiling(proc_name, m=flops)
end subroutine arnoldi_r_2d

subroutine arnoldi_c_2d(A, V, k, B, H)
  complex(wp), intent(in)                 :: A(:,:)
  complex(wp), intent(in)                 :: V(:,:)
  integer, intent(in)                     :: k
  complex(wp), allocatable, intent(inout) :: B(:,:,:)
  complex(wp), allocatable, intent(inout) :: H(:,:,:)
  !local variables
  integer                                 :: kk, i, j, n, m, p, flops
  complex(wp), allocatable                :: Av(:,:)
  character(len=*), parameter             :: proc_name = "arnoldi"

  if (profile_code) call start_profiling(proc_name)

  n = size(V, 1)
  m = size(V, 2)
  kk = min(n, k+1)

  if (allocated(B)) deallocate(B)
  if (allocated(H)) deallocate(H)
  allocate(B(n,m,kk), H(m,kk,kk), Av(n,m))

  B = zeroc
  H = zeroc

  do i = 1, m
    B(:,i,1) = V(:,i) / sqrt(real(dot_product(V(:,i),V(:,i)), kind=wp))
  end do

  do j = 1, kk
    if (profile_code) call start_profiling("arnoldi_gemm")
    call gemm("n", "n", n, m, n, onec, A, n, B(:,:,j), n, zeroc, Av, n)
    if (profile_code) call end_profiling("arnoldi_gemm", n, m, n, "c")
    do i = 1, j
      do p = 1, m
        H(p,i,j) = dot_product(B(:,p,i), Av(:,p))
        Av(:,p) = Av(:,p) - H(p,i,j) * B(:,p,i)
      end do
    end do
    if (j /= kk) then
      do p = 1, m
        H(p,j+1,j) = sqrt(real(dot_product(Av(:,p), Av(:,p)), kind=wp))
        B(:,p,j+1) = Av(:,p) / H(p,j+1,j)
      end do
    end if
  end do

  flops = 8*n*m*kk*(n+ (kk+2)/2)
  if (profile_code) call end_profiling(proc_name, m=flops)
end subroutine arnoldi_c_2d

subroutine arnoldi_rc_2d(A, V, k, B, H)
  real(wp), intent(in)                    :: A(:,:)
  complex(wp), intent(in)                 :: V(:,:)
  integer, intent(in)                     :: k
  complex(wp), allocatable, intent(inout) :: B(:,:,:)
  complex(wp), allocatable, intent(inout) :: H(:,:,:)
  !local variables
  integer                                 :: kk, i, j, n, m, p, flops
  complex(wp), allocatable                :: Av(:,:)
  character(len=*), parameter             :: proc_name = "arnoldi"

  if (profile_code) call start_profiling(proc_name)

  n = size(V, 1)
  m = size(V, 2)
  kk = min(n, k+1)

  if (allocated(B)) deallocate(B)
  if (allocated(H)) deallocate(H)
  allocate(B(n,m,kk), H(m,kk,kk), Av(n,m))

  B = zeroc
  H = zeroc

  do i = 1, m
    B(:,i,1) = V(:,i) / sqrt(real(dot_product(V(:,i),V(:,i)), kind=wp))
  end do

  do j = 1, kk
    call start_profiling("arnoldi")
    call gemm("n", "n", n, m, n, onec, A, n, B(:,:,j), n, zeroc, Av, n)
    call end_profiling("arnoldi", n, m, n, "rc")
    do i = 1, j
      do p = 1, m
        H(p,i,j) = dot_product(B(:,p,i), Av(:,p))
        Av(:,p) = Av(:,p) - H(p,i,j) * B(:,p,i)
      end do
    end do
    if (j /= kk) then
      do p = 1, m
        H(p,j+1,j) = sqrt(real(dot_product(Av(:,p), Av(:,p)), kind=wp))
        B(:,p,j+1) = Av(:,p) / H(p,j+1,j)
      end do
    end if
  end do

  flops = 4*n*n*m*kk + 4*n*m*kk*(kk+2)
  if (profile_code) call end_profiling(proc_name, m=flops)
end subroutine arnoldi_rc_2d

subroutine arnoldi_cr_2d(A, V, k, B, H)
  complex(wp), intent(in)                 :: A(:,:)
  real(wp), intent(in)                    :: V(:,:)
  integer, intent(in)                     :: k
  complex(wp), allocatable, intent(inout) :: B(:,:,:)
  complex(wp), allocatable, intent(inout) :: H(:,:,:)
  !local variables
  complex(wp), allocatable                :: V_(:,:)

  allocate(V_(size(V,1), size(V,2)))
  V_ = cmplx(V, 0.0_wp, kind=wp)

  call arnoldi(A, V_, k, B, H)
end subroutine arnoldi_cr_2d

subroutine block_arnoldi_r(A, V, k, B, H)
  real(wp), intent(in)                 :: A(:,:)
  real(wp), intent(in)                 :: V(:,:)
  integer, intent(in)                  :: k
  real(wp), allocatable, intent(inout) :: B(:,:)
  real(wp), allocatable, intent(inout) :: H(:,:)
  !local variables
  integer                              :: kk, i, i1, i2, j, j1, j2, n, m, p, flops
  real(wp), allocatable                :: Av(:,:), R(:,:)
  character(len=*), parameter          :: proc_name = "block_arnoldi"

  if (profile_code) call start_profiling(proc_name)

  n = size(V, 1)
  m = size(V, 2)
  kk = min(n/m, k+1)

  if (allocated(B)) deallocate(B)
  if (allocated(H)) deallocate(H)
  allocate(B(n,m*kk), H(m*kk,m*kk), Av(n,m), R(m,m))

  B = zeror
  H = zeror

  B(:,1:m) = V
  call getqr(B(:,1:m), R)

  do j = 1, kk
    j1 = (j-1)*m + 1
    j2 = j*m
    if (profile_code) call start_profiling("arnoldi_gemm")
    call gemm("n", "n", n, m, n, oner, A, n, B(:,j1:j2), n, zeror, Av, n)
    if (profile_code) call end_profiling("arnoldi_gemm", n, m, n, "r")
    do i = 1, j
      i1 = (i-1)*m + 1
      i2 = i*m
      if (profile_code) call start_profiling("arnoldi_orth")
      call gemm("t", "n", m, m, n, oner, B(:,i1:i2), n, Av, n, zeror, H(i1:i2,j1:j2), m)
      call gemm("n", "n", n, m, m, -oner, B(:,i1:i2), n, H(i1:i2,j1:j2), m, oner, Av, n)
      if (profile_code) call end_profiling("arnoldi_orth", m=4*n*m*m)
      !calculate   V(:,i1:i1)^{\dagger} Av = H(i1:i2,j1:j2)
      !update      Av = Av - V(:,i1:i2)H(i1:i2,j1:j2)
    end do
    if (j /= kk) then
      i1 = j*m + 1
      i2 = (j+1) * m
      call getqr(Av, H(i1:i2,j1:j2))
      B(:,i1:i2) = Av
    end if
  end do

  flops = 2*n*m*kk*(n+ 2*m) + 2*(n-m/3)*m*m*kk
  if (profile_code) call end_profiling(proc_name, m=flops)
end subroutine block_arnoldi_r

subroutine block_arnoldi_c(A, V, k, B, H)
  complex(wp), intent(in)                 :: A(:,:)
  complex(wp), intent(in)                 :: V(:,:)
  integer, intent(in)                     :: k
  complex(wp), allocatable, intent(inout) :: B(:,:)
  complex(wp), allocatable, intent(inout) :: H(:,:)
  !local variables
  integer                                 :: kk, i, i1, i2, j, j1, j2, n, m, p, flops
  complex(wp), allocatable                :: Av(:,:), R(:,:)
  character(len=*), parameter             :: proc_name = "block_arnoldi"

  if (profile_code) call start_profiling(proc_name)

  n = size(V, 1)
  m = size(V, 2)
  kk = min(n/m, k+1)

  if (allocated(B)) deallocate(B)
  if (allocated(H)) deallocate(H)
  allocate(B(n,m*kk), H(m*kk,m*kk), Av(n,m), R(m,m))

  B = zeror
  H = zeror

  B(:,1:m) = V
  call getqr(B(:,1:m), R)

  do j = 1, kk
    j1 = (j-1)*m + 1
    j2 = j*m
    if (profile_code) call start_profiling("arnoldi_gemm")
    call gemm("n", "n", n, m, n, onec, A, n, B(:,j1:j2), n, zeroc, Av, n)
    if (profile_code) call end_profiling("arnoldi_gemm", n, m, n, "c")
    do i = 1, j
      i1 = (i-1)*m + 1
      i2 = i*m
      if (profile_code) call start_profiling("arnoldi_orth")
      call gemm("c", "n", m, m, n, onec, B(:,i1:i2), n, Av, n, zeroc, H(i1:i2,j1:j2), m)
      call gemm("n", "n", n, m, m, -onec, B(:,i1:i2), n, H(i1:i2,j1:j2), m, onec, Av, n)
      if (profile_code) call end_profiling("arnoldi_orth", m=16*n*m*m)
      !calculate   V(:,i1:i1)^{\dagger} Av = H(i1:i2,j1:j2)
      !update      Av = Av - V(:,i1:i2)H(i1:i2,j1:j2)
    end do
    if (j /= kk) then
      i1 = j*m + 1
      i2 = (j+1) * m
      call getqr(Av, H(i1:i2,j1:j2))
      B(:,i1:i2) = Av
    end if
  end do

  flops = 8*n*m*kk*(n+ 2*m) + 8*(n-m/3)*m*m*kk
  if (profile_code) call end_profiling(proc_name, m=flops)
end subroutine block_arnoldi_c

subroutine block_arnoldi_rc(A, V, k, B, H)
  real(wp), intent(in)                    :: A(:,:)
  complex(wp), intent(in)                 :: V(:,:)
  integer, intent(in)                     :: k
  complex(wp), allocatable, intent(inout) :: B(:,:)
  complex(wp), allocatable, intent(inout) :: H(:,:)
  !local variables
  integer                                 :: kk, i, i1, i2, j, j1, j2, n, m, p, flops
  complex(wp), allocatable                :: Av(:,:), R(:,:)
  character(len=*), parameter             :: proc_name = "block_arnoldi"

  if (profile_code) call start_profiling(proc_name)

  n = size(V, 1)
  m = size(V, 2)
  kk = min(n/m, k+1)

  if (allocated(B)) deallocate(B)
  if (allocated(H)) deallocate(H)
  allocate(B(n,m*kk), H(m*kk,m*kk), Av(n,m), R(m,m))

  B = zeror
  H = zeror

  B(:,1:m) = V
  call getqr(B(:,1:m), R)

  do j = 1, kk
    j1 = (j-1)*m + 1
    j2 = j*m
    if (profile_code) call start_profiling("arnoldi_gemm")
    call gemm("n", "n", n, m, n, onec, A, n, B(:,j1:j2), n, zeroc, Av, n)
    if (profile_code) call end_profiling("arnoldi_gemm", n, m, n, "rc")
    do i = 1, j
      i1 = (i-1)*m + 1
      i2 = i*m
      if (profile_code) call start_profiling("arnoldi_orth")
      call gemm("c", "n", m, m, n, onec, B(:,i1:i2), n, Av, n, zeroc, H(i1:i2,j1:j2), m)
      call gemm("n", "n", n, m, m, -onec, B(:,i1:i2), n, H(i1:i2,j1:j2), m, onec, Av, n)
      if (profile_code) call end_profiling("arnoldi_orth", m=16*n*m*m)
      !calculate   V(:,i1:i1)^{\dagger} Av = H(i1:i2,j1:j2)
      !update      Av = Av - V(:,i1:i2)H(i1:i2,j1:j2)
    end do
    if (j /= kk) then
      i1 = j*m + 1
      i2 = (j+1) * m
      call getqr(Av, H(i1:i2,j1:j2))
      B(:,i1:i2) = Av
    end if
  end do

  flops = 4*n*m*kk*(n+ 4*m) + 8*(n-m/3)*m*m*kk
  if (profile_code) call end_profiling(proc_name, m=flops)
end subroutine block_arnoldi_rc

subroutine block_arnoldi_cr(A, V, k, B, H)
  complex(wp), intent(in)                 :: A(:,:)
  real(wp), intent(in)                    :: V(:,:)
  integer, intent(in)                     :: k
  complex(wp), allocatable, intent(inout) :: B(:,:)
  complex(wp), allocatable, intent(inout) :: H(:,:)
  !local variables
  complex(wp), allocatable                :: V_(:,:)

  allocate(V_(size(V,1), size(V,2)))
  V_ = cmplx(V, 0.0_wp, kind=wp)

  call block_arnoldi(A, V_, k, B, H)
end subroutine block_arnoldi_cr

end module krylov