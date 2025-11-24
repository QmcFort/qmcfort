! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module gauss_integrals

#include "../preproc.inc"
use constants
use profiling
use qmcfort_pos
use lapack 
use mc_murchie_davidson
use gauss_base
use orth

implicit none 

contains

!******************************************************************************** 
!
! Collect one-electron integrals: 
!
!    overlap matrix        - S_ab
!    kinetic energy matrix - T_ab
!    external potential    - V_ab
!
!******************************************************************************** 
subroutine collect_int_1e(bs, struc, S_ab, T_ab, V_ab, ispin)
  type(gauss_basis_set), intent(in)    :: bs
  type(structure), intent(in)          :: struc
  real(wp), allocatable, intent(inout):: S_ab(:,:,:)
  real(wp), allocatable, intent(inout) :: T_ab(:,:,:)
  real(wp), allocatable, intent(inout) :: V_ab(:,:,:)
  integer, intent(in)                  :: ispin
  !local variables
  integer                              :: a, b
  integer                              :: a_start, a_end
  integer                              :: b_start, b_end
  integer                              :: n, spin
  real(wp), allocatable                :: S_ab_(:,:), T_ab_(:,:), V_ab_(:,:)

  n = bs%get_norbitals()

  if (.not. allocated(S_ab)) allocate(S_ab(n,n,ispin))
  if (.not. allocated(T_ab)) allocate(T_ab(n,n,ispin))
  if (.not. allocated(V_ab)) allocate(V_ab(n,n,ispin))

  do b = 1, size(bs%shell) 
    b_start = bs%get_start_indx(b) + 1
    b_end = bs%get_start_indx(b) + bs%shell(b)%get_shell_size()

    do a = 1, size(bs%shell) !b
      a_start = bs%get_start_indx(a) + 1
      a_end = bs%get_start_indx(a) + bs%shell(a)%get_shell_size()

      call gauss_int_1e(bs%shell(a), bs%shell(b), struc, S_ab_, T_ab_, V_ab_)

      S_ab(a_start:a_end, b_start:b_end, 1) = S_ab_
      T_ab(a_start:a_end, b_start:b_end, 1) = T_ab_
      V_ab(a_start:a_end, b_start:b_end, 1) = V_ab_
    end do
  end do

  !debug:
  !!call trans_symm(S_ab(:,:,1), "u")
  !!call trans_symm(T_ab(:,:,1), "u")
  !!call trans_symm(V_ab(:,:,1), "u")

  do spin = 2, ispin
    S_ab(:,:,spin) = S_ab(:,:,1)
    T_ab(:,:,spin) = T_ab(:,:,1)
    V_ab(:,:,spin) = V_ab(:,:,1)
  end do
end subroutine collect_int_1e


!******************************************************************************** 
!
! Collect two-electron integrals: 
!
!    Coulomb integrals     - g_abcd
!
!******************************************************************************** 
subroutine collect_int_2e(bs, g_abcd, ispin)
  type(gauss_basis_set), intent(in)    :: bs
  real(wp), allocatable, intent(inout) :: g_abcd(:,:,:,:,:)
  integer, intent(in)                  :: ispin
  !local variables
  integer                              :: a, b, c, d
  integer                              :: a_start, a_end
  integer                              :: b_start, b_end
  integer                              :: c_start, c_end
  integer                              :: d_start, d_end
  integer                              :: n, spin
  real(wp), allocatable                :: g_abcd_(:,:,:,:)
  character(len=*), parameter          :: proc_name = "collect_int_2e"

  if (use_profiling) call start_profiling(proc_name)

  n = bs%get_norbitals()

  if (.not. allocated(g_abcd)) allocate(g_abcd(n,n,n,n,ispin))

  !$omp parallel do private(a,a_start,a_end,b,b_start,b_end,c,c_start,c_end,d,d_start,d_end,g_abcd_) schedule(dynamic) collapse(4)
  do d = 1, size(bs%shell)
    do c = 1, size(bs%shell)
      do b = 1, size(bs%shell)
        do a = 1, size(bs%shell)

          if (c > d) cycle
          if (b > d) cycle
          if (a > b) cycle

          a_start = bs%get_start_indx(a) + 1
          a_end = bs%get_start_indx(a) + bs%shell(a)%get_shell_size()
          b_start = bs%get_start_indx(b) + 1
          b_end = bs%get_start_indx(b) + bs%shell(b)%get_shell_size()
          c_start = bs%get_start_indx(c) + 1
          c_end = bs%get_start_indx(c) + bs%shell(c)%get_shell_size()
          d_start = bs%get_start_indx(d) + 1
          d_end = bs%get_start_indx(d) + bs%shell(d)%get_shell_size()
    
          call gauss_int_2e(bs%shell(a), bs%shell(b), bs%shell(c), bs%shell(d), g_abcd_)
    
          g_abcd(a_start:a_end, b_start:b_end, c_start:c_end, d_start:d_end, 1) = g_abcd_
        end do
      end do
    end do
  end do
  !$omp end parallel do
    
  call trans_symm(g_abcd(:,:,:,:,1))

  do spin = 2, ispin
    g_abcd(:,:,:,:,spin) = g_abcd(:,:,:,:,1)
  end do

  if (use_profiling) call end_profiling(proc_name)
end subroutine collect_int_2e


!******************************************************************************** 
!
! Collect two-electron integrals for fixed 2 shells c, and d 
!
!    Coulomb integrals     - g_abcd
!
!******************************************************************************** 
subroutine collect_int_2e_2sh(bs, c, d, V_cd)
  type(gauss_basis_set), intent(in)    :: bs
  integer, intent(in)                  :: c
  integer, intent(in)                  :: d
  real(wp), allocatable, intent(inout) :: V_cd(:,:)
  !local variables
  integer                              :: a, b
  integer                              :: ab_start
  integer                              :: n, n_ab, n_cd
  real(wp), allocatable                :: g_abcd(:,:,:,:)
  character(len=*), parameter          :: proc_name = "collect_int_2e_sh"

  if (profile_code) call start_profiling(proc_name)

  n = bs%get_norbitals()
  n_cd = bs%shell(c)%get_shell_size() * bs%shell(d)%get_shell_size()

  if (allocated(V_cd)) deallocate(V_cd)
  allocate(V_cd(n*n,n_cd))

  !$omp parallel do private(a,n_ab,ab_start,g_abcd) schedule(dynamic) collapse(2)
  do b = 1, size(bs%shell)
    do a = 1, size(bs%shell)
      if (b > a) cycle

      call gauss_int_2e(bs%shell(a), bs%shell(b), bs%shell(c), bs%shell(d), g_abcd)
      n_ab = bs%shell(a)%get_shell_size() * bs%shell(b)%get_shell_size()

      ab_start = bs%get_start_indx(a, b) 
      V_cd(ab_start+1:ab_start+n_ab,:) = reshape(g_abcd, shape=[n_ab, n_cd])

      ab_start = bs%get_start_indx(b, a)
      call transpose_g_abcd(g_abcd)
      V_cd(ab_start+1:ab_start+n_ab,:) = reshape(g_abcd, shape=[n_ab, n_cd])
    end do
  end do
  !$omp end parallel do

  if (profile_code) call end_profiling(proc_name)
contains
  subroutine transpose_g_abcd(g_abcd)
    real(wp), allocatable, intent(inout) :: g_abcd(:,:,:,:)
    !local variables
    integer                              :: c, d
    real(wp), allocatable                :: g_abcd_(:,:,:,:)
    
    if (size(g_abcd,1)==1 .or. size(g_abcd,2)==1) return
    allocate(g_abcd_(size(g_abcd,2), size(g_abcd,1), size(g_abcd,3), size(g_abcd,4)))

    do d = 1, size(g_abcd, 4)
      do c = 1, size(g_abcd, 3)
        g_abcd_(:,:,c,d) = transpose(g_abcd(:,:,c,d))
      end do
    end do

    call move_alloc(g_abcd_, g_abcd)
  end subroutine transpose_g_abcd
end subroutine collect_int_2e_2sh


!******************************************************************************** 
!
! Collect diagonal elemets of two-electron integrals
!
!    Diagonal Coulomb matrix elements - g_abab ==> V_AA
!
!******************************************************************************** 
subroutine collect_int_2e_diag(bs, V_diag)
  type(gauss_basis_set), intent(in)    :: bs
  real(wp), allocatable, intent(inout) :: V_diag(:)
  !local variables
  integer                              :: a, b, ab_start, aa, bb
  integer                              :: n, n_a, n_b, i
  real(wp), allocatable                :: g_abab(:,:,:,:)

  n = bs%get_norbitals()

  if (allocated(V_diag)) deallocate(V_diag)
  allocate(V_diag(n*n))

  !$omp parallel do private(a,aa,bb,n_a,n_b,ab_start,i,g_abab) schedule(dynamic) collapse(2)
  do b = 1, size(bs%shell)
    do a = 1, size(bs%shell)
      n_a = bs%shell(a)%get_shell_size()
      n_b = bs%shell(b)%get_shell_size()
      ab_start = bs%get_start_indx(a,b)

      call gauss_int_2e(bs%shell(a), bs%shell(b), bs%shell(a), bs%shell(b), g_abab)
  
      do bb = 1, n_b
        do aa = 1, n_a
          i = ab_start + (bb-1)*n_a + aa
          V_diag(i) = g_abab(aa,bb,aa,bb)
        end do
      end do
    end do
  end do
  !$omp end parallel do
end subroutine collect_int_2e_diag


!******************************************************************************** 
!
! Calculates one-electron Gaussian integrals for given Gaussian shells
!    (overlap, kinetic energy, and external potential)
!
!    V_ab = (sh_a | sh_b)
!
! output matrix size: (na x nb), where na = 2la+1, nb = 2lb+1
!
!******************************************************************************** 
subroutine gauss_int_1e(shell_a, shell_b, struc, S_ab, T_ab, V_ab)
  type(gauss_shell), intent(in)        :: shell_a
  type(gauss_shell), intent(in)        :: shell_b
  type(structure), intent(in)          :: struc
  real(wp), allocatable, intent(inout) :: S_ab(:,:)
  real(wp), allocatable, intent(inout) :: T_ab(:,:)
  real(wp), allocatable, intent(inout) :: V_ab(:,:)
  !local variables
  integer                            :: va, vb, cont_len_a, cont_len_b, cont_len, indx
  integer                            :: na_c, nb_c, n_c
  real(wp), allocatable              :: S_temp(:,:), T_temp(:,:), V_temp(:,:)
  real(wp), allocatable              :: cont(:), temp(:)
  character(len=*), parameter        :: proc_name = "gauss_int_1e"

  if (use_profiling) call start_profiling(proc_name)
  
  na_c = (shell_a%l+1) * (shell_a%l+2) / 2
  nb_c = (shell_b%l+1) * (shell_b%l+2) / 2
  n_c = na_c * nb_c

  cont_len_a = size(shell_a%cont_coeff)
  cont_len_b = size(shell_b%cont_coeff)
  cont_len = cont_len_a * cont_len_b

  if (allocated(S_ab)) deallocate(S_ab)
  if (allocated(T_ab)) deallocate(T_ab)
  if (allocated(V_ab)) deallocate(V_ab)
  allocate(S_ab(na_c,nb_c), T_ab(na_c,nb_c), V_ab(na_c,nb_c))
  allocate(cont(cont_len), temp(n_c))
  allocate(S_temp(n_c, cont_len))
  allocate(T_temp(n_c, cont_len))
  allocate(V_temp(n_c, cont_len))

  !primitive contractions
  do vb = 1, cont_len_b
    do va = 1, cont_len_a
      indx = (vb-1)*cont_len_a + va
      cont(indx) = shell_a%cont_coeff(va) * shell_b%cont_coeff(vb)
      call mc_murchie_cart_overlap(shell_a%l, shell_b%l, shell_a%cont_exp(va), shell_b%cont_exp(vb), &
                                   shell_a%cite, shell_b%cite, temp)
      S_temp(:,indx) = temp

      call mc_murchie_cart_kin(shell_a%l, shell_b%l, shell_a%cont_exp(va), shell_b%cont_exp(vb), &
                               shell_a%cite, shell_b%cite, temp)
      T_temp(:,indx) = temp

      call mc_murchie_cart_coulomb_1e(shell_a%l, shell_b%l, shell_a%cont_exp(va), shell_b%cont_exp(vb), &
                                      shell_a%cite, shell_b%cite, struc%Zeff, struc%coords, temp)
      V_temp(:,indx) = temp
    end do
  end do
  
  call gemv("n", n_c, cont_len, 1.0_wp, S_temp, n_c, cont, 1, 0.0_wp, temp, 1)
  S_ab = reshape(temp, shape=[na_c,nb_c])

  call gemv("n", n_c, cont_len, 1.0_wp, T_temp, n_c, cont, 1, 0.0_wp, temp, 1)
  T_ab = reshape(temp, shape=[na_c,nb_c])

  call gemv("n", n_c, cont_len, 1.0_wp, V_temp, n_c, cont, 1, 0.0_wp, temp, 1)
  V_ab = reshape(temp, shape=[na_c,nb_c])

  !cartesian - solid harmonics transformation
  call transform_1e_ham(S_ab, shell_a%trafo, shell_b%trafo)
  call transform_1e_ham(T_ab, shell_a%trafo, shell_b%trafo)
  call transform_1e_ham(V_ab, shell_a%trafo, shell_b%trafo)

  if (use_profiling) call end_profiling(proc_name)
end subroutine gauss_int_1e


!******************************************************************************** 
!
! Calculates two-electron Coulomb Gaussian integrals for given Gaussian shells
!
!    g_abcd = (sh_a sh_b | sh_c sh_d)
!
! output array size: (na x nb x nc x nd), where nx = 2lx+1
!
!******************************************************************************** 
subroutine gauss_int_2e(shell_a, shell_b, shell_c, shell_d, g_abcd)
  use standalone, only: tensor_product
  type(gauss_shell), intent(in)      :: shell_a
  type(gauss_shell), intent(in)      :: shell_b
  type(gauss_shell), intent(in)      :: shell_c
  type(gauss_shell), intent(in)      :: shell_d
  real(wp), allocatable, intent(out) :: g_abcd(:,:,:,:)
  !local variables
  integer                            :: va, vb, vc, vd
  integer                            :: na_c, nb_c, nc_c, nd_c, n_c, indx_ab, indx_cd
  integer                            :: cont_len_a, cont_len_b, cont_len_c, cont_len_d, cont_len_ab, cont_len_cd
  real(wp), allocatable              :: g_temp_ab(:,:), g_temp_cd(:,:), temp(:), cont_ab(:), cont_cd(:)
  character(len=*), parameter        :: proc_name = "gauss_int_2e"
  character(len=*), parameter        :: timer_name = "prim_cont_2e"

  if (use_profiling) call start_profiling(proc_name)

  na_c = (shell_a%l+1) * (shell_a%l+2) / 2
  nb_c = (shell_b%l+1) * (shell_b%l+2) / 2
  nc_c = (shell_c%l+1) * (shell_c%l+2) / 2
  nd_c = (shell_d%l+1) * (shell_d%l+2) / 2

  n_c = na_c * nb_c * nc_c * nd_c

  cont_len_a = size(shell_a%cont_coeff)
  cont_len_b = size(shell_b%cont_coeff)
  cont_len_c = size(shell_c%cont_coeff)
  cont_len_d = size(shell_d%cont_coeff)
  cont_len_ab = cont_len_a * cont_len_b
  cont_len_cd = cont_len_c * cont_len_d

  if (allocated(g_abcd)) deallocate(g_abcd)
  allocate(g_abcd(na_c, nb_c, nc_c, nd_c), temp(n_c))
  allocate(g_temp_ab(n_c, cont_len_ab), g_temp_cd(n_c, cont_len_cd))
  allocate(cont_ab(cont_len_ab), cont_cd(cont_len_cd))

  cont_ab = tensor_product(shell_a%cont_coeff, shell_b%cont_coeff)
  cont_cd = tensor_product(shell_c%cont_coeff, shell_d%cont_coeff)

  if (use_profiling) call start_profiling(timer_name)
  do vd = 1, cont_len_d
    do vc = 1, cont_len_c
      
      do vb = 1, cont_len_b
        do va = 1, cont_len_a
          call mc_murchie_cart_coulomb_2e(shell_a%l, shell_b%l, shell_c%l, shell_d%l, &
               shell_a%cont_exp(va), shell_b%cont_exp(vb), shell_c%cont_exp(vc), shell_d%cont_exp(vd), &
               shell_a%cite, shell_b%cite, shell_c%cite, shell_d%cite, temp)
          indx_ab = (vb-1)*cont_len_a + va
          g_temp_ab(:,indx_ab) = temp
        end do
      end do

      indx_cd = (vd-1)*cont_len_c + vc
!debug:
      !call gemv("n", n_c, cont_len_ab, 1.0_wp, g_temp_ab, n_c, cont_ab, 1, 0.0_wp, temp, 1)
      !call gemv("n", n2, m2, fac, F_cd_wyz, n2, R_1d, 1, 0.0_wp, g_cd_tuv(:,indx), 1)
temp = matmul(g_temp_ab, cont_ab)
      g_temp_cd(:,indx_cd) = temp
    end do
  end do

!debug:
  !call gemv("n", n_c, cont_len_cd, 1.0_wp, g_temp_cd, n_c, cont_cd, 1, 0.0_wp, temp, 1)
temp = matmul(g_temp_cd, cont_cd)
  g_abcd = reshape(temp, shape=[na_c, nb_c, nc_c, nd_c])
  if (use_profiling) call end_profiling(timer_name)

  call transform_2e_ham(g_abcd, shell_a%trafo, shell_b%trafo, shell_c%trafo, shell_d%trafo)

  if (use_profiling) call end_profiling(proc_name)
end subroutine gauss_int_2e

end module gauss_integrals