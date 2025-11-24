! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module mc_murchie_davidson

#include "../preproc.inc"
use constants
use profiling
use lapack
use boys_functions
use gauss_base

implicit none

real(wp), parameter :: hfac_ = 2.0_wp * pi**(5.0_wp/2.0_wp)

contains

!******************************************************************************** 
!
! Computes cartesian overlap integrals within shell
!
!    S_ab = (pi/p)^{3/2} E_{000}^{ab}
!
! with
!
!    p = a+b
!
!******************************************************************************** 
subroutine mc_murchie_cart_overlap(la, lb, a, b, Ra, Rb, S_ab)
  integer, intent(in)                  :: la
  integer, intent(in)                  :: lb
  real(wp), intent(in)                 :: a
  real(wp), intent(in)                 :: b
  real(wp), intent(in)                 :: Ra(3)
  real(wp), intent(in)                 :: Rb(3)
  real(wp), allocatable, intent(inout) :: S_ab(:)
  !local variables
  integer                              :: na, nb, a_indx, b_indx, i, j, k, l, m, n, indx
  real(wp), allocatable                :: E_abt(:,:,:), E_abu(:,:,:), E_abv(:,:,:)
  character(len=*), parameter          :: proc_name = "mc_murchie_overlap"

  if (profile_code) call start_profiling(proc_name)

  na = (la+1) * (la+2) / 2
  nb = (lb+1) * (lb+2) / 2

  call overlap_exp_coeff(la, lb, a, b, Ra(1), Rb(1), E_abt)
  call overlap_exp_coeff(la, lb, a, b, Ra(2), Rb(2), E_abu)
  call overlap_exp_coeff(la, lb, a, b, Ra(3), Rb(3), E_abv)

  call null_overlap_exp_coeff(E_abt, E_abu, E_abv, S_ab)
  S_ab = S_ab * (pi/(a+b))**(1.5_wp)

  if (profile_code) call end_profiling(proc_name)
end subroutine mc_murchie_cart_overlap


!******************************************************************************** 
!
! Computes cartesian kinetic integrals within shell
!
!    T_ab =  T_ij E_{0}^{kl} E_{0}^{mn} + 
!            E_{0}^{ij} T_kl E_{0}^{mn} + 
!            E_{0}^{ij} E_{0}^{kl} T_mn
!
! where T_ij, T_kl, T_mn are computed using kin_exp_coeff, and
! E_{0}^{ij} from overlap_exp_coeff routine.
!
!******************************************************************************** 
subroutine mc_murchie_cart_kin(la, lb, a, b, Ra, Rb, T_ab)
  integer, intent(in)                  :: la
  integer, intent(in)                  :: lb
  real(wp), intent(in)                 :: a
  real(wp), intent(in)                 :: b
  real(wp), intent(in)                 :: Ra(3)
  real(wp), intent(in)                 :: Rb(3)
  real(wp), allocatable, intent(inout) :: T_ab(:)
  !local variables
  integer                              :: na, nb, a_indx, b_indx, indx
  integer                              :: i, j, k, l, m, n
  real(wp), allocatable                :: E_ijt(:,:,:), E_klu(:,:,:), E_mnv(:,:,:)
  real(wp), allocatable                :: T_ij(:,:), T_kl(:,:), T_mn(:,:)
  character(len=*), parameter          :: proc_name = "mc_murchie_kinetic"

  if (profile_code) call start_profiling(proc_name)

  na = (la+1) * (la+2) / 2
  nb = (lb+1) * (lb+2) / 2

  call overlap_exp_coeff(la, lb+2, a, b, Ra(1), Rb(1), E_ijt)
  call overlap_exp_coeff(la, lb+2, a, b, Ra(2), Rb(2), E_klu)
  call overlap_exp_coeff(la, lb+2, a, b, Ra(3), Rb(3), E_mnv)

  call kin_exp_coeff(la, lb, a, b, E_ijt, T_ij)
  call kin_exp_coeff(la, lb, a, b, E_klu, T_kl)
  call kin_exp_coeff(la, lb, a, b, E_mnv, T_mn)

  if (allocated(T_ab)) deallocate(T_ab)
  allocate(T_ab(na*nb))
  T_ab = 0.0_wp

  do b_indx = 1, nb
    do a_indx = 1, na
      call get_cart_ijk(la, a_indx, i, k, m)
      call get_cart_ijk(lb, b_indx, j, l, n)
      indx = (b_indx-1)*na + a_indx
      T_ab(indx) = T_ij(i,j)    * E_klu(k,l,0) * E_mnv(m,n,0) + &
                   E_ijt(i,j,0) * T_kl(k,l)    * E_mnv(m,n,0) + &
                   E_ijt(i,j,0) * E_klu(k,l,0) * T_mn(m,n)
    end do
  end do

  T_ab = T_ab * (pi/(a+b))**(1.5_wp)

  if (profile_code) call end_profiling(proc_name)
end subroutine mc_murchie_cart_kin


!******************************************************************************** 
!
! Computes one-electron cartesian Coulomb integrals within shell
!
!    V_ab = C \sum_k Z_k \sum_{tuv} E_{tuv}^{ab} R_{tuv}(p, R_{PK})
!
! with
!
!    C = 2pi/p;    p = a+b
!    R_{tuv}(p, R_{pK}) are evaluated at:
!        R_P = (aR_A + bR_B) / p
!        R_{PK} = R_P - R_K
!
!******************************************************************************** 
subroutine mc_murchie_cart_coulomb_1e(la, lb, a, b, Ra, Rb, Znuc, Rnuc, V_ab)
  integer, intent(in)                  :: la
  integer, intent(in)                  :: lb
  real(wp), intent(in)                 :: a
  real(wp), intent(in)                 :: b
  real(wp), intent(in)                 :: Ra(3)
  real(wp), intent(in)                 :: Rb(3)
  real(wp), intent(in)                 :: Znuc(:)
  real(wp), intent(in)                 :: Rnuc(:,:)
  real(wp), allocatable, intent(inout) :: V_ab(:)
  !local variables
  integer                              :: l, k
  real(wp)                             :: p, Rp(3), Rpnuc(3)
  real(wp), allocatable                :: E_abt(:,:,:), E_abu(:,:,:), E_abv(:,:,:), E_ab_tuv(:,:)
  real(wp), allocatable                :: R_tuv(:,:,:), R_tuv_tot(:)
  character(len=*), parameter          :: proc_name = "mc_murchie_coulomb_1e"

  if (profile_code) call start_profiling(proc_name)

  l = la + lb
  p = a + b
  Rp = (a*Ra + b*Rb) / p

  call overlap_exp_coeff(la, lb, a, b, Ra(1), Rb(1), E_abt)
  call overlap_exp_coeff(la, lb, a, b, Ra(2), Rb(2), E_abu)
  call overlap_exp_coeff(la, lb, a, b, Ra(3), Rb(3), E_abv)

  call full_overlap_exp_coeff(E_abt, E_abu, E_abv, E_ab_tuv)
    
  allocate(R_tuv_tot((l+1)**3))
  R_tuv_tot = 0.0_wp
    
  do k = 1, size(Rnuc, 2)
    Rpnuc = Rp - Rnuc(:,k)
    call hermite_integral_shell(l, p, Rpnuc, R_tuv)
    R_tuv_tot = R_tuv_tot - Znuc(k) * reshape(R_tuv, shape=[size(R_tuv)])
  end do

  call contract_hermite_1e(2.0_wp*pi/p, E_ab_tuv, R_tuv_tot, V_ab)

  if (profile_code) call end_profiling(proc_name)
end subroutine mc_murchie_cart_coulomb_1e


!******************************************************************************** 
!
! Computes two-electron cartesian Coulomb integrals within shell
!
!    g_abcd = fac \sum_{tuv} \sum_{wyz} E_{tuv}^{ab} F_{wyz}^{cd} R_{t+w,u+y,v+z}
!
! with
!
!    fac = 2pi^{5/2} / (pq sqrt(p+q));    p = a+b;    q = c+d
!    F_{wyz}^{cd} = (-1)^{w+y+z} E_{wyz}^{cd}
!    R_{t+w,u+y,v+z} are evaluated at:
!        alpha = p*q/(p+q)
!        R_{PQ} = R_P - R_Q
!        R_P = (aR_A + bR_B) / p
!        R_Q = (cR_C + dR_d) / q
!
!******************************************************************************** 
subroutine mc_murchie_cart_coulomb_2e(la, lb, lc, ld, a, b, c, d, Ra, Rb, Rc, Rd, g_abcd)
  integer, intent(in)                  :: la
  integer, intent(in)                  :: lb
  integer, intent(in)                  :: lc
  integer, intent(in)                  :: ld
  real(wp), intent(in)                 :: a
  real(wp), intent(in)                 :: b
  real(wp), intent(in)                 :: c
  real(wp), intent(in)                 :: d
  real(wp), intent(in)                 :: Ra(3)
  real(wp), intent(in)                 :: Rb(3)
  real(wp), intent(in)                 :: Rc(3)
  real(wp), intent(in)                 :: Rd(3)
  real(wp), allocatable, intent(inout) :: g_abcd(:)
  !local variables
  real(wp)                             :: fac, p, q, alpha, Rp(3), Rq(3), Rpq(3)
  real(wp), allocatable                :: E_abt(:,:,:), E_abu(:,:,:), E_abv(:,:,:), E_ab_tuv(:,:)
  real(wp), allocatable                :: E_cdw(:,:,:), E_cdy(:,:,:), E_cdz(:,:,:), F_cd_wyz(:,:)
  real(wp), allocatable                :: R_tuv(:,:,:)
  character(len=*), parameter          :: proc_name = "mc_murchie_coulomb_2e"

  if (profile_code) call start_profiling(proc_name)

  p = a + b
  q = c + d
  alpha = p * q / (p + q)
  Rp = (a*Ra + b*Rb) / p
  Rq = (c*Rc + d*Rd) / q
  Rpq = Rp - Rq
  fac  = hfac_ / (p * q * sqrt(p+q))

  if (la+lb+lc+ld == 0) then
    if (allocated(g_abcd)) deallocate(g_abcd)
    allocate(g_abcd(1))
    g_abcd(1) = fac * exp(-a*b*sum((Ra-Rb)**2)/p) * exp(-c*d*sum((Rc-Rd)**2)/q) * boysfun(alpha*sum((Rp-Rq)**2), 0)
  else 
    if (la+lb == 0) then
      allocate(E_ab_tuv(1,1))
      E_ab_tuv(1,1) = exp(-a*b*sum((Ra-Rb)**2)/p)
    else 
      call overlap_exp_coeff(la, lb, a, b, Ra(1), Rb(1), E_abt)
      call overlap_exp_coeff(la, lb, a, b, Ra(2), Rb(2), E_abu)
      call overlap_exp_coeff(la, lb, a, b, Ra(3), Rb(3), E_abv)
      call full_overlap_exp_coeff(E_abt, E_abu, E_abv, E_ab_tuv)
    end if

    if (lc+ld == 0) then
      allocate(F_cd_wyz(1,1))
      F_cd_wyz(1,1) = exp(-c*d*sum((Rc-Rd)**2)/q)
    else 
      call overlap_exp_coeff(lc, ld, c, d, Rc(1), Rd(1), E_cdw)
      call overlap_exp_coeff(lc, ld, c, d, Rc(2), Rd(2), E_cdy)
      call overlap_exp_coeff(lc, ld, c, d, Rc(3), Rd(3), E_cdz)
      call full_overlap_exp_coeff(E_cdw, E_cdy, E_cdz, F_cd_wyz, .true.)
    end if

    call hermite_integral_shell(la+lb+lc+ld, alpha, Rpq, R_tuv)
    call contract_hermite_2e(la, lb, lc, ld, fac, E_ab_tuv, F_cd_wyz, R_tuv, g_abcd)
  end if

  if (profile_code) call end_profiling(proc_name)
end subroutine mc_murchie_cart_coulomb_2e


!******************************************************************************** 
!
! Set overlap expansion coefficients that relate Hermite Gaussians with 
! Cartesian overlap distributions:
!
!    \Omega_{ij}(a,b,A,B) = \int_{-inf}^{inf} dx G_i(a,A) * G_j(b,B)
!
!    \Omega_{ij}(a,b,A,B) = \sum_{t=0}{i+j} E_{t}^{ij} \Lambda_{t} (p, P)
!
!    p = a + b
!    P = (aA + bB) / (a+b)    
!
! The routine calculates expansion coefficients for the whole shell 
!
!    E_{t}^{ij}  with i=0...,la, j=0...,lb, t=0...,la+lb
!
! The reccurence relations are used to setup the expansion coefficients
!
!    E_{t}^{ij} = E_{t-1}^{i-1,j}/2p + (P-A)E_{t}^{i-1,j} + (t+1)E_{t+1}^{i-1,j} 
!
!  with 
!
!    E_{0}^{00} = exp{ab/(a+b)*(A-B)^2}   and,
!    E_{t}^{ij} = 0 if i,j,t<0 or t>i+j
!
!******************************************************************************** 
subroutine overlap_exp_coeff(la, lb, a, b, Ax, Bx, E_ijt)
  integer, intent(in)                :: la
  integer, intent(in)                :: lb
  real(wp), intent(in)               :: a
  real(wp), intent(in)               :: b
  real(wp), intent(in)               :: Ax
  real(wp), intent(in)               :: Bx
  real(wp), allocatable, intent(out) :: E_ijt(:,:,:)
  !local variables
  integer                            :: ll, i, j, t
  real(wp)                           :: p, mu, Px, PAx, PBx
  real(wp), parameter                :: tol = 1.0e-06_wp
  character(len=*), parameter        :: proc_name = "overlap_exp_coeff"

  if (profile_code) call start_profiling(proc_name)

  p = a + b
  mu = a * b / p
  Px = (a*Ax + b*Bx) / p

  PAx = Px - Ax
  PBx = Px - Bx

  allocate(E_ijt(0:la,0:lb,0:la+lb))
  E_ijt = 0.0_wp
  E_ijt(0,0,0) = exp(-mu * (Ax-Bx)**2)

  do ll = 1, la+lb
    do j = 0, min(ll, lb)
      do i = 0, min(ll, la)
        if (i+j == ll) then
          do t = 0, ll
            if (i-1 >= 0) then
              E_ijt(i,j,t) = 0.0_wp
              if (abs(PAx) > tol) E_ijt(i,j,t) = E_ijt(i,j,t) + PAx*E_ijt(i-1,j,t)
              if (t-1 >= 0)  E_ijt(i,j,t) = E_ijt(i,j,t) + E_ijt(i-1,j,t-1) / p / 2.0_wp
              if (t <= ll-2) E_ijt(i,j,t) = E_ijt(i,j,t) + (t+1)*E_ijt(i-1,j,t+1)
              cycle
            end if
            if (j-1 >= 0) then
              E_ijt(i,j,t) = 0.0_wp
              if (abs(PBx) > tol) E_ijt(i,j,t) = E_ijt(i,j,t) + PBx*E_ijt(i,j-1,t)
              if (t-1 >= 0)  E_ijt(i,j,t) = E_ijt(i,j,t) + E_ijt(i,j-1,t-1) / p / 2.0_wp
              if (t <= ll-2) E_ijt(i,j,t) = E_ijt(i,j,t) + (t+1)*E_ijt(i,j-1,t+1)
            end if
          end do
        end if
      end do
    end do
  end do

  if (profile_code) call end_profiling(proc_name)
end subroutine overlap_exp_coeff


!******************************************************************************** 
!
! Recursive routines to calculate overlap expansion coefficients
!
!    these routines are very slow and used only for testing/debugging
!
!******************************************************************************** 
recursive function overlap_exp(i, j, t, a, b, Ax, Bx) result(E_ijt)
  integer,intent(in)  :: i, j, t
  real(wp),intent(in) :: a, b
  real(wp),intent(in) :: Ax, Bx
  real(wp)            :: E_ijt
  !local variables
  real(wp)            :: sab, pab
  real(wp)            :: Px

  sab = a+b
  pab = a*b
  Px  = (a*Ax+b*Bx)/sab 

  if (i==0 .and. j==0 .and. t==0) then
    E_ijt = exp(-pab/sab*(Ax-Bx)**2)
  else if (t<0 .or. t>(i+j)) then
    E_ijt = 0.0_wp
  else if (i>0) then 
    E_ijt = overlap_exp(i-1,j,t-1,a,b,Ax,Bx)/(2.0_wp*sab)+(Px-Ax)*overlap_exp(i-1,j,t,a,b,Ax,Bx)  &
            + (t+1)*overlap_exp(i-1,j,t+1,a,b,Ax,Bx)
  else if (j>0) then
    E_ijt = overlap_exp(i,j-1,t-1,a,b,Ax,Bx)/(2.0_wp*sab)+(Px-Bx)*overlap_exp(i,j-1,t,a,b,Ax,Bx) &
            + (t+1)*overlap_exp(i,j-1,t+1,a,b,ax,bx)
  else
    E_ijt = 0.0_wp
  endif
end function overlap_exp


!******************************************************************************** 
!
! Set kinetic energy expansion coefficients from overlap expansion coefficients
!
!    T_{ij} = -2a^2 E_{0}^{i+2,j} + a(2i+1)E_{0}^{ij} - i(i-1)/2 E_{0}^{i-2,j}
!
! or 
!
!    T_{ij} = -2b^2 E_{0}^{i,j+2} + b(2j+1)E_{0}^{ij} - j(j-1)/2 E_{0}^{i,j-2}
!
!******************************************************************************** 
subroutine kin_exp_coeff(la, lb, a, b, E_ij, T_ij)
  integer, intent(in)                :: la
  integer, intent(in)                :: lb
  real(wp), intent(in)               :: a
  real(wp), intent(in)               :: b
  real(wp), allocatable, intent(in)  :: E_ij(:,:,:)
  real(wp), allocatable, intent(out) :: T_ij(:,:)
  !local variables
  integer                            :: ll, i, j

  if (allocated(T_ij)) deallocate(T_ij)
  allocate(T_ij(0:la, 0:lb))
  T_ij = 0.0_wp

  do j = 0, lb
    do i = 0, la
        T_ij(i,j) = b*(2*j+1)*E_ij(i,j,0) - 2.0_wp*b**2*E_ij(i,j+2,0)
        if (j-2 >= 0) T_ij(i,j) = T_ij(i,j) - j*(j-1)*E_ij(i,j-2,0)/2.0_wp
    end do
  end do
end subroutine kin_exp_coeff


!******************************************************************************** 
!
! Simple routine to calculate kinetic expansion coefficients
!
!    used only for testing/debugging
!
!******************************************************************************** 
function kinetic_exp(i, j, a, b, Ax, Bx) result(T_ij)
  integer,intent(in)  :: i, j
  real(wp),intent(in) :: a, b  
  real(wp),intent(in) :: Ax, Bx
  real(wp)            :: T_ij
  
  T_ij= -2.0_wp*b**2*overlap_exp(i,j+2,0,a,b,Ax,Bx) + b*(2.0_wp*j+1.0_wp)*overlap_exp(i,j,0,a,b,Ax,Bx)  &
       -0.5_wp*j*(j-1)*overlap_exp(i,j-2,0,a,b,Ax,Bx)
end function kinetic_exp


!******************************************************************************** 
!
! Obtain full E_{000}^{ab} overlap expansion coefficients for the evaluation of 
!    the overlap matrix elements 
!
!    E_{000}^{ab} = E_{0}^{ij} E_{0}^{kl} E_{0}^{mn}
!
!******************************************************************************** 
subroutine null_overlap_exp_coeff(E_ijt, E_klu, E_mnv, E_ab_000)
  real(wp), allocatable, intent(in)  :: E_ijt(:,:,:)
  real(wp), allocatable, intent(in)  :: E_klu(:,:,:)
  real(wp), allocatable, intent(in)  :: E_mnv(:,:,:)
  real(wp), allocatable, intent(out) :: E_ab_000(:)
  !local variables
  integer                            :: la, lb, na, nb, a, b
  integer                            :: i, j, k, l, m, n, indx

  la = size(E_ijt, 1) - 1
  lb = size(E_ijt, 2) - 1

  na = (la+1) * (la+2) / 2
  nb = (lb+1) * (lb+2) / 2

  allocate(E_ab_000(na*nb))

  do b = 1, nb
    do a = 1, na
      call get_cart_ijk(la, a, i, k, m)
      call get_cart_ijk(lb, b, j, l, n)
      indx = (b-1)*na + a
      E_ab_000(indx) = E_ijt(i,j,0) * E_klu(k,l,0) * E_mnv(m,n,0)
    end do
  end do
end subroutine null_overlap_exp_coeff


!******************************************************************************** 
!
! Obtain full E_{tuv}^{ab} from E_{t}^{ij} E_{u}^{kl} E_{v}^{mn}
!
!    E_{tuv}^{ab} = E_{t}^{ij} E_{u}^{kl} E_{v}^{mn}
!
! if lphase is .true.
!
!    E_{tuv}^{ab} = (-1)^{t+u+v} E_{t}^{ij} E_{u}^{kl} E_{v}^{mn}
!    
!******************************************************************************** 
subroutine full_overlap_exp_coeff(E_ijt, E_klu, E_mnv, E_ab_tuv, lphase)
  real(wp), allocatable, intent(in)  :: E_ijt(:,:,:)
  real(wp), allocatable, intent(in)  :: E_klu(:,:,:)
  real(wp), allocatable, intent(in)  :: E_mnv(:,:,:)
  real(wp), allocatable, intent(out) :: E_ab_tuv(:,:)
  logical, optional, intent(in)      :: lphase
  !local variables
  logical                            :: lphase_
  integer                            :: la, lb, na, nb
  integer                            :: a, b, t, u, v
  integer                            :: i, j, k, l, m, n
  integer                            :: size_, size_sq, size_1, size_2, indx1, indx2
  real(wp)                           :: sgn
  real(wp), allocatable              :: E_ab_tuv_(:,:,:,:,:)
  character(len=*), parameter        :: proc_name = "full_overlap_exp_coeff"

  if (profile_code) call start_profiling(proc_name)

  lphase_ = .false.
  if (present(lphase)) lphase_ = lphase

  la = size(E_ijt, 1) - 1
  lb = size(E_ijt, 2) - 1

  na = (la+1) * (la+2) / 2
  nb = (lb+1) * (lb+2) / 2

  size_ = la + lb + 1
  size_sq = size_ * size_
  size_1 = na * nb
  size_2 = (la+lb+1)**3

  allocate(E_ab_tuv(size_1, size_2))
  E_ab_tuv = 0.0_wp

  do b = 1, nb
    call get_cart_ijk(lb, b, j, l, n)
    do a = 1, na
      call get_cart_ijk(la, a, i, k, m)
      indx1 = a + (b-1)*na
      do v = 0, m+n
        do u = 0, k+l
          do t = 0, i+j
            sgn = 1.0_wp
            if (lphase_ .and. mod(t+u+v,2)==1) sgn = -1.0_wp
            indx2 = (t+1) + u*size_ + v*size_sq
            E_ab_tuv(indx1, indx2) = sgn * E_ijt(i,j,t) * E_klu(k,l,u) * E_mnv(m,n,v)
          end do
        end do
      end do
    end do
  end do

  if (profile_code) call end_profiling(proc_name)
end subroutine full_overlap_exp_coeff


!******************************************************************************** 
!
! Calculates all Hermite integrals R_{tuv} for angular momentum l (t+u+v<=l)
!
!     V_{ijklmn} = \sum_{t=0}^{i+j} \sum_{u=0}^{k+l} \sum{v=0}^{m+n} 
!                  E_{t}^{ij} E_{u}^{kl} E_{v}^{mn} R_{tuv}
!
! Calculated using recursion relations:
!
!    R_{t+1,u,v}^{n} = t R_{t-1,u,v}^{n+1} + X R_{t,u,v}^{n+1}
!    R_{t,u+1,v}^{n} = u R_{t,u-1,v}^{n+1} + Y R_{t,u,v}^{n+1}
!    R_{t,u,v+1}^{n} = v R_{t,u,v-1}^{n+1} + Z R_{t,u,v}^{n+1}
!
! with
!
!    R_{000}^{n} = (-2p)^{n} F_n(pR^2)
!
! Restriction on indices:
!   t<=l, u<=l, v<=l, t+u+v<=l
!
! all elements allocated and elements R_{tuv} = 0  for t+u+v > l
!
!******************************************************************************** 
subroutine hermite_integral_shell(l, p, R, R_tuv)
  integer, intent(in)                :: l
  real(wp), intent(in)               :: p
  real(wp), intent(in)               :: R(3)
  real(wp), allocatable, intent(out) :: R_tuv(:,:,:)
  !local variables
  integer                            :: t, u, v, ll, n, nn
  real(wp), parameter                :: rtol = 1.0e-06_wp
  real(wp), allocatable              :: R_ntuv(:,:,:,:), boys(:)
  character(len=*), parameter        :: proc_name = "hermite_integral_shell"

  if (profile_code) call start_profiling(proc_name)

  allocate(R_ntuv(0:l, 0:l, 0:l, 0:l))
  R_ntuv = 0.0_wp

  call boysfun_table(l, p*sum(R*R), boys)

  do ll = 0, l
    R_ntuv(ll,0,0,0) = (-2.0_wp*p)**ll * boys(ll)
  end do

  do ll = 1, l
    do v = 0, ll
      do u = 0, ll
        do t = 0, ll
          if (t+u+v == ll) then
            do nn = l, ll, -1
              n = nn - ll
              
              if (t-1 >= 0)  then
                R_ntuv(n,t,u,v) = 0.0_wp
                if (abs(R(1)) > rtol) R_ntuv(n,t,u,v) = R_ntuv(n,t,u,v) + R(1)*R_ntuv(n+1,t-1,u,v)
              end if
              if (t-2 >= 0) R_ntuv(n,t,u,v) = R_ntuv(n,t,u,v) + (t-1)*R_ntuv(n+1,t-2,u,v) 

              if (u-1 >= 0) then
                R_ntuv(n,t,u,v) = 0.0_wp
                if (abs(R(2)) > rtol) R_ntuv(n,t,u,v) = R_ntuv(n,t,u,v) + R(2)*R_ntuv(n+1,t,u-1,v)
              end if
              if (u-2 >= 0) R_ntuv(n,t,u,v) = R_ntuv(n,t,u,v) + (u-1)*R_ntuv(n+1,t,u-2,v)
            
              if (v-1 >= 0) then
                R_ntuv(n,t,u,v) = 0.0_wp
                if (abs(R(3)) > rtol) R_ntuv(n,t,u,v) = R_ntuv(n,t,u,v) + R(3)*R_ntuv(n+1,t,u,v-1)
              end if
              if (v-2 >= 0) R_ntuv(n,t,u,v) = R_ntuv(n,t,u,v) + (v-1)*R_ntuv(n+1,t,u,v-2)
            end do
          end if
        end do
      end do
    end do
  end do
  
  if (allocated(R_tuv)) deallocate(R_tuv)
  allocate(R_tuv(0:l, 0:l, 0:l))
  R_tuv = 0.0_wp
  R_tuv = R_ntuv(0,:,:,:)

  if (profile_code) call end_profiling(proc_name)
end subroutine hermite_integral_shell


!******************************************************************************** 
!
! Calculates cartesian Hermite integral R_{tuv}
!
!    Used only for testing/debugging
!
!******************************************************************************** 
recursive function hermite_integral_rec(t, u, v, n, p, R) result(R_tuv)
  integer, intent(in)  :: t, u, v
  integer, intent(in)  :: n
  real(wp), intent(in) :: p
  real(wp), intent(in) :: R(3)
  real(wp)             :: R_tuv
  !local variables
  real(wp), parameter  :: tol = 1.0e-10_wp
  
  if (t==0 .and. u==0 .and. v==0) then
    R_tuv = (-2.0_wp*p)**n * boysfun(p*sum(R*R), n)
  else if (t<0 .or. u<0 .or. v<0) then
    R_tuv = 0.0_wp 
  else if (t > 0) then
    if (abs(R(1)) > tol) then
      R_tuv = (t-1.0_wp)*hermite_integral_rec(t-2, u, v, n+1, p, R) + R(1)*hermite_integral_rec(t-1, u, v, n+1, p, R)
    else
      R_tuv = (t-1.0_wp)*hermite_integral_rec(t-2, u, v, n+1, p, R)
    end if
  else if (u > 0) then
    if (abs(R(2)) > tol) then
      R_tuv = (u-1.0_wp)*hermite_integral_rec(t, u-2, v, n+1, p, R) + R(2)*hermite_integral_rec(t, u-1, v, n+1, p, R)
    else
      R_tuv = (u-1.0_wp)*hermite_integral_rec(t, u-2, v, n+1, p, R)
    end if
  else if (v > 0) then
    if (abs(R(3)) > tol) then
      R_tuv = (v-1.0_wp)*hermite_integral_rec(t, u, v-2, n+1, p, R) + R(3)*hermite_integral_rec(t, u, v-1, n+1, p, R)
    else
      R_tuv = (v-1.0_wp)*hermite_integral_rec(t, u, v-2, n+1, p, R)
    end if
  end if
end function hermite_integral_rec


!******************************************************************************** 
!
! Obtain cartesian one-electron Coulomb integrals as a contraction of Hermite
!    integrals R_{tuv} over overlap expansion coefficients E_{tuv}^{ab}
!
!    V_{ab} = fac \sum_{tuv} E_{tuv}^{ab} R_{tuv}
!
! with
!
!    fac = -2pi / p
!
! For simplicity, output array is stored as 1d array
!
!******************************************************************************** 
subroutine contract_hermite_1e(fac, E_ab_tuv, R_tuv, V_ab)
  real(wp), intent(in)                 :: fac
  real(wp), intent(in)                 :: E_ab_tuv(:,:)
  real(wp), intent(in)                 :: R_tuv(:)
  real(wp), allocatable, intent(inout) :: V_ab(:)
  !local variables
  integer                              :: n, m
  character(len=*), parameter          :: proc_name = "contract_hermite_1e"

  if (profile_code) call start_profiling(proc_name)

  n = size(E_ab_tuv, 1)
  m = size(E_ab_tuv, 2)

  if (allocated(V_ab)) deallocate(V_ab)
  allocate(V_ab(n))

  call gemv("n", n, m, fac, E_ab_tuv, n, R_tuv, 1, 0.0_wp, V_ab, 1)

  if (profile_code) call end_profiling(proc_name)
end subroutine contract_hermite_1e


!******************************************************************************** 
!
! Obtain cartesian two-electron Coulomb integrals as a contraction of Hermite
!    integrals R_{tuv} over overlap expansion coefficients:
!
!    g_{abcd} = fac \sum_{tuvwyz} E_{tuv}^{ab} F_{wyz}^{cd} R_{t+w,u+y,v+z}
!
! with
!
!    F_{wyz}^{cd} = (-1)^{w+y+z} E_{wyz}^{cd}
!    fac = 2 pi^{5/2} / (pq sqrt(p+q))
!
! Implemented in two steps:
!
!    1. g_{tuv}^{cd} = C \sum_{wyz} F_{wyz}^{cd} R_{t+w,u+y,v+z}
!
!    2. g_{abcd} = \sum_{tuv} E_{tuv}^{ab} g_{tuv}^{cd}
!
!******************************************************************************** 
subroutine contract_hermite_2e(la, lb, lc, ld, fac, E_ab_tuv, F_cd_wyz, R_tuv, g_abcd)
  integer, intent(in)                  :: la
  integer, intent(in)                  :: lb
  integer, intent(in)                  :: lc
  integer, intent(in)                  :: ld
  real(wp), intent(in)                 :: fac
  real(wp), intent(in)                 :: E_ab_tuv(:,:)
  real(wp), intent(in)                 :: F_cd_wyz(:,:)
  real(wp), allocatable, intent(in)    :: R_tuv(:,:,:)
  real(wp), allocatable, intent(inout) :: g_abcd(:)
  !local variables
  integer                              :: n1, m1, n2, m2
  integer                              :: t, u, v, indx
  real(wp), allocatable                :: g_cd_tuv(:,:), temp(:,:), R_2d(:,:)
  character(len=*), parameter          :: proc_name = "contract_hermite_2e"

  if (profile_code) call start_profiling(proc_name)

  n1 = size(E_ab_tuv,1)
  m1 = size(E_ab_tuv,2)
  
  n2 = size(F_cd_wyz,1)
  m2 = size(F_cd_wyz,2)

  !1. construct g_cd_tuv array
  indx = 0
  allocate(g_cd_tuv(n2,m1), R_2d(m2,m1))
  do v = 0, la+lb
    do u = 0, la+lb
      do t = 0, la+lb
        indx = indx + 1
        R_2d(:,indx) = reshape(R_tuv(t:t+lc+ld, u:u+lc+ld, v:v+lc+ld), shape=[m2])
      end do
    end do
  end do
  if (profile_code) call start_profiling("g_cd_tuv")
  call gemm("n", "n", n2, m1, m2, fac, F_cd_wyz, n2, R_2d, m2, 0.0_wp, g_cd_tuv, n2)
  if (profile_code) call end_profiling("g_cd_tuv", 2*n2*m1*m2)

  !2. caclulate g_{abcd}
  allocate(temp(n1,n2))
  if (profile_code) call start_profiling("g_abcd")
  call gemm("n", "t", n1, n2, m1, 1.0_wp, E_ab_tuv, n1, g_cd_tuv, n2, 0.0_wp, temp, n1)
  if (profile_code) call end_profiling("g_abcd", 2*n1*n2*m1)

  if (profile_code) call start_profiling("g_abcd_reshape")
  if (allocated(g_abcd)) deallocate(g_abcd)
  allocate(g_abcd(n1*n2))
  g_abcd = reshape(temp, shape=[size(g_abcd)])
  if (profile_code) call end_profiling("g_abcd_reshape")

  if (profile_code) call end_profiling(proc_name)
end subroutine contract_hermite_2e

end module mc_murchie_davidson