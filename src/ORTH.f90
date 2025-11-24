! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module orth

#include "preproc.inc"
use constants
use mpi, only: comm_world
use profiling
use qmcfort_io
use lapack
use standalone
use hamilton_vars

implicit none

interface get_integral
  module procedure get_integral_trafo, get_integral_dec, get_integral_dec_trafo
end interface

interface transform_1e_ham
  module procedure transform_1e_ham_1, transform_1e_ham_2
end interface transform_1e_ham

interface transform_2e_ham
  module procedure transform_2e_ham_1, transform_2e_ham_2, transform_2e_ham_4
end interface transform_2e_ham

type basis_trafo_descriptor
  integer                         :: orth_mode = 3
  character(len=charlen)          :: basis_orig
  character(len=charlen)          :: basis
  logical                         :: basis_changed
contains
  procedure                       :: init            => init_trafo_des
  procedure                       :: reader          => trafo_des_reader
  procedure                       :: print_header    => print_basis_trafo_header
  procedure                       :: print_footer    => print_basis_trafo_footer
end type basis_trafo_descriptor

contains

!********************************************************************************
!
!       Copy information from ham_des to trafo_des and read additional info
!
!******************************************************************************** 
subroutine init_trafo_des(trafo_des, ham_des, basis)
  class(basis_trafo_descriptor), intent(out) :: trafo_des
  type(hamil_descriptor), intent(inout)      :: ham_des
  character(len=*), optional, intent(in)     :: basis

  trafo_des%basis_orig = ham_des%basis
  trafo_des%basis = ham_des%basis
  call trafo_des%reader()

  if (present(basis)) trafo_des%basis = basis
  trafo_des%basis_changed = trim(ham_des%basis) /= trim(trafo_des%basis)
  ham_des%basis = trafo_des%basis
end subroutine init_trafo_des


!******************************************************************************** 
!
!       Small internal reader for basis transformation
!
!******************************************************************************** 
subroutine trafo_des_reader(trafo_des)
  use qmcfort_in, only: add_input
  class(basis_trafo_descriptor), intent(inout) :: trafo_des
  
  call add_input("basis"       , trafo_des%basis)
  call add_input("orth"        , trafo_des%orth_mode)
  call add_input("orth_mode"   , trafo_des%orth_mode)
end subroutine trafo_des_reader


!******************************************************************************** 
!
!       Print header for the basis transformation
!
!******************************************************************************** 
subroutine print_basis_trafo_header(trafo_des)
  class(basis_trafo_descriptor), intent(in) :: trafo_des
  !local variables
  character(len=charlen)                    :: intfor, logfor, realfor, charfor
  !
  intfor  = '(1X,T5,A,T50,"= ",10I8)'
  logfor  = '(1X,T5,A,T50,"= ",10L)'
  realfor = '(1X,T5,A,T50,"= ",10F12.6)'
  charfor = '(1X,T5,A,T50,"= ",A)'
  !
  if (comm_world%mpirank == 0) then
    write(io%screen%funit,*)       "BASIS TRANSFORMATION"
    write(io%screen%funit,*)       starshort
    write(io%screen%funit,*)       "basis trafo descriptor:"
    write(io%screen%funit,charfor) "from basis"          , trim(trafo_des%basis_orig)
    write(io%screen%funit,charfor) "to basis"            , trim(trafo_des%basis)
    write(io%screen%funit,intfor)  "ortgonalization mode", trafo_des%orth_mode
    write(io%qmcfort_out%funit,*)       "BASIS TRANSFORMATION"
    write(io%qmcfort_out%funit,*)       starshort
    write(io%qmcfort_out%funit,*)       "basis trafo descriptor:"
    write(io%qmcfort_out%funit,charfor) "from basis"          , trim(trafo_des%basis_orig)
    write(io%qmcfort_out%funit,charfor) "to basis"            , trim(trafo_des%basis)
    write(io%qmcfort_out%funit,intfor)  "ortgonalization mode", trafo_des%orth_mode
  end if
end subroutine print_basis_trafo_header  


!******************************************************************************** 
!
!       Print information after the basis transformation
!
!******************************************************************************** 
subroutine print_basis_trafo_footer(trafo_des)
  class(basis_trafo_descriptor), intent(in) :: trafo_des
  
  if (comm_world%mpirank == 0) then
    write(io%screen%funit,*) starlong
    write(io%qmcfort_out%funit,*) starlong
  end if
end subroutine print_basis_trafo_footer


!******************************************************************************** 
!
!       Calculate transformation matrix trafo
!
!******************************************************************************** 
subroutine calculate_trafo_matrix(trafo_des, overlap, phi, trafo)
  type(basis_trafo_descriptor), intent(in)  :: trafo_des
  real(wp), intent(in)                      :: overlap(:,:,:)
  real(wp), intent(in)                      :: phi(:,:,:)
  real(wp), allocatable, intent(out)        :: trafo(:,:,:,:)
  !
  if (trim(trafo_des%basis) == "mo") then
    call mean_field_orthogonalization(phi, trafo)
  else if (trim(trafo_des%basis) == "ao_orth") then
    if (trafo_des%orth_mode == 1) then
      call symmetric_orthogonalization(overlap, trafo)
    else if (trafo_des%orth_mode == 2) then
      call canonical_orthogonalization(overlap, trafo)
    else if (trafo_des%orth_mode == 3) then
      call cholesky_orthogonalization(overlap, trafo)
    else 
      if (comm_world%mpirank == 0) write(*,*) "Orthogonalization mode not recognized - Program stops now"
      stop 
    end if
  end if 
end subroutine calculate_trafo_matrix


!******************************************************************************** 
!
!       Mean Field orthogonalization
!
!         Forward transformation X      = X(:,:,1,spin) = phi(:,:,spin)
!         Back transformation    X^{-1} = X(:,:,2,spin) = phi(:,:,spin)^{-1}
!
!******************************************************************************** 
subroutine mean_field_orthogonalization(phi, trafo)
  real(wp), intent(in)               :: phi(:,:,:)
  real(wp), allocatable, intent(out) :: trafo(:,:,:,:)
  !local variables
  integer                            :: n, spin, ispin
  !
  n = size(phi, 1)
  ispin = size(phi, 3)
  if (.not. allocated(trafo)) allocate(trafo(n,n,2,ispin))
  !
  do spin = 1, ispin
    trafo(:,:,1,spin) = phi(:,:,spin)
    trafo(:,:,2,spin) = phi(:,:,spin)
    call inverse(trafo(:,:,2,spin))
  end do
end subroutine mean_field_orthogonalization


!******************************************************************************** 
!
!       Symmetric orthogonalization scheme:
!
!         Forward transformation X      = X(:,:,spin,1) = S^{-1/2}(:,:,spin)
!         Back transformation    X^{-1} = X(:,:,spin,2) = S^{1/2}(:,:,spin)
!
!******************************************************************************** 
subroutine symmetric_orthogonalization(overlap, trafo)
  real(wp), intent(in)               :: overlap(:,:,:)
  real(wp), allocatable, intent(out) :: trafo(:,:,:,:)
  !local variables
  integer                            :: i, n, spin, ispin
  real(wp), allocatable              :: s(:), U(:,:), help(:,:)
  !
  n = size(overlap, 1)
  ispin = size(overlap, 3)
  if (.not. allocated(trafo)) allocate(trafo(n,n,2,ispin))
  allocate(s(n), U(n,n), help(n,n))
  !
  do spin = 1, ispin
    U = overlap(:,:,spin)
    call syev(U, s)             !diagonalize overlap matrix
    !
    do i = 1, n
      s(i) = sqrt(s(i))
    end do
    !
    do i = 1, n                 !Calculate S^{1/2}
      help(:,i) = U(:,i) * s(i)
    end do
    !
    call gemm("n", "t", n, n, n, 1.0_wp, help, n, U, n, 0.0_wp, trafo(:,:,2,spin), n)
    !
    do i = 1, n                 !Calculate S^{-1/2} 
      help(:,i) = U(:,i) / s(i)
    end do
    !
    call gemm("n", "t", n, n, n, 1.0_wp, help, n, U, n, 0.0_wp, trafo(:,:,1,spin), n)
  end do
end subroutine symmetric_orthogonalization


!******************************************************************************** 
!
!       Canonical orthogonalization scheme:
!
!         Forward transformation  X      = X(:,:,1) = Us^{-1/2}
!         Back transforamtion     X^{-1} = X(:,:,2) = U^{+}s^{1/2}
!
!******************************************************************************** 
subroutine canonical_orthogonalization(overlap, trafo)
  real(wp), intent(in)               :: overlap(:,:,:)
  real(wp), allocatable, intent(out) :: trafo(:,:,:,:)
  !local variables
  integer                            :: i, n, spin, ispin
  real(wp), allocatable              :: s(:), U(:,:), help(:,:)
  !
  n = size(overlap, 1)
  ispin = size(overlap, 3)
  if (.not. allocated(trafo)) allocate(trafo(n,n,2,ispin))
  allocate(s(n), U(n,n), help(n,n))
  !
  do spin = 1, ispin
    U = overlap(:,:,spin)
    call syev(U, s)             !diagonalize overlap matrix
    !
    do i = 1, n
      s(i) = sqrt(s(i))
    end do
    !
    do i = 1, n                 !Calculate X^{-1} = s^{1/2}U^T
      trafo(i,:,2,spin) = s(i) * U(:,i)
    end do
    !
    do i = 1, n                 !Calculate X = Us^{-1/2}
      trafo(:,i,1,spin) = U(:,i) / s(i)
    end do
  end do
end subroutine canonical_orthogonalization


!******************************************************************************** 
!
!       Cholesky orthogonalization scheme:
!
!         Forward transformation  X      = X(:,:,1) = L^{-1}
!         Back transforamtion     X^{-1} = X(:,:,2) = L
!
!******************************************************************************** 
subroutine cholesky_orthogonalization(overlap, trafo)
  real(wp), intent(in)               :: overlap(:,:,:)
  real(wp), allocatable, intent(out) :: trafo(:,:,:,:)
  !local variables
  integer                            :: i, j, k, n, spin, ispin
  real(wp), allocatable              :: s(:), U(:,:), help(:,:)
  !
  n = size(overlap, 1)
  ispin = size(overlap, 3)
  if (.not. allocated(trafo)) allocate(trafo(n,n,2,ispin))
  allocate(s(n), U(n,n), help(n,n))
  !
  do spin = 1, ispin
    U = overlap(:,:,spin)
    call potrf(U, "u")
    trafo(:,:,2,spin) = U
    !
    call inverse_tri(U, "u", "n")
    trafo(:,:,1,spin) = U
    !Fill lower triangular part with zeros, because it could be unreferenced
    do k = 1, 2
      do i = 1, n
        do j = 1, i-1
          trafo(i,j,k,spin) = 0.0_wp
        end do
      end do
    end do
  end do
end subroutine cholesky_orthogonalization


!******************************************************************************** 
!       Transforms one-electron Hamiltonian - applies for ham1 and overlap
!               h --> X^{+} h X
!******************************************************************************** 
subroutine transform_ham1(trafo, ham1)
  real(wp), intent(in)        :: trafo(:,:,:,:)
  real(wp), intent(inout)     :: ham1(:,:,:)
  !local variables
  integer                     :: spin
  character(len=*), parameter :: proc_name = "transform_ham1"
  
  if (use_profiling) call start_profiling(proc_name)

  do spin = 1, size(ham1, 3)
    call transform_1e_ham(ham1(:,:,spin), trafo(:,:,1,spin))
  end do
  
  if (comm_world%mpirank == 0) write(io%screen%funit,*) "ham1 transformed"

  if (use_profiling) call end_profiling(proc_name)
end subroutine transform_ham1


!******************************************************************************** 
!
!       Transform two-electron Hamiltonian
!
!******************************************************************************** 
subroutine transform_ham2(trafo, ham2, ham2_dec)
  real(wp), intent(in)                 :: trafo(:,:,:,:)
  real(wp), allocatable, intent(inout) :: ham2(:,:,:,:,:)
  real(wp), allocatable, intent(inout) :: ham2_dec(:,:,:)
  !local variables
  character(len=*), parameter          :: proc_name = "transform_ham2"
  
  if (use_profiling) call start_profiling(proc_name)

  if (allocated(ham2)) call transform_ham2_full(trafo, ham2)
  if (allocated(ham2_dec)) call transform_ham2_dec(trafo, ham2_dec)

  if (use_profiling) call end_profiling(proc_name)
end subroutine transform_ham2


!******************************************************************************** 
!
!       Transform full two-body Hamiltonian
!
!******************************************************************************** 
subroutine transform_ham2_full(trafo, ham2)
  real(wp), intent(in)                     :: trafo(:,:,:,:)
  real(wp), allocatable, intent(inout)     :: ham2(:,:,:,:,:)
  !local variables
  integer                                  :: spin
  
  do spin = 1, size(trafo, 4)
    call transform_2e_ham(ham2(:,:,:,:,spin), trafo(:,:,1,spin))
  end do
  
  if (size(ham2, 5) == 3) call transform_2e_ham(ham2(:,:,:,:,3), trafo(:,:,1,1), trafo(:,:,1,2))
  
  if (comm_world%mpirank == 0) write(io%screen%funit,*) "ham2 transformed"
end subroutine transform_ham2_full


!******************************************************************************** 
!       
!       Transform decomposed two-body Hamiltonian
!
!******************************************************************************** 
subroutine transform_ham2_dec(trafo, ham2_dec)
  real(wp), intent(in)                     :: trafo(:,:,:,:)
  real(wp), allocatable, intent(inout)     :: ham2_dec(:,:,:)
  !local variables
  integer                                  :: spin, n, size_, ng, g
  real(wp), allocatable                    :: temp(:,:)
  
  n = size(trafo, 1)
  size_ = size(ham2_dec, 1)
  ng = size(ham2_dec, 2)

  allocate(temp(n,n))

  do spin = 1, size(ham2_dec, 3)
    do g = 1, ng
      temp = reshape(ham2_dec(:,g,spin), shape=[n,n])
      call transform_1e_ham(temp, trafo(:,:,1,spin))
      ham2_dec(:,g,spin) = reshape(temp, shape=[size_])
    end do
  end do
  
  if (comm_world%mpirank == 0) write(io%screen%funit,*) "ham2_dec transformed"
end subroutine transform_ham2_dec


!******************************************************************************** 
!
!       Transforms orbitals according to x^{-1}
!
!******************************************************************************** 
subroutine transform_orbitals(trafo, phi)
  real(wp), intent(in)        :: trafo(:,:,:,:)
  real(wp), intent(inout)     :: phi(:,:,:)
  !local variables
  integer                     :: spin, sp, n
  real(wp)                    :: phi_temp(size(phi,1), size(phi,2))
  character(len=*), parameter :: proc_name = "transform_orbitals"
  
  if (use_profiling) call start_profiling(proc_name)

  n = size(phi, 1)
  do spin = 1, size(phi, 3)
    sp = spin
    if (size(trafo,4) == 1) sp = 1
    phi_temp = phi(:,:,spin)
    call gemm("n", "n", n, n, n, 1.0_wp, trafo(:,:,2,sp), n, phi_temp, n, 0.0_wp, phi(:,:,spin), n)
  end do
  
  if (comm_world%mpirank == 0) write(io%screen%funit,*) "orbitals transformed"

  if (use_profiling) call end_profiling(proc_name)
end subroutine transform_orbitals


!******************************************************************************** 
!
! Basis transformation of the one-body integrals
!
!    \bar H = T_1^{\dagger} H T_2
!
! transform_1e_ham_1(ham1, trafo1) - square matrices
! transform_1e_ham_2(ham1, trafo1, trafo2)
!
!******************************************************************************** 
subroutine transform_1e_ham_1(ham1, trafo1)
  real(wp), intent(inout) :: ham1(:,:)
  real(wp), intent(in)    :: trafo1(:,:)
  !local variables
  real(wp), allocatable   :: help(:,:), temp(:,:)
  integer                 :: n
  
  n = size(trafo1, 1)

  allocate(help, mold=ham1)
  
  call gemm("n", "n", n, n, n, 1.0_wp, ham1, n, trafo1, n, 0.0_wp, help, n)
  call gemm("t", "n", n, n, n, 1.0_wp, trafo1, n, help, n, 0.0_wp, ham1, n)
end subroutine transform_1e_ham_1

subroutine transform_1e_ham_2(ham1, trafo1, trafo2)
  real(wp), allocatable, intent(inout) :: ham1(:,:)
  real(wp), intent(in)                 :: trafo1(:,:)
  real(wp), intent(in)                 :: trafo2(:,:)
  !local variables
  real(wp), allocatable                :: help(:,:), temp(:,:)
  integer                              :: n1, n2, m1, m2
  
  n1 = size(trafo1, 1)
  m1 = size(trafo1, 2)
  n2 = size(trafo2, 1)
  m2 = size(trafo2, 2)

  allocate(help(n1,m2), temp(m1,m2))
  
  call gemm("n", "n", n1, m2, n2, 1.0_wp, ham1, n1, trafo2, n2, 0.0_wp, help, n1)
  call gemm("t", "n", m1, m2, n1, 1.0_wp, trafo1, n1, help, n1, 0.0_wp, temp, m1)

  call move_alloc(temp, ham1)
end subroutine transform_1e_ham_2


!********************************************************************************
!
! Basis transformation of the two-body Coulomb integrals
!
! transform_2e_ham_1(ham2, trafo1) - assumes square matrices
! transform_2e_ham_2(ham2, trafo1, trafo2) - assumes square matrices &
!    trafo1 transform first two indices and trafo2 last two indices
! transform_2e_ham_2(ham2, trafo1, trafo2, trafo3, trafo4) - general trafo
!
!********************************************************************************
subroutine transform_2e_ham_1(ham2, trafo1)
  real(wp), intent(inout) :: ham2(:,:,:,:)
  real(wp), intent(in)    :: trafo1(:,:)
  !local variables
  integer                 :: i, j, k, n
  real(wp), allocatable   :: tempvec(:)

  n = size(trafo1, 1)

  allocate(tempvec(n))
  
  do k = 1, n
    do j = 1, n
      do i = 1, n
         call mv_prod(trafo1, ham2(:,i,j,k), tempvec)
         ham2(:,i,j,k) = tempvec
       end do
     end do
   end do
   
   do k = 1, n
     do j = 1, n
       do i = 1, n
         call mv_prod(trafo1, ham2(i,:,j,k), tempvec)
         ham2(i,:,j,k) = tempvec
       end do
     end do
   end do
   
   do k = 1, n
     do j = 1, n
       do i = 1, n
         call mv_prod(trafo1, ham2(i,j,:,k), tempvec)
         ham2(i,j,:,k) = tempvec
       end do
     end do
   end do
   
   do k = 1, n
     do j = 1, n
       do i = 1, n
         call mv_prod(trafo1, ham2(i,j,k,:), tempvec)
         ham2(i,j,k,:) = tempvec
       end do
     end do
   end do
end subroutine transform_2e_ham_1

subroutine transform_2e_ham_2(ham2, trafo1, trafo2)
  real(wp), intent(inout) :: ham2(:,:,:,:)
  real(wp), intent(in)    :: trafo1(:,:)
  real(wp), intent(in)    :: trafo2(:,:)
  !local variables
  integer                 :: i, j, k, n
  real(wp), allocatable   :: tempvec(:)

  n = size(trafo1, 1)

  allocate(tempvec(n))
  
  do k = 1, n
    do j = 1, n
      do i = 1, n
         call mv_prod(trafo1, ham2(:,i,j,k), tempvec)
         ham2(:,i,j,k) = tempvec
       end do
     end do
   end do
   
   do k = 1, n
     do j = 1, n
       do i = 1, n
         call mv_prod(trafo1, ham2(i,:,j,k), tempvec)
         ham2(i,:,j,k) = tempvec
       end do
     end do
   end do
   
   do k = 1, n
     do j = 1, n
       do i = 1, n
         call mv_prod(trafo2, ham2(i,j,:,k), tempvec)
         ham2(i,j,:,k) = tempvec
       end do
     end do
   end do
   
   do k = 1, n
     do j = 1, n
       do i = 1, n
         call mv_prod(trafo2, ham2(i,j,k,:), tempvec)
         ham2(i,j,k,:) = tempvec
       end do
     end do
   end do
end subroutine transform_2e_ham_2

subroutine transform_2e_ham_4(ham2, trafo1, trafo2, trafo3, trafo4)
  real(wp), allocatable, intent(inout) :: ham2(:,:,:,:)
  real(wp), intent(in)                 :: trafo1(:,:)
  real(wp), intent(in)                 :: trafo2(:,:)
  real(wp), intent(in)                 :: trafo3(:,:)
  real(wp), intent(in)                 :: trafo4(:,:)
  !local variables
  integer                              :: n1, n2, n3, n4
  integer                              :: m1, m2, m3, m4
  integer                              :: i, j, k
  real(wp), allocatable                :: temp1(:,:,:,:), tempvec1(:)
  real(wp), allocatable                :: temp2(:,:,:,:), tempvec2(:)
  real(wp), allocatable                :: temp3(:,:,:,:), tempvec3(:)
  real(wp), allocatable                :: temp4(:,:,:,:), tempvec4(:)
  character(len=*), parameter          :: proc_name = "transform_ham2"

  if (use_profiling) call start_profiling(proc_name)

  n1 = size(trafo1, 1)
  n2 = size(trafo2, 1)
  n3 = size(trafo3, 1)
  n4 = size(trafo4, 1)

  m1 = size(trafo1, 2)
  m2 = size(trafo2, 2)
  m3 = size(trafo3, 2)
  m4 = size(trafo4, 2)

  allocate(temp1(m1,n2,n3,n4), tempvec1(m1))
  do k = 1, n4
    do j = 1, n3
      do i = 1, n2
        call mv_prod(trafo1, ham2(:,i,j,k), tempvec1)
        temp1(:,i,j,k) = tempvec1
      end do
    end do
  end do
   
  allocate(temp2(m1,m2,n3,n4), tempvec2(m2))
  do k = 1, n4
    do j = 1, n3
      do i = 1, m1
        call mv_prod(trafo2, temp1(i,:,j,k), tempvec2)
        temp2(i,:,j,k) = tempvec2
      end do
    end do
  end do
  
  allocate(temp3(m1,m2,m3,n4), tempvec3(m3))
  do k = 1, n4
    do j = 1, m2
      do i = 1, m1
        call mv_prod(trafo3, temp2(i,j,:,k), tempvec3)
        temp3(i,j,:,k) = tempvec3
      end do
    end do
  end do
  
  allocate(temp4(m1,m2,m3,m4), tempvec4(m4))
  do k = 1, m3
    do j = 1, m2
      do i = 1, m1
        call mv_prod(trafo4, temp3(i,j,k,:), tempvec4)
        temp4(i,j,k,:) = tempvec4
      end do
    end do
  end do
  
  call move_alloc(temp4, ham2)

  if (use_profiling) call end_profiling(proc_name)
end subroutine transform_2e_ham_4


!******************************************************************************** 
!
! Help routine for transformation of two-body integrals
!
!******************************************************************************** 
subroutine mv_prod(A, x, y)
  real(wp), intent(in)    :: A(:,:)
  real(wp), intent(in)    :: x(:)
  real(wp), intent(inout) :: y(:)
  !local variables
  integer                 :: n, m

  n = size(A, 1)
  m = size(A, 2)

  call gemv("t", n, m, 1.0_wp, A, n, x, 1, 0.0_wp, y, 1)
end subroutine mv_prod


!******************************************************************************** 
!
!       Transform one particular Coulomb integral with indices p q r s
!
!       get_integral_trafo     - T[C1 C2, V]
!       get_integral_dec       - [ij|kl] = sum_g Lg(ij) Lg(kl)
!       get_integral_dec_trafo - T[C1 C2, L]
!
!******************************************************************************** 
function get_integral_trafo(ham2, trafo1, trafo2, p, q, r, s) result(integral)
  real(wp), intent(in) :: ham2(:,:,:,:,:)
  real(wp), intent(in) :: trafo1(:,:)
  real(wp), intent(in) :: trafo2(:,:)
  integer, intent(in)  :: p, q, r, s
  real(wp)             :: integral
  !local variables
  integer              :: i, j, k, l, n
  ! 
  integral = 0.0_wp
  n = size(ham2, 1)
  !
  do l = 1, n
    do k = 1, n
      do j = 1, n
        do i = 1, n
            integral = integral + trafo1(i,p)*trafo1(j,q)*trafo2(k,r)*trafo2(l,s)*ham2(i,j,k,l,1)
        end do
      end do
    end do
  end do
end function get_integral_trafo

function get_integral_dec(ham2_dec, p, q, spin1, spin2) result(integral)
  real(wp), intent(in)         :: ham2_dec(:,:,:)
  integer, intent(in)          :: p, q
  integer, optional,intent(in) :: spin1, spin2
  real(wp)                     :: integral
  !local variables
  integer                      :: spin1_, spin2_, g
  
  spin1_ = 1
  spin2_ = 1
  if (present(spin1)) spin1_ = spin1
  if (present(spin2)) spin2_ = spin2
  integral = 0.0_wp
  
  do g = 1, size(ham2_dec, 2)
   integral = integral + ham2_dec(p,g,spin1_) * ham2_dec(q,g,spin2_) 
  end do
end function get_integral_dec

function get_integral_dec_trafo(ham2_dec, trafo1, trafo2, p, q) result(integral)
  real(wp), intent(in)         :: ham2_dec(:,:,:)
  real(wp), intent(in)         :: trafo1(:,:)
  real(wp), intent(in)         :: trafo2(:,:)
  integer, intent(in)          :: p, q
  real(wp)                     :: integral
  !local variables
  integer                      :: a, b, c, d
  integer                      :: i, j, g, n, ng
  real(wp), allocatable        :: temp1(:), temp2(:), tempmat(:,:)
  
  n = size(trafo1, 1)
  ng = size(ham2_dec, 2)
  allocate(temp1(ng), temp2(ng), tempmat(n,n))
  temp1 = 0.0_wp
  temp2 = 0.0_wp
  
  call get_index(n, n, a, b, p, 2)
  call get_index(n, n, c, d, q, 2)

  do g = 1, ng
    tempmat = reshape(ham2_dec(:,g,1), shape=[n,n])
    do j = 1, n
      do i = 1, n
        temp1(g) = temp1(g) + trafo1(i,a) * trafo1(j,b) * tempmat(i,j)
        temp2(g) = temp2(g) + trafo2(i,c) * trafo2(j,d) * tempmat(i,j)
      end do
    end do
  end do
  
  integral = sum(temp1 * temp2)
end function get_integral_dec_trafo


!********************************************************************************
!       subroutine that transform one index of  canonical integrals to obtain       
!               integrals, where one index is local.                     
!               INPUT:          N  - Number of canonical orbitals                        
!                               M  - Number of local orbitals with M<=N                 
!               inout:          HAM_1(1:M,:)  - Part corresponding to one local index   
!                               INT_2(1:M,:,:,:) - same for two-integals                
!               The last index of an array will be transformed                          
!********************************************************************************
subroutine GEN_TRAFO_INT_1(N,M,PHI,HAM_1,INT_2)
  integer,intent(IN):: N, M
  real(wp),DIMENSION(M,N),intent(IN):: PHI
  real(wp),DIMENSION(N,N),intent(inout):: HAM_1
  real(wp),DIMENSION(N,N,N,N),intent(inout):: INT_2

  !local variables
  integer:: i,j,k
  real(wp),DIMENSION(M):: tempvec

  !Transform one-electron integrals
  !  j -> beta
  do i=1,N
    call mv_prod(N,M,PHI,HAM_1(i,:),tempvec)
    HAM_1(i,1:M) = tempvec
  end do

  !Transform two-electron integrals
  !  l -> delta
  do k=1,N
    do j=1,N
      do i=1,N
        call mv_prod(N,M,PHI,INT_2(i,j,k,:),tempvec)
        INT_2(i,j,k,1:M) = tempvec
      end do
    end do
  end do

contains
  subroutine mv_prod(n, m, A, x, y)
    integer, intent(in)     :: n
    integer, intent(in)     :: m
    real(wp), intent(in)    :: A(:,:)
    real(wp), intent(in)    :: x(:)
    real(wp), intent(inout) :: y(:)
    !local variables
    real(wp)                :: alpha, beta
    !
    alpha = 1.0_wp
    beta = 0.0_wp
    call gemv("n", m, n, alpha, A, m, x, 1, beta, y, 1)
  end subroutine mv_prod
end subroutine GEN_TRAFO_INT_1

end module orth
