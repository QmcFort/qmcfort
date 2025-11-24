! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module gauss_base

use constants
use file_handle, only: FileHandle
use standalone
use qmcfort_io
use qmcfort_pos
use hamilton_vars
use string, only: lowercase
use mpi, only: comm_world

implicit none

!******************************************************************************** 
!
! Stores Gaussian basis set for specific geometry
!
!******************************************************************************** 
type gauss_basis_set
  type(gauss_shell), allocatable :: shell(:)
contains
  generic, public                :: get_start_indx => get_start_indx_1, get_start_indx_2
  generic, public                :: get_shells_from_indx  => get_shells_from_indx_1, get_shells_from_indx_2
  procedure                      :: get_nshells 
  procedure                      :: get_norbitals
  procedure                      :: read => read_bs
  procedure                      :: read_mpi => read_bs_mpi
  procedure                      :: add_shell => add_basis_set_shell
  procedure, private             :: get_start_indx_1, get_start_indx_2
  procedure, private             :: get_shells_from_indx_1, get_shells_from_indx_2
end type gauss_basis_set


!******************************************************************************** 
!
! Stores one Gaussian shell that is uniquely determined by:
!
!    n          - principal quantum number
!    l          - angular quantum number
!    cite       - 3D center of the orbitals
!
! Contains (2l+1) orbitals lm with m = -l, -l+1 ... 0 ... l-1, l
!
! Orbital representation:
!
!    \chi_{lm}(exp, cite) = \sum_{d}^{cont_len} cont_coeff_{d} G_{lm}(exp_d, cite)
!
! where 
! 
!    G_{lm}  is spherical-harmonic GTO function
!
!******************************************************************************** 
type gauss_shell
  character(len=charlen) :: name
  character(len=charlen) :: atom
  integer                :: l
  real(wp)               :: cite(3)
  real(wp), allocatable  :: cont_coeff(:)
  real(wp), allocatable  :: cont_exp(:)
  real(wp), allocatable  :: trafo(:,:)
contains
  generic, public        :: assignment(=) => assign_gauss_shell
  procedure              :: get_shell_size
  procedure, private     :: assign_gauss_shell
end type gauss_shell

integer, allocatable :: cart_indx_matrix(:,:,:)
integer, allocatable :: cart_ijk_matrix(:,:,:)

contains

!******************************************************************************** 
!
! Determine number of shells in the basis set
!
!******************************************************************************** 
function get_nshells(bs) result(nshells)
  class(gauss_basis_set), intent(in) :: bs
  integer                            :: nshells

  nshells = 0
  if (allocated(bs%shell)) nshells = size(bs%shell)
end function get_nshells


!******************************************************************************** 
!
! Determine number of orbitals in the basis set
!
!******************************************************************************** 
function get_norbitals(bs) result(norbitals)
  class(gauss_basis_set), intent(in) :: bs
  integer                            :: norbitals
  !local variables
  integer                            :: i

  norbitals = 0

  do i = 1, size(bs%shell)
    norbitals = norbitals + bs%shell(i)%get_shell_size()
  end do
end function get_norbitals


!******************************************************************************** 
!
! Determine start orbital index for the a-th shell in basis set
!
!******************************************************************************** 
function get_start_indx_1(bs, a) result(indx)
  class(gauss_basis_set), intent(in) :: bs
  integer, intent(in)                :: a
  integer                            :: indx
  !local variables
  integer                            :: i

  indx = 0

  do i = 1, a-1
    indx = indx + bs%shell(i)%get_shell_size()
  end do
end function get_start_indx_1


!******************************************************************************** 
!
! Determine start orbital index for the pair shells (shell_a shell_b | ... )
!
!******************************************************************************** 
function get_start_indx_2(bs, a, b) result(indx)
  class(gauss_basis_set), intent(in) :: bs
  integer, intent(in)                :: a
  integer, intent(in)                :: b
  integer                            :: indx
  !local variables
  integer                            :: i, n, n_a, n_b

  n = bs%get_norbitals()

  n_a = 0
  do i = 1, a-1
    n_a = n_a + bs%shell(i)%get_shell_size()
  end do

  n_b = 0
  do i = 1, b-1
    n_b = n_b + bs%shell(i)%get_shell_size()
  end do

  indx = n*n_b + n_a*bs%shell(b)%get_shell_size()
end function get_start_indx_2


!******************************************************************************** 
!
! Obtain shell indices a, b for the orbital pair index in range 1..N^2
!
!******************************************************************************** 
subroutine get_shells_from_indx_1(bs, indx, a, b)
  class(gauss_basis_set), intent(in) :: bs
  integer, intent(in)                :: indx
  integer, intent(out)               :: a
  integer, intent(out)               :: b
  !local variables
  integer                            :: aa, bb, size_a, size_b, indx_
  
  indx_ = 0 

  do bb = 1, bs%get_nshells()
    size_b = bs%shell(bb)%get_shell_size()
    do aa = 1, bs%get_nshells()
      size_a = bs%shell(aa)%get_shell_size()
      indx_ = indx_ + size_a*size_b
      if (indx_ >= indx) then
        a = aa
        b = bb
        return
      end if
    end do
  end do
end subroutine get_shells_from_indx_1

!******************************************************************************** 
!
! Obtain index within shell form the full orbital pair index in range 1..N^2
!
!******************************************************************************** 
subroutine get_shells_from_indx_2(bs, indx, indx_sh)
  class(gauss_basis_set), intent(in) :: bs
  integer, intent(in)                :: indx
  integer, intent(out)               :: indx_sh
  !local variables
  integer                            :: aa, bb, size_a, size_b, indx_
  
  indx_ = 0 

  do bb = 1, bs%get_nshells()
    size_b = bs%shell(bb)%get_shell_size()
    do aa = 1, bs%get_nshells()
      size_a = bs%shell(aa)%get_shell_size()
      indx_ = indx_ + size_a*size_b
      if (indx_ >= indx) then
        indx_sh = indx - (indx_ - size_a*size_b)
        return
      end if
    end do
  end do
end subroutine get_shells_from_indx_2


!******************************************************************************** 
!
! Reorder Cholesky vectors from shell-shell ordering to orbital-orbital ordering
!
! Example:
!    Consider basis set that contains 2 shells and 4 orbitals in total:
!        [AB] := [A1,A2,B1,B2]
!
!  shell-shell:        orbital-orbital:
!     A1 A1                 A1 A1
!     A2 A1                 A2 A1
!     A1 A2                 B1 A1
!     A2 A2                 B2 A1
!     B1 A1                 A1 A2
!     B2 A1                 A2 A2
!     B1 A2                 B1 A2
!     B2 A2                 B2 A2
!     ...                   ...
!
!******************************************************************************** 
subroutine shell_reorder(bs, h2_gres)
  type(gauss_basis_set), intent(in) :: bs
  real(wp), intent(inout)           :: h2_gres(:,:)
  !local variables
  integer                           :: i, g
  integer, allocatable              :: indx_array(:)
  real(wp), allocatable             :: temp(:)

  call get_indx_array(bs, indx_array)
  allocate(temp(size(h2_gres,1)))

  do g = 1, size(h2_gres,2)
    do i = 1, size(h2_gres,1)
      temp(i) = h2_gres(indx_array(i),g)
    end do
    h2_gres(:,g) = temp
  end do
contains 
  subroutine get_indx_array(bs, indx_array)
    type(gauss_basis_set), intent(in) :: bs
    integer, allocatable, intent(out) :: indx_array(:)
    !local variables
    integer                           :: a, b, aa, bb, a_full, b_full, a_start, b_start
    integer                           :: n, indx, indx_

    n = bs%get_norbitals()
    allocate(indx_array(n*n))

    indx_ = 0
    do b = 1, size(bs%shell)
      b_start = bs%get_start_indx(b)
      do a = 1, size(bs%shell)
        a_start = bs%get_start_indx(a)
        do bb = 1, bs%shell(b)%get_shell_size()
          b_full = b_start + bb
          do aa = 1, bs%shell(a)%get_shell_size()
            a_full = a_start + aa
            indx = (b_full-1)*n + a_full
            indx_ = indx_ + 1
            indx_array(indx) = indx_
          end do
        end do
      end do
    end do
  end subroutine get_indx_array
end subroutine


!******************************************************************************** 
!
! MPI wrapper to read Gaussian basis sets
!
!******************************************************************************** 
subroutine read_bs_mpi(basis_set, ham_des, struc)
  class(gauss_basis_set), intent(inout) :: basis_set
  type(hamil_descriptor), intent(inout) :: ham_des
  type(structure), intent(in)           :: struc
  !local variables
  integer                               :: proc
  
  do proc = 0, comm_world%mpisize-1
    if (proc == comm_world%mpirank) call basis_set%read(ham_des, struc)
    call comm_world%barrier()
  end do
end subroutine read_bs_mpi


!******************************************************************************** 
!
! Read Gaussian basis set
!
!******************************************************************************** 
subroutine read_bs(basis_set, ham_des, struc)
  class(gauss_basis_set), intent(inout) :: basis_set
  type(hamil_descriptor), intent(inout) :: ham_des
  type(Structure), intent(in)           :: struc
  !local variables
  integer                               :: n, q, v, ierror
  integer                               :: nspecies, nshells, cont_len
  character(len=atomc)                  :: atom_symb
  character(len=3)                      :: shell_name
  character(len=linec)                  :: line
  real(wp)                              :: zeta
  real(wp), allocatable                 :: cont_exp(:), cont_coeff(:), cont_coeff_(:)
  real(wp), allocatable                 :: coords(:,:)
  type(FileHandle)                      :: fh

  fh = FileHandle("basis_set")
  call fh%open(status="old", action="read")

  nspecies = 0
  do
    read(fh%funit,"(a)",iostat=ierror) line
    if (ierror /= 0) exit

    !cycle if inappropriate line
    if (len_trim(line) <= 1) cycle
    if (line(1:1)=="!" .or. line(1:1)=="#" .or. line(1:1)=="*") cycle
    
    read(line,*,iostat=ierror) shell_name , cont_len , zeta
    if (ierror /= 0) then
      read(line,*,iostat=ierror) atom_symb, nshells
      if (ierror == 0) nspecies = nspecies + 1
      cycle
    else 
      if (allocated(cont_exp))  deallocate(cont_exp)
      if (allocated(cont_coeff))  deallocate(cont_coeff)
      if (allocated(cont_coeff_))  deallocate(cont_coeff_)
      allocate(cont_exp(cont_len), cont_coeff(cont_len), cont_coeff_(cont_len))

      v = 1
      do 
        read(fh%funit,"(a)",iostat=ierror) line

        !cycle if inappropriate line
        if (len_trim(line) <= 1) cycle
        if (line(1:1)=="!" .or. line(1:1)=="#" .or. line(1:1)=="*") cycle

        read(line,*,iostat=ierror) cont_exp(v) , cont_coeff(v), cont_coeff_(v)
        if (ierror == 0) then
          v = v + 1
        else 
          read(line,*,iostat=ierror) cont_exp(v) , cont_coeff(v)
          if (ierror == 0) v = v + 1
        end if

        if (v-1 == cont_len) exit
      end do

      coords = struc%get_coords_for_species(atom_symb)
      do q = 1, size(coords,2)
        if (lowercase(shell_name) == "sp") then 
          call basis_set%add_shell("s", atom_symb, coords(:,q), zeta, cont_exp, cont_coeff)
          call basis_set%add_shell("p", atom_symb, coords(:,q), zeta, cont_exp, cont_coeff_)
        else 
          call basis_set%add_shell(shell_name, atom_symb, coords(:,q), zeta, cont_exp, cont_coeff)
        end if
      end do
    end if
  end do

  if (nspecies /= struc%get_nspecies()) then
    if (comm_world%mpirank == 0) write(*,*) "basis set not read for all species found in Structure object"
    stop
  end if

  n = basis_set%get_norbitals()
  ham_des%n = n
  if (n == 0) then
    if (comm_world%mpirank == 0) write(*,*) "fatal error: number of read gaussians is 0 - PROGRAM STOPS NOW"
    stop
  end if

  if (comm_world%mpirank == 0) write(*,*) "file BASIS_SET read successfully"
  call fh%close()
end subroutine read_bs


!******************************************************************************** 
!
! Generate basis_set file
!
!******************************************************************************** 
subroutine generate_basis_set_file(basis_set, species)
  character(len=*), intent(in) :: basis_set
  character(len=*), intent(in) :: species
  !local
  character(len=120)           :: command_str

  command_str = TOOLS_DIR // "/makebasis.sh" // " " // trim(basis_set) // " " // trim(species)

  if (comm_world%mpirank == 0) call system(command_str)
end subroutine generate_basis_set_file


!******************************************************************************** 
!
! Add additional shell into the basis_set
!
! Shel is described via:
!    shell_name     -     S, P, D, F, ...
!    atom_symb      -     species name
!    cite           -     Gaussian center
!    zeta           -     scaling parameter for primitive contraction exponents
!    cont_exp       -     contraction exponents from the file 
!    cont_coeff     -     contract
!
!******************************************************************************** 
subroutine add_basis_set_shell(basis_set, shell_name, atom_symb, cite, zeta, cont_exp, cont_coeff)
  class(gauss_basis_set), intent(inout) :: basis_set
  character(len=*), intent(in)          :: shell_name
  character(len=*), intent(in)          :: atom_symb
  real(wp), intent(in)                  :: cite(:)
  real(wp), intent(in)                  :: zeta
  real(wp), intent(in)                  :: cont_exp(:)
  real(wp), intent(in)                  :: cont_coeff(:)
  !local
  integer                               :: i, n
  type(gauss_shell)                     :: shell
  type(gauss_shell), allocatable        :: new_shells(:)

  n = basis_set%get_nshells()
  allocate(new_shells(n+1))
  do i = 1, n
    new_shells(i) = basis_set%shell(i)
  end do

  shell%l = get_l_number(shell_name)
  shell%name = shell_name
  shell%atom = atom_symb
  shell%cite = cite
  allocate(shell%cont_exp(size(cont_exp)))
  allocate(shell%cont_coeff(size(cont_coeff)))
  shell%cont_exp = cont_exp * zeta**2
  shell%cont_coeff = cont_coeff * gto_norm(shell%cont_exp, shell%l)
  shell%trafo = harm_cart_trafo(shell%l)

  new_shells(n+1) = shell
  call move_alloc(new_shells, basis_set%shell)
end subroutine add_basis_set_shell


!******************************************************************************** 
!
! Assignment operator for the Gaussian shell
!
!******************************************************************************** 
subroutine assign_gauss_shell(to, from)
  class(gauss_shell), intent(inout) :: to
  type(gauss_shell), intent(in)     :: from

  to%l = from%l
  to%name = from%name
  to%atom = from%atom
  to%cite = from%cite

  allocate(to%cont_coeff, source=from%cont_coeff)
  allocate(to%cont_exp, source=from%cont_exp)
  allocate(to%trafo, source=from%trafo)
end subroutine assign_gauss_shell


!******************************************************************************** 
!
! Determine shell size (number of orbitals in the shell)
!
!******************************************************************************** 
function get_shell_size(shell) result(shell_size)
  class(gauss_shell), intent(in) :: shell
  integer                        :: shell_size

  shell_size = 2*shell%l + 1
end function get_shell_size


!******************************************************************************** 
!
! Norm of the spherical-harmonic Gaussian orbital with quantum number l
!
!******************************************************************************** 
elemental real(wp) function gto_norm(alpha, l)
  real(wp), intent(in) :: alpha
  integer, intent(in)  :: l

  gto_norm = (2.0*alpha/pi)**(0.75) * (4.0*alpha)**(l/2.0) / sqrt(real(double_factorial(2*l-1),wp))
end function gto_norm


!******************************************************************************** 
!
! Transformation matrix from cartesian to spherical-harmonics Gaussian basis
!
!
!  G_{lm} = N_{lm} \sum_{t=0}^{tmax} \sum_{u=0}^{t} \sum_{v=vm}^{vmax}  x
!           C_{tuv}^{lm} G_{xyz}
!
!    vm = 0 if m>=0; 0.5 if m<0
!    tmax = floor((l-|m|)/2)
!    vmax = floor(|m|/2 - vm) + vm
!    x = 2t + |m| = 2(u+v)
!    y = 2(u+v)
!    z = l - 2t - |m|
!
!    N_{lm} = 1/(2^{|m|} l!) sqrt(2 (l+|m|)! (l-|m|)! / 2**{delta_{0m}})
!    C_{tuv}^{lm} = (-1)^{t+v-vm} (0.25)^{t} bin(l,t) bin(l-t,|m|+t) bin(t,u) bin(|m|, 2v)
!
!******************************************************************************** 
function harm_cart_trafo(l) result(trafo)
  integer, intent(in)   :: l
  real(wp), allocatable :: trafo(:,:)
  !local variables
  integer               :: m, harm_size, cart_size, harm_indx, cart_indx
  integer               :: t, u, v, tmax, vmax, x, y, z
  real(wp)              :: vm, vv, norm

  harm_size = 2*l + 1
  cart_size = (l+1) * (l+2) / 2
  allocate(trafo(cart_size, harm_size))
  trafo = 0.0_wp

  do m = l, -l, -1
    harm_indx = get_harm_indx(l, m)

    if (m >= 0) then
      vm = 0.0_wp
    else 
      vm = 0.5_wp
    end if

    tmax = (l-abs(m)) / 2
    vmax = floor(abs(m)/2.0 - vm)

    norm = Nlm(l, m)

    do t = 0, tmax
      do u = 0, t
        do v = 0, vmax
          vv = v + vm

          x = 2*t + abs(m) - floor(2*(u+vv))
          y = floor(2*(u+vv))
          z = l - 2*t - abs(m)

          cart_indx = get_cart_indx(x, y, z)

          trafo(cart_indx, harm_indx) = norm * Clm_tuv(l, m, t, u, vv)
        end do
      end do
    end do
  end do
contains 
  function Nlm(l, m) result(norm)
    integer, intent(in) :: l
    integer, intent(in) :: m
    real(wp)            :: norm

    norm = sqrt(2.0_wp*factorial(l+abs(m))*factorial(l-abs(m))) / (2**(abs(m))*factorial(l))
    if (m==0) norm = norm / sqrt(2.0_wp)
  end function Nlm

  function Clm_tuv(l, m, t, u, v) result(coeff)
    integer, intent(in)  :: l
    integer, intent(in)  :: m
    integer, intent(in)  :: t
    integer, intent(in)  :: u
    real(wp), intent(in) :: v
    real(wp)             :: coeff
    !local variables
    real(wp)             :: vm

    if (m >= 0) then
      vm = 0.0_wp
    else 
      vm = 0.5_wp
    end if 

    coeff = 0.25_wp**(t) * binom(l, t) * binom(l-t, abs(m)+t) * binom(t, u) * binom(abs(m), floor(2*v))
    if (mod(floor(t+v-vm), 2) == 1) coeff = -coeff
  end function Clm_tuv
end function harm_cart_trafo


!******************************************************************************** 
!
! Transform shell name (s, p, d, f...) to the angular quantum number l
!
!******************************************************************************** 
function get_l_number(shell_name) result(l)
  character(len=*), intent(in) :: shell_name
  integer                      :: l

  select case (trim(shell_name))
    case ("s", "S")
      l = 0
    case ("p", "P")
      l = 1
    case ("d", "D")
      l = 2
    case ("f", "F")
      l = 3
    case ("g", "G")
      l = 4
    case ("h", "H")
      l = 5
    case ("i", "I")
      l = 6
    case ("j", "J")
      l = 7
    case ("k", "K")
      l = 8
    case ("l", "L")
      l = 9
    case ("m", "M")
      l = 10
    case default
      if (comm_world%mpirank == 0) write(*,*)"Shell is not recognised - reading of the BASIS_SET failed - PROGRAM STOPS NOW"
      call exit
  end select
end function get_l_number
  

!******************************************************************************** 
!
! Determines the ordering of spherical-harmonics GTOs within l shell
!
!******************************************************************************** 
function get_harm_indx(l, m) result(indx)
  integer, intent(in) :: l
  integer, intent(in) :: m
  integer             :: indx

  indx = m + l + 1
  indx = 2*l + 2 - indx
end function get_harm_indx


!******************************************************************************** 
!
! Determines the ordering of cartesian GTOs within l shell
!
!******************************************************************************** 
function get_cart_indx(i, j, k) result(indx)
  integer, intent(in) :: i
  integer, intent(in) :: j
  integer, intent(in) :: k
  integer             :: indx

  if (.not. allocated(cart_indx_matrix)) call get_cart_indx_matrix() 
  if (i+j+k > size(cart_indx_matrix,1)-1) call get_cart_indx_matrix(i+j+k)

  indx = cart_indx_matrix(i,j,k)
end function get_cart_indx


!******************************************************************************** 
!
! Determine ijk of the Cartesian Gaussian G_{ijk} from the ordering index indx
!
!******************************************************************************** 
subroutine get_cart_ijk(l, indx, i, j, k)
  integer, intent(in)  :: l
  integer, intent(in)  :: indx
  integer, intent(out) :: i
  integer, intent(out) :: j
  integer, intent(out) :: k

  if (.not. allocated(cart_ijk_matrix)) call get_cart_ijk_matrix() 
  if (l > size(cart_ijk_matrix,1)-1) call get_cart_ijk_matrix(l)

  i = cart_ijk_matrix(l, indx, 1)
  j = cart_ijk_matrix(l, indx, 2)
  k = cart_ijk_matrix(l, indx, 3)
end subroutine get_cart_ijk


!******************************************************************************** 
!
! cart_indx_matrix yields ordering index of the cartesian Gaussian G_{ijk}
!
!    indx = cart_order_matrix(i,j,k)
!
! indx  is obtained via get_cart_indx(i,j,k) function
!
!******************************************************************************** 
subroutine get_cart_indx_matrix(lmax)
  integer, optional, intent(in) :: lmax
  !local variables
  integer                       :: l, lmax_
  integer                       :: k, p, q, r

  !$omp single
  lmax_ = 14
  if (present(lmax)) lmax_ = lmax
  if (allocated(cart_indx_matrix)) deallocate(cart_indx_matrix)
  allocate(cart_indx_matrix(0:lmax_, 0:lmax_, 0:lmax_))
  cart_indx_matrix = 0

  do l = 0, lmax_
    k = 0
    do r = 0, l
      do q = 0, l
        do p = 0, l
          if (p+q+r == l) then
            k = k + 1
            cart_indx_matrix(p,q,r) = k 
          end if
        end do
      end do
    end do
  end do
  !$omp end single
end subroutine get_cart_indx_matrix


!******************************************************************************** 
!
! cart_ijk_matrix is used to obtain ijk indices from G_{ijk} notation for a 
! given ordering indx (ordering indx ==> ca)
!
!******************************************************************************** 
subroutine get_cart_ijk_matrix(lmax)
  integer, optional, intent(in) :: lmax
  !local variables
  integer                       :: l, lmax_
  integer                       :: k, kmax
  integer                       :: p, q, r

  !$omp single
  lmax_ = 14
  if (present(lmax)) lmax_ = lmax
  kmax = (lmax_+1) * (lmax_+2) / 2
  if (allocated(cart_ijk_matrix)) deallocate(cart_ijk_matrix)
  allocate(cart_ijk_matrix(0:lmax_, kmax, 3))
  cart_ijk_matrix = 0

  do l = 0, lmax_
    k = 0
    do r = 0, l
      do q = 0, l
        do p = 0, l
          if (p+q+r == l) then
            k = k + 1
            cart_ijk_matrix(l,k,:) = [p, q, r]
          end if
        end do
      end do
    end do
  end do
  !$omp end single
end subroutine get_cart_ijk_matrix

end module gauss_base
