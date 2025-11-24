! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module slater
!******************************************************************************** 
!
! This module implements basic types for the manipulation with  orthogonal 
! Slater determinants
!
!    slater_str type
!    slater_det type
!
!******************************************************************************** 

use constants
use string

integer, parameter :: strkind = i8
integer, parameter :: strsize = i8size

type slater_str
  integer(strkind) :: str
contains
  generic, public    :: excite => excite_str_1, excite_str_2 
  procedure, public  :: len => str_len
  procedure, public  :: nel => str_nel
  procedure, public  :: nmax => str_nmax
  procedure, public  :: occ => str_occupation
  procedure, public  :: exc_deg => excitation_degree_string
  procedure, public  :: content => str_content
  procedure, public  :: diff => str_diff
  procedure, public  :: ovlp => str_ovlp
  procedure, public  :: phase => str_phase
  procedure, private :: excite_str_1, excite_str_2
end type slater_str

type slater_det
  integer(strkind) :: str(2)
contains
  generic, public    :: excite => excite_det_1, excite_det_2
  procedure, public  :: len => det_len
  procedure, public  :: nel => det_nel
  procedure, public  :: nmax => det_nmax
  procedure, public  :: occ => det_occupation
  procedure, public  :: exc_deg => excitation_degree_det
  procedure, public  :: content => det_content
  procedure, public  :: diff => det_diff
  procedure, public  :: phase => det_phase
  procedure, public  :: ovlp => det_ovlp
  procedure, public  :: read_from_char => read_det_from_char
  procedure, private :: excite_det_1, excite_det_2
end type slater_det

interface slater_str
  procedure :: init_slater_str
end interface 

interface slater_det
  procedure :: init_slater_det_int, init_slater_det_str
end interface 

contains

!******************************************************************************** 
!
! Slater string constructor
!
!******************************************************************************** 
type(slater_str) function init_slater_str(intstr)
  integer(strkind), intent(in) :: intstr
  
  init_slater_str%str = intstr
end function init_slater_str


!******************************************************************************** 
!
! Returns number of available bits, i.e. maximal number of electrons 
! that can be stored in the Slater string
!
!******************************************************************************** 
function str_len(str) result(len)
  class(slater_str), intent(in) :: str
  integer                       :: len

  len = storage_size(str%str)
end function str_len


!******************************************************************************** 
!
! Counts number of electrons stored in Slater determinant
!
!******************************************************************************** 
function str_nel(str) result(nel)
  class(slater_str), intent(in) :: str
  integer                       :: nel

  nel = popcnt(str%str)
end function str_nel


!******************************************************************************** 
!
! Returs highest occupied orbital in Slater string
!
!******************************************************************************** 
function str_nmax(str) result(n)
  class(slater_str), intent(in) :: str
  integer                       :: n 

  n = str%len() - leadz(str%str)
end function str_nmax


!******************************************************************************** 
!
! Get occupation of the orbital in Slater string
!
! returns either 0 or 1
! valid oribtals index: from 1 to strsize
!
!******************************************************************************** 
function str_occupation(str, i) result(occ)
  class(slater_str), intent(in) :: str
  integer, intent(in)           :: i
  integer                       :: occ

  if (btest(str%str, i-1)) then
    occ = 1
  else 
    occ = 0
  end if
end function str_occupation


!******************************************************************************** 
!
! Determine excitation degree for two Slater strings
!
! exciation degrer = number of single excitation operators needed to apply
! to transform str1 into str2:
!
!    str2 = T str1,  with T = E_{i}^{a} E_{j}^{b} ...
!
!******************************************************************************** 
function excitation_degree_string(str1, str2) result(ndiff)
  class(slater_str), intent(in) :: str1
  type(slater_str), intent(in)  :: str2
  integer                       :: ndiff
  !local variables
  integer(strkind)              :: diff_

  diff_ = ieor(str1%str, str2%str)
  ndiff = popcnt(diff_) / 2
end function excitation_degree_string


!******************************************************************************** 
!
! Apply single excitation operator E_{i}^{a} onto Slater string str 
!
!******************************************************************************** 
function excite_str_1(str, i, a) result(str_)
  class(slater_str), intent(in) :: str
  integer, intent(in)           :: i
  integer, intent(in)           :: a
  type(slater_str)              :: str_
  
  str_ = str
  str_%str = ibset(ibclr(str%str, i-1), a-1)
end function excite_str_1

!******************************************************************************** 
!
! Apply double excitation operator E_{ij}^{ab} onto Slater string str 
!
!******************************************************************************** 
function excite_str_2(str, i, a, j, b) result(str_)
  class(slater_str), intent(in) :: str
  integer, intent(in)           :: i
  integer, intent(in)           :: j
  integer, intent(in)           :: a
  integer, intent(in)           :: b
  type(slater_str)              :: str_
  
  str_ = str%excite(i, a)
  str_ = str_%excite(j, b)
end function excite_str_2


!******************************************************************************** 
!
! Get list of the occupied orbitals in the Slater string
!
!    example:
!        integer rep :  199
!        bit rep     :  1 1 1 0 0 0 1 1 0 
!        content     :  [1 2 3 7 8]
!   
!******************************************************************************** 
function str_content(str) result(ind)
  class(slater_str), intent(in) :: str
  integer, allocatable          :: ind(:)
  !local variables
  integer                       :: i, indx, spin, nelect
  
  nelect = str%nel()

  allocate(ind(nelect))
  
  indx = 0
  do i = 1, str%nmax()
    if (btest(str%str, i-1)) then
      indx = indx + 1
      ind(indx) = i
    end if
  end do
end function str_content


!******************************************************************************** 
!
! Obtain indices of excitation operator E_{ijkl..}^{abcd...}
! that transforms str1 to str2
!
! Following format is used: diff(2,ndiff)
!
!    ndiff = i j k l ....
!            a b c d ....           
!
!******************************************************************************** 
function str_diff(str1, str2) result(diff)
  class(slater_str), intent(in) :: str1
  type(slater_str), intent(in)  :: str2
  integer, allocatable          :: diff(:,:)
  !local variables
  integer                       :: i, a, k, ndiff
  integer(kind=strkind)         :: diff_

  ndiff = str1%exc_deg(str2)

  if (ndiff == 0) then
    if (allocated(diff)) deallocate(diff)
    allocate(diff(2,ndiff))
    return
  else if (ndiff > 2) then
    if (allocated(diff)) deallocate(diff)
    allocate(diff(2,0))
    return
  else
    if (allocated(diff)) deallocate(diff)
    allocate(diff(2,ndiff))
  end if
  
  i = 0
  a = 0
  diff_ = ieor(str1%str, str2%str)
  do k = 1, strsize
    if (btest(diff_, k-1)) then
      if (btest(str1%str, k-1)) then
        i = i + 1
        diff(1, i) = k
      else
        a = a + 1
        diff(2, a) = k
      end if
    end if
  end do
end function str_diff


!******************************************************************************** 
!
! Estimate phase factor between two Slater strings 
!
!    phase = (-1)^{nperm}
!
! where nperm is number of anticommutation relations needed to bring str1 into str2
!
! if diff array (str_diff) is passed as input: it will be used
! if not: the diff array will be created on the fly
!
! diff array has following form:
!
!    i  j ...
!    a  b ...
!
!******************************************************************************** 
function str_phase(str1, str2, diff) result(phase)
  class(slater_str), intent(in)              :: str1
  type(slater_str), intent(in)               :: str2
  integer, allocatable, optional, intent(in) :: diff(:,:)
  integer                                    :: phase
  !local variables
  integer                                    :: p, exc, ndiff, nperm
  integer, allocatable                       :: diff_(:,:)

  if (present(diff)) then
    allocate(diff_, source=diff)
  else 
    diff_ = str1%diff(str2)
  end if

  nperm = 0
  ndiff = size(diff_, 2)

  do exc = 1, ndiff
    if (diff_(2,exc) > diff_(1,exc)) then
      do p =  diff_(1,exc)+1, diff_(2,exc)-1
        if (btest(str1%str, p-1)) nperm = nperm + 1
      end do
    else 
      do p =  diff_(2,exc)+1, diff_(1,exc)-1
        if (btest(str1%str, p-1)) nperm = nperm + 1
      end do
    end if
  end do

  if (ndiff == 2) then
    if (diff_(2,1)>diff_(1,2) .or. diff_(1,1)>diff_(2,2)) nperm = nperm + 1
  end if

  phase = (-1)**(nperm)
end function str_phase


!******************************************************************************** 
!
! Calculate overlpa between two Slater strings
!
!    <Phi_i|Phi_j> = delta_{ij}
!
!******************************************************************************** 
function str_ovlp(str1, str2) result(ovlp)
  class(slater_str), intent(in) :: str1
  type(slater_str), intent(in)  :: str2
  real(wp)                      :: ovlp

  ovlp = 0.0_wp
  if (str1%str==str2%str) ovlp = 1.0_wp
end function str_ovlp


!******************************************************************************** 
!
! Slater determinant constructor using two integers
!
!******************************************************************************** 
type(slater_det) function init_slater_det_int(intstr_a, intstr_b)
  integer(strkind), intent(in) :: intstr_a
  integer(strkind), intent(in) :: intstr_b
  
  init_slater_det_int%str(1) = intstr_a
  init_slater_det_int%str(2) = intstr_b
end function init_slater_det_int


!******************************************************************************** 
!
! Slater determinant constructor using slater strings
!
!******************************************************************************** 
type(slater_det) function init_slater_det_str(str_a, str_b)
  type(slater_str), intent(in) :: str_a
  type(slater_str), intent(in) :: str_b
  
  init_slater_det_str%str(1) = str_a%str
  init_slater_det_str%str(2) = str_b%str
end function init_slater_det_str


!******************************************************************************** 
!
! Determine excitation degree for two Slater determinants
!
! exciation degree = number of single excitation operators needed to apply
! to transform det1 into det2:
!
!    det2 = T det1,  with T = E_{i}^{a} E_{j}^{b} ...
!
!******************************************************************************** 
function excitation_degree_det(det1, det2) result(ndiff)
  class(slater_det), intent(in) :: det1
  type(slater_det), intent(in)  :: det2
  integer                       :: ndiff
  !local variables
  integer(strkind)              :: diff_

  diff_ = ieor(det1%str(1), det2%str(1))
  ndiff = popcnt(diff_) / 2

  diff_ = ieor(det1%str(2), det2%str(2))
  ndiff = ndiff + popcnt(diff_) / 2
end function excitation_degree_det


!******************************************************************************** 
!
! Returns number of available bits, i.e. maximal number of electrons 
! that can be stored in the Slater determinant
!
!******************************************************************************** 
function det_len(det) result(len)
  class(slater_det), intent(in) :: det
  integer                       :: len

  len = storage_size(det%str(1))
end function det_len


!******************************************************************************** 
!
! Counts number of electrons stored in Slater determinant
!
!******************************************************************************** 
function det_nel(det) result(nel)
  class(slater_det), intent(in) :: det
  integer                       :: nel(2)

  nel(1) = popcnt(det%str(1))
  nel(2) = popcnt(det%str(2))
end function det_nel


!******************************************************************************** 
!
! Returs highest occupied orbital in Slater determinant
!
!******************************************************************************** 
function det_nmax(det) result(n)
  class(slater_det), intent(in) :: det
  integer                       :: n 

  n = det%len() - min(leadz(det%str(1)), leadz(det%str(2)))
end function det_nmax


!******************************************************************************** 
!
! Get occupation of the orbital in Slater determinant
!
! if spin variable is supplied - spin occupation is returned (0 or 1)
! if spin variable is not used - total occupation is returned (0, 1, or 2)
!
! valid oribtals index: from 1 to strsize
! valid spin index: 1 or 2
!
!******************************************************************************** 
function det_occupation(det, i, spin) result(occ)
  class(slater_det), intent(in) :: det
  integer, intent(in)           :: i
  integer, optional, intent(in) :: spin
  integer                       :: occ

  occ = 0
  if (present(spin)) then
    if (btest(det%str(spin), i-1)) occ = occ + 1
  else 
    if (btest(det%str(1), i-1)) occ = occ + 1
    if (btest(det%str(2), i-1)) occ = occ + 1
  end if
end function det_occupation


!******************************************************************************** 
!
! Apply single excitation operator E_{i}^{a}(sigma) onto Slater determinant
!
!******************************************************************************* 
function excite_det_1(det, i, a, spin) result(det_)
  class(slater_det), intent(in) :: det
  integer, intent(in)           :: i
  integer, intent(in)           :: a
  integer, intent(in)           :: spin
  type(slater_det)              :: det_
  
  det_ = det
  det_%str(spin) = ibset(ibclr(det_%str(spin), i-1), a-1)
end function excite_det_1

!******************************************************************************** 
!
! Apply double excitation operator onto Slater determinant
!
!    double exc. operator:  E_{i}^{a}(sigma) E_{j}^{b}(tau)
!
!******************************************************************************** 
function excite_det_2(det, i, a, spin1, j, b, spin2) result(det_)
  class(slater_det), intent(in) :: det
  integer, intent(in)           :: i
  integer, intent(in)           :: j
  integer, intent(in)           :: a
  integer, intent(in)           :: b
  integer, intent(in)           :: spin1
  integer, intent(in)           :: spin2
  type(slater_det)              :: det_
  
  det_ = det%excite(i, a, spin1)
  det_ = det_%excite(j, b, spin2)
end function excite_det_2


!******************************************************************************** 
!
! Get list of the occupied orbitals in the Slater determinant
!
!    example:
!        integer rep    :  199 31
!        bit rep        :  1 1 1 0 0 0 1 1 0 
!                       :  1 1 1 1 1 0 0 0 0
!        content (nex2) :  [1 2 3 7 8 1 2 3 4 5]
!                          [1 1 1 1 1 2 2 2 2 2]
!   
!******************************************************************************** 
function det_content(det) result(ind)
  class(slater_det), intent(in) :: det
  integer, allocatable          :: ind(:,:)
  !local variables
  integer                       :: i, indx, spin, nel(2), nelect
  
  nel = det%nel()
  nelect = sum(nel)

  allocate(ind(nelect, 2))
  
  indx = 0
  do spin = 1, 2
    do i = 1, det%nmax()
      if (btest(det%str(spin), i-1)) then
        indx = indx + 1
        ind(indx,1) = i
        ind(indx,2) = spin
      end if
    end do
  end do
end function det_content


!******************************************************************************** 
!
! Obtain indices of excitation operator E_{ijkl..}^{abcd...}(s1s2s3s4....)
! that transforms det1 to det2
!
! Following format is used: diff(3,ndiff)
!
!    diff =  s1 s2 s3 s4 ...
!            i  j  k  l  ...
!            a  b  c  d  ...           
!
!******************************************************************************** 
function det_diff(det1, det2) result(diff)
  class(slater_det), intent(in) :: det1
  type(slater_det), intent(in)  :: det2
  integer, allocatable          :: diff(:,:)
  !local variables
  integer                       :: i, a, k, ndiff, spin
  integer(kind=strkind)         :: diff_

  ndiff = det1%exc_deg(det2)

  if (ndiff == 0) then
    if (allocated(diff)) deallocate(diff)
    allocate(diff(3,ndiff))
    return
  else if (ndiff > 2) then
    if (allocated(diff)) deallocate(diff)
    allocate(diff(3,0))
    return
  else
    if (allocated(diff)) deallocate(diff)
    allocate(diff(3,ndiff))
  end if
  
  i = 0
  a = 0
  do spin = 1, 2
    diff_ = ieor(det1%str(spin), det2%str(spin))
    do k = 1, strsize
      if (btest(diff_, k-1)) then
        if (btest(det1%str(spin), k-1)) then
          i = i + 1
          diff(1, i) = k
          diff(3, i) = spin
        else
          a = a + 1
          diff(2, a) = k
        end if
      end if
    end do
  end do
end function det_diff


!******************************************************************************** 
!
! Estimate phase factor between two Slater determinants 
!
!    phase = (-1)^{nperm}
!
! where nperm is number of anticomm. relations needed to bring det1 into det2
!
! if diff array (det_diff) is passed as input: it will be used
! if not: the diff array will be created on the fly
!
! diff array has following form:
!
!    i  j  ...
!    a  b  ...
!    s1 s2 ...
!
!******************************************************************************** 
function det_phase(det1, det2, diff) result(phase)
  class(slater_det), intent(in)              :: det1
  type(slater_det), intent(in)               :: det2
  integer, allocatable, optional, intent(in) :: diff(:,:)
  integer                                    :: phase
  !local variables
  integer                                    :: p, exc, ndiff, nperm, spin
  integer, allocatable                       :: diff_(:,:)

  if (present(diff)) then
    allocate(diff_, source=diff)
  else 
    diff_ = det1%diff(det2)
  end if

  nperm = 0
  ndiff = size(diff_, 2)

  do exc = 1, ndiff
    spin = diff_(3,exc)
    if (diff_(2,exc) > diff_(1,exc)) then
      do p =  diff_(1,exc)+1, diff_(2,exc)-1
        if (btest(det1%str(spin), p-1)) nperm = nperm + 1
      end do
    else 
      do p =  diff_(2,exc)+1, diff_(1,exc)-1
        if (btest(det1%str(spin), p-1)) nperm = nperm + 1
      end do
    end if
  end do

  !add permutation for double excitation with crossing (a>j or i>b)
  if (ndiff == 2) then
    if(diff_(3,1) == diff_(3,2)) then
      if (diff_(2,1)>diff_(1,2) .or. diff_(1,1)>diff_(2,2)) nperm = nperm + 1
    end if
  end if

  phase = (-1)**(nperm)
end function det_phase


!******************************************************************************** 
!
! Calculate overlpa between two Slater determinants 
!
!    <Phi_i|Phi_j> = delta_{ij}
!
!******************************************************************************** 
function det_ovlp(det1, det2) result(ovlp)
  class(slater_det), intent(in) :: det1
  type(slater_det), intent(in)  :: det2
  real(wp)                      :: ovlp

  ovlp = 0.0_wp
  if (det1%str(1)==det2%str(1) .and. det1%str(2)==det2%str(2)) ovlp = 1.0_wp
end function det_ovlp


!******************************************************************************** 
!
! Read determinant from character string of the format
!
!    [ i j k l m n ]  [ ii jj kk ll ]  ci_coeff
!      alpha_string     beta_string    corresponding coefficient
!
!******************************************************************************** 
function read_det_from_char(det, char_det) result(coeff)
  class(slater_det), intent(out) :: det
  character(len=*), intent(in)   :: char_det
  complex(wp)                    :: coeff
  !local variables 
  integer                        :: i, posl, posr, ierror
  integer, allocatable           :: left_occ(:), right_occ(:)
  real(wp)                       :: coeff_r, coeff_i
  character(len=:), allocatable  :: left_str, right_str, coeff_str

  posl = scan(char_det, "[")
  posr = scan(char_det, "]")
  left_str = char_det(posl+1:posr-1)

  right_str =  char_det(posr+1:)
  posl = scan(right_str, "[")
  posr = scan(right_str, "]")
  coeff_str = adjustl(right_str(posr+1:))
  right_str = right_str(posl+1:posr-1)

  left_occ = int_arr_from_string(left_str, del=" ")
  det%str(1) = 0
  do i = 1, size(left_occ)
    det%str(1) = ibset(det%str(1), left_occ(i))
  end do
  !det%str(1) = sum(2**left_occ)

  right_occ = int_arr_from_string(right_str, del=" ")
  det%str(2) = 0
  do i = 1, size(right_occ)
    det%str(2) = ibset(det%str(2), right_occ(i))
  end do
  !det%str(2) = sum(2**right_occ)

  coeff_r = 0.0_wp
  coeff_i = 0.0_wp

  read(coeff_str, *, iostat=ierror) coeff_r, coeff_i
  if (ierror /= 0) then
    read(coeff_str, *, iostat=ierror) coeff_r
    coeff_i = 0.0_wp
  end if

  coeff = cmplx(coeff_r, coeff_i, kind=dp)
end function read_det_from_char

end module slater