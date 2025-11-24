! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module test_energy_types

use energy_types

use constants, only: sp, dp, wp, spprec
use fruit, only: assert_equals, run_test_case, test_set_summary

implicit none

private
public :: test_energy_types_set

contains

subroutine test_energy_types_set
  call run_test_case(test_en_assign_simple, "Test simple assignment of the energy type")
  call run_test_case(test_en_reset, "Test the reset feature of the energy type")
  call run_test_case(test_en_add, "Test the addition of two energy objects")
  call run_test_case(test_en_add_mult, "Test the addition and multiplication of energy objects with scalars")
  call run_test_case(test_en_add_div, "Test the addition and division of energy objects with scalars")

  call test_set_summary("src/energy_types.f90")
end subroutine test_energy_types_set

subroutine test_en_assign_simple
  type(Energy) :: e1, e2
  
  e1 = Energy()
  e1%enuc = 0.0_wp
  e1%e1 = 1.0_wp
  e1%eh = 2.0_wp
  e1%ex = 3.0_wp

  e2 = e1
  
  call assert_equals(e1%e1, e2%e1)
  call assert_equals(e1%eh, e2%eh)
  call assert_equals(e1%ex, e2%ex)
end subroutine test_en_assign_simple

subroutine test_en_reset
  type(CEnergy)       :: en
  complex(wp)         :: earr(en_len)
  real(wp), parameter :: enuc = 10.0_wp, e0=0.0_wp
  
  earr = [1.0_wp, 2.0_wp, 3.0_wp]

  en = CEnergy(enuc=enuc, e0=e0, energy_array=earr)
  call en%reset()

  call assert_equals(enuc, en%enuc)
end subroutine test_en_reset

subroutine test_en_add
  type(CEnergy)       :: en1, en2, en, en_
  real(wp), parameter :: enuc = 1.0_wp, e0 = 0.0_wp
  complex(wp)         :: earr(en_len) = [(1.0_wp,1.0_wp), (2.0_wp,2.0_wp), (3.0_wp,3.0_wp)]
  complex(wp)         :: earr2(en_len) = [(2.0_wp,2.0_wp), (4.0_wp,4.0_wp), (6.0_wp,6.0_wp)]
  
  en1 = CEnergy(enuc=enuc, e0=e0, energy_array=earr)
  en2 = CEnergy(enuc=enuc, e0=e0, energy_array=earr)


  en = CEnergy(enuc=enuc, e0=e0, energy_array=earr2)
  en_ = en1 + en2

  call assert_equals(en%enuc, en_%enuc)
  call assert_equals(en%ex, en%ex)
  call assert_equals(en%total_energy(), en_%total_energy())
end subroutine test_en_add

subroutine test_en_add_mult
  type(CEnergy)       :: en1, en2, en, en_
  real(wp), parameter :: enuc = 1.0_wp, e0=-4.0_wp, factor = 2.0_wp
  complex(wp)         :: earr(en_len) = [(1.0_wp,1.0_wp), (2.0_wp,2.0_wp), (3.0_wp,3.0_wp)]
  complex(wp)         :: earr2(en_len) = [(3.0_wp,3.0_wp), (6.0_wp,6.0_wp), (9.0_wp,9.0_wp)]
  
  en1 = CEnergy(enuc=enuc, e0=e0, energy_array=earr)
  en2 = CEnergy(enuc=enuc, e0=e0, energy_array=earr)

  en = CEnergy(enuc=enuc, e0=e0, energy_array=earr2)
  en_ = en2 + factor * en1
  
  call assert_equals(en%e0, en_%e0)
  call assert_equals(en%e1, en_%e1)
  call assert_equals(en%electron_energy(), en_%electron_energy())
end subroutine test_en_add_mult

subroutine test_en_add_div
  type(CEnergy)       :: en1, en2, en, en_
  real(wp), parameter :: enuc = 1.0_wp, e0=-4.0_wp, factor = 2.0_wp
  complex(wp)         :: earr(en_len) = [(1.0_wp,1.0_wp), (2.0_wp,2.0_wp), (3.0_wp,3.0_wp)]
  
  en1 = CEnergy(enuc=enuc, e0=e0, energy_array=earr)
  en2 = CEnergy(enuc=enuc, e0=e0, energy_array=earr)
  
  en = CEnergy(enuc=enuc, e0=e0, energy_array=earr)
  en_ = (en1 + en2) / factor
  
  call assert_equals(en%enuc, en_%enuc)
  call assert_equals(en%electron_energy(), en_%electron_energy())
  call assert_equals(en%total_energy(), en_%total_energy())
end subroutine test_en_add_div

end module test_energy_types
