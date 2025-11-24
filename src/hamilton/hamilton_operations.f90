! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module hamilton_operations

#include "../preproc.inc"
use constants
use mpi
use lapack
use spin_utils
use standalone, only: trace2

implicit none

public

interface sum_gres_potential
  module procedure sum_gres_potential_r, sum_gres_potential_c
end interface sum_gres_potential

contains

!********************************************************************************
!
! Calculate Hartree potential from ERIs in the canonical basis
!
!    J_{pq}^{sigma} = \sum_{tau} \sum_{i} [pq|ii]^{sigma,tau}
!
! Note:
!    Summation over index i runs from first to last
!    potential calculated for states p,q from first_state to last_state
!
!********************************************************************************
subroutine calculate_hartree_potential(h2, first, last, ispin, rspin, hpot, first_state, last_state)
  real(wp), target, intent(in)       :: h2(:,:,:,:,:)
  integer, intent(in)                :: first(2)
  integer, intent(in)                :: last(2)
  integer, intent(in)                :: ispin
  integer, intent(in)                :: rspin
  real(wp), allocatable, intent(out) :: hpot(:,:,:)
  integer, optional, intent(in)      :: first_state
  integer, optional, intent(in)      :: last_state
  !local variables
  integer                            :: ntot, n, i
  integer                            :: spin, spin1, spin2, sp, ispin2
  integer                            :: first_state_, last_state_
  real(wp), pointer                  :: mat(:,:)

  ntot = size(h2, 1)
  ispin2 = size(h2, 5)

  if (present(first_state)) then
    first_state_ = first_state
  else
    first_state_ = 1
  end if

  if (present(last_state)) then
    last_state_ = last_state
  else
    last_state_ = ntot
  end if

  n = last_state_ - first_state_ + 1

  allocate(hpot(n,n,ispin))
  hpot = 0.0_wp

  do spin1 = 1, ispin
    do spin2 = 1, ispin
      sp = get_spin2_channel(spin1, spin2, ispin2)

      do i = first(spin2), last(spin2)
        if (sp==3 .and. spin1==2) then
          mat => h2(i, i, first_state_:last_state_, first_state_:last_state_, sp)
        else 
          mat => h2(first_state_:last_state_, first_state_:last_state_, i, i, sp)
        end if 

        hpot(:,:,spin) = hpot(:,:,spin) + rspin * mat
      end do
    end do
  end do
end subroutine calculate_hartree_potential


!********************************************************************************
!
! Calculate Hartree potential from Cholesky vectors in the canonical basis
!
!    J_{pq}^{sigma} = \sum_g l_g L_{g,pq}^{sigma}      with
!
!    l_g = \sum_{tau} \sum_i L_{g,ii}^{tau}
!
! Note:
!    Summation over index i runs from first to last
!    potential calculated for states p,q from first_state to last_state
!
!********************************************************************************
subroutine calculate_hartree_potential_gres(h2_gres, first, last, ispin, rspin, hpot, first_state, last_state)
  real(wp), target, intent(in)       :: h2_gres(:,:,:,:)
  integer, intent(in)                :: first(2)
  integer, intent(in)                :: last(2)
  integer, intent(in)                :: ispin
  integer, intent(in)                :: rspin
  real(wp), allocatable, intent(out) :: hpot(:,:,:)
  integer, optional, intent(in)      :: first_state
  integer, optional, intent(in)      :: last_state
  !local variables
  integer                            :: ntot, n, g, ng, i
  integer                            :: spin, sp, ispin1
  integer                            :: first_state_, last_state_
  real(wp), allocatable              :: lgmean(:)
  real(wp), pointer                  :: mat(:,:)

  ntot = size(h2_gres, 1)
  ng = size(h2_gres, 3)
  ispin1 = size(h2_gres, 4)

  if (present(first_state)) then
    first_state_ = first_state
  else
    first_state_ = 1
  end if

  if (present(last_state)) then
    last_state_ = last_state
  else
    last_state_ = ntot
  end if

  n = last_state_ - first_state_ + 1

  allocate(lgmean(ng))
  lgmean = 0.0_wp

  do spin = 1, ispin
    sp = get_spin1_channel(spin, ispin1)
    do g = 1, ng
      do i = first(spin), last(spin)
        lgmean(g) = lgmean(g) + rspin * h2_gres(i,i,g,sp)
      end do
    end do
  end do

  allocate(hpot(n,n,ispin))
  hpot = 0.0_wp

  do spin = 1, ispin
    sp = get_spin1_channel(spin, ispin1)
    do g = 1, ng
      mat => h2_gres(first_state_:last_state_, first_state_:last_state_, g, sp)
      hpot(:,:,spin) = hpot(:,:,spin) + lgmean(g) * mat
    end do
  end do
end subroutine calculate_hartree_potential_gres


!********************************************************************************
!
! Calculate g-resolved Hartree potential from Cholesky vectors 
! in the canonical basis
!
!    J_{g,pq}^{sigma} = l_g L_{g,pq}^{sigma}      with
!
!    l_g = \sum_{tau} \sum_i L_{g,ii}^{tau}
!
! Note:
!    Summation over index i runs from first to last
!    potential calculated for states p,q from first_state to last_state
!
!********************************************************************************
subroutine calculate_gres_hartree_potential(h2_gres, first, last, ispin, rspin, hpot_gres, first_state, last_state)
  real(wp), target, intent(in)       :: h2_gres(:,:,:,:)
  integer, intent(in)                :: first(2)
  integer, intent(in)                :: last(2)
  integer, intent(in)                :: ispin
  integer, intent(in)                :: rspin
  real(wp), allocatable, intent(out) :: hpot_gres(:,:,:,:)
  integer, optional, intent(in)      :: first_state
  integer, optional, intent(in)      :: last_state
  !local variables
  integer                            :: ntot, n, g, ng, i
  integer                            :: spin, sp, ispin1
  integer                            :: first_state_, last_state_
  real(wp), allocatable              :: lgmean(:)
  real(wp), pointer                  :: mat(:,:)

  ntot = size(h2_gres, 1)
  ng = size(h2_gres, 3)
  ispin1 = size(h2_gres, 4)

  if (present(first_state)) then
    first_state_ = first_state
  else
    first_state_ = 1
  end if

  if (present(last_state)) then
    last_state_ = last_state
  else
    last_state_ = ntot
  end if

  n = last_state_ - first_state_ + 1

  allocate(lgmean(ng))
  lgmean = 0.0_wp

  do spin = 1, ispin
    sp = get_spin1_channel(spin, ispin1)
    do g = 1, ng
      do i = first(spin), last(spin)
        lgmean(g) = lgmean(g) + rspin * h2_gres(i,i,g,sp)
      end do
    end do
  end do

  allocate(hpot_gres(n,n,ng,ispin))
  hpot_gres = 0.0_wp

  do spin = 1, ispin
    sp = get_spin1_channel(spin, ispin1)

    do g = 1, ng
      mat => h2_gres(first_state_:last_state_, first_state_:last_state_, g, sp)
      hpot_gres(:,:,g,spin) = lgmean(g) * mat   
    end do
  end do
end subroutine calculate_gres_hartree_potential


!********************************************************************************
!
! Calculate exchange potential from ERIs in the canonical basis
!
!    K_{pq}^{sigma} = \sum_{i} [pi|iq]^{sigma,sigma}
!
! Note:
!    Summation over index i runs from first to last
!    potential calculated for states p,q from first_state to last_state
!
!********************************************************************************
subroutine calculate_exchange_potential(h2, first, last, ispin, xpot, first_state, last_state)
  real(wp), target, intent(in)       :: h2(:,:,:,:,:)
  integer, intent(in)                :: first(2)
  integer, intent(in)                :: last(2)
  integer, intent(in)                :: ispin
  real(wp), allocatable, intent(out) :: xpot(:,:,:)
  integer, optional, intent(in)      :: first_state
  integer, optional, intent(in)      :: last_state
  !local variables
  integer                            :: ntot, n, i
  integer                            :: spin, sp, ispin2
  integer                            :: first_state_, last_state_
  real(wp), pointer                  :: mat(:,:)

  ntot = size(h2, 1)
  ispin2 = size(h2, 5)

  if (present(first_state)) then
    first_state_ = first_state
  else
    first_state_ = 1
  end if

  if (present(last_state)) then
    last_state_ = last_state
  else
    last_state_ = ntot
  end if

  n = last_state_ - first_state_ + 1

  allocate(xpot(n,n,ispin))
  xpot = 0.0_wp

  do spin = 1, ispin
    sp = get_spin1_channel(spin, ispin2)
    do i = first(spin), last(spin)
      !debug: may be good to store to the array not to the pointer
      mat => h2(first_state_:last_state_, i, i, first_state_:last_state_, sp)
      xpot(:,:,spin) = xpot(:,:,spin) + mat
    end do
  end do
end subroutine calculate_exchange_potential


!********************************************************************************
!
! Calculate exchange potential from Cholesky vectors in the canonical basis
!
!    K{pq}^{sigma} = \sum_g \sum_{i} L_{g,pi}^{sigma} L_{g,iq}^{sigma}
!
! Note:
!    Summation over index i runs from first to last
!    potential calculated for states p,q from first_state to last_state
!
!********************************************************************************
subroutine calculate_exchange_potential_gres(h2_gres, first, last, ispin, xpot, first_state, last_state)
  real(wp), target, intent(in)       :: h2_gres(:,:,:,:)
  integer, intent(in)                :: first(2)
  integer, intent(in)                :: last(2)
  integer, intent(in)                :: ispin
  real(wp), allocatable, intent(out) :: xpot(:,:,:)
  integer, optional, intent(in)      :: first_state
  integer, optional, intent(in)      :: last_state
  !local variables
  integer                            :: ntot, n, m, g, ng
  integer                            :: spin, sp, ispin1
  integer                            :: first_state_, last_state_
  real(wp), pointer                  :: mat(:,:)

  ntot = size(h2_gres, 1)
  ng = size(h2_gres, 3)
  ispin1 = size(h2_gres, 4)

  if (present(first_state)) then
    first_state_ = first_state
  else
    first_state_ = 1
  end if

  if (present(last_state)) then
    last_state_ = last_state
  else
    last_state_ = ntot
  end if

  n = last_state_ - first_state_ + 1

  allocate(xpot(n,n,ispin))
  xpot = 0.0_wp

  do spin = 1, ispin
    sp = get_spin1_channel(spin, ispin1)
    m = last(spin) - first(spin) + 1
    if (m <= 0) cycle

    do g = 1, ng
      mat => h2_gres(first_state_:last_state_, first(spin):last(spin), g, sp)
      call gemm("n", "t", n, n, m, 1.0_wp, mat, n, mat, n, 1.0_wp, xpot(:,:,spin), n)
    end do
  end do
end subroutine calculate_exchange_potential_gres


!********************************************************************************
!
! Calculate g-resolved exchange potential from Cholesky vectors
! in the canonical basis
!
!    K_{g,pq}^{sigma} = \sum_{i} L_{g,pi}^{sigma} L_{g,iq}^{sigma}
!
! Note:
!    Summation over index i runs from first to last
!    potential calculated for states p,q from first_state to last_state
!
!********************************************************************************
subroutine calculate_gres_exchange_potential(h2_gres, first, last, ispin, xpot_gres, first_state, last_state)
  real(wp), target, intent(in)       :: h2_gres(:,:,:,:)
  integer, intent(in)                :: first(2)
  integer, intent(in)                :: last(2)
  integer, intent(in)                :: ispin
  real(wp), allocatable, intent(out) :: xpot_gres(:,:,:,:)
  integer, optional, intent(in)      :: first_state
  integer, optional, intent(in)      :: last_state
  !local variables
  integer                            :: ntot, g, ng, n, m
  integer                            :: spin, sp, ispin1
  integer                            :: first_state_, last_state_
  real(wp), pointer                  :: mat(:,:)

  ntot = size(h2_gres, 1)
  ng = size(h2_gres, 3)
  ispin1 = size(h2_gres, 4)

  if (present(first_state)) then
    first_state_ = first_state
  else
    first_state_ = 1
  end if

  if (present(last_state)) then
    last_state_ = last_state
  else
    last_state_ = ntot
  end if

  n = last_state_ - first_state_ + 1

  allocate(xpot_gres(n,n,ng,ispin))
  xpot_gres = 0.0_wp

  do spin = 1, ispin
    sp = get_spin1_channel(spin, ispin1)
    m = last(spin) - first(spin) + 1
    if (m <= 0) cycle

    do g = 1, ng
      mat => h2_gres(first_state_:last_state_, first(spin):last(spin), g, sp)
      call gemm("n", "t", n, n, m, 1.0_wp, mat, n, mat, n, 0.0_wp, xpot_gres(:,:,g,spin), n)
    end do
  end do
end subroutine calculate_gres_exchange_potential


!********************************************************************************
!
! Calculate Hartree potential from ERIs using RDM
!
!    J_{pq}^{sigma} = \sum_{tau} \sum_{rs} [pq|rs]^{sigma,tau} R_{rs}^{sigma}
!
! Note:
!    potential calculated for states p,q from first_state to last_state
!
!********************************************************************************
subroutine calculate_hartree_potential_rdm(h2, rdm, ispin, rspin, hpot, first_state, last_state)
  real(wp), target, intent(in)       :: h2(:,:,:,:,:)
  real(wp), intent(in)               :: rdm(:,:,:)
  integer, intent(in)                :: ispin
  integer, intent(in)                :: rspin
  real(wp), allocatable, intent(out) :: hpot(:,:,:)
  integer, optional, intent(in)      :: first_state
  integer, optional, intent(in)      :: last_state
  !local variables
  integer                            :: ntot, n, p, q, r, s, pp, qq
  integer                            :: spin1, spin2, sp, ispin2
  integer                            :: first_state_, last_state_
  real(wp)                           :: val

  ntot = size(h2, 1)
  ispin2 = size(h2, 5)

  if (present(first_state)) then
    first_state_ = first_state
  else
    first_state_ = 1
  end if

  if (present(last_state)) then
    last_state_ = last_state
  else
    last_state_ = ntot
  end if

  n = last_state_ - first_state_ + 1

  allocate(hpot(n,n,ispin))
  hpot = 0.0_wp

  do spin1 = 1, ispin
    do q = first_state_, last_state_
      qq = q - first_state_ + 1
      do p = first_state_, last_state_
        pp = p - first_state_ + 1
        val = 0.0_wp 
        do spin2 = 1, ispin
          sp = get_spin2_channel(spin1, spin2, ispin2)
          do s = 1, ntot
            do r = 1, ntot
              if (sp==3 .and. spin1==2) then
                val = val + rdm(r,s,spin2) * h2(r,s,p,q,sp) * rspin
              else 
                val = val + rdm(r,s,spin2) * h2(p,q,r,s,sp) * rspin
              end if 
            end do
          end do
        end do
        hpot(pp,qq,spin1) = hpot(pp,qq,spin1) + val
      end do
    end do
  end do
end subroutine calculate_hartree_potential_rdm


!********************************************************************************
!
! Calculate Hartree potential from Cholesky vectors using RDM
!
!    J_{pq}^{sigma} = \sum_g l_g L_{g,pq}^{sigma}      with
!
!    l_g = \sum_{tau} \sum_{rs} L_{g,rs}^{tau} R_{rs}^{tau}
!
! Note:
!    potential calculated for states p,q from first_state to last_state
!
!********************************************************************************
subroutine calculate_hartree_potential_rdm_gres(h2_gres, rdm, ispin, rspin, hpot, first_state, last_state)
  real(wp), target, intent(in)       :: h2_gres(:,:,:,:)
  real(wp), intent(in)               :: rdm(:,:,:)
  integer, intent(in)                :: ispin
  integer, intent(in)                :: rspin
  real(wp), allocatable, intent(out) :: hpot(:,:,:)
  integer, optional, intent(in)      :: first_state
  integer, optional, intent(in)      :: last_state
  !local variables
  integer                            :: ntot, n, g, ng
  integer                            :: spin, sp, ispin1
  integer                            :: first_state_, last_state_
  real(wp), allocatable              :: lgmean(:)
  real(wp), pointer                  :: mat(:,:)

  ntot = size(h2_gres, 1)
  ng = size(h2_gres, 3)
  ispin1 = size(h2_gres, 4)

  if (present(first_state)) then
    first_state_ = first_state
  else
    first_state_ = 1
  end if

  if (present(last_state)) then
    last_state_ = last_state
  else
    last_state_ = ntot
  end if

  n = last_state_ - first_state_ + 1

  allocate(lgmean(ng))
  lgmean = 0.0_wp

  do spin = 1, ispin
    sp = get_spin1_channel(spin, ispin1)
    do g = 1, ng
      lgmean(g) = lgmean(g) + rspin * trace2(h2_gres(:,:,g,sp), rdm(:,:,spin))
    end do
  end do

  allocate(hpot(n,n,ispin))
  hpot = 0.0_wp

  do spin = 1, ispin
    sp = get_spin1_channel(spin, ispin1)
    do g = 1, ng
      mat => h2_gres(first_state_:last_state_, first_state_:last_state_, g, sp)
      hpot(:,:,spin) = hpot(:,:,spin) + lgmean(g) * mat
    end do
  end do
end subroutine calculate_hartree_potential_rdm_gres


!********************************************************************************
!
! Calculate g-resolved Hartree potential from Cholesky vectors uding RDM
!
!    J_{g,pq}^{sigma} = l_g L_{g,pq}^{sigma}      with
!
!    l_g = \sum_{tau} \sum_{rs} \sum_i L_{g,rs}^{tau} R_{rs}^{tau}
!
! Note:
!    potential calculated for states p,q from first_state to last_state
!
!********************************************************************************
subroutine calculate_gres_hartree_potential_rdm(h2_gres, rdm, ispin, rspin, hpot_gres, first_state, last_state)
  real(wp), target, intent(in)       :: h2_gres(:,:,:,:)
  real(wp), intent(in)               :: rdm(:,:,:)
  integer, intent(in)                :: ispin
  integer, intent(in)                :: rspin
  real(wp), allocatable, intent(out) :: hpot_gres(:,:,:,:)
  integer, optional, intent(in)      :: first_state
  integer, optional, intent(in)      :: last_state
  !local variables
  integer                            :: ntot, n, g, ng
  integer                            :: spin, sp, ispin1
  integer                            :: first_state_, last_state_
  real(wp), allocatable              :: lgmean(:)
  real(wp), pointer                  :: mat(:,:)

  ntot = size(h2_gres, 1)
  ng = size(h2_gres, 3)
  ispin1 = size(h2_gres, 4)

  if (present(first_state)) then
    first_state_ = first_state
  else
    first_state_ = 1
  end if

  if (present(last_state)) then
    last_state_ = last_state
  else
    last_state_ = ntot
  end if

  n = last_state_ - first_state_ + 1

  allocate(lgmean(ng))
  lgmean = 0.0_wp

  do spin = 1, ispin
    sp = get_spin1_channel(spin, ispin1)
    do g = 1, ng
      lgmean(g) = lgmean(g) + rspin * trace2(h2_gres(:,:,g,sp), rdm(:,:,spin))
    end do
  end do

  allocate(hpot_gres(n,n,ng,ispin))
  hpot_gres = 0.0_wp

  do spin = 1, ispin
    sp = get_spin1_channel(spin, ispin1)

    do g = 1, ng
      mat => h2_gres(first_state_:last_state_, first_state_:last_state_, g, sp)
      hpot_gres(:,:,g,spin) = lgmean(g) * mat   
    end do
  end do
end subroutine calculate_gres_hartree_potential_rdm


!********************************************************************************
!
! Calculate exchange potential from ERIs using RDM
!
!    K_{pq}^{sigma} = \sum_{rs} [pr|sq]^{sigma,sigma} R_{rs}^{sigma}
!
! Note:
!    potential calculated for states p,q from first_state to last_state
!
!********************************************************************************
subroutine calculate_exchange_potential_rdm(h2, rdm, ispin, xpot, first_state, last_state)
  real(wp), target, intent(in)       :: h2(:,:,:,:,:)
  real(wp), intent(in)               :: rdm(:,:,:)
  integer, intent(in)                :: ispin
  real(wp), allocatable, intent(out) :: xpot(:,:,:)
  integer, optional, intent(in)      :: first_state
  integer, optional, intent(in)      :: last_state
  !local variables
  integer                            :: ntot, n, p, q, r, s, pp, qq
  integer                            :: spin, sp, ispin2
  integer                            :: first_state_, last_state_
  real(wp)                           :: val
  real(wp), pointer                  :: mat(:,:)

  ntot = size(h2, 1)
  ispin2 = size(h2, 5)

  if (present(first_state)) then
    first_state_ = first_state
  else
    first_state_ = 1
  end if

  if (present(last_state)) then
    last_state_ = last_state
  else
    last_state_ = ntot
  end if

  n = last_state_ - first_state_ + 1

  allocate(xpot(n,n,ispin))
  xpot = 0.0_wp

  do spin = 1, ispin
    sp = get_spin1_channel(spin, ispin2)
    do q = first_state_, last_state_
      qq = q - first_state_ + 1
      do p = 1, ntot
        pp = p - first_state_ + 1
        val = 0.0_wp 
        do s = 1, ntot
          do r = 1, ntot
            val = val + rdm(r,s,spin) * h2(p,r,s,q,sp)
          end do
        end do
        xpot(pp,qq,spin) = val
      end do
    end do
  end do
end subroutine calculate_exchange_potential_rdm


!********************************************************************************
!
! Calculate exchange potential from Cholesky vectors in the orthonormal basis
!
!    K_{pq}^{sigma} = \sum_{g,rs} L_{g,pr}^{sigma} R_{rs}^{sigma} L_{g,sq}^{sigma}
!
! Note:
!    potential calculated for states p,q from first_state to last_state
!
!********************************************************************************
subroutine calculate_exchange_potential_rdm_gres(h2_gres, rdm, ispin, xpot, first_state, last_state)
  real(wp), target, intent(in)       :: h2_gres(:,:,:,:)
  real(wp), intent(in)               :: rdm(:,:,:)
  integer, intent(in)                :: ispin
  real(wp), allocatable, intent(out) :: xpot(:,:,:)
  integer, optional, intent(in)      :: first_state
  integer, optional, intent(in)      :: last_state
  !local variables
  integer                            :: ntot, n, g, ng
  integer                            :: spin, sp, ispin1
  integer                            :: first_state_, last_state_
  real(wp), allocatable              :: temp(:,:)
  real(wp), pointer                  :: mat(:,:)

  ntot = size(h2_gres, 1)
  ng = size(h2_gres, 3)
  ispin1 = size(h2_gres, 4)

  if (present(first_state)) then
    first_state_ = first_state
  else
    first_state_ = 1
  end if

  if (present(last_state)) then
    last_state_ = last_state
  else
    last_state_ = ntot
  end if

  n = last_state_ - first_state_ + 1

  allocate(xpot(n,n,ispin))
  xpot = 0.0_wp

  allocate(temp(ntot,n))

  do spin = 1, ispin
    sp = get_spin1_channel(spin, ispin1)

    do g = 1, ng
      mat => h2_gres(first_state_:last_state_, :, g, sp)
      call gemm("n", "t", ntot, n, ntot, 1.0_wp, rdm(:,:,spin), ntot, mat, n, 0.0_wp, temp, ntot)
      call gemm("n", "n", n, n, ntot, 1.0_wp, mat, n, temp, ntot, 1.0_wp, xpot(:,:,spin), n)
    end do
  end do
end subroutine calculate_exchange_potential_rdm_gres


!********************************************************************************
!
! Calculate g-resolved exchange potential from Cholesky vectors
! in the orthonormal basis
!
!    K_{g,pq}^{sigma} = \sum_{rs} L_{g,pr}^{sigma} R_{rs,sigma} L_{g,sq}^{sigma}
!
! Note:
!    potential calculated for states p,q from first_state to last_state
!
!********************************************************************************
subroutine calculate_gres_exchange_potential_rdm(h2_gres, rdm, ispin, xpot_gres, first_state, last_state)
  real(wp), target, intent(in)       :: h2_gres(:,:,:,:)
  real(wp), intent(in)               :: rdm(:,:,:)
  integer, intent(in)                :: ispin
  real(wp), allocatable, intent(out) :: xpot_gres(:,:,:,:)
  integer, optional, intent(in)      :: first_state
  integer, optional, intent(in)      :: last_state
  !local variables
  integer                            :: ntot, n, g, ng
  integer                            :: spin, sp, ispin1
  integer                            :: first_state_, last_state_
  real(wp), allocatable              :: temp(:,:)
  real(wp), pointer                  :: mat(:,:)

  ntot = size(h2_gres, 1)
  ng = size(h2_gres, 3)
  ispin1 = size(h2_gres, 4)

  if (present(first_state)) then
    first_state_ = first_state
  else
    first_state_ = 1
  end if

  if (present(last_state)) then
    last_state_ = last_state
  else
    last_state_ = ntot
  end if

  n = last_state_ - first_state_ + 1

  allocate(xpot_gres(n,n,ng,ispin))
  xpot_gres = 0.0_wp

  allocate(temp(ntot,n))

  do spin = 1, ispin
    sp = get_spin1_channel(spin, ispin1)

    do g = 1, ng
      mat => h2_gres(first_state_:last_state_, :, g, sp)
      call gemm("n", "t", ntot, n, ntot, 1.0_wp, rdm(:,:,spin), ntot, mat, n, 0.0_wp, temp, ntot)
      call gemm("n", "n", n, n, ntot, 1.0_wp, mat, n, temp, ntot, 1.0_wp, xpot_gres(:,:,g,spin), n)
    end do
  end do
end subroutine calculate_gres_exchange_potential_rdm


!********************************************************************************
!
! Sum up real-valued g-resolved potential
!
!********************************************************************************
subroutine sum_gres_potential_r(pot_gres, pot)
  real(wp), target, contiguous, intent(in)   :: pot_gres(:,:,:,:)
  real(wp), allocatable, target, intent(out) :: pot(:,:,:)
  !local variables
  integer                                    :: n, m, nm, ng, spin, ispin
  real(wp), allocatable                      :: vec(:)
  real(wp), pointer                          :: pot_gres_ptr(:,:), pot_ptr(:)

  n = size(pot_gres, 1)
  m = size(pot_gres, 2)
  ng = size(pot_gres, 3)
  ispin = size(pot_gres, 4)
  nm = n * m

  allocate(pot(n,m,ispin))

  allocate(vec(ng))
  vec = oner

  do spin = 1, ispin
    pot_gres_ptr(1:nm, 1:ng) => pot_gres(:,:,:,spin)
    pot_ptr(1:nm) => pot(:,:,spin)
    call gemv("n", nm, ng, oner, pot_gres_ptr, nm, vec, 1, zeror, pot_ptr, 1)
  end do
end subroutine sum_gres_potential_r


!********************************************************************************
!
! Sum up complex-valued g-resolved potential
!
!********************************************************************************
subroutine sum_gres_potential_c(pot_gres, pot)
  complex(wp), target, contiguous, intent(in)   :: pot_gres(:,:,:,:)
  complex(wp), allocatable, target, intent(out) :: pot(:,:,:)
  !local variables
  integer                                    :: n, m, nm, ng, spin, ispin
  complex(wp), allocatable                      :: vec(:)
  complex(wp), pointer                          :: pot_gres_ptr(:,:), pot_ptr(:)

  n = size(pot_gres, 1)
  m = size(pot_gres, 2)
  ng = size(pot_gres, 3)
  ispin = size(pot_gres, 4)
  nm = n * m

  allocate(pot(n,m,ispin))

  allocate(vec(ng))
  vec = onec

  do spin = 1, ispin
    pot_gres_ptr(1:nm, 1:ng) => pot_gres(:,:,:,spin)
    pot_ptr(1:nm) => pot(:,:,spin)
    call gemv("n", nm, ng, onec, pot_gres_ptr, nm, vec, 1, zeroc, pot_ptr, 1)
  end do
end subroutine sum_gres_potential_c

end module hamilton_operations