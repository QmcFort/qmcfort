! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module wave_trial

!******************************************************************************** 
!
! This module defines an abstract type WaveTrial that provides an interface 
! for the trial wave functions mainly for the AFQMC calculations
!
! The usage of the WaveTrial objects is split into three categories
!
!    1. Initialization: used to construct the concrete WaveTrial objects
!       It is not a part of the interface since each concrete type is allowed 
!       to have its own set of data for the initialization.
!
!           call init_wave_trial_...(..., wave_trial...)
!
!    2. Setup: calculate arrays that are needed for the efficient energy 
!       evaluation. It is splitted from the initialization since, it is
!       quite expensive
!
!           call wave_trial...%setup(lself)
!
!    3. Energy evaluation: After the initialization and setup, WaveTrial objects
!       can be used for the local energy evaluation
! 
!           call wave_trial%local_energy(...)
!
! Note: In contrast to local energy evaluation, overlap to other determinants
!       can be calculated directly after the initialization
!
!******************************************************************************** 

use constants
use file_handle, only: FileHandle
use qmcfort_io
use mpi, only: comm_world
use profiling
use energy_types
use lapack
use hamilton_vars, only: hamil_descriptor
use hamilton_layout
use hamilton_type
use nosd
use wave_trial_des, only: WaveTrialDes

implicit none

private
public :: WaveTrial

type, abstract :: WaveTrial
  type(WaveTrialDes)      :: des
  type(Hamilton), pointer :: ham
contains
  procedure                                  :: base_assign => base_assign_wave_trial
  procedure                                  :: base_init => base_init_wave_trial
  procedure                                  :: pack_energy => wave_trial_pack_energy
  procedure                                  :: pack_energy_gres => wave_trial_pack_energy_gres

  generic                                    :: assignment(=) => assign_wave_trial
  generic                                    :: trial_energy => trial_energy_c, trial_energy_full
  generic                                    :: initial_energy => initial_energy_c, initial_energy_full
  generic                                    :: get_ovlp => get_ovlp_c, get_ovlp_batch
  generic                                    :: h1_mean => h1_mean_c, h1_mean_batch
  generic                                    :: h2_gres_mean => h2_gres_mean_c, h2_gres_mean_batch
  generic                                    :: exchange_energy => exchange_energy_c, exchange_energy_batch
  generic                                    :: exchange_energy_gres => exchange_energy_c_gres, exchange_energy_batch_gres
  generic                                    :: self_energy => self_energy_c, self_energy_batch
  generic                                    :: self_energy_gres => self_energy_c_gres, self_energy_batch_gres
  generic                                    :: local_energy => local_energy_c_packed, local_energy_batch_packed, &
                                                                local_energy_c, local_energy_batch
  generic                                    :: local_energy_gres => local_energy_c_gres, local_energy_batch_gres

  procedure, private                         :: trial_energy_c, initial_energy_c
  procedure, private                         :: get_ovlp_c, h1_mean_c, h2_gres_mean_c
  procedure, private                         :: exchange_energy_c, exchange_energy_c_gres
  procedure, private                         :: self_energy_c, self_energy_c_gres
  procedure, private                         :: local_energy_c,local_energy_c_packed, local_energy_batch_packed
  procedure, private                         :: local_energy_c_gres

  !deferred routines
  procedure(Iassign_wave_trial), deferred    :: assign_wave_trial
  procedure(Isetup), deferred                :: setup
  procedure(Ireport), deferred               :: report
  procedure(Isizeof), deferred               :: sizeof

  procedure(Itrial_energy), deferred         :: trial_energy_full
  procedure(Itrial_h2_gres_mean), deferred   :: trial_h2_gres_mean
  procedure(Itrial_energy), deferred         :: initial_energy_full
  procedure(Itrial_h2_gres_mean), deferred   :: initial_h2_gres_mean

  procedure(Iget_ovlp), deferred             :: get_ovlp_batch
  procedure(Ih1_mean), deferred              :: h1_mean_batch
  procedure(Ih2_gres_mean), deferred         :: h2_gres_mean_batch
  procedure(Iexchange_energy), deferred      :: exchange_energy_batch
  procedure(Iexchange_energy_gres), deferred :: exchange_energy_batch_gres
  procedure(Iself_energy), deferred          :: self_energy_batch
  procedure(Iself_energy_gres), deferred     :: self_energy_batch_gres
  procedure(Ilocal_energy), deferred         :: local_energy_batch
  procedure(Ilocal_energy_gres), deferred    :: local_energy_batch_gres
end type WaveTrial

abstract interface
  subroutine Iassign_wave_trial(to, from)
    import WaveTrial
    class(WaveTrial), intent(inout) :: to
    class(WaveTrial), intent(in)    :: from
  end subroutine Iassign_wave_trial

  subroutine Isetup(self, lself)
    import WaveTrial
    class(WaveTrial), intent(inout) :: self
    logical, optional, intent(in)   :: lself
  end subroutine Isetup

  subroutine Ireport(self, comm, fh)
    import WaveTrial, mpi_communicator, FileHandle
    class(WaveTrial), intent(inout)    :: self
    type(mpi_communicator), intent(in) :: comm
    type(FileHandle), intent(in)       :: fh
  end subroutine Ireport

  function Isizeof(self) result(mem)
    import WaveTrial, dp
    class(WaveTrial), intent(in) :: self
    real(dp)                     :: mem
  end function Isizeof

  subroutine Itrial_energy(self, e1, lgmean, eh, ex, eself)
    import WaveTrial, wp
    class(WaveTrial), intent(inout)       :: self
    complex(wp), intent(out)              :: e1
    complex(wp), allocatable, intent(out) :: lgmean(:)
    complex(wp), allocatable, intent(out) :: eh(:)
    complex(wp), allocatable, intent(out) :: ex(:)
    complex(wp), allocatable, intent(out) :: eself(:)
  end subroutine Itrial_energy

  subroutine Itrial_h2_gres_mean(self, lgmean)
    import WaveTrial, wp
    class(WaveTrial), intent(inout)       :: self
    complex(wp), allocatable, intent(out) :: lgmean(:)
  end subroutine Itrial_h2_gres_mean

  subroutine Iget_ovlp(self, coeff, ovlp)
    import WaveTrial, wp
    class(WaveTrial), intent(inout)       :: self
    complex(wp), intent(in)               :: coeff(:,:)
    complex(wp), allocatable, intent(out) :: ovlp(:)
  end subroutine Iget_ovlp

  subroutine Ih1_mean(self, coeff, e1)
    import WaveTrial, wp
    class(WaveTrial), intent(inout)       :: self
    complex(wp), intent(in)               :: coeff(:,:)
    complex(wp), allocatable, intent(out) :: e1(:)
  end subroutine Ih1_mean

  subroutine Ih2_gres_mean(self, coeff, lgmean, eh)
    import WaveTrial, wp
    class(WaveTrial), intent(inout)       :: self
    complex(wp), intent(in)               :: coeff(:,:)
    complex(wp), allocatable, intent(out) :: lgmean(:,:)
    complex(wp), allocatable, intent(out) :: eh(:,:)
  end subroutine Ih2_gres_mean

  subroutine Iexchange_energy(self, coeff, ex)
    import WaveTrial, wp
    class(WaveTrial), intent(inout)       :: self
    complex(wp), intent(in)               :: coeff(:,:)
    complex(wp), allocatable, intent(out) :: ex(:)
  end subroutine Iexchange_energy

  subroutine Iexchange_energy_gres(self, coeff, ex)
    import WaveTrial, wp
    class(WaveTrial), intent(inout)       :: self
    complex(wp), intent(in)               :: coeff(:,:)
    complex(wp), allocatable, intent(out) :: ex(:,:)
  end subroutine Iexchange_energy_gres

  subroutine Iself_energy(self, coeff, eself)
    import WaveTrial, wp
    class(WaveTrial), intent(inout)       :: self
    complex(wp), intent(in)               :: coeff(:,:)
    complex(wp), allocatable, intent(out) :: eself(:)
  end subroutine Iself_energy

  subroutine Iself_energy_gres(self, coeff, eself)
    import WaveTrial, wp
    class(WaveTrial), intent(inout)       :: self
    complex(wp), intent(in)               :: coeff(:,:)
    complex(wp), allocatable, intent(out) :: eself(:,:)
  end subroutine Iself_energy_gres

  subroutine Ilocal_energy(self, coeff, e1, lgmean, eh, ex, eself, lsamplex, lself, en)
    import WaveTrial, CEnergy, wp
    class(WaveTrial), intent(inout)                   :: self
    complex(wp), intent(in)                           :: coeff(:,:)
    complex(wp), allocatable, intent(out)             :: e1(:)
    complex(wp), allocatable, intent(out)             :: lgmean(:,:)
    complex(wp), allocatable, intent(out)             :: eh(:)
    complex(wp), allocatable, intent(out)             :: ex(:)
    complex(wp), allocatable, intent(out)             :: eself(:)
    logical, optional, intent(in)                     :: lsamplex
    logical, optional, intent(in)                     :: lself
    type(CEnergy), allocatable, optional, intent(out) :: en(:)
  end subroutine Ilocal_energy

  subroutine Ilocal_energy_gres(self, coeff, e1, lgmean, eh, ex, eself, lsamplex, lself, en)
    import WaveTrial, CEnergy, wp
    class(WaveTrial), intent(inout)                   :: self
    complex(wp), intent(in)                           :: coeff(:,:)
    complex(wp), allocatable, intent(out)             :: e1(:)
    complex(wp), allocatable, intent(out)             :: lgmean(:,:)
    complex(wp), allocatable, intent(out)             :: eh(:,:)
    complex(wp), allocatable, intent(out)             :: ex(:,:)
    complex(wp), allocatable, intent(out)             :: eself(:,:)
    logical, optional, intent(in)                     :: lsamplex
    logical, optional, intent(in)                     :: lself
    type(CEnergy), allocatable, optional, intent(out) :: en(:)
  end subroutine Ilocal_energy_gres
end interface

contains

!******************************************************************************** 
!
! Base Assignment for all WaveTrial extensions
!
!******************************************************************************** 
subroutine base_assign_wave_trial(to, from)
  class(WaveTrial), intent(inout) :: to
  class(WaveTrial), intent(in)    :: from

  to%des = from%des
  to%ham => from%ham
end subroutine base_assign_wave_trial


!******************************************************************************** 
!
! Base setup of the WaveTrial components
!
!******************************************************************************** 
subroutine base_init_wave_trial(self, ham)
  class(WaveTrial), intent(inout)    :: self
  type(Hamilton), target, intent(in) :: ham

  self%ham => ham
end subroutine base_init_wave_trial


!debug: wave_trial : this feature temporarily removed
!******************************************************************************** 
!
! Difference in exchange enrgy between compressed and uncompressed Cholesky vectors
!
!     Delta_x = E_x(Phi_init; tol=0) - E_x(Phi_init; tol)
!
! then:
!
!     E_x(tol=0) \approx E_x(tol) + Delta_x
!
!******************************************************************************** 
!!subroutine estimate_delta_x(self)
!!  class(WaveTrial), intent(inout) :: self
!!  !local
!!  logical                         :: fast_return
!!  complex(wp)                     :: e1, e1_
!!  complex(wp), allocatable        :: lgmean(:), eh(:), ex(:), eself(:)
!!  complex(wp), allocatable        :: lgmean_(:), eh_(:), ex_(:), eself_(:)
!!  
!!
!!  self%des%delta_x = 0.0_wp
!!
!!  if (.not. self%des%compress_x) return
!!  fast_return = .not. (self%des%exchange_mode == "cholesky")
!!  if (fast_return) return
!!
!!  self%des%compress_x = .false.
!!  call self%setup_exchange()
!!  call self%initial_energy(e1, lgmean, eh, ex, eself)
!!
!!  self%des%compress_x = .true.
!!  call self%setup_exchange()
!!  call self%initial_energy(e1_, lgmean_, eh_, ex_, eself_)
!!
!!  self%des%delta_x = real(sum(ex) - sum(ex_), wp)
!!end subroutine estimate_delta_x


!******************************************************************************** 
!
! Pack the result of the local_energy routine to the CEnergy object
!
!******************************************************************************** 
function wave_trial_pack_energy(self, e1, eh, ex) result(en)
  class(WaveTrial), intent(in) :: self
  complex(wp), intent(in)      :: e1
  complex(wp), intent(in)      :: eh
  complex(wp), intent(in)      :: ex
  type(CEnergy)                :: en

  en = CEnergy(self%ham%enuc, e0=self%ham%h0, energy_array=[e1, eh, ex])
end function wave_trial_pack_energy


!******************************************************************************** 
!
! Pack the result of the local_energy_gres routine to the CEnergy object
!
!******************************************************************************** 
function wave_trial_pack_energy_gres(self, e1, eh, ex) result(en)
  class(WaveTrial), intent(in) :: self
  complex(wp), intent(in)      :: e1
  complex(wp), intent(in)      :: eh(:)
  complex(wp), intent(in)      :: ex(:)
  type(CEnergy)                :: en

  en = CEnergy(self%ham%enuc, e0=self%ham%h0, energy_array=[e1, sum(eh), sum(ex)])
end function wave_trial_pack_energy_gres


!******************************************************************************** 
!  
! Pack trial energy directly to the CEnergy object
!
!******************************************************************************** 
subroutine trial_energy_c(self, en)
  class(WaveTrial), intent(inout) :: self
  type(CEnergy), intent(out)      :: en
  !local
  complex(wp)                     :: e1
  complex(wp), allocatable        :: lgmean(:), eh(:), ex(:), eself(:)

  call self%trial_energy(e1, lgmean, eh, ex, eself)
  en = self%pack_energy_gres(e1, eh, ex)
end subroutine trial_energy_c


!******************************************************************************** 
!  
! Pack initial energy directly to the CEnergy object
!
!******************************************************************************** 
subroutine initial_energy_c(self, en)
  class(WaveTrial), intent(inout) :: self
  type(CEnergy), intent(out)      :: en
  !local
  complex(wp)                     :: e1
  complex(wp), allocatable        :: lgmean(:), eh(:), ex(:), eself(:)

  call self%initial_energy(e1, lgmean, eh, ex, eself)
  en = self%pack_energy_gres(e1, eh, ex)
end subroutine initial_energy_c


!******************************************************************************** 
!
! Calculate overlap for a single walker
!
!    overlap = <Psi_T | Phi>
!
! The actual implementation in one of WaveTrial... modules
!
!******************************************************************************** 
subroutine get_ovlp_c(self, coeff, ovlp)
  class(WaveTrial), intent(inout) :: self
  complex(wp), intent(in)         :: coeff(:,:)
  complex(wp), intent(out)        :: ovlp
  !local
  complex(wp), allocatable        :: ovlp_(:)

  call self%get_ovlp(coeff, ovlp_)
  ovlp = ovlp_(1)
end subroutine get_ovlp_c


!******************************************************************************** 
!
! Calculate one-electron energy for a single walker
!
!    E1 = <Psi_T | H1 | Phi>
!
! The actual implementation in one of WaveTrial... modules
!
!******************************************************************************** 
subroutine h1_mean_c(self, coeff, mean)
  class(WaveTrial), intent(inout) :: self
  complex(wp), intent(in)         :: coeff(:,:)
  complex(wp), intent(out)        :: mean
  !local
  complex(wp), allocatable        :: mean_(:)
  
  call self%h1_mean(coeff, mean_)
  mean = mean_(1)
end subroutine h1_mean_c


!******************************************************************************** 
!
! Calculate expectation value of the Cholesky vectors for a single walker
!
!    l_g = <Psi_T | L_g | Phi>
!
! The actual implementation in one of WaveTrial... modules
!
!******************************************************************************** 
subroutine h2_gres_mean_c(self, coeff, lgeman, eh)
  class(WaveTrial), intent(inout)       :: self
  complex(wp), intent(in)               :: coeff(:,:)
  complex(wp), allocatable, intent(out) :: lgeman(:)
  complex(wp), allocatable, intent(out) :: eh(:)
  !local
  complex(wp), allocatable              :: lgeman_(:,:), eh_(:,:)
  
  call self%h2_gres_mean(coeff, lgeman_, eh_)

  allocate(lgeman(size(lgeman_, 1)))
  lgeman = lgeman_(:,1)

  allocate(eh(size(eh_, 1)))
  eh = eh_(:,1)
end subroutine h2_gres_mean_c


!******************************************************************************** 
!
! Calculate exchange energy for a single walker
!
!    Ex = <Psi_T| H_x | Phi> 
!
! The actual implementation in one of WaveTrial... modules
!
!******************************************************************************** 
subroutine exchange_energy_c(self, coeff, ex)
  class(WaveTrial), intent(inout) :: self
  complex(wp), intent(in)         :: coeff(:,:)
  complex(wp), intent(out)        :: ex
  !local
  complex(wp), allocatable        :: ex_(:)
  
  call self%exchange_energy(coeff, ex_)
  ex = ex_(1)
end subroutine exchange_energy_c

!******************************************************************************** 
!
! Calculate g-resolved exchange energy for a single walker
!
!    Ex_g = <Psi_T| H_x | Phi> 
!
! The actual implementation in one of WaveTrial... modules
!
!******************************************************************************** 
subroutine exchange_energy_c_gres(self, coeff, ex)
  class(WaveTrial), intent(inout)       :: self
  complex(wp), intent(in)               :: coeff(:,:)
  complex(wp), allocatable, intent(out) :: ex(:)
  !local
  complex(wp), allocatable              :: ex_(:,:)
  
  call self%exchange_energy_gres(coeff, ex_)
  allocate(ex(size(ex_,1)))
  ex = ex_(:,1)
end subroutine exchange_energy_c_gres


!******************************************************************************** 
!
! Calculate self energy for a single walker
!
!    Eself = <Psi_T| Hself |Phi>
!
! The actual implementation in one of WaveTrial... modules
!
!******************************************************************************** 
subroutine self_energy_c(self, coeff, eself)
  class(WaveTrial), intent(inout) :: self
  complex(wp), intent(in)         :: coeff(:,:)
  complex(wp), intent(out)        :: eself
  !local
  complex(wp), allocatable        :: eself_(:)
  
  call self%self_energy(coeff, eself_)
  eself = eself_(1)
end subroutine self_energy_c


!******************************************************************************** 
!
! Calculate g-resolved self energy for a single walker
!
!    Eself_g = <Psi_T | L_g^2 |Phi>
!
! The actual implementation in one of WaveTrial... modules
!
!******************************************************************************** 
subroutine self_energy_c_gres(self, coeff, eself)
  class(WaveTrial), intent(inout)       :: self
  complex(wp), intent(in)               :: coeff(:,:)
  complex(wp), allocatable, intent(out) :: eself(:)
  !local
  complex(wp), allocatable              :: eself_(:,:)
  
  call self%self_energy_gres(coeff, eself_)
  allocate(eself(size(eself_, 1)))
  eself = eself_(:,1)
end subroutine self_energy_c_gres


!******************************************************************************** 
!
! Calculate local energy contributions for a single walker
!
!    E1 = <Psi_T | H1 | Phi>
!    l_g = <Psi_T | L_g | Phi>
!    Eh = <Psi_T | Hh | Phi>
!    Ex = <Psi_T | Hx | Phi>
!    Eself = <Psi_T | Hself | Phi>
!
! The actual implementation in one of WaveTrial... modules
!
!******************************************************************************** 
subroutine local_energy_c(self, coeff, e1, lgmean, eh, ex, eself, lsamplex, lself)
  class(WaveTrial), intent(inout)       :: self
  complex(wp), intent(in)               :: coeff(:,:)
  complex(wp), intent(out)              :: e1
  complex(wp), allocatable, intent(out) :: lgmean(:)
  complex(wp), intent(out)              :: eh
  complex(wp), intent(out)              :: ex
  complex(wp), intent(out)              :: eself
  logical, optional, intent(in)         :: lsamplex
  logical, optional, intent(in)         :: lself
  !local
  integer                               :: ng, nw
  complex(wp), allocatable              :: lgmean_(:,:), e1_(:), eh_(:), ex_(:), eself_(:)
  
  call self%local_energy(coeff, e1_, lgmean_, eh_, ex_, eself_, lsamplex, lself)

  ng = size(lgmean_, 1)
  nw = size(lgmean_, 2)

  e1 = e1_(1)
  eh = eh_(1)
  ex = ex_(1)
  eself = eself_(1)

  allocate(lgmean(ng))
  lgmean = lgmean_(:,1)
end subroutine local_energy_c


!******************************************************************************** 
!
! Pack local energy contributions of a single walker to the CEnergy variable
!
!******************************************************************************** 
subroutine local_energy_c_packed(self, coeff, en, lsamplex, lself)
  class(WaveTrial), intent(inout)       :: self
  complex(wp), intent(in)               :: coeff(:,:)
  type(CEnergy), intent(out)            :: en
  logical, optional, intent(in)         :: lsamplex
  logical, optional, intent(in)         :: lself
  !local
  complex(wp)                           :: e1, eh, ex, eself
  complex(wp), allocatable              :: lgmean(:)

  call self%local_energy(coeff, e1, lgmean, eh, ex, eself, lsamplex, lself)
  en = self%pack_energy(e1, eh, ex)
end subroutine local_energy_c_packed


!******************************************************************************** 
!
! Pack local energy contributions of a batch of walkers to the CEnergy variables
!
!******************************************************************************** 
subroutine local_energy_batch_packed(self, coeff, en, lsamplex, lself)
  class(WaveTrial), intent(inout)         :: self
  complex(wp), intent(in)                 :: coeff(:,:)
  type(CEnergy), allocatable, intent(out) :: en(:)
  logical, optional, intent(in)           :: lsamplex
  logical, optional, intent(in)           :: lself
  !local
  integer                                 :: w
  complex(wp), allocatable                :: e1(:), eh(:), ex(:), eself(:), lgmean(:,:)

  call self%local_energy(coeff, e1, lgmean, eh, ex, eself, lsamplex, lself)

  allocate(en(size(e1)))
  do w = 1, size(e1)
    en(w) = self%pack_energy(e1(w), eh(w), ex(w))
  end do
end subroutine local_energy_batch_packed


!******************************************************************************** 
!
! Calculate g-resolved local energy contributions for a single walker
!
!    E1 = <Psi_T | H1 | Phi>
!    l_g = <Psi_T | L_g | Phi>
!    Eh = <Psi_T | Hh | Phi>
!    Ex_g = <Psi_T | Hx_g | Phi>
!    Eself_g = <Psi_T | L_g^2 | Phi>
!
! The actual implementation in one of WaveTrial... modules
!
!******************************************************************************** 
subroutine local_energy_c_gres(self, coeff, e1, lgmean, eh, ex, eself, lsamplex, lself)
  class(WaveTrial), intent(inout)       :: self
  complex(wp), intent(in)               :: coeff(:,:)
  complex(wp), intent(out)              :: e1
  complex(wp), allocatable, intent(out) :: lgmean(:)
  complex(wp), allocatable, intent(out) :: eh(:)
  complex(wp), allocatable, intent(out) :: ex(:)
  complex(wp), allocatable, intent(out) :: eself(:)
  logical, optional, intent(in)         :: lsamplex
  logical, optional, intent(in)         :: lself
  !local
  integer                               :: ng, nw
  complex(wp), allocatable              :: e1_(:), lgmean_(:,:), eh_(:,:), ex_(:,:), eself_(:,:)
  
  call self%local_energy_gres(coeff, e1_, lgmean_, eh_, ex_, eself_, lsamplex, lself)

  ng = size(lgmean_, 1)
  nw = size(lgmean_, 2)

  e1 = e1_(1)

  allocate(lgmean(ng))
  lgmean = lgmean_(:,1)

  allocate(eh(ng))
  eh = eh_(:,1)

  allocate(ex(ng))
  ex = ex_(:,1)

  allocate(eself(ng))
  eself = eself_(:,1)
end subroutine local_energy_c_gres

end module wave_trial