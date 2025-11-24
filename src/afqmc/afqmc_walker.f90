! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module afqmc_walker
!******************************************************************************** 
!
! AfqmcWalker type implementation
!
!******************************************************************************** 

#include "../preproc.inc"
use constants
use file_handle, only: FileHandle
use mpi, only: mpi_communicator
use profiling
use standalone, only: phase
use lapack, only: determinant
use energy_types
use hamilton_vars, only: hamil_descriptor
use nosd, only: get_ovlp_matrix
use wave_trial, only: WaveTrial

implicit none

type AfqmcWalker
  integer :: id
  real(wp)                 :: weight
  real(wp)                 :: bare_weight
  real(wp)                 :: phase

  real(wp)                 :: rew_w, rew_p
  real(wp)                 :: dphase
  complex(wp)              :: acc

  type(CEnergy)            :: energy
  complex(wp)              :: overlap, overlap_old, doverlap
  complex(wp)              :: energyl, energyl_old
  complex(wp)              :: energyh, energyh_old

  real(wp), allocatable    :: random_field(:)
  complex(wp), allocatable :: shifted_field(:)
  complex(wp), allocatable :: shift(:)
  complex(wp), allocatable :: scale(:)

  complex(wp), allocatable :: lgmean(:)
  complex(wp), allocatable :: eh(:)
  complex(wp), allocatable :: ex(:)
  complex(wp), allocatable :: eself(:)
  complex(wp), allocatable :: coeff(:,:)

  complex(wp)              :: isf
  complex(wp)              :: isf_mf
  
  logical                  :: move
  integer                  :: life
  integer                  :: revlife
contains
  generic                  :: assignment(=) => assign_afqmc_walker
  procedure                :: alloc => alloc_afqmc_walker
  procedure                :: dealloc => dealloc_afqmc_walker
  procedure                :: hybrid_energy => calculate_hybrid_energy
  procedure                :: get_ovlp => get_walker_overlap
  procedure                :: sizeof => sizeof_afqmc_walker
  procedure                :: read_  => read_afqmc_walker
  procedure                :: write_ => write_afqmc_walker
  procedure, private       :: assign_afqmc_walker
end type AfqmcWalker

interface AfqmcWalker
  module procedure finit_afqmc_walker
end interface AfqmcWalker

interface alloc_walker
  module procedure alloc_afqmc_walkers_1d, alloc_afqmc_walkers_2d
end interface alloc_walker

interface dealloc_walker
  module procedure dealloc_afqmc_walkers_1d, dealloc_afqmc_walkers_2d
end interface dealloc_walker

interface copy_afqmc_walkers
  module procedure copy_afqmc_walkers_1d, copy_afqmc_walkers_2d
end interface

interface batch_walkers
  module procedure batch_walkers_1d, batch_walkers_2d
end interface batch_walkers

interface unbatch_walkers
  module procedure unbatch_walkers_1d, unbatch_walkers_2d
end interface unbatch_walkers

interface unbatch_walkers_gres
  module procedure unbatch_walkers_gres_1d, unbatch_walkers_gres_2d
end interface unbatch_walkers_gres

contains

!******************************************************************************** 
!
! Initialization of the AFQMCWalker object
!
!******************************************************************************** 
subroutine init_afqmc_walker(phi_t, coeff, self, sp_proj)
  class(WaveTrial), target, intent(inout) :: phi_t
  complex(wp), intent(in)                 :: coeff(:,:,:)
  type(AfqmcWalker), intent(out)          :: self
  logical, optional, intent(in)           :: sp_proj
  !local variables
  integer                                 :: spin, w1, w2, nel
  logical                                 :: sp_proj_
  complex(wp)                             :: e1
  complex(wp), allocatable                :: lgmean(:), eh(:), ex(:), eself(:)
  type(hamil_descriptor), pointer         :: hdes

  sp_proj_ = .false.
  if (present(sp_proj)) sp_proj_ = sp_proj

  hdes => phi_t%ham%des

  call self%alloc(hdes)

  do spin = 1, hdes%ispin
    nel = hdes%nel(spin)
    if (nel == 0) cycle
    w1 = hdes%nell(spin) + 1
    w2 = hdes%nell(spin) + nel
    
    if (sp_proj_) then
      self%coeff(:,w1:w2) = coeff(:,1:nel,1)
    else
      self%coeff(:,w1:w2) = coeff(:,1:nel,spin)
    end if
  end do

  call phi_t%get_ovlp(self%coeff, self%overlap)
  self%overlap_old = self%overlap

  self%id = 1
  self%weight = abs(self%overlap)
  self%bare_weight = abs(self%overlap)
  self%phase = phase(self%overlap)
  self%rew_w = 1.0_wp
  self%rew_p = 0.0_wp
  self%dphase = 0.0_wp
  self%acc = onec

  call phi_t%initial_energy(e1, lgmean, eh, ex, eself)
  self%lgmean = lgmean
  self%eh = eh
  self%ex = ex
  self%eself = eself

  self%energy = phi_t%pack_energy_gres(e1, eh, ex)
  self%energyl = self%energy%electron_energy()
  self%energyl_old = self%energy%electron_energy()
  self%energyh = self%energy%electron_energy()
  self%energyh_old = self%energy%electron_energy()

  self%isf = zeroc
  self%isf_mf = zeroc

  self%life = 0
  self%revlife = 0
  self%move = .true.
end subroutine init_afqmc_walker


!******************************************************************************** 
!
! AfqmcWalker constructor
!
!******************************************************************************** 
function finit_afqmc_walker(phi_t, coeff, sp_proj) result(self)
  class(WaveTrial), intent(inout) :: phi_t
  complex(wp), intent(in)         :: coeff(:,:,:)
  logical, optional, intent(in)   :: sp_proj
  type(AfqmcWalker)               :: self

  call init_afqmc_walker(phi_t, coeff, self, sp_proj)
end function finit_afqmc_walker


!******************************************************************************** 
!
! Assignment operator for the AfqmcWalker type
!	
!******************************************************************************** 
subroutine assign_afqmc_walker(to, from)
  class(AfqmcWalker), intent(inout) :: to
  type(AfqmcWalker), intent(in)     :: from

  to%id = from%id
  to%weight = from%weight
  to%bare_weight = from%bare_weight
  to%phase = from%phase

  to%rew_w = from%rew_w
  to%rew_p = from%rew_p
  to%dphase = from%dphase
  to%acc = from%acc

  to%overlap = from%overlap
  to%overlap_old = from%overlap_old
  to%doverlap = from%doverlap

  to%energy = from%energy
  to%energyl = from%energyl
  to%energyl_old = from%energyl_old
  to%energyh = from%energyh
  to%energyh_old = from%energyh_old

  to%isf = from%isf
  to%isf_mf = from%isf_mf

  to%life = from%life
  to%revlife = from%revlife
  to%move = from%move

  if (allocated(from%random_field)) then 
    if (allocated(to%random_field)) deallocate(to%random_field)
    allocate(to%random_field, source=from%random_field)
  end if

  if (allocated(from%shifted_field)) then 
    if (allocated(to%shifted_field)) deallocate(to%shifted_field)
    allocate(to%shifted_field, source=from%shifted_field)
  end if

  if (allocated(from%shift)) then
    if (allocated(to%shift)) deallocate(to%shift)
    allocate(to%shift, source=from%shift)
  end if

  if (allocated(from%scale)) then 
    if (allocated(to%scale)) deallocate(to%scale)
    allocate(to%scale, source=from%scale)
  end if

  if (allocated(from%lgmean)) then
    if (allocated(to%lgmean)) deallocate(to%lgmean)
    allocate(to%lgmean, source=from%lgmean)
  end if

  if (allocated(from%eh)) then 
    if (allocated(to%eh)) deallocate(to%eh)
    allocate(to%eh, source=from%eh)
  end if

  if (allocated(from%ex)) then 
    if (allocated(to%ex)) deallocate(to%ex)
    allocate(to%ex, source=from%ex)
  end if

  if (allocated(from%eself)) then
    if (allocated(to%eself)) deallocate(to%eself)
    allocate(to%eself, source=from%eself)
  end if

  if (allocated(from%coeff)) then 
    if (allocated(to%coeff)) deallocate(to%coeff)
    allocate(to%coeff, source=from%coeff)
  end if
end subroutine assign_afqmc_walker


!******************************************************************************** 
!
! Allocate ararys in AfqmcWalker variable
!
!******************************************************************************** 
subroutine alloc_afqmc_walker(self, hdes)
  class(AfqmcWalker), intent(inout)  :: self
  type(hamil_descriptor), intent(in) :: hdes
  
  call self%dealloc()
  allocate(self%lgmean(hdes%ng)); self%lgmean = zeroc
  allocate(self%eh(hdes%ng)); self%eh = zeroc
  allocate(self%ex(hdes%ng)); self%ex = zeroc
  allocate(self%eself(hdes%ng)); self%eself = zeroc
  allocate(self%coeff(hdes%n, hdes%nocc)); self%coeff = zeroc
end subroutine alloc_afqmc_walker


!******************************************************************************** 
!
! Deallocate ararys in AfqmcWalker variable
!
!******************************************************************************** 
subroutine dealloc_afqmc_walker(self)
  class(AfqmcWalker), intent(inout) :: self
  
  if (allocated(self%random_field)) deallocate(self%random_field)
  if (allocated(self%shifted_field)) deallocate(self%shifted_field)
  if (allocated(self%shift)) deallocate(self%shift)
  if (allocated(self%scale)) deallocate(self%scale)
  if (allocated(self%lgmean)) deallocate(self%lgmean)
  if (allocated(self%eh)) deallocate(self%ex)
  if (allocated(self%ex)) deallocate(self%ex)
  if (allocated(self%eself)) deallocate(self%eself)
  if (allocated(self%coeff)) deallocate(self%coeff)
end subroutine dealloc_afqmc_walker


!******************************************************************************** 
!
! Calculate hybrid energy of the AFQMC walker
!
!    E_h = - log(S'/S + isf_mf + isf) / tau
!
!******************************************************************************** 
function calculate_hybrid_energy(walker, tau) result(hybrid_energy)
  class(AfqmcWalker), intent(in) :: walker
  real(wp), intent(in)           :: tau
  complex(wp)                    :: hybrid_energy

  hybrid_energy = -(log(walker%overlap / walker%overlap_old) + walker%isf_mf + walker%isf) / tau
end function calculate_hybrid_energy


!******************************************************************************** 
!
! Calculate AFqmcWalker overlap <Phi|Phi>
!
!******************************************************************************** 
function get_walker_overlap(self, hdes) result(det)
  class(AfqmcWalker), intent(inout)  :: self
  type(hamil_descriptor), intent(in) :: hdes
  complex(wp)                        :: det
  !local variables
  integer                            :: spin, nocc, w1, w2
  complex(wp)                        :: det_
  complex(wp), allocatable           :: ovlp(:,:)

  det = onec

  do spin = 1, hdes%ispin
    nocc = hdes%nel(spin)
    if (nocc == 0) cycle
    w1 = hdes%nell(spin) + 1
    w2 = hdes%nell(spin) + nocc
    call get_ovlp_matrix(self%coeff(:,w1:w2), self%coeff(:,w1:w2), ovlp)
    call determinant(ovlp, det_)
    det = det * det_**hdes%rspin
  end do
end function get_walker_overlap


!******************************************************************************** 
!
! Memory size of the AfqmcWalker instnce
!
!******************************************************************************** 
function sizeof_afqmc_walker(self) result(mem)
  class(AfqmcWalker), intent(in) :: self
  real(wp)                       :: mem

  mem = sizeof(self) +  sizeof(self%random_field) + sizeof(self%shifted_field) + size(self%shift) + &
        sizeof(self%scale) + sizeof(self%lgmean) + sizeof(self%eh) + sizeof(self%ex) + sizeof(self%eself) + sizeof(self%coeff)
end function sizeof_afqmc_walker


!******************************************************************************** 
!
! Read the AfqmcWalker obejct from the file
!
!******************************************************************************** 
subroutine read_afqmc_walker(self, fh)
  class(AfqmcWalker), intent(inout) :: self
  type(FileHandle), intent(in)      :: fh
  
  read(fh%funit) self%id
  read(fh%funit) self%lgmean
  read(fh%funit) self%eh
  read(fh%funit) self%ex
  read(fh%funit) self%eself
  read(fh%funit) self%coeff

  read(fh%funit) self%overlap, self%overlap_old
  read(fh%funit) self%isf_mf, self%isf
  read(fh%funit) self%weight, self%bare_weight, self%rew_w, self%rew_p, self%phase, self%dphase, self%acc
  read(fh%funit) self%energyl, self%energyl_old
  read(fh%funit) self%energyh, self%energyh_old
  read(fh%funit) self%energy%enuc, self%energy%e0
  read(fh%funit) self%energy%e1, self%energy%eh
  read(fh%funit) self%energy%ex

  read(fh%funit) self%life
  read(fh%funit) self%revlife
  read(fh%funit) self%move
end subroutine read_afqmc_walker


!******************************************************************************** 
!
! Write the AfqmcWalker obejct to the file
!
!******************************************************************************** 
subroutine write_afqmc_walker(self, fh)
  class(AfqmcWalker), intent(in) :: self
  type(FileHandle), intent(in)   :: fh

  write(fh%funit) self%id
  write(fh%funit) self%lgmean
  write(fh%funit) self%eh
  write(fh%funit) self%ex
  write(fh%funit) self%eself
  write(fh%funit) self%coeff
  
  write(fh%funit) self%overlap, self%overlap_old
  write(fh%funit) self%isf_mf, self%isf
  write(fh%funit) self%weight, self%bare_weight, self%rew_w, self%rew_p, self%phase, self%dphase, self%acc
  write(fh%funit) self%energyl, self%energyl_old
  write(fh%funit) self%energyh, self%energyh_old
  write(fh%funit) self%energy%enuc, self%energy%e0
  write(fh%funit) self%energy%e1, self%energy%eh
  write(fh%funit) self%energy%ex

  write(fh%funit) self%life
  write(fh%funit) self%revlife
  write(fh%funit) self%move
end subroutine write_afqmc_walker


!******************************************************************************** 
!
! Allocate arary of the AfqmcWalker objects
!
!******************************************************************************** 
subroutine alloc_afqmc_walkers_1d(walkers, hdes, nw)
  type(AfqmcWalker), allocatable, intent(inout) :: walkers(:)
  type(hamil_descriptor), intent(in)            :: hdes
  integer, intent(in)                           :: nw
  !local variables
  integer                                       :: w
  
  if (allocated(walkers)) call dealloc_walker(walkers)
  allocate(walkers(nw))
  
  do w = 1, nw
    call walkers(w)%alloc(hdes)
  end do
end subroutine alloc_afqmc_walkers_1d

subroutine alloc_afqmc_walkers_2d(walkers, hdes, nw, nt)
  type(AfqmcWalker), allocatable, intent(inout) :: walkers(:,:)
  type(hamil_descriptor), intent(in)            :: hdes
  integer, intent(in)                           :: nw
  integer, intent(in)                           :: nt
  !local variables
  integer                                       :: w, i
  
  if (allocated(walkers)) call dealloc_walker(walkers)
  allocate(walkers(nw,nt))
  
  do i = 1, nt
    do w = 1, nw
      call walkers(w,i)%alloc(hdes)
    end do
  end do
end subroutine alloc_afqmc_walkers_2d


!******************************************************************************** 
!
! Deallocate arary of the AfqmcWalker objects
!
!******************************************************************************** 
subroutine dealloc_afqmc_walkers_1d(walkers)
  type(AfqmcWalker), allocatable, intent(inout) :: walkers(:)
  !local variables
  integer                                       :: w
  
  if (allocated(walkers)) then
    do w = 1, size(walkers)
      call walkers(w)%dealloc()
    end do
    deallocate(walkers)
  end if
end subroutine dealloc_afqmc_walkers_1d

subroutine dealloc_afqmc_walkers_2d(walkers)
  type(AfqmcWalker), allocatable, intent(inout) :: walkers(:,:)
  !local variables
  integer                                       :: w, i
  
  if (allocated(walkers)) then
    do i = 1, size(walkers, 2)
      do w = 1, size(walkers, 1)
        call walkers(w,i)%dealloc()
      end do
    end do
    deallocate(walkers)
  end if
end subroutine dealloc_afqmc_walkers_2d


!******************************************************************************** 
!
! Copy an array of the AfqmcWalker objects
!
!******************************************************************************** 
subroutine copy_afqmc_walkers_1d(from, to)
  type(AfqmcWalker), allocatable, intent(in)    :: from(:)
  type(AfqmcWalker), allocatable, intent(inout) :: to(:)
  !local variables
  integer                                       :: i

  if (.not. allocated(to)) allocate(to, source=from)

  do i = 1, size(from)
    to(i) = from(i)
  end do
end subroutine copy_afqmc_walkers_1d

subroutine copy_afqmc_walkers_2d(from, to)
  type(AfqmcWalker), allocatable, intent(in)    :: from(:,:)
  type(AfqmcWalker), allocatable, intent(inout) :: to(:,:)
  !local variables
  integer                                       :: i, j

  if (.not. allocated(to)) allocate(to, source=from)

  do j = 1, size(from, 2)
    do i = 1, size(from, 1)
      to(i,j) = from(i,j)
    end do
  end do
end subroutine copy_afqmc_walkers_2d


!******************************************************************************** 
!
!       Pack array of walkers into walker_block
!
!               Input:
!                   ham_des  - Hamiltonian descriptor
!                   walkers  - array of afqmc walkers
!                   mode     - determines the data layout - see below
!
!               Output:
!                   coeff    - packed orbitals
!
!   default - mode not present
!       -------------------------------------------------------------
!       |       |       |     |       |       |       |     |       |
!       |       |       |     |       |       |       |     |       |
!       | phi_1 | phi_2 | ... | phi_n | phi_1 | phi_2 | ... | phi_n |
!       |  _up  |  _up  |     |  _up  |  _dn  |  _dn  |     |  _dn  |
!       |       |       |     |       |       |       |     |       |
!       -------------------------------------------------------------
!
!   if present mode  
!       -------------------------------------------------------
!       |       |       |       |       |     |       |       |
!       |       |       |       |       |     |       |       |
!       | phi_1 | phi_1 | phi_2 | phi_2 | ... | phi_n | phi_n |
!       |  _up  |  _dn  |  _up  |  _dn  |     |  _up  |  _dn  |
!       |       |       |       |       |     |       |       |
!       -------------------------------------------------------
!   
!******************************************************************************** 
subroutine batch_walkers_1d(ham_des, walkers, coeff, mode)
  type(hamil_descriptor), intent(in)      :: ham_des
  type(AfqmcWalker), intent(in)           :: walkers(:)
  complex(wp), allocatable, intent(inout) :: coeff(:,:) 
  integer, optional, intent(in)           :: mode
  !local variables
  integer                                 :: nw, w, w1, w2, w1_, w2_, spin, nocc
  character(len=*), parameter             :: proc_name = "batch_walkers"
  
  if (use_profiling) call start_profiling(proc_name)

  nw = size(walkers)
  if (.not. allocated(coeff)) allocate(coeff(ham_des%n, ham_des%nocc*nw))
 
  if (present(mode)) then
    w1_ = 0
    w2_ = 0
    do w = 1, nw
      do spin = 1, ham_des%ispin
        nocc = ham_des%nel(spin)
        w1 = ham_des%nell(spin) + 1
        w2 = ham_des%nell(spin) + nocc
        w1_ = w2_ + 1
        w2_ = w2_ + nocc
        coeff(:,w1_:w2_) = walkers(w)%coeff(:,w1:w2)
      end do
    end do
  else
    w1_ = 0
    w2_ = 0
    do spin = 1, ham_des%ispin
      nocc = ham_des%nel(spin)
      w1 = ham_des%nell(spin) + 1
      w2 = ham_des%nell(spin) + nocc
      do w = 1, nw
        w1_ = w2_ + 1
        w2_ = w2_ + nocc
        coeff(:,w1_:w2_) = walkers(w)%coeff(:,w1:w2)
      end do
    end do
  end if

  if (use_profiling) call end_profiling(proc_name)
end subroutine batch_walkers_1d

!******************************************************************************** 
!
! Batch orbitals of the AfqmcWalker array into a single array of orbitals
!
! Data layout:
!
!   default - mode not present - used for exchange energy evaluation: 
!       -------------------------------------------------------------
!       |       |       |     |       |       |       |     |       |
!       |       |       |     |       |       |       |     |       |
!       | phi_1 | phi_2 | ... | phi_n | phi_1 | phi_2 | ... | phi_n |
!       |  _up  |  _up  |     |  _up  |  _dn  |  _dn  |     |  _dn  |
!       |       |       |     |       |       |       |     |       |
!       -------------------------------------------------------------
!
!   if present mode - used for force bias / hartree energy calculation       
!       -------------------------------------------------------
!       |       |       |       |       |     |       |       |
!       |       |       |       |       |     |       |       |
!       | phi_1 | phi_1 | phi_2 | phi_2 | ... | phi_n | phi_n |
!       |  _up  |  _dn  |  _up  |  _dn  |     |  _up  |  _dn  |
!       |       |       |       |       |     |       |       |
!       -------------------------------------------------------
!   
!******************************************************************************** 
subroutine batch_walkers_2d(hdes, walkers, coeff, mode)
  type(hamil_descriptor), intent(in)      :: hdes
  type(AfqmcWalker), intent(in)           :: walkers(:,:)
  complex(wp), allocatable, intent(inout) :: coeff(:,:) 
  integer, optional, intent(in)           :: mode
  !local variables
  integer                                 :: i, w, w1, w2, w1_, w2_, spin, nocc
  character(len=*), parameter             :: proc_name = "batch_walkers"
  
  if (use_profiling) call start_profiling(proc_name)

  if (.not. allocated(coeff)) allocate(coeff(hdes%n, hdes%nocc*size(walkers)))
 
  if (present(mode)) then
    w1_ = 0
    w2_ = 0
    do i = 1, size(walkers, 2)
      do w = 1, size(walkers, 1)
        do spin = 1, hdes%ispin
          nocc = hdes%nel(spin)
          w1 = hdes%nell(spin) + 1
          w2 = hdes%nell(spin) + nocc
          w1_ = w2_ + 1
          w2_ = w2_ + nocc
          coeff(:,w1_:w2_) = walkers(w,i)%coeff(:,w1:w2)
        end do
      end do
    end do
  else
    w1_ = 0
    w2_ = 0
    do spin = 1, hdes%ispin
      nocc = hdes%nel(spin)
      w1 = hdes%nell(spin) + 1
      w2 = hdes%nell(spin) + nocc
      do i = 1, size(walkers, 2)
        do w = 1, size(walkers, 1)
          w1_ = w2_ + 1
          w2_ = w2_ + nocc
          coeff(:,w1_:w2_) = walkers(w,i)%coeff(:,w1:w2)
        end do
      end do
    end do
  end if

  if (use_profiling) call end_profiling(proc_name)
end subroutine batch_walkers_2d


!******************************************************************************** 
!
! Unpack AfqmcWalker batch into the array of AfqmcWalker
!
!******************************************************************************** 
subroutine unbatch_walkers_1d(walkers, ovlp, lgmean, en)
  type(AfqmcWalker), intent(inout)    :: walkers(:)
  complex(wp), intent(in)             :: ovlp(:)
  complex(wp), intent(in)             :: lgmean(:,:)
  type(CEnergy), optional, intent(in) :: en(:)
  !local variables
  integer                             :: w, nw
  character(len=*), parameter         :: proc_name = "unbatch_walkers"
  
  if (use_profiling) call start_profiling(proc_name)

  nw = size(walkers)
  
  do w = 1, nw
    walkers(w)%overlap_old = walkers(w)%overlap
    walkers(w)%overlap = ovlp(w)
    walkers(w)%lgmean = lgmean(:,w)
  end do

  if (present(en)) then
    do w = 1, nw
      walkers(w)%energy = en(w)
    end do
  end if

  if (use_profiling) call end_profiling(proc_name)
end subroutine unbatch_walkers_1d

subroutine unbatch_walkers_2d(walkers, ovlp, lgmean, en)
  type(AfqmcWalker), intent(inout)    :: walkers(:,:)
  complex(wp), intent(in)             :: ovlp(:)
  complex(wp), intent(in)             :: lgmean(:,:)
  type(CEnergy), optional, intent(in) :: en(:)
  !local variables
  integer                             :: i, w, w_
  character(len=*), parameter         :: proc_name = "unbatch_walkers"
  
  if (use_profiling) call start_profiling(proc_name)
  
  do i = 1, size(walkers, 2)
    do w = 1, size(walkers, 1)
      w_ = (i-1) * size(walkers, 1) + w
      walkers(w,i)%overlap_old = walkers(w,i)%overlap
      walkers(w,i)%overlap = ovlp(w_)
      walkers(w,i)%lgmean = lgmean(:,w_)
    end do
  end do

  if (present(en)) then
    do i = 1, size(walkers, 2)
      do w = 1, size(walkers, 1)
        w_ = (i-1) * size(walkers, 1) + w
        walkers(w,i)%energy = en(w_)
      end do
    end do
  end if

  if (use_profiling) call end_profiling(proc_name)
end subroutine unbatch_walkers_2d


!******************************************************************************** 
!
! Unbatch AfqmcWalker batch into the array of AfqmcWalker (g-resolved version)
!
!******************************************************************************** 
subroutine unbatch_walkers_gres_1d(walkers, ovlp, lgmean, eh, ex, eself, en)
  type(AfqmcWalker), intent(inout)    :: walkers(:)
  complex(wp), intent(in)             :: ovlp(:)
  complex(wp), intent(in)             :: lgmean(:,:)
  complex(wp), intent(in)             :: eh(:,:)
  complex(wp), intent(in)             :: ex(:,:)
  complex(wp), intent(in)             :: eself(:,:)
  type(CEnergy), optional, intent(in) :: en(:)
  !local variables
  integer                             :: w, nw
  character(len=*), parameter         :: proc_name = "unbatch_walkers"
  
  if (use_profiling) call start_profiling(proc_name)

  nw = size(walkers)
  
  do w = 1, nw
    walkers(w)%overlap_old = walkers(w)%overlap
    walkers(w)%overlap = ovlp(w)
    walkers(w)%lgmean = lgmean(:,w)
    walkers(w)%eh = eh(:,w)
    walkers(w)%ex = ex(:,w)
    walkers(w)%eself = eself(:,w)
  end do

  if (present(en)) then
    do w = 1, nw
      walkers(w)%energy = en(w)
    end do
  end if

  if (use_profiling) call end_profiling(proc_name)
end subroutine unbatch_walkers_gres_1d

subroutine unbatch_walkers_gres_2d(walkers, ovlp, lgmean, eh, ex, eself, en)
  type(AfqmcWalker), intent(inout)    :: walkers(:,:)
  complex(wp), intent(in)             :: ovlp(:)
  complex(wp), intent(in)             :: lgmean(:,:)
  complex(wp), intent(in)             :: eh(:,:)
  complex(wp), intent(in)             :: ex(:,:)
  complex(wp), intent(in)             :: eself(:,:)
  type(CEnergy), optional, intent(in) :: en(:)
  !local variables
  integer                             :: i, w, w_
  character(len=*), parameter         :: proc_name = "unbatch_walkers"
  
  if (use_profiling) call start_profiling(proc_name)
  
  do i = 1, size(walkers, 2)
    do w = 1, size(walkers, 1)
      w_ = (i-1) * size(walkers, 1) + w
      walkers(w,i)%overlap_old = walkers(w,i)%overlap
      walkers(w,i)%overlap = ovlp(w_)
      walkers(w,i)%lgmean = lgmean(:,w_)
      walkers(w,i)%eh = eh(:,w_)
      walkers(w,i)%ex = ex(:,w_)
      walkers(w,i)%eself = eself(:,w_)
    end do
  end do

  if (present(en)) then
    do i = 1, size(walkers, 2)
      do w = 1, size(walkers, 1)
        w_ = (i-1) * size(walkers, 1) + w
        walkers(w,i)%energy = en(w_)
      end do
    end do
  end if

  if (use_profiling) call end_profiling(proc_name)
end subroutine unbatch_walkers_gres_2d


!******************************************************************************** 
!
! MPI Send for the AfqmcWalker object
!
!debug: cleanup needed to see what is really needed to be send around
!******************************************************************************** 
subroutine send_afqmc_walker(walker, comm, dest, tag)
  type(AfqmcWalker), intent(in)      :: walker
  type(mpi_communicator), intent(in) :: comm
  integer, intent(in)                :: dest
  integer, intent(in)                :: tag
  
  call comm%send(walker%id, dest, tag)
  call comm%send(walker%lgmean, dest, tag)
  call comm%send(walker%eh, dest, tag)
  call comm%send(walker%ex, dest, tag)
  call comm%send(walker%eself, dest, tag)
  call comm%send(walker%coeff, dest, tag)

  call comm%send(walker%energy%enuc, dest, tag)
  call comm%send(walker%energy%e0, dest, tag)
  call comm%send(walker%energy%e1, dest, tag)
  call comm%send(walker%energy%eh, dest, tag)
  call comm%send(walker%energy%ex, dest, tag)
  call comm%send(walker%overlap, dest, tag)
  call comm%send(walker%weight, dest, tag)
  call comm%send(walker%bare_weight, dest, tag)
  call comm%send(walker%phase, dest, tag)
  call comm%send(walker%energyl, dest, tag)
  call comm%send(walker%energyh, dest, tag)
  call comm%send(walker%life, dest, tag)
  call comm%send(walker%revlife, dest, tag)
  call comm%send(walker%move, dest, tag)
  call comm%send(walker%acc, dest, tag)

  !probably not necessary for communication
  call comm%send(walker%doverlap, dest, tag)
  call comm%send(walker%rew_w, dest, tag)
  call comm%send(walker%rew_p, dest, tag)
  call comm%send(walker%overlap_old, dest, tag)
  call comm%send(walker%energyl_old, dest, tag)
  call comm%send(walker%energyh_old, dest, tag)
  call comm%send(walker%isf_mf, dest, tag)
  call comm%send(walker%isf, dest, tag)
  call comm%send(walker%dphase, dest, tag)
end subroutine send_afqmc_walker


!******************************************************************************** 
!
! MPI Recv for the AfqmcWalker object
!
!******************************************************************************** 
subroutine recv_afqmc_walker(walker, comm, source, tag)
  type(AfqmcWalker), intent(inout)   :: walker
  type(mpi_communicator), intent(in) :: comm
  integer, intent(in)                :: source
  integer, intent(in)                :: tag
  
  call comm%recv(walker%id, source, tag)
  call comm%recv(walker%lgmean, source, tag)
  call comm%recv(walker%eh, source, tag)
  call comm%recv(walker%ex, source, tag)
  call comm%recv(walker%eself, source, tag)
  call comm%recv(walker%coeff, source, tag)

  call comm%recv(walker%energy%enuc, source, tag)
  call comm%recv(walker%energy%e0, source, tag)
  call comm%recv(walker%energy%e1, source, tag)
  call comm%recv(walker%energy%eh, source, tag)
  call comm%recv(walker%energy%ex, source, tag)
  call comm%recv(walker%overlap, source, tag)
  call comm%recv(walker%weight, source, tag)
  call comm%recv(walker%bare_weight, source, tag)
  call comm%recv(walker%phase, source, tag)
  call comm%recv(walker%energyl, source, tag)
  call comm%recv(walker%energyh, source, tag)
  call comm%recv(walker%life, source, tag)
  call comm%recv(walker%revlife, source, tag)
  call comm%recv(walker%move, source, tag)
  call comm%recv(walker%acc, source, tag)

  !probably not necessary for communication
  call comm%recv(walker%doverlap, source, tag)
  call comm%recv(walker%rew_w, source, tag)
  call comm%recv(walker%rew_p, source, tag)
  call comm%recv(walker%overlap_old, source, tag)
  call comm%recv(walker%energyl_old, source, tag)
  call comm%recv(walker%energyh_old, source, tag)
  call comm%recv(walker%isf_mf, source, tag)
  call comm%recv(walker%isf, source, tag)
  call comm%recv(walker%dphase, source, tag)
end subroutine recv_afqmc_walker


!******************************************************************************** 
!
! Broadcast AfqmcWalker object
!
!******************************************************************************** 
subroutine broadcast_afqmc_walker(walker, comm, root)
  type(AfqmcWalker), intent(inout)   :: walker
  type(mpi_communicator), intent(in) :: comm
  integer, intent(in)                :: root
  !local variables
  integer                            :: i

  if (comm%mpirank == root) then
    do i = 0, comm%mpisize-1
      if (i == root) cycle
      call send_afqmc_walker(walker, comm, i, i)
    end do
  else
    call recv_afqmc_walker(walker, comm, root, comm%mpirank)
  end if
end subroutine broadcast_afqmc_walker


!******************************************************************************** 
!
! Gather AfqmcWalker array
!
! Note:
!    Number of walkers may differ per MPI process
!
!******************************************************************************** 
subroutine gather_walkers(comm, hdes, from, to, root)
  type(mpi_communicator), intent(in)          :: comm
  type(hamil_descriptor), intent(in)          :: hdes
  type(AfqmcWalker), intent(in)               :: from(:)
  type(AfqmcWalker), allocatable, intent(out) :: to(:)
  integer, optional, intent(in)               :: root
  !local variables
  integer                                     :: w, w_, root_, source
  integer, allocatable                        :: nw_per_proc(:)

  if (present(root)) then
    root_ = root
  else
    root_ = 0
  end if

  call comm%gather(size(from), nw_per_proc, root_)

  if (comm%mpirank == root_) then
    call alloc_walker(to, hdes, sum(nw_per_proc))

    w_ = 0
    do source = 1, size(nw_per_proc)
      do w = 1, nw_per_proc(source)
        w_ = w_ + 1
        if (source-1 == root_) then
          to(w_) = from(w)
        else
          call recv_afqmc_walker(to(w_), comm, source-1, 1000)
        end if
      end do
    end do
  else
    do w = 1, size(from)
      call send_afqmc_walker(from(w), comm, root_, 1000)
    end do
  end if
end subroutine gather_walkers


!******************************************************************************** 
!
! Scatter AfqmcWalker array
!
!******************************************************************************** 
subroutine scatter_walkers(comm, from, to)
  type(mpi_communicator), intent(in) :: comm
  type(AfqmcWalker), intent(in)      :: from(:)
  type(AfqmcWalker), intent(inout)   :: to(:)
  !local variables
  integer                            :: w, w_, dest, tag

  if (comm%mpirank == 0) then
    do w = 1, size(from)
      dest = (w-1) / size(to)
      w_  = w - dest*size(to)
      if (dest == 0) then
        to(w_) = from(w)
      else
        tag = 1000 + dest
        call send_afqmc_walker(from(w), comm, dest, tag)
      end if
    end do
  else
    do w = 1, size(to)
      !debug: very strange i'm not sure how it worked
      tag = 1000 + comm%mpirank
      call recv_afqmc_walker(to(w), comm, 0, tag)
    end do
  end if
end subroutine scatter_walkers

end module afqmc_walker