! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-FileCopyrightText: Copyright (c) 2025 Martin Schlipf
! SPDX-License-Identifier: Apache-2.0

module afqmc_rebalance
!******************************************************************************** 
!
! Rebalance routines for AFQMC walkers
!
!******************************************************************************** 

#include "../preproc.inc"
use constants
use mpi
use profiling
use standalone, only: modul
use brng_factories, only: brng_glb
use afqmc_walker

implicit none

integer, parameter, private :: free = -1

type pop_controll
  logical                                     :: lnorm = .false.
  character(len=charlen)                      :: method = "comb"
  real(wp)                                    :: min_weight = 0.1_wp
  real(wp)                                    :: max_weight = 3.0_wp
  integer                                     :: freq = unset
  integer                                     :: count_ 
  integer                                     :: sent_walkers, max_walkers
  integer                                     :: sent_walkers_tot, max_walkers_tot
  integer                                     :: nwalkers
  integer, allocatable                        :: map(:,:)
  procedure(get_instances_interface), pointer :: get_instances => null()
contains
  procedure                                   :: init => init_pop_controll
  procedure                                   :: reader => pop_controll_reader
  procedure                                   :: average => pop_controll_average
  procedure                                   :: check_weights
  procedure                                   :: rebalance
  procedure                                   :: rebalance_walkers
end type pop_controll

abstract interface
  subroutine get_instances_interface(pop, weights, instances, bias)
    import wp, pop_controll
    class(pop_controll), intent(in)   :: pop
    real(wp), intent(in)              :: weights(:,:)
    integer, allocatable, intent(out) :: instances(:,:)
    real(wp), optional, intent(in)    :: bias
  end subroutine get_instances_interface
end interface

contains

!******************************************************************************** 
!
!       Population control - interface routine used in afqmc module
!
!******************************************************************************** 
subroutine rebalance_walkers(pop, walkers, step, lpop)
  class(pop_controll), intent(inout) :: pop
  type(AfqmcWalker), intent(inout)   :: walkers(:)
  integer, intent(in)                :: step
  logical, intent(out)               :: lpop
  
  lpop = .false.

  if (pop%freq == 0) return

  if (pop%freq > 0) then
    if (modul(step, pop%freq)) lpop = .true.
  else
    call pop%check_weights(walkers%weight, lpop)
  end if

  if (.not. lpop) return
  pop%count_ = pop%count_ + 1
  
  call pop%rebalance(comm_world, walkers)
  
  if (comm_world%mpirank == 0) then
    pop%sent_walkers_tot = pop%sent_walkers_tot + pop%sent_walkers
    pop%max_walkers_tot = max(pop%max_walkers_tot, pop%max_walkers)
  end if
end subroutine rebalance_walkers


!******************************************************************************** 
!
!       Same as above, but MPI version for full walker ensemble
!
!       Implemented by Z. Sukruma and M. Schlipf
!
!******************************************************************************** 
subroutine rebalance(pop, comm, walkers, bias)
  class(pop_controll), intent(inout) :: pop
  type(mpi_communicator), intent(in) :: comm
  type(AfqmcWalker), intent(inout)   :: walkers(:)
  real(wp), intent(in), optional     :: bias
  !local variables
  real(wp)                           :: norm
  real(wp), allocatable              :: weights(:,:)
  integer, allocatable               :: instances(:,:), map(:,:)
  character(len=*), parameter        :: proc_name = "rebalance_walkers"
  
  if (use_profiling) call start_profiling(proc_name)

  call gather_weights(comm, walkers%weight, weights)

  if (comm%mpirank == 0) then
    call normalize_weights(weights, norm)
    call pop%get_instances(weights, instances, bias)
    call determine_map(instances, map)
    pop%sent_walkers = sum(abs(instances-1)) / 2 
    pop%max_walkers = maxval(instances-1)   
  end if

  call broadcast_map(comm, size(walkers), norm, map)
  call copy_walkers(comm, map, walkers)

  if (pop%lnorm) then
    walkers%weight = norm
  else
    walkers%weight = 1.0_wp 
  end if

  if (comm%mpirank == 0) then
    if (.not. allocated(pop%map)) allocate(pop%map, mold=map)
    pop%map = map
  end if

  if (use_profiling) call end_profiling(proc_name)
end subroutine rebalance


!******************************************************************************** 
!
!       Collect all weigths on the master rank
!
!******************************************************************************** 
subroutine gather_weights(comm, local_weights, global_weights)
  type(mpi_communicator), intent(in) :: comm
  real(wp), intent(in)               :: local_weights(:)
  real(wp), allocatable, intent(out) :: global_weights(:,:)
  !local
  real(wp), allocatable              :: local_weights_2d(:,:)
  
  if (comm%mpirank == 0) then
    allocate(global_weights(size(local_weights), comm%mpisize))
  else 
    allocate(global_weights(1,1))
  end if

  allocate(local_weights_2d(size(local_weights), 1))
  local_weights_2d(:,1) = local_weights

  call comm_world%gather(local_weights_2d, global_weights, 0)
end subroutine gather_weights


!******************************************************************************** 
!
!       Calculate normalized weights
!
!******************************************************************************** 
subroutine normalize_weights(weights, norm)
  real(wp), intent(inout) :: weights(:,:)
  real(wp), intent(out)   :: norm
  !
  norm = sum(weights) / size(weights)
  weights = weights / norm
end subroutine normalize_weights


!******************************************************************************** 
!
!       Number of replicas for each walker in std method
!
!******************************************************************************** 
subroutine get_instances_std(pop, weights, instances, bias)
  class(pop_controll), intent(in)   :: pop
  real(wp), intent(in)              :: weights(:,:)
  integer, allocatable, intent(out) :: instances(:,:)
  real(wp), intent(in), optional    :: bias
  !local variables
  real(wp), allocatable             :: cumulative_sum(:), r(:)
  real(wp)                          :: sum_
  integer                           :: walker, image, i, w
  
  allocate(instances(size(weights, 1), size(weights, 2)))
  allocate(cumulative_sum(0:size(weights)))
  allocate(r(size(weights)))
  instances = 0
  cumulative_sum = 0.0_wp
  sum_ = 0.0_wp
  do image = 1, size(weights, 2)
    do walker = 1, size(weights, 1)
      i = (image-1)*size(weights,1) + walker
      sum_ = sum_ + weights(walker, image)
      cumulative_sum(i) = sum_
    end do
  end do
  
  call brng_glb%rand(r)
  do w = 1, size(r)
    do i = 1, size(r)
      if (r(w)>=cumulative_sum(i-1) .and. r(w)<cumulative_sum(i)) then
        image = (i - 1) / size(weights, 1) + 1
        walker = mod(i - 1, size(weights, 1)) + 1
        instances(walker, image) = instances(walker, image) + 1
        exit 
      end if
    end do
  end do
end subroutine get_instances_std


!******************************************************************************** 
!
!       Number of replicas for each walker in comb method
!
!******************************************************************************** 
subroutine get_instances_comb(pop, weights, instances, bias)
  class(pop_controll), intent(in)   :: pop
  real(wp), intent(in)              :: weights(:,:)
  integer, allocatable, intent(out) :: instances(:,:)
  real(wp), intent(in), optional    :: bias
  !local variables
  real(wp)                          :: cumulative_sum
  integer                           :: walker, image, current, previous
 
  allocate(instances(size(weights, 1), size(weights, 2)))
  cumulative_sum = determine_bias(bias)
  previous = 0

  do image = 1, size(weights, 2)
    do walker = 1, size(weights, 1)
      cumulative_sum = cumulative_sum + weights(walker, image)
      current = ceiling(cumulative_sum)
      instances(walker, image) = current - previous
      previous = current
    end do
  end do
end subroutine get_instances_comb


!******************************************************************************** 
!
!       Bias for the first walker in series
!
!******************************************************************************** 
real(wp) function determine_bias(bias) result (bias_)
  real(wp), intent(in), optional :: bias
  
  if (present(bias)) then
    bias_ = bias
  else
    call brng_glb%rand(bias_)
    bias_ = -1.0_wp * bias_
  end if
end function determine_bias


!******************************************************************************** 
!
!       Determine the id-s of new walkers form old walkers
!
!                     walkers_new = walkers[map]
!
!******************************************************************************** 
subroutine determine_map(instances, map)
  integer, allocatable, intent(in)  :: instances(:,:)
  integer, allocatable, intent(out) :: map(:,:)
  !local variables
  integer, allocatable              :: remaining(:), free_slots(:), new_slots(:)
  integer                           :: image
  
  allocate(map, mold=instances)
  do image = 1, size(instances, 2)
    call determine_local_map(instances(:,image), map(:,image), image)
    free_slots = get_free_slots(instances(:,image))
    new_slots = get_new_slots(instances(:, image), image)
    call merge_free_new_slots(map(:,image), remaining, free_slots, new_slots)
    !call merge_new_slots(remaining, new_slots)
  end do
  if (size(remaining) > 0) map = unpack(remaining, map==free, map)
end subroutine determine_map


!******************************************************************************** 
!
!       Determine the map indices within the node
!
!******************************************************************************** 
subroutine determine_local_map(instances, local_map, image)
  integer, intent(in)    :: instances(:)
  integer, intent(inout) :: local_map(:)
  integer, intent(in)    :: image
  !local variables
  integer                :: w 
  
  local_map = free
  do w = 1, size(instances)
    if (instances(w) > 0) local_map(w) = w
  end do
  where (local_map /= free)
    local_map = local_map + (image - 1) * size(instances)
  end where
end subroutine determine_local_map


!******************************************************************************** 
!       
!       Determines the array of walkers which are killed - each entry
!       represents the id of a walker which has to be replaced by some 
!                           other walker in new_slots
!
!******************************************************************************** 
function get_free_slots(instances) result(free_slots)
    integer, intent(in)  :: instances(:)
    integer, allocatable :: free_slots(:)
    !local variables
    integer              :: w, i
    !
    if (allocated(free_slots)) deallocate(free_slots)
    allocate(free_slots(count(instances == 0)))
    !
    i = 0
    do w = 1, size(instances)
      if (instances(w) == 0) then
        i = i + 1
        free_slots(i) = w
      end if
    end do
end function get_free_slots


!******************************************************************************** 
!       
!       Determines the array of walkers which are replicated - each entry
!                            represents the id of a walker
!
!******************************************************************************** 
function get_new_slots(instances, image) result(new_slots)
  integer, intent(in)  :: instances(:)
  integer, intent(in)  :: image
  integer, allocatable :: new_slots(:)
  !local variables
  integer              :: w, i, ii
  !
  if (allocated(new_slots)) deallocate(new_slots)
  allocate(new_slots(sum(pack(instances-1, instances-1>0))))
  !
  ii = 0
  do w = 1, size(instances)
    do i = 2, instances(w)
      ii = ii + 1
      new_slots(ii) = w
    end do
  end do
  new_slots = new_slots + (image - 1) * size(instances)
end function get_new_slots


!******************************************************************************** 
!
!       Fills the empty slots within the local map with replicated walkers from
!       the same image. The rest of the free slots remains free and rest of the 
!       rest of the replicated slots goes to the remianing.
!       At the end, full global map is filled by entries from remaining 
!
!******************************************************************************** 
subroutine merge_free_new_slots(map, remaining, free_slots, new_slots)
  integer, intent(inout)              :: map(:)
  integer, allocatable, intent(inout) :: remaining(:)
  integer, intent(in)                 :: free_slots(:)
  integer, intent(in)                 :: new_slots(:)
  !local variables
  integer, allocatable                :: work(:)
  integer                             :: i, ii
  !
  do i = 1, min(size(free_slots), size(new_slots))
    map(free_slots(i)) = new_slots(i)
  end do
  !
  if (size(free_slots) >= size(new_slots)) return
  if (.not. allocated(remaining)) allocate(remaining(0))
  call move_alloc(remaining, work)
  allocate(remaining(size(work) + size(new_slots) - size(free_slots)))
  remaining(1:size(work)) = work
  !
  do i = size(free_slots)+1, size(new_slots)
    ii = size(work) + i - size(free_slots)
    remaining(ii) = new_slots(i)
  end do
end subroutine merge_free_new_slots


!******************************************************************************** 
!
!       Copies new slots from the local map to the remainging
!
!******************************************************************************** 
subroutine merge_new_slots(remaining, new_slots)
  integer, allocatable, intent(inout) :: remaining(:)
  integer, intent(in)                 :: new_slots(:)
  !local variables
  integer, allocatable                :: work(:)
  integer                             :: i, ii
  !
  call move_alloc(remaining, work)
  allocate(remaining(size(work) + size(new_slots)))
  remaining(1:size(work)) = work
  do i = 1, size(new_slots) 
    ii = size(work) + i
    remaining(ii) = new_slots(i)
  end do
end subroutine merge_new_slots


!******************************************************************************** 
!
!       Broadcast map to all nodes
!
!******************************************************************************** 
subroutine broadcast_map(comm, num_walkers, norm, map)
  type(mpi_communicator), intent(in)  :: comm
  integer, intent(in)                 :: num_walkers
  real(wp), intent(inout)             :: norm
  integer, allocatable, intent(inout) :: map(:,:)
  
  if (comm%mpirank /= 0) allocate(map(num_walkers, comm%mpisize))

  call comm%bcast(map, 0)
  call comm%bcast(norm, 0)
end subroutine broadcast_map

!******************************************************************************** 
!
!       Walker copies according to the map
!
!******************************************************************************** 
subroutine copy_walkers(comm, map, walkers)
  type(mpi_communicator), intent(in) :: comm
  integer, intent(in)                :: map(:,:)
  type(AfqmcWalker), intent(inout)   :: walkers(:)
  !local variables
  integer, allocatable               :: from_image(:), from_walker(:), map_image(:)
  
  call determine_from_where(map(:,comm%mpinode), from_image, from_walker)
  call copy_locally(from_image == comm%mpinode, from_walker, walkers)
  call receive_from_other_images(comm, from_image, walkers)
  call send_to_other_images(comm, map, walkers)
end subroutine copy_walkers


!******************************************************************************** 
!
!       Determine who sends the walker       
!
!******************************************************************************** 
subroutine determine_from_where(map_image, from_image, from_walker)
  integer, intent(in)               :: map_image(:)
  integer, allocatable, intent(out) :: from_image(:)
  integer, allocatable, intent(out) :: from_walker(:)
  !
  allocate(from_image(size(map_image)), from_walker(size(map_image)))
  from_image = (map_image - 1) / size(map_image) + 1
  from_walker = mod(map_image - 1, size(map_image)) + 1
end subroutine determine_from_where


!******************************************************************************** 
!
!       Copy walker withing the node
!       
!******************************************************************************** 
subroutine copy_locally(local_copy, from_walker, walkers)
  logical, intent(in)              :: local_copy(:)
  integer, intent(in)              :: from_walker(:)
  type(AfqmcWalker), intent(inout) :: walkers(:)
  !local variables
  integer                           :: walker
  !
  do walker = 1, size(walkers)
    if (.not. local_copy(walker)) cycle
    if (walker == from_walker(walker)) cycle
    walkers(walker) = walkers(from_walker(walker))
  end do
end subroutine copy_locally


!******************************************************************************** 
!
!       Routine to receive walker variables to other nodes
!
!******************************************************************************** 
subroutine receive_from_other_images(comm, from_image, walkers)
  type(mpi_communicator), intent(in) :: comm
  integer, intent(in)                :: from_image(:)
  type(AfqmcWalker), intent(inout)   :: walkers(:)
  !local variables
  integer                            :: walker
  !
  do walker = 1, size(walkers)
    if (from_image(walker) == comm%mpinode) cycle
    call recv_afqmc_walker(walkers(walker), comm, from_image(walker)-1, 200)
  end do
end subroutine receive_from_other_images


!******************************************************************************** 
!
!       Routine to send walker variables to other nodes
!
!******************************************************************************** 
subroutine send_to_other_images(comm, map, walkers)
  type(mpi_communicator), intent(in) :: comm
  integer, intent(in)                :: map(:,:)
  type(AfqmcWalker), intent(inout)   :: walkers(:)
  !local variables
  integer, allocatable               :: from_image(:), from_walker(:), map_image(:)
  integer                            :: image, walker, lb, ub
  
  do image = 1, comm%mpisize
    if (image == comm%mpinode) cycle
    call determine_from_where(map(:,image), from_image, from_walker)
    do walker = 1, size(walkers)
      if (from_image(walker) /= comm%mpinode) cycle
      call send_afqmc_walker(walkers(from_walker(walker)), comm, image-1, 200)
    end do
  end do
end subroutine send_to_other_images


!******************************************************************************** 
!
!       Initialize pop_controll type
!
!******************************************************************************** 
subroutine init_pop_controll(pop, nwalkers)
  class(pop_controll), intent(inout) :: pop
  integer, intent(in)                :: nwalkers
  
  call pop%reader
  
  call associate_get_instances(pop)
  
  pop%nwalkers = nwalkers
  pop%count_ = 0
  pop%sent_walkers = 0
  pop%sent_walkers_tot = 0
  pop%max_walkers = 0
  pop%max_walkers_tot = 0
end subroutine init_pop_controll


!******************************************************************************** 
!
!       Small routine to read pop_controll variables from qmcfort_in file
!
!******************************************************************************** 
subroutine pop_controll_reader(pop)
  use qmcfort_in, only: add_input
  class(pop_controll), intent(inout) :: pop
  !
  call add_input("pop_method",       pop%method)
  call add_input("rebalance_method", pop%method)
  call add_input("rebalance",        pop%freq)
  call add_input("pop",              pop%freq) 
  call add_input("pop_norm",         pop%lnorm)
  call add_input("rebalance_norm",   pop%lnorm)
  call add_input("min_weight",       pop%min_weight)
  call add_input("max_weight",       pop%max_weight)
  call add_input("pop_min_weight",   pop%min_weight)
  call add_input("pop_max_weight",   pop%max_weight)
end subroutine pop_controll_reader


!******************************************************************************** 
!
!       Choose specified get_instances method
!
!******************************************************************************** 
subroutine associate_get_instances(pop)
  type(pop_controll), intent(inout) :: pop
  
  select case (trim(pop%method))
    case ("comb")
      pop%get_instances => get_instances_comb
    case ("std")
      pop%get_instances => get_instances_std
    case default
      pop%method = "comb"
      pop%get_instances => get_instances_comb
  end select
end subroutine associate_get_instances


!******************************************************************************** 
!
!   Caluclate average number of walkers sent around per call
!
!******************************************************************************** 
integer function pop_controll_average(pop)
  class(pop_controll), intent(in) :: pop
  if (pop%count_ == 0) then
    pop_controll_average = 0.0_wp
  else
    pop_controll_average = pop%sent_walkers_tot / pop%count_
  end if
end function pop_controll_average


!******************************************************************************** 
!
!       Check weights
!
!******************************************************************************** 
subroutine check_weights(pop, weights, lpop)
  class(pop_controll), intent(inout) :: pop
  real(wp), intent(in)               :: weights(:)
  logical, intent(inout)             :: lpop
  !local variables
  real(wp)                           :: min_weight, max_weight, av_weight
  integer                            :: ierror
  character(len=*), parameter        :: proc_name = "check_weights"
  
  if (use_profiling) call start_profiling(proc_name)
  
  min_weight = minval(weights)
  max_weight = maxval(weights)
  av_weight = sum(weights)

  call comm_world%reduce(min_weight, 0, MPI_Min)
  call comm_world%reduce(max_weight, 0, MPI_Max)
  call comm_world%mpisum(av_weight, 0)
  
  av_weight = av_weight / pop%nwalkers
  lpop = max_weight>pop%max_weight .or. min_weight<pop%min_weight

  call comm_world%bcast(lpop, 0)

  if (use_profiling) call end_profiling(proc_name)
end subroutine check_weights

end module afqmc_rebalance
