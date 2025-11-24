! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module profiling
!******************************************************************************** 
!       This module implemetns small routines for time profiling
!       Use:
!           start_profiling(name)
!           end_profiling(name)
!
!       Uses module variables
!
!       Linked (tree) list structure for structured timers:
!               Routine name with its content - displays total time and number of calls
!
!       Linked list for total timers:
!               Routine name    total time      total number of calls 
!
!******************************************************************************** 

#include "preproc.inc"
!$ use omp_lib
#ifdef MPI
use mpi_f08
#endif
use constants
use mpi, only: comm_world
use qmcfort_io

implicit none

type timer 
  character(len=50)    :: name_
  integer              :: ncalls
  integer              :: ncalls_total
  integer              :: nest_level

  real(dp)             :: start
  real(dp)             :: start_cpu
  real(dp)             :: time
  real(dp)             :: time_self
  real(dp)             :: cputime
  real(dp)             :: cputime_self
  real(dp)             :: flops
  real(dp)             :: flops_total

  type(timer), pointer :: parent
  type(timer), pointer :: child
  type(timer), pointer :: next
end type timer

type(timer), pointer :: head
type(timer), pointer :: head_flat
type(timer), pointer :: last
type(timer), pointer :: self_timer
type(timer), pointer :: root_timer

contains

subroutine start_profiling(name_)
  character(len=*), intent(in)           :: name_
  !local variables
  real(dp)                            :: time1, time2
  
  if (.not. associated(root_timer)) return
      
  !$omp master
  call snap(time1)
  call find_new_timer(last, name_)
  last%start = time1
  call snap(time2)
  self_timer%time = self_timer%time + time2 - time1
  self_timer%cputime = self_timer%cputime + time2 - time1
  !$omp end master
end subroutine start_profiling


subroutine end_profiling(name_, m, n, p, opt)
  character(len=*), intent(in)           :: name_
  integer, optional, intent(in)          :: m
  integer, optional, intent(in)          :: n
  integer, optional, intent(in)          :: p
  character(len=*), optional, intent(in) :: opt
  !local variables
  real(dp)                            :: time1, time2
  integer                                :: nthreads = 1
  
  if (.not. associated(root_timer)) return
#ifdef _OPENMP
  nthreads = omp_get_num_threads()
#endif

  !$omp master
  call snap(time1)
  if (trim(last%name_) /= trim(name_)) write(*,*) "Caution: something strange happend with timers!"
  last%time = last%time + time1 - last%start
  last%cputime = last%cputime + (time1 - last%start)*nthreads
  last%ncalls = last%ncalls + 1
  last%ncalls_total = last%ncalls_total + nthreads 
  call update_flops(last, nthreads, m, n, p, opt)
  if (associated(last%parent)) last => last%parent
  call snap(time2)
  self_timer%time = self_timer%time + time2 - time1
  self_timer%cputime = self_timer%cputime + (time2-time1)*nthreads
  !$omp end master
end subroutine end_profiling


!******************************************************************************** 
!
!       Initializes root_timer, parent of all other timers and parent of head
!
!******************************************************************************** 
subroutine initialize_timers()
  real(dp) :: tim1, tim1_cpu

  if (.not. use_profiling) return

  !$omp master
  call snap(tim1)
  call snap_cpu(tim1_cpu)
  allocate(root_timer)
  root_timer%name_ = "QMCFORT"
  root_timer%ncalls = 0
  root_timer%nest_level = 0
  root_timer%flops = 0.0_dp
  root_timer%flops_total = 0.0_dp
  root_timer%start = tim1
  root_timer%start_cpu = tim1_cpu
  root_timer%parent => null()
  root_timer%child  => null()
  root_timer%next   => null()
  !$omp end master
end subroutine initialize_timers


!******************************************************************************** 
!
!       Finalizes root_timer, full program timing determined - print only in flat
!
!********************************************************************************
subroutine finalize_timers()
  real(dp) :: tim2, tim2_cpu
  integer     :: nthreads = 1

  if (.not. associated(root_timer)) return

#ifdef _OPENMP
  nthreads = omp_get_num_threads()
#endif

  !$omp master
  call snap(tim2)
  call snap_cpu(tim2_cpu)
  root_timer%ncalls = root_timer%ncalls + 1
  root_timer%ncalls_total = root_timer%ncalls_total + nthreads
  root_timer%time = tim2 - root_timer%start
  root_timer%cputime = tim2_cpu - root_timer%start_cpu
  !$omp end master

  call dump_timings
end subroutine finalize_timers


!******************************************************************************** 
!
!       Initializes new timer
!               note: %start is initialized separately
!
!******************************************************************************** 
subroutine initialize_timer(timer_new, name_, parent, child, next)
  type(timer), pointer           :: timer_new
  character(len=*)               :: name_
  type(timer), pointer, optional :: parent
  type(timer), pointer, optional :: child
  type(timer), pointer, optional :: next
  !local variables
  type(timer), pointer           :: my_parent, my_child, my_next
  integer                        :: nest_level
  
  if (present(parent)) then
    my_parent => parent
    nest_level = parent%nest_level + 1
  else
    my_parent => null()
    nest_level = 0
  end if
  
  if (present(child)) then
    my_child => child
  else
    my_child => null()
  end if
  
  if (present(next)) then
    my_next => next
  else
    my_next => null()
  end if
  
  if (.not. associated(timer_new)) allocate(timer_new)

  timer_new%name_ = name_
  timer_new%ncalls = 0
  timer_new%ncalls_total = 0
  timer_new%nest_level = nest_level
  timer_new%time = 0.0_dp
  timer_new%time_self = 0.0_dp
  timer_new%cputime = 0.0_dp
  timer_new%cputime_self = 0.0_dp
  timer_new%flops = 0.0_dp
  timer_new%flops_total = 0.0_dp
  timer_new%parent => my_parent
  timer_new%child => my_child
  timer_new%next => my_next
end subroutine initialize_timer


!********************************************************************************  
!
!       Move pointer last to appropriate place after start_profiling call:
!
!               1. allocates head timer - first timer
!               2. allocates first child of the pointer last
!               3. allocates next child last%children
!               4. Reactivates one child among last%children
!
!******************************************************************************** 
subroutine find_new_timer(last, name_)
  type(timer), pointer :: last
  character(len=*), intent(in) :: name_
  !local variables
  type(timer), pointer :: child => null()
  
  if (.not. associated(head)) then
    !allocate first timer and let head show to it
    allocate(head)
    call  initialize_timer(head, name_, parent=root_timer)
    root_timer%child => head
    last => head
    !allocate self timer, too
    allocate(self_timer)
    call initialize_timer(self_timer, "self_timer")
  else
    if (associated(last%child)) then
      !Check for the first child
      child => last%child
      if (trim(child%name_) == trim(name_)) then
        last => child
      else
        do
          !Search among other children
          if (associated(child%next)) then
            child => child%next
            if (trim(child%name_) == trim(name_)) then
              last => child
              exit
            end if
          else
            !if matching timer is not found, allocate new (next) child
            allocate(child%next)
            child => child%next
            call initialize_timer(child, name_, parent=last)
            last => child
            exit
          end if
        end do
      end if
    else
      !allocate first child, if it still doesn't exist
      allocate(last%child)
      child => last%child
      call initialize_timer(child, name_, parent=last)
      last => child
    end if
  end if
end subroutine find_new_timer


!******************************************************************************** 
!
!       Updates the nummber of FLOPS of the specific timer
!
!   Either adds passedd number of FLOPS or assumes the FLOPS of the 
!   assumed matrix-matrix multiplication (A[m,n] * B[n,p])
!
!******************************************************************************** 
subroutine update_flops(self, nthreads, m, n, p, opt)
  type(timer), pointer, intent(inout)    :: self
  integer, intent(in)                    :: nthreads
  integer, optional, intent(in)          :: m
  integer, optional, intent(in)          :: n
  integer, optional, intent(in)          :: p
  character(len=*), optional, intent(in) :: opt
  !local variables
  real(dp)                            :: flops
  
  if (present(m)) then
    if (present(n) .and. present(p)) then
      flops = 1.0_dp * m * n * p * get_prefactor(opt)
    else
      flops = 1.0_dp * m
    end if
    self%flops = self%flops + flops
    self%flops_total = self%flops_total + flops * nthreads
  else
    if (self%flops > 0.1_dp) write(*,*) "CAUTION: Performance measurements could be wrong"
  end if
contains
  real(dp) function get_prefactor(opt)
    use string, only: lowercase
    character(len=*), optional, intent(in) :: opt
    
    select case (trim(lowercase(opt)))
      case ("r")
        get_prefactor = 2.0_dp
      case ("cr", "rc")
        get_prefactor = 4.0_dp
      case("c")
        get_prefactor = 8.0_dp
      case default
        get_prefactor = 2.0_dp
    end select
  end function get_prefactor
end subroutine update_flops


!******************************************************************************** 
!
!       Main output routine from this module
!
!               1. writes structured timers
!               2. writes flat timers
!
!******************************************************************************** 
subroutine dump_timings
  real(dp)          :: time1, time2
  
  call snap(time1)

  if (.not. associated(head)) then
    write(*,*) "Problem with time profiling occured - timers are not allocated!"
    return
  end if

!debug: NOCI - since NOCI is done only on master node and calls time profiling
!reduce_timers will fail
  !call reduce_timers(root_timer)
  call get_self_timers(root_timer)
  call decorate_self_timer

  if (comm_world%mpirank /= 0) return

  call create_flat_timers(root_timer, head_flat)
  call sort_timers(head_flat)
  
  call io%timing%open(status="replace", action="write")

  call graph_timer_header
  call write_graph_timers(head, root_timer)
  write(io%timing%funit,*) starlong

  call graph_frac_header
  call write_graph_fracs(head, root_timer)
  write(io%timing%funit,*) starlong

  call flat_timer_header
  call write_flat_timers(head_flat)
  write(io%timing%funit,*) starlong

  call flops_header
  call write_flops_items(head_flat)
  write(io%timing%funit,*) starlong

  call snap(time2)
  call timing_summary(time2-time1)

  call io%timing%close()
end subroutine dump_timings


!******************************************************************************** 
!
!       Short summary of the timings
!
!******************************************************************************** 
subroutine timing_summary(t)
  real(dp), intent(in) :: t

  write(io%timing%funit,*)
  write(io%timing%funit,102) "=================================================="
  write(io%timing%funit,*)
  write(io%timing%funit,100) "     Total elapsed time         =", write_time(root_timer%time)
  write(io%timing%funit,100) "     Total cpu time             =", write_time(root_timer%cputime)
  write(io%timing%funit,101) "     Ratio cpu / elapsed time   =", root_timer%cputime/root_timer%time
  write(io%timing%funit,*)
  write(io%timing%funit,100) "     Time for measuring timers  =", write_time(self_timer%time)
  write(io%timing%funit,100) "     Time for writing timings   =", write_time(t)
  write(io%timing%funit,102) "=================================================="
  write(io%timing%funit,*)
  write(io%timing%funit,*) starlong 

  100 format(1x,t20,a,t55,a)
  101 format(1x,t20,a,t55,f8.2)
  102 format(1x,t20,a)
end subroutine timing_summary


!******************************************************************************** 
!
!       writes timers structured - called from wrapper function write_
!
!               1. Looks for child and calls routine for ptr%child
!               2. Looks for brother and calls routine for ptr%next
!
!******************************************************************************** 
recursive subroutine write_graph_timers(ptr, parent)
  type(timer), pointer            :: ptr
  type(timer), pointer            :: parent

  call write_graph_timer(ptr, parent)
  if (associated(ptr%child)) call write_graph_timers(ptr%child, ptr)
  if (associated(ptr%next))  call write_graph_timers(ptr%next, parent)
end subroutine write_graph_timers


!******************************************************************************** 
!
!       Displays one particular timer for structured timer output:  
!
!               identation      name     #calls         total time
!               Identation is important to keep structure transparent
!
!******************************************************************************** 
subroutine write_graph_timer(ptr, parent)
  type(timer), pointer  :: ptr
  type(timer), pointer  :: parent
  !local variables
  character(len=2)      :: space = "  "
  character(len=50)     :: name_, time, time_self
  integer               :: ncalls
  
  name_  = repeat(space, ptr%nest_level) // trim(ptr%name_)
  ncalls = ptr%ncalls
  time = repeat(space, ptr%nest_level) // adjustl(write_time(ptr%time))
  time_self = repeat(space, ptr%nest_level) // adjustl(write_time(ptr%time_self))

  write(io%timing%funit,100) trim(name_), ncalls, time_self, time
  100 format(1x,a,t45,i10,t55,a,t85,a)
end subroutine write_graph_timer


!******************************************************************************** 
!
!       Header for graph timers
!
!******************************************************************************** 
subroutine graph_timer_header
  write(io%timing%funit,*)
  write(io%timing%funit,*)
  write(io%timing%funit,*) "     Time Profiling - graph - elapsed time:"
  write(io%timing%funit,*) starshort
  write(io%timing%funit,100) "       routine Name      ", "    ncalls", & 
                             "  self time                  ", "  acc time                   "
  write(io%timing%funit,100) "       ============      ", "    ======", & 
                             "  =========                  ", "  ========                   "
  100 format(1x,a,t45,a,t55,a,t85,a)
end subroutine graph_timer_header


!******************************************************************************** 
!
!       writes percentages structured
!
!               1. Looks for child and calls routine for ptr%child
!               2. Looks for brother and calls routine for ptr%next
!
!******************************************************************************** 
recursive subroutine write_graph_fracs(ptr, parent)
  type(timer), pointer            :: ptr
  type(timer), pointer            :: parent

  call write_graph_frac(ptr, parent)
  if (associated(ptr%child)) call write_graph_fracs(ptr%child, ptr)
  if (associated(ptr%next))  call write_graph_fracs(ptr%next, parent)
end subroutine write_graph_fracs


!******************************************************************************** 
!
!       Displays percentages of one partiuclar timer:  
!
!         identation      name     #calls    self %    total %    
!
!******************************************************************************** 
subroutine write_graph_frac(ptr, parent)
  type(timer), pointer  :: ptr
  type(timer), pointer  :: parent
  !local variables
  character(len=2)      :: space = "  "
  character(len=50)     :: name_, time, time_self
  character(len=50)     :: perc, perc_self
  integer               :: ncalls
  
  name_  = repeat(space, ptr%nest_level) // trim(ptr%name_)
  ncalls = ptr%ncalls
  perc = repeat(space, ptr%nest_level) // adjustl(write_perc(ptr%time/parent%time))
  perc_self = repeat(space, ptr%nest_level) // adjustl(write_perc(ptr%time_self/ptr%time))

  write(io%timing%funit,100) trim(name_), ncalls, perc_self, perc
  100 format(1x,a,t45,i10,t55,a,t85,a)
end subroutine write_graph_frac


!******************************************************************************** 
!
!       Header for graph fracs
!
!******************************************************************************** 
subroutine graph_frac_header
  write(io%timing%funit,*)
  write(io%timing%funit,*)
  write(io%timing%funit,*) "     Time Profiling - graph - percentage :"
  write(io%timing%funit,*) starshort
  write(io%timing%funit,100) "       routine Name      ", "    ncalls", & 
                             "  self time %                 ", "  acc time %                  "
  write(io%timing%funit,100) "       ============      ", "    ======", & 
                             "  ===========                 ", "  ==========                  "
  100 format(1x,a,t45,a,t55,a,t85,a)
end subroutine graph_frac_header


!******************************************************************************** 
!
!       Runs over timers and writes flat timers
!
!******************************************************************************** 
subroutine write_flat_timers(head)
  type(timer), pointer            :: head
  !local variables
  type(timer), pointer            :: ptr
  integer                         :: count_
  
  ptr => head
  count_ = 1
  call write_flat_timer(count_, ptr)

  do
    if (associated(ptr%next)) then
      ptr => ptr%next
      count_ = count_ + 1
      call write_flat_timer(count_, ptr)
    else
      exit
    end if
  end do
end subroutine write_flat_timers


!******************************************************************************** 
!
!       Displays one particular timer for flat timer output:  
!
!       counter    name    #calls    self time    % self    acc time    % accu
!
!                                  acc time/call
!
!******************************************************************************** 
subroutine write_flat_timer(count_, ptr)
  integer, intent(in)   :: count_
  type(timer), pointer  :: ptr
  !local variables
  character(len=50)     :: name_, time, time_self, perc, perc_self, time_rel
  integer               :: ncalls
  
  name_ = trim(ptr%name_)
  ncalls = ptr%ncalls
  time = trim(write_time(ptr%time))
  time_self = trim(write_time(ptr%time_self))
  time_rel = trim(write_time(ptr%time/ptr%ncalls))
  perc = trim(write_perc(ptr%time/root_timer%time))
  perc_self = trim(write_perc(ptr%time_self/root_timer%time))

  write(io%timing%funit,100) count_, name_, ncalls, time_self, perc_self, &
                             time, perc, time_rel
  100 format(1x,i5,t10,a,t40,i10,t55,a,t70,a,t80,a,t95,a,t105,a)
end subroutine write_flat_timer


!******************************************************************************** 
!
!       Header for flat timers
!
!******************************************************************************** 
subroutine flat_timer_header
  write(io%timing%funit,*)
  write(io%timing%funit,*)
  write(io%timing%funit,*) "     Time Profiling - flat - elapsed time:"
  write(io%timing%funit,*) starshort

  write(io%timing%funit,100) "  cnt  ", "       routine name           ", "    ncalls     ", "   self time   ", "  % self  ", &
                             "  acc time      ", "  % acc   ", "  atime/call   "
  write(io%timing%funit,100) "  ===  ", "       ============           ", "    ======     ", "   =========   ", "  ======  ", & 
                             "  ========      ", "  =====   ", "  ==========   "
  100 format(1x,a,t10,a,t40,a,t55,a,t70,a,t80,a,t95,a,t105,a)
end subroutine flat_timer_header


!******************************************************************************** 
!
!       Runs over timers and writes flops
!
!******************************************************************************** 
subroutine write_flops_items(head)
  type(timer), pointer            :: head
  !local variables
  real(dp)                     :: tot_flops, tot_time
  type(timer), pointer            :: ptr
  integer                         :: count_
  
  ptr => head
  count_ = 0

  tot_flops = 0.0_dp
  tot_time = 0.0_dp

  do
    if (associated(ptr%next)) then
      if (ptr%flops > 0.1_dp) then
        count_ = count_ + 1
        call write_flops_item(count_, ptr)
        tot_flops = tot_flops + ptr%flops
        tot_time = tot_time + ptr%time
      end if
      ptr => ptr%next
    else
      exit
    end if
  end do
  
  write(io%timing%funit,*)
  write(io%timing%funit,100) "average FLOPS/core", trim(write_flops(tot_flops, tot_time))
  100 format(1x,t5,a,t85,a)
end subroutine write_flops_items


!******************************************************************************** 
!
!       Displays flops of the particular timer 
!
!       counter    name    #calls    real time    time/call    FLOPS
!
!******************************************************************************** 
subroutine write_flops_item(count_, ptr)
  integer, intent(in)  :: count_
  type(timer), pointer :: ptr
  !local variables
  character(len=50)    :: name_, time, perc, flops, time_rel
  integer              :: ncalls
  
  name_ = trim(ptr%name_)
  flops = trim(write_flops(ptr%flops, ptr%time))

  ncalls = ptr%ncalls
  time = trim(write_time(ptr%time))
  time_rel = trim(write_time(ptr%time/ptr%ncalls))
  perc = trim(write_perc(ptr%time/root_timer%time))

  write(io%timing%funit,100) count_, name_, ncalls, time, time_rel, flops 
  100 format(1x,i5,t10,a,t40,i10,t55,a,t70,a,t85,a)
end subroutine write_flops_item


!******************************************************************************** 
!
!       Header for FLOPS writing
!
!******************************************************************************** 
subroutine flops_header
  write(io%timing%funit,*)
  write(io%timing%funit,*)
  write(io%timing%funit,*) "     Time Profiling - FLOPS report  "
  write(io%timing%funit,*) starshort
  write(io%timing%funit,100) "  cnt  ", "  routine name      ", "    ncalls", & 
                             " acc time      ", " atime/call ", "   FLOPS/core  " 
  write(io%timing%funit,100) "  ===  ", "  ============      ", "    ======", & 
                             " ========      ", " ========== ", "   ==========  " 
  100 format(1x,a,t10,a,t40,a,t55,a,t70,a,t85,a)
end subroutine flops_header


!******************************************************************************** 
!
!       Sum timers over all processes to obtain total CPU times
!
!******************************************************************************** 
recursive subroutine reduce_timers(ptr)
  type(timer), pointer :: ptr

  call reduce_timer(ptr)
  if (associated(ptr%child)) call reduce_timers(ptr%child)
  if (associated(ptr%next))  call reduce_timers(ptr%next)
end subroutine reduce_timers


!******************************************************************************** 
!
!       Reduce specific timers - obtains total cpu times and total flops
!
!******************************************************************************** 
subroutine reduce_timer(ptr)
  type(timer), pointer :: ptr

  call comm_world%mpisum(ptr%ncalls_total)
  call comm_world%mpisum(ptr%cputime)
  call comm_world%mpisum(ptr%flops_total)
end subroutine reduce_timer


!******************************************************************************** 
!
!       Obtain self times - subtract times spent in child routines
!
!******************************************************************************** 
recursive subroutine get_self_timers(ptr)
  type(timer), pointer :: ptr

  call get_self_timer(ptr)
  if (associated(ptr%child)) call get_self_timers(ptr%child)
  if (associated(ptr%next))  call get_self_timers(ptr%next)
end subroutine get_self_timers


!******************************************************************************** 
!
!       Obtain self time of the specific timer
!
!******************************************************************************** 
subroutine get_self_timer(ptr)
  type(timer), pointer :: ptr
  !local variables
  type(timer), pointer :: child
  real(dp)          :: time_self, cputime_self
  
  time_self = 0.0_dp
  cputime_self = 0.0_dp

  child => ptr%child
  do 
    if (associated(child)) then
      time_self = time_self + child%time
      cputime_self = cputime_self + child%cputime
      child => child%next
    else
      exit
    end if
  end do

  ptr%time_self = ptr%time - time_self
  ptr%cputime_self = ptr%cputime - cputime_self

  if (ptr%name_ == "QMCFORT") then
    ptr%cputime_self = ptr%time_self
  end if
end subroutine get_self_timer


!******************************************************************************** 
!
!       Transforms graph structured timers into flat timers
!
!******************************************************************************** 
recursive subroutine create_flat_timers(ptr, head_flat)
  type(timer), pointer :: ptr
  type(timer), pointer :: head_flat
  
  call check_flat_timer(ptr, head_flat)
  if (associated(ptr%child)) call create_flat_timers(ptr%child, head_flat)
  if (associated(ptr%next))  call create_flat_timers(ptr%next, head_flat)
end subroutine create_flat_timers


!******************************************************************************** 
!
!       Adds particular timer from structured to flat format
!
!******************************************************************************** 
subroutine check_flat_timer(ptr, head_flat)
  type(timer), pointer :: ptr
  type(timer), pointer :: head_flat
  !local variables
  type(timer), pointer :: ptr2
  
  if (.not. associated(head_flat)) then
    allocate(head_flat)
    head_flat%name_ = trim(ptr%name_)
    head_flat%ncalls = ptr%ncalls
    head_flat%ncalls_total = ptr%ncalls_total
    head_flat%time = ptr%time
    head_flat%time_self = ptr%time_self
    head_flat%cputime = ptr%cputime
    head_flat%cputime_self = ptr%cputime_self
    head_flat%flops = ptr%flops
    head_flat%flops_total = ptr%flops_total
    head_flat%next => null()
  else
    ptr2 => head_flat
    do
      if (trim(ptr%name_) == trim(ptr2%name_)) then
        ptr2%ncalls = ptr2%ncalls + ptr%ncalls
        ptr2%ncalls_total = ptr2%ncalls_total + ptr%ncalls_total
        ptr2%time = ptr2%time + ptr%time
        ptr2%time_self = ptr2%time_self + ptr%time_self
        ptr2%cputime = ptr2%cputime + ptr%cputime
        ptr2%cputime_self = ptr2%cputime_self + ptr%cputime_self
        ptr2%flops = ptr2%flops + ptr%flops
        ptr2%flops_total = ptr2%flops_total + ptr%flops_total
        exit
      else
        if (associated(ptr2%next)) then
          ptr2 => ptr2%next
        else
          allocate(ptr2%next)
          ptr2 => ptr2%next
          ptr2%name_ = trim(ptr%name_)
          ptr2%ncalls = ptr%ncalls
          ptr2%ncalls_total = ptr%ncalls_total
          ptr2%time = ptr%time
          ptr2%time_self = ptr%time_self
          ptr2%cputime = ptr%cputime
          ptr2%cputime_self = ptr%cputime_self
          ptr2%flops = ptr%flops
          ptr2%flops_total = ptr%flops_total
          ptr2%next => null()
          exit
        end if
      end if
    end do
  end if
end subroutine check_flat_timer


!******************************************************************************** 
!
!       Sort timers according to timer%time_self
!
!******************************************************************************** 
subroutine sort_timers(head)
  type(timer), pointer :: head
  !local variables
  type(timer), pointer :: ptr1, ptr2

  ptr1 => head

  do
    if (associated(ptr1%next)) then
      ptr2 => ptr1%next
    else
      exit
    end if

    do
      if (ptr2%time_self > ptr1%time_self) then
        call swap_timers(ptr1, ptr2)
      else
        if (associated(ptr2%next)) then
          ptr2 => ptr2%next
        else
          exit
        end if
      end if
    end do

    ptr1 => ptr1%next
  end do
end subroutine sort_timers


!******************************************************************************** 
!
!       Swap two timer pointers
!
!******************************************************************************** 
subroutine swap_timers(ptr1, ptr2)
  type(timer), pointer :: ptr1, ptr2
  !local variables
  type(timer)          :: temp

  temp%name_ = ptr1%name_
  temp%ncalls = ptr1%ncalls
  temp%ncalls_total = ptr1%ncalls_total
  temp%nest_level = ptr1%nest_level
  temp%time = ptr1%time
  temp%time_self = ptr1%time_self
  temp%cputime = ptr1%cputime
  temp%cputime_self = ptr1%cputime_self
  temp%flops = ptr1%flops
  temp%flops_total = ptr1%flops_total

  ptr1%name_ = ptr2%name_
  ptr1%ncalls = ptr2%ncalls
  ptr1%ncalls_total = ptr2%ncalls_total
  ptr1%nest_level = ptr2%nest_level
  ptr1%time = ptr2%time
  ptr1%time_self = ptr2%time_self
  ptr1%cputime = ptr2%cputime
  ptr1%cputime_self = ptr2%cputime_self
  ptr1%flops = ptr2%flops
  ptr1%flops_total = ptr2%flops_total

  ptr2%name_ = temp%name_
  ptr2%ncalls = temp%ncalls
  ptr2%ncalls_total = temp%ncalls_total
  ptr2%nest_level = temp%nest_level
  ptr2%time = temp%time
  ptr2%time_self = temp%time_self
  ptr2%cputime = temp%cputime
  ptr2%cputime_self = temp%cputime_self
  ptr2%flops = temp%flops
  ptr2%flops_total = temp%flops_total
end subroutine swap_timers


!******************************************************************************** 
!
!       Manualy set self_timer, since it is not called by start/end_profiling
!
!******************************************************************************** 
subroutine decorate_self_timer
  self_timer%ncalls = 1
  self_timer%ncalls_total = 1
  self_timer%time_self = self_timer%time
  self_timer%cputime_self = self_timer%cputime
end subroutine decorate_self_timer


!******************************************************************************** 
!
!       Count elapsed time
!
!******************************************************************************** 
subroutine snap(time)
  use iso_fortran_env, only: int64
  real(dp), intent(out) :: time
  !logical variables
  integer(int64)           :: time_count, clock_rate
 
#ifdef MPI
  time = MPI_Wtime()
#else
  call system_clock(time_count, clock_rate)
  time = real(time_count, dp) / real(clock_rate, dp)
#endif
end subroutine snap


!******************************************************************************** 
!
!       Count CPU time
!
!******************************************************************************** 
subroutine snap_cpu(time)
  real(dp), intent(out) :: time
 
   call cpu_time(time)
   !call snap(time)
end subroutine snap_cpu

end module profiling
