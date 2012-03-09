program main

  implicit none

  ! print header information

  ! initialize problem
  call initialize()

  ! run problem
  call run_problem()

  ! finalize problem 
  call finalize()

  ! terminate program
  stop

contains

!===============================================================================
! INTIALIZE
! Doxygen comment
!===============================================================================

  subroutine initialize()

    use hdf5
    use global,    only: seed,allocate_problem,mat,tal,emax,emin
    use input,     only: read_input
    use materials, only: compute_totxs
    use output,    only: print_heading
    use tally,     only: setup_tallies

    ! local variables
    integer :: error ! hdf5 error
    real(8) :: rn    ! initial random number

    ! initialize the fortran hdf5 interface
    call h5open_f(error)

    ! print heading information
    call print_heading()

    ! read input
    call read_input()

    ! initalize random number generator
    rn = rand(seed)

    ! allocate problem
    call allocate_problem()

    ! precompute total cross section of materials
    call compute_totxs(mat)

    ! set up tally structure
    call setup_tallies(tal,size(tal),emax,emin)

  end subroutine initialize

!===============================================================================
! RUN_PROBLEM
! Doxygen comment
!===============================================================================

  subroutine run_problem()

    use global,    only: nhistories,mat,neut,tal,eidx,emin
    use particle,  only: init_particle
    use physics,   only: sample_source,perform_physics,get_eidx
    use tally,     only: add_to_tally,bank_tally

    ! local variables
    integer :: i  ! iteration counter

    ! begin loop over histories
    do i = 1,nhistories

      ! intialize history
      call init_particle(neut)

      ! sample source energy
      call sample_source()

      ! begin transport of neutron
      do while (neut%alive)

        ! call index routine for first tally
        eidx = get_eidx(neut%E)

        ! record collision temp tally
        call add_to_tally(tal,size(tal),sum(mat%totalxs(eidx,:)),neut%E)

        ! perform physics
        call perform_physics()

        ! check for energy cutoff
        if (neut%E < emin) neut%alive = .FALSE.

      end do

      ! neutron is dead if out of transport loop (ecut or absorb) --> bank tally
      call bank_tally(tal,size(tal))

    end do

do i = 1,size(tal(1)%sum)
  write(101,*) tal(1)%sum(i)
end do
 
  end subroutine run_problem

!===============================================================================
! FINALIZE
! Doxygen comment
!===============================================================================

  subroutine finalize()

    use hdf5

    ! local variables
    integer :: error ! hdf5 error

    ! write output

    ! deallocate problem

    ! close the fortran interface
    call h5close_f(error)

  end subroutine finalize

end program main
