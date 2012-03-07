program main

  implicit none

  ! print header information

  ! initialize problem

  ! run problem

  ! finalize problem 

  ! terminate program
  stop

contains

!===============================================================================
! INTIALIZE
! Doxygen comment
!===============================================================================

  subroutine initialize()

    ! read input

    ! initalize random number generator

    ! allocate problem

    ! precompute total cross section of materials

  end subroutine initialize

!===============================================================================
! RUN_PROBLEM
! Doxygen comment
!===============================================================================

  subroutine run_problem()

    ! begin loop over histories

      ! sample source energy

      ! begin transport of neutron

        ! bank record collision temp tally

        ! perform physics

        ! check for energy cutoff

      ! neutron is dead if out of transport loop (ecut or absorb) --> bank tally
 
  end subroutine run_problem

!===============================================================================
! FINALIZE
! Doxygen comment
!===============================================================================

  subroutine finalize()

    ! write output

    ! deallocate problem

  end subroutine finalize

end program main
