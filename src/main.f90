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

    use global,  only: nhistories,mat,neut,tal,ntals,eidx
    use physics, only: perform_physics
    use tally,   only: add_to_tally,bank_tally

    ! local variables
    integer :: i  ! iteration counter

    ! begin loop over histories
    do i = 1,nhistories

      ! sample source energy

      ! begin transport of neutron
      do while (neut%alive)

        ! record collision temp tally
        call add_to_tally(tal,ntals,sum(mat%totalxs(eidx,:)))

        ! perform physics
        call perform_physics()

        ! check for energy cutoff
        if (neut%E < 1e-11) neut%alive = .FALSE.

      end do

      ! neutron is dead if out of transport loop (ecut or absorb) --> bank tally
      call bank_tally(tal,ntals)

    end do
 
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
