program main
!==============================================================================!
!> @mainpage SlowMC: Slowing Down Monte Carlo 
!>
!> @section Overview
!>
!> This program solves the slowing down neutron transport equation in either
!> infinite medium or effective two-region collision probability theory. It 
!> models parts of the same physics performed by the NJOY data processing code.
!> This code is for strictly academic purposes and allows the user to see the
!> relative impact of physics in the generation of multigroup cross sections
!> and on flux spectra. This code currently uses the following external 
!> libraries:
!>  - HDF5 v1.8.# 
!>
!> The package HDF5 can be downloaded from http://www.hdfgroup.org/HDF5/ 
!>
!> @section Compiling
!>
!> Compiling is as easy as running the Makefile with:
!>
!> @verbatim
!>   make xml-fortran
!>   make slowmc
!> @endverbatim
!>
!> @section Running
!>
!> To run SlowMC, execute the following:
!>
!> @verbatim
!>   slowmc
!> @endverbatim
!>
!==============================================================================!

  implicit none

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
!> @brief high level routine for intializing problem 
!===============================================================================

  subroutine initialize()

    use hdf5
    use global,    only: seed,allocate_problem,mat,tal,emax,emin
    use input,     only: read_input
    use materials, only: compute_macroxs
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

    ! precompute macroscopic cross section of materials
    call compute_macroxs(mat)

    ! set up tally structure
    call setup_tallies(tal,size(tal),emax,emin)

  end subroutine initialize

!===============================================================================
! RUN_PROBLEM
!> @brief main routine for executing the transport calculation
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

      ! print update to user
      if (mod(i,nhistories/10) == 0) then
        write(*,'(/A,1X,I0,1X,A)') 'Simulated',i,'neutrons...'
      end if

    end do

do i = 1,size(tal(1)%sum)
  write(101,*) tal(1)%sum(i)
end do
 
  end subroutine run_problem

!===============================================================================
! FINALIZE
!> @brief routine that finalizes the problem
!===============================================================================

  subroutine finalize()

    use global, only: deallocate_problem
    use hdf5

    ! local variables
    integer :: error ! hdf5 error

    ! write output

    ! deallocate problem
    call deallocate_problem()

    ! close the fortran interface
    call h5close_f(error)

  end subroutine finalize

end program main
