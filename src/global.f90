!==============================================================================!
! MODULE: global 
!
!> @author Bryan Herman
!>
!> @brief Contains all of the global variables
!==============================================================================!

module global

  use materials, only: material_type
  use particle,  only: particle_type
  use tally,     only: tally_type

  implicit none
  save

  ! version information
  integer :: VERSION_MAJOR   = 0
  integer :: VERSION_MINOR   = 1
  integer :: VERSION_RELEASE = 1

  ! list all types
  type(material_type)           :: mat
  type(particle_type)           :: neut
  type(tally_type), allocatable :: tal(:)

  ! list history input information
  integer :: nhistories
  integer :: seed
  integer :: source_type

  ! list global vars that are set during run
  integer :: eidx   ! energy index for cross sections

  ! set max and min energy
  real(8) :: emin = 1e-11_8
  real(8) :: emax = 20.0_8

  ! kT value base on 300K
  real(8) :: kT = 8.6173324e-5_8*300*1.0e-6_8

contains

!==============================================================================
! ALLOCATE_PROBLEM
!> @brief allocates global variables for calculation
!==============================================================================

  subroutine allocate_problem()

    ! allocate tallies
    if(.not.allocated(tal)) allocate(tal(1))

  end subroutine allocate_problem

end module global 
