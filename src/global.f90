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
  use timing,    only: Timer

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
  integer :: eidx      ! energy index for cross sections
  integer :: n_tallies ! number of tallies

  ! set max and min energy
  real(8) :: emin = 1e-11_8
  real(8) :: emax = 20.0_8

  ! kT value base on 300K
  real(8) :: kT = 8.6173324e-5_8*300*1.0e-6_8

  ! timers
  type(Timer) :: time_init
  type(Timer) :: time_run

contains

!==============================================================================
! ALLOCATE_PROBLEM
!> @brief allocates global variables for calculation
!==============================================================================

  subroutine allocate_problem(n)

    ! formal variables
    integer :: n ! size of tallies

    ! allocate tallies
    if (.not.allocated(tal)) allocate(tal(n))

  end subroutine allocate_problem

!===============================================================================
! DEALLOCATE_PROBLEM
!> @brief deallocates global variables
!===============================================================================

  subroutine deallocate_problem()

    use materials, only: deallocate_material
    use tally,     only: deallocate_tally

    ! local variables
    integer :: i  ! loop counter

    ! deallocate material
    call deallocate_material(mat)

    ! deallocate within tallies
    do i = 1,n_tallies

      ! deallocate tally
      call deallocate_tally(tal(i))

    end do

    ! deallocate tally variable
    if (allocated(tal)) deallocate(tal)

  end subroutine deallocate_problem

!===============================================================================
! ADD_TO_TALLIES
!> @brief routine that adds temporary value to tallies
!===============================================================================

  subroutine add_to_tallies()

    use tally, only: add_to_tally

    ! local variables
    integer :: i            ! loop counter
    real(8) :: fact = 1.0_8 ! multiplier factor
    real(8) :: totxs        ! total macroscopic xs of material

    ! compute macroscopic cross section
    totxs = sum(mat%totalxs(eidx,:))

    ! begin loop over tallies
    do i = 1,n_tallies

      ! set multiplier
      select case(tal(i)%react_type)

        ! flux only
        case(0)
          fact = 1.0_8

        ! absorption
        case(1)
          fact = sum(mat%absorxs(eidx,:))

        ! scattering
        case(2)
          fact = sum(mat%scattxs(eidx,:))

        ! micro capture
        case(3)
          fact = mat%isotopes(tal(i)%isotope)%xs_capt(eidx)
        case DEFAULT
          fact = 1.0_8

      end select

      ! call routine to add tally
      call add_to_tally(tal(i),fact,totxs,neut%E)

    end do

  end subroutine add_to_tallies

!===============================================================================
! BANK_TALLIES
!> @brief routine that record temporary history information in tallies 
!===============================================================================

  subroutine bank_tallies()

    use tally, only: bank_tally

    ! local variables
    integer :: i  ! loop counter

    ! begin loop over tallies
    do i = 1,n_tallies

      ! call routine to bank tally
      call bank_tally(tal(i))

    end do

  end subroutine bank_tallies

end module global 
