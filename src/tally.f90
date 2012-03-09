!==============================================================================!
! MODULE: tally 
!
!> @author Bryan Herman
!>
!> @brief Contains information about tallying quantities 
!==============================================================================!

module tally

  implicit none
  private
  public :: setup_tallies,add_to_tally,bank_tally

  type, public :: tally_type

    real(8), allocatable :: E(:)      ! user defined energy structure
    real(8), allocatable :: val(:)    ! the temporary value
    real(8), allocatable :: sum(:)    ! the sum for the mean and var
    real(8), allocatable :: sum_sq(:) ! the sum for the variable
    logical :: flux_tally             ! is this the flux tally
    integer :: nbins                  ! number of tally regions
    real(8) :: width                  ! the uniform width
    real(8) :: emax                   ! max e
    real(8) :: emin                   ! min e

  end type tally_type

contains

!===============================================================================
! SETUP_TALLIES
!> @brief routine to initialize all tallies 
!===============================================================================

  subroutine setup_tallies(this,n,emax,emin)

    ! formal variables
    integer          :: n       ! number of tally arrays
    type(tally_type) :: this(n) ! array of tallies
    real(8)          :: emax    ! max e
    real(8)          :: emin    ! min e

    ! set up automatic flux tally
    this(1)%flux_tally = .true.
    this(1)%nbins = 1000
    this(1)%emax = emax
    this(1)%emin = emin
    this(1)%width = (log10(emax) - log10(emin))/dble(this(1)%nbins)

    ! preallocate vectors
    if(.not.allocated(this(1)%val)) allocate(this(1)%val(1000)) 
    if(.not.allocated(this(1)%sum)) allocate(this(1)%sum(1000))
    if(.not.allocated(this(1)%sum_sq)) allocate(this(1)%sum_sq(1000)) 

    ! zero out tallies
    this(1)%val = 0.0_8
    this(1)%sum = 0.0_8
    this(1)%sum_sq = 0.0_8

  end subroutine setup_tallies

!===============================================================================
! ADD_TO_TALLY
!> @brief routine to add quantities during transport of a particle 
!===============================================================================

  subroutine add_to_tally(this,n,totxs,E)

    ! formal variables
    integer          :: n        ! size of tallies
    type(tally_type) :: this(n)  ! array of tallies
    real(8)          :: totxs    ! totalxs
    real(8)          :: E        ! neutron energy

    ! local variables
    integer :: i    ! iteration counter
    integer :: idx  ! index in tally grid

    ! begin loop over tallies
    do i = 1,n

      ! use uniform grid sampling if flux tally
      if (this(i)%flux_tally) then

        ! calculate index
        idx = ceiling((log10(E) - log10(this(i)%emin))/this(i)%width)

      else
        write(*,*) 'Something is wrong right now'
        stop
      end if

      ! add to tally
      this(i)%val(idx) = this(i)%val(idx) + 1.0_8/totxs

    end do

  end subroutine add_to_tally

!===============================================================================
! BANK_TALLY
!> @brief routine to bank a histories tallies
!===============================================================================

  subroutine bank_tally(this,n)

    ! formal variables
    integer          :: n       ! size of tallies
    type(tally_type) :: this(n) ! array of tallies 

    ! local variables
    integer :: i    ! iteration counter

    ! begin loop over tallies
    do i = 1,n

      ! record to sums
      this(i)%sum    = this(i)%sum    + this(i)%val
      this(i)%sum_sq = this(i)%sum_sq + this(i)%val**2

    end do

 end subroutine bank_tally

end module tally
