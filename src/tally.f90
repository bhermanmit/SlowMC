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
  public :: set_spectrum_tally,add_to_tally,bank_tally,deallocate_tally,       &
 &          set_user_tally

  type, public :: tally_type

    real(8), allocatable :: E(:)      ! user defined energy structure
    real(8), allocatable :: val(:)    ! the temporary value
    real(8), allocatable :: sum(:)    ! the sum for the mean and var
    real(8), allocatable :: sum_sq(:) ! the sum for the variable
    logical :: flux_tally = .false.   ! is this the flux tally
    integer :: nbins                  ! number of tally regions
    real(8) :: width                  ! the uniform width
    real(8) :: emax                   ! max e
    real(8) :: emin                   ! min e
    integer :: react_type             ! reaction type id

  end type tally_type

contains

!===============================================================================
! SET_USER_TALLY
!> @brief routine to intialize user-defined tallies
!===============================================================================

  subroutine set_user_tally(this,Ebins,n,react_type)

    ! formal variables
    type(tally_type) :: this        ! a tally
    integer          :: n           ! size of Ebins
    integer          :: react_type  ! reaction type
    real(8)          :: Ebins(n)    ! vector of energy bins
    
    ! preallocate user-defined energy structure
    if (.not.allocated(this%E)) allocate(this%E(n))

    ! set energy structure
    this%E = Ebins

    ! set reaction type
    this%react_type = react_type

    ! preallocate vectors
    if(.not.allocated(this%val)) allocate(this%val(n-1))
    if(.not.allocated(this%sum)) allocate(this%sum(n-1))
    if(.not.allocated(this%sum_sq)) allocate(this%sum_sq(n-1))

    ! zero out tallies
    this%val = 0.0_8
    this%sum = 0.0_8
    this%sum_sq = 0.0_8


  end subroutine set_user_tally

!===============================================================================
! SET_SPECTRUM_TALLY
!> @brief routine to initialize all tallies 
!===============================================================================

  subroutine set_spectrum_tally(this,emax,emin)

    ! formal variables
    type(tally_type) :: this    ! array of tallies
    real(8)          :: emax    ! max e
    real(8)          :: emin    ! min e

    ! set up automatic flux tally
    this%flux_tally = .true.
    this%nbins = 1000
    this%emax = emax
    this%emin = emin
    this%width = (log10(emax) - log10(emin))/dble(this%nbins)

    ! preallocate vectors
    if(.not.allocated(this%val)) allocate(this%val(1000)) 
    if(.not.allocated(this%sum)) allocate(this%sum(1000))
    if(.not.allocated(this%sum_sq)) allocate(this%sum_sq(1000)) 

    ! zero out tallies
    this%val = 0.0_8
    this%sum = 0.0_8
    this%sum_sq = 0.0_8

  end subroutine set_spectrum_tally

!===============================================================================
! ADD_TO_TALLY
!> @brief routine to add quantities during transport of a particle 
!===============================================================================

  subroutine add_to_tally(this,fact,totxs,E)

    ! formal variables
    type(tally_type) :: this     ! array of tallies
    real(8)          :: fact     ! multiplier for tally
    real(8)          :: totxs    ! totalxs
    real(8)          :: E        ! neutron energy

    ! local variables
    integer :: i    ! iteration counter
    integer :: idx  ! index in tally grid

    ! use uniform grid sampling if flux tally
    if (this%flux_tally) then

      ! calculate index
      idx = ceiling((log10(E) - log10(this%emin))/this%width)

    else

      ! begin loop around energy vector to get index
      do i = 1,size(this%E)
        if (E < this%E(i)) then
          idx = i - 1
          exit
        end if
      end do

    end if

    ! add to tally
    if (idx /= 0) this%val(idx) = this%val(idx) + fact/totxs

  end subroutine add_to_tally

!===============================================================================
! BANK_TALLY
!> @brief routine to bank a histories tallies
!===============================================================================

  subroutine bank_tally(this)

    ! formal variables
    type(tally_type) :: this ! array of tallies 

    ! record to sums
    this%sum    = this%sum    + this%val
    this%sum_sq = this%sum_sq + this%val**2

    ! zero out temp value
    this%val = 0.0_8

  end subroutine bank_tally

!===============================================================================
! DEALLOCATE_TALLY
!> @brief routine to deallocate tally types
!===============================================================================

  subroutine deallocate_tally(this)

    ! formal variables
    type(tally_type) :: this ! array of tallies 

    ! deallocate all
    if (allocated(this%E)) deallocate(this%E)
    if (allocated(this%val)) deallocate(this%val)
    if (allocated(this%sum)) deallocate(this%sum)
    if (allocated(this%sum_sq)) deallocate(this%sum_sq)

  end subroutine deallocate_tally

end module tally
