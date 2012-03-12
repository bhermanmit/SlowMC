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
 &          set_user_tally,calculate_statistics

  type, public :: tally_type

    real(8), allocatable :: E(:)      ! user defined energy structure
    real(8), allocatable :: val(:,:)    ! the temporary value
    real(8), allocatable :: sum(:,:)    ! the sum for the mean and var
    real(8), allocatable :: sum_sq(:,:) ! the sum for the variable
    real(8), allocatable :: mean(:,:)   ! mean of tallies
    real(8), allocatable :: std(:,:)    ! standard deviation of tallies
    logical :: flux_tally = .false.   ! is this the flux tally
    integer :: nbins                  ! number of tally regions
    real(8) :: width                  ! the uniform width
    real(8) :: emax                   ! max e
    real(8) :: emin                   ! min e
    integer :: react_type             ! reaction type id
    integer :: isotope                ! isotope number

  end type tally_type

contains

!===============================================================================
! SET_USER_TALLY
!> @brief routine to intialize user-defined tallies
!===============================================================================

  subroutine set_user_tally(this,Ebins,n,react_type,isotope,n_materials)

    ! formal variables
    type(tally_type) :: this        ! a tally
    integer          :: n           ! size of Ebins
    integer          :: react_type  ! reaction type
    integer          :: isotope     ! isotope for multiplier
    integer          :: n_materials ! number of material tally regions
    real(8)          :: Ebins(n)    ! vector of energy bins
    
    ! preallocate user-defined energy structure
    if (.not.allocated(this%E)) allocate(this%E(n))

    ! set energy structure
    this%E = Ebins

    ! set reaction type
    this%react_type = react_type

    ! set isotope
    this%isotope = isotope

    ! preallocate vectors
    if(.not.allocated(this%val)) allocate(this%val(n-1,n_materials))
    if(.not.allocated(this%sum)) allocate(this%sum(n-1,n_materials))
    if(.not.allocated(this%sum_sq)) allocate(this%sum_sq(n-1,n_materials))

    ! preallocate mean and stdev
    if (.not.allocated(this%mean)) allocate(this%mean(size(this%sum),          &
   &                                                               n_materials))
    if (.not.allocated(this%std))  allocate(this%std(size(this%sum),           &
   &                                                               n_materials))

    ! zero out tallies
    this%val = 0.0_8
    this%sum = 0.0_8
    this%sum_sq = 0.0_8

  end subroutine set_user_tally

!===============================================================================
! SET_SPECTRUM_TALLY
!> @brief routine to initialize all tallies 
!===============================================================================

  subroutine set_spectrum_tally(this,emax,emin,n_materials)

    ! formal variables
    type(tally_type) :: this         ! a tally
    integer          :: n_materials  ! number of materials
    real(8)          :: emax         ! max e
    real(8)          :: emin         ! min e

    ! set up automatic flux tally
    this%flux_tally = .true.
    this%nbins = 1000
    this%emax = emax
    this%emin = emin
    this%width = (log10(emax) - log10(emin))/dble(this%nbins)

    ! preallocate vectors
    if(.not.allocated(this%val)) allocate(this%val(1000,n_materials))
    if(.not.allocated(this%sum)) allocate(this%sum(1000,n_materials))
    if(.not.allocated(this%sum_sq)) allocate(this%sum_sq(1000,n_materials))

    ! preallocate mean and stdev
    if (.not.allocated(this%mean)) allocate(this%mean(size(this%sum),          &
   &                                                               n_materials))
    if (.not.allocated(this%std))  allocate(this%std(size(this%sum),           &
   &                                                               n_materials))

    ! zero out tallies
    this%val = 0.0_8
    this%sum = 0.0_8
    this%sum_sq = 0.0_8

  end subroutine set_spectrum_tally

!===============================================================================
! ADD_TO_TALLY
!> @brief routine to add quantities during transport of a particle 
!===============================================================================

  subroutine add_to_tally(this,fact,totxs,E,region)

    ! formal variables
    type(tally_type) :: this     ! a tally 
    integer          :: region   ! region id
    real(8)          :: fact     ! multiplier for tally
    real(8)          :: totxs    ! totalxs
    real(8)          :: E        ! neutron energy

    ! local variables
    integer :: i      ! iteration counter
    integer :: idx=0  ! index in tally grid

    ! use uniform grid sampling if flux tally
    if (this%flux_tally) then

      ! calculate index
      idx = ceiling((log10(E) - log10(this%emin))/this%width)

    else

      ! check for output bounds
      if (E < minval(this%E) .or. E > maxval(this%E)) return

      ! begin loop around energy vector to get index
      do i = 1,size(this%E)
        if (E < this%E(i)) then
          idx = i - 1
          exit
        end if
      end do

    end if

    ! add to tally
    if (idx /= 0) this%val(idx,region) = this%val(idx,region) + fact/totxs

  end subroutine add_to_tally

!===============================================================================
! BANK_TALLY
!> @brief routine to bank a histories tallies
!===============================================================================

  subroutine bank_tally(this)

    ! formal variables
    type(tally_type) :: this ! a tally

    ! record to sums
    this%sum    = this%sum    + this%val
    this%sum_sq = this%sum_sq + this%val**2

    ! zero out temp value
    this%val = 0.0_8

  end subroutine bank_tally

!===============================================================================
! CALCULATE_STATISTICS
!> @brief routine to compute mean and standard deviation of tallies
!===============================================================================

  subroutine calculate_statistics(this,n)

    ! formal variables
    type(tally_type) :: this ! a tally
    integer          :: n    ! number of histories run

    ! compute mean
    this%mean =  this%sum / dble(n)

    ! compute standard deviation
    this%std = sqrt((this%sum_sq/dble(n) - this%mean**2)/dble(n))

  end subroutine calculate_statistics

!===============================================================================
! DEALLOCATE_TALLY
!> @brief routine to deallocate tally types
!===============================================================================

  subroutine deallocate_tally(this)

    ! formal variables
    type(tally_type) :: this ! a tally

    ! deallocate all
    if (allocated(this%E)) deallocate(this%E)
    if (allocated(this%val)) deallocate(this%val)
    if (allocated(this%sum)) deallocate(this%sum)
    if (allocated(this%sum_sq)) deallocate(this%sum_sq)

  end subroutine deallocate_tally

end module tally
