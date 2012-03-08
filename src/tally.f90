module tally

  implicit none
  private
  public :: add_to_tally,bank_tally

  type, public :: tally_type

    real(8), allocatable :: E(:)      ! user defined energy structure
    real(8), allocatable :: val(:)    ! the temporary value
    real(8), allocatable :: sum(:)    ! the sum for the mean and var
    real(8), allocatable :: sum_sq(:) ! the sum for the variable
    logical :: flux_tally             ! is this the flux tally
    real(8) :: width                  ! the uniform width

  end type tally_type

contains

!===============================================================================
! ADD_TO_TALLY
! Doxygen comment
!===============================================================================

  subroutine add_to_tally(this,n,totxs)

    ! formal variables
    type(tally_type) :: this(n)  ! array of  tallies
    real(8)          :: totxs    ! totalxs
    integer          :: n        ! size of tallies

    ! local variables
    integer :: i    ! iteration counter
    integer :: idx  ! index in tally grid

    ! begin loop over tallies
    do i = 1,n

      ! use uniform grid sampling if flux tally
      if (this(i)%flux_tally) then

        ! calculate index
        idx = ceiling(abs(log10(20.0_8) - log10(1.0e-11_8))/this(i)%width);

      end if

      ! add to tally
      this(i)%val(idx) = this(i)%val(idx) + 1/totxs

    end do

  end subroutine add_to_tally

!===============================================================================
! BANK_TALLY
! Doxygen comment
!===============================================================================

  subroutine bank_tally(this,n)

    ! formal variables
    type(tally_type) :: this(n) ! array of tallies 
    integer          :: n       ! size of tallies

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
