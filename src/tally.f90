module tally

  implicit none
  private

  type, public :: tally_type

    real(8), allocatable :: E(:)
    real(8), allocatable :: val(:)
    real(8), allocatable :: sum(:)
    real(8), allocatable :: sum_sq(:)

  end type tally_type

contains

end module tally
