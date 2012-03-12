!==============================================================================!
! MODULE: particle
!
!> @author Bryan Herman
!>
!> @brief Contains information about the particle that is transporting
!==============================================================================!

module particle

  implicit none
  private
  public :: init_particle

  type, public :: particle_type

    real(8) :: E       ! particle's energy
    logical :: alive   ! am i alive?
    integer :: region  ! material location
    integer :: isoidx  ! isotope index in region
    integer :: reactid ! reaction id

  end type particle_type

contains

!===============================================================================
! INIT_PARTICLE
!> @brief routine to initialize a particle
!===============================================================================

  subroutine init_particle(this)

    ! formal variables
    type(particle_type) :: this ! a particle

    ! initialize
    this%E = 0.0_8
    this%alive = .true.
    this%region = 1
    this%isoidx = 0
    this%reactid = 0

  end subroutine init_particle

end module particle
