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

    real(8) :: E      ! particle's energy
    logical :: alive  ! am i alive?

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

  end subroutine init_particle

end module particle
